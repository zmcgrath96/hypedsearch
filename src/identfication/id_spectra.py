from src.file_io import mzML
from src.identfication import search
from src.types.database import Database
from src.types.aligner import Aligner
from src.types.objects import Spectrum, MassSequence
from src.scoring import mass_comparisons
from src.spectra.gen_spectra import gen_spectrum

from collections import defaultdict
import math 
from bisect import bisect

############## Constants ##############
TOP_N = 5
#######################################

def make_all_base_mers(database: Database, base_mer: int) -> defaultdict:
    '''
    Create the list of all the base mers from 1 to base_mer with mass, protein, start position information
    
    Inputs:
        database:    (Database) source of the sequences
        base_mer:    (int) the base k-mer to make up to
    Outputs:
        list of MassSequence for all singly, doubly b and y masses
    '''
    allbasemers = defaultdict(list)
    database.set_kmer_size(base_mer)
    database.index()
    md = database.metadata
    for mer in md:
        mer_spec = gen_spectrum(mer)['spectrum']
        for mass in mer_spec:
            mass_key = math.floor(mass)
            allbasemers[mass_key].append(MassSequence(mass, mer))

    return allbasemers

def search_base_kmers(observed: Spectrum, allbasemers: dict, tolerance: float) -> list:
    '''
    Search through all of the base kmers and find those that gave us good hits
    
    Inputs:
        spectrum:    (Spectrum) what to sequence
        allbasemers: (dict of list of MassSequence) all of the basemers made from the function 'make_all_base_mers_hash'
        tolerance:   (float) the ppm tolerance to accept for each mass
    Outputs:
        list of MassSequence for all masses that were in the acceptable range of an observed mass
    '''
    hits = []
    for mass in observed.spectrum:
        tol = mass_comparisons.ppm_opt(mass, tolerance)
        lb_mass = mass - tol
        ub_mass = mass + tol
        lb_mass_key = math.floor(lb_mass)
        ub_mass_key = math.floor(ub_mass)
        
        hits += [x.sequence for x in allbasemers[ub_mass_key] if lb_mass <= x.mass <= ub_mass]
        if lb_mass_key != ub_mass_key:
            hits += [x.sequence for x in allbasemers[ub_mass_key] if lb_mass <= x.mass <= ub_mass]
            
    return list(set(hits))




def id_spectrum(spectrum: Spectrum, database: Database, base_kmers: list, n=TOP_N, min_peptide_len=5, max_peptide_len=20) -> Aligner:
    '''
    Run a scoring and alignment on a specific spectrum 

    Inputs:
        spectrum:   Spectrum namedtuple instance
        database:   Database class instance
        base_kmers: list of top scoring candiates for extension
    kwargs:
        n:                  number of alignments to return for each spectrum. Default=3
        min_peptide_len:     minimum peptide length to consider. Default=5
        max_peptide_len:    maximum peptide length to consider. Default=20
    Outputs:
        Aligner
    '''
    # get the top b and top y scores from the database
    top_b, top_y = search.search_database(spectrum, database, base_kmers, TOP_N)
    # top b and y are dictionaries with rank as key and values are Scores objects
    # try and find some alignment from these scores
    aligner = Aligner(spectrum)
    aligner.add_scores([x for _, x in top_b.items()], 'b')
    aligner.add_scores([x for _, x in top_y.items()], 'y')
    aligner.make_alignments()
        
    return aligner


def id_spectra(spectra_files: list, database_file: str, verbose=True, min_peptide_len=5, max_peptide_len=20, ppm_tolerance=20) -> dict:
    '''
    Run a scoring and alignment on each spectra passed in

    Inputs:
        spectra_files:      (list of strings) of file names of spectra
        database_file:      (string) full path to a .fasta database
    kwargs: 
        verbose:            (bool) whether or not to print messages. Default=True
        min_peptide_len:    (int) minimum length sequence to consider for alignment. Default=5
        max_peptide_len:    (int) maximum length sequence to consider for alignemtn. Default=20
        ppm_tolerance:      (float) tolerance for ppm to include in search. Default=20
    Outputs:
        dict containing the results. All information is keyed by the spectrum file name with scan number appended and the values are Aligner objects
    '''
    verbose and print('Loading database...')
    database = Database(database_file, is_uniprot=True, kmer_size=min_peptide_len)
    verbose and print('Done. Indexing database...')
    base_k_mer_masses = make_all_base_mers(database, min_peptide_len)
    verbose and print('Done.')
    verbose and print('Number of {}-mers found in the database: {}'.format(min_peptide_len, len(database.metadata.keys())))
    results = {}
    # go through all of the mzml files
    for i, spectrum_file in enumerate(spectra_files):
        print('Analyzing spectra file {}/{}[{}%]\n'.format(i + 1, len(spectra_files), int(float(i)/float(len(spectra_files)) * 100)))

        spectra = mzML.read(spectrum_file)
        # go through each spectrum in the mzml file
        for j, spec in enumerate(spectra):
            print('Analyzing spectrum {}/{}[{}%]\r'.format(j + 1, len(spectra), int(float(j)/float(len(spectra)) * 100)), end='')
            # make a Spectrum namedtuple object
            spectrum = Spectrum(spec['spectrum'], spec['abundance'], spec['level'], spec['scan_no'], spec['precursor_mass'], spectrum_file)
            # get the good kmers
            interesting_kmers = search_base_kmers(spectrum, base_k_mer_masses, ppm_tolerance)
            aligned_spectrum = id_spectrum(spectrum, database, interesting_kmers, min_peptide_len=min_peptide_len, max_peptide_len=max_peptide_len)
            entry_name = '{}_{}'.format(spectrum_file, spectrum.scan_number)
            results[entry_name] = aligned_spectrum

    return results
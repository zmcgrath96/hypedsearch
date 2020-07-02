from src.file_io import mzML
from src.identfication.search import search_kmers_hash
from src.identfication.alignment import attempt_alignment
from src.types.database import Database
from src.types.objects import Spectrum, MassSequence, KmerMasses, KmerMassesResults, Alignments, DatabaseEntry
from src.scoring import mass_comparisons
from src.sequence.gen_spectra import gen_spectrum

from collections import defaultdict
import math 
from bisect import bisect

###################################################################################
#                   PRE-SEARCH DATABASE FUNCTIONS
###################################################################################

def build_kmermasses(
    database: Database, 
    min_peptide_len: int, 
    max_peptide_len: int, 
    verbose=False
) -> KmerMasses:
    '''
    Build a KmerMasses object from a database. The entries to the KmerMasses object 
    are dictionaries where keys are integer values of masses and the entries are 
    lists of MassSequence objects. Hashing integer masses lets us more quickly search
    each mass.
    
    Inputs: 
        database:          (Database) object used for reading in proteins, creating the tree, and indexing
        min_peptide_len:   (int) the minimum length peptide to consider. NOTE: this is also the minimum length 
                                 any protein can contribute to a hybrid peptide. 
                                 Example:
                                     protein 1: ABCDEFGHIJK, protein 2: LMNOPQRSTUV
                                     true hybrid sequence: IJK-LMNOPQ
                                     min_peptide_len should be set to 3
        max_peptide_len:   (int) the maximum peptide length to consider
    Outputs:
        KmerMasses object
    '''
    # defaultdicts to hash integer value of masses into 
    bs = defaultdict(list)
    bd = defaultdict(list)
    ys = defaultdict(list)
    yd = defaultdict(list)
    
    # keep track of what kmers we've seen to avoid re-analyzing and 
    # inserting kmers we've seen before
    kmer_tracker = defaultdict(str)

    prot_entry: DatabaseEntry
    
    for i, prot_entry in enumerate(database):

        verbose and print(f'Looking at protein {i + 1}/{len(database)}\r', end='')

        # go through the kmers for the b sequences
        for j in range(len(prot_entry.sequence) - min_peptide_len):

            # make a kmer sequence. Do the max (to generate the kmer spec once) then 
            # just iterate through it
            kmer_len = max_peptide_len if j + max_peptide_len <= len(prot_entry.sequence) \
                else len(prot_entry.sequence) - j
            kmer = prot_entry.sequence[j:j+kmer_len]

            # generate the singly and doubly b spectra
            kmer_spec_b_s = gen_spectrum(kmer, ion='b', charge=1)['spectrum']
            kmer_spec_b_d = gen_spectrum(kmer, ion='b', charge=2)['spectrum']

            # iterate through the spectra and add the entry to the table
            for k in range(min_peptide_len, kmer_len):

                if 'b' in kmer_tracker[kmer[:k]]:
                    continue

                kmer_tracker[kmer[:k]] += 'b'

                bs[math.floor(kmer_spec_b_s[k-1])].append(MassSequence(kmer_spec_b_s[k-1], kmer[:k]))
                bd[math.floor(kmer_spec_b_d[k-1])].append(MassSequence(kmer_spec_b_d[k-1], kmer[:k]))
             
        # go through the kmers for the b sequences
        for j in range(len(prot_entry.sequence) - min_peptide_len):

            # make a kmer sequence. Do the max (to generate the kmer spec once) then 
            # just iterate through it
            kmer_len = max_peptide_len if j + max_peptide_len <= len(prot_entry.sequence) \
                else len(prot_entry.sequence) - j

            kmer = prot_entry.sequence[-j - kmer_len: -j] if j != 0 else prot_entry.sequence[-kmer_len:]

            # generate the singly and doubly b spectra
            kmer_spec_y_s = gen_spectrum(kmer, ion='y', charge=1)['spectrum']
            kmer_spec_y_d = gen_spectrum(kmer, ion='y', charge=2)['spectrum']

            # iterate through the spectra and add the entry to the tayle
            for k in range(min_peptide_len, kmer_len):

                if 'y' in kmer_tracker[kmer[-k:]]:
                    continue

                kmer_tracker[kmer[-k:]] += 'y'

                ys[math.floor(kmer_spec_y_s[k-1])].append(MassSequence(kmer_spec_y_s[k-1], kmer[-k:]))
                yd[math.floor(kmer_spec_y_d[k-1])].append(MassSequence(kmer_spec_y_d[k-1], kmer[-k:]))
 
    return KmerMasses(bs, bd, ys, yd)

###################################################################################
#                   /PRE-SEARCH DATABASE FUNCTIONS
###################################################################################

def id_spectrum(
    spectrum: Spectrum, 
    db: Database, 
    kmermasses: KmerMasses, 
    min_peptide_len: int, 
    n=3, 
    ppm_tolerance=20, 
    precursor_tolerance=1,
    scoring_alg='ibb', 
    DEBUG=False
) -> Alignments:
    '''
    Run an alignemnt in the form of an Amino Acid sequence with a score to the caller for this spectrum.
    The top n results are returned. If no alignments can be created, None is returned.

    Inputs:
        spectrum:               (Spectrum) the spectrum to id
        kmermasses:             (KmerMasses) hash tables with al the necessary dictionaries
    kwargs:
        n:                      (int) number of alignments to return for each spectrum. Default=3
        ppm_tolerance:          (int) the ppm tolerance to allow when searching. Default=20
        precursor_tolerance:    (float) the tolerance to allow when matching precusor masses. Default=1
        scoring_alg:            (str) the name of the scoring algoirhtm to use. Options are 'bb', 'ion', 'ibb'. Default=bb
    Outputs:
        Alignments namedtuple or None
    '''
    # search the mass tables
    bs, bd, ys, yd = search_kmers_hash(spectrum, kmermasses.bs, 20), \
                    search_kmers_hash(spectrum, kmermasses.bd, 20), \
                    search_kmers_hash(spectrum, kmermasses.ys, 20), \
                    search_kmers_hash(spectrum, kmermasses.yd, 20)
    
    # put the results into a structrue
    hits = KmerMassesResults(bs, bd, ys, yd)

    # if we get no hits whatsoever, return None
    if all([len(x) == 0 for x in kmermasses]):
        return None
        
    # attempt alignments
    a = attempt_alignment(
        spectrum, 
        db, 
        hits, 
        min_peptide_len, 
        ppm_tolerance=ppm_tolerance, 
        scoring_alg=scoring_alg, 
        precursor_tolerance=precursor_tolerance, 
        DEBUG=DEBUG
    )
    
    # return the alignments in a structure
    return Alignments(spectrum, a)

def id_spectra(
    spectra_files: list, 
    database_file: str, 
    verbose=True, 
    min_peptide_len=5, 
    max_peptide_len=20, 
    result_count=3, 
    peak_filter=0, 
    relative_abundance_filter=0.0,
    ppm_tolerance=20, 
    precursor_tolerance=1, 
    scoring_alg='ibb', 
    DEBUG=False
) -> dict:
    '''
    Run a scoring and alignment on each spectra passed in and give spectra a sequence of 
    Amino Acids that best describes the spectrum

    Inputs:
        spectra_files:          (list of strings) of file names of spectra
        database_file:          (string) full path to a .fasta database
    kwargs: 
        verbose:                (bool) whether or not to print messages. Default=True
        min_peptide_len:        (int) minimum length sequence to consider for alignment. Default=5
        max_peptide_len:        (int) maximum length sequence to consider for alignemtn. Default=20
        ppm_tolerance:          (int) tolerance for ppm to include in search. Default=20
        precursor_tolerance:    (float) the tolerance to allow when matching precusor masses. Default=1
        scoring_alg:            (str) the scoring algorithm to use. Either 'bb' or 'ion'. Default=bb
    Outputs:
        dict containing the results. 
        All information is keyed by the spectrum file name with scan number appended 
        and the values are list of SequenceAligment objects
    '''

    # load the database into memory
    verbose and print('Loading database...')
    db = Database(database_file,min_len=min_peptide_len, max_len=max_peptide_len , verbose=verbose)
    verbose and print('\nDone.')

    # build the kmermasses namedtuple object once for fast search
    verbose and print('Building hashes for kmers...')
    kmermasses = build_kmermasses(db, min_peptide_len, max_peptide_len, verbose=verbose)
    verbose and print(f'\nDone.')

    # keep track of the results
    results = {}

    # go through all of the spectra files
    for i, spectrum_file in enumerate(spectra_files):

        verbose and print('Analyzing spectra file {}/{}[{}%]\n'.format(i + 1, len(spectra_files), int(float(i)/float(len(spectra_files)) * 100)))

        # load the spectra into memory
        spectra = mzML.read(spectrum_file, peak_filter=peak_filter, relative_abundance_filter=relative_abundance_filter)

        # go through each spectrum in the mzml file
        for j, spec in enumerate(spectra):

            verbose and print('Analyzing spectrum {}/{}[{}%]\r'.format(j + 1, len(spectra), int(float(j)/float(len(spectra)) * 100)), end='')            
            
            # align this spectrum
            aligned_spectrum = id_spectrum(
                spec,
                 db,
                 kmermasses,
                 min_peptide_len,
                 result_count,
                 ppm_tolerance=ppm_tolerance,
                 precursor_tolerance=precursor_tolerance,
                 scoring_alg=scoring_alg,
                DEBUG=DEBUG
            )
            
            # if an alignment cannot be created, continue
            if aligned_spectrum is None:
                continue

            # save the results in the dictionary for now
            entry_name = '{}_{}'.format(spectrum_file, spec.scan_number)
            results[entry_name] = aligned_spectrum

    return results
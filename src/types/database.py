from suffix_tree import Tree
import string
from src.types.objects import KmerMetaData, DatabaseEntry, KmerMasses, MassSequence, Spectrum
from src.sequence.gen_spectra import gen_spectrum
from collections import defaultdict
from pyteomics import fasta

from src.utils import ppm_to_da

import math
    
class Database: 
    '''
    Database

    Container for multiple database operations. Packages 
    '''
    def __init__(
        self, 
        fasta_file_name='', 
        is_uniprot=False, 
        min_len=3, 
        max_len=20, 
        verbose=False
    ) -> None:
        '''
        Init the database with a fasta file name

        Inputs:
        kwargs:
            fasta_file_name:    str name of fasta file to read
            is_uniprot:         bool Set to True if the database is from UnitPort. If True, extra informtion is kept
            min_len:          int size of kmer to use for indexing the database
            verbose:        (bool) extra printing. Default=False

        Outputs: 
            None
        '''
        self.fasta_file: str = fasta_file_name
        self.proteins: dict = self.__read_fasta(fasta_file_name) if '.fasta' in fasta_file_name else {}
        self.min_len: int = min_len
        self.max_len: int = max_len
        self.metadata: dict = None
        self.verbose: bool = verbose
        self.tree: Tree = None
        self.kmer_masses: KmerMasses = None

    ##################### Overloaded operators ################

    def __iter__(self):
        for _, entry in self.proteins.items():
            yield entry

    def __len__(self):
        return len(self.proteins)
    
    ##################### Private Methods #####################

    def __read_fasta(self, fasta_file: str, is_uniprot=False) -> dict:
        '''
        Read proteins into memory from fasta file
        
        Inputs:
            fasta_file: str path to fasta file
        kwargs:
            is_uniprot: bool adds attribute 'human_readable_name' to dictionary if True. Default=False
        Outputs:
            dictionary of Entry class keyed by the protein name 
        '''
        prots = {}
        e: DatabaseEntry

        # split the name on the OS value if it exists
        get_name = lambda name: name[:name.index('OS=')-1] if 'OS=' in name else name

        # go through each entry in the fasta and put it in memory
        for i, entry in enumerate(fasta.read(fasta_file)):

            # take the description without the 'sp' value
            desc = entry.description.split('|')[1:] if '|' in entry.description else entry.description

            # if the id is in the description, take it
            if len(desc) > 1:
                id_ = desc[0]
                name = get_name(desc[1])
                
            # make the id just the number
            else:
                id_ = i
                name = get_name(desc[0])
                
            # get the sequence
            seq = entry.sequence

            # make the entry and add it to prots
            e = DatabaseEntry(name, seq, id_, name)
            prots[name] = e

        return prots

    ##################### Public Methods #####################
    def set_max_len(self, k: int) -> None:
        '''
        Change the max size of the database for indexing

        Inputs:
            k:      (int) new max size 
        Outputs:
            None
        '''
        self.max_len = k
        
    def get_entry_by_name(self, name: str) -> DatabaseEntry:
        '''
        Get a protein entry from the database

        Inputs: 
            name:       str name of the protein to get
        Outputs:
            Entry class of the protein. Empty Entry if not found
        '''
        return self.proteins.get(name, DatabaseEntry())

    def add_entry(
        self, 
        protein_name: str, 
        protein_sequnece: str, 
        protein_id='', 
        human_readable_name=''
    ) -> bool:
        '''
        Add an entry to the database. If it cannot be added to the database, False is returned

        Inputs:
            protein_name:       (str) name of the protein being added
            protein_sequence:   (str) string of amino acids describing the protein
        kwargs:
            protein_id:         (str) identifier of the protein. Default=''
            human_readable_name:(str) name that is easy for people to read. Default=''
        Outputs:
            bool True if added, False if it can't be
        '''
        if protein_name in self.proteins:
            return False

        e = DatabaseEntry(protein_name, protein_sequnece, protein_id, human_readable_name)
        self.proteins[protein_name] = e 

        # add it to the tree
        if not self.tree:
            self.tree = Tree()

        self.tree.add(protein_name, protein_sequnece)

        return True

    def build(
        self, 
        min_peptide_len: int, 
        max_peptide_len: int, 
        # missed_cleavages: int,
        # digest: str,
        verbose=False
    ) -> KmerMasses:
        '''
        Build a KmerMasses object from a database. The entries to the KmerMasses object 
        are dictionaries where keys are integer values of masses and the entries are 
        lists of MassSequence objects. Hashing integer masses lets us more quickly search
        each mass.
        
        Inputs: 
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
        self.tree = Tree()

        # defaultdicts to hash integer value of masses into 
        bs = defaultdict(list)
        bd = defaultdict(list)
        ys = defaultdict(list)
        yd = defaultdict(list)
        
        # keep track of what kmers we've seen to avoid re-analyzing and 
        # inserting kmers we've seen before
        kmer_tracker = defaultdict(str)

        
        for i, prot_name in enumerate(self.proteins):

            verbose and print(f'Looking at protein {i + 1}/{len(self.proteins)}\r', end='')

            prot_entry: DatabaseEntry = self.proteins[prot_name]

            # add the protein to the tree
            self.tree.add(prot_entry.name, prot_entry.sequence)

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
    
        self.kmer_masses = KmerMasses(bs, bd, ys, yd)


    def search(self, observed: Spectrum, kmers: str, tolerance: float) -> list:
        '''
        Search through all masses and saved kmers to find masses that are within our tolerance
        
        Inputs:
            spectrum:    (Spectrum) what to sequence
            allbasemers: (dict of list of MassSequence) all of the basemers made from the function 'make_all_base_mers_hash'
            tolerance:   (float) the ppm tolerance to accept for each mass
        Outputs:
            list of MassSequence for all masses that were in the acceptable range of an observed mass
        '''
        kmer_dict = self.kmer_masses.bs if kmers == 'bs' else (self.kmer_masses.bd if kmers == 'bd' \
                    else (self.kmer_masses.ys if kmers == 'ys' else self.kmer_masses.yd))
        hits = []
        for mass in observed.spectrum:
            tol = ppm_to_da(mass, tolerance)
            lb_mass = mass - tol
            ub_mass = mass + tol
            lb_mass_key = math.floor(lb_mass)
            ub_mass_key = math.floor(ub_mass)
            hits += [x.sequence for x in kmer_dict[ub_mass_key] if lb_mass <= x.mass <= ub_mass]
            if lb_mass_key != ub_mass_key:
                hits += [x.sequence for x in kmer_dict[ub_mass_key] if lb_mass <= x.mass <= ub_mass]
                
        return hits
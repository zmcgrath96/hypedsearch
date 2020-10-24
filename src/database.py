from pyteomics import fasta
from collections import namedtuple, defaultdict
from math import ceil
from typing import Iterable

from src import gen_spectra
from src.utils import predicted_len, ppm_to_da, overlap_intervals

import array as arr

BATCH_SIZE = 300

class Database:
    '''
    Class containing database data and methods for speeding 
    up searches/reducing search space

    ...

    Attributes
    ----------
    fasta_file : str
        the string of the path to the source protein database
    proteins : dict
        dictionary of proteins keyed by names and values of list of
        sequences (1 sequence for no digest, 1+ for digested)
    kmer_to_prots : dict
        dictionary that maps kmers found to source proteins
    tagged_b_ions : dict
        dictionary mapping masses to sequences
    tagged_y_ions : dict
        dictionary mapping masses to sequences

    Methods
    -------
    get_proteins_with_subsequence(subsequence)
        Gets all proteins with a certain subsequence

    get_protein(prot_name)
        Get all sequences associated with a protein

    reduce_search_space(masses)
        Index the database against a set of masses for faster lookups later

    get_mass_tags(mass, tolerance)
        Get all amino acid tags assicated with a mass. NOTE: reduce_search_space
        must be run first in order to run this
    '''

    fasta_file = ''
    proteins = {}
    kmer_to_prots = defaultdict(list)
    tagged_b_ions = defaultdict(list)
    tagged_y_ions = defaultdict(list)
    _mz_to_boundaries = {}

    def __init__(self, fasta_file: str): 
        '''
        Inputs:
            fasta_file: (str)   fasta database file
        '''
        self.fasta_file = fasta_file

        # build the databse immediately
        self._build()

    #----------------- Private functions --------------------#
    def _build(self):
        # pull the name out
        get_name = lambda x: x.split('|')[-1].split()[0]

        for entry in fasta.read(self.fasta_file):
            p_name = get_name(entry.description)
            self.proteins[p_name] = [entry.sequence]

    def _hashable_boundaries(self, boundaries: list) -> str:
        return '-'.join([str(x) for x in boundaries])

    def _merge(self, P: Iterable, indices: Iterable, kmers: Iterable, boundaries: Iterable):
        b_i, p_i = 0, 0

        matched_masses = defaultdict(list)

        while b_i < len(boundaries) and p_i < len(P):

            # if P[p_i] is in the boundary, keep track of it increment p_i
            if boundaries[b_i][0] <= P[p_i] <= boundaries[b_i][1]:
                matched_masses[self._hashable_boundaries(boundaries[b_i])] += kmers[indices[p_i - 1]:indices[p_i]]
                p_i += 1

            # if the upper boundary is less than the p_i, increment b_i
            elif P[p_i] > boundaries[b_i][1]:
                b_i += 1

            # if the lower boundary is greater than the p_i, increment p_i
            elif P[p_i] < boundaries[b_i][0]:
                p_i += 1

        return matched_masses

    def _make_database_set(self, proteins: list, max_len: int) -> (list, list, list):
        '''
        Create parallel lists of (masses, index_maps, kmers) for the merge sort operation

        Inputs:
            proteins:   (list) protein entries of shape (name, entry) where entry has a .sequence attribute
            max_len:    (int) the max length kmer to generate
        Ouputs:
            db_list_b:      (array) all b ion masses created in iteration through the proteins
            index_list_b:   (array) the index mapping. Equal size to db_list_b where each entry maps 
                                    idx(db_list_b) -> range of indices in kmer_list_b of kmers corresponding to masses
            kmer_list_b:    (list) the list of kmers associated to db_list_b through index_list_b

            db_list_y:      (array) all y ion masses created in iteration through the proteins
            index_list_y:   (array) he index mapping. Equal size to db_list_y where each entry maps 
                                    idx(db_list_y) -> range of indices in kmer_list_y of kmers corresponding to masses
            kmer_list_y:    (list) the list of kmers associated to db_list_b through index_list_y

            kmer_set:       (dict) mappings of kmers -> source proteins

        '''
        db_dict_b = defaultdict(set)
        db_dict_y = defaultdict(set)
        kmer_set = defaultdict(list)

        # function to add all masses of b+, b++, y+, y++ ions to the correct dict for sorting
        # and keep track of the source proteins
        def add_all(kmer, prot_name):
            for ion in 'by':
                for charge in [1, 2]:
                    spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)['spectrum']

                    for i, mz in enumerate(spec):
                        kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:]
                        r_d = db_dict_b if ion == 'b' else db_dict_y
                        r_d[mz].add(kmer_to_add)
                        kmer_set[kmer_to_add].append(prot_name)

        plen = len(proteins)

        # go through each protein and add all kmers to the correct dictionary for later sorting
        for i, (prot_name, prot_seq) in enumerate(proteins):
            
            print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')

            for j in range(1, max_len):
                kmer = prot_seq[:j]
                add_all(kmer, prot_name)
                
            for j in range(len(prot_seq) - max_len):
                kmer = prot_seq[j:j+max_len]
                add_all(kmer, prot_name)

            for j in range(len(prot_seq) - max_len, len(prot_seq)):
                kmer = prot_seq[j:]
                add_all(kmer, prot_name)

        print('\nSorting the set of protein masses...')
        
        # arrays take less space than lists
        db_list_b, index_list_b, kmer_list_b = arr.array('f'), arr.array('i'), []
        db_list_y, index_list_y, kmer_list_y = arr.array('f'), arr.array('i'), []

        # sort the keys to make sure that we are going through masses the correct way
        sorted_keys = sorted(db_dict_b.keys())
        for mz in sorted_keys:

            # get all kmers associated with this mass, append it the b ion list, and keep track 
            # of the kmers in the kmer list
            kmers = db_dict_b[mz]
            db_list_b.append(mz)
            offset = 0 if not len(index_list_b) else index_list_b[-1]
            index_list_b.append(len(kmers) + offset)
            kmer_list_b += kmers

        sorted_keys = sorted(db_dict_y.keys())
        for mz in sorted_keys:

            # get all kmers associated with this mass, append it the y ion list, and keep track 
            # of the kmers in the kmer list
            kmers = db_dict_y[mz]
            db_list_y.append(mz)
            offset = 0 if not len(index_list_y) else index_list_y[-1]
            index_list_y.append(len(kmers) + offset)
            kmer_list_y += kmers

        return db_list_b, index_list_b, kmer_list_b, db_list_y, index_list_y, kmer_list_y, kmer_set


    def _find_mz_to_boundaries(self, masses: list, boundaries: list):
        '''
        Find the mapping of what masses map to which boundaries 

        Inputs:
            masses:     (list) original input masses
            boundaires: (list) list of [lower_bound, upper_bound]
        '''

        # keep track of where at in each list we are
        m_i = 0
        b_i = 0

        while m_i < len(masses) and b_i < len(boundaries):

            # see if m_i is in between the current boundaries
            if boundaries[b_i][0] <= masses[m_i] <= boundaries[b_i][1]:
                self._mz_to_boundaries[masses[m_i]] = self._hashable_boundaries(boundaries[b_i])
                m_i += 1
            
            # otherwise if the lower bound is greater than mass, increase mass
            elif masses[m_i] <= boundaries[b_i][0]:
                m_i += 1

            # otherwise upper bound is less than m so increase b_i
            else:
                b_i += 1

    #---------------- Public functions --------------------#
    def get_proteins_with_subsequence(self, subsequence: str) -> list:
        '''
        Get a list of protien names that contain the subsequence

        Inputs:
            subsequence: (str) the subsequence to look for
        Ouptus:
            (list) names of proteins that contain the subseqeunce
        '''
        prots = []

        # first check the kmers dictionary
        if subsequence in self.kmer_to_prots:
            return self.kmer_to_prots[subsequence]

        # if its not in the kmers dict, look through every prot
        for prot_name, sequences in self.proteins.items():
            found = False
            for seq in sequences:
                if subsequence in seq:
                    found = True
                    prots.append(prot_name)
                    break 

            if found: 
                break 
        
        return prots

    def get_protein(self, prot_name: str) -> list:
        '''
        Get all sequences associated with a protein

        Inputs:
            prot_name:  (str) the name of the protein to get
        Outputs:
            (list) all sequences associated with the protein
        '''
        if prot_name in self.proteins:
            return self.proteins[prot_name]

        return []


    def reduce_search_space(self, masses: list, tolerance: int):
        '''
        Index the database against a set of masses for faster lookups 
        by using get_mass_tags

        Inputs:
            masses:     (list)  floats of all masses in an input run
            tolerance:  (int)   the tolerance (in ppm) to allow when searching
        '''
        print('Indexing database against spectra for faster searches...')

        # create the boundaries and remove overlaps
        def make_boundaries(mz):
            da_tol = ppm_to_da(mz, tolerance)
            return [mz - da_tol, mz + da_tol]

        boundaries = [make_boundaries(mz) for mz in masses]

        # make overlapped boundaries larger boundaries
        boundaries = overlap_intervals(boundaries)

        # create the mappings from mz to boundaries
        self._find_mz_to_boundaries(masses, boundaries)

        # estimate the max len
        max_len = predicted_len(boundaries[-1][1])

        # calc the number of batches needed
        num_batches = ceil(len(self.proteins) / BATCH_SIZE)

        # create batches of proteins in the form of (prot name, prot entry)
        kv_prots = [(k, v) for k, v in self.proteins.items()]
        batched_prots = [kv_prots[i*BATCH_SIZE:(i+1)*BATCH_SIZE] for i in range(num_batches)]

        # go through each batch set and create the list representation, merge, and keep good prots
        for batch_num, batch_set in enumerate(batched_prots):

            print(f'On batch {batch_num + 1}/{num_batches}\n', end='')

            # expand batch set ot be prot_name, seq for each of the sequnces in our v
            expanded_batch_set = [(k, seq) for k, v in batch_set for seq in v]

            # create our list representation
            batch_b_list, index_list_b, batch_kmer_b, batch_y_list, index_list_y, batch_kmer_y, batch_kmer_set = self._make_database_set(expanded_batch_set, max_len)

            # find tha batch matched masses for both b and y ions
            matched_masses_b_batch = self._merge(batch_b_list, index_list_b, batch_kmer_b, boundaries)
            matched_masses_y_batch = self._merge(batch_y_list, index_list_y, batch_kmer_y, boundaries)

            # add these these hits to our function scoped set of masses and to the kmer set
            for k, v in matched_masses_b_batch.items():
                self.tagged_b_ions[k] += v
                
                # add all kmers to kmer set
                for kmer in v:
                    self.kmer_to_prots[kmer] += batch_kmer_set[kmer]

            for k, v in matched_masses_y_batch.items():
                self.tagged_y_ions[k] += v 

                # add all kmers to kmer set
                for kmer in v:
                    self.kmer_to_prots[kmer] += batch_kmer_set[kmer]

        print('Done')

    def get_mass_tags(self, mass: float) -> (list, list):
        '''
        Get all amino acid sequences associated with a mass

        Inputs:
            mass:   (float) the mass to search for
        Outputs:
            (list, list) b-ion matches, y-ion matches
        '''
        # first get the boundary map
        boundary_map = self._mz_to_boundaries[mass]

        b_matches = []
        y_matches = [] 

        # check if its in each set of ion mappings and add to our return
        if boundary_map in self.tagged_b_ions:
            b_matches = self.tagged_b_ions[boundary_map]

        if boundary_map in self.tagged_y_ions:
            y_matches = self.tagged_y_ions[boundary_map]

        return (b_matches, y_matches)
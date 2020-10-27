from src.utils import hashable_boundaries, predicted_len
from src.objects import Database
from src.cppModules import gen_spectra

from collections import defaultdict
from typing import Iterable
from math import ceil

import array as arr

BATCH_SIZE = 300

def merge(P: Iterable, indices: Iterable, kmers: Iterable, boundaries: Iterable):
    b_i, p_i = 0, 0

    matched_masses = defaultdict(list)

    while b_i < len(boundaries) and p_i < len(P):

        # if P[p_i] is in the boundary, keep track of it increment p_i
        if boundaries[b_i][0] <= P[p_i] <= boundaries[b_i][1]:
            matched_masses[hashable_boundaries(boundaries[b_i])] += kmers[indices[p_i - 1]:indices[p_i]]
            p_i += 1

        # if the upper boundary is less than the p_i, increment b_i
        elif P[p_i] > boundaries[b_i][1]:
            b_i += 1

        # if the lower boundary is greater than the p_i, increment p_i
        elif P[p_i] < boundaries[b_i][0]:
            p_i += 1

    return matched_masses

def make_database_set(proteins: list, max_len: int) -> (list, list, list):
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
                spec = gen_spectra.gen_spectrum(kmer, ion=ion, charge=charge)

                for i, mz in enumerate(spec):
                    kmer_to_add = kmer[:i+1] if ion == 'b' else kmer[-i-1:]
                    r_d = db_dict_b if ion == 'b' else db_dict_y
                    r_d[mz].add(kmer_to_add)
                    kmer_set[kmer_to_add].append(prot_name)

    plen = len(proteins)

    # go through each protein and add all kmers to the correct dictionary for later sorting
    for i, (prot_name, prot_entry) in enumerate(proteins):
        
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')

        for j in range(1, max_len):
            kmer = prot_entry.sequence[:j]
            add_all(kmer, prot_name)
            
        for j in range(len(prot_entry.sequence) - max_len):
            kmer = prot_entry.sequence[j:j+max_len]
            add_all(kmer, prot_name)

        for j in range(len(prot_entry.sequence) - max_len, len(prot_entry.sequence)):
            kmer = prot_entry.sequence[j:]
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


def match_masses(spectra_boundaries: list, db: Database, max_pep_len=30) -> (dict, dict, Database):
    '''
    Take in a list of boundaries from observed spectra and return a b and y
    dictionary that maps boundaries -> kmers

    Inputs:
        spectra_boundaries: (list) lists with [lower_bound, upper_bound]
        db:                 (Database) entries to look for kmers in
    kwargs:
        max_pep_len:        (int) maximum peptide length. Default=30
    Outputs:
        (dict, dict) mappings from boundaries -> kmers for (b, y) ions
    '''

    '''
    We need to create mappings from boundaries -> kmers for both b and y ions
    Emprically we found that ~300 proteins worth of kmers and floats takes ~4GB RAM
    So in order to keep in down we need to batch it. Heres what we'll do
    
      1. Define a constant above for the number of proteins to have per batch
      2. For the number of proteins in a batch
        1. Pass those proteins to make_database_set (returns a list of b, y kmers with masses and indices)
        2. Pass those lists ot merge
        3. The results of that are mappings of boundaries -> kmers, add those to some constant mapping we have
        4. Delete those old lists to clear memory
        5. Repeat
    
    Adding each batch's results to our mapping will create a large mapping. We will return
    that mapping and the updated database

    '''
    # keep track of all of the good mass matches and kmers
    matched_masses_b, matched_masses_y, kmer_set = defaultdict(list), defaultdict(list), defaultdict(list)

    # estimate the max len
    max_len = min(predicted_len(spectra_boundaries[-1][1]), max_pep_len)

    # calc the number of batches needed
    num_batches = ceil(len(db.proteins) / BATCH_SIZE)

    # create batches of proteins in the form of (prot name, prot entry)
    kv_prots = [(k, v) for k, v in db.proteins.items()]
    batched_prots = [kv_prots[i*BATCH_SIZE:(i+1)*BATCH_SIZE] for i in range(num_batches)]

    # go through each batch set and create the list representation, merge, and keep good prots
    for batch_num, batch_set in enumerate(batched_prots):

        print(f'On batch {batch_num + 1}/{num_batches}\n', end='')

        extended_batch_set = [(k, entry) for (k, v) in batch_set for entry in v]

        # create our list representation
        batch_b_list, index_list_b, batch_kmer_b, batch_y_list, index_list_y, batch_kmer_y, batch_kmer_set = make_database_set(extended_batch_set, max_len)

        # find tha batch matched masses for both b and y ions
        matched_masses_b_batch = merge(batch_b_list, index_list_b, batch_kmer_b, spectra_boundaries)
        matched_masses_y_batch = merge(batch_y_list, index_list_y, batch_kmer_y, spectra_boundaries)

        # add these these hits to our function scoped set of masses and to the kmer set
        for k, v in matched_masses_b_batch.items():
            matched_masses_b[k] += v 
            
            # add all kmers to kmer set
            for kmer in v:
                kmer_set[kmer] += batch_kmer_set[kmer]

        for k, v in matched_masses_y_batch.items():
            matched_masses_y[k] += v 

            # add all kmers to kmer set
            for kmer in v:
                kmer_set[kmer] += batch_kmer_set[kmer]

    #update kmers in db to the kmer set
    db = db._replace(kmers=kmer_set)

    return (matched_masses_b, matched_masses_y, db)
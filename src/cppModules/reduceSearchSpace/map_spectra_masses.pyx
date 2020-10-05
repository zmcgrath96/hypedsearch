# distutils: language = c++

from libcpp.string cimport string 
from libcpp.vector cimport vector
from cpython cimport array

from mapSpectraMasses cimport mappings, protein, boundary, mapBoundaries 

def map_boundaries(boundaries: list, proteins: list, max_kmer_length: int):
    '''
    Perform a merge search of boundaries on the list of input proteins and keep 
    the mappings from boundaries -> kmers for b and y ions as well as a mapping
    from kmer -> source proteins

    Inputs: 
        boundaries: (list) list of iterables of size two with (lowerbound, upperbound)
        proteins:   (list) list of objects with attributes .name and .sequence
        max_kmer_length:    (int) the longest kmer to allow

    Outputs:
        (dict, dict, dict) dictonaries of 
        (range->bionmatches, range->yionmatches, kmer->source proteins)
    '''
    # turn boundaries into a vector of boundary objects
    cdef vector[boundary] bs 
    for lBuB in boundaries:
        bs.push_back(boundary(lBuB[0], lBuB[1]))

    # turn proteins in to protein classes
    cdef vector[protein] ps 
    for entry in proteins:
        ps.push_back(protein(str.encode(entry.name), str.encode(entry.sequence)))
 
    # call map boundaries 
    cdef mappings * m = mapBoundaries(bs, ps, max_kmer_length)

    # return them all in dict form
    return (dict(m.matchedBMasses), dict(m.matchedYMasses), dict(m.kmerToProts))
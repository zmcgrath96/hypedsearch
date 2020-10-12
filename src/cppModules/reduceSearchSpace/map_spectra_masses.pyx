# distutils: language = c++

from libcpp.string cimport string 
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from cython.operator import dereference, postincrement

from mapSpectraMasses cimport mappings, protein, boundary, mapBoundaries 

boundary_delta = .002

def hashable_boundaries(boundary: list) -> str:
    return '-'.join([str(x) for x in boundary])

def map_cpp_boundaries_to_inputs(cpp_boundaries: list, input_boundaries: list) -> dict:
    '''
    Since c++ treats floats slightyly different than python, we need to convert
    the boundary hashes from the cpp floats back to the input floats from python

    We will do a linear scan of both and see if both the first and last 
    are less than  .002 different and that is our hash mapping
    '''
    c_i, b_i = 0, 0

    cpp_to_py = {}

    hash_to_boundaries = lambda x: [float(x.split('-')[0]), float(x.split('-')[1])]

    while c_i < len(cpp_boundaries) and b_i < len(input_boundaries):
        cpp_b = hash_to_boundaries(cpp_boundaries[c_i])
        inp_b = input_boundaries[b_i]

        # check to see if they are close enough
        if abs(cpp_b[0] - inp_b[0]) < boundary_delta and abs(cpp_b[1] - inp_b[1]) < boundary_delta:
            cpp_to_py[cpp_boundaries[c_i]] = hashable_boundaries(input_boundaries[b_i])

            # increment just c_i
            c_i += 1

        # if cpp is lower than input, increment c_pp
        elif (cpp_b[0] - inp_b[0] + cpp_b[1] - inp_b[1]) < 0:
            c_i += 1

        # otherwise increment b_i
        else: 
            b_i += 1

    return cpp_to_py

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
    print('Converting python objects to c++...')
    # turn boundaries into a vector of boundary objects
    cdef vector[boundary] bs 
    cdef boundary b
    for lBuB in boundaries:
        b = boundary(lBuB[0], lBuB[1])
        bs.push_back(b)

    # turn proteins in to protein classes
    cdef vector[protein] ps 
    cdef protein p
    for entry in proteins:
        p = protein(str.encode(entry.name), str.encode(entry.sequence))
        ps.push_back(p)
    print('Done')
    # call map boundaries 
    cdef mappings * m = mapBoundaries(bs, ps, max_kmer_length)

    print('\nConverting from c++ objects to python...')
    # convert from unordered_maps to dictionaries
    matched_b_masses, matched_y_masses, kmer_to_prots = {}, {}, {}

    # first do b masses
    cdef unordered_map[string, vector[string]].iterator bIt = m.matchedBMasses.begin()
    # keep track of the keys in order to map back to python
    cpp_b_masses_keys = []

    while(bIt != m.matchedBMasses.end()):
        # get the key
        key = dereference(bIt).first.decode()
        cpp_b_masses_keys.append(key)

        # get the value
        value = [x.decode() for x in dereference(bIt).second]

        # put it in my return dictionary
        matched_b_masses[key] = value

        # increment to the next element
        postincrement(bIt)
 
    # next y masses
    cdef unordered_map[string, vector[string]].iterator yIt = m.matchedYMasses.begin()
    # keep track of the keys in order to map back to python
    cpp_y_masses_keys = []

    while(yIt != m.matchedYMasses.end()):
        # get the key
        key = dereference(yIt).first.decode()
        cpp_y_masses_keys.append(key)

        # get the value
        value = [x.decode() for x in dereference(yIt).second]

        # put it in my return dictionary
        matched_y_masses[key] = value

        # increment to the next element
        postincrement(yIt)

    # finally kmer to protein
    cdef unordered_map[string, unordered_set[string]].iterator kIt = m.kmerToProts.begin()
    
    while(kIt != m.kmerToProts.end()):
        # get the key
        key = dereference(kIt).first.decode()

        # get the value
        value = [x.decode() for x in dereference(kIt).second]

        # put it in my return dictionary
        kmer_to_prots[key] = list(set(value))

        # increment to the next element
        postincrement(kIt)

    # get our key mapping
    b_cpp_to_py = map_cpp_boundaries_to_inputs(
        sorted(cpp_b_masses_keys, key=lambda x: float(x.split('-')[0])), 
        boundaries
    )
    y_cpp_to_py = map_cpp_boundaries_to_inputs(
        sorted(cpp_y_masses_keys, key=lambda x: float(x.split('-')[0])), 
        boundaries
    )

    # remap to a NEW dictionary
    py_matched_b_masses, py_matched_y_masses = {}, {}

    for k, v in matched_b_masses.items():
        if k not in b_cpp_to_py:
            continue
        py_key = b_cpp_to_py[k]
        py_matched_b_masses[py_key] = list(set(v))

    del matched_b_masses

    for k, v in matched_y_masses.items():
        if k not in y_cpp_to_py:
            continue
        py_key = y_cpp_to_py[k]
        py_matched_y_masses[py_key] = list(set(v))

    print('Done')

    return (py_matched_b_masses, py_matched_y_masses, kmer_to_prots)
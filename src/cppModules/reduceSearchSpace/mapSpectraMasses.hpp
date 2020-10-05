#ifndef MAPSPECTRAMASSES_H
#define MAPSPECTRAMASSES_H

#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include<array>

#include "massSearch.hpp"


/**
 * Outward facing function for the caller. Does the merge search and mappings to map boundaries to kmers
 * 
 * @param boundaries    std::vector<boundary>   vector of float arrays of length 2 of the form [lowerbound, upperbound]
 * @param proteins      std::vector<protein>    vector of protein objects needed for indexing
 * @param maxKmerLength int                     the maximum length kmer to allow for
 * 
 * @return mappings
*/
mappings * mapBoundaries(std::vector<boundary> boundaries, std::vector<protein> proteins, int maxKmerLength);


#endif
#ifndef MAPSPECTRAMASSES_H
#define MAPSPECTRAMASSES_H

#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>

#include "massSearch.hpp"

/**
 * Mappings that we need to pass back to caller in order to go from ranges ('lowerbound-upperbound') -> kmers
 * for both b and y ions as well as a mapping of the important kmers -> source proteins
*/
class mappings {
public:
    std::unordered_map<std::string, std::vector<std::string>> matchedBMasses;
    std::unordered_map<std::string, std::vector<std::string>> matchedYMasses;
    std::unordered_map<std::string, std::vector<std::string>> kmerToProts;

    mappings();
    ~mappings();
};

/**
 * Class for called to turn proteins into so its easier for us to handle
*/
class protein {
public:
    std::string name;
    std::string sequence;
};

/**
 * Internal mappings that we need to generate that we do NOT pass back to caller. This is just or us 
 * to be able to "index" the database then call our merge operations
*/
class _internalMappings {
public:
    std::vector<float> bIonRefMasses;
    std::vector<int> bMassesToKmerIndices; 
    std::vector<std::string> bKmers; 
    std::vector<float> yIonRefMasses;
    std::vector<int> yMassesToKmerIndices; 
    std::vector<std::string> yKmers; 
    std::unordered_map<std::string, std::vector<std::string>> kmerToProts;

    _internalMappings();
    ~_internalMappings();
};

/**
 * Outward facing function for the caller. Does the merge search and mappings to map boundaries to kmers
 * 
 * @param boundaries    std::vector<float []>   vector of float arrays of length 2 of the form [lowerbound, upperbound]
 * @param proteins      std::vector<protein>    vector of protein objects needed for indexing
 * @param maxKmerLength int                     the maximum length kmer to allow for
 * 
 * @return mappings
*/
mappings * mapBoundaries(std::vector<float []> boundaries, std::vector<protein> proteins, int maxKmerLength);


#endif
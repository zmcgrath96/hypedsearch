#ifndef MASSSEARCH_H
#define MASSSEARCH_H

#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <iostream>

#include "mapSpectraMasses.hpp"

/**
 * Merge search a vector of boundaries against a vector of masses and keep track 
 * of the masses that were found in the returned map. The map is of the form
 * 'lowerbound-upperbound' and the value is a vector of resulting strings
 * 
 * @param masses        std::vector<float>      the vector of database masses 
 * @param indices       std::vector<float>      a parallel vector to masses that maps index(mass) -> set of kmers
 * @param kmers         std::vector<std::string>a vector with ranges for kmers that match a mass
 * @param boundaries    std::vector<float []>   a vector where every entry has an array of length 2 for [lowerbound, upperbound]
 * 
 * @return  std::unordered_map<std::string, std::vector<std::string>> maps the 'lowerbound-upperbound' to matched kmers
*/
std::unordered_map<std::string, std::vector<std::string>> mergeSearch(
    std::vector<float> masses, 
    std::vector<int> indices, 
    std::vector<std::string> kmers, 
    std::vector<float []> boundaries
);

/**
 * Create all the internal mappings needed for mergeSearch 
 * 
 * @param boundaries    std::vector<float []>   vector of float arrays of the form [lowerbound, upperbound]
 * @param proteins      std::vector<protein>    vector of protein objects to index
 * @param maxKmerLength int                     max kmer length to consider
 * 
 * @return _internalMappings *
*/
_internalMappings * indexProteins(std::vector<float []> boundaries, std::vector<protein> proteins, int maxKmerLength);

/**
 * Turn a boundary set into a string to be able to hash
 * 
 * @param boundary  float []    an array of length 2 of [lowerbound, upperbound]
 * 
 * @return std::string of the form 'lowerbound-upperbound'
*/
std::string hashifyBoundary(float boundary []);


#endif
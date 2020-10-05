#ifndef MASSSEARCH_H
#define MASSSEARCH_H

#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <iostream>

/**
 * Mappings that we need to pass back to caller in order to go from ranges ('lowerbound-upperbound') -> kmers
 * for both b and y ions as well as a mapping of the important kmers -> source proteins
*/
class mappings {
public:
    std::unordered_map<std::string, std::vector<std::string>> matchedBMasses;
    std::unordered_map<std::string, std::vector<std::string>> matchedYMasses;
    std::unordered_map<std::string, std::vector<std::string>> kmerToProts;

    mappings(){}
    ~mappings(){}
};

/**
 * Class for called to turn proteins into so its easier for us to handle
*/
class protein {
public:
    std::string name;
    std::string sequence;

    protein(std::string name, std::string sequence);
    protein(){}
    ~protein(){}
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

    _internalMappings(){}
    ~_internalMappings(){}
};

/**
 * Class for us to be able to handle boundaries. Trying to wrap c++ std::array
 * is finicky for cython, so its much easier to just do this
*/
class boundary{
public:
    float lowerBound;
    float upperBound;

    boundary(float lowerBound, float upperBound);
    boundary(){};
    ~boundary(){};
};

/**
 * Merge search a vector of boundaries against a vector of masses and keep track 
 * of the masses that were found in the returned map. The map is of the form
 * 'lowerbound-upperbound' and the value is a vector of resulting strings
 * 
 * @param masses        std::vector<float>      the vector of database masses 
 * @param indices       std::vector<float>      a parallel vector to masses that maps index(mass) -> set of kmers
 * @param kmers         std::vector<std::string>a vector with ranges for kmers that match a mass
 * @param boundaries    std::vector<boundary>   a vector where every entry has an array of length 2 for [lowerbound, upperbound]
 * 
 * @return  std::unordered_map<std::string, std::vector<std::string>> maps the 'lowerbound-upperbound' to matched kmers
*/
std::unordered_map<std::string, std::vector<std::string>> mergeSearch(
    std::vector<float> masses, 
    std::vector<int> indices, 
    std::vector<std::string> kmers, 
    std::vector<boundary> boundaries
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
_internalMappings * indexProteins(std::vector<boundary> boundaries, std::vector<protein> proteins, int maxKmerLength);

/**
 * Turn a boundary set into a string to be able to hash
 * 
 * @param boundary  boundary    an array of length 2 of [lowerbound, upperbound]
 * 
 * @return std::string of the form 'lowerbound-upperbound'
*/
std::string hashifyBoundary(boundary boundary);


#endif
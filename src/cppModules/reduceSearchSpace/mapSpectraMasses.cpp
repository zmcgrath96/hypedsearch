#include "mapSpectraMasses.hpp"

/**
 * Outward facing function for the caller. Does the merge search and mappings to map boundaries to kmers
 * 
 * @param boundaries    std::vector<float []>   vector of float arrays of length 2 of the form [lowerbound, upperbound]
 * @param proteins      std::vector<protein>    vector of protein objects needed for indexing
 * @param maxKmerLength int                     the maximum length kmer to allow for
 * 
 * @return mappings
*/
mappings * mapBoundaries(std::vector<float []> boundaries, std::vector<protein> proteins, int maxKmerLength){
    // create the internal mappings needed
    _internalMappings * iMaps = indexProteins(boundaries, proteins, maxKmerLength);

    // create the return type
    mappings * returnMappings = new mappings();

    // do the b mappings first
    returnMappings->matchedBMasses = mergeSearch(iMaps->bIonRefMasses, iMaps->bMassesToKmerIndices, iMaps->bKmers, boundaries);

    // delete the stuff from the iMap that we don't need
    iMaps->bIonRefMasses.clear();
    iMaps->bMassesToKmerIndices.clear();

    // do the y mappings second
    returnMappings->matchedBMasses = mergeSearch(iMaps->bIonRefMasses, iMaps->bMassesToKmerIndices, iMaps->bKmers, boundaries);

    // delete teh y ioin stuff we don't need
    iMaps->yIonRefMasses.clear();
    iMaps->yMassesToKmerIndices.clear();

    // keep the kmers that we need from our kmerToProt mapping
    for (auto const entry: returnMappings->matchedBMasses){
        // the second entry is our list of kmers
        std::vector<std::string> kmers = entry.second;

        // got through these, find the source prot, and put it in the returnMappings kmer mapping
        for (std::string kmer: kmers){
            std::vector<std::string> sourceProts = iMaps->kmerToProts[kmer];

            // put these in return val
            returnMappings->kmerToProts[kmer].insert(returnMappings->kmerToProts[kmer].end(), sourceProts.begin(), sourceProts.end());
        }
    }

    // keep the kmers that we need from our kmerToProt mapping
    for (auto const entry: returnMappings->matchedYMasses){
        // the second entry is our list of kmers
        std::vector<std::string> kmers = entry.second;

        // got through these, find the source prot, and put it in the returnMappings kmer mapping
        for (std::string kmer: kmers){
            std::vector<std::string> sourceProts = iMaps->kmerToProts[kmer];

            // put these in return val
            returnMappings->kmerToProts[kmer].insert(returnMappings->kmerToProts[kmer].end(), sourceProts.begin(), sourceProts.end());
        }
    }

    delete iMaps;
    return returnMappings;
}


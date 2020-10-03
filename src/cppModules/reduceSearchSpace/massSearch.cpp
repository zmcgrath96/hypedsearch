#include "massSearch.hpp"
#include "../spectraGeneration/genSpectra.hpp"

/**
 * Turn a boundary set into a string to be able to hash
 * 
 * @param boundary  float []    an array of length 2 of [lowerbound, upperbound]
 * 
 * @return std::string of the form 'lowerbound-upperbound'
*/
std::string hashifyBoundary(float boundary[]){
    return std::to_string(boundary[0]) + "-" + std::to_string(boundary[1]);
}

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
){
    // index into masses and boundaries respectively
    int m_i, b_i = 0;

    // keep track of masses that we have matched
    std::unordered_map<std::string, std::vector<std::string>> matchedMasses;

    while (m_i < masses.size() && b_i < boundaries.size()){

        // if masses at m_i is inside our current boundaries[b_i], keep track of the 
        // kmers at this position and increment m_i
        if (boundaries[b_i][0] <= masses[m_i] && masses[m_i] <= boundaries[b_i][1]){
            // get the string representation of the boundary for hashing
            std::string boundaryHash = hashifyBoundary(boundaries[b_i]);

            // get the kmers that need to be added from the kmer list
            std::vector<std::string> kmersToAdd;

            // start index for getting kmers
            int startIndex = m_i == 0? 0: indices[m_i - 1];

            for (int i = startIndex; i < indices[m_i]; i ++){
                kmersToAdd.push_back(kmers[i]);
            }

            // append our set of kmers to the set of kmers in the map
            // vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
            matchedMasses[boundaryHash].insert(matchedMasses[boundaryHash].end(), kmersToAdd.begin(), kmersToAdd.end());
        
            // finally increment m_i
            m_i ++;
        }

        // if the upper bound is less than masses[m_i], increment b_i
        else if (masses[m_i] > boundaries[b_i][1]) b_i ++;

        // the lower bound is greater than p_i, incrment p_i
        else m_i ++;
    }

    // finally return the map
    return matchedMasses;
}

/**
 * Create all the internal mappings needed for mergeSearch 
 * 
 * @param boundaries    std::vector<float []>   vector of float arrays of the form [lowerbound, upperbound]
 * @param proteins      std::vector<protein>    vector of protein objects to index
 * @param maxKmerLength int                     max kmer length to consider
 * 
 * @return _internalMappings *
*/
_internalMappings * indexProteins(std::vector<float []> boundaries, std::vector<protein> proteins, int maxKmerLength){
    // init  our return value
    _internalMappings * returnVal = new _internalMappings();

    // keep track of the kmers being assigned to masses
    std::map<float, std::vector<std::string>> bIonMassKmers;
    std::map<float, std::vector<std::string>> yIonMassKmers;

    auto addKmerLambda = [&bIonMassKmers, &yIonMassKmers, &returnVal](std::string kmer, std::string proteinName){
        for (std::string ion: {"b", "y"}){
            for (int charge: {1, 2}){
                // create a spectrum
                std::vector<float> spectrum = genSpectrum(kmer, ion, charge, true);

                // go through the masses and add all masses with their respective kmer
                // in. EX: kmer=ABCD, ion='b'. Spectrum would be {100, 250, 340, 500}
                // so then we would iteratively add pairs of 
                // A, 100   AB, 250    ABC, 340    ABCD, 500
                // to the bionmasskmers
                for (int i = 0; i < spectrum.size(); i ++){

                    // get the correct substring from the kmer
                    std::string kmerToAdd = ion == "b"? kmer.substr(0, i+1) : kmer.substr(kmer.size()-i-1, i+1);

                    // get the mass from the spectrum
                    float mz = spectrum[i];

                    // if its a b ion, we want to add to the bIonMassKmers, otherwise the yIonMassKmers
                    if (ion == "b") bIonMassKmers[mz].push_back(kmerToAdd);
                    else yIonMassKmers[mz].push_back(kmerToAdd);
                    
                    // finally add the kmer to the kmer set 
                    returnVal->kmerToProts[kmerToAdd].push_back(proteinName);
                }
                
            }
        }
    }; // end of lambda

    int numProts = proteins.size();

    // go through each protein and add all kmers to dictionaries
    for (int protNum = 0; protNum < numProts; protNum ++){
        std::cout << "\rOn protein " << protNum + 1 << "/" << numProts << " [" << (int)(protNum + 1) * 100 / numProts << "%]";

        // get the current prot to avoid calling it more than once
        protein prot = proteins[protNum];

        // do the first bit 
        int edgeIterRange = maxKmerLength > prot.sequence.length() ?
                                prot.sequence.length() : maxKmerLength;

        // do the edges of the protein first
        for (int j = 0; j < edgeIterRange; j ++){
            // get the kmer and call addKmerLambda on it
            std::string frontKmer = prot.sequence.substr(0, j + 1);
            std::string backKmer = prot.sequence.substr(prot.sequence.length() - j - 1, j+1);
            addKmerLambda(frontKmer, prot.name);
            addKmerLambda(backKmer, prot.name);
        }

        // if the length of the protein is less than the max kmer length, continue
        if (maxKmerLength > prot.sequence.length()) continue;

        // otherwise do all the kmers in the middle
        for (int j = 1; j < prot.sequence.length() - maxKmerLength; j ++){
            std::string kmer = prot.sequence.substr(j, maxKmerLength);
            addKmerLambda(kmer, prot.name);
        }
    
    }

    // all the keys in our b and y ion maps are sorted (http://www.cplusplus.com/reference/map/map/)
    // so we can just iterate through the values normally to create our parallel lists

    for (auto const entry: bIonMassKmers){
        float mass = entry.first;
        std::vector<std::string> kmers = entry.second;

        // append the mass to the interal returnVal b ref masses
        returnVal->bIonRefMasses.push_back(mass);

        // get the offset needed for the index list
        int offset = returnVal->bMassesToKmerIndices.size() == 0 ? 0 : returnVal->bMassesToKmerIndices.back();

        // append the next range of kmers to the massesToKmerIndices
        returnVal->bMassesToKmerIndices.push_back(kmers.size() + offset);

        // append these kmers the the internal kmers of b returnval
        returnVal->bKmers.insert(returnVal->bKmers.end(), kmers.begin(), kmers.end());
    }

    // same thing for the y ion
    for (auto const entry: yIonMassKmers){
        float mass = entry.first;
        std::vector<std::string> kmers = entry.second;

        // append the mass to the interal returnVal y ref masses
        returnVal->yIonRefMasses.push_back(mass);

        // get the offset needed for the index list
        int offset = returnVal->yMassesToKmerIndices.size() == 0 ? 0 : returnVal->yMassesToKmerIndices.back();

        // append the next range of kmers to the massesToKmerIndices
        returnVal->yMassesToKmerIndices.push_back(kmers.size() + offset);

        // append these kmers the the internal kmers of b returnval
        returnVal->yKmers.insert(returnVal->yKmers.end(), kmers.begin(), kmers.end());
    }

}



#include "loadFasta.cpp"
#include "genSpectra.hpp"
#include "set"

int main(){
    std::string fastaFile = "/Users/zacharymcgrath/Desktop/nod2 data/all data/NOD2_mouse_database.fasta";
    std::cout << "Loading fasta...\n";
    std::vector<std::string> prots = readFasta(fastaFile);
    std::cout << "Done.\n";
    std::vector<float> databaseMassSet;

    for (int i = 0; i < prots.size(); i++){
        std::string prot = prots[i];
        std::cout << "\rOn protein " << i + 1 << "/" << prots.size();

        if (prot.size() < 30) continue;

        for (int j = 0; j < 29; j ++){
            std::string kmer = prot.substr(0, j);
            std::vector<float> thisSpec = genSpectrum(kmer, "", -1, false);
            for (float mass: thisSpec) databaseMassSet.push_back(mass);
        }

        // break into kmers of size 30
        for (int j = 0; j < prot.length() - 30; j ++){
            std::string kmer = prot.substr(j, 30);
            std::vector<float> thisSpec = genSpectrum(kmer, "", -1, false);
            for (float mass: thisSpec) databaseMassSet.push_back(mass);
        }

        for (int j = prot.length() - 29; j < prot.length() - 1; j ++){
            std::string kmer = prot.substr(j, prot.length() - j);
            std::vector<float> thisSpec = genSpectrum(kmer, "", -1, false);
            for (float mass: thisSpec) databaseMassSet.push_back(mass);
        }
    
    }

    std::cout << "\nNumber of elements in the list: " << databaseMassSet.size() << "\n";

    return 1;
}
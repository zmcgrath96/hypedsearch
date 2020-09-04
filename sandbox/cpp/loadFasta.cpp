#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <regex>

std::vector<std::string> readFasta(std::string fastaFile){
    std::string line;
    std::ifstream database (fastaFile);
    std::string prot;
    std::vector<std::string> prots;

    std::regex newlines_re("\n+");

    if (database.is_open()){
        while (getline (database, line))
        {
            if (line.substr(0, 1) == ">"){
                // add prot to the list
                prots.push_back(prot);

                // clear prot
                prot = "";
            }
            else {
                // remove newline and append to prot
                std::string result = std::regex_replace(line, newlines_re, "");
                prot += result;
            }
        }
        prots.push_back(prot);
        database.close();
    }
    return prots;
}
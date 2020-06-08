# hypedsearch
## **hy**brid **pe**pti**d**e **search**ing tool

hypedsearch is a tool for identifying both hybrid and non-hybrid protein sequences from mass spectrometry data. hypedsearch takes MS/MS files (`mzML` only currently) and `fasta` files as inputs and identifies peptides, both hybrid and nonhybrid.  

## Usage

**hypedsearch** uses python3, thus to run the program you must have python3 installed. If you have not installed it, instructions to do so are [here](https://www.python.org/). 
### Intsallation
First clone the repository
```bash
$> git clone https://github.com/zmcgrath96/hypedsearch.git
```
### Install dependencies
```bash
$> cd hypedsearch
$/hypedsearch> pip3 install -r requirements.txt
```
### Run tests
### Testing
To run all unit tests run the following:
```bash
$> cd hypedsearch
$hypedsearch> python3 -m unittest -v
```

### Examples

## References
### pyteomics
Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717
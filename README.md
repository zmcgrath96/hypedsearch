# hypedsearch
## **hy**brid **pe**pti**d**e **search**ing tool

`hypedsearch` is a tool for identifying both hybrid and non-hybrid protein sequences from mass spectrometry data. `hypedsearch` takes MS/MS files (`mzML` only currently) and `fasta` files as inputs and identifies peptides, both hybrid and nonhybrid.  

`hypedsearch` identifies sequences by using a method called "k-mer extension". If you're not familiar, a "k-mer" is a k long string, or in this case, k long sequence of amino acids. The process, at a high level, works like this:
1. Pre-processing of the protein database to identify all k-mers for k in the range `min_peptide_len` to `max_peptide_len`. These k-mers are generated from both the `N` and `C` termini
2. Using the k-mers generated from the `N` termninus side, attempt to identify a sequence of amino acids that describe the `b` ions in the observed spectrum.
3. Repeat step 2, but from the `C` terminus side and try to describe the `y` ions in the observed spectrum
4. Filter out the poor scoring sequences
5. For the rest of the sequences, attempt to align and overlap the two sequences to the spectrum
6. If two sequences have no overlap, or do overlap but are from different proteins, report the alignment as hybrid
7. Save all alignments

Currently `hypedsearch` works on Unix systems and has yet to be tested on Windows.

## Setup

**hypedsearch** uses python3. If you don't have python3, you can find instructions to download [here](https://www.python.org/). 
### Intsallation
First clone the repository
```bash
$> git clone https://github.com/zmcgrath96/hypedsearch.git
```
### Install dependencies
```bash
$> cd hypedsearch
$hypedsearch> pip3 install -r requirements.txt
```
### Run tests
### Testing
To run all unit tests run the following:
```bash
$> cd hypedsearch
$hypedsearch> python3 -m unittest -v
```

## Usage
`hypedsearch` is fairly strightforward. The only necessary inputs are 
1. A directory where all of your spectra are saved in `mzML` file formats
2. A `fasta` file protein database

Other optional parameters are available, however these are the only two needed to get started. If you run with the default options, output is saved in the `~/` directory. If you would like to see all of the available parameters, simply run `python3 -m src.main --help` from the `hypedsearch` directory. You will then see the following help menu

```
hypedsearch> python3 -m src.main --help
usage: main.py [-h] [--output-dir SAVE_DIR]
               [--min-peptide-len MIN_PEPTIDE_LEN]
               [--max-peptide-len MAX_PEPTIDE_LEN] [--tolerance TOLERANCE]
               [--verbose VERBOSE] [--score SCORING_ALG]
               SF DB

Tool for identifying proteins, both hybrid and non hybrid from MS/MS data

positional arguments:
  SF                    Path to folder containing spectra files.
  DB                    Path to .fasta file containing proteins

optional arguments:
  -h, --help            show this help message and exit
  --output-dir SAVE_DIR
                        Directory to save all figures. Default=~/
  --min-peptide-len MIN_PEPTIDE_LEN
                        Minimum peptide length to consider. Default=5
  --max-peptide-len MAX_PEPTIDE_LEN
                        Maximum peptide length to consider. Default=20
  --tolerance TOLERANCE
                        ppm tolerance to allow in search. Deafult=20
  --score SCORING_ALG   Scoring algorithm to use. Options are [bb, ion, ibb]
                        for backbone, ion, and ion backbone respectively.
                        Default=bb
  --verbose VERBOSE     Extra printing to console during run. Default=True

```


### Parameters
* __output-dir__: The directory to save all output data in. Default is the home directory `~/`.
* __min-peptide-len__: This is the minimum peptide length to consider. This is also the minimum length to consider for half of a hybrid peptide. For example, if you want to search for hybrid peptides that may only have 3 amino acids on a side, this parameter should be set to 3. 
* __max-peptide-len__: The maximum peptide length to consider. This value is used for pre-processing on the protein database. Longer peptide lengths, especially on larger databases, increases both the time and memory needed to run. 
* __tolerance__: The parts per million (ppm) tolerance to allow when trying to match an experimental mass to a theoretical mass. This value is converted to Da. The resulting range allowed is a theoretical mass +/- the mass. For example, a 20 ppm tolerance would be a 40 ppm range (both + and - 20 ppm).
* __score__: The scoring algorithm to use when creating an alignment. Three options exists:
    1. Backbone score (bb): Score is a function not of number of ions matched, but what they describe. For example, 4 ions that describe 4 bond sites scores better than 4 matched ions that describe 2 bonds. 
    2. Ion backbone score (ibb)(Default): Same as the backbone score, but only matches to `b` or `y` ions when coming building a kmer from the `N` (left) or `C` (right) terminus.
    3. Ion matching (ion): Score is blind to where the ions are, but rather looks for the number of ions being matched. 
* __verbose__: Printing to the console for constant updates. Recommended to leave as True (the default)


## References
### pyteomics
*Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717
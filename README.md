# hypedsearch
## **hy**brid **pe**pti**d**e **search**ing tool

For more information about hypedsearch, installation and usage instructions, 
please see the documentation [found here](https://hypedsearch.readthedocs.io/en/latest/). 

---

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
## References
### pyteomics
*Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717
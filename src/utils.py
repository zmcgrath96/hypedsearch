import os, gzip, shutil, copy, math
from typing import Iterable, Any
from itertools import product
import numpy as np
from collections import namedtuple

from src.objects import Spectrum
from src import gen_spectra

import math

def make_valid_dir_string(dir_path: str) -> str:
    '''
    Add / character to end of directory string to make valid directory path

    Inputs:
        dir_path:   string name of directory to check
    Outputs:
        string      dir_path with / as final character
    '''
    return dir_path + '/' if dir_path[-1] != '/' else dir_path

def make_dir(dir_path: str) -> bool:
    '''
    Check directory path for existing directory or make one with the name given

    Inputs:
        dir_path:   string full path to directory to create
    Ouputs:
        bool        True if successful False otherwise
    '''
    dir_path = make_valid_dir_string(dir_path)
    try:
        if not os.path.exists(dir_path): 
            os.makedirs(dir_path)
        return True 
    except:
        return False

def make_valid_text_file(file_name: str) -> str:
    '''
    Ensure some string path has .txt appended to it for appropriate .txt extension

    Inputs:
        file_name:  string file name to make valid txt file name
    Outputs:
        string      file name guaranteed to have .txt at the end of it
    '''
    file_name = file_name + '.txt' if '.txt' not in file_name else file_name
    return file_name

def make_valid_json_file(file_name: str) -> str:
    '''
    Ensure some string path has .json appended to it for appropriate .json extension

    Inputs:
        file_name:  string file name to make valid json file name
    Outputs:
        string      file name guaranteed to have .json at the end of it
    '''
    file_name = file_name + '.json' if '.json' not in file_name else file_name
    return file_name

def make_valid_csv_file(file_name: str) -> str:
    '''
    Ensure some string path has .csv appended to it for appropriate .csv extension

    Inputs:
        file_name:  string file name to make valid csv file name
    Outputs:
        string      file name guaranteed to have .csv at the end of it
    '''
    file_name = file_name + '.csv' if '.csv' not in file_name else file_name
    return file_name

def make_valid_fasta_file(file_name: str) -> str:
    '''
    Ensure some string path has .fasta appended to it for appropriate .fasta extension

    Inputs:
        file_name:  string file name to make valid fasta file name
    Outputs:
        string      file name guaranteed to have .fasta at the end of it
    '''
    file_name = file_name + '.fasta' if '.fasta' not in file_name else file_name
    return file_name

def file_exists(file_name: str) -> bool:
    '''
    Determine if a file exists

    Inputs:
        file_name:  string path to the file to check for
    Outputs:
        bool        True if file exists False otherwise
    '''
    return os.path.isfile(file_name)

def gzip_file(file_name: str, delete_old=True) -> str:
    '''
    Compress a file with gzip compression

    Inputs:
        file_name:  string path to the file name of the file to compress
    kwargs:
        delete_old: bool delete the unzipped file. Default=True
    Outputs:
        string      file path of the compressed file
    '''
    compressed_file_name = file_name + '.gz'
    with open(file_name, 'rb') as f_in:
        with gzip.open(compressed_file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    delete_old and os.remove(file_name)
    return compressed_file_name

def gunzip_file(compressed_file_name: str, delete_old=True) -> str:
    '''
    Decompress a gzipped file

    Inputs:
        compressed_file_name:   string path to the gzipped file
    kwargs: 
        delete_old:             bool delete the compressed file. Default=True
    Ouputs:
        string                  path to the unzipped file
    '''
    file_name = compressed_file_name if '.gz' not in compressed_file_name else compressed_file_name.replace('.gz', '')
    with gzip.open(compressed_file_name, 'rb') as f_in:
        with open(file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    delete_old and os.remove(compressed_file_name)
    return file_name

def is_gzipped(file_name: str) -> bool:
    '''
    Determines if a file has been compressed with gzip

    Inputs:
        file_name:  string path to file to test
    Outputs:
        bool        True if the folder is gzipped false otherwise
    '''
    return '.gz' == file_name[-3:]

def gzip_dir(d: str, delete_old=True) -> str:
    '''
    Compress a directory with gzip compression

    Inputs:
        d:          string path to the directory
    kwargs:
        delete_old: bool delete the original directory. Default=True
    Ouptuts:
        str         path to new zipped directory
    '''
    root = '/'.join(d.split('/')[:-1])
    shutil.make_archive(d, 'zip', root)
    delete_old and shutil.rmtree(d)
    return d + '.zip'

def is_json(file: str) -> bool:
    '''
    Determine if a file is a json file

    Inputs:
        file:   string file name to test
    Ouputs:
        bool    True if is json false otherwise
    '''
    return True if '.json' in file else False

def is_fasta(file: str) -> bool:
    '''
    Determine if a file is a fasta file 

    Inputs:
        file:   string file to determine if it is fasta
    Outputs:
        bool    True if is fasta False otherwise
    '''
    return True if '.fasta' in file else False

def is_dir(dir_path: str) -> bool:
    '''
    Determine if a path is a valid path to a directory

    Inputs:
        dir_path:   string full path to directory
    Outputs:
        bool        True if directory exists False otherwise
    '''
    return os.path.isdir(dir_path)

def is_file(file: str) -> bool:
    '''
    Determine if a file exists

    Inputs:
        file:   string full path to the file
    Outputs:   
        bool    True if file exists False otherwise
    '''
    return os.path.isfile(file)

def insort_by_key(value: Any, a: Iterable, key: str) -> Iterable:
    '''
    Insert value into list in order. List a given should be sorted prior to insort
    
    Inputs:
        value:    (Any) thing to be inserted
        a:        (Iterable) the list or list-like value to insert into
        key:      (str) the key to use for comparison of two values
    Ouputs:
        Iterable a with value inserted in order by key
    '''
    if len(a) == 0:
        return [value]
    elif len(a) == 1:
        return a + [value] if a[0][key] <= value[key] else [value] + a
    
    mid = math.floor(len(a)/2)
    if a[mid-1][key] <= value[key] <= a[mid][key]:
        return a[:mid] + [value] + a[mid:]
    elif a[mid][key] > value[key]:
        return insort_by_key(value, a[:mid], key) + a[mid:]
    else:
        return a[:mid] + insort_by_key(value, a[mid:], key)

def insort_by_index(value: Any, a: Iterable, index: int) -> Iterable:
    '''
    Insert value into list in order. List a given should be sorted prior to insort
    
    Inputs:
        value:    (Any) thing to be inserted
        a:        (Iterable) the list or list-like value to insert into
        index:    (int) the index to use for comparison of two values
    Ouputs:
        Iterable a with value inserted in order by index
    '''
    if len(a) == 0:
        return [value]
    elif len(a) == 1:
        return a + [value] if a[0][index] <= value[index] else [value] + a
    
    mid = math.floor(len(a)/2)
    if a[mid-1][index] <= value[index] <= a[mid][index]:
        return a[:mid] + [value] + a[mid:]
    elif a[mid][index] > value[index]:
        return insort_by_index(value, a[:mid], index) + a[mid:]
    else:
        return a[:mid] + insort_by_index(value, a[mid:], index)

def insort_by_func(value: Any, a: Iterable, func: callable) -> Iterable:
    '''
    Insort value into a using func as the the callable (on 2 inputs to be compared)

    Example:
        value:  12.34
        a:      [1.90, 5.90, 12.33, 45]

        When func is called, it is called on 2 inputs (func(12.34, 12.33)) and returns
        True if the first input is smaller than the second, otherwise False

    Inputs:
        value:  (Any) value to be added to the list
        a:      (Iterable) the iterable to add Any into
        func:   (Callable) function with inputs (value, value2) and outputs True if value2 > value
    Outputs:
        (Iterable) updated value of a
    '''

    if len(a) == 0:
        return [value]
    elif len(a) == 1:
        return a + [value] if func(a[0], value) else [value] + a
    
    mid = math.floor(len(a)/2)
    if func(a[mid-1], value) and func(value, a[mid]):
        # check to see if the nieghboring values are the same, if not don't add the new one
        if a[mid-1] == value or a[mid] == value:
            return a
        return a[:mid] + [value] + a[mid:]
    elif func(value, a[mid]):
        return insort_by_func(value, a[:mid], func) + a[mid:]
    else:
        return a[:mid] + insort_by_func(value, a[mid:], func)

def all_perms_of_s(s: str, keyletters: str) -> list:
    '''
    Find all permutations of a string that has values 'keyletters' in them
    
    Inputs:
        s:          (str) the string to evaluate
        keyletters: (str) the letters to permutate
    Outputs:
        list of all the permutations of keyletters
    '''

    # Convert input string into a list so we can easily substitute letters
    seq = list(s)
    
    perms = []

    # Find indices of key letters in seq
    indices = [ i for i, c in enumerate(seq) if c in keyletters ]

    # Generate key letter combinations & place them into the list
    for t in product(keyletters, repeat=len(indices)):
        for i, c in zip(indices, t):
            seq[i] = c
        perms.append(''.join(seq))
    return perms

def ppm_to_da(mass: float, ppm_tolerance: float) -> float:
    '''
    Calculate the absolute boundary value from an mass mass. Value returned is in Da
    
    Inputs:
        mass:           (float) the mass mass used to calculate the tolerances
        ppm_tolerance:  (float or int) the tolerance in ppm to use to convert to Da
    Outputs:
        float value in Da 
    '''
    return abs((ppm_tolerance / 1000000)*mass)

def make_sparse_array(spectrum: list, width: float, value=50) -> np.ndarray:
    '''
    Make a spectrum (a list of floats) into a sparsely populated array for xcorr 
    calculation. Indices are calculated by
    
    idx = int(m/w), m is mass, w is bin width
    
    width is the tolerance in Da to allow when calculating scores. All peaks with some value
    are given a new value of 50.
    
    Inputs:
        spectrum:   (list) float mass values of peaks
        width:      (float) mass tolerance to accept to make bin width
    kwargs:
        value:      (number) value to give peaks at the new index. Default=50
    Outputs:
        (np.ndarray) sparesly populated list
    '''
    # find the largest mass and make that the length of the array
    list_size = int(max(spectrum)//width) + 1
    
    sparse = np.zeros(list_size)
    
    # populate sparse at the index for each mass
    for m in spectrum:
        sparse[int(m // width)] = value

    return sparse

def overlap_intervals(intervals: list) -> list:
    '''
    Take a list of intervals and turn it into a smaller list by finding any 
    overlapping intervals and making it a larger interval
    Inputs:
        intervals:  (list) intervals (in the form of lists) of upper and lower bounds (inclusive)
    Outputs:
        (list) intervals of lists
    '''
    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged  

def to_percent(index, total):
    return int(100 * (index + 1)/total)

def predicted_len(precursor_mass: float, precursor_charge: int) -> int:
    '''
    The predicted length of a spectrum based on its maximum mass

    Inputs:
        precursor_mass:   float   the maximum mass of the sequence
        precursor_charge: int     the charge of the precursor
    Outputs:
        (int) predicted length
    '''
    return math.ceil(precursor_mass / gen_spectra.get_precursor('G', precursor_charge))

def predicted_len_precursor(spectrum: Spectrum, sequence: str) -> int:
    '''
    Make a prediction of the peptide length give a spectrum and the current 
    sequence. 

    Inputs:
        spectrum: (Spectrum)
        sequence: (str)
    Outputs:
        int
    '''
    # first get the theoretical precursor of the sequence
    theoretical_prec = gen_spectra.get_precursor(sequence, spectrum.precursor_charge)

    # now we can estimate length   spec/seq = real/theory
    estimated_len = math.ceil(len(sequence) * (spectrum.precursor_mass / theoretical_prec))

    return estimated_len

def hashable_boundaries(boundaries: list) -> str:
    '''
    Turn a lower and upper bound into a string in order 
    to hash

    Inputs:
        boundaries:     (list) [lower_bound: float, upper_bound: float]
    Outputs:
        (str) <lower_bound>-<upper_bound>
    '''
    return '-'.join([str(x) for x in boundaries])

def cosine_similarity(a: list, b: list) -> float:
    '''
    Calculate the cosine similarity of two vectors
    
    Inputs:
        a:   (list)
        b:   (list)
    Outputs:
        (float) cosine similarity (a dot b)/(norm(a) * norm(b))
    '''
    # pad with zeros if not the same size
    if len(a) > len(b):
        b = list(b) + [0 for _ in range(len(a) - len(b))]
        
    if len(b) > len(a):
        a = list(a) + [0 for _ in range(len(b) - len(a))]

    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

def DEV_contains_truth_parts(truth_seq: str, hybrid: bool, b_seqs: list, y_seqs: list) -> bool:
    '''
    DEV FUNCTION ONLY

    Determines if a set of b and y sequences can potentially create the truth. If so, 
    True is returned, otherwise False

    Inputs:
        truth_seq:  (str)   the truth sequence trying to find
        hybrid:     (bool)  True if the truth is a hybrid
        b_seqs:     (list)  k-mers identified from the b ion score
        y_seqs:     (list)  k-mers identified from the y ion score
    Outputs:
        True if the truth could be constructructed from current set of seqs, False otherwise
    '''
    has_left = False
    has_right = False

    left_half = ''
    right_half = ''

    # if the sequence is hybrid, replace all I and L with B
    if hybrid:
        
        # split into left and right halves for later 
        if '-' in truth_seq:
            left_half = truth_seq.split('-')[0]
            right_half = truth_seq.split('-')[1]

        elif '(' in truth_seq and ')' in truth_seq:
            left_half = truth_seq.split(')')[0].replace('(', '')
            right_half = truth_seq.split('(')[1].replace(')', '')

        else: 
            left_half = truth_seq[:2]
            right_half = truth_seq[-2:]

        truth_seq = truth_seq.replace('I', 'B').replace('L', 'B').replace('-', '').replace('(', '').replace(')', '')

        b_seqs = [x.replace('I', 'B').replace('L', 'B') for x in b_seqs]
        y_seqs = [x.replace('I', 'B').replace('L', 'B') for x in y_seqs]

    # see if any of the b seqs match the first bit of the sequence
    has_left = any(
        [x == truth_seq[:len(x)] for x in b_seqs if len(x) > 1]
    ) or any(
        [truth_seq == x[:len(truth_seq)] for x in b_seqs if len(x) > 1]
    )

    has_right = any(
        [x == truth_seq[-len(x):] for x in y_seqs if len(x) > 1]
    ) or any(
        [truth_seq == x[-len(truth_seq):] for x in y_seqs if len(x) > 1]
    )

    # if its a hybrid, both left and right must be true, otherwise just one will do 
    if hybrid:
        
        # if we dont' have complete matches, thats ok. If we can get the halves 
        # of the "truth" in some longer b or y seq, thats good enough
        if not has_left:
            left_half = left_half.replace('I', 'B').replace('L', 'B')
            has_left = any([left_half == x[:len(left_half)] for x in b_seqs])

        if not has_right:
            right_half = right_half.replace('I', 'B').replace('L', 'B')
            has_right = any([right_half == x[-len(right_half):] for x in y_seqs])

        return has_left and has_right

    return has_left or has_right


def DEV_contains_truth_exact(truth_seq: str, hybrid: bool, seqs: list) -> bool:
    '''
    DEV FUNCTION ONLY

    Determines of the truth sequence is held in the list of sequences. If so, 
    True is returned, otherwise False

    Inputs:
        truth_seq:  (str)   the truth sequence trying to find
        hybrid:     (bool)  True if the truth is a hybrid
        seqs:       (list)  sequences to look through
    Outputs:
        True if the truth is in current set of seqs, False otherwise
    '''

    # if hybrid, replace all L and I with B, and also remove any - or ()
    if hybrid:
        truth_seq = truth_seq \
                        .replace('I', 'B') \
                        .replace('L', 'B') \
                        .replace('-', '') \
                        .replace('(', '') \
                        .replace(')', '')

        seqs = [
            x \
                .replace('I', 'B') \
                .replace('L', 'B') \
                .replace('-', '') \
                .replace('(', '') \
                .replace(')', '')
            for x in seqs
        ]

    contains_exact = any([x == truth_seq for x in seqs])

    return contains_exact
import os, gzip, shutil, copy, math
from typing import Iterable, Any
from itertools import product

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
        return a + [value] if getattr(a[0], key) <= getattr(value, key) else [value] + a
    
    mid = math.floor(len(a)/2)
    if getattr(a[mid-1], key) <= getattr(value, key) <= getattr(a[mid], key):
        return a[:mid] + [value] + a[mid:]
    elif getattr(a[mid], key) > getattr(value, key):
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
import os, gzip, shutil, copy

def get_related_files(files: list, sub: str, not_sub=None) -> list:
    '''
    Find the subset of files passed in with the substring without the exclusive substring

    Inputs:
        files:      list of file names (strings)
        sub:        string substring to search for in each file name
    kwargs:
        not_sub:    string substring to make sure is not in each file name
    Outputs:
        list        list of strings. Subset of files that have sub and don't have not_sub
    '''
    if not_sub is not None and not_sub != '':
        return [x for x in files if sub in x and not_sub not in x]
    return [x for x in files if sub in x]

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


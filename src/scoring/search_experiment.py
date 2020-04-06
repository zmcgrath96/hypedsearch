import pyopenms
import os
import shutil
from subprocess import call
from utils.utils import make_dir, make_valid_dir_string, is_gzipped, gunzip_file, gzip_file, gzip_dir
from scoring import search
from file_io import fasta, mzML

crux_to_rm = ['tide-search.decoy.txt', 'tide-search.log.txt', 'tide-search.params.txt']
scoring_functions = ['crux', 'custom']

#######################################################################################
#                       BEGIN "PRIVATE" FUNCTIONS
#######################################################################################

def __parse_spectrum_name(spec_name):
    return str(spec_name.split('/')[-1]).lower().replace('.mzml', '')

def __parse_db_name(db_name):
    return  str(db_name.split('/')[-1]).replace('.fasta', '')

def __index_db_files(path_to_crux_cmd: str, db_files: list) -> list:
    '''
    Do some pre indesing when running to crux search to avoid doing it while running search

    Inputs:
        path_to_crux_cmd:   string absolute path to the crux bin command
        db_files:           list of strings of database file paths to index
    Outputs:
        list of indexed file names (absolute path)
    '''
    idx_names = []
    num_dbs = len(db_files)

    for i, db_file in enumerate(db_files):
        print('On database {}/{}[{}%]\r'.format(i+1, num_dbs, int(((i+1)/num_dbs) * 100)), end="")
        
        this_output_dir = '/'.join(str(db_file).split('/')[:-1])
        this_output_dir = make_valid_dir_string(this_output_dir) + 'indexed/'
        make_dir(this_output_dir)
        idx_name = this_output_dir + str(db_file).replace('.fasta', '_index').split('/')[-1]

        indx_cmd = [
            path_to_crux_cmd, 
            'tide-index', 
            db_file, 
            idx_name, 
            '--min-length', '2', 
            '--min-mass', '50', 
            '--output-dir', this_output_dir, 
            '--overwrite', 'T', 
            '--min-peaks', '2', 
            '--precursor-window', '1000000000',
            '--enzyme', 'no-enzyme', 
            '--verbosity', '0'
        ]
        call(indx_cmd)
        # remove extra output files
        os.remove(this_output_dir + 'tide-index.params.txt')
        os.remove(this_output_dir + 'tide-index.log.txt')

        idx_names.append(idx_name)
    return idx_names

def __remove_indices(index_file):
    if isinstance(index_file, list):
        index_file = index_file[0]
    rm_dir = '/'.join(index_file.split('/')[:-1])
    shutil.rmtree(rm_dir)

def __crux_search(spectra_files: list, database_files: list, path_to_crux_cmd: str, output_dir: str, compress=True):
    '''
    Generate scores from the crux scoring algorihm

    Inputs:
        spectra_files: list of str paths to all the spectra (.mzML) files
        database_files: list of str paths to all the database (.fasta) files
        path_to_crux_cmd: str path to the executable for crux.
        output_dir: str path to the directory to save files
    kwargs:
        compress: bool compress the output result. Default=True
    Outputs:
        list of str of output files
    '''
    output_dir = make_valid_dir_string(output_dir) + 'search_output/'
    make_dir(output_dir)
    spec_dir = '/'.join(spectra_files[0].split('/')[:-1])

    is_compressed = is_gzipped(spectra_files[0])
    print('Pre-indexing database files...')
    indexed_db_files = __index_db_files(path_to_crux_cmd, database_files)
    print('\nDone. Scoring..')

    output_count = 0
    output_files = []
    num_dbs = len(indexed_db_files)
    num_specs = len(spectra_files)

    for i, spec_file in enumerate(spectra_files):
        spec_file = spec_file if not is_compressed else gunzip_file(spec_file)
        for j, database_file in enumerate(indexed_db_files):
            this_db_name = __parse_db_name(database_file)
            print('On spectrum {}/{}[{}%]\tOn database file {}/{}[{}%]\r'.format(i+1, num_specs, int(((i+1)/num_specs) * 100), j+1, num_dbs, int(((j+1)/num_dbs)*100)), end="")
            this_output_dir = output_dir + '{}_vs_{}'.format(__parse_spectrum_name(spec_file), this_db_name)
            search_cmd = [
                path_to_crux_cmd, 
                'tide-search', 
                spec_file, 
                database_file, 
                '--min-length', '2', 
                '--min-mass', '50', 
                '--output-dir', this_output_dir, 
                '--overwrite', 'T', 
                '--min-peaks', '2', 
                '--precursor-window', '1000000000',
                '--enzyme', 'no-enzyme', 
                '--verbosity', '0'
                ] 
            call(search_cmd)
            output_count += 1
            o = this_output_dir + '/tide-search.target.txt' if not compress else gzip_file(this_output_dir + '/tide-search.target.txt')
            o_tsv = o.replace('.txt', '.tsv')
            os.rename(o, o_tsv)
            output_files.append(o_tsv)

            # saving space, so should remove the extra stuff
            if is_compressed:
                for rm in crux_to_rm:
                    os.remove(this_output_dir + '/' + rm)

        # if the files were compressed, we were trying to save disk space so just remove the mzml files
        is_compressed and os.remove(spec_file)

    is_compressed and os.rmdir(spec_dir)
    __remove_indices(indexed_db_files)
    return output_files

def __load_databases(database_files: list) -> dict:
    '''
    Load the .fasta databases into memory

    Inputs:
        database_files: list of strings of full paths to .fasta files
    Outputs:
        dictionary of lists of the form
        {
            '<database file name>': [{protein info}]
        }
    '''
    dbs = {}
    for db_name in database_files:
        dbs[db_name] = fasta.read(db_name)
    return dbs

def __load_spectra(spectra_files: list) -> dict:
    '''
    Load the .mzML files into memory

    Inputs:
        spectra_files: list of strings of full paths to .mzML files
    Ouputs:
        dictionary of lists of the form
        {
            '<spectra file name>': [spectra]
        }
    '''
    spectra = {}
    for spectra_file in spectra_files:
        spectra_file = spectra_file if not is_gzipped(spectra_file) else gunzip_file(spectra_file)
        spectra[spectra_file] = mzML.read(spectra_file)
    return spectra

def __custom_search(spectra_files: list, database_files: list, output_dir: str, compress=True) -> list:
    '''
    Use a custom scoring function to score spectra

    Inputs:
        spectra_files: list of strings paths to all .mzML files
        database_files: list of strings paths to all .fasta files
        output_dir: str path to directory to save results under
    OPTONAL:
        compress: bool whether or not to compress results. Default=True
    Outputs:
        list of str of output files
    '''
    output_dir = make_valid_dir_string(output_dir)
    make_dir(output_dir)
    output_files = []

    no_spec = len(spectra_files)
    no_db = len(database_files)

    # preload stuff into memory to reduce file io
    print('Loading databases into memory...')
    databases = __load_databases(database_files)
    print('Done')
    print('Loading spectra into memory...')
    spectra = __load_spectra(spectra_files)
    print('Done')

    # running out of room so gzip those directories
    # db_dir = '/'.join(database_files[0].split('/')[:-1])
    # spect_dir = '/'.join(spectra_files[0].split('/')[:-1])
    # gzip_dir(db_dir)
    # gzip_dir(spect_dir)

    for spec_no, spectra_file in enumerate(spectra):
        for db_no, db in enumerate(databases):
            print('On spectrum: {}/{} [{}%]   On database: {}/{} [{}%]\r'.format(spec_no, no_spec, int(float(spec_no) / float(no_spec) * 100), db_no, no_db, int(float(db_no) / float(no_db) * 100)), end='')
            output_name = output_dir + 'search_output/' + '{}_vs_{}'.format(__parse_spectrum_name(spectra_file), __parse_db_name(db))
            output_file = search.search_database(spectra[spectra_file], databases[db], output_name)
            output_file = output_file if not compress else gzip_file(output_file)
            output_files.append(output_file)
    return output_files
#######################################################################################
#                        END "PRIVATE" FUNCTIONS
#######################################################################################

def score_peptides(spectra_files: list, database_files: list, output_dir: str, compress=True, score_func='custom', path_to_crux_cmd='') -> list:
    '''
    Use a scoring function or scoring tool to generate scores of spectra vs databases
  
    Inputs:
        spectra_files: list of str paths to all the spectra (.mzML) files
        database_files: list of str paths to all the database (.fasta) files
        output_dir: str path to the directory to save files
    kwargs:
        compress: bool compress the output result. Default=True
        score_func: str determine which scoring function to use. Default='custom'
        path_to_crux_cmd: str path to the executable for crux. Default=''
    Outputs:
        list of str of output files
    '''
    score_func = 'custom' if score_func.lower() not in scoring_functions else score_func.lower()
    if score_func == 'crux':
        return __crux_search(spectra_files, database_files, path_to_crux_cmd, output_dir, compress=compress)
    else:
        return __custom_search(spectra_files, database_files, output_dir, compress=compress)

from src.objects import MassSequence, Spectrum, Database
from src.utils import ppm_to_da, make_dir, is_file
from src.sequence.gen_spectra import max_mass

from collections import defaultdict
from pyteomics import fasta

import pandas as pd
import string
import math
import pickle

########################## Private Functions ##########################

def __make_db_file(db: Database) -> str:
    '''
    Make the name and directories to save the datbase into. The subdirectories
    are made, but no file is made. The structure is as follows:

        databases/<fasta_file_name>/f<min len>t<max len>.db

    If no fasta file is made, then "noname" is used as the intermediate folder name

    Inputs:
        db:     (Databse) the database for which the file name is made for
    Outputs:
        (str) the name of the file (full path) for which the database is saved
    '''
    # make the database directory
    db_file = './databases/'
    make_dir(db_file)

    # make the intermediate folder
    db_file += db.fasta_file.split('/')[-1].replace('.fasta', '') + '/' \
        if db.fasta_file is not None and db.fasta_file != '' else 'noname/'
    make_dir(db_file)

    # add the f t clause
    db_file += f'f{db.min_len}t{db.max_len}.db'

    return db_file


def __saved_db_exists(db: Database) -> bool:
    '''
    Check to see if a saved database exists

    Inputs:
        db:     (Database) the namedtuple with fasta file and min, max lens set to see if exists
    Outputs:
        (bool) True if the database exists, false otherwise
    '''
    db_file = __make_db_file(db)
    
    # if noname is the middle directory, return false
    if 'noname' in db_file:
        return False

    print(f'{db_file} file exists: {is_file(db_file)}')
    # return if the file exists
    return is_file(db_file)

def __read_fasta(db: Database, fasta_file: str) -> dict:
    '''
    Read proteins into memory from fasta file
    
    Inputs:
        db:         (Database) object to insert proteins into
        fasta_file: (str) path to fasta file
    Outputs:
        (Database) updated with the proteins
    '''
    prots = []
    db = db._replace(fasta_file=fasta_file)

    db.verbose and print('Loading fasta file into memory...')

    # split the name on the OS value if it exists
    get_name = lambda name: name[:name.index('OS=')-1] if 'OS=' in name else name

    # go through each entry in the fasta and put it in memory
    for i, entry in enumerate(fasta.read(fasta_file)):

        # take the description without the 'sp' value
        desc = entry.description.split('|')[1:] if '|' in entry.description else entry.description

        # if the id is in the description, take it
        if len(desc) > 1:
            id_ = desc[0]
            name = get_name(desc[1])
            
        # make the id just the number
        else:
            id_ = i
            name = get_name(desc[0])
            
        # get the sequence
        seq = entry.sequence

        # make the entry and add it to prots
        prots.append({'name': name, 'id': id_, 'sequence': seq})

    # throw it into a pandas dataframe
    proteins = pd.DataFrame(prots)

    # update db's proteins attribute
    db = db._replace(proteins=proteins)

    db.verbose and print('Done')

    return db

########################## /Private Functions ##########################

########################## Public Functions ##########################

def save_db(db: Database) -> bool:
    '''
    Pickle and dump the database to use again later. 

    Inputs:
        db:     (Database) the database to save
    Outputs:
        (bool) True if the database could be saved successfuly, False otherwise
    '''
    db_file = __make_db_file(db)

    # try:
    print(f'Saving to: {db_file}')
    db = db._replace(tree=bytearray(db.tree))
    pickle.dump(db, open(db_file, 'wb'), )
    # except:
    #     return False

    return True


def build_or_load_db(db: Database) -> Database:
    '''
    Load the database if it exists, otherwise build it

    Inputs:
        db:     (Database) the database tuple to store the new database into
    Outputs:
        (Database) the updated database
    '''
    if False:#saved_db_exists(db):
        print('Loading saved database...')
        db = pickle.load(open(__make_db_file(db), 'rb'))
        print('Done')

    else: 
        db = __read_fasta(db, db.fasta_file)
        db = build(db)

        # # save it 
        # if not save_db(db):
        #     print('WARNING: Could not save built database to file')

    return db


def get_proteins_with_subsequence(db: Database, subsequence: str) -> list:
    '''
    Find all proteins that have the subsequence provided

    Inputs:
        db:             (Database) the database with all of the proteins
        subsequence:    (str) the subsequence to look for
    Outputs:
        (list) names of the proteins that contain the subsequence
    '''
    return list(db.proteins[db.proteins['sequence'].str.contains(subsequence)]['name'])


def get_entry_by_name(db: Database, name: str) -> dict:
    '''
    Get a protein entry from the database

    Inputs: 
        db:         (Database) the database to search
        name:       (str) name of the protein to get
    Outputs:
        (dict) Entry of the protein. Empty Entry if not found
    '''
    res = db.proteins[db.proteins['name'] == name].to_dict('r')

    if len(res):
        return res[0]
    return None


def build(
    db: Database
) -> Database:
    '''
    Build the database. In order to build, the fasta file should have been read and the database should
    contain the proteins dictionary. Call "read_fasta" first. This function builds the internal 
    KmerMasses object and the internal Tree 
    
    Inputs: 
        db:         (Database) the database to build internal structures for
    Outputs:
        (Database) updated databse
    '''

    # breakdown a sequence s into all its possible subsequences from min to max len
    def breakdown(s: str) -> list:
        kmers = []

        # go through each starting position of the sequence
        for j in range(len(s) - db.min_len):

            # make a kmer sequence. Do the max (to generate the kmer spec once) then 
            # just iterate through it
            kmer_len = db.max_len if j + db.max_len <= len(s) else len(s) - j

            for k in range(db.min_len, kmer_len):
                kmers.append(s[j:j+k])
            
        return kmers

    # get the singly and doubly b and y masses for this subsequence
    def spectrify(s: str) -> dict:
        f = {}
        f['bs'] = max_mass(s, 'b', 1)
        f['bd'] = max_mass(s, 'b', 2)
        f['ys'] = max_mass(s, 'y', 1)
        f['yd'] = max_mass(s, 'y', 2)
        f['sequence'] = s
        return f

    # for each protein sequence:
    #   1. break it down to all the subsequences via breakdown function (list output)
    #   2. make it a dataframe (DataFrame output)
    #   3. make each column's list entry (that is a list) into its own row (DataFrame output)
    #   4. remove any duplicates
    a = db.proteins['sequence']    \
        .apply(breakdown)       \
        .to_frame()             \
        .explode('sequence')    \
        .drop_duplicates('sequence')

    # for each of the subsequences that we created in a, spectrify it 
    # which creates a bs, bd, ys, yd column for each subsequence
    b = pd.DataFrame(list(a['sequence'].apply(spectrify)))

    # update db kmer_masses to be the dataframe
    db = db._replace(kmer_masses=b)

    return db


def search(db: Database, observed: Spectrum, kmers: str, tolerance: float) -> list:
    '''
    Search through all masses and saved kmers to find masses that are within our tolerance
    
    Inputs:
        db:         (Database) the database containing the kmermasses dictionaries to search
        spectrum:   (Spectrum) what to sequence
        kmers:      (str) which ion type we are searching. Options: ['bs', 'bd', 'ys', 'yd']
        tolerance:  (float) the ppm tolerance to accept for each mass
    Outputs:
        list of MassSequence for all masses that were in the acceptable range of an observed mass
    '''
    
    bounds = []
    hits = []

    # get all of the bounds of a spectrum
    for mz in observed.spectrum:
        tol = ppm_to_da(mz, tolerance)
        bounds.append((mz - tol, mz + tol))

    # go through each bound pair and search the dataframe
    for b in bounds:
        hits += list(db.kmer_masses[db.kmer_masses[kmers].between(b[0], b[1])]['sequence'])
            
    return hits

########################## /Public Functions ##########################
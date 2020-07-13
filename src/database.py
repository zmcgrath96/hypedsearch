
from src.objects import MassSequence, Spectrum, Database
from src.utils import ppm_to_da, make_dir, is_file
from src.sequence.gen_spectra import max_mass

from collections import defaultdict
from pyteomics import fasta
from more_itertools import flatten

import src.sqlite_interface as sql
import pandas as pd
import string
import math
import pickle
import psutil

########################## Private Functions ##########################

def __make_db_dir(db: Database) -> str:
    '''
    Make the directories to save the datbase into. The subdirectories
    are made. The structure is as follows:

        databases/<fasta_file_name>/

    If no fasta file is made, then "noname" is used as the intermediate folder name

    Inputs:
        db:     (Databse) the database for which the file name is made for
    Outputs:
        (str) the name of the full path for which the database is saved
    '''
    # make the database directory
    db_dir = './databases/'
    make_dir(db_dir)

    # make the intermediate folder
    db_dir += db.fasta_file.split('/')[-1].replace('.fasta', '') + '/' \
        if db.fasta_file is not None and db.fasta_file != '' else 'noname/'
    make_dir(db_dir)

    return db_dir


def __make_db_file(db: Database) -> str:
    '''
    Make the name and directories to save the datbase into. The subdirectories
    are made, but no file is made. The structure is as follows:

        databases/<fasta_file_name>/f<min len>t<max len>.pkl

    If no fasta file is made, then "noname" is used as the intermediate folder name

    Inputs:
        db:     (Databse) the database for which the file name is made for
    Outputs:
        (str) the name of the file (full path) for which the database is saved
    '''
    # make the database directory
    db_dir = __make_db_dir(db)

    # add the f t clause
    db_file = db_dir + f'f{db.min_len}t{db.max_len}.pkl'

    return db_file


def __make_sql_db_file(db: Database) -> str:
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
    # get the db dir
    db_dir = __make_db_dir(db)

    # make the .db file 
    db_file = db_dir + f'f{db.min_len}t{db.max_len}.db'

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

    # return if the file exists
    return is_file(db_file)


def __save_db(db: Database) -> bool:
    '''
    Pickle and dump the database to use again later. 

    Inputs:
        db:     (Database) the database to save
    Outputs:
        (bool) True if the database could be saved successfuly, False otherwise
    '''
    db_file = __make_db_file(db)

    try:
        print(f'Saving to: {db_file}')
        
        # temporarily remove the connection as its not pickleable
        conn = db.conn
        db = db._replace(conn=None)
        pickle.dump(db, open(db_file, 'wb'))

        # replace db conn
        db = db._replace(conn=conn)
    except:
        return False

    return True


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

    # create the connection and table for proteins
    conn = sql.create_connection(__make_sql_db_file(db))
    db = db._replace(conn=conn)
    sql.create_table(
        db.conn,
        '''
        CREATE TABLE IF NOT EXISTS proteins (
            name TEXT, 
            id TEXT PRIMARY KEY, 
            sequence TEXT
        )
        '''
    )

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

    # add these sequences to the table
    sql.execute_indel_many(
        db.conn, 
        '''
        INSERT INTO proteins (name, id, sequence) VALUES (?, ?, ?)
        ''', 
        [(x['name'], x['id'], x['sequence']) for x in prots]
    )

    db.verbose and print('Done')

    return db


def __build(
    db: Database
) -> Database:
    '''
    Build the database. In order to build, the fasta file should have been read and the database should
    contain the proteins dictionary. Call "read_fasta" first. This function builds pandas dataframes 
    
    Inputs: 
        db:         (Database) the database to build internal structures for
    Outputs:
        (Database) updated databse
    '''

    # breakdown a sequence s into all its possible subsequences from min to max len
    def breakdown(s: str) -> list:
        kmers = []

        # go through each starting position of the sequence
        for j in range(len(s) - db.min_len + 1):

            # make a kmer sequence. Do the max (to generate the kmer spec once) then 
            # just iterate through it
            kmer_len = db.max_len + 1 if j + db.max_len <= len(s) else len(s) - j + 1

            for k in range(db.min_len, kmer_len):
                kmers.append(s[j:j+k])

        return kmers

    # get the singly and doubly b and y masses for this subsequence
    def spectrify(s: str) -> dict:
        return (
            max_mass(s, 'b', 1),
            max_mass(s, 'b', 2),
            max_mass(s, 'y', 1),
            max_mass(s, 'y', 2),
            s
        )
       
    # create the kmers database in sqlite
    sql.create_table(
        db.conn, 
        '''
        CREATE TABLE IF NOT EXISTS kmers (
            bs REAL, 
            bd REAL, 
            ys REAL, 
            yd REAL, 
            sequence TEXT
        )
        '''
    )

    db.verbose and print(f'Finding all subsequences in range {db.min_len} to {db.max_len}')

    # in batches of 1000 proteins, apply the breakdown and spectrify functions
    #   1. get the number of proteins in the database
    #   2. calculate the number of batches to run
    #   3. for each batch   
    #       a. get the sequences of proteins in batch
    #       b. apply breakdown to each sequence
    #       c. apply spectrify to all sequences
    #   4. add the results of spectrify to the table
    batch_size = 1000

    # get the protein sequences
    sequences = [x[0] for x in sql.execute_return(
        db.conn, 
        '''
        SELECT sequence FROM proteins
        '''
    ).fetchall()]
    num_prots = len(sequences)

    # calculate the number of batches
    num_batches = math.ceil(num_prots / batch_size)

    # go through each batch
    for batch in range(num_batches):
        print(f'Indexing sequences for faster searches... {int(batch * 100 / num_batches)}%\r', end='')
        
        # the tuples to add to kmers
        kmer_additions = []

        # the sequences in the batch
        batch_sequences = sequences[batch*batch_size:(batch+1)*batch_size]
        
        for sequence in batch_sequences:
            kmer_additions += [
                spectrify(kmer) for kmer in breakdown(sequence)
            ]
        
        # add it to the table
        sql.execute_indel_many(
            db.conn,
            '''
            INSERT INTO kmers (bs, bd, ys, yd, sequence) VALUES (?, ?, ?, ?, ?)
            ''',
            kmer_additions
        )

    db.verbose and print('Done.')

    return db

########################## /Private Functions ##########################

########################## Public Functions ##########################


def build_or_load_db(db: Database) -> Database:
    '''
    Load the database if it exists, otherwise build it

    Inputs:
        db:     (Database) the database tuple to store the new database into
    Outputs:
        (Database) the updated database
    '''
    if __saved_db_exists(db):
        print('Loading saved database...')
        db = pickle.load(open(__make_db_file(db), 'rb'))

        # establish the sqlite connection
        db = db._replace(conn=sql.create_connection(__make_sql_db_file(db)))
        print('Done')

    else: 
        db = __read_fasta(db, db.fasta_file)
        db = __build(db)

        # save it 
        if not __save_db(db):
            print('WARNING: Could not save built database to file')

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
    # we take x[0] because the results come back as tuples in the form ('name', )
    return [x[0] for x in 
        sql.execute_return(
            db.conn, 
            f'''
            SELECT name FROM proteins
            WHERE sequence LIKE "%{subsequence}%"
            '''
        ).fetchall()
    ]


def get_entry_by_name(db: Database, name: str) -> dict:
    '''
    Get a protein entry from the database

    Inputs: 
        db:         (Database) the database to search
        name:       (str) name of the protein to get
    Outputs:
        (dict) Entry of the protein. Empty Entry if not found
    '''
    a = list(sql.execute_return(
        db.conn, 
        f'''
        SELECT * FROM proteins
        WHERE name == "{name}"
        '''
    ).fetchall())[0]

    if len(a) == 0:
        return {}

    return({
        'name': a[0],
        'id': a[1],
        'sequence': a[2]
    })

def cache_spectra(db: Database, spectra: list, ppm_tolerance=0, da_tolerance=0) -> None:
    '''
    Caching the hits of as many spectra peaks as possible. One
    of the two tolerances should be set in order to perform a search. If non is set, then default
    is 20 ppm.
    
    Inputs:
        db:         (Database) the database containing the kmermasses dictionaries to search
        spectra:    (list) spectra (each of Spectrum type) to cache
    kwargs:
        ppm_tolerance:  (float) the ppm tolerance to accept for each mass. Default=0
        da_tolerance:   (int) the tolerance to accept for each mass. Default=0
    Outputs:
        None
    '''


def search(db: Database, observed: Spectrum, kmers: str, ppm_tolerance=0, da_tolerance=0) -> list:
    '''
    Search through all masses and saved kmers to find masses that are within our tolerance. One
    of the two tolerances should be set in order to perform a search. If non is set, then default
    is 20 ppm.
    
    Inputs:
        db:         (Database) the database containing the kmermasses dictionaries to search
        spectrum:   (Spectrum) what to sequence
        kmers:      (str) which ion type we are searching. Options: ['bs', 'bd', 'ys', 'yd']
    kwargs:
        ppm_tolerance:  (float) the ppm tolerance to accept for each mass. Default=0
        da_tolerance:   (int) the tolerance to accept for each mass. Default=0
    Outputs:
        list of MassSequence for all masses that were in the acceptable range of an observed mass
    '''

    # create a tuple of a lower and upper bound
    def get_bounds(mz):
        tol = ppm_to_da(mz, ppm_tolerance) if ppm_tolerance > 0 \
            else (da_tolerance if da_tolerance > 0 else ppm_to_da(mz, 20))

        return (mz - tol, mz + tol)

    # get all bounds for all mz values in the observed
    bounds = [get_bounds(mz) for mz in observed.spectrum]

    # create a query as large as bounds / 2
    between_query = [f'{kmers} BETWEEN {b1} and {b2}' for b1, b2 in bounds]

    # create the full select query to run
    full_query = f'SELECT sequence FROM kmers WHERE {" or ".join(between_query)}'

    # run the query. we take x[0] because entries come as tuples ('sequence',)
    return [x[0] for x in list(sql.execute_return(
        db.conn, 
        full_query
    ).fetchall())]
            

########################## /Public Functions ##########################
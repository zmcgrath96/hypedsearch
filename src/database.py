
from src.objects import DatabaseEntry, KmerMasses, MassSequence, Spectrum, Database
from src.utils import ppm_to_da, make_dir, is_file
from src.sequence.gen_spectra import gen_spectrum

from collections import defaultdict
from pyteomics import fasta
from suffix_tree import Tree

import string
import math
import pickle

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


def saved_db_exists(db: Database) -> bool:
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
        db = read_fasta(db, db.fasta_file)
        db = build(db)

        # # save it 
        # if not save_db(db):
        #     print('WARNING: Could not save built database to file')

    return db

def read_fasta(db: Database, fasta_file: str, is_uniprot=False) -> dict:
    '''
    Read proteins into memory from fasta file
    
    Inputs:
        db:         (Database) object to insert proteins into
        fasta_file: (str) path to fasta file
    kwargs:
        is_uniprot: (bool) adds attribute 'human_readable_name' to dictionary if True. Default=False
    Outputs:
        (Database) updated with the proteins
    '''
    prots = {}
    e: DatabaseEntry
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
        e = DatabaseEntry(name, seq, id_, name)
        prots[name] = e

    db = db._replace(proteins=prots)

    db.verbose and print('Done')

    return db


def get_entry_by_name(db: Database, name: str) -> DatabaseEntry:
    '''
    Get a protein entry from the database

    Inputs: 
        db:         (Database) the database to search
        name:       (str) name of the protein to get
    Outputs:
        (DatbaseEntry) Entry of the protein. Empty Entry if not found
    '''
    return db.proteins.get(name, DatabaseEntry())


def add_entry(
    db: Database, 
    protein_name: str, 
    protein_sequnece: str, 
    protein_id='', 
    human_readable_name=''
) -> Database:
    '''
    Add an entry to the database. If it cannot be added to the database, False is returned

    Inputs:
        db:                 (Database) the database to add entries to
        protein_name:       (str) name of the protein being added
        protein_sequence:   (str) string of amino acids describing the protein
    kwargs:
        protein_id:         (str) identifier of the protein. Default=''
        human_readable_name:(str) name that is easy for people to read. Default=''
    Outputs:
        (Database) the updated database
    '''
    if protein_name in db.proteins:
        return db

    e = DatabaseEntry(protein_name, protein_sequnece, protein_id, human_readable_name)
    db.proteins[protein_name] = e 

    # add it to the tree
    if not db.tree:
        db = db._replace(tree=Tree())
        db.tree = Tree()

    db.tree.add(protein_name, protein_sequnece)

    return db


def build(
    db: Database
) -> KmerMasses:
    '''
    Build the database. In order to build, the fasta file should have been read and the database should
    contain the proteins dictionary. Call "read_fasta" first. This function builds the internal 
    KmerMasses object and the internal Tree 
    
    Inputs: 
        db:         (Database) the database to build internal structures for
    Outputs:
        (Database) updated databse
    '''
    db = db._replace(tree=Tree())

    # defaultdicts to hash integer value of masses into 
    bs = defaultdict(list)
    bd = defaultdict(list)
    ys = defaultdict(list)
    yd = defaultdict(list)
    
    # keep track of what kmers we've seen to avoid re-analyzing and 
    # inserting kmers we've seen before
    kmer_tracker = defaultdict(str)
    
    for i, prot_name in enumerate(db.proteins):

        db.verbose and print(f'Looking at protein {i + 1}/{len(db.proteins)}\r', end='')

        prot_entry: DatabaseEntry = db.proteins[prot_name]

        # add the protein to the tree
        db.tree.add(prot_entry.name, prot_entry.sequence)

        # go through the kmers for the b sequences
        for j in range(len(prot_entry.sequence) - db.min_len):

            # make a kmer sequence. Do the max (to generate the kmer spec once) then 
            # just iterate through it
            kmer_len = db.max_len if j + db.max_len <= len(prot_entry.sequence) \
                else len(prot_entry.sequence) - j
            kmer = prot_entry.sequence[j:j+kmer_len]

            # generate the singly and doubly b spectra
            kmer_spec_b_s = gen_spectrum(kmer, ion='b', charge=1)['spectrum']
            kmer_spec_b_d = gen_spectrum(kmer, ion='b', charge=2)['spectrum']

            # iterate through the spectra and add the entry to the table
            for k in range(db.min_len, kmer_len):

                if 'b' in kmer_tracker[kmer[:k]]:
                    continue

                kmer_tracker[kmer[:k]] += 'b'

                bs[math.floor(kmer_spec_b_s[k-1])].append(MassSequence(kmer_spec_b_s[k-1], kmer[:k]))
                bd[math.floor(kmer_spec_b_d[k-1])].append(MassSequence(kmer_spec_b_d[k-1], kmer[:k]))
            
        # go through the kmers for the b sequences
        for j in range(len(prot_entry.sequence) - db.min_len):

            # make a kmer sequence. Do the max (to generate the kmer spec once) then 
            # just iterate through it
            kmer_len = db.max_len if j + db.max_len <= len(prot_entry.sequence) \
                else len(prot_entry.sequence) - j

            kmer = prot_entry.sequence[-j - kmer_len: -j] if j != 0 else prot_entry.sequence[-kmer_len:]

            # generate the singly and doubly b spectra
            kmer_spec_y_s = gen_spectrum(kmer, ion='y', charge=1)['spectrum']
            kmer_spec_y_d = gen_spectrum(kmer, ion='y', charge=2)['spectrum']

            # iterate through the spectra and add the entry to the tayle
            for k in range(db.min_len, kmer_len):

                if 'y' in kmer_tracker[kmer[-k:]]:
                    continue

                kmer_tracker[kmer[-k:]] += 'y'

                ys[math.floor(kmer_spec_y_s[k-1])].append(MassSequence(kmer_spec_y_s[k-1], kmer[-k:]))
                yd[math.floor(kmer_spec_y_d[k-1])].append(MassSequence(kmer_spec_y_d[k-1], kmer[-k:]))

    db = db._replace(kmer_masses=KmerMasses(bs, bd, ys, yd))
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
    # get the dictionary to search
    kmer_dict = db.kmer_masses.bs if kmers == 'bs' else (db.kmer_masses.bd if kmers == 'bd' \
                else (db.kmer_masses.ys if kmers == 'ys' else db.kmer_masses.yd))

    # go through each mass in the observed and get any mass in the tolerance
    hits = []
    for mass in observed.spectrum:

        # get the tolerance in Da
        tol = ppm_to_da(mass, tolerance)

        # get upper and lower bound keys and masses
        lb_mass = mass - tol
        ub_mass = mass + tol
        lb_mass_key = math.floor(lb_mass)
        ub_mass_key = math.floor(ub_mass)

        # go through all the values in thelower bound key and collect hits that are in our tolerance
        hits += [x.sequence for x in kmer_dict[lb_mass_key] if lb_mass <= x.mass <= ub_mass]

        # if our upper bound key is different, also search the upper bound key
        if lb_mass_key != ub_mass_key:
            hits += [x.sequence for x in kmer_dict[ub_mass_key] if lb_mass <= x.mass <= ub_mass]
            
    return hits
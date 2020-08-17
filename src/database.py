from pyteomics import fasta
from collections import namedtuple

from src.objects import Database

def extract_protein_name(prot_entry: namedtuple) -> str:
    '''
    Extract the protein name from a protein entry namedtuple from pyteomics 
    fasta read. 
    
    Inputs:
        prot_entry:     (namedtuple) the namedtuple with the entry 'description'
    Outputs:
        (str) the name of the protein to extract
    '''
    if '|' in prot_entry.description:
        return prot_entry.description.split('|')[-1].split(' ')[0]

    return prot_entry.description.split(' ')[0]

def build(fasta_file: str) -> Database:
    '''
    Create a Database namedtuple from a fasta file

    Inputs:
        fasta_file:     (str) the full path to a source database
    Outputs:
        (Database) a Database namedtuple with the fasta_file and proteins
                    feilds filled and an initialized tree
    '''
    db = Database(fasta_file)

    for entry in fasta.read(fasta_file):
        db.proteins[entry.description] = entry

    return db


def get_proteins_with_subsequence(db: Database, sequence: str) -> list:
    '''
    Find the name of all proteins that have the subsequence provided. A 
    list of these names are returned

    Inputs:
        db:         (Database) holder of the sequences
        sequence:   (str) the subsequence to look for
    Outputs:
        (list) string names of the source proteins
    '''
    return db.tree.search(sequence)

def get_entry_by_name(db: Database, name: str) -> namedtuple:
    '''
    Get a namedtuple of the protein entry from the database. 

    Inputs:
        db:     (Database) the database with the proteins
        name:   (str) the name of the protein to find
    Ouputs:
        (namedtuple) has entries 'description' and 'sequence'
    '''
    return db.proteins[name]
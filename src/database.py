from pyteomics import fasta
from collections import namedtuple, defaultdict

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

    prots = defaultdict(list)

    # pull the name out
    get_name = lambda x: x.split('|')[-1].split()[0]

    for entry in fasta.read(fasta_file):
        p_name = get_name(entry.description)
        prots[p_name].append(entry)

    db = db._replace(proteins=prots)
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
    return list(set(db.kmers[sequence]))

def get_proteins_with_subsequence_ion(db: Database, sequence: str, ion: str) -> list:
    '''
    Find all protein names that have the subsequence. Recursivley search
    if the full sequence is not found immediately

    Inputs:
        db:         (Database) holder of the sequences
        sequence:   (str) the subsequence to look for
        ion:        (str) the ion type. {'b', 'y'}
    Outputs:
        (list) string names of the source proteins
    '''

    hits = []
    subseq = sequence

    while len(hits) == 0 and len(subseq) > 0:

        # get some hits
        hs = get_proteins_with_subsequence(db, subseq)

        # make sure that these hits HAVE the full sequence
        for h in hs:
            for entry in get_entry_by_name(db, h):
                if sequence in entry.sequence:
                    hits.append(h)

        # if no hits, reduce the subseq according to the ion
        if len(hits) == 0:
            subseq = subseq[1:] if ion == 'y' else subseq[:-1]

    return hits

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
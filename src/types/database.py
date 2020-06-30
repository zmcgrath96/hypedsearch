from datrie import Trie
import string
from src.types.objects import KmerMetaData, DatabaseEntry
from collections import defaultdict
from pyteomics import fasta
    
class Database: 
    '''
    Database

    Container for multiple database operations. Packages 
    '''
    def __init__(
        self, 
        fasta_file_name='', 
        is_uniprot=False, 
        min_len=3, 
        max_len=20, 
        verbose=False
    ) -> None:
        '''
        Init the database with a fasta file name

        Inputs:
        kwargs:
            fasta_file_name:    str name of fasta file to read
            is_uniprot:         bool Set to True if the database is from UnitPort. If True, extra informtion is kept
            min_len:          int size of kmer to use for indexing the database
            verbose:        (bool) extra printing. Default=False

        Outputs: 
            None
        '''
        self.fasta_file: str = fasta_file_name
        self.proteins: dict = self.__read_fasta(fasta_file_name) if '.fasta' in fasta_file_name else {}
        self.min_len: int = min_len
        self.max_len: int = max_len
        self.metadata: dict = None
        self.verbose: bool = verbose
        self.tree: Trie = self.__build_tree()

    ##################### Overloaded operators ################

    def __iter__(self):
        for _, entry in self.proteins.items():
            yield entry

    def __len__(self):
        return len(self.proteins)
    
    ##################### Private Methods #####################
    def __build_tree(self) -> Trie:
        '''
        Add a prefix tree to the database
        '''
        t = Trie(string.ascii_uppercase)
        plen = len(self.proteins)
        i = 0
        for key, value in self.proteins.items():
            self.verbose and print(f'Adding protein {i + 1}/{plen} to tree\r', end='')
            i += 1

            # add each subsequence to the tree
            for j in range(len(value.sequence) - self.min_len):
                subseqlen = self.max_len if j + self.max_len < len(value.sequence) - 1 else len(value.sequence) - j
                t[value.sequence[j:j+subseqlen]] = key
    
        return t


    def __read_fasta(self, fasta_file: str, is_uniprot=False) -> dict:
        '''
        Read proteins into memory from fasta file
        
        Inputs:
            fasta_file: str path to fasta file
        kwargs:
            is_uniprot: bool adds attribute 'human_readable_name' to dictionary if True. Default=False
        Outputs:
            dictionary of Entry class keyed by the protein name 
        '''
        prots = {}
        e: DatabaseEntry

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

        return prots

    ##################### Public Methods #####################
    def set_max_len(self, k: int) -> None:
        '''
        Change the max size of the database for indexing

        Inputs:
            k:      (int) new max size 
        Outputs:
            None
        '''
        self.max_len = k
        
    def get_entry_by_name(self, name: str) -> DatabaseEntry:
        '''
        Get a protein entry from the database

        Inputs: 
            name:       str name of the protein to get
        Outputs:
            Entry class of the protein. Empty Entry if not found
        '''
        return self.proteins.get(name, DatabaseEntry())

    def add_entry(
        self, 
        protein_name: str, 
        protein_sequnece: str, 
        protein_id='', 
        human_readable_name=''
    ) -> bool:
        '''
        Add an entry to the database. If it cannot be added to the database, False is returned

        Inputs:
            protein_name:       (str) name of the protein being added
            protein_sequence:   (str) string of amino acids describing the protein
        kwargs:
            protein_id:         (str) identifier of the protein. Default=''
            human_readable_name:(str) name that is easy for people to read. Default=''
        Outputs:
            bool True if added, False if it can't be
        '''
        if protein_name in self.proteins:
            return False

        e = DatabaseEntry(protein_name, protein_sequnece, protein_id, human_readable_name)
        self.proteins[protein_name] = e 

        # add it to the tree
        if not self.tree:
            self.tree = self.__build_tree()

        else:
            for i in range(len(protein_sequnece) - self.min_len):
                subseqlen = self.max_len if i + self.max_len < len(protein_sequnece) -1 else len(protein_sequnece) - i
                self.tree[protein_sequnece[i:i+subseqlen]] = protein_name

        return True
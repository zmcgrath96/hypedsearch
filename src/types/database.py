from suffix_tree import Tree
from src.types.objects import KmerMetaData

class Entry(object):
    '''
    Class to contain protein entry information
    '''
    def __init__(self, name: str, sequence: str, readable_name='', id='', kmer_size=3):
        '''
        init

        Inputs:
            name:       string name of the protein
            sequence:   string the amino acid sequence of the protein
        kwargs:
            readable_name:  string the easy to read name of the protein
            id:             id of the protein from the source
            kmer_size:      int the size of kmer to use for indexing
        '''
        self.name: str = name
        self.sequence: str = sequence
        self.readable_name: str = readable_name
        self.id: str = id
        self.kmers: list = None
        self.kmer_size: int = None

    ##################### Setters #####################
    def set_kmer_size(self, k: int):
        self.kmer_size = k

    ##################### Public Methods #####################
    def index(self):
        '''
        Index a protein by generating pairs from the protein of (kmer, start_position, end_position (inclusive))
        
        Outputs:
            list of tuples in the form (kmer, start_position, end_position)
        '''
        kmers = []
        for i in range(len(self.sequence) + self.kmer_size -1):
            entry = (self.sequence[i: i+self.kmer_size], i, i + self.kmer_size-1)
            kmers.append(entry)
        self.kmers = kmers
        return self.kmers

class Database(object): 
    '''
    Database

    Container for multiple database operations. Packages 
    '''
    def __init__(self, fasta_file_name='', is_uniprot=False, kmer_size=3, verbose=False) -> None:
        '''
        Init the database with a fasta file name

        Inputs:
        kwargs:
            fasta_file_name:    str name of fasta file to read
            is_uniprot:         bool Set to True if the database is from UnitPort. If True, extra informtion is kept
            kmer_size:          int size of kmer to use for indexing the database
            verbose:        (bool) extra printing. Default=False

        Outputs: 
            None
        '''
        self.fasta_file: str = fasta_file_name
        self.proteins: dict = self.__read_fasta(fasta_file_name) if '.fasta' in fasta_file_name else {}
        self.kmer_size: int = kmer_size
        self.metadata: dict = None
        self.verbose: bool = verbose
        self.tree: Tree = self.__build_tree()
    
    ##################### Private Methods #####################
    def __build_tree(self) -> Tree:
        '''
        Add a suffix tree to the database
        '''
        t = Tree()
        plen = len(self.proteins)
        printskiplen = plen // 100
        printskipc = 0
        i = 0
        for key, value in self.proteins.items():
            if printskipc == printskiplen:
                printskipc = 0
                self.verbose and print(f'Adding protein {i + 1}/{plen} to tree\r', end='')
            
            i += 1
            printskipc += 1

            t.add(key, value.sequence)
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
        with open(fasta_file, 'r') as i:
            name = None 
            seq = '' 
            identifier = ''
            hmn_rdble_name = ''
            for line in i:
                if '>' in line: #name line

                    # add the last thing to the list
                    if not ((name is None or name == '') and (seq is None or seq == '')):
                        entry = {
                            'sequence': seq,
                            'identifier': identifier
                        }
                        entry['human_readable_name'] = hmn_rdble_name if is_uniprot else ''
                        e = Entry(name, entry['sequence'], entry['human_readable_name'], entry['identifier'])
                        prots[name] = e

                    seq = '' 
                    name = str(str(line.split('|')[2]).split(' ')[0]).replace('\n', '')
                    identifier = str(line.split('|')[1])
                    if is_uniprot:
                        after_bar = str(line.split('|')[2])
                        hmn_rdble_name = str(' '.join(after_bar.split(' ')[1:]).split('OS=')[0]).strip()
                else:
                    seq += line.replace('\n', '')
            # add the last one
            entry = {
                'sequence': seq,
                'identifier': identifier
            }
            entry['human_readable_name'] = hmn_rdble_name if is_uniprot else ''
            e = Entry(name, entry['sequence'], entry['human_readable_name'], entry['identifier'])
            prots[name] = e
        return prots

    ##################### Public Methods #####################
    def set_kmer_size(self, k: int) -> None:
        '''
        Change the k-mer size of the database for indexing

        Inputs:
            k:      (int) new kmer size 
        Outputs:
            None
        '''
        self.kmer_size = k

    def index(self) -> None:
        '''
        Create an indexing of a database by the kmer size. The indexing information 
        is saved in the self.metadata attribute as Kmer namedtuple instances

        Inputs:
            None
        Outputs: 
            None
        '''
        kmers = {}
        for name, e in self.proteins.items():
            # e is an entry class
            e.set_kmer_size(self.kmer_size)
            # index the entry and add it to my index
            mers = e.index()
            for i in range(len(mers)):
                mer, start_pos, end_pos = mers[i]
                if mer not in kmers:
                    kmers[mer] = []
                pairing = KmerMetaData(name, start_pos, end_pos)
                kmers[mer].append(pairing)
        # set my metadata and entries
        self.verbose and print('{} unique kmers'.format(len(kmers.keys())))
        self.metadata = kmers
        
    def get_entry_by_name(self, name: str) -> Entry:
        '''
        Get a protein entry from the database

        Inputs: 
            name:       str name of the protein to get
        Outputs:
            Entry class of the protein. Empty Entry if not found
        '''
        if name in self.proteins: 
            return self.proteins[name]
        return Entry('', '')

    def add_entry(self, protein_name: str, protein_sequnece: str, protein_id='', human_readable_name='') -> bool:
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

        e = Entry(protein_name, protein_sequnece, human_readable_name, protein_id, self.kmer_size)
        self.proteins[protein_name] = e 

        if self.metadata:
            # we need to index this bad boi
            mers = e.index()
            for i in range(len(mers)):
                mer, start_pos, end_pos = mers[i]
                if mer not in self.metadata:
                    self.metadata[mer] = []
                pairing = (protein_name, start_pos, end_pos)
                self.metadata[mer].append(pairing)

        # add it to the tree
        if self.tree:
            self.tree.add(protein_name, protein_sequnece)
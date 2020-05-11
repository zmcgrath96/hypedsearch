from src.database.entry import Entry

class Database(object): 
    '''
    Database

    Container for multiple database operations. Packages 
    '''
    def __init__(self, fasta_file_name: str, is_uniprot=False, kmer_size=3) -> None:
        '''
        Init the database with a fasta file name

        Inputs:
            fasta_file_name:    str name of fasta file to read
        kwargs:
            is_uniprot:         bool Set to True if the database is from UnitPort. If True, extra informtion is kept
            kmer_size:          int size of kmer to use for indexing the database
        Outputs: 
            None
        '''
        self.fasta_file = fasta_file_name
        self.proteins = self.__read_fasta(fasta_file_name)
        self.kmer_size = kmer_size
        self.metadata = None
    
    ##################### Private Methods #####################
    def __read_fasta(self, fasta_file: str, is_uniprot=False) -> list:
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
    def index(self):
        '''
        Create an indexing of a database by the kmer size. The indexing information 
        is saved in the self.metadata attribute

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
                pairing = (name, start_pos, end_pos)
                kmers[mer].append(pairing)
        # set my metadata and entries
        print('{} unique kmers'.format(len(kmers.keys())))
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

    def get_metadata_of_mer(self, mer: str) -> list:
        '''
        Get the metadata associated with a mer. This data will be a list of tuples of the form
            (protein name: str, starting_position: int, ending_position: int)

        Inputs: 
            mer:    str the subsequence 
        Outputs:
            list fo tuples. If the mer is not found, and empty list is returned
        '''
        if mer not in self.metadata:
            return []
        return self.metadata[mer]
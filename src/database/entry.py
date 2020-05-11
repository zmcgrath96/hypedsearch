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
        self.name = name
        self.sequence = sequence
        self.readable_name = readable_name
        self.id = id
        self.kmers = None
        self.kmer_size = None


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
            entry = (self.sequence[i: i+self.kmer_size], i, i + self.kmer_size - 1)
            kmers.append(entry)
        self.kmers = kmers
        return self.kmers
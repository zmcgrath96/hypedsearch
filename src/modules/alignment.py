from src.modules.spectrum import Spectrum
from src.modules.kmer import Kmer 
from src.modules.scores import Scores

class Aligned:
    '''
    Class that holds one possible alignement. No spectral information is saved here

    Public Methods:

    Properties:
        b_alignment:    Scores object used to make predicted alignment on the left
        y_alignment:    Scores object used to make predicted alignment on the right
        confidence:     (float) score of the confidence scoring algorithm
        hybrid:         (bool) determines whether an alignment predicts a hybrid peptide
    '''
    def __init__(self, b_score: Scores, y_score: Scores, spectrum: Spectrum):
        self.b_alignment = b_score
        self.y_alignment = y_score
        self.confidence = self.__alignment_score(b_score, y_score, spectrum)
        self.hybrid = b_score.kmer.protein != y_score.kmer.protein

    ################################### Private Methods ###################################
    def __alignment_score(self, b_entry: Scores, y_entry: Scores, spectrum: Spectrum) -> float:
        '''
        CREATED 11 MAY 2020
        Recreate the spectrum by adding the two sequences and scoring it on the spectrum

        Inputs:
            (b_entry, y_entry):     Scores object
            spectrum:               Spectrum object of the spectrum investigated
        Outputs:
            float the score calculated by the above formula
        '''
        spliced_sequence = Kmer(b_entry.kmer.sequence + y_entry.kmer.sequence, '', 0)
        spliced_score = Scores(spectrum, spliced_sequence)
        return spliced_score.b_score + spliced_score.y_score


class Aligner:
    '''
    Class that holds alignment information on a spectrum. 

    Public Methods:
        make_alignment:     make alignment predictions
        add_scores:         add score objects to make the alignment

    Properties
        spectrum:           (Spectrum) Spectrum object of the spectrum being identified
        alignments:         (list) of Aligned objects
        b_scores:           (list) Score objects used to create alignemnts from the left
        y_scores:           (list) Score objects used to create alignments from the right

    '''
    
    def __init__(self, spectrum: Spectrum):
        self.spectrum = spectrum
        self.alignments = []
        self.b_scores = []
        self.y_scores = []

    ################################ Public Methods ################################
    def make_alignment(self, n=3) -> list:
        '''
        Try and align all of the b and y scores given to the Aligner. Returns top n results

        kwargs: 
            n:      (int) number of alignments to create. Default=3
        Outputs:
            (list) list of Aligned objects
        '''
        alignments = []
        for bscore in self.b_scores:
            for yscore in self.y_scores:
                alignments.append(Aligned(bscore, yscore, self.spectrum))
        alignments.sort(key=lambda x: x.confidence, reverse=True)
        self.alignments = alignments
        return_limit = min(n, len(alignments))
        return alignments[:return_limit]

    def add_scores(self, scores: list, ion: str, clear=False) -> None:
        '''
        Add scores to the object

        Inputs:
            scores:     (list) of Scores object to add
            ion:        (str) name of ion these scores are for. Options are ['b', 'y']
        kwargs:
            clear:      (bool) erase the current contents of scores. Default=False
        Outputs:
            None
        '''
        if not len(scores) > 0:
            return
        if type(scores[0]) != type(Scores):
            print('ERROR: Scores should be a list of Scores objects')
            return

        if ion.lower() == 'b':
            if clear:
                self.b_scores = scores 
            else:
                self.b_scores += scores

            self.b_scores.sort(key=lambda x: x.b_score, reverse=True)

        else: 
            if clear:
                self.y_scores = scores
            else:
                self.y_scores = scores

            self.y_scores.sort(key=lambda x: x.y_score, reverse=True)

    ################################ Private Methods ################################
    def __find_protein_pairings(self) -> dict:
        '''
        Find any pairings of the b and y score by parent proteins.

        Inputs:
            (b_score, y_score):    dict of top ion scores with the form 
                {
                    0: {starting_position: int, ending_position: int, k: int, b_score: float, y_score: float, sequence: str, protein_name: str}
                    ...
                    n-1: {...}
                }
        Outputs:
            dict        dictionary of the pairings of the form
                {
                    protein_name: [
                        (b_rank, y_rank, b_score, y_score)  
                    ]
                }
        '''
        pairings = {}

        b_prots = [x.kmer.protein for x in self.b_scores]
        y_prots = [x.kmer.protein for x in self.y_scores]
        # see if any pairing does exist
        if not (any([True if x in y_prots else False for x in b_prots])):
            return pairings

        # at least 1 pairing exists
        for rank, bscore in enumerate(self.b_scores):
            prot = bscore.kmer.protein
            # if the protein is not found in y proteins, skip
            if prot not in y_prots:
                continue
            # if the protein is not in pairings, add it
            if prot not in pairings: 
                pairings[prot] = []
            # find the pairing we know exists and add it
            for yrank, yscore in enumerate(self.y_scores):
                if yscore.kmer.protein == prot:
                    pairings[prot].append((rank, yrank, bscore, yscore))

        return pairings

    ################################ Overridden Methods ################################
    def __iter__(self):
        '''
        Iterates through the aligned object in alignments
        '''
        for x in self.alignments:
            yield x
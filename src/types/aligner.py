from src.types.objects import AlignedScoredKmers, Kmer, ScoredKmer, Spectrum
from src.scoring.scoring import score_subsequence    
from collections import namedtuple

class Aligner(object):
    '''
    Class that holds alignment information on a spectrum. 
    Public Methods:
        make_alignments:    make alignment predictions
        add_scores:         add score objects to make the alignment
        to_json:            function to turn object into a dict for json dumping
    Properties
        spectrum:           (Spectrum) Spectrum namedtuple instance of the spectrum being identified
        alignments:         (list) of AlignedScoredKmer namedtuple instances of aligned kmers
        b_scores:           (list) ScoredKmer namedtuple instances used to create alignments from the left
        y_scores:           (list) ScoredKmer namedtuple instances used to create alignments from the left
    '''
    
    def __init__(self, spectrum: Spectrum):
        self.spectrum = spectrum
        self.alignments = []
        self.b_scores = []
        self.y_scores = []

    ################################ Public Methods ################################
    def make_alignments(self, n=3) -> list:
        '''
        Try and align all of the b and y scores given to the Aligner. Returns top n results
        kwargs: 
            n:      (int) number of alignments to create. Default=3
        Outputs:
            (list) list of AlignedScoredKmers objects
        '''
        alignments = []
        if not self.b_scores or not self.y_scores:
            return []
        for bscore in self.b_scores:
            for yscore in self.y_scores:
                # run an alignment score on the b and y score objects
                aligned = self.__make_alignment(bscore, yscore)
                alignments.append(aligned)
        alignments.sort(key=lambda x: x.alignment_score, reverse=True)
        return_limit = min(n, len(alignments))
        alignments = alignments[:return_limit]
        self.alignments = alignments
        return alignments

    def add_scores(self, scores: list, ion: str, clear=False) -> None:
        '''
        Add scores to the object
        Inputs:
            scores:     (list) of ScoredKmer namedtuple instances
            ion:        (str) name of ion these scores are for. Options are ['b', 'y']
        kwargs:
            clear:      (bool) erase the current contents of scores. Default=False
        Outputs:
            None
        '''
        if not len(scores) > 0:
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

    def to_json(self) -> dict:
        '''
        Turns object into a dictionary recursively
        Inputs:
            None
        Outputs:
            dict
        '''
        return {
            'spectrum': self.__jsonify_namedtuple(self.spectrum),
            'alignments': [self.__jsonify_namedtuple(x) for x in self.alignments],
            'b_scores': [self.__jsonify_namedtuple(x) for x in self.b_scores],
            'y_scores': [self.__jsonify_namedtuple(x) for x in self.y_scores]
        }

    ################################ Private Methods ################################
    def __alignment_score(self, predicted_sequence: str) -> float:
        '''
        CREATED 11 MAY 2020
        Recreate the spectrum by adding the two sequences and scoring it on the spectrum
        Inputs:
            predicted_sequence:     (str) the predicted sequence of amino acids that align to the spectrum
        Outputs:
            (float) sum of the predicted sequence's b and y scores
        '''
        b_score, y_score = score_subsequence(self.spectrum.spectrum, predicted_sequence)
        return b_score + y_score
    
    def __make_alignment(self, b_entry: ScoredKmer, y_entry: ScoredKmer,  missing_overlap=5) -> AlignedScoredKmers:
        '''
        CREATED 12 MAY 2020
        Attempt an alignment between b and y entries

        Inputs:
            b_entry:            (ScoredKmer) the entry to contribute for the left alignment
            y_entry:            (ScoredKmer) the entry to contribute for the right alignment
            spectrum:           (Spectrum) the spectrum this is being aligned to
        kwargs:
            missing_overlap:    (int) the number of missing amino acids between left an right alignment to consider for a hybrid. Defualt=5
        Outputs:
            AlignedScoredKmer
        '''
        hybrid_seq = ''
        # if the same protein
        if b_entry.kmer.protein == y_entry.kmer.protein:
            # check to see if y is left of b
            if y_entry.kmer.start_position < b_entry.kmer.start_position and y_entry.kmer.end_position < b_entry.kmer.end_position:
                spliced_name = '{}~{}~hybrid'.format(b_entry.kmer.protein, y_entry.kmer.protein)
            elif (y_entry.kmer.start_position - b_entry.kmer.end_position) > missing_overlap:
                spliced_name = '{}~{}~hybrid'.format(b_entry.kmer.protein, y_entry.kmer.protein)
            else:
                spliced_name = b_entry.kmer.protein
        # either the proteins aren't the same or they're too far away
        else:
            spliced_name = '{}~{}~hybrid'.format(b_entry.kmer.protein, y_entry.kmer.protein)
        if 'hybrid' not in spliced_name:
            # find the overlap
            spliced_seq = self.__align_overlap(b_entry, y_entry)
        else:
            spliced_seq = b_entry.kmer.sequence + y_entry.kmer.sequence
            hybrid_seq = b_entry.kmer.sequence + '-' + y_entry.kmer.sequence
        
        # get the alignment score
        alignment_score = self.__alignment_score(spliced_seq)
        hybrid = 'hybrid' in spliced_name
        ask = AlignedScoredKmers(b_entry, y_entry, self.spectrum, spliced_name, alignment_score, spliced_seq, hybrid, hybrid_seq)
        return ask
        

    def __align_overlap(self, b_entry: ScoredKmer, y_entry: ScoredKmer) -> str:
        '''
        Align 2 peptide sequences 

        Inputs:
            b_entry:    (ScoredKmer) the left contributer
            y_entry:    (ScoredKmer) the right contributer
        Outputs:
            (str) the aligned subsequence
        '''
        # if there's a gap, fill it with X's for now
        gap_len = y_entry.kmer.start_position - b_entry.kmer.end_position 
        if gap_len > 0:
            return b_entry.kmer.sequence + 'X'*gap_len + y_entry.kmer.sequence
        # check to see if either fully encompasses the other
        if b_entry.kmer.start_position <= y_entry.kmer.start_position and b_entry.kmer.end_position >= y_entry.kmer.end_position:
            return b_entry.kmer.sequence
        elif y_entry.kmer.start_position <= b_entry.kmer.start_position and y_entry.kmer.end_position >= b_entry.kmer.end_position:
            return y_entry.kmer.sequence
        # align the overlap
        else:
            for i in range(len(b_entry.kmer.sequence)):
                if b_entry.kmer.sequence[i:] in y_entry.kmer.sequence:
                    overlapped_seq = b_entry.kmer.sequence[i:]
                    b_cont = b_entry.kmer.sequence[:i]
                    y_cont = y_entry.kmer.sequence[len(overlapped_seq):]
                    return b_cont + overlapped_seq + y_cont

    def __jsonify_namedtuple(self, tup: namedtuple) -> dict:
        '''
        Recursive method to make any nested namedtuple json friendly

        Inputs:
            tup:    (namedtuple) thing to be flattened
        Outputs:
            (dict) json friendly namedtuple
        '''
        td = dict(tup._asdict())
        for key, value in td.items():
            if any([x in str(type(value)) for x in ['ScoredKmer', 'Kmer', 'AlignedScoredKmers', 'Specturm']]):
                td[key] = self.__jsonify_namedtuple(value)
        return td


    ################################ Overridden Methods ################################
    def __iter__(self):
        '''
        Iterates through the aligned object in alignments
        '''
        for x in self.alignments:
            yield x
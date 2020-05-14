from collections import namedtuple

'''
Named tuples for lighter-weight object like interaction
'''

'''
Spectrum:
    Holds information regarding an MS or MS/MS spectrum

    Properties:
        spectrum:       (list) m/z float values from an MS run
        abundance:      (list) floats that describe abundance of each peak value
        ms_level:       (int) MS experiment level
        scan_number:    (int) scan number of the spectrum
        precursor_mass: (float) precursor mass of the spectrum 
        file_name:      (string) name of the file that the spectrum was taken from
'''
Spectrum = namedtuple(
    'Spectrum', 
    ['spectrum', 'abundance', 'ms_level', 'scan_number', 'precursor_mass', 'file_name'],
    defaults=[[], [], 0, -1, 0, '']
)

'''
Kmer:
    holds kmer information for indexing databases and incrementing sequences

    Properties:
        k:              (int) the length of the kmer
        sequence:       (str) the actual kmer
        protein:        (str) the name of the protein the kmer is taken from
        start_position: (int) the starting position from the protein sequence the kmer is taken from
        end_position:   (int) the ending position from the protein sequence the kmer is taken from
'''
Kmer = namedtuple(
    'Kmer', 
    ['k', 'sequence', 'protein', 'start_position', 'end_position'], 
    defaults=[[0, '', '', -1, -1]]
)

'''
BasicScoredKmer:
    holds basic scoring info of 2 scores and a string

    Properties:
        b_score:    (float) the b-ion score of the kmer 
        y_score:    (float) the y-ion score of the kmer
        kmer:       (str) the kmer sequence
'''
BasicScoredKmer = namedtuple(
    'BasicScoredKmer', 
    ['b_score', 'y_score', 'kmer'],
    defaults=[0, 0, '']
)

'''
ScoredKmer:
    holds information regarding a scored kmer

    Properties:
        b_score:    (float) the b-ion score of the kmer 
        y_score:    (float) the y-ion score of the kmer
        kmer:       (Kmer) a Kmer namedtuple instance
'''
ScoredKmer = namedtuple(
    'ScoredKmer', 
    ['b_score', 'y_score', 'kmer'],
    defaults=[0, 0, Kmer(0, '', '', -1, -1)]
)

'''
AlignedScoredKmers:
    holds score information of two scored kmers and an attempted alignment

    Properties:
        b_alignment:    (Kmer) the Kmer namedtuple instance associated with the left alignment
        y_alignment:    (Kmer) the Kmer namedtuple instance associated with the right alignment
        spectrum:       (Spectrum) the Spectrum namedtuple instance this alignment was created for
        protein:        (str) the name of the protein this is assumed to be from. If hybrid the name will be <left protein>-<right protein>-hybrid
        b_score:        (float) the b_score associated with the sequence compared to some spectrum
        y_score:        (float) the y_score associated with the sequence compared to some spectrum
        sequence:       (str) the predicted amino acid sequence 
        hybrid:         (bool) True if predicted sequence is a hybrid false otherwise
        hybrid_sequence:(str) if the sequence is a hybrid, a sequence is shown here with a - where the junction is
'''
AlignedScoredKmers = namedtuple(
    'AlignedScoreKmers',
    ['b_alignment', 'y_alignment', 'spectrum', 'protein', 'alignment_score', 'sequence', 'hybrid', 'hybrid_sequence'],
    defaults=[ScoredKmer(0, 0,  Kmer(0, '', '', -1, -1)), ScoredKmer(0, 0,  Kmer(0, '', '', -1, -1)), Spectrum([], [], 0, -1, 0, ''), '', 0, '', False, '']
)
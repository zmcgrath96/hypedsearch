from collections import namedtuple

'''
Named tuples for lighter-weight object like interaction
'''

'''
KmerMassesResults:
    Holds the hits from a hash on the entries of KmerMassesResults.

    Properties:
        (bs, bd, ys, yd):   (list) MassSequence hits 
'''
KmerMassesResults = namedtuple(
    'KmerMassesResults', 
    ['bs', 'bd', 'ys', 'yd'], 
    defaults=[[], [], [], []]
)

'''
__cache:
    Holds protein and kmer hits that were in the sqlite database. 

    Properties:
        proteins:           (DataFrame) cached proteins 
        kmers:              (DataFrame) cached proteins
        max_mass:           (float) the highest m/z value we cached
        proteins_cached:    (bool) proteins were cached into memory
'''
cache__ = namedtuple(
    'cache__',
    ['proteins', 'kmers', 'max_mass', 'proteins_cached'], 
    defaults=[None, None, 0, False]
)

'''
Database:
    Holds proteins, fasta file, protein tree, and kmer masses 

    Properties:
        fasta_file:     (str) the name of the input fasta file
        min_len:        (int) the minimum length peptide to consider
        max_len:        (int) the maximum length peptide to consider
        verbose:        (bool) extra printing
        conn:           (sqlite3.Connection) the connection to the sqlite database
        __cache:        (__cache) data to be held in memory
'''
Database = namedtuple(
    'Database', 
    ['fasta_file', 'min_len', 'max_len', 'verbose', 'conn', 'cache__'], 
    defaults=['', 0, 0, True, None, None]
)

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
    ['spectrum', 'abundance', 'total_intensity', 'ms_level', 'scan_number', 'precursor_mass', 'file_name', 'id'],
    defaults=[[], [], 0, 0, -1, 0, '', '']
)

'''
MassSequence:
    Basic container for holding some mass associated with a sequence

    Properties: 
        mass:       (float) some mass associated with a sequence
        sequence:   (str) the sequence associated with the mass
'''
MassSequence = namedtuple(
    'MassSequence', 
    ['mass', 'sequence'], 
    defaults=[0.0, '']
)

'''
SequenceAlignment:
    Information on a nonhybrid peptide alignment. 

    Properties:
        proteins:       (list) proteins that this sequence is found in
        sequence:       (str) Amino acid sequence tagged
        b_score:        (float) b ion score of the sequence
        y_score:        (float) y ion score of the sequence
        total_score:    (float) score given to the sequence
'''
SequenceAlignment = namedtuple(
    'SequenceAlignment', 
    ['proteins', 'sequence', 'b_score', 'y_score', 'total_score', 'precursor_distance'],
    defaults=[[], '', 0.0, 0.0, 0.0, 100]
)

'''
HybridSequenceAlignment
    Information on a hybrid sequence alignment

    Properties:
        left_proteins:      (list) proteins that contain the sequence of amino acids that contribute
                                   to the left side of the hybrid peptide
        right_proteins:     (list) proteins that contain the sequence of amino acids that contribute 
                                   to the right side of the hybrid peptide
        sequence:           (str) Amino acid sequence tagged
        hybrid_sequence:    (str) Amino acid sequence tagged with indicators for the junction area
        b_score:            (float) b ion score of the sequence
        y_score:            (float) y ion score of the sequence
        total_score:        (float) the score given to the sequence
'''
HybridSequenceAlignment = namedtuple(
    'HybridSequenceAlignment', 
    ['left_proteins', 'right_proteins', 'sequence', 'hybrid_sequence', 
        'b_score', 'y_score', 'total_score', 'precursor_distance'],
    defaults=[[], [], '', '', 0.0, 0.0, 0.0, 100]
)

'''
Alignments
    Contains the spectrum with SequenceAlignments and HybridSequenceAlignments

    Properties:
        spectrum:       (Spectrum) The spectrum being aligned
        alignments:     (list) contains both HybridSequenceAlignment and SequenceAlignment 
'''
Alignments = namedtuple(
    'Alignments', 
    ['spectrum', 'alignments'], 
    defaults=[Spectrum([], [], 0, 0, 0.0, ''), []]
)
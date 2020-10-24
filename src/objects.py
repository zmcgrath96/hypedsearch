from collections import namedtuple

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

'''
Digest
    contains basic information of what a digest does

    Properties:
        start:      (list) tuples of the form (char, char) whre the first entry is what 
                    amino acid could start a valid digest and the second is either
                    L to cut to the left of it or R to cut to the right of it
        end:        (list) tuples of the form (char, char) whre the first entry is what 
                    amino acid could start a valid digest and the second is either
                    L to cut to the left of it or R to cut to the right of it
'''
digest = namedtuple(
    'digest',
    ['start', 'end'],
    defaults=[[], []]
)
import unittest
from src.scoring import scoring
from src.spectra.gen_spectra import gen_spectrum

class test_scoring(unittest.TestCase):
    def setUp(self):
        pass 

    def test_score_subsequence(self):
        seq1 = 'MALWA'
        seq2 = 'MALVR'
        seq3 = 'WAQSST'

        spec1 = gen_spectrum(seq1)['spectrum']

        ss11 = scoring.score_subsequence(spec1, seq1)
        ss12 = scoring.score_subsequence(spec1, seq2)
        ss13 = scoring.score_subsequence(spec1, seq3)

        self.assertTrue(ss11 > ss12, 'scoring a sequence against its own spectrum should be greater than an incorrect alignment scoring')
        self.assertTrue(ss11 > ss13, 'scoring a sequence against its own spectrum should be greater than a jumble of AAs')
        self.assertTrue(ss12 > ss13, 'scoring a partial alignment should be greater than a jumble of AAs')
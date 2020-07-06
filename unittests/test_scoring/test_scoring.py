import unittest
from src.scoring import scoring
from src.sequence.gen_spectra import gen_spectrum
from src.objects import Spectrum

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

    def test_backbone_score(self):
        seq = 'MALWARM'
        
        # full coverage from 1 ion should give 100 score
        fullbion = gen_spectrum(seq, ion='b', charge=1)['spectrum']
        fullbionscore = scoring.backbone_score(Spectrum(fullbion), seq, 10)
        self.assertEqual(fullbionscore, 100, 'full coverage of one ion type should give a score of 100')

        # full coverage from 2 ions with no overlap should give 100 score
        halfbs = gen_spectrum(seq[:4], ion='b', charge=1)['spectrum']
        halfys = gen_spectrum(seq[5:], ion='y', charge=1)['spectrum']
        halfbshalfysscore = scoring.backbone_score(Spectrum(halfbs + halfys), seq, 10)
        self.assertEqual(halfbshalfysscore, 100, 'full coverage from 2 ions without overlap should give a score of 100')

        # partial coverage by 1 ion should give a score represented by the integer value of percentage coverage
        partialbs = gen_spectrum(seq[:4], ion='b', charge=1)['spectrum']
        partialbsscore = scoring.backbone_score(Spectrum(partialbs), seq, 10)
        self.assertEqual(partialbsscore, int(100 * len(partialbs)/(len(seq)-1)), 'partial coverage by 1 ion should give a score represented by the integer value of percentage coverage')

        # partial coverage by 2 ions should give a score represented by the integer value of percentage coverage
        partialbs = gen_spectrum(seq[:2], ion='b', charge=1)['spectrum']
        partialys = gen_spectrum(seq[5:], ion='y', charge=1)['spectrum']
        partialbspartialysscore = scoring.backbone_score(Spectrum(partialbs + partialys), seq, 10)
        self.assertEqual(partialbspartialysscore, int(100*(len(partialbs) + len(partialys))/(len(seq)-1)), 'partial coverage by 2 ions should give a score represented by the integer value of percentage coverage')

        # full coverage for 2 ions should give score 100 + number bonds
        fullbs = gen_spectrum(seq, ion='b', charge=1)['spectrum']
        fullys = gen_spectrum(seq, ion='y', charge=1)['spectrum']
        fullbsfullysscore = scoring.backbone_score(Spectrum(fullbs + fullys), seq, 10)
        self.assertEqual(fullbsfullysscore, 100+len(seq)-1, 'full coverage for 2 ions should give score 100 + number bonds')

        # full coverage by all 4 ions should give screo 100 + 3*number bonds
        allspec = gen_spectrum(seq)['spectrum']
        allspecscore = scoring.backbone_score(Spectrum(allspec), seq, 10)
        self.assertEqual(allspecscore, 100 + 3*(len(seq)-1), 'full coverage by all 4 ions should give screo 100 + 3*number bonds')

        # empty list should return score 0
        zerospec = scoring.backbone_score(Spectrum(), seq, 10)
        self.assertEqual(zerospec, 0, 'empty list should return score 0')

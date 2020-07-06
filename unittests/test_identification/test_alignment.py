import unittest
from src.identfication import alignment
from src.database import Database

class test_alignment(unittest.TestCase):
    def setUp(self):
        self.db = Database()
        self.db.add_entry('prot1', 'MALWARMQSSRTVMLK')
        self.db.add_entry('prot2', 'VVQQERLLYWMPYNVVS')

    def test_align_overlaps(self):
        seq1 = 'ABCDE'
        seq2 = 'DEFGH'
        seq3 = 'ABCDEFGH'
        seq4 = 'LMNOP'

        self.assertEqual(alignment.align_overlaps(seq1, seq2), seq3, 'ABCDE and DEFGH should alignm to make ABCDEFGH')
        self.assertEqual(alignment.align_overlaps(seq3, seq1), seq3, 'ABCDEFGH and ABCDE should return the one that contains the subsequence')
        self.assertEqual(alignment.align_overlaps(seq3, seq2), seq3, 'ABCDEFGH and DEFGH should return the one that contains the subsequence')
        self.assertEqual(alignment.align_overlaps(seq1, seq4), seq1 + '-' + seq4, 'Non overlapping sequences should be concatenated')

        self.assertEqual(alignment.align_overlaps(seq2, seq1), 'DEFGH-ABCDE', 'DEFGH and ABCDE should just append the latter to the former')

    def test_hybrid_alignemnt(self):
        seq1 = 'ABCDE'
        seq2 = 'DEFGH'
        seq3 = 'LMNOP'

        self.assertEqual(alignment.hybrid_alignment(seq1, seq2)[1], 'ABC(DE)FGH', 'hybrid alignment with an overlap should put ambiguous section in parenthesis')
        self.assertEqual(alignment.hybrid_alignment(seq1, seq3)[1], seq1 + '-' + seq3, 'hybrid alignment with non overlapping sequences should concatenate with a -')

    def test_align_b_y(self):
        blist = ['MALWARM', 'QSSRT', 'YNVVS']
        ylist = ['MALWARM', 'SRTVMLK', 'ERLL']
        shouldbe = ['MALWARM', 'QSSRTVMLK', 'YNVVS-ERLL']

        self.assertTrue(all([x in alignment.align_b_y(blist, ylist, self.db) for x in shouldbe]), f'{shouldbe} all in alignments from b and y lists')

    def test_get_parents(self):
        seq1 = 'MALWARM'
        seq2 = 'SSRTV-PYNVVS'
        seq3 = 'TVM(L)LYWMP'

        self.assertEqual(alignment.get_parents(seq1, self.db), (['prot1'], None), 'nonhybrid sequence should have 1 correct protein as the parent')
        self.assertEqual(alignment.get_parents(seq2, self.db), (['prot1'], ['prot2']), 'hybrid sequence should return the two correct parents')
        self.assertEqual(alignment.get_parents(seq3, self.db), (['prot1'], ['prot2']), 'ambiguous hybrid sequence should return the two correct parents')

    def test_replace_ambiguous_hybrids(self):
        hybs = ['SSRTV-PYNVVS', 'MQSS(RT)VMLK' ]
        updated_hybs = alignment.replace_ambiguous_hybrids(hybs, self.db)
        
        self.assertEqual(hybs[0], updated_hybs[0], 'True hybrid alignment should not be replaced')
        self.assertEqual(hybs[1].replace('(', '').replace(')', ''), updated_hybs[1], 'hybrid sequence found in database should be replaced with nonhybrid')
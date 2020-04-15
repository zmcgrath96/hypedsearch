import unittest
from src.alignment import search
from src.spectra import gen_spectra

class test_search(unittest.TestCase):
    def setUp(self):
        self.prots = [{
            'sequence': 'MALWARMQRST',
            'name': 'prot1',
            'id': 'id1'
        },
        {
            'sequence': 'QSVYMALCTR', 
            'name': 'prot2',
            'id': 'id2'
        }]
        self.unk_seq = 'ARMQRS'
        self.spectrum = gen_spectra.gen_spectrum(self.unk_seq)['spectrum']

    def test_merge_and_sort(self):
        old_b = {
            0: { 'b_score': 4.4, 'sequence': 'MALWA' },
            1: { 'b_score': 3.7, 'sequence': 'MALW' },
            2: { 'b_score': 1.1, 'sequence': 'AL' }
        }
        new_b = {
            0: { 'b_score': 4.8, 'sequence': 'MALWA' },
            1: { 'b_score': 1.1, 'sequence': 'AFK' },
            2: { 'b_score': 0.1, 'sequence': 'Q' }
        }
        old_y = {
            0: { 'y_score': 4.8, 'sequence': 'AWLAM' },
            1: { 'y_score': 1.0, 'sequence': 'M' },
            2: { 'y_score': 0.2, 'sequence': 'F' }
        }
        new_y = {
            0: { 'y_score': 5.1, 'sequence': 'SAWLAM' },
            1: { 'y_score': 4.9, 'sequence': 'AFKLA' },
            2: { 'y_score': 0.3, 'sequence': 'A' }
        }
        
        updated_b, updated_y = search.merge_and_sort(old_b, old_y, new_b, new_y)

        self.assertEqual([new_b[0], old_b[0], old_b[1]], [x for _, x in updated_b.items()], 'the merge and sort should return the top 3 scores')
        self.assertEqual([new_y[0], new_y[1], old_y[0]], [x for _, x in updated_y.items()], 'the merge and sort should return the top 3 scores')

    def test_search_protein(self): 
        b_results, y_results = search.search_protein({'spectrum': self.spectrum}, self.prots[0])
        self.assertEqual(b_results[0]['sequence'], self.unk_seq, 'top result for b result and ideal spectrum should be the correct sequence')
        self.assertEqual(y_results[0]['sequence'], self.unk_seq, 'top result for y result and ideal specturm should be the correct sequence')

        bad_b_res, bad_y_res = search.search_protein({'spectrum': self.spectrum}, self.prots[1])
        self.assertTrue(b_results[0]['b_score'] > bad_b_res[0]['b_score'], 'the b score for the correct protein should be better than the wrong protein')
        self.assertTrue(y_results[0]['y_score'] > bad_y_res[0]['y_score'], 'the y score for the correct protein should be better than the wrong protein')

    def test_serach_proteins(self):
        b_results, y_results = search.search_proteins({'spectrum': self.spectrum}, self.prots)
        self.assertEqual(b_results[0]['sequence'], self.unk_seq, 'top b result for all proteins searched and ideal spectrum should be the correct sequence')
        self.assertEqual(y_results[0]['sequence'], self.unk_seq, 'top y result for all proteins searched and ideal spectrum should be the correct sequence')

if __name__ == '__main__':
    unittest.main()
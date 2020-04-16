import unittest
from src.alignment import aligners

class test_aligners(unittest.TestCase):
    def setUp(self): 
        self.b_scores = {
            0: {'starting_position': 5, 'ending_position': 10, 'b_score': 6, 'y_score': 0, 'protein_name': 'protein1'},
            1: {'starting_position': 30, 'ending_position': 33, 'b_score': 2, 'y_score': 0, 'protein_name': 'protein1' },
            2: {'starting_position': 100, 'ending_position': 102, 'b_score': 1, 'y_score': 0, 'protein_name': 'protein2'}
        }
        self.y_scores = {
            0: {'starting_position': 33, 'ending_position': 45, 'b_score': 0, 'y_score': 12, 'protein_name': 'protein2'},
            1: {'starting_position': 5, 'ending_position': 10, 'b_score': 0, 'y_score': 6, 'protein_name': 'protein1'},
            2: {'starting_position': 101, 'ending_position': 102, 'b_score': 0, 'y_score': 1, 'protein_name': 'protein2'}
        }
        self.no_intersection = {
            0: {'starting_position': 0, 'ending_position': 5, 'b_score': 9, 'y_score': 9, 'protein_name': 'protein3'},
            1: {'starting_position': 10, 'ending_position': 15, 'b_score': 9, 'y_score': 9, 'protein_name': 'protein3'},
            2: {'starting_position': 100, 'ending_position': 105, 'b_score': 9, 'y_score': 9, 'protein_name': 'protein3'}
        } 

    def test_find_protein_pairings(self):
        pairings = aligners.find_protein_pairings(self.b_scores, self.y_scores)
        self.assertEqual(len(pairings), 2, 'there should be 4 parirings with the chosen proteins')

        pairings = aligners.find_protein_pairings(self.b_scores, self.no_intersection)
        self.assertEqual(len(pairings), 0, 'there should be no pairings when no proteins are shared')

    def test_align_spectrum_by_protein_ions(self):
        aligned = aligners.align_spectrum_by_protein_ions([1, 2, 3], self.b_scores, self.y_scores)
        self.assertEqual(aligned[0]['confidence'], 100.0, 'the confidence of the best match should be 100 for this case')
        self.assertEqual(aligned[0]['protein_name'], 'protein1', 'the best match should have the protein "protein1"')

        aligned = aligners.align_spectrum_by_protein_ions([1, 2, 3], self.b_scores, self.no_intersection)
        self.assertEqual(len(aligned), 0, 'with no shared proteins, no alignment should be returned')
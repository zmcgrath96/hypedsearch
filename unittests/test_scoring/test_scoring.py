import unittest
from src.scoring import scoring

class test_scoring(unittest.TestCase):
    def setUp(self):
        pass 

    def test_confidence(self):
        b_entry = {
            'b_score':  3,
            'starting_position': 3,
            'ending_position': 5,
        }
        y_entry = {
            'y_score': 3, 
            'starting_position': 3,
            'ending_position': 5
        }
        #((# overlapped_b_ions * b_score) + (# overlapped_y_ions * y_score) + # total ions identified) - (# non overlapped_b_ions + # non overlapped_y_ions)
        confidence_perfect = 9 + 9 + 6
        confidence_score = scoring.confidence(b_entry, y_entry)
        self.assertEqual(confidence_perfect, confidence_score, 'The perfect alignment should give the perfect confidence score')

        b_entry['ending_position'] = 4
        y_entry['starting_position'] = 4
        confidence_slightly_off = (1 * 3) + (1 * 3) + 4 - 2
        confidence_score = scoring.confidence(b_entry, y_entry)
        self.assertEqual(confidence_slightly_off, confidence_score, 'the slighty off alignment should give a slightly lower score')

        y_entry['starting_position'] = 6
        y_entry['ending_position'] = 8
        confidence_way_off = 0 + 0 + 5 - 5
        confidence_score = scoring.confidence(b_entry, y_entry)
        self.assertEqual(confidence_way_off, confidence_score, 'a completely missaligned b and y score should give a score of 0')

    def test_confidence_simple(self): 
        b_entry = {
            'b_score':  3,
            'starting_position': 3,
            'ending_position': 5,
        }
        y_entry = {
            'y_score': 3, 
            'starting_position': 3,
            'ending_position': 5
        }
        #((# overlapped_b_ions * b_score) + (# overlapped_y_ions * y_score) + # total ions identified) - (# non overlapped_b_ions + # non overlapped_y_ions)
        confidence_perfect = 100.0
        confidence_score = scoring.confidence_simple(b_entry, y_entry)
        self.assertEqual(confidence_perfect, confidence_score, 'The perfect alignment should give the perfect confidence score')

        b_entry['ending_position'] = 4
        y_entry['starting_position'] = 4
        confidence_slightly_off = 50.0
        confidence_score = scoring.confidence_simple(b_entry, y_entry)
        self.assertEqual(confidence_slightly_off, confidence_score, 'the slighty off alignment should give a slightly lower score')

        y_entry['starting_position'] = 6
        y_entry['ending_position'] = 8
        confidence_way_off = 0.0
        confidence_score = scoring.confidence_simple(b_entry, y_entry)
        self.assertEqual(confidence_way_off, confidence_score, 'a completely missaligned b and y score should give a score of 0')

import unittest
from src.alignment import comparisons
class test_comapisons(unittest.TestCase):
    def setUp(self):
        self.spectrum = [1, 2, 3, 4, 7]
        self.reference = [1, 2, 3, 4, 5]

    def test_compare_masses(self):
        '''
        equation is (# found in spec + streak #)/(len(reference)/2)
        ref ref = (5 + 5)/(2.5) 
        spec ref = (4 + 4)/2.5
        '''
        spec_ref = comparisons.compare_masses(self.spectrum, self.reference)
        ref_ref = comparisons.compare_masses(self.reference, self.reference)
        empty_ref = comparisons.compare_masses([], self.reference)

        self.assertEqual(spec_ref, 16/5, 'comparison between spectrum and reference should be 16/5')
        self.assertEqual(ref_ref, 4, 'comparison of reference to reference should give 4')
        self.assertEqual(empty_ref, 0, 'comparison of empty to reference should return 0')

    # don't test the other functions as they use compare masses or spectra generation

if __name__ == '__main__':
    unittest.main()
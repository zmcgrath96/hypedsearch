import unittest
from src.scoring import mass_comparisons
class test_comparisons(unittest.TestCase):
    def setUp(self):
        self.spectrum = [1, 2, 3, 4, 7]
        self.reference = [1, 2, 3, 4, 5]

    def test_compare_masses(self):
        '''
        equation is (# found in spec + streak #)/(len(reference))
        ref ref = (5 + 5)/(5) 
        spec ref = (4 + 4)/5
        '''
        spec_ref = mass_comparisons.compare_masses(self.spectrum, self.reference)
        ref_ref = mass_comparisons.compare_masses(self.reference, self.reference)
        empty_ref = mass_comparisons.compare_masses([], self.reference)

        self.assertEqual(spec_ref, 8/5, 'comparison between spectrum and reference should be 8/5')
        self.assertEqual(ref_ref, 2, 'comparison of reference to reference should give 2')
        self.assertEqual(empty_ref, 0, 'comparison of empty to reference should return 0')

        spec_ref_opt = mass_comparisons.optimized_compare_masses(self.spectrum, self.reference)
        ref_ref_opt = mass_comparisons.optimized_compare_masses(self.reference, self.reference)
        empty_ref_opt = mass_comparisons.optimized_compare_masses([], self.reference)

        self.assertEqual(spec_ref_opt, 8/5, 'optimized compare masses of reference and observed should be 8/5')
        self.assertEqual(ref_ref_opt, 2, 'optimized compare masses of reference and reference should be 2')
        self.assertEqual(empty_ref_opt, 0, 'optimized comparsion of empty to reference should return 0')

    # don't test the other functions as they use compare masses or spectra generation

if __name__ == '__main__':
    unittest.main()
import unittest
from src.identfication import filtering
from src.types.objects import Spectrum
from src.sequence.gen_spectra import gen_spectrum
from collections import namedtuple

class test_filtering(unittest.TestCase):
    def setUp(self):
        pass

    def test_slope_filtering(self):
        scores = [10, 9, 8, 7, 6.5, 6, 5.5, 5.1, 4.7, 4.5, 4.4, 4.35, 4.31, 4.299, 4.287]
        self.assertTrue(len(filtering.slope_filtering(scores)) < len(scores)/2, 'filtered slope scores should be reduced by more than half')

    def test_result_filterings(self):
        spec = gen_spectrum('MALWAR')
        s = Spectrum(spec['spectrum'], [100 for _ in range(len(spec['spectrum']))], 2, 0, spec['precursor_mass'], '')
        
        hits = namedtuple('hits', ['b', 'y'])
        h = hits(['MAL', 'MALW', 'MALWAR', 'MLA', 'LMA', 'LAM', 'AML', 'ALM'], ['WAR', 'LWAR', 'MALWAR', 'WRA', 'RAW', 'RWA', 'ARW', 'AWR']) 

        bres, yres = filtering.result_filtering(s, h, 3)
        self.assertEqual(bres[0], 'MALWAR', 'top hit for b side should be the correct sequence')
        self.assertEqual(yres[0], 'MALWAR', 'top hit for y side should be the correct sequence')
    
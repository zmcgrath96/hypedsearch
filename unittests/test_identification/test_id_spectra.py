import unittest
from src.identfication import id_spectra
from src.database import Database
from src.sequence.gen_spectra import gen_spectrum
import math

class test_id_spectra(unittest.TestCase):
    def setUp(self):
        self.prot1seq = 'MALWARM'
        self.prot2seq = 'WQSSRT'
        self.db = Database()
        self.db.add_entry('prot1', self.prot1seq)
        self.db.add_entry('prot2', self.prot2seq)

    def test_build_kmermasses(self):
        # there should be an entry for each mass in the kmermasses entry (setting min to 1)
        spec1 = gen_spectrum(self.prot1seq)['spectrum']
        spec2 = gen_spectrum(self.prot2seq)['spectrum']

        kmermasses = id_spectra.build_kmermasses(self.db, 1, 10)
    
        spec1hits = []
        for mass in spec1:
            hit = False
            for _, entrydict in kmermasses._asdict().items():
                masskey = math.floor(mass)
                if masskey in entrydict:
                    hit = any([mass == x[0] for x in entrydict[masskey]])
                    if hit:
                        break
            spec1hits.append(hit)
        self.assertTrue(all(spec1hits), 'all masses should have been found in the kmermasses object')

        spec2hits = []
        for mass in spec2:
            hit = False
            for _, entrydict in kmermasses._asdict().items():
                masskey = math.floor(mass)
                if masskey in entrydict:
                    hit = any([mass == x[0] for x in entrydict[masskey]])
                    if hit:
                        break
            spec2hits.append(hit)
        self.assertTrue(all(spec2hits), 'all masses should have been found in the kmermasses object')

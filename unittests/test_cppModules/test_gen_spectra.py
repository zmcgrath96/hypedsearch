import unittest
from src.cppModules import gen_spectra
import src.gen_spectra as py_gen_spectra

class test_gen_spectra(unittest.TestCase):
    def setUp(self):
        self.sequence = 'MALWARSMQT'

    def test_b_spectra(self):

        # generate the python one
        py_bs_spec = py_gen_spectra.gen_spectrum(self.sequence, 1, 'b')['spectrum']

        # generate the cpp one
        bs_spec = gen_spectra.gen_spectrum(self.sequence, 'b', 1)
        
        # check that each element is ~equal (2 decimal places)
        for (x1, x2) in zip(sorted(py_bs_spec), sorted(bs_spec)):
            self.assertAlmostEqual(x1, x2, 2, 'Singly b sepctra should be the same')

        # generate the python one
        py_bd_spec = py_gen_spectra.gen_spectrum(self.sequence, 2, 'b')['spectrum']

        # generate the cpp one
        bd_spec = gen_spectra.gen_spectrum(self.sequence, 'b', 2)

        # check that each element is ~equal (2 decimal places)
        for (x1, x2) in zip(sorted(py_bd_spec), sorted(bd_spec)):
            self.assertAlmostEqual(x1, x2, 2, 'Doubly b sepctra should be the same')

    def test_y_spectra(self):

        # generate the python one
        py_ys_spec = py_gen_spectra.gen_spectrum(self.sequence, 1, 'y')['spectrum']

        # generate the cpp one
        ys_spec = gen_spectra.gen_spectrum(self.sequence, 'y', 1)

        # check that each element is ~equal (2 decimal places)
        for (x1, x2) in zip(sorted(py_ys_spec), sorted(ys_spec)):
            self.assertAlmostEqual(x1, x2, 2, 'Singly y sepctra should be the same')

        # generate the python one
        py_yd_spec = py_gen_spectra.gen_spectrum(self.sequence, 2, 'y')['spectrum']

        # generate the cpp one
        yd_spec = gen_spectra.gen_spectrum(self.sequence, 'y', 2)

        # check that each element is ~equal (2 decimal places)
        for (x1, x2) in zip(sorted(py_yd_spec), sorted(yd_spec)):
            self.assertAlmostEqual(x1, x2, 2, 'Doubly y sepctra should be the same')
import unittest
import os
from src import gen_spectra

class test_gen_spectra(unittest.TestCase): 
    def setUp(self): 
        self.sequence = 'MALWAR'
    
    def test_b_ions(self):
        #test the b_ions function with the sequence, "MALWAR".
        #The expected sequence of b ions is [131.040485, 202.077599, 315.161663, 501.240976, 572.27809, 728.379201]
        self.assertEqual(gen_spectra.b_ions(self.sequence), [131.040485, 202.077599, 315.161663, 501.240976, 572.27809, 728.379201])
    
    def test_y_ions(self):
        #Test the y_ions function with the sequence, 'MALWAR".
        #The expected sequence of y ions is [156.101111, 227.138225, 340.222289, 526.301602, 597.338716, 728.379201]
        self.assertEqual(gen_spectra.y_ions(self.sequence), [156.101111, 227.138225, 340.222289, 526.301602, 597.338716, 728.379201])

    def test_calc_masses(self):
        #Test the calc masses function with the sequence, "MALWAR".
        #The masses should be all of the b and y ion masses in one list.
        self.assertEqual(sorted(gen_spectra.calc_masses(self.sequence)), (sorted((gen_spectra.b_ions(self.sequence)) 
            + (gen_spectra.y_ions(self.sequence)))))
    
    def test_max_mass(self):
        #Test the max mass function with the sequence, "MALWAR"
        Expected_max = 728.379201
        charge = 1
        charge2 = 2
        self.assertEqual(gen_spectra.max_mass(self.sequence, 'b', charge), Expected_max)
        self.assertEqual(gen_spectra.max_mass(self.sequence, 'y', charge2), Expected_max / 2)
    
    def test_get_precursor(self):
        #Test the get_precursory function with the sequence, "MALWAR"
        Expected_prec = 728.379201
        self.assertEqual(gen_spectra.get_precursor(self.sequence), 728.379201)
    
    def test_gen_spectrum(self):
        #Test the gen_spectrum function with the sequence, "MALWAR"
        Expected_spectrum = gen_spectra.b_ions(self.sequence) + gen_spectra.y_ions(self.sequence)
        self.assertEqual(gen_spectra.gen_spectrum(self.sequence), Expected_spectrum)
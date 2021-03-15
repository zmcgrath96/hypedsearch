import unittest
import os
from src import gen_spectra

class test_gen_spectra(unittest.TestCase): 
    def setUp(self): 
        sequence = 'MALWAR'
    
    def test_b_ions(self):
        #test the b_ions function with the sequence, "MALWAR".
        #The expected sequence of b ions is [131.040485, 202.077599, 315.161663, 501.240976, 572.27809, 728.379201]
        self.assertEqual(gen_spectra.b_ions(sequence), [131.040485, 202.077599, 315.161663, 501.240976, 572.27809, 728.379201])
    
    def test_y_ions(self):
        #Test the y_ions function with the sequence, 'MALWAR".
        #The expected sequence of y ions is [156.101111, 227.138225, 340.222289, 526.301602, 597.338716, 728.379201]
        self.assertEqual(gen_spectra.y_ions(spectrum), [156.101111, 227.138225, 340.222289, 526.301602, 597.338716, 728.379201])

    def test_calc_masses(self):
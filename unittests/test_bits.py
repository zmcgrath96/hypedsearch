import unittest
import random

from src import bits

class test_bits(unittest.TestCase):
    def setUp(self):
        self.alphabet = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'U', 'W', 'Y', 'V', 'X', 'B', 'Z']
        pass 

    def testEncodeDecodeScale(self):
        
        strings = []
        # generate 1000 strings and check that we can encode and decode correctly
        for i in range(1000):
            str_len = random.randint(3, 30)
            string = ''
            for j in range(str_len):
                string += self.alphabet[random.randint(0, len(self.alphabet)-1)]

        # ensure encoding and decoding works
        are_same = []
        for s in strings:
            e = bits.str_to_bits(s)
            d = bits.bits_to_str(e)

            are_same.append(e == s)

        # check all are true
        self.assertTrue(all(are_same), 'All encoding/decoding returns the same value')
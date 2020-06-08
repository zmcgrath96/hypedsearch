import unittest
import os
from src import utils

class test_utils(unittest.TestCase):
    def setUp(self):
        self.tosortkey = [
            {'a': 'b', 'b': 2}, 
            {'a': 'e', 'b': 18},
            {'a': 'd', 'b': 11},
            {'a': 'a', 'b': -1}, 
            {'a': 'c', 'b': 3}    
        ]
        self.tosortindex = [
            ('b', 2), ('e', 18), ('d', 11), ('a', -1), ('c', 3)
        ]

    def test_insorts(self):
        a = []
        for i in range(len(self.tosortkey)):
            a = utils.insort_by_key(self.tosortkey[i], a, 'b')
        
        self.assertEqual('abcde', ''.join([x['a'] for x in a]), 'insort by key should be sorted correctly')

        b = []
        for i in range(len(self.tosortindex)):
            b = utils.insort_by_index(self.tosortindex[i], b, 1)
        
        self.assertEqual('abcde', ''.join([x[0] for x in b]), 'insort by index should be sorted correctly')

    def test_all_perms_of(self):
        s = 'AABBCC'
        self.assertEqual(16, len(list(set(utils.all_perms_of_s(s, 'AB')))), 'All permutations of A and B should return 16 unique sequences')

    def test_ppm_to_da(self):
        mass = 100
        tol = 20
        self.assertEqual(utils.ppm_to_da(mass, tol), .002, '20 ppm of 100 should be .002')
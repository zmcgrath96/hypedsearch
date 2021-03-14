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

    def test_ppm_to_da(self):
        #Tests the conver ppm to da function with a mass of 100 and a tolerance of 20. 20 ppm of 100 should be .002
        mass = 100
        tol = 20
        self.assertEqual(utils.ppm_to_da(mass, tol), .002, '20 ppm of 100 should be .002')

    def test_file_exists(self)
        #Run the utils.file_exists function with two sample files.
        filename = 'requirements.txt'
        filename2 = 'thisfiledoesnotexist.txt'
        self.assertEqual(utils.file_exists(filename), True) "This file exists in the directory"
        self.assertEqual(utils.file_exists(filename2), True) "This file does not exist in the directory"

    def test_make_valid_dir(self)
        #Run the utils.make_valid_dir_string function with path which includes os seperator character and
        # one path which does not include an os seperator character
        path1 = 'C:\Users'
        path2 = 'C:\Users\'
        self.assertEqual(utils.make_valid_dir_string(path1), 'C:\Users\')
        self.assertEqual(utils.make_valid_dir_string(path2), 'C\Users\')
    
    def test_make_dir(self)
        #Run the utils.make_dir function with a directory which exists and a directory which does not exist
        #Case1: Directory does not exist. Should make a new directory
        dir = os.path.join('foo', 'bar')
        self.assertFalse(os.path.isdir(dir))
        make_dir(dir)
        self.assertTrue(os.path.isdir(dir))

        #Case2: Directory already exists. Nothing should happen
        self.assertTrue(os.path.isdir(dir))
        make_dir(dir)
        self.assertTrue(os.path.isdir(dir))
        shutil.rmtree(dir)
    
    def test_make_valid_text_file(self)
        #Run the utils.make_valid_text_file function with a file which is a text file and a file which is not a text file
        filename = 'requirements.txt'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_text_file(filename), 'requirements.txt')
        self.assertEqual(utils.make_valid_text_file(filename2), 'testfile.yaml.txt')
    
    def test_make_valid_json_file(self)
        #Run the utils.make_valid_text_file function with a file which is a json file and a file which is not a json file
        filename = 'test.json'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_text_file(filename), 'test.json')
        self.assertEqual(utils.make_valid_text_file(filename2), 'testfile.yaml.json')

    def test_make_valid_csv_file(self)
        #Run the utils.make_valid_text_file function with a file which is a json file and a file which is not a json file
        filename = 'test.csv'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_text_file(filename), 'test.csv')
        self.assertEqual(utils.make_valid_text_file(filename2), 'testfile.yaml.csv')

    def test_make_valid_fasta_file(self)
        #Run the utils.make_valid_text_file function with a file which is a json file and a file which is not a json file
        filename = 'test.fasta'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_text_file(filename), 'test.fasta')
        self.assertEqual(utils.make_valid_text_file(filename2), 'testfile.yaml.fasta')
    
    def test_is_json(self)
        #Run the utils.is_json function with a file which is a json and a file which is not a json file
        filename = 'test.json'
        filename2 = 'test.csv'
        self.assertEqual(utils.is_json(filename), True)
        self.assertEqual(utils.is_json(filename2), False)
    
    def test_is_fasta(self)
        #Run the utils.is_fasta function with a file which is a fasta and a file which is not a fasta file
        filename = 'test.fasta'
        filename2 = 'test.csv'
        self.assertEqual(utils.is_fasta(filename), True)
        self.assertEqual(utils.is_fasta(filename2), False)
    
    def test_is_dir(self) #TODO
        #Run the utils.is_dir function with a path which is a valid path to a directory and a path which isn't
        filename = ''
        filename2 = ''
        self.assertEqual(utils.is_dir(filename), True)
        self.assertEqual(utils.is_dir(filename2), False)

    def test_is_fle(self)
        #Run the utils.is_file function with a file and a file which does not exist
        filename = 'requirements.txt'
        filename2 = 'thisfiledoesnotexist.txt'
        self.assertEqual(utils.is_file(filename), True)
        self.assertEqual(utils.is_file(filename2), False)

    def test_all_perms_of_s(self) 
        #Run the utils.all_perms_of_s function with two strings and varying keywords

        string1 = 'LMWHOMP'
        keyletter1 = 'L,J,I' #expected result: ['LMWHOMP', 'JMWHOMP', 'IMWHOMP']
        string2 = 'MALWAR MZHL'
        keyletter2 = 'L,H' #expected result: ['MALWAR MZHL', 'MAHWAR MZHL', 'MALWAR MZLL', 'MALWAR MZHH', 'MAHWAR MZLL', 'MAHWAR MZHH', 'MAHWAR MZLH']
        self.assertEqual(utils.all_perms_of_s(string1, keyletter1), ['LMWHOMP', 'JMWHOMP', 'IMWHOMP'])
        self.assertEqual(utils.all_perms_of_s(string2, keyletter2), ['MALWAR MZHL', 'MAHWAR MZHL', 'MALWAR MZLL', 'MALWAR MZHH', 'MAHWAR MZLL', 'MAHWAR MZHH', 'MAHWAR MZLH'])
    
    def test_make_sparse_array(self) #TODO
        #Run the utils.make_sparse_array function with two sample spectrums and bin widthss
        spectrum1 = [1.5,3,2.3,5,2.2,35,5,16]
        spectrum2 = [3.5,62,5,6,7.8,1.1,2,4.556,34]

    def test_overlap_intervals
        #Run the utils.overlap_intervals function with two different intervals. One set will overlap and one won't
        interval1 = [0,3]
        interval2 = [2,5] #Expected to return [0,5]
        interval3 = [-1,4]
        interval4 = [6,15] #Expected not to return anything
        self.assertEqual(Utils.overlap_intervals(interval1, interval2), [0,5])
        self.assertEqual(Utils.overlap_intervals(interval3, interval4))
    
    def test_to_percent
        #Run the utils.to_percent function with two different values to convert to a percent.
        value1 = 53/100 #Expected: 53%
        value2 = 234/3456 #Expected: 7%
        self.assertEqual(utils.to_percent(value1, total1), 53)
        self.assertEqual(utils.to_percent(value2, total2), 7)
    
    def test_predicted_len
        
        
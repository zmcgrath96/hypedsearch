import unittest
import os
from src import utils
from src import gen_spectra
from src.objects import Spectrum

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

    def test_file_exists(self):
        #Run the utils.file_exists function with two sample files.
        filename = 'requirements.txt'
        filename2 = 'thisfiledoesnotexist.txt'
        self.assertTrue(utils.file_exists(filename)) # This file exists in the directory
        self.assertFalse(utils.file_exists(filename2)) # This file does not exist in the directory

    def test_make_valid_dir(self):
        #Run the utils.make_valid_dir_string function with path which includes os seperator character and
        # one path which does not include an os seperator character
        path1 = 'C:\\Users'
        path2 = "C:\\Users\\"
        self.assertEqual(utils.make_valid_dir_string(path1), 'C:\\Users\\')
        self.assertEqual(utils.make_valid_dir_string(path2), 'C\\Users\\')
    
    def test_make_dir(self):
        #Run the utils.make_dir function with a directory which exists and a directory which does not exist
        #Case1: Directory does not exist. Should make a new directory
        dir = os.path.join('foo', 'bar')
        self.assertFalse(os.path.isdir(dir))
        utils.make_dir(dir)
        self.assertTrue(os.path.isdir(dir))

        #Case2: Directory already exists. Nothing should happen
        self.assertTrue(os.path.isdir(dir))
        utils.make_dir(dir)
        self.assertTrue(os.path.isdir(dir))
        shutil.rmtree(dir)
    
    def test_make_valid_text_file(self):
        #Run the utils.make_valid_text_file function with a file which is a text file and a file which is not a text file
        filename = 'requirements.txt'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_text_file(filename), 'requirements.txt')
        self.assertEqual(utils.make_valid_text_file(filename2), 'testfile.yaml.txt')
    
    def test_make_valid_json_file(self):
        #Run the utils.make_valid_json_file function with a file which is a json file and a file which is not a json file
        filename = 'test.json'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_json_file(filename), 'test.json')
        self.assertEqual(utils.make_valid_json_file(filename2), 'testfile.yaml.json')

    def test_make_valid_csv_file(self):
        #Run the utils.make_valid_csv_file function with a file which is a json file and a file which is not a json file
        filename = 'test.csv'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_csv_file(filename), 'test.csv')
        self.assertEqual(utils.make_valid_csv_file(filename2), 'testfile.yaml.csv')

    def test_make_valid_fasta_file(self):
        #Run the utils.make_valid_fasta_file function with a file which is a json file and a file which is not a json file
        filename = 'test.fasta'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_fasta_file(filename), 'test.fasta')
        self.assertEqual(utils.make_valid_fasta_file(filename2), 'testfile.yaml.fasta')
    
    def test_is_json(self):
        #Run the utils.is_json function with a file which is a json and a file which is not a json file
        filename = 'test.json'
        filename2 = 'test.csv'
        self.assertEqual(utils.is_json(filename), True)
        self.assertEqual(utils.is_json(filename2), False)
    
    def test_is_fasta(self):
        #Run the utils.is_fasta function with a file which is a fasta and a file which is not a fasta file
        filename = 'test.fasta'
        filename2 = 'test.csv'
        self.assertEqual(utils.is_fasta(filename), True)
        self.assertEqual(utils.is_fasta(filename2), False)
    
    def test_is_dir(self): 
        #Run the utils.is_dir function with a path which is a valid path to a directory and a path which isn't
        filename = os.path.join('foo', 'bar')
        self.assertFalse(utils.is_dir(filename))
        utils.make_dir(filename)
        self.assertTrue(utils.is_dir(filename))

    def test_is_fle(self):
        #Run the utils.is_file function with a file and a file which does not exist
        filename = 'requirements.txt'
        filename2 = 'thisfiledoesnotexist.txt'
        self.assertEqual(utils.is_file(filename), True)
        self.assertEqual(utils.is_file(filename2), False)

    def test_all_perms_of_s(self): 
        #Run the utils.all_perms_of_s function with two strings and varying keywords

        string1 = 'LMWHOMP'
        keyletter1 = 'LJI' #expected result: ['LMWHOMP', 'JMWHOMP', 'IMWHOMP']
        string2 = 'MALWAR MZHL'
        keyletter2 = 'LH' #expected result: ['MALWAR MZHL', 'MAHWAR MZHL', 'MALWAR MZLL', 'MALWAR MZHH', 'MAHWAR MZLL', 
        #'MAHWAR MZHH', 'MAHWAR MZLH', 'MALWAR MZLH']
        self.assertEqual(sorted(utils.all_perms_of_s(string1, keyletter1)), sorted(['LMWHOMP', 'JMWHOMP', 'IMWHOMP']))
        self.assertEqual(sorted(utils.all_perms_of_s(string2, keyletter2)), sorted(['MALWAR MZHL', 'MAHWAR MZHL', 
            'MALWAR MZLL', 'MALWAR MZHH', 'MAHWAR MZLL', 'MAHWAR MZHH', 'MAHWAR MZLH', 'MALWAR MZLH']))
    
    def test_make_sparse_array(self): #TODO
        #Run the utils.make_sparse_array function with two sample spectrums and bin widths
        spectrum1 = [1.5, 3, 2.3, 5, 2.2, 35, 5, 16]
        spectrum2 = [3.5, 62, 5, 6, 7.8, 1.1, 2, 4.556, 34]

    def test_overlap_intervals(self):
        #Run the utils.overlap_intervals function with two different intervals. One set will overlap and one won't
        intervals1 = [[0,3], [2,5], [3,7]] #Expected to return [0,7]
        intervals2 = [[-1,4], [6,15]] #Expected not to return anything
        self.assertEqual(utils.overlap_intervals(intervals1), [0,5])
        self.assertEqual(utils.overlap_intervals(intervals2))
    
    def test_to_percent(self):
        #Run the utils.to_percent function with two different values to convert to a percent.
        value1 = 53
        total1 = 100 #Expected: 53%
        value2 = 234
        total2 = 3456 #Expected: 7%
        self.assertEqual(utils.to_percent(value1, total1), 53)
        self.assertEqual(utils.to_percent(value2, total2), 7)
    
    def test_predicted_len(self):
        #Run the utils.predicted_len function with two different precursor masses.
        precursor = 240 #Expected len would be 5
        charge = 1
        self.assertEqual(utils.predicted_len(precursor, charge), 5)
        precursor = precursor * 5
        charge = charge + 1 #Expected len would be 11
        self.assertEqual(utils.predicted_len(precursor, charge), 11)
    
    def test_predicted_len_precursor(self): 
        sequence = 'MALWAR'
        spectrum = Spectrum(gen_spectra.gen_spectrum(sequence), gen_spectra.get_precursor(sequence))
        expected_length = (len(sequence) * (spectrum.precursor_mass / 728.379201))
        self.assertEqual(utils.predicted_len_precursor(spectrum, sequence), expected_length)
    
    def test_hashable_boundaries(self):
        #run the hashable_boundaries function with two lists of boundaries
        boundaries2 = [2,10] #Expected to return 2-10
        boundaries3 = [3,5,7] #Expected to return nothing
        self.assertEqual(utils.hashable_boundaries(boundaries2), '2-10')
        self.assertEqual(utils.hashable_boundaries(boundaries3), None)

    def test_split_hybrid(self):  
        #run the split_hybrid function with two samples sequence
        Sequence1 = 'MAL-WAR' #Expected left is "MAL" and expected right is "WAR"
        Sequence2 = 'DLTQTL-B' #Expected left is "DLTQTL" and expected right is "B"
        self.assertEqual(utils.split_hybrid(Sequence1), ('MAL', 'WAR'))
        self.assertEqual(utils.split_hybrid(Sequence2), ('DLTQTL', 'B'))


import unittest
from src.alignment import mass_comparisons, kmer_extension
from src.spectra import gen_spectra
class test_kmer_extension(unittest.TestCase):
    def setUp(self):
        self.pepsequence = 'MALWARM'
        self.protsequence = 'TWSKMALWARMQVCE'
        self.pepspectra = gen_spectra.gen_spectrum(self.pepsequence)['spectrum']

    def test_new_entry(self):
        old_entry_valid_b = {
            'k': 3, 
            'sequence': 'MAL',
            'starting_position': 4,
            'ending_position': 6,
            'b_score': 2, 
            'y_score': 0
        }
        old_entry_invalid_b = {
            'k': None, 
            'sequence': None,
            'starting_position': 0,
            'ending_position': 4,
            'b_score': None, 
            'y_score': None
        }
        old_entry_valid_y = {
            'k': 3, 
            'sequence': 'ARM',
            'starting_position': 8,
            'ending_position': 10,
            'b_score': 0,
            'y_score': 2
        }
        old_entry_invalid_y = {
            'k': None,
            'sequence': None,
            'starting_position': 10,
            'ending_position': 14,
            'b_score': None,
            'y_score': None
        }

        new_b_entry = kmer_extension.new_entry(old_entry_valid_b, self.protsequence, self.pepspectra)
        new_y_entry = kmer_extension.new_entry(old_entry_valid_y, self.protsequence, self.pepspectra, 'y')

        self.assertEqual(4, new_b_entry['k'], 'Valid old b entry should be extended to k=4')
        self.assertEqual('MALW', new_b_entry['sequence'], 'Valid old b entry should extend sequnce to MALW')
        self.assertEqual(old_entry_invalid_b, kmer_extension.new_entry(old_entry_invalid_b, self.protsequence, self.pepsequence), 'Invalid old b entry should return the same entry')

        self.assertEqual(4, new_y_entry['k'], 'Valid old y entry should be extended to k=4')
        self.assertEqual('WARM', new_y_entry['sequence'], 'Valid old y entry should be extended to WARM')
        self.assertEqual(old_entry_invalid_y, kmer_extension.new_entry(old_entry_invalid_y, self.protsequence, self.pepspectra, 'y'), 'Invalid old y entry should return same entry')

    def test_kmer_extend(self): 
        score_b = lambda seq, spec: mass_comparisons.compare_masses(spec, gen_spectra.gen_spectrum(seq, ion='b')['spectrum'])
        score_y = lambda seq, spec: mass_comparisons.compare_masses(spec, gen_spectra.gen_spectrum(seq, ion='y')['spectrum'])
        kmers = [
            {
                'k': 3, 
                'starting_position': i, 
                'ending_position': i+2, 
                'sequence': self.protsequence[i: i+3], 
                'b_score': score_b(self.protsequence[i: i+3], self.pepspectra),
                'y_score': score_y(self.protsequence[i: i+3], self.pepspectra)
            } for i in range(len(self.pepspectra) - 2)
        ]

        b_list, y_list = kmer_extension.kmer_extend(self.pepspectra, self.protsequence, kmers, kmers)
        self.assertEqual('MALWARM', b_list[0]['sequence'], 'For ideal spectra, the best result should be the full peptide sequence')
        self.assertEqual('MALWARM', y_list[0]['sequence'], 'For ideal spectra, the best result shoudl be the full peptide sequence')



if __name__ == '__main__':
    unittest.main()
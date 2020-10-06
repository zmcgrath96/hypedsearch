import unittest
from src.utils import predicted_len
from src.cppModules import map_spectra_masses

from collections import namedtuple

BasicProtein = namedtuple('BasicProtein', ['name', 'sequence'], defaults=['', ''])

class test_map_masses(unittest.TestCase):
    def setUp(self):

        # create a set of proteins
        self.basic_prots = [
            BasicProtein('Protein 1', 'MALWARM'), 
            BasicProtein('Protein 2', 'GGGQQUVRS'), 
            BasicProtein('Protein 3', 'NRGQWEVE')
        ]

        # create a set of boundaries
        self.boundaries = [
            [58.0276, 58.0299],
            [106.0477, 106.052],
            [115.0479, 115.0525],
            [132.0451, 132.0504],
            [148.0575, 148.0634],
            [150.0553, 150.0613],
            [172.0682, 172.0751],
            [203.0808, 203.0889],
            [247.1239, 247.1338],
            [262.1457, 262.1562],
            [271.1459, 271.1567],
            [300.1242, 300.1362],
            [306.1533, 306.1656],
            [316.1626, 316.1753],
            [328.1662, 328.1793],
            [361.2122, 361.2266],
            [376.1639, 376.179],
            [377.189, 377.2041],
            [428.1803, 428.1974],
            [456.2222, 456.2405],
            [502.2382, 502.2583],
            [512.1628, 512.1833],
            [562.2395, 562.262],
            [563.2646, 563.2871],
            [573.2739, 573.2968],
            [579.1309, 579.154],
            [640.2188, 640.2444],
            [642.2978, 642.3235],
            [676.3464, 676.3735],
            [678.1973, 678.2244],
            [690.2955, 690.3231],
            [729.3719, 729.4011],
            [747.3158, 747.412],
            [768.2748, 768.3055],
            [771.3378, 771.3687],
            [825.2951, 825.3282],
            [834.2953, 834.3287],
            [860.4098, 860.4442],
            [870.4043, 870.4391],
            [878.42, 878.4551],
            [882.3155, 882.3508],
            [903.4138, 903.45],
            [921.3256, 921.3624],
            [939.3358, 939.3734],
            [999.4443, 999.4843],
            [1017.4545, 1017.4952]
        ]

        self.correct_mappings = {
            '58.0276-58.0299': ['G', 'GG', 'N'],
            '106.0477-106.052': ['S'],
            '115.0479-115.0525': ['GG', 'N'],
            '132.0451-132.0504': ['M'],
            '148.0575-148.0634': ['E'],
            '150.0553-150.0613': ['M'],
            '172.0682-172.0751': ['GGG'],
            '203.0808-203.0889': ['MA'],
            '247.1239-247.1338': ['VE', 'EV'],
            '262.1457-262.1562': ['RS'],
            '271.1459-271.1567': ['NR'],
            '300.1242-300.1362': ['GGGQ'],
            '306.1533-306.1656': ['RM'],
            '316.1626-316.1753': ['MAL'],
            '328.1662-328.1793': ['NRG'],
            '361.2122-361.2266': ['VRS'],
            '376.1639-376.179': ['EVE'],
            '377.189-377.2041': ['ARM'],
            '428.1803-428.1974': ['GGGQQ'],
            '456.2222-456.2405': ['NRGQ'],
            '502.2382-502.2583': ['MALW'],
            '512.1628-512.1833': ['UVRS'],
            '562.2395-562.262': ['WEVE'],
            '563.2646-563.2871': ['WARM'],
            '573.2739-573.2968': ['MALWA'],
            '579.1309-579.154': ['GGGQQU'],
            '640.2188-640.2444': ['QUVRS'],
            '642.2978-642.3235': ['NRGQW'],
            '676.3464-676.3735': ['LWARM'],
            '678.1973-678.2244': ['GGGQQUV'],
            '690.2955-690.3231': ['QWEVE'],
            '729.3719-729.4011': ['MALWAR', 'ALWARM'],
            '747.3158-747.412': ['GQWEVE', 'ALWARM', 'MALWAR'],
            '768.2748-768.3055': ['QQUVRS'],
            '771.3378-771.3687': ['NRGQWE'],
            '825.2951-825.3282': ['GQQUVRS'],
            '834.2953-834.3287': ['GGGQQUVR'],
            '860.4098-860.4442': ['MALWARM'],
            '870.4043-870.4391': ['NRGQWEV'],
            '878.42-878.4551': ['MALWARM'],
            '882.3155-882.3508': ['GGQQUVRS'],
            '903.4138-903.45': ['RGQWEVE'],
            '921.3256-921.3624': ['GGGQQUVRS'],
            '939.3358-939.3734': ['GGGQQUVRS'],
            '999.4443-999.4843': ['NRGQWEVE'],
            '1017.4545-1017.4952': ['NRGQWEVE']
        }

        self.correct_kmer_mappings = {
            'M': ['Protein 1'],
            'A': ['Protein 1'],
            'L': ['Protein 1'],
            'W': ['Protein 3', 'Protein 1'],
            'R': ['Protein 2', 'Protein 3', 'Protein 1'],
            'MA': ['Protein 1'],
            'AL': ['Protein 1'],
            'LW': ['Protein 1'],
            'WA': ['Protein 1'],
            'AR': ['Protein 1'],
            'RM': ['Protein 1'],
            'MAL': ['Protein 1'],
            'ALW': ['Protein 1'],
            'LWA': ['Protein 1'],
            'WAR': ['Protein 1'],
            'ARM': ['Protein 1'],
            'MALW': ['Protein 1'],
            'ALWA': ['Protein 1'],
            'LWAR': ['Protein 1'],
            'WARM': ['Protein 1'],
            'MALWA': ['Protein 1'],
            'ALWAR': ['Protein 1'],
            'LWARM': ['Protein 1'],
            'MALWAR': ['Protein 1'],
            'ALWARM': ['Protein 1'],
            'MALWARM': ['Protein 1'],
            'G': ['Protein 2', 'Protein 3'],
            'Q': ['Protein 2', 'Protein 3'],
            'U': ['Protein 2'],
            'V': ['Protein 2', 'Protein 3'],
            'S': ['Protein 2'],
            'GG': ['Protein 2'],
            'GQ': ['Protein 2', 'Protein 3'],
            'QQ': ['Protein 2'],
            'QU': ['Protein 2'],
            'UV': ['Protein 2'],
            'VR': ['Protein 2'],
            'RS': ['Protein 2'],
            'GGG': ['Protein 2'],
            'GGQ': ['Protein 2'],
            'GQQ': ['Protein 2'],
            'QQU': ['Protein 2'],
            'QUV': ['Protein 2'],
            'UVR': ['Protein 2'],
            'VRS': ['Protein 2'],
            'GGGQ': ['Protein 2'],
            'GGQQ': ['Protein 2'],
            'GQQU': ['Protein 2'],
            'QQUV': ['Protein 2'],
            'QUVR': ['Protein 2'],
            'UVRS': ['Protein 2'],
            'GGGQQ': ['Protein 2'],
            'GGQQU': ['Protein 2'],
            'GQQUV': ['Protein 2'],
            'QQUVR': ['Protein 2'],
            'QUVRS': ['Protein 2'],
            'GGGQQU': ['Protein 2'],
            'GGQQUV': ['Protein 2'],
            'GQQUVR': ['Protein 2'],
            'QQUVRS': ['Protein 2'],
            'GGGQQUV': ['Protein 2'],
            'GGQQUVR': ['Protein 2'],
            'GQQUVRS': ['Protein 2'],
            'GGGQQUVR': ['Protein 2'],
            'GGQQUVRS': ['Protein 2'],
            'GGGQQUVRS': ['Protein 2'],
            'N': ['Protein 3'],
            'E': ['Protein 3'],
            'NR': ['Protein 3'],
            'RG': ['Protein 3'],
            'QW': ['Protein 3'],
            'WE': ['Protein 3'],
            'EV': ['Protein 3'],
            'VE': ['Protein 3'],
            'NRG': ['Protein 3'],
            'RGQ': ['Protein 3'],
            'GQW': ['Protein 3'],
            'QWE': ['Protein 3'],
            'WEV': ['Protein 3'],
            'EVE': ['Protein 3'],
            'NRGQ': ['Protein 3'],
            'RGQW': ['Protein 3'],
            'GQWE': ['Protein 3'],
            'QWEV': ['Protein 3'],
            'WEVE': ['Protein 3'],
            'NRGQW': ['Protein 3'],
            'RGQWE': ['Protein 3'],
            'GQWEV': ['Protein 3'],
            'QWEVE': ['Protein 3'],
            'NRGQWE': ['Protein 3'],
            'RGQWEV': ['Protein 3'],
            'GQWEVE': ['Protein 3'],
            'NRGQWEV': ['Protein 3'],
            'RGQWEVE': ['Protein 3'],
            'NRGQWEVE': ['Protein 3']
        }

        # calculate the max length
        self.est_length = predicted_len(self.boundaries[-1][1])

    def test_map_spectra_masses(self):
        matched_b_masses, matched_y_masses, kmer_to_prots = map_spectra_masses.map_boundaries(self.boundaries, self.basic_prots, self.est_length)

        # see if all of the keys are in correct mappings
        for k, v in matched_b_masses.items():
            self.assertIn(k, self.correct_mappings.keys(), 'Each key should be in the correct mappings')
            self.assertEqual(sorted(v), sorted(self.correct_mappings[k]), 'mappings should have the same kmers')

        # see if all of the keys are in correct mappings
        for k, v in matched_y_masses.items():
            self.assertIn(k, self.correct_mappings.keys(), 'Each key should be in the correct mappings')
            self.assertEqual(sorted(v), sorted(self.correct_mappings[k]), 'mappings should have the same kmers')

        # check that all the kmer mappings are present
        for k, v in kmer_to_prots.items():
            self.assertIn(k, self.correct_kmer_mappings.keys(), 'Every kmer should be present in the correct mappings')
            self.assertEqual(sorted(v), sorted(self.correct_kmer_mappings[k]), 'Kmers should map to the correct source proteins')

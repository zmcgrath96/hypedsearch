import unittest
from src.identfication import search
from src.types.objects import Spectrum, MassSequence
from src.spectra.gen_spectra import gen_spectrum
from collections import defaultdict
import math

class test_search(unittest.TestCase):
    def setUp(self):
        pass

    def test_search_kmers_hash(self):
        seq = 'MALWAR'
        spec = gen_spectrum(seq)

        s = Spectrum(spec['spectrum'], [100 for _ in range(len(spec['spectrum']))], 2, 0, spec['precursor_mass'], '')
        d = defaultdict(list)
        for i in range(1, len(seq) + 1):
            bseq = seq[:i]
            yseq = seq[len(seq)-i:]

            bmass = max(gen_spectrum(bseq, ion='b', charge=1)['spectrum'])
            ymass = max(gen_spectrum(yseq, ion='y', charge=1)['spectrum'])

            d[math.floor(bmass)].append(MassSequence(bmass, bseq))
            d[math.floor(ymass)].append(MassSequence(ymass, yseq))

        hits = search.search_kmers_hash(s, d, 10)

        self.assertEqual(len(hits), 12, 'for singly charged ions, 2xlen(sequence) is the expected hit count')
            

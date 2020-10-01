from src import bits 
from pyteomics import fasta 
from src.params import DATABASE_FILE
from src.sequence import gen_spectra

max_kmer = 30

kmers = []
all_spec = []
for i, entry in enumerate(fasta.read(DATABASE_FILE)):
    print(f'\rOn protein {i+1}/~17000', end='')

    seq = entry.sequence

    for j in range(min(max_kmer, len(seq))):
        # kmers.append(bits.str_to_bits(seq[:j]))
        # kmers.append(seq[:j])
        all_spec.append(gen_spectra.gen_spectrum(seq[:j])['spectrum'])

    for j in range(min(max_kmer, len(seq))):
        # kmers.append(bits.str_to_bits(seq[-j-1:]))
        # kmers.append(seq[-j-1:])
        all_spec.append(gen_spectra.gen_spectrum(seq[-j-1:])['spectrum'])

    if max_kmer > len(seq):
        continue

    for j in range(1, len(seq) - max_kmer):
        # kmers.append(bits.str_to_bits(seq[j:j+max_kmer]))
        # kmers.append(seq[j:j+max_kmer])
        all_spec.append(gen_spectra.gen_spectrum(seq[j:j+max_kmer])['spectrum'])
        
print(f'\nNumber of elements in the list: {len(kmers)}')
print()
# spin so that I can see the memory
# for i in range(len(kmers)):
#     print(f'\r{kmers[i]}/1000000000', end='')
# print('\nDone')
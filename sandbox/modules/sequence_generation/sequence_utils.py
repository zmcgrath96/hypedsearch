import numpy as np

'''length_dist

DESC:
    return a list of values that are the lengths peptides should be
Inputs:
    num_peptides: int number of lengths to generate
kwargs:
    dist: str the type of distribution to use. Default=beta (based on experimental data)
            other values: random
    min_length: int minium length peptide to generate. Default=4
    max_length: int maximum length peptide to generate. Default=35
    a: float alpha value. Default=2.3 (based on experimental data)
    b: float beta value. Default=7 (based on experimental data)
Outputs:
    list of int lengths drawn from distribution
'''
def length_dist(num_peptides, dist='beta', min_length=4, max_length=35, a=2.3, b=7):
    if dist.lower() == 'beta':
        y = max_length
        x = min_length
        return x + (np.random.beta(a, b, size=(num_peptides)) * (y - x))
    else: # assume random
        return np.random.random_integers(min_length, high=max_length, size=num_peptides)


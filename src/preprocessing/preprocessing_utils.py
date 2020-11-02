from src.file_io import mzML
from src.utils import ppm_to_da, overlap_intervals

def load_spectra(spectra_files: list, ppm_tol: int, peak_filter=0, relative_abundance_filter=0) -> (list, list, dict):
    '''
    Load all the spectra files into memory and merge all 
    spectra into one massive list for reduction of the search space

    indicesnputs:
        spectra_files:  (list) full paths to all spectra files
    kwargs:
        peak_filter:    (int) the top X most abundant spectra to keep. Default=0
        relative_abundance_filter:  (float) the percentage (0, 1) of the total abundance 
                                    a peak must have to be considered significant. Default=0
    Outputs:
        (list, dict, list) all spectra in memory of boundaries namedtuples, a list of floats of all masses
    '''
    # the single dimension of all the rounded spectra to use for reducing the search space
    linear_spectra = []

    # a list of all boundaries namedtuples
    all_spectra = []

    # go through each spectra file and load them into memory
    for spectra_file in spectra_files:
        these_spectra = mzML.read(
            spectra_file, 
            peak_filter=peak_filter, 
            relative_abundance_filter=relative_abundance_filter
        )

        all_spectra += these_spectra

        # go through each mass of each s, load it into memory, 
        # round the numbers to 3 decimal places for easier, and append to linear_spectra
        linear_spectra += list(set([
            x for spectrum in these_spectra for x in spectrum.spectrum
        ]))

    # sort the linear spectra
    linear_spectra.sort()

    # turn the all spectra list into a list of boundaries
    def make_boundaries(mz):
        da_tol = ppm_to_da(mz, ppm_tol)
        return [mz - da_tol, mz + da_tol]

    boundaries = [make_boundaries(mz) for mz in linear_spectra]

    # make overlapped boundaries larger boundaries
    boundaries = overlap_intervals(boundaries)

    # make a mapping for mz -> boundaries
    b_i, s_i = 0, 0
    mz_mapping = {}
    while s_i < len(linear_spectra):
        
        # if the current boundary encapsulates s_i, add to list
        if boundaries[b_i][0] <= linear_spectra[s_i] <= boundaries[b_i][1]:
            mz_mapping[linear_spectra[s_i]] = b_i 
            s_i += 1

        # else if the s_i < boundary, increment s_i
        elif linear_spectra[s_i] < boundaries[b_i][0]:
            s_i += 1

        # else if s_i > boundary, incrment b_i
        elif linear_spectra[s_i] > boundaries[b_i][1]:
            b_i += 1

    return (all_spectra, boundaries, mz_mapping)
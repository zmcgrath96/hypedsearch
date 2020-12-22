def relative_abundance_filtering(masses: list, abundances: list, percentage: float) -> (list, list):
    '''
    Take all peaks from the spectrum who's abundance is at least percentage of the total abundances. 
    It is assumed that the masses and abundances lists share ordering

    Inputs:
        masses:     (list) floats of masses
        abundances: (list) floats of the abundances of the corresponding masses
        percentage: (float) the minimum percentage of abundance a peak must have. Must be in range (0, 1)
    Outputs:
        (list, list) (masses, abundances) pair (pairing is the same as the input) that
                    pass the filter
    '''
    # total intensity
    ti = sum(abundances)

    # find the filter value
    min_value = ti * percentage

    # zip them up and take values that pass the filter
    filtered_mass_abundances = [x for x in zip(masses, abundances) if x[1] >= min_value]

    # split them off and return
    masses = [float(x) for x, _ in filtered_mass_abundances]
    abundances = [float(x) for _, x in filtered_mass_abundances]

    return (masses, abundances)


def peak_filtering(masses: list, abundances: list, num_peaks: int) -> (list, list):
    '''
    Take the most abundant peaks and return the sorted masses with the abundances.
    It is assumed that the masses and abundances lists share ordering

    Inputs:
        masses:     (list) floats of masses
        abundances: (list) floats of the abundances of the corresponding masses
        num_peaks:  (int) the number of peaks to retain
    Outputs:
        (list, list) (masses, abundances) pair (pairing is the same as the input) that
                    pass the filter
    '''
    # zip the abundance and the m/z values together
    mass_abundances = zip(masses, abundances)
    
    # sort by key 1, the abundance, and take the top peak filter results
    mass_abundances = sorted(mass_abundances, key=lambda x: x[1], reverse=True)[:num_peaks]

    # sort them now by the value of m/z
    mass_abundances.sort(key=lambda x: x[0])

    # seperate them
    masses = [float(x) for x, _ in mass_abundances]
    abundances = [float(x) for _, x in mass_abundances]

    return (masses, abundances)
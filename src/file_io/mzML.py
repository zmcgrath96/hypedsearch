from src.utils import file_exists
from src.types.objects import Spectrum
from pyteomics import mzml

def read(filename: str, peak_filter=50) -> list:
    '''
    read an .mzML file into memory

    Inputs:
        filename:       (str) path to the file to import
    kwargs:
        peak_filter:    (int) the top number of peaks to keep. Default=50
    Outputs:
        (list) Spectrum namedtuple instances
    '''
    if not file_exists(filename):
        print('File {} not found. Please make sure that this file exists'.format(filename))
        return
    else: 
        spectra = []
        
        filecontents = mzml.read(filename)
        for content in filecontents:

            # zip the abundance and the m/z values together
            mass_abundances = zip(list(content['m/z array']), list(content['intensity array']))
            
            # sort by key 1, the abundance, and take the top peak filter results
            mass_abundances = sorted(mass_abundances, key=lambda x: x[1], reverse=True)[:peak_filter]

            # sort them now by the value of m/z
            mass_abundances.sort(key=lambda x: x[0])

            # seperate them
            masses = [float(x) for x, _ in mass_abundances]
            abundances = [float(x) for _, x in mass_abundances]

            # get the precursor
            precursor = None
            precursor_list = content['precursorList']['precursor']
            for p in precursor_list:
                for selected_ion in p['selectedIonList']['selectedIon']:
                    if int(selected_ion['charge state']) == 2:
                        precursor = int(selected_ion['charge state'])

            precursor = precursor if precursor is not None else max(masses)/2
          
            spectra.append(Spectrum(
                masses,
                abundances,
                int(content['ms level']),
                int(content['index']),
                precursor,
                filename
            ))

        return spectra
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
            

            # sort the spectra 
            sorted_labeled_spec = sorted([(i, float(spec)) for i, spec in list(enumerate(content['m/z array']))], key=lambda x: x[1])[:peak_filter]
            sorted_keys = [x[0] for x in sorted_labeled_spec]
            sorted_spec = [x[1] for x in sorted_labeled_spec]
            sorted_abundances = [float(content['intensity array'][i]) for i in sorted_keys]

            spectra.append(Spectrum(
                sorted_spec,
                sorted_abundances,
                int(content['ms level']),
                int(content['index']),
                int(max(content['m/z array'])),
                filename
            ))

        return spectra
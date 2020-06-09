from src.utils import file_exists
from src.types.objects import Spectrum
from pyteomics import mzml

def read(filename: str) -> list:
    '''
    read an .mzML file into memory

    Inputs:
        filename:   (str) path to the file to import
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
            spectra.append(Spectrum(
                [float(x) for x in content['m/z array']],
                [float(x) for x in content['intensity array']],
                int(content['ms level']),
                int(content['index']),
                int(max(content['m/z array'])),
                filename
            ))

        return spectra
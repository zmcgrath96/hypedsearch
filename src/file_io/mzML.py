from pyopenms import MSExperiment, MzMLFile
from utils.utils import file_exists

def read(file):
    '''
    read an .mzML file into memory

    Inputs:
        file: str path to the file to import
    Outputs:
        NONE if file is not found
    '''
    if not file_exists(file):
        print('File {} not found. Please make sure that this file exists'.format(file))
        return
    else: 
        spectra = []
        exp = MSExperiment()
        MzMLFile().load(file, exp)

        for s in exp.getSpectra():
            spectra.append({
                'level': s.getMSLevel(),
                'spectrum': list(list(s.get_peaks())[0]), 
                'scan_no': int(str(s.getNativeID()).split('=')[-1].replace("'", ''))
            })

        return spectra
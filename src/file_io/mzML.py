from pyopenms import MSExperiment, MzMLFile
from src.utils.utils import file_exists

def read(file: str) -> list:
    '''
    read an .mzML file into memory

    Inputs:
        file: str path to the file to import
    Outputs:
        list of dictionaries with spectra information of the form
        {
            'level': int, 
            'spectrum': list of floats,
            'abundance': list of floats,
            'scan_no': int
        }
    '''
    if not file_exists(file):
        print('File {} not found. Please make sure that this file exists'.format(file))
        return
    else: 
        spectra = []
        exp = MSExperiment()
        MzMLFile().load(file, exp)

        for s in exp.getSpectra():
            # get_peaks returns two lists. 0 is the mz and the second is the abundance
            spec = [float(x) for x in list(s.get_peaks())[0]]
            abundance = [float(x) for x in list(s.get_peaks())[1]]
            if not len(spec) > 0:
                continue
            precursor = max(spec)
            spectra.append({
                'level': int(s.getMSLevel()),
                'spectrum': spec, 
                'abundance': abundance, 
                'scan_no': int(str(s.getNativeID()).split('=')[-1].replace("'", '')),
                'precursor_mass': precursor
            })

        return spectra
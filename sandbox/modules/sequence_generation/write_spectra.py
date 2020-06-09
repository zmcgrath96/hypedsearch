from pyopenms import MSExperiment, MSSpectrum, MzMLFile, Peak1D, Precursor
import xml.etree.ElementTree as ET

import struct
import base64

ET.register_namespace('', "http://psi.hupo.org/ms/mzml")
'''write_mzml

DESC:
    create a mass spectrum file in mzml of sequences
Inputs:
    file_name: str name to save the file in
    spectra: list of dictionaries of the form [{spectrum: list[floats], precursor_mass: float}]
             This data is written to file (the spectrum and the precursor)
kwargs:
    title_prefix: str name to give as prefix to the name of each spectrum. Default=Spectrum
    output_dir: str name of the directory to save files to. Default=./
Outputs:
    list of strings of file paths
'''
def write_mzml(file_name, spectra, title_prefix='Spectrum ', output_dir='./'):
    if '.mzml' not in file_name.lower():
        file_name += '.mzML'

    output_file = output_dir + file_name

    exp = MSExperiment()
    sp_count = 0

    tic = []
    for spectrum in spectra:
        spec = MSSpectrum()
        spec.setMSLevel(2)
        name = str.encode(title_prefix + str(sp_count))
        spec.setName(name)
        sp_count += 1
        if 'abundance' not in spectrum:
            print('oops')
            i = [500 for _ in spectrum['spectrum']]
        else:
            i = spectrum['abundance']
        tic.append(sum(i))
        spec.set_peaks([spectrum['spectrum'], i])
        spec.setMSLevel(2)
        prec = Precursor()
        prec.setCharge(2)
        prec.setMZ(spectrum['precursor_mass'])
        spec.setPrecursors([prec])
        spec.sortByPosition()
        exp.addSpectrum(spec)

    MzMLFile().store(output_file, exp)

    # load it back and fix it
    tree = ET.parse(output_file)
    mzml = tree.find('{http://psi.hupo.org/ms/mzml}mzML')
    # add an id to the mzml
    mzml.set('id', 'testSpectraFileFixed')
    mzml.set('xmlns', "http://psi.hupo.org/ms/mzml")
    mzml.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    # remove some info from the mzml header [accession, version]
    toDelMzml = ['accession']
    for todel in toDelMzml:
        del mzml.attrib[todel]
        
    # add file description stuff to the fileDescription tag
    fileDescription = mzml.find('{http://psi.hupo.org/ms/mzml}fileDescription')


    sourceFileElement = ET.fromstring('<sourceFileList count="1"> <sourceFile id="NOD2_E3.mzXML" name="NOD2_E3.mzXML" location="file:///C:/Users/zachmcgrath/Downloads"> <cvParam cvRef="MS" accession="MS:1000776" name="scan number only nativeID format" value=""/> <cvParam cvRef="MS" accession="MS:1000566" name="ISB mzXML format" value=""/> <cvParam cvRef="MS" accession="MS:1000569" name="SHA-1" value="89ecb0dd31ca3a2fdf5ef2c4f5341f6e5e9f06f0"/> </sourceFile> </sourceFileList>')
    fileDescription.append(sourceFileElement)

    instrumentConfigurationList = mzml.find('{http://psi.hupo.org/ms/mzml}instrumentConfigurationList')

    instrumentConfiguration = ET.fromstring(' <instrumentConfiguration id="IC1"> <componentList count="3"> <source order="1"> <userParam name="msIonisation" value="HPLC-Chip/MS"/>  </source> <analyzer order="1"> <userParam name="msMassAnalyzer" value="Q-TOF"/>  </analyzer> <detector order="1"> <userParam name="msDetector" value="ADC"/> </detector> </componentList> </instrumentConfiguration>')
    for child in instrumentConfigurationList:
        del child

    instrumentConfigurationList.append(instrumentConfiguration)

    run = mzml.find('{http://psi.hupo.org/ms/mzml}run')
    spectrumList = run.find('{http://psi.hupo.org/ms/mzml}spectrumList')
    print(len(spectrumList))
    for spectrumElement in spectrumList:
        # need to add this
        # <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
        centroidElement = ET.Element('{http://psi.hupo.org/ms/mzml}cvParam', attrib={'name': 'centroid spectrum', 'accession':'MS:1000127', 'value': ''})
        pl = spectrumElement.find('{http://psi.hupo.org/ms/mzml}precursorList')
        bdal = spectrumElement.find('{http://psi.hupo.org/ms/mzml}binaryDataArrayList')
        if bdal is None or pl is None:
            print(spectrumElement)
            spectrumList.remove(spectrumElement)
            
        # add the range. Not good for MSConvert for whatever reason
        # scanListEl = spectrumElement.find('{http://psi.hupo.org/ms/mzml}scanList')
        # scanEl = scanListEl.find('{http://psi.hupo.org/ms/mzml}scan')[0]
        # scanWindowListElement = ET.fromstring('<scanWindowList count="1"> <scanWindow> <cvParam cvRef="MS" accession="MS:1000501" value="0" name="scan window lower limit" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" /> <cvParam cvRef="MS" accession="MS:1000500" value="10000" name="scan window upper limit" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" /> </scanWindow></scanWindowList>')
        # scanEl.append(scanWindowListElement)
        
        spectrumElement.append(centroidElement)

    timearray = base64.encodebytes(struct.pack(f'<{len(tic)}f', *[i+1 for i in range(len(tic))]))
    ticarray = base64.encodebytes(struct.pack(f'<{len(tic)}f', *tic))


    chromatogramListElmentString = f' <chromatogramList count="1" defaultDataProcessingRef="mzLibProcessing"> <chromatogram id="TIC" index="0" defaultArrayLength="1025" dataProcessingRef="mzLibProcessing"> <cvParam cvRef="MS" accession="MS:1000235" value="" name="total ion current chromatogram" /> <binaryDataArrayList count="2"> <binaryDataArray encodedLength="10936"> <cvParam cvRef="MS" accession="MS:1000523" value="" name="64-bit float" /> <cvParam cvRef="MS" accession="MS:1000576" value="" name="no compression" /> <cvParam cvRef="MS" accession="MS:1000595" value="" name="time array" unitAccession="UO:0000031" unitName="Minutes" unitCvRef="UO" /> <binary>{timearray}</binary> </binaryDataArray> <binaryDataArray encodedLength="10936"> <cvParam cvRef="MS" accession="MS:1000523" value="" name="64-bit float" /> <cvParam cvRef="MS" accession="MS:1000576" value="" name="no compression" /> <cvParam cvRef="MS" accession="MS:1000515" value="" name="intensity array" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS" /> <binary>{ticarray}</binary> </binaryDataArray> </binaryDataArrayList> </chromatogram> </chromatogramList>'
    chromatogramListElment = ET.fromstring(chromatogramListElmentString)
    run.append(chromatogramListElment)

    tree.write(output_file, encoding='utf-8', xml_declaration=True)

    return output_file
    

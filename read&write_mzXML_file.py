#This code reads a mzXML file, transforms it to a mzML file and writes on it using pyopensms

import numpy as np
import matplotlib as plt
from pyopenms import MSExperiment, MzXMLFile

def correct_oscilations (spectrum):
    """
    Corrects the oscilations of the MS 
    Args:
        spectrum()
        
        CREATE THE CODE FOR THIS FUNCTION
        1st-> extract all data from file
            decompress zlib
            decode base 64 
        2nd-> apply filtering of that data
        
    """
def plot_spectra(original_spectrum, corrected_spectrum):
    """
    Plots the original data obtained from the mzXML file and the new representation
    of the data with corrected oscillations on the same graph with different colors.
    
    Args:
        original_spectra: List of spectra before correction.
        corrected_spectra: List of spectra after correction.
    """
    plt.figure(figsize=(12, 6))
    
    for orig, corr in zip(original_spectrum, corrected_spectrum):
        mz_orig, intensity_orig = orig.get_peaks()
        mz_corr, intensity_corr = corr.get_peaks()
        plt.plot(mz_orig, intensity_orig, label='Original Spectrum', color='blue', alpha=0.6)
        plt.plot(mz_corr, intensity_corr, label='Corrected Spectrum', color='red', alpha=0.8)
    
    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.legend()
    plt.title('Original vs Corrected Spectra')
    plt.show()
          
        
def process_mzxml(file_path, save_as):    
    """
    Reads a mzXML file and writes and saves the modified/corrected file 
    with no oscilations 

    Args:
        file_path ([string]): input file path
        save_as ([string]): output file name

    Returns:
        [none]
    """
    input_map=MSExperiment()
    MzXMLFile.load(file_path, input_map)
    
    #Plot both 
    original_spectra = []
    corrected_spectra = []
   
    for spectrum in input_map:
       original_spectra.append(spectrum.clone())  #original
       corrected_spectra.append(correct_oscilations(spectrum))#corrected
   
    plot_spectra(original_spectra, corrected_spectra)
    
    #Save the corrected spectrum
    # if no other name is given for the file-> se sobreescribe
    output_file=save_as
    MzXMLFile.store(output_file, input_map)
    print(f"Processed file saved as {save_as}")
    
    if __name__ == "__main__":
        file_path="C:/Users/maipa/repos/oscilations_TFG_repository/MS_files"
        input_mzxml = "/CTRL_103_01_c_afterreboot.mzXML"
        save_as= "/filtered_CTRL_103_01_c_afterreboot.mzXML"
        process_mzxml(input_mzxml, save_as)
        
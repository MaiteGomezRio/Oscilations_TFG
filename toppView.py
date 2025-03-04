#This code reads a mzXML file, transforms it to a mzML file and writes on it using pyopensms

import numpy as np
import matplotlib.pyplot as plt
import base64
import zlib
import pyopenms as oms


#pyopnems does this automatically
def decompress_peaks(encoded_data):
    """
    Extraction of the data from file
    Decodes base64 and decompresses encoded (with zlib) intensity values from mzXML.
    """
    decoded = base64.b64decode(encoded_data)
    decompressed = zlib.decompress(decoded)
    return np.frombuffer(decompressed, dtype=np.float64)

def correct_oscilations (spectrum):
    """
    Corrects the oscilations of the MS 
    TODO: apply filtering of that data  
    
    Args:
        original spectrum from file
    Returns: 
        corrected/filtered spectrum
    """
    #Obtain the peaks (data)
    mz, intensity = spectrum.get_peaks()
    
    
    #I can also do directly the set_peaks but just in case:
    #Create instance of an object spectrum to store the corrected data
    corrected_spectrum = oms.MSSpectrum()
    corrected_spectrum.set_peaks((mz, intensity))
    
    #In TOPPView- RT was not correctly saved
    corrected_spectrum.setRT(spectrum.getRT())
    
    #Apply correction-TODO: LOOK FOR FILTERING/SMOOTHING THAT WORKS BEST FOR METABOLITES 
    gf=oms.GaussFilter()
    param = gf.getParameters()
    param.setValue("gaussian_width", 1.0)  # needs wider width?
    gf.setParameters(param)
    gf.filter(corrected_spectrum)
    

    return corrected_spectrum
    

def plot_spectra(original_spectra, corrected_spectra, max_spectra=5):
    """
    Plots the original data obtained from the mzXML file and the new representation
    of the data with corrected oscilations on the same graph with different colors.

    Args:
        original_spectra: List of spectra before correction.
        corrected_spectra: List of spectra after correction.
        max_spectra: Number of spectra to plot (default=5).
    """
    plt.figure(figsize=(12, 6))

    plotted = 0  # Counter for how many spectra we've plotted

    for i, (original, corrected) in enumerate(zip(original_spectra, corrected_spectra)):
        mz_original, intensity_original = original.get_peaks()
        mz_corrected, intensity_corrected = corrected.get_peaks()

        # Check if the spectrum contains peaks and skip empty spectra
        if len(mz_original) == 0 or len(mz_corrected) == 0:
            print(f"Skipping empty spectrum {i}") 
            continue  

        # Plot only non-empty spectra
        plt.plot(mz_original, intensity_original, color='red', alpha=0.6)
        plt.plot(mz_corrected, intensity_corrected, color='blue', alpha=0.8)

        plotted += 1
        if plotted >= max_spectra:  
            break

    if plotted == 0:
        print("No non-empty spectra found to plot!")
    else:
        plt.plot([], [], color='red', label='Original Spectrum')
        plt.plot([], [], color='blue', label='Corrected Spectrum')
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.legend()
        plt.title(f'Original vs Corrected Spectra')#with this file only five spectra are shown
        plt.show()
              
def process_mzxml(file_path, save_as_mzml):   
        
    """
    Reads a mzXML file and writes and saves the modified/corrected file 
    with no oscilations 

    Args:
        file_path ([string]): input file path
        save_as ([string]): output file name

    Returns:
        [none]
    """
    input_map = oms.MSExperiment()
    
    # Load mzXML file
    oms.MzXMLFile().load(file_path, input_map)
    
    # Apply corrections
    corrected_spectra = []
    for spectrum in input_map:
        corrected_spectrum = correct_oscilations(spectrum)
        corrected_spectra.append(corrected_spectrum)
    
    # Save corrected spectra as mzML
    output_map = oms.MSExperiment()
    output_map.setSpectra(corrected_spectra)
    
    oms.MzMLFile().store(save_as_mzml, output_map)
    print(f"Processed mzML file saved as {save_as_mzml}")
    
    
if __name__ == "__main__":
    file_path="C:/Users/maipa/repos/oscilations_TFG_repository/MS_files"
    input_mzxml = "/CTRL_103_01_c_afterreboot.mzXML"
    save_as= "/filtered_CTRL_103_01_c_afterreboot.mzML"
    process_mzxml(file_path+input_mzxml, file_path+save_as)#When I do the final processing 
    #I should put save_as=None so the file is overwritten and corrected 
    # Store as mzML

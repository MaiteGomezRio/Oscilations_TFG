import numpy as np
import matplotlib.pyplot as plt
import pyopenms as oms
from scipy.signal import savgol_filter

def correct_oscillations(spectrum, window_length, filter_order):
    """
    Applies Savitzky-Golay filter to an MS spectrum.
    
    Args:
        spectrum (MSSpectrum): Original mass spectrum.
        window_length (int): Window size for the filter (odd num)
        filter_order (int): order of the filter 
    
    Returns:
        MSSpectrum: Smoothed/filtered spectrum.
    """
    mz, intensity = spectrum.get_peaks()
    
    #Savitzky-Golay filter
    if len(intensity) > window_length: #to ensure that intensity is larger that window length
        smoothed_intensity=savgol_filter(intensity, window_length, filter_order)
    else:
        smoothed_intensity=intensity #no smoothing
        
    # Copy the spectrum to store the corrected data-- correction done because it showed Rt=-1
    corrected_spectrum = oms.MSSpectrum()
    corrected_spectrum.set_peaks((mz, smoothed_intensity))
    corrected_spectrum.setRT(spectrum.getRT())  # Keep retention time! 
    

    return corrected_spectrum

def plot_spectra(original_spectra, corrected_spectra, max_spectra=5):#see if I have to change this max_spectra
    """
    Plots original and Gaussian-filtered spectra for comparison.
    """
    plt.figure(figsize=(12, 6))
    plotted = 0  

    for i, (original, corrected) in enumerate(zip(original_spectra, corrected_spectra)):
        mz_original, intensity_original = original.get_peaks()
        mz_corrected, intensity_corrected = corrected.get_peaks()

        plt.plot(mz_original, intensity_original, color='black', alpha=0.6, label='Original' if plotted == 0 else "")
        plt.plot(mz_corrected, intensity_corrected, color='red', alpha=0.8, label='Corrected' if plotted == 0 else "")

        plotted += 1
        if plotted >= max_spectra:
            break

    if plotted > 0:
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.legend()
        plt.title('Original vs Corrected')
        plt.show()
    else:
        print("No spectra found to plot.")

def process_mzxml(file_path, save_as, window_length, filter_order):
    """
    Reads an mzXML file, applies Gaussian filtering, and saves as mzML.
    """
    input_map = oms.MSExperiment()
    output_map=oms.MSExperiment()
    
    # Load mzXML file
    oms.MzXMLFile().load(file_path, input_map)

    original_spectra = []
    corrected_spectra = []

    for spectrum in input_map:
        original_spectra.append(spectrum)
        corrected_spectrum = correct_oscillations(spectrum, window_length, filter_order)
        corrected_spectra.append(corrected_spectrum)

    plot_spectra(original_spectra, corrected_spectra)

    # Replace spectra with corrected versions
    #input_map.clear(False) if I want to change the original file I will have to remove the input spectra
    input_map.setSpectra(corrected_spectra)

    # Save as mzML
    oms.MzMLFile().store(save_as, input_map)
    print(f"Processed file saved as {save_as}")

if __name__ == "__main__":
    file_path = "C:/Users/maipa/repos/oscilations_TFG_repository/MS_files/CTRL_103_01_c_afterreboot.mzXML"
    save_as = "C:/Users/maipa/repos/oscilations_TFG_repository/MS_files/filtered2_CTRL_103_01_c_afterreboot.mzML"
    process_mzxml(file_path, save_as, window_length=11, filter_order=3)

#Read and write pyteomics

import numpy as np
import matplotlib.pyplot as plt
import base64
import zlib
import pyopenms as oms
from pyteomics import mzxml

# Decompress function (if needed)
def decompress_peaks(encoded_data):
    """
    Extraction of the data from file.
    Decodes base64 and decompresses encoded (with zlib) intensity values from mzXML.
    """
    decoded = base64.b64decode(encoded_data)
    decompressed = zlib.decompress(decoded)
    return np.frombuffer(decompressed, dtype=np.float64)

# Function to correct oscilations
def correct_oscilations(spectrum):
    """
    Corrects the oscillations of the MS using a Gaussian filter.
    
    Args:
        spectrum (MSSpectrum): A spectrum from mzXML.
    
    Returns: 
        corrected spectrum (MSSpectrum)
    """
    # Obtain peaks (m/z and intensity)
    mz, intensity = spectrum.get_peaks()
    
    # Create a new spectrum object to store corrected data
    corrected_spectrum = oms.MSSpectrum()
    corrected_spectrum.set_peaks((mz, intensity))
    
    # Apply Gaussian filter
    gf = oms.GaussFilter()
    param = gf.getParameters()
    param.setValue("gaussian_width", 1.0)  # Adjust this parameter if needed
    gf.setParameters(param)
    gf.filter(corrected_spectrum)

    return corrected_spectrum

# Function to write mzXML using Pyteomics
def write_mzxml_directly(output_file, corrected_spectra):
    """
    Writes spectra data into an mzXML file using Pyteomics.
    
    Args:
        output_file (str): Path to save the mzXML file.
        corrected_spectra (list): List of corrected spectra dictionaries.
    """
    with mzxml.write(output_file, overwrite=True) as writer:
        for spectrum in corrected_spectra:
            writer.write(spectrum)

# Plot function to compare original vs corrected spectra
def plot_spectra(original_spectra, corrected_spectra, max_spectra=5):
    """
    Plots the original spectra obtained from mzXML and the corrected spectra.

    Args:
        original_spectra (list): Spectra before correction.
        corrected_spectra (list): Spectra after correction.
        max_spectra (int): Number of spectra to plot (default=5).
    """
    plt.figure(figsize=(12, 6))

    plotted = 0  # Counter

    for i, (original, corrected) in enumerate(zip(original_spectra, corrected_spectra)):
        mz_original, intensity_original = original.get_peaks()
        mz_corrected, intensity_corrected = corrected.get_peaks()

        if len(mz_original) == 0 or len(mz_corrected) == 0:
            print(f"Skipping empty spectrum {i}") 
            continue  

        plt.plot(mz_original, intensity_original, color='red', alpha=0.6, label="Original Spectrum" if plotted == 0 else "")
        plt.plot(mz_corrected, intensity_corrected, color='blue', alpha=0.8, label="Corrected Spectrum" if plotted == 0 else "")

        plotted += 1
        if plotted >= max_spectra:  
            break

    if plotted == 0:
        print("No non-empty spectra found to plot!")
    else:
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.legend()
        plt.title(f'Original vs Corrected Spectra')
        plt.show()
          
# Main processing function
def process_mzxml(file_path, save_as):    
    """
    Reads an mzXML file, applies corrections, and writes back the corrected file.
    
    Args:
        file_path (str): Input mzXML file path.
        save_as (str): Output mzXML file path.
    
    Returns:
        None
    """
    input_map = oms.MSExperiment()
    
    # Load mzXML file using OpenMS
    oms.MzXMLFile().load(file_path, input_map)
    
    # Prepare lists to store spectra
    original_spectra = []
    corrected_spectra = []

    # Process each spectrum
    for spectrum in input_map:
        original_spectra.append(spectrum)
        corrected_spectrum = correct_oscilations(spectrum)
        corrected_spectra.append(corrected_spectrum)

    # Plot spectra comparison
    plot_spectra(original_spectra, corrected_spectra)

    # Convert OpenMS spectra into Pyteomics-compatible format
    corrected_spectra_pyteomics = []
    for spectrum in corrected_spectra:
        mz, intensity = spectrum.get_peaks()
        corrected_spectra_pyteomics.append({
            'm/z array': np.array(mz),
            'intensity array': np.array(intensity),
            'scan start time': spectrum.getRT(),
            'ms level': 1  # Adjust if necessary
        })

    # Save the corrected spectrum back to mzXML using Pyteomics
    write_mzxml_directly(save_as, corrected_spectra_pyteomics)
    print(f"Processed file saved as {save_as}")

# Run script
if __name__ == "__main__":
    file_path = "C:/Users/maipa/repos/oscilations_TFG_repository/MS_files"
    input_mzxml = "/CTRL_103_01_c_afterreboot.mzXML"
    save_as = "/filtered_pyteomics_CTRL_103_01_c_afterreboot.mzXML"
    process_mzxml(file_path + input_mzxml, file_path + save_as)



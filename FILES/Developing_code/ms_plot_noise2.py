# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from pyopenms import *

def plot_ms2_spectrum(mzML_file_path):
    """
    Loads the mzML file and plots the first spectrum MS2 found.
    
    Args:
        mzML_file_path (str): File path with metabolite data
    Returns:
        null: Shows MS graph with noise   
    """

    # Load the mzML file
    exp=MSExperiment()
    MzMLFile().load(mzML_file_path, exp)

    # Extract MS2 spectra (fragment ions) for better visualization
    ms2_spectra = [s for s in exp if s.getMSLevel() == 2]

    if not ms2_spectra:
        print("No MS2 spectra found in the file.")
    else:
        spectrum = ms2_spectra[0]  # Take the first MS2 spectrum
        mz_values, intensities = spectrum.get_peaks()

    # Plot the spectrum
    plt.figure(figsize=(10, 6))
    plt.stem(mz_values, intensities, basefmt=" ", use_line_collection=True)
    plt.xlabel("m/z (Mass-to-Charge Ratio)")
    plt.ylabel("Intensity")
    plt.title("Mass Spectrum (MS2)")
    plt.xlim(min(mz_values) - 10, max(mz_values) + 10)
    plt.ylim(0, max(intensities) * 1.1)
    plt.grid(True)

    # Show the plot
    plt.show()

def main(): 
    mzML_file_path=r"C:\Users\maipa\.spyder-py3\160528_FHS_Eicosanoid_Plate_13_03.mzML"
    plot_ms2_spectrum(mzML_file_path)
if __name__ == "__main__":
    main()
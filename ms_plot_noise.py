# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from pyopenms import *
from pyteomics import mzxml

def plot_ms2_spectrum(mzXML_file_path):
    """
    Loads the mzXML file and plots the first MS2 spectrum found.

    Args:
        mzXML_file_path (str): File path with metabolite data.
    Returns:
        null: Shows MS graph with noise   
    """

    # Load the mzXML file using pyteomics
    with mzxml.read(mzXML_file_path) as spectra:
        ms2_spectra = [s for s in spectra if s['msLevel'] == 2]

        if not ms2_spectra:
            print("No MS2 spectra found in the file.")
            return
        
        spectrum = ms2_spectra[0]  # Take the first MS2 spectrum
        mz_values = spectrum['m/z array']
        intensities = spectrum['intensity array']

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

if __name__ == "__main__":
    mzXML_file_path = r"C:\Users\maipa\OneDrive - Fundación Universitaria San Pablo CEU\Documentos\CURSOS ACADÉMICOS\4ºIBM\TFG\CTRL_103_01_c_afterreboot.mzXML"
    plot_ms2_spectrum(mzXML_file_path)

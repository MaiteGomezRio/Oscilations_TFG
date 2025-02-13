# -*- coding: utf-8 -*-

#código para leer un archivo de tipo mzXML

from pyteomics import mzxml
import matplotlib.pyplot as plt


# Ruta del archivo .mzXML
mzxml_file = r"C:\Users\maipa\repos\oscilations_TFG_repository\MS_files\CTRL_103_01_c_afterreboot.mzXML"

# Leer el archivo
with mzxml.read(mzxml_file) as spectra:
    for spectrum in spectra:
        print(f"Scan: {spectrum['num']}, m/z array: {spectrum['m/z array'][:5]}, Intensity: {spectrum['intensity array'][:5]}")
        break  # Solo imprime el primer espectro para verificar que funciona

#Plot of the MS

def plot_mass_spectrum(mz_values, intensity_values, scan_number=None):
    """
    Plots a mass spectrum given a list of m/z values and corresponding intensity values.

    Parameters:
    - mz_values (list or np.array): Mass-to-charge ratio (m/z) values.
    - intensity_values (list or np.array): Corresponding intensity values.
    - scan_number (int, optional): The scan number (for labeling the plot).
    """
    if not mz_values or not intensity_values:
        print("No data to plot.")
        return
    
    plt.figure(figsize=(10, 5))
    plt.stem(mz_values, intensity_values, basefmt=" ")

    plt.xlabel("m/z (Mass-to-Charge Ratio)")
    plt.ylabel("Intensity")
    plt.title(f"Mass Spectrum" + (f" (Scan {scan_number})" if scan_number else ""))
    plt.grid(True, linestyle="--", alpha=0.5)

    plt.show()


if __name__ == "__main__":
    mz_values = [60.07272471, 61.06897566, 62.09066769, 63.08439054, 63.14925878]
    intensity_values = [7621.23, 1735.00, 21293.96, 4955.96, 814.26]
    
    plot_mass_spectrum(mz_values, intensity_values, scan_number=3077)

















plt.show()

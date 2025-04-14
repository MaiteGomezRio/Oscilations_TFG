
#This code reads the mzXML file and directly writes on it

from lxml import etree
import numpy as np

def process_and_save_mzxml(input_file, output_file, intensity_percentile=10):
    """
    Reads an mzXML file, applies noise filtering, plots spectra, and saves a filtered version.
    """
    # Parse the mzXML file as an XML tree
    tree = etree.parse(input_file)
    root = tree.getroot()

    # Find all spectrum entries
    for scan in root.findall("entries"):
        #save intensity and mz values
        mz_values = []
        intensity_values = []

        # Extract mz and intensity values from peaks {http://sashimi.sourceforge.net/schema/mzXML_3.2}peaks
        peaks = scan.find("entries")

        # Apply noise filtering
        filtered_mz, filtered_intensity, threshold = filter_noise(mz_values, intensity_values, percentile=intensity_percentile)

        # Plot the scan
        plot_mass_spectra(mz_values, intensity_values, filtered_mz, filtered_intensity, scan_to_plot)

        # Convert filtered values back to binary string
        #COMPLETAR CON CÓDIGO

    # Save the modified mzXML file
    tree.write(output_file)
    print(f"Filtered mzXML saved")
          
          
          
          
          
          
          
          
          
          
          
          d as: {output_file}")
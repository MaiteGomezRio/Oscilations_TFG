import numpy as np
import matplotlib.pyplot as plt#plot of the ms
import pyopenms as oms
from scipy.signal import savgol_filter#smoothing for freq
from scipy.fftpack import fft#for obtaining freq
import os#to identify file extension 
import subprocess#allows to run system commands and interact w/ external programs (proteoWizard)
import time

def convert_mzxml_2_mzml(file_path):#works
    """
        Since pyopenms cannot overwrite a mzxml file, it transforms it to mzml w/ proteoWizard
        
        Args:
            file_path (str): path to the mzXML file
        
        Returns:
        str: Path to the converted mzML file.
    """
    #Msconvert from ProteoWizard
    output_folder = os.path.dirname(file_path)  # Carpeta del archivo original
    mzml_file_path = os.path.join(output_folder, os.path.basename(file_path).replace(".mzXML", ".mzML"))
    try:
        '''subprocess.run([
            "msconvert", file_path, "--mzML", "--outdir", output_folder  # <-- Asegura que se guarda en la misma carpeta
        ], check=True)'''
        subprocess.run([
            "msconvert", file_path, "--mzML", "--outdir", output_folder,"--64", "--zlib",
            "--filter", "peakPicking true 1-"  # Ensures centroiding but keeps all peaks
        ], check=True)

        if os.path.exists(mzml_file_path):
            #print("File correctly generated.")
            return mzml_file_path
        else:
            raise RuntimeError(f"MsConvert could not generate the expected file: {mzml_file_path}")

    except subprocess.CalledProcessError as e:
        print(f"Error executing msconvert: {e}")
        raise RuntimeError("ProteoWizard (msconvert) failed conversion.")

def apply_savgol_filter(intensity, window_length, filter_order):
    """
    Applies a Savitzsky-Golay filter for smoothing of the signal
  
    Args:
        intensity : 
        window_length : 
        filter_order :

    Returns:
        Smoothed intensity after applying the filter

    """
    #print(f"Before filtering: {intensity[:10]} ...")
    if len(intensity) > window_length: #to ensure that intensity is larger that window length
        smoothed_intensity=savgol_filter(intensity, window_length, filter_order)
        #print(f"After filtering: {smoothed_intensity[:10]} ...")
    else:
        smoothed_intensity=intensity #no smoothing
        
    return smoothed_intensity

def calculate_freq(intensity):
    """
    Computes the frequency components of a signal using FFT and from that 
    returns the main freq (highest) --but it should all be constant?
    
    Args:
        intensity: Intensity values of the MS
    Returns:
        frequency: Returns the frequency values after applying FFT
        magnitude: Returns the magnitude of each freq
    """
    #FFT computation
    fft_result=fft(intensity)#how much of each freq there is
    freqs=np.fft.fftfreq(len(intensity))#Frequency
    
    #Obtaining the frequencies (From one because I eliminate the dc component which is always 0)
    fft_freqs=freqs[1:len(freqs) // 2]
    fft_magnitude= np.abs(fft_result[1:len(freqs) // 2])#how strong each freq is (like amplitude kind of)
    
    #main_index = np.argmax(fft_magnitude)
    #main_freq = abs(fft_freqs[main_index])#freq at wich signal oscillates with the highest amplitude
    
    return fft_freqs, fft_magnitude

def obtain_amplitudes(mz, intensity, bin_size):
    """
    Applies binning by m/z values to obtain the amplitudes of the signal
    
    Args:
        mz : m/z values
        intensity : filtered intensity at that point in the spectrum
        bin_size : width of the m/z bin performed

    Returns:
        amplitudes: amplitude per bin
        bin_centers: center of each bin

    """
    min_mz=np.min(mz)
    max_mz=np.max(mz)
    #calculo el bin size de cada spectrum
    bins=np.arange(min_mz,max_mz+bin_size,bin_size)
    amplitudes, bin_edges=np.histogram(mz,bins,weights=intensity)

    #amplitud media por bin?
    bin_centers=(bins[:-1]+bins[1:])/2#halla el centro 
    
    return amplitudes, bin_centers
    
    
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
    
    #1. We apply the filter
    smoothed_intensity=apply_savgol_filter(intensity, window_length, filter_order)
    #CÓDIGO PARA PROBAR EL FILTRO--funciona
    #smoothed_spectrum = oms.MSSpectrum()
    #smoothed_spectrum.set_peaks((mz, smoothed_intensity))
    #smoothed_spectrum.setRT(spectrum.getRT())  # Keep retention time! 
    
    #2. We perform the FFT to obtain the freq
    fft_freqs, fft_magnitude=calculate_freq(smoothed_intensity)

    #3. We obtain amplitude for each m/z values by binning
    amplitudes,bin_centers=obtain_amplitudes(mz, smoothed_intensity, bin_size=0.01)
    
    #4. TODO--We make the corrections
    #corrected_intensity=....
    
    #5.We apply the corrections to the spectrum
    # Copy the spectrum to store the corrected data-- correction done because it showed Rt=-1
    corrected_spectrum=oms.MSSpectrum()
    #corrected_spectrum.set_peaks((mz,corrected_intensity))
    corrected_spectrum.setRT(spectrum.getRT())
    
    #solo devolver lo necesario--borrar las utilizadas para probar código
    return corrected_spectrum, fft_freqs, fft_magnitude, amplitudes, bin_centers 


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



def process_file(file_path, save_as, window_length, filter_order):
    """
    Reads an mzXML/mzML file, applies Savitzsky-Golay filter, and saves as mzML.
    """
    input_map = oms.MSExperiment()
    #output_map=oms.MSExperiment()
    
    #To detect if it is mzML, mzXML
    file_extension = os.path.splitext(file_path)[1]
    file_extension = file_extension.lower()
    if file_extension == ".mzxml":
        # Load mzXML file
        #oms.MzXMLFile().load(file_path, input_map)
        mzml_file_path=convert_mzxml_2_mzml(file_path)
        print("Converting to mzML...")
        time.sleep(3)
        if not os.path.exists(mzml_file_path):
           raise RuntimeError(f"Error: file not found")
        else:
            print("Loading mzML file...")
            oms.MzMLFile().load(mzml_file_path, input_map)
    else:
        # Load mzML file
        oms.MzMLFile().load(file_path, input_map)

    original_spectra = []
    smoothed_spectra = []
    main_freqs=[]

    for spectrum in input_map:
        original_spectra.append(spectrum)
        corrected_spectrum, fft_freqs, fft_magnitude, amplitudes, bin_centers = correct_oscillations(spectrum, window_length, filter_order)

        smoothed_spectra.append(corrected_spectrum)
        #Main freq of each spectrum
        main_freq_index = np.argmax(fft_magnitude)
        main_freq = abs(fft_freqs[main_freq_index])#the freq with the max magnitude
        main_freqs.append(main_freq)#all main freqs for that spectrum
        
    #pARA COMPROBAR QUE EL CÓDIGO TIENE SENTIDO
    mean_freq=np.mean(main_freqs)
    std_freqs=np.std(main_freqs)
    print(f"Overall mean freq: {mean_freq:.10f}")
    print(f"Overall std freq: {std_freqs:.10f}")
    plt.hist(main_freqs, bins='auto')
    plt.xlabel("Main Frequency")
    plt.ylabel("Count")
    plt.title("Distribution of Main Frequencies")
    plt.show()
    #PARA COMPROBAR--- BORRAR AL FINAL
    plt.bar(bin_centers, amplitudes, 0.01)
    plt.xlabel("m/z")
    plt.ylabel("Amplitude (binned)")
    plt.title("Binned Amplitudes by m/z")
    plt.show()
    
    #SAVE THE CHANGES 
    input_map.setSpectra(smoothed_spectra)
    # Replace spectra with corrected versions
    #input_map.clear(False) if I want to change the original file I will have to remove the input spectra
    #Plot of the spectra
    plot_spectra(original_spectra, smoothed_spectra)
    # Save as mzML
    oms.MzMLFile().store(mzml_file_path, input_map)
    print(f"Processed file saved as {save_as}")

if __name__ == "__main__":
    file_path = "C:/Users/maipa/repos/oscilations_TFG_repository/TESTS/CTRL_103_01_c_afterreboot.mzXML"
    save_as = "C:/Users/maipa/repos/oscilations_TFG_repository/TESTS/CTRL_103_01_c_afterreboot.mzML"
    process_file(file_path, save_as, window_length=5, filter_order=3)
    
    
        

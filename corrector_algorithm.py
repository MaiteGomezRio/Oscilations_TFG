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

def apply_savgol_filter(intensities, window_length, filter_order):
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
    if len(intensities) > window_length: #to ensure that intensity is larger that window length
        smoothed_intensities=savgol_filter(intensities, window_length, filter_order)
        #print(f"After filtering: {smoothed_intensity[:10]} ...")
    else:
        smoothed_intensities=intensities #no smoothing
        
    return smoothed_intensities

def calculate_freq(intensities):
    """
    Computes the frequency components of a signal using FFT and from that 
    returns the main freq 
    
    Args:
        intensity: Intensity values of the MS
    Returns:
        frequency: Returns the frequency values after applying FFT
        magnitude: Returns the magnitude of each freq
    """
    #FFT computation
    fft_result=fft(intensities)#how much of each freq there is
    freqs=np.fft.fftfreq(len(intensities))#Frequency
    
    #Obtaining the frequencies (From one because I eliminate the dc component which is always 0)
    fft_freqs=freqs[1:len(freqs) // 2]
    fft_magnitude= np.abs(fft_result[1:len(freqs) // 2])#how strong each freq is (like amplitude kind of)
    
    #main_index = np.argmax(fft_magnitude)
    #main_freq = abs(fft_freqs[main_index])#freq at wich signal oscillates with the highest amplitude
    
    return fft_freqs, fft_magnitude

def obtain_amplitudes(mzs, intensities, bin_size):
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
    min_mz=np.min(mzs)
    max_mz=np.max(mzs)
    #calculo el bin size de cada spectrum
    bins=np.arange(min_mz,max_mz+bin_size,bin_size)
    amplitudes, bin_edges=np.histogram(mzs,bins,weights=intensities)

    #amplitud media por bin?
    bin_centers=(bins[:-1]+bins[1:])/2#halla el centro 
    
    return amplitudes, bin_centers
        
def correct_oscillations(spectrum, window_length, filter_order):
    """
    Corrects the oscillations in the spectrum.
    
    Args:
        spectrum (MSSpectrum): Original mass spectrum.
        window_length (int): Window size for the filter (odd num)
        filter_order (int): order of the filter 
    
    Returns:
        MSSpectrum: Smoothed/filtered spectrum.
    """
    
    """leer las intensidades y Establecer una baseline en el mínimo de los valores, 
    "mover la señal y establecer esa misma linea para la señal reconstruida para que al restar no se 
    me quede ningún valor negativo, """

    mzs, intensities = spectrum.get_peaks()
    
    #1. We apply the filter
    smoothed_intensities=apply_savgol_filter(intensities, window_length, filter_order)
    #CÓDIGO PARA PROBAR EL FILTRO--funciona
    #smoothed_spectrum = oms.MSSpectrum()
    #smoothed_spectrum.set_peaks((mzs, smoothed_intensity))
    #smoothed_spectrum.setRT(spectrum.getRT())  # Keep retention time! 
    
    #2. We perform the FFT to obtain the freq
    fft_freqs, fft_magnitude=calculate_freq(smoothed_intensities)
    main_freq_index = np.argmax(fft_magnitude)
    main_freq = abs(fft_freqs[main_freq_index])
    #plt.hist(fft_freqs, bins='auto')-> para comprobar que la freq se ha reducido
    
    #3. We obtain amplitude for each m/z values by binning
    amplitudes, bin_centers=obtain_amplitudes(mzs, smoothed_intensities, bin_size=0.01)
    
    #TODO--borrar si funciona
    #count=0
    #print("\nAmplitudes por bin (para cada centro de bin):")
    #for center, amp in zip(bin_centers, amplitudes):
        #if amp!=0:
            #print(f"  m/z: {center:.4f}, Amplitud: {amp:.2f}")
            #count += 1
            #if count >= 100:
                #break
    
    #4. We make the corrections
    #como tengo una amplitud asociada a cada bin pero no a cada valor de mz tengo interpolar
    #con coherencia para poder construir la señal correctamente 
    
    interpolated_amplitudes=np.interp(mzs, bin_centers, amplitudes)
    #max_amp = np.max(interpolated_amplitudes)
    
    #TODO--comprobar que está bien y tiene sentido
    oscillation = interpolated_amplitudes * np.sin(2 * np.pi * main_freq * np.arange(len(mzs)))
    #we subtract the oscillation signal created to the original one w/ oscillations
    #So the corrected intensities do not give negative values
    limited_oscillation = np.minimum(oscillation, smoothed_intensities)
    
    corrected_intensities=smoothed_intensities-limited_oscillation
    
    
    #5.We apply the corrections to the spectrum
    # Copy some of the spectrum info to store the corrected data-- correction done because it showed Rt=-1
    corrected_spectrum=oms.MSSpectrum()
    corrected_spectrum.set_peaks((mzs,corrected_intensities))#for m/z range
    corrected_spectrum.setInstrumentSettings(spectrum.getInstrumentSettings())#for polarity
    
    #For mzMine to recognize them    
    corrected_spectrum.setRT(spectrum.getRT())
    corrected_spectrum.setMSLevel(spectrum.getMSLevel())
    corrected_spectrum.setDriftTime(spectrum.getDriftTime())
    corrected_spectrum.setAcquisitionInfo(spectrum.getAcquisitionInfo())
    corrected_spectrum.setNativeID(spectrum.getNativeID())
    
    #ACQUISITION INFO
    info = spectrum.getAcquisitionInfo()
    print("ADQ INFO:", info.size())
    
    #corrected_spectrum.setMSLevel(1)
    #corrected_spectrum.setMetaValue("scan definition", "Full")
    #corrected_spectrum.setMetaValue("filter string", "Full scan")
    #corrected_spectrum.setMetaValue("Polarity", "+")
    
    if spectrum.metaValueExists("scan definition"):
        corrected_spectrum.setMetaValue("scan definition", spectrum.getMetaValue("scan definition"))
    if spectrum.metaValueExists("filter string"):
        corrected_spectrum.setMetaValue("filter string", spectrum.getMetaValue("filter string"))
    if spectrum.metaValueExists("Polarity"):
        corrected_spectrum.setMetaValue("Polarity", spectrum.getMetaValue("Polarity"))
    
    return corrected_spectrum 

def plot_spectra(original_spectra, corrected_spectra, max_spectra):
    """
    Plots a specified number of original and corrected spectra for comparison.

    Args:
        original_spectra (list): List of original spectra (MSSpectrum).
        corrected_spectra (list): List of corrected spectra (MSSpectrum).
        max_spectra (int): Number of spectra to plot.
    """
    total_spectra = min(len(original_spectra), len(corrected_spectra), max_spectra)

    for i in range(total_spectra):
        mz_original, intensity_original = original_spectra[i].get_peaks()
        mz_corrected, intensity_corrected = corrected_spectra[i].get_peaks()

        plt.figure(figsize=(10, 4))
        plt.plot(mz_original, intensity_original, label="Original", alpha=0.6, color="black")
        plt.plot(mz_corrected, intensity_corrected, label="Corrected", alpha=0.8, color="red")
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.title(f"Spectrum {i+1} - Original vs Corrected")
        plt.legend()
        plt.tight_layout()
        plt.show()

def plot_last_spectra(original_spectra, corrected_spectra, n_last_spectra):
    """
    Plots the last n original and corrected spectra for comparison.

    Args:
        original_spectra (list): List of original MSSpectrum objects.
        corrected_spectra (list): List of corrected MSSpectrum objects.
        n_last_spectra (int): Number of last spectra to plot.
    """
    total_spectra = min(len(original_spectra), len(corrected_spectra), n_last_spectra)
    start_idx = -total_spectra  # Start from the end

    for i in range(start_idx, 0):
        mz_original, intensity_original = original_spectra[i].get_peaks()
        mz_corrected, intensity_corrected = corrected_spectra[i].get_peaks()

        plt.figure(figsize=(10, 4))
        plt.plot(mz_original, intensity_original, label="Original", alpha=0.6, color="black")
        plt.plot(mz_corrected, intensity_corrected, label="Corrected", alpha=0.8, color="red")
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.title(f"Spectrum {len(original_spectra) + i + 1} (last {total_spectra})")
        plt.legend()
        plt.tight_layout()
        plt.show()

def plot_and_save_individual_spectra(original_spectra, corrected_spectra, indices_to_plot, output_folder="spectra_plots"):
    """
    Plots and saves selected spectra comparing original and corrected versions.

    Args:
        original_spectra (list): List of original spectra.
        corrected_spectra (list): List of corrected spectra.
        indices_to_plot (list of int): Indices of spectra to plot.
        output_folder (str): Folder where plots are saved.
    """
    import os
    os.makedirs(output_folder, exist_ok=True)

    for i in indices_to_plot:
        mz_orig, int_orig = original_spectra[i].get_peaks()
        mz_corr, int_corr = corrected_spectra[i].get_peaks()

        plt.figure(figsize=(10, 4))
        plt.plot(mz_orig, int_orig, label="Original", color="black", alpha=0.6)
        plt.plot(mz_corr, int_corr, label="Corrected", color="red", alpha=0.7)
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.title(f"Spectrum {i} - Original vs Corrected")
        plt.legend()
        plt.tight_layout()

        # Save as PNG
        plot_path = os.path.join(output_folder, f"spectrum_{i}_comparison.png")
        plt.savefig(plot_path)
        plt.close()
        print(f"Saved: {plot_path}")
        
def process_file(file_path, save_as, window_length, filter_order):
    """
    Reads an mzXML/mzML file, applies Savitzsky-Golay filter, and saves as mzML.
    """
    input_map = oms.MSExperiment()
    #output_map=oms.MSExperiment()
    
    #To detect if it is mzML, mzXML
    file_extension = os.path.splitext(file_path)[1].lower()
    
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
        #mzml_file_path=file_path
        oms.MzMLFile().load(file_path, input_map)

    original_spectra = []
    corrected_spectra = []
    #main_freqs=[]

    for spectrum in input_map:
        original_spectra.append(spectrum)
        corrected_spectrum = correct_oscillations(spectrum, window_length, filter_order)
        if corrected_spectrum is not None:
            corrected_spectra.append(corrected_spectrum)
        #Main freq of each spectrum-CÓDIGO DE COMPROBACIÓN
        #main_freq_index = np.argmax(fft_magnitude)
        #main_freq = abs(fft_freqs[main_freq_index])#the freq with the max magnitude
        #main_freqs.append(main_freq)#all main freqs for that spectrum
        
    
    #PARA COMPROBAR--- BORRAR AL FINAL
    #plt.bar(bin_centers, amplitudes, 0.01)
    #plt.xlabel("m/z")
    #plt.ylabel("Amplitude (binned)")
    #plt.title("Binned Amplitudes by m/z")
    #plt.show()
    
    #SAVE THE CHANGES 
    input_map.setSpectra(corrected_spectra)
    # Replace spectra with corrected versions
    #input_map.clear(False) if I want to change the original file I will have to remove the input spectra
    
    #Plot of the spectra
    #plot_last_spectra(original_spectra, corrected_spectra, n_last_spectra=5)  
    #plot_spectra(original_spectra, corrected_spectra, max_spectra=10)
    #plot_and_save_individual_spectra(original_spectra, corrected_spectra, indices_to_plot=[0, 100, 1000, 1500, 3665])
    # Save as mzML
    oms.MzMLFile().store(save_as, input_map)
    print(f"Corrected file saved: {save_as}")
    
if __name__ == "__main__":
    file_path = "C:/Users/maipa/repos/oscilations_TFG_repository/TESTS/MSCONVERT_CTRL_103_01_c_afterreboot.mzML"
    save_as = "C:/Users/maipa/repos/oscilations_TFG_repository/TESTS/CTRL_103_01_c_afterreboot2.mzML"
    process_file(file_path, save_as, window_length=5, filter_order=3)
    
    
        

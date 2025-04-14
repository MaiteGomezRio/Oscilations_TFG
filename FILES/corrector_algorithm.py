import numpy as np
import matplotlib.pyplot as plt#plot of the ms
import pyopenms as oms
from scipy.signal import savgol_filter#smoothing for freq
from scipy.fftpack import fft#for obtaining freq
import os#to identify file extension 
import subprocess#allows to run system commands and interact w/ external programs (proteoWizard)
import time
from collections import Counter

def validate_fft_spectrum(intensities, sampling_interval, rt=None):
    freqs, magnitude, main_freq = calculate_freq(intensities, sampling_interval)
    
    plt.figure(figsize=(8, 3))
    plt.plot(freqs, magnitude)
    plt.axvline(main_freq, color='red', linestyle='--', label=f"Main freq: {main_freq:.4f}")
    plt.title(f"FFT espectro - RT {rt:.2f} min" if rt else "FFT espectro")
    plt.xlabel("Frecuencia (ciclos/minuto)")
    plt.ylabel("Magnitud")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def validate_sine_vs_signal(mzs, intensities, main_freq):
    x = np.arange(len(mzs))
    pure_sine = np.sin(2 * np.pi * main_freq * x)

    plt.figure(figsize=(10, 3))
    plt.plot(intensities / np.max(intensities), label='Señal normalizada', alpha=0.6)
    plt.plot(pure_sine, label=f"Seno (freq: {main_freq:.4f})", alpha=0.8)
    plt.title("Comparación señal vs seno detectado")
    plt.legend()
    plt.tight_layout()
    plt.show()
    

def validate_calculated_freqs(freqs, rts, tic, bins=100):
    """
    Valida que calculate_freq haya estimado correctamente la frecuencia oscilatoria.

    - Muestra un histograma de las frecuencias obtenidas por espectro.
    - Reconstruye una oscilación con la frecuencia más común.
    - La compara visualmente contra el TIC real.

    Args:
        freqs (list or np.array): Frecuencias estimadas por calculate_freq().
        rts (list or np.array): Tiempo de retención por espectro.
        tic (list or np.array): Total ion current (suma de intensidades por espectro).
        bins (int): Número de bins para el histograma.
    """
    freqs = np.array(freqs)
    rts = np.array(rts)
    tic = np.array(tic)

    # Histograma de frecuencias
    rounded_freqs = np.round(freqs, 5)

    plt.figure(figsize=(8, 4))
    plt.hist(rounded_freqs, bins=bins, color="skyblue", edgecolor="black")
    plt.xlabel("Frecuencia estimada (ciclos/minuto)")
    plt.ylabel("Número de espectros")
    plt.title("Distribución de frecuencias por `calculate_freq()`")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

    # Encuentra la frecuencia dominante
    most_common_freq, count = Counter(rounded_freqs).most_common(1)[0]
    print(f" Frecuencia más común: {most_common_freq:.5f} ciclos/minuto ({count} espectros)")

    # Reconstrucción de oscilación con la frecuencia dominante
    amplitude = np.max(tic) * 0.3  # Ajustable
    #reconstructed = amplitude * np.sin(2 * np.pi * most_common_freq * rts)
    
    oscillation = np.array([
        amplitude * np.sin(2 * np.pi * f * t)
        for f, t in zip(freqs, rts)
    ])

    # Normalización segura
    def safe_normalize(arr):
        min_val, max_val = np.min(arr), np.max(arr)
        return (arr - min_val) / (max_val - min_val) if max_val != min_val else arr * 0

    tic_norm = safe_normalize(tic)
    osc_norm = safe_normalize(oscillation)

    # Plot de comparación
    plt.figure(figsize=(12, 4))
    plt.plot(rts, tic_norm, label="TIC Real (normalizado)", color="orange")
    plt.plot(rts, osc_norm, label="Oscilación reconstruida", color="blue", alpha=0.7)
    plt.xlabel("Tiempo de Retención (RT)")
    plt.ylabel("Intensidad Normalizada")
    plt.title("Validación de frecuencia dominante según `calculate_freq()`")
    plt.legend()
    plt.tight_layout()
    plt.show()

def validate_freqs(freqs, rts, tic):
    """
    Valida las frecuencias detectadas vs tiempo de retención y TIC.
    Args:
        freqs: Lista de frecuencias principales por espectro
        rts: Lista de tiempos de retención
        tic: Lista de TIC (total ion current) por espectro
    """
    plt.figure(figsize=(12, 4))
    
    # Frecuencia vs RT
    plt.subplot(1, 2, 1)
    plt.plot(rts, freqs, 'o-', label="Frecuencia principal", color="blue")
    plt.xlabel("Tiempo de retención (min)")
    plt.ylabel("Frecuencia (ciclos/minuto)")
    plt.title("Frecuencia principal vs RT")
    plt.grid(True)
    
    # Frecuencia vs TIC
    plt.subplot(1, 2, 2)
    plt.plot(tic, freqs, 'x-', label="Frecuencia principal", color="green")
    plt.xlabel("TIC")
    plt.ylabel("Frecuencia (ciclos/minuto)")
    plt.title("Frecuencia principal vs TIC")
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

def plot_freq_histogram(freqs):
    """
    Histograma de frecuencias principales detectadas.
    """
    plt.figure(figsize=(6, 4))
    plt.hist(freqs, bins=30, color='skyblue', edgecolor='black')
    plt.xlabel("Frecuencia principal")
    plt.ylabel("Número de espectros")
    plt.title("Distribución de Frecuencias Principales")
    
    print("Stats:")
    print(f" - Frecuencia mínima: {np.min(freqs):.4f}")
    print(f" - Media: {np.mean(freqs):.4f}")
    print(f" - Mediana: {np.median(freqs):.4f}")
    print(f" - Máxima: {np.max(freqs):.4f}")
    plt.tight_layout()
    plt.show()

def validate_amplitudes(mzs, intensities, bin_size=0.01):
    """
    Visualiza los bins de amplitud generados y la señal original.
    """
    amplitudes, bin_centers = obtain_amplitudes(mzs, intensities, bin_size)

    plt.figure(figsize=(10, 4))
    plt.plot(mzs, intensities, label="Señal original", alpha=0.6)
    plt.stem(bin_centers, amplitudes, linefmt='r-', markerfmt='ro', basefmt=' ', label='Amplitudes por bin')
    plt.xlabel("m/z")
    plt.ylabel("Intensidad")
    plt.title("Validación de Amplitudes Binned")
    plt.legend()
    plt.tight_layout()
    plt.show()
    
def validate_baseline(mzs, intensities, smoothed_intensities, oscillation, corrected_intensities, baseline):
    """
    Plotea la señal con la baseline y oscilación reconstruida para validar su efecto.
    """
    plt.figure(figsize=(12, 4))
    plt.plot(mzs, intensities, label="Original", color="orange", alpha=0.5)
    plt.plot(mzs, smoothed_intensities, label="Smoothed", color="green")
    plt.plot(mzs, corrected_intensities, label="Corrected", color="black")
    plt.plot(mzs, oscillation, label="Oscillation", color="blue", linestyle='--')
    plt.axhline(y=baseline, color='grey', linestyle='--', label="Baseline")
    plt.title("Validación de baseline y oscilación reconstruida")
    plt.xlabel("m/z")
    plt.ylabel("Intensidad")
    plt.legend()
    plt.tight_layout()
    plt.show()
    
def validate_global_baseline(rts, tic, global_baseline):
    
    plt.figure(figsize=(12, 4))
    plt.plot(rts, tic, label="TIC original", color='black')
    plt.axhline(y=global_baseline, linestyle='--', color='grey', label='Baseline global (3º percentil TIC)')
    plt.xlabel("Tiempo de retención (min)")
    plt.ylabel("TIC")
    plt.title("TIC vs RT con baseline global")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
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

def calculate_freq(intensities, sampling_interval=1.0, plot_spectrum=False):
    """
    Computes the frequency components of a signal using FFT and from that 
    returns the main freq 
    
    Args:
        intensity: Intensity values of the MS
    Returns:
        frequency: Returns the frequency values after applying FFT
        magnitude: Returns the magnitude of each freq
    """
    
    centered_signal = intensities - np.mean(intensities)
    
    # FFT
    fft_result = fft(centered_signal)
    freqs = np.fft.fftfreq(len(centered_signal), d=sampling_interval)

    # Solo nos quedamos con frecuencias positivas
    pos_mask = freqs > 0
    fft_freqs = freqs[pos_mask]
    fft_magnitude = np.abs(fft_result[pos_mask])
    
    #la freq principal
    main_freq = fft_freqs[np.argmax(fft_magnitude)]
    
    
    # Plot opcional del espectro de frecuencias
    if plot_spectrum:
        plt.figure(figsize=(10, 4))
        plt.plot(fft_freqs, fft_magnitude, color='darkgreen')
        plt.title("Espectro de Frecuencia (FFT)")
        plt.xlabel("Frecuencia (ciclos/minuto)")
        plt.ylabel("Magnitud")
        plt.tight_layout()
        plt.show()
    
    return fft_freqs, fft_magnitude, main_freq

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
       
def correct_oscillations(spectrum, window_length, filter_order, global_baseline):
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
    
    if len(mzs) < window_length:
        return spectrum  # Skip correction if spectrum too small
    
    #1. We apply the filter
    smoothed_intensities = apply_savgol_filter(intensities, window_length, filter_order)
    smoothed_intensities = np.clip(smoothed_intensities, 0, None)#porque salen algunas negativas

    
    #2. We perform the FFT to obtain the freq
    sampling_interval=np.mean(np.diff(mzs))
    fft_freqs, fft_magnitude, main_freq=calculate_freq(smoothed_intensities, sampling_interval)
    
    
    
    #3. We obtain amplitude for each m/z values by binning
    amplitudes, bin_centers=obtain_amplitudes(mzs, smoothed_intensities, bin_size=0.01)
    #como tengo una amplitud asociada a cada bin pero no a cada valor de mz tengo interpolar
    #con coherencia para poder construir la señal correctamente 
    interpolated_amplitudes=np.interp(mzs, bin_centers, amplitudes)
    
    
    #4. We make the corrections
    
    #Creo la señal oscilatoria
    
    x=np.arange(len(mzs))
    oscillation = interpolated_amplitudes * np.sin(2 * np.pi * main_freq * x)
    
    
    #we establish a baseline so we can substract the signal correctly and to avoid negative values 
    baseline = np.quantile(smoothed_intensities, 0.03)  # w/ np.min(smoothed_intensities) coge la mínima pero claro la mínima se sale de lo que es "la media"
    #corrected_values = smoothed_intensities - oscillation
    #corrected_intensities = np.clip(corrected_values, global_baseline, None)  # Para que no baje de la línea base
    oscillation_correction = smoothed_intensities - np.abs(oscillation)
    corrected_intensities = np.clip(oscillation_correction, baseline, None)
    
    
    
    #intensity_threshold = np.percentile(intensities, 98)#si los picos son extremadamente altos-> raro no? el pico ese que sale-> con esto se arregla
    #si es un pico que no es real se tiene que corregir (true) y si no no se tiene que corregir (false)
    #mask = intensities < intensity_threshold  # true o false-> solo se corrige si no es un pico real porque si no corrige zonas de la señal que no son
    #corrected_intensities=smoothed_intensities.copy()
    #mask2 = rts < 15
    #corrected_values = smoothed_intensities[mask] - oscillation[mask]
    #corrected_intensities[mask] = np.clip(corrected_values, baseline, None)#para que solo se corrijan las que tienen que corregirse


    
    #5.We apply the corrections to the spectrum
    # Copy some of the spectrum info to store the corrected data-- correction done because it showed Rt=-1
    corrected_spectrum=oms.MSSpectrum()
    corrected_spectrum.set_peaks((mzs,corrected_intensities))#for m/z range
    
    #For mzMine to recognize them    
    corrected_spectrum.setRT(spectrum.getRT())
    corrected_spectrum.setMSLevel(spectrum.getMSLevel())
    corrected_spectrum.setDriftTime(spectrum.getDriftTime())
    corrected_spectrum.setAcquisitionInfo(spectrum.getAcquisitionInfo())
    corrected_spectrum.setNativeID(spectrum.getNativeID())
    corrected_spectrum.setInstrumentSettings(spectrum.getInstrumentSettings())#for polarity
    
    #Force metadata to appear in mzMine
    for meta in ["scan definition", "filter string", "Polarity"]:
        if spectrum.metaValueExists(meta):
            corrected_spectrum.setMetaValue(meta, spectrum.getMetaValue(meta))
            
    global plot_count#debug variable
    
    # VALIDACIÓN PARA LOS PRIMEROS 10 ESPECTROS
    global plot_count
    if plot_count < 10:
        print(f"--- Validación #{plot_count + 1} ---")
        validate_amplitudes(mzs, smoothed_intensities)
        validate_baseline(mzs, intensities, smoothed_intensities, oscillation, corrected_intensities, baseline)
        
        # VALIDACIÓN EXTRA: forma del seno que estoy generando
        x = np.arange(len(mzs))
        pure_sine = np.sin(2 * np.pi * main_freq * x)
        plt.figure(figsize=(10, 3))
        plt.plot(x, pure_sine, color='purple')
        plt.title(f"Seno puro con frecuencia principal detectada: {main_freq:.4f}")
        plt.xlabel("Índice de punto")
        plt.ylabel("Valor de seno")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    
    plot_count += 1
    
    return corrected_spectrum, main_freq 

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
    freqs=[]#para poder comprobar la freq
    rts=[]#para comprobar
    tic=[]#para comprobar
    global_baseline = None

    for spectrum in input_map:
        original_spectra.append(spectrum)
        
        mzs, intensities=spectrum.get_peaks()#debug
        rt=spectrum.getRT()#debug
        rts.append(rt)#debug
        tic.append(np.sum(intensities))
        global_baseline = np.quantile(tic, 0.03)
        corrected_spectrum, freq = correct_oscillations(spectrum, window_length, filter_order,global_baseline)
        
        if corrected_spectrum is not None and freq is not None:
            corrected_spectra.append(corrected_spectrum)
            freqs.append(freq)#debug
    
    #SAVE THE CHANGES 
    input_map.setSpectra(corrected_spectra)
    # Replace spectra with corrected versions
    #input_map.clear(False) if I want to change the original file I will have to remove the input spectra
    # Save as mzML
    oms.MzMLFile().store(save_as, input_map)
    print(f"Corrected file saved: {save_as}")
    
    # VALIDACIÓN DE FRECUENCIAS 
    print("\n--- VALIDACIÓN DE FRECUENCIAS DETECTADAS ---")

    # Visualización espectro por espectro (los 5 primeros)
    for i in range(min(5, len(original_spectra))):
        print(f"\n[Espectro #{i+1}] RT: {rts[i]:.2f} min | Freq: {freqs[i]:.4f}")
        mzs, intensities = original_spectra[i].get_peaks()
        sampling_interval = np.mean(np.diff(mzs))
        validate_fft_spectrum(intensities, sampling_interval, rts[i])
        validate_sine_vs_signal(mzs, intensities, freqs[i])
    
    # Histograma de frecuencias
    plot_freq_histogram(freqs)
    #validate_freqs(freqs, rts, tic)
    
    
    #VALIDACIÓN DEL BASELINE 
    validate_global_baseline(rts, tic, global_baseline)
    
    
    
    
  
if __name__ == "__main__":
    plot_count=0
    file_path = "C:/Users/maipa/repos/oscilations_TFG_repository/TESTS/MSCONVERT_CTRL_103_01_c_afterreboot_original.mzML"
    save_as = "C:/Users/maipa/repos/oscilations_TFG_repository/TESTS/CTRL_103_01_c_afterreboot_corrected_sinfase.mzML"
    if os.path.exists(save_as):
        print(f"El archivo {save_as} ya existe. Eliminándolo para reemplazarlo...")
        os.remove(save_as)
    process_file(file_path, save_as, window_length=11, filter_order=5)#he cambiado esto 08-04-25
    
    
        

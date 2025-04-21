# processing/corrector.py
import numpy as np
import pyopenms as oms
from utils.filters import apply_savgol_filter
from utils.fft_utils import calculate_freq, obtain_amplitudes, estimate_local_freqs
from validation.signal_validator import validate_amplitudes, validate_baseline
import matplotlib.pyplot as plt


def correct_oscillations(spectrum, window_length, filter_order):
    global plot_count
    
    #Obtengo los mzs e intensidades de la señal original
    mzs, intensities = spectrum.get_peaks()
    if len(mzs) < window_length:
        return spectrum, None
    
    #1. Aplico el filtro a la señal para hacer un smothing
    smoothed_intensities = apply_savgol_filter(intensities, window_length, filter_order)
    smoothed_intensities = np.clip(smoothed_intensities, 0, None)#para que sean todas positivas--no haría falta

    #2. Obtengo las frecuencias para cada una de las intensidades del espectro y las guardo en un array
    #sampling_interval = np.mean(np.diff(mzs))
    #fft_freqs, fft_magnitude = calculate_freq(smoothed_intensities, sampling_interval)
    local_freqs=estimate_local_freqs(smoothed_intensities, mzs)
    
    #3. Obtengo las amplitudes para cada uno de los mz de ese espectro
    amplitudes, bin_centers = obtain_amplitudes(mzs, smoothed_intensities, bin_size=0.01)
    interpolated_amplitudes = np.interp(mzs, bin_centers, amplitudes)
    
    #4. Me creo una nueva "señal oscilatoria" de intensidades y le resto a cada smooth_intensity (de la original)
    #la frecuencia en esa intensidad por la amplitud para esa intensidad en ese mz
    #a esa señal creada le hago clip en el baseline
    #Inew = I_i - f(i) * A_i--> 
    #Ii es la intensidad original (suavizada). f€[0,1]
    #f(I_i) = f(i) es la frecuencia local en el punto i.
    #Ai = amplitud interpolada en el punto i
    new_intensities = smoothed_intensities - local_freqs * interpolated_amplitudes

    #5. Tengo las intensidades de la señal original y les resto las intensidades nuevas de "mi oscilación"
    corrected_intensities=smoothed_intensities-new_intensities
    
    # Establezco el baseline (cogiendo la señal suavizada)
    baseline=np.min(smoothed_intensities)
    #baseline = np.quantile(smoothed_intensities, 0.03)
    #corrected_baseline = np.quantile(corrected_intensities, 0.03)
    #corrected_intensities = np.clip(corrected_intensities, 0, None)
    
    #print(f"Corrected intensities: {corrected_intensities}")
    corrected_intensities+=baseline
    
    #6. Guardo los cambios del espectro, fuerzo los parámetros para que no se pierda información y devuelvo el espectro corregido
    corrected_spectrum = oms.MSSpectrum()
    corrected_spectrum.set_peaks((mzs, corrected_intensities))
    
    corrected_spectrum.setRT(spectrum.getRT())
    corrected_spectrum.setMSLevel(spectrum.getMSLevel())
    corrected_spectrum.setDriftTime(spectrum.getDriftTime())
    corrected_spectrum.setAcquisitionInfo(spectrum.getAcquisitionInfo())
    corrected_spectrum.setNativeID(spectrum.getNativeID())
    corrected_spectrum.setInstrumentSettings(spectrum.getInstrumentSettings())

    for meta in ["scan definition", "filter string", "Polarity"]:
        if spectrum.metaValueExists(meta):
            corrected_spectrum.setMetaValue(meta, spectrum.getMetaValue(meta))
    
    count=0
    #if count<=10:
        #print(f"[DEBUG] RT={spectrum.getRT():.2f} | baseline={baseline:.2e} | corrected_baseline={corrected_baseline:.2e} | offset={offset:.2e}")
        #mzs, spectrum_intensities=corrected_spectrum.get_peaks()
        #print(f"Intesities: {spectrum_intensities[count]}")
        
    
    return corrected_spectrum



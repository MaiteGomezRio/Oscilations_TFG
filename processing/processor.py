
# processing/processor.py
import pyopenms as oms
import os
import time
import numpy as np
from collections import Counter
from processing.corrector import correct_oscillations
from utils.io_utils import convert_mzxml_2_mzml
from validation.signal_validator import validate_fft_spectrum, plot_freq_histogram, validate_sine_vs_signal, validate_global_baseline

plot_count = 0  # Global for debugging

def process_file(file_path, save_as, window_length, filter_order):
    input_map = oms.MSExperiment()

    file_extension = os.path.splitext(file_path)[1].lower()

    if file_extension == ".mzxml":
        mzml_file_path = convert_mzxml_2_mzml(file_path)
        print("Converting to mzML...")
        time.sleep(3)
        if not os.path.exists(mzml_file_path):
            raise RuntimeError(f"Error: file not found")
        else:
            print("Loading mzML file...")
            oms.MzMLFile().load(mzml_file_path, input_map)
    else:
        oms.MzMLFile().load(file_path, input_map)

    original_spectra = []
    corrected_spectra = []
   
    rts = []
    tic = []
    global_baseline = None
    print("Processing file....")
    print(">>> Empezando corrección")
    for spectrum in input_map:
        original_spectra.append(spectrum)
        mzs, intensities = spectrum.get_peaks()
        rt = spectrum.getRT()
        rts.append(rt)
        tic.append(np.sum(intensities))
        global_baseline = np.quantile(tic, 0.03)
        corrected_spectrum = correct_oscillations(spectrum, window_length, filter_order, global_baseline)

        if corrected_spectrum is not None:
            corrected_spectra.append(corrected_spectrum)
    print("<<< Corrección terminada") 

    '''# 1. Recolectar todas las intensidades corregidas
    all_corrected_intensities = []
   
    for spectrum in corrected_spectra:
       _, intensities = spectrum.get_peaks()
       all_corrected_intensities.extend(intensities)

    # 2. Calcular baseline actual de los espectros corregidos
    corrected_baseline = np.quantile(all_corrected_intensities, 0.03)

    # 3. Calcular cuánto hay que subir
    offset = global_baseline - corrected_baseline
    print(f"> Offset aplicado a todos los espectros: {offset:.2e}")

    # 4. Aplicar desplazamiento a todos los espectros corregidos
    for spectrum in corrected_spectra:
        mzs, intensities = spectrum.get_peaks()
        shifted_intensities = intensities + offset
        spectrum.set_peaks((mzs, shifted_intensities))'''       

    input_map.setSpectra(corrected_spectra)
    oms.MzMLFile().store(save_as, input_map)
    print(f"Corrected file saved: {save_as}")


    

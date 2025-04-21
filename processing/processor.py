
# processing/processor.py
import pyopenms as oms
import os
import time
import numpy as np
from collections import Counter
from processing.corrector import correct_oscillations
from utils.io_utils import convert_mzxml_2_mzml
from validation.signal_validator import validate_fft_spectrum, plot_freq_histogram, validate_sine_vs_signal, validate_global_baseline
import time
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
    print(">>> Corrigiendo")
    start_time=time.time()
    for spectrum in input_map:
        original_spectra.append(spectrum)
        mzs, intensities = spectrum.get_peaks()
        rt = spectrum.getRT()
        rts.append(rt)
        tic.append(np.sum(intensities))
        
    #global_baseline = np.quantile(tic, 0.03)#por qué hago quantile? calcula el percentil 3 porque si lo hago con min coge un pico demasiado pequeño
    global_baseline=np.min(tic)
    print(f"[DEBUG] global_baseline calculado: {global_baseline:.2e}")
    

    for spectrum in input_map:
        corrected_spectrum = correct_oscillations(spectrum, window_length, filter_order, global_baseline)
        if corrected_spectrum is not None:
            corrected_spectra.append(corrected_spectrum)
    end_time=time.time()
    elapsed_time=end_time-start_time
    print("<<< Corrección terminada") 
    print(f"El programa ha tardado {elapsed_time:.2f} secs")#3mins más o menos     

    input_map.setSpectra(corrected_spectra)
    oms.MzMLFile().store(save_as, input_map)
    print(f"Corrected file saved: {save_as}")


    

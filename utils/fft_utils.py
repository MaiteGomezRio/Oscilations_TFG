# utils/fft_utils.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft

def calculate_freq(intensities, sampling_interval=1.0, plot_spectrum=False):
    centered_signal = intensities - np.mean(intensities)
    fft_result = fft(centered_signal)
    freqs = np.fft.fftfreq(len(centered_signal), d=sampling_interval)

    pos_mask = freqs > 0
    fft_freqs = freqs[pos_mask]
    fft_magnitude = np.abs(fft_result[pos_mask])
    

    if plot_spectrum:
        plt.figure(figsize=(10, 4))
        plt.plot(fft_freqs, fft_magnitude, color='darkgreen')
        plt.title("Espectro de Frecuencia (FFT)")
        plt.xlabel("Frecuencia (ciclos/minuto)")
        plt.ylabel("Magnitud")
        plt.tight_layout()
        plt.show()
        
    main_freq = freqs[np.argmax(fft_magnitude)]

    return fft_freqs, fft_magnitude, main_freq

def estimate_local_freqs(intensities, mzs, window_size=21):
    """
    Para cada punto i, toma una pequeña ventana de intensidades y sus m/z reales.
    Calcula el sampling_interval solo para esa ventana.
    Aplica calculate_freq usando ese intervalo.
    Guarda la frecuencia dominante local f(i).
    
    Args:
        
    Return:    
        
    """
    half_w = window_size // 2
    local_freqs = []

    for i in range(len(intensities)):
        # Define ventana centrada en i
        start = max(0, i - half_w)
        end = min(len(intensities), i + half_w + 1)

        window_intensities = intensities[start:end]
        window_mzs = mzs[start:end]

        # Asegurar tamaño uniforme rellenando si hace falta
        if len(window_intensities) < window_size:
            pad_len = window_size - len(window_intensities)
            window_intensities = np.pad(window_intensities, (0, pad_len), mode='constant')
            window_mzs = np.pad(window_mzs, (0, pad_len), mode='edge')  # repetir bordes

        # Calcular intervalo de muestreo local
        local_sampling_interval = np.mean(np.diff(window_mzs))

        # Calcular frecuencia dominante en esa ventana
        _, _, main_freq = calculate_freq(window_intensities, local_sampling_interval)
        local_freqs.append(main_freq)

    return np.array(local_freqs)


def obtain_amplitudes(mzs, intensities, bin_size):
    min_mz = np.min(mzs)
    max_mz = np.max(mzs)
    bins = np.arange(min_mz, max_mz + bin_size, bin_size)
    amplitudes, _ = np.histogram(mzs, bins, weights=intensities)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    return amplitudes, bin_centers


# validation/signal_validator.py
import numpy as np
import matplotlib.pyplot as plt
from utils.fft_utils import calculate_freq

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

def validate_amplitudes(mzs, intensities, bin_size=0.01):
    from utils.fft_utils import obtain_amplitudes

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


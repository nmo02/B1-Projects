import numpy as np
import matplotlib.pyplot as plt
import wfdb
import numpy as np
from scipy.fft import fft, fftfreq

def load_ecg_data(record_name):
    record = wfdb.rdrecord(record_name)
    signal = record.p_signal[:, 0]  # Extract first signal channel
    sampling_rate = record.fs
    return signal, sampling_rate

def smooth_signal(signal, sigma=5):
    N = len(signal)
    x = np.linspace(-N//2, N//2, N)
    gaussian_window = np.exp(-0.5 * (x / sigma)**2)
    gaussian_window /= np.sum(gaussian_window)  # Normalize window
    return np.convolve(signal, gaussian_window, mode='same')

def ai_fft(signal, sampling_rate):
    N = len(signal)
    T = 1.0 / sampling_rate
    fft_output = fft(signal)
    frequencies = fftfreq(N, T)[:N // 2]
    amplitudes = 2.0 / N * np.abs(fft_output[:N // 2])
    return frequencies, amplitudes

def manual_fft(signal, sampling_rate):
    N = len(signal)
    T = 1.0 / sampling_rate
    if N & (N - 1) != 0:
        next_pow2 = 2**int(np.ceil(np.log2(N)))
        padded_signal = np.zeros(next_pow2)
        padded_signal[:N] = signal
        signal = padded_signal
        N = len(signal)
    
    fft_output = fft(signal)
    frequencies = np.arange(N) / (N * T)
    amplitudes = 2.0 / N * np.abs(fft_output[:N // 2])
    return frequencies[:N // 2], amplitudes

def main():
    patient_records = ["20800001", "20800002", "20800003", "20800004", "20800005"]
    fig, axs = plt.subplots(3, 2, figsize=(16, 10))
    axs = axs.flatten()[:-1]
    
    for idx, record_name in enumerate(patient_records):
        ecg_signal, sampling_rate = load_ecg_data(record_name)
        ecg_segment = ecg_signal[:sampling_rate * 10]
        smoothed_signal = smooth_signal(ecg_segment)
        
        frequencies_ai, amplitudes_ai = ai_fft(smoothed_signal, sampling_rate)
        frequencies_manual, amplitudes_manual = manual_fft(smoothed_signal, sampling_rate)
        
        axs[idx].plot(frequencies_ai, amplitudes_ai, label="AI FFT (Smoothed)")
        axs[idx].plot(frequencies_manual, amplitudes_manual, linestyle='dashed', label="Self FFT (Smoothed)")
        axs[idx].set_title(f"FFT Comparison on Record {record_name}")
        axs[idx].set_xlabel("Frequency (Hz)")
        axs[idx].set_ylabel("Amplitude")
        axs[idx].legend()
        axs[idx].grid()
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

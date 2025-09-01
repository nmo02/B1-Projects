import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# Recursive FFT implementation
def self_FFT(signal):
    N = len(signal)
    
    # Pad the signal to the next power of 2 if needed
    if N & (N - 1) != 0:  # Check if N is not a power of 2
        next_power_of_2 = 2 ** (N - 1).bit_length()
        padded_signal = np.zeros(next_power_of_2, dtype=complex)
        padded_signal[:N] = signal
        signal = padded_signal
        N = len(signal)

    if N <= 1:
        return np.array(signal, dtype=complex)  # Base case: return the signal itself

    # Recursive calls for even and odd indices
    X_even = self_FFT(signal[::2])
    X_odd = self_FFT(signal[1::2])

    # Combine results with twiddle factors
    factor = np.exp(-2j * np.pi * np.arange(N) / N)
    X = np.zeros(N, dtype=complex)
    X[:N // 2] = X_even + factor[:N // 2] * X_odd
    X[N // 2:] = X_even - factor[:N // 2] * X_odd
    return X

# Wrapper to compute frequencies and FFT
def manual_fft(signal, sampling_rate):
    N = len(signal)  # Length of the signal
    T = 1 / sampling_rate  # Sampling interval

    # Compute FFT using self_FFT
    fft_output = self_FFT(signal)

    # Generate frequency components
    frequencies = np.arange(N) / (N * T)
    half_N = N // 2

    # Normalize the FFT output (divide by N once)
    normalized_fft_output = 2.0 / N * np.abs(fft_output[:half_N])

    return frequencies[:half_N], normalized_fft_output

# AI FFT using SciPy
def ai_FFT(signal, sampling_rate):
    N = len(signal)
    T = 1 / sampling_rate
    fft_output = fft(signal)
    frequencies = fftfreq(N, T)[:N//2]
    return frequencies, 2.0 / N * abs(fft_output[:N//2])

# Apply smoothing to reduce noise in the amplitude spectrum
def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

# Main function
# Main function with updated smoothing
def main():
    # Generate a synthetic signal
    sampling_rate = 100  # Hz
    duration = 10  # seconds
    t = np.linspace(0, duration, sampling_rate * duration, endpoint=False)
    signal = np.sin(2 * np.pi * 5 * t) + np.sin(2 * np.pi * 15 * t)

    # Manual FFT
    xf_manual, yf_manual = manual_fft(signal, sampling_rate)

    # AI FFT
    xf_ai, yf_ai = ai_FFT(signal, sampling_rate)

    # Apply Gaussian smoothing
    smoothing_sigma = 1  # Standard deviation of the Gaussian kernel
    yf_manual_smooth = smooth(yf_manual, smoothing_sigma)
    yf_ai_smooth = smooth(yf_ai, smoothing_sigma)

    # Plot results for comparison
    plt.plot(xf_ai, yf_ai_smooth, label="AI FFT (Smoothed)")
    plt.plot(xf_manual, yf_manual_smooth, linestyle='dashed', label="Self FFT (Smoothed)")
    plt.legend()
    plt.title("FFT Comparison of a Synthetic Signal")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude")
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()

import numpy as np
from astropy.io import fits
import os

def create_dummy_data():
    # Parameters
    width = 50
    height = 50
    n_frames_A = 100
    n_frames_B = 120
    framerate_A = 10.0 # 10 Hz
    framerate_B = 15.0 # 15 Hz

    true_offset = 0.87 # seconds.
    # We want A(t) matches B(t + true_offset)
    # A(t) = S(t)
    # B(t') = S(t' - true_offset)
    # Check: B(t + true_offset) = S(t + true_offset - true_offset) = S(t) = A(t). Correct.

    # Create a matrix for the quadratic relationship
    n_pixels = height * width
    intermediate_dim = 100  # d from the user's example
    np.random.seed(42)  # for reproducibility
    M_linear = np.random.randn(height , width)
    Mquadra = np.random.randn(height , width)


    # Generate random spatial patterns (eigenimages)
    n_modes = 5
    spatial_modes = np.random.randn(n_modes, height, width)

    # Generate time series for modes - Gaussian Pulse
    t_end = max(n_frames_A/framerate_A, n_frames_B/framerate_B + true_offset) + 5
    t_fine = np.linspace(0, t_end, 10000)

    temporal_modes = np.zeros((n_modes, len(t_fine)))
    center_time = 5.0
    for i in range(n_modes):
        # Distinct pulses
        sigma = 1.0
        temporal_modes[i] = np.exp(-((t_fine - center_time - i*0.5)**2) / (2 * sigma**2))

    # Sample for Dataset A
    times_A = np.arange(n_frames_A) / framerate_A
    data_A = np.zeros((n_frames_A, height, width))

    for i, t in enumerate(times_A):
        idx = np.searchsorted(t_fine, t)
        if idx >= len(t_fine): idx = len(t_fine) - 1

        frame = np.zeros((height, width))
        for m in range(n_modes):
            frame += spatial_modes[m] * temporal_modes[m, idx]

        data_A[i] = frame + 0.01 * np.random.randn(height, width)

    # Sample for Dataset B
    # B(t_rec) = Signal(t_rec - true_offset)

    times_B = np.arange(n_frames_B) / framerate_B
    # Shift B start time so we have overlap.
    # Pulse is around t=5.
    # A sees 0..10. A sees pulse.
    # B needs to see pulse.
    # Pulse is at t=5 in physical time.
    # B records it at t_rec = 5 + true_offset = 7.5.
    # B range 0..8 (approx).
    # So B might miss the pulse if true_offset is 2.5?
    # Wait, B(t_rec) = S(t_rec - 2.5).
    # Peak at t=5.
    # S(t_rec - 2.5) peaks when t_rec - 2.5 = 5 => t_rec = 7.5.
    # If B goes up to 8s (120 frames / 15 Hz = 8s), it barely sees it.
    # Let's start B later or make B longer.
    times_B = times_B + 1.0 # Start B at 1.0s (rec time). Ends at 9.0s.
    # Rec range [1.0, 9.0]. 7.5 is inside.

    data_B = np.zeros((n_frames_B, height, width))

    for i, t_rec in enumerate(times_B):
        t_signal = t_rec - true_offset

        idx = np.searchsorted(t_fine, t_signal)
        if idx >= len(t_fine): idx = len(t_fine) - 1
        if idx < 0: idx = 0

        frame = np.zeros((height, width))
        for m in range(n_modes):
            frame += spatial_modes[m] * temporal_modes[m, idx]

        # Introduce a matrix-based non-linear (quadratic) relationship
        frame_squared = np.square(frame)
        
        # Combine linear and non-linear parts
        frame_combined = frame @ M_linear  + frame_squared @ Mquadra

        # Normalize to prevent extreme values from blowing up the signal
        if np.max(np.abs(frame_combined)) > 0:
            frame_combined /= np.max(np.abs(frame_combined))

        data_B[i] = frame_combined + 0.01 * np.random.randn(height, width)

    # Save files
    hdu_A = fits.PrimaryHDU(data_A)
    hdu_A.writeto('dataA.fits', overwrite=True)
    np.savetxt('dataA.txt', times_A, fmt='%.6f')

    hdu_B = fits.PrimaryHDU(data_B)
    hdu_B.writeto('dataB.fits', overwrite=True)
    np.savetxt('dataB.txt', times_B, fmt='%.6f')

    print(f"Generated data with true offset = {true_offset}")

if __name__ == "__main__":
    create_dummy_data()
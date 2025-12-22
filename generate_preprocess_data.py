import os
import sys
import numpy as np
from astropy.io import fits
import shutil

def generate_stream_data(stream_name, start_time, duration, interval=0.001):
    os.makedirs(stream_name, exist_ok=True)

    # Parameters for FITS
    width = 10
    height = 10
    n_frames_per_file = 10

    current_time = start_time
    end_time = start_time + duration

    global_index = 0

    while current_time < end_time:
        import datetime
        dt = datetime.datetime.fromtimestamp(current_time)
        timestamp_str = dt.strftime("%H:%M:%S.%f") + "000"

        filename_base = f"{stream_name}_{timestamp_str}"
        txt_filename = os.path.join(stream_name, f"{filename_base}.txt")
        fits_filename = os.path.join(stream_name, f"{filename_base}.fits")

        # Create FITS file
        # data[t, y, x] = global_index + t + x/100
        # To make diff meaningful, let's make it change over time
        data = np.zeros((n_frames_per_file, height, width))
        for i in range(n_frames_per_file):
            t_val = (global_index + i) % 10 # Periodic
            data[i, :, :] = t_val

        hdu = fits.PrimaryHDU(data)
        hdu.writeto(fits_filename, overwrite=True)

        with open(txt_filename, 'w') as f:
            f.write("# Header\n")
            for i in range(n_frames_per_file):
                t_acq = current_time + i * interval
                line = f"{i:>10} {global_index:>10} {0.0:>16.9f} {t_acq:>20.9f} {t_acq:>20.6f} {global_index:>10} {0:>10}\n"
                f.write(line)
                global_index += 1

        current_time += n_frames_per_file * interval

def generate_mask(stream_name):
    # Mask center pixel (5,5)
    mask = np.zeros((10, 10))
    mask[5, 5] = 1.0

    filename = f"{stream_name}.mask.fits"
    hdu = fits.PrimaryHDU(mask)
    hdu.writeto(filename, overwrite=True)

if __name__ == "__main__":
    t0 = 1000.0
    if os.path.exists("streamA"): shutil.rmtree("streamA")
    if os.path.exists("streamB"): shutil.rmtree("streamB")

    # Generate A: 10 frames @ 0.1s (10 Hz).
    generate_stream_data("streamA", t0, 1.0, interval=0.1)
    generate_mask("streamA")

    # Generate B: 10 frames @ 0.1s. Same.
    generate_stream_data("streamB", t0, 1.0, interval=0.1)
    generate_mask("streamB")

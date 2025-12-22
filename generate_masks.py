import numpy as np
from astropy.io import fits
import sys

def create_masks():
    width = 50
    height = 50

    # Create maskA: 1s everywhere except a corner
    maskA = np.ones((height, width))
    maskA[0:10, 0:10] = 0.0

    # Create maskB: 1s everywhere
    maskB = np.ones((height, width))

    hdu_A = fits.PrimaryHDU(maskA)
    hdu_A.writeto('maskA.fits', overwrite=True)

    hdu_B = fits.PrimaryHDU(maskB)
    hdu_B.writeto('maskB.fits', overwrite=True)

    print("Generated maskA.fits and maskB.fits")

def create_bad_mask():
    # Wrong dimensions
    mask = np.ones((20, 20))
    hdu = fits.PrimaryHDU(mask)
    hdu.writeto('mask_bad.fits', overwrite=True)
    print("Generated mask_bad.fits")

if __name__ == "__main__":
    create_masks()
    create_bad_mask()

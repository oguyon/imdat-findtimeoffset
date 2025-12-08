# Time Offset Finder for 2D Image Time Series

This program finds the time offset between two time series of 2D images (stored as 3D FITS files) that maximizes the correlation between them. It is written in C for performance, utilizing CFITSIO for file I/O and OpenBLAS/LAPACK for linear algebra operations.

## Problem Description

Given two datasets, `DataA` and `DataB`, each consisting of a sequence of images and associated timestamps:
- `dataA.fits` (3D FITS: Time x Y x X)
- `dataA.txt` (Timestamps for frames in A)
- `dataB.fits` (3D FITS: Time x Y x X)
- `dataB.txt` (Timestamps for frames in B)

The program seeks a time offset $\Delta t$ such that the correlation between $A(t)$ and $B(t + \Delta t)$ is maximized.

## Method

The high-dimensionality of image data makes direct correlation computationally expensive and potentially noisy. This program uses a two-step dimensionality reduction and correlation approach:

1.  **Principal Component Analysis (PCA)**:
    -   Each dataset is reduced to its top $k$ principal components (spatial modes and temporal scores).
    -   This captures the dominant variance in the data while significantly reducing the problem size.
    -   Implemented using `cblas_dgemm` for the Gram matrix and `LAPACKE_dsyev` for eigen-decomposition.

2.  **Canonical Correlation Analysis (CCA)**:
    -   For a candidate offset $\Delta t$, the temporal scores of `DataB` are interpolated to match the time points of `DataA` shifted by $\Delta t$.
    -   CCA finds the linear combinations of the principal components of A and B that are maximally correlated.
    -   Implemented using QR decomposition (`LAPACKE_dgeqrf`/`dorgqr`) and SVD (`LAPACKE_dgesvd`).

3.  **Optimization**:
    -   The program scans a range of time offsets to find the global maximum of the canonical correlation coefficient.

## Prerequisites

-   **CMake** (>= 3.10)
-   **GCC**
-   **CFITSIO** (`libcfitsio-dev`)
-   **OpenBLAS** (`libopenblas-dev`)
-   **LAPACKE** (`liblapacke-dev`)

### Installation (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install cmake gcc libcfitsio-dev libopenblas-dev liblapacke-dev
```

## Build Instructions

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

1.  Ensure your input files are named `dataA.fits`, `dataA.txt`, `dataB.fits`, and `dataB.txt` and are located in the working directory.
    -   **FITS files**: 3D arrays where the 3rd axis is time.
    -   **Text files**: ASCII files with one floating-point timestamp per line, corresponding to the frames in the FITS file.

2.  Run the executable:

```bash
./find_offset
```

The program will output the scanning progress and finally the estimated best time offset and the maximum correlation found.

## Testing

A Python script `generate_data.py` is provided to generate synthetic data with a known time offset for testing purposes.

```bash
# Install dependencies
pip install numpy astropy

# Generate data
python3 generate_data.py

# Run program
./build/find_offset
```

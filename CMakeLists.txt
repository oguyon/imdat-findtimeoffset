cmake_minimum_required(VERSION 3.10)
project(TimeOffsetFinder C)

set(CMAKE_C_STANDARD 99)

# Find CFITSIO
find_path(CFITSIO_INCLUDE_DIR fitsio.h PATHS /usr/local/include /opt/homebrew/include)
find_library(CFITSIO_LIBRARY cfitsio PATHS /usr/local/lib /opt/homebrew/lib)

# Find OpenBLAS
find_path(CBLAS_INCLUDE_DIR cblas.h PATHS /usr/include/x86_64-linux-gnu /usr/include/openblas /usr/include /usr/local/include /opt/homebrew/include)
find_library(OPENBLAS_LIBRARY openblas PATHS /usr/local/lib /opt/homebrew/lib)

# Find LAPACKE
find_path(LAPACKE_INCLUDE_DIR lapacke.h PATHS /usr/local/include /opt/homebrew/include)
find_library(LAPACKE_LIBRARY lapacke PATHS /usr/local/lib /opt/homebrew/lib)
find_library(LAPACK_LIBRARY lapack)

include_directories(${CFITSIO_INCLUDE_DIR} ${CBLAS_INCLUDE_DIR} ${LAPACKE_INCLUDE_DIR})

add_executable(find_offset main.c)
target_link_libraries(find_offset ${CFITSIO_LIBRARY} ${OPENBLAS_LIBRARY} ${LAPACKE_LIBRARY} ${LAPACK_LIBRARY} m)

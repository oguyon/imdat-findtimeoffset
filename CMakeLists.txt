cmake_minimum_required(VERSION 3.10)
project(TimeOffsetFinder C)

set(CMAKE_C_STANDARD 99)

find_package(CFITSIO REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(find_offset main.c)

# CFITSIO
if (TARGET CFITSIO::CFITSIO)
    target_link_libraries(find_offset PRIVATE CFITSIO::CFITSIO)
else()
    target_include_directories(find_offset PRIVATE ${CFITSIO_INCLUDE_DIRS})
    target_link_libraries(find_offset PRIVATE ${CFITSIO_LIBRARIES})
endif()

# Link LAPACK + LAPACKE
target_link_libraries(find_offset PRIVATE
    ${LAPACK_LIBRARIES}                   # For LAPACK & LAPACKE
    /usr/local/opt/lapack/lib/liblapacke.dylib  # Explicit LAPACKE
)
# Includes for lapacke.h and cblas.h
target_include_directories(find_offset
    PRIVATE
    ${LAPACK_INCLUDE_DIRS}
    /usr/local/opt/openblas/include
)

# Link CFITSIO
find_library(CFITSIO_LIBRARY cfitsio PATHS /usr/local/opt/cfitsio/lib)
target_link_libraries(find_offset PRIVATE ${CFITSIO_LIBRARY})

# Add math library if needed
target_link_libraries(find_offset PRIVATE m)


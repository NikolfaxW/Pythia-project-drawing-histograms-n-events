cmake_minimum_required(VERSION 3.10)
project(Pythia-project-drawing-histograms)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# PYTHIA configuration
set(PYTHIA_PREFIX "/home/nikolfax/pythia8311")
set(PYTHIA_INCLUDE ${PYTHIA_PREFIX}/include)
set(PYTHIA_LIB ${PYTHIA_PREFIX}/lib)

# LHAPDF configuration
set(LHAPDF_PREFIX "/home/nikolfax/LHPDF")
set(LHAPDF_INCLUDE ${LHAPDF_PREFIX}/include)
set(LHAPDF_LIB ${LHAPDF_PREFIX}/lib)

# FastJet configuration
set(FASTJET_PREFIX "/home/nikolfax/fastjet-install")
set(FASTJET_INCLUDE ${FASTJET_PREFIX}/include)
set(FASTJET_LIB ${FASTJET_PREFIX}/lib)

# Execute root-config to get ROOT configuration
execute_process(COMMAND /snap/root-framework/931/usr/local/bin/root-config --cflags
        OUTPUT_VARIABLE ROOT_CFLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND /snap/root-framework/931/usr/local/bin/root-config --libs
        OUTPUT_VARIABLE ROOT_LIBS
        OUTPUT_STRIP_TRAILING_WHITESPACE)

# Include directories
include_directories(
        ${PYTHIA_INCLUDE}
        ${LHAPDF_INCLUDE}
        ${FASTJET_INCLUDE}
        ${ROOT_INCLUDE_DIRS} # Use ROOT_INCLUDE_DIRS obtained from root-config
)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CFLAGS} -O2 -std=c++17 -pedantic -W -Wall -Wshadow -fPIC -pthread")

# Add source files to your project
add_executable(Pythia-project-drawing-histograms drawF.h drawF.cpp main.cpp
)

# Link the necessary libraries
target_link_libraries(Pythia-project-drawing-histograms PRIVATE
        ${PYTHIA_LIB}/libpythia8.so
        ${LHAPDF_LIB}/libLHAPDF.so
        ${FASTJET_LIB}/libfastjet.so
        ${ROOT_LIBS} # Use ROOT_LIBS obtained from root-config
        # Add other ROOT libraries here as needed
)
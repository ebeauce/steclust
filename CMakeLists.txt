cmake_minimum_required(VERSION 3.24)
project(steclust LANGUAGES C CXX)

# 1. Initialize Python Extension Build Support
find_package(Python REQUIRED COMPONENTS Development.Module)

# 2. Handle OpenMP Configuration Natively across Linux and macOS
find_package(OpenMP REQUIRED)

# -----------------------------------------------
#        CPU Beamforming Engine
# -----------------------------------------------
python_add_library(space_time_energy_cpu MODULE steclust/src/space_time_energy.c)

target_compile_options(space_time_energy_cpu PRIVATE 
    -O3 
    -fPIC 
    -ftree-vectorize 
    -march=native
)

# Connect OpenMP variables natively found by CMake (handles Homebrew/venvs/Linux automatically)
target_link_libraries(space_time_energy_cpu PRIVATE OpenMP::OpenMP_C)

# Install into the output layout inside the destination module's lib directory
install(TARGETS space_time_energy_cpu DESTINATION steclust/lib)

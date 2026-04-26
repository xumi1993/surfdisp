cmake_minimum_required(VERSION 3.14)
project(surfdisp LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -O2 for both to get comparable optimisation levels
set(CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS}     -O2")
# gfortran: allow lines longer than 72 chars in fixed-form source (.f)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -ffixed-line-length-none")

include_directories(${PROJECT_SOURCE_DIR}/include)

# ---- Fortran library (original surfdisp96 + C wrapper) ----
add_library(surfdisp_f STATIC
    src/surfdisp96.f
    src/surfdisp_fwrapper.f90
)

# ---- C++ library ----
add_library(surfdisp_cpp STATIC
    src/surfdisp.cpp
)
target_include_directories(surfdisp_cpp PUBLIC include)

# ---- Test / comparison executable ----
add_executable(test_compare test/test_compare.cpp)
target_include_directories(test_compare PRIVATE include)
target_link_libraries(test_compare surfdisp_cpp surfdisp_f)

# ---- Benchmark executable ----
add_executable(benchmark test/benchmark.cpp)
target_include_directories(benchmark PRIVATE include)
target_link_libraries(benchmark surfdisp_cpp surfdisp_f)

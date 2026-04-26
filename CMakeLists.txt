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

# ---- constants Fortran module (shared by slegn and sregn) ----
add_library(constants_f STATIC src/constants.f90)

# ---- slegn C++ library ----
add_library(slegn_cpp STATIC src/slegn.cpp)
target_include_directories(slegn_cpp PUBLIC include)

# ---- slegn Fortran library ----
add_library(slegn_f STATIC src/slegn96.f90)
target_link_libraries(slegn_f constants_f)

# ---- slegn comparison test ----
add_executable(test_slegn test/test_slegn.cpp)
target_include_directories(test_slegn PRIVATE include)
target_link_libraries(test_slegn slegn_cpp slegn_f constants_f surfdisp_cpp surfdisp_f)

# ---- slegn benchmark ----
add_executable(benchmark_slegn test/benchmark_slegn.cpp)
target_include_directories(benchmark_slegn PRIVATE include)
target_link_libraries(benchmark_slegn slegn_cpp slegn_f constants_f surfdisp_cpp surfdisp_f)

# ---- sregn C++ library ----
add_library(sregn_cpp STATIC src/sregn.cpp)
target_include_directories(sregn_cpp PUBLIC include)

# ---- sregn Fortran library ----
add_library(sregn_f STATIC src/sregn96.f90 src/sregn96_fwrapper.f90)
target_link_libraries(sregn_f constants_f)

# ---- sregn comparison test ----
add_executable(test_sregn test/test_sregn.cpp)
target_include_directories(test_sregn PRIVATE include)
target_link_libraries(test_sregn sregn_cpp sregn_f constants_f surfdisp_cpp surfdisp_f)

# ---- sregn benchmark ----
add_executable(benchmark_sregn test/benchmark_sregn.cpp)
target_include_directories(benchmark_sregn PRIVATE include)
target_link_libraries(benchmark_sregn sregn_cpp sregn_f constants_f surfdisp_cpp surfdisp_f)

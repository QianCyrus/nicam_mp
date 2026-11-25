cmake_minimum_required(VERSION 3.20)
project(NICAM_NSW6_KOKKOS_BENCH LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Kokkos REQUIRED)

add_executable(bench_kokkos_full bench_kokkos_full.cpp)
target_include_directories(bench_kokkos_full PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../nicam_mp/cuda_kernels)
target_link_libraries(bench_kokkos_full PRIVATE Kokkos::kokkos)

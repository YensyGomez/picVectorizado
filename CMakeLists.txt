cmake_minimum_required(VERSION 2.8)
project(PIC2D_CPU)
find_package(FFTW REQUIRED)
add_executable(PIC2D_CPU pic2dAplanado.cpp)
target_link_libraries(PIC2D_CPU m fftw3)

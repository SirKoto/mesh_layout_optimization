cmake_minimum_required(VERSION 3.14)

add_library(eigen INTERFACE)

target_include_directories(eigen INTERFACE eigen/)

add_library(spectra INTERFACE)

target_include_directories(spectra INTERFACE spectra/include)


# tinyply
add_library(tinyply STATIC tinyply/source/tinyply.cpp)
target_include_directories(tinyply PUBLIC tinyply/source/)
target_compile_features(tinyply PRIVATE cxx_std_11)

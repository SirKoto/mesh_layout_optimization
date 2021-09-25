cmake_minimum_required(VERSION 3.20)

add_library(eigen INTERFACE)

target_include_directories(eigen INTERFACE eigen/)
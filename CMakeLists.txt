cmake_minimum_required(VERSION 3.20)

project(mesh_layout_opt
    VERSION 0.0.1
    LANGUAGES CXX C)

# Libraries
add_subdirectory(libs)

# Builds executable
add_subdirectory(src)


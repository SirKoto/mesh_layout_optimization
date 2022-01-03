cmake_minimum_required(VERSION 3.14)


add_executable(${PROJECT_NAME}
    main.cpp
    TriangleMesh.cpp TriangleMesh.hpp
    Args.cpp Args.hpp 
    LayoutMaker.cpp LayoutMaker.hpp
    LayoutOptimizer.cpp LayoutOptimizer.hpp)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_link_libraries(${PROJECT_NAME} PRIVATE eigen tinyply spectra)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()
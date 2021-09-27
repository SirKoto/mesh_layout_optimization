cmake_minimum_required(VERSION 3.14)


add_executable(${PROJECT_NAME}
    main.cpp
    TriangleMesh.cpp TriangleMesh.hpp
    Args.cpp Args.hpp)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_link_libraries(${PROJECT_NAME} PRIVATE eigen tinyply)
cmake_minimum_required(VERSION 3.20)


add_executable(${PROJECT_NAME}
    main.cpp)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_link_libraries(${PROJECT_NAME} PRIVATE eigen)
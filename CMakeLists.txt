cmake_minimum_required(VERSION 3.10)
project(Cloud)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(
    Cloud
    main.cpp

    Grid_Construction/vector3D.cpp

    Surface/ply.c
    Surface/MarchingCubes.cpp

)


target_include_directories(
    Cloud
    PUBLIC ${CMAKE_SOURCE_DIR}
)
cmake_minimum_required(VERSION 3.10)
project(Cloud)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(
    Cloud
    main.cpp
    Advection/advection.cpp

    Force/General_Force/Buoyancy.cpp
    Force/General_Force/GravityForce.cpp

    Grid_Construction/vector3D.cpp

    Solver/solver.cpp
    
    Surface/MarchingCubes.cpp
    Surface/ply.c


)


target_include_directories(
    Cloud
    PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/Surface
    ${CMAKE_SOURCE_DIR}/NOISE
    ${CMAKE_SOURCE_DIR}/Grid_Construction
)
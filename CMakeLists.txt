cmake_minimum_required(VERSION 3.10)
project(sph)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set (src src/sph.cpp main.cpp)

add_executable(sph ${src})
target_include_directories(sph PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include/>
)
target_link_libraries(sph Eigen3::Eigen GL GLU glut)
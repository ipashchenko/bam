cmake_minimum_required(VERSION 3.11)
project(bam)

include_directories("/usr/inlcude/eigen3", "/home/ilya/github/DNest4/code", "src")


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++17 -O3 -DNDEBUG -w -DEIGEN_MPL2_ONLY")

set(SOURCE_FILES src/main.cpp src/Component.cpp src/Data.cpp src/SkyModel.cpp src/DNestModel.cpp src/utils.h)
add_executable(bam ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(bam /home/ilya/github/DNest4/code/libdnest4.a)
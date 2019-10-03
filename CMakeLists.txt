cmake_minimum_required(VERSION 3.9)
project(bam)

include_directories("/usr/include/eigen3", "/home/ilya/github/DNest4/code", "src")


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11 -O3 -DNDEBUG -w -DEIGEN_MPL2_ONLY")

set(SOURCE_FILES src/main.cpp src/Component.cpp src/Data.cpp src/MyConditionalPrior.h src/MyConditionalPrior.cpp src/DNestModel.cpp)
add_executable(bam ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(bam /home/ilya/github/DNest4/code/libdnest4.a)
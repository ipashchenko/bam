cmake_minimum_required(VERSION 3.11)
project(bam)

include_directories("/usr/include/eigen3", "/home/ilya/github/DNest4/code", "src")


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11 -march=native -O3 -DNDEBUG -w -DEIGEN_MPL2_ONLY")

include(FetchContent)
FetchContent_Declare(
        rapidcsv
        GIT_REPOSITORY "https://github.com/d99kris/rapidcsv.git"
        GIT_TAG        "v8.80"
)
FetchContent_MakeAvailable(rapidcsv)

set(SOURCE_FILES src/main.cpp src/Data.cpp src/MyConditionalPrior.cpp src/DNestModel.cpp src/LogNormal.h src/LogNormal.cpp src/Gaussian2D.h src/Gaussian2D.cpp src/Fixed.h src/Fixed.cpp src/Utils.h)
add_executable(bam ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(bam /home/ilya/github/DNest4/code/libdnest4.a)
TARGET_LINK_LIBRARIES(bam rapidcsv)

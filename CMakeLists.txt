cmake_minimum_required(VERSION 3.23)

# set(CMAKE_C_COMPILER "D:/ProgramData/MinGW/mingw64/bin/gcc.exe")
# set(CMAKE_CXX_COMPILER "D:/ProgramData/MinGW/mingw64/bin/g++.exe")

project(SINGLEX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Debug)
# You have to set these variables as Windows environment variables:
# EIGEN3_INCLUDE_DIR <- %EIGEN3_ROOT%
#
SET( EIGEN3_INCLUDE_DIR "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3" )
SET( TCLAP_INCLUDE_DIR "/opt/cLibs/tclap-1.2.2/include" )

find_package( Eigen3 3.3 REQUIRED )
# include_directories is needed for the compiler to know where looking for Eigen3 header files to be included 
include_directories( ${EIGEN3_INCLUDE_DIR} ${TCLAP_INCLUDE_DIR} )

add_executable(SXCpp main.cpp Slip.cpp Grain.cpp disvel.cpp func.cpp global.cpp load.cpp read_input.cpp write_output.cpp singleX.h)

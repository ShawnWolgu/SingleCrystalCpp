cmake_minimum_required(VERSION 3.23)
project(SINGLEX)
set(CMAKE_CXX_STANDARD 17)

# You have to set these variables as Windows environment variables:
# EIGEN3_INCLUDE_DIR <- %EIGEN3_ROOT%
# EIGEN3_DIR <- %EIGEN3_ROOT%\cmake
#
# EIGEN3_INCLUDE_DIR: variable needed for file %EIGEN3_ROOT%/cmake/FindEigen3.cmake
#
# CMAKE_MODULE_PATH: Search path for the module Eigen3 to be loaded by find_package
#
SET( EIGEN3_INCLUDE_DIR "E:/eigen-3.4.0/" )
SET( EIGEN3U_INCLUDE_DIR "E:/eigen-3.4.0/unsupported/" )
SET( CMAKE_MODULE_PATH "E:/eigen-3.4.0/cmake" )

find_package( Eigen3 3.3 REQUIRED )
# include_directories is needed for the compiler to know where looking for Eigen3 header files to be included 
include_directories( ${EIGEN3_INCLUDE_DIR} )
include_directories( ${EIGEN3U_INCLUDE_DIR} )

add_executable(cmake_app main.cpp Slip.cpp Grain.cpp func.cpp)
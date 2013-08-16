cmake_minimum_required(VERSION 2.6)
set(OpenCV_DIR "/usr/local/share/OpenCV")
project (MTC)

find_package(OpenCV REQUIRED core highgui imgproc gpufilters gpuarithm )
#include("/usr/local/include")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -g3 -Wall  ")
set(CMAKE_CXX_FLAGS_RELEASE   "-std=c++11 -Ofast -Wall ")

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo
MinSizeRel."
       FORCE )
endif()



add_subdirectory(modules)
add_subdirectory(data)



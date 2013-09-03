cmake_minimum_required(VERSION 2.8)
set(OpenCV_DIR "/usr/local/share/OpenCV")
project (MTC)

find_package(OpenCV REQUIRED core highgui imgproc gpufilters gpuarithm gpuimgproc ml)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.54.0 COMPONENTS regex)
if(Boost_FOUND)
message("${Boost_LIBRARIES}")
endif()
find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
message("find eigen3" )
message("${EIGEN3_INCLUDE_DIR}" )
endif()

#include("/usr/local/include")
#set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -g3 -Wall  ")
set(CMAKE_CXX_FLAGS_RELEASE   "-std=c++11 -Ofast -Wall ")

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo
MinSizeRel."
       FORCE )
endif()



add_subdirectory(modules)
#add_subdirectory(data)




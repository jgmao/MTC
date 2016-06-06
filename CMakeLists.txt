cmake_minimum_required(VERSION 2.8)
project (MTC)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(OS "MacOSX")
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(OS "Linux")
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "WINDOWS")
  set(OS "WIN32")
endif()
message("${OS}")

cmake_policy(SET CMP0017 OLD)

if ("${OS}" MATCHES "Linux")
set(OpenCV_DIR "/usr/local/share/OpenCV")
else()
set(OpenCV_DIR "/opt/local/lib/cmake")
endif()

message (GCC_MAJOR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake-2.8/Modules")
message("${CMAKE_MODULE_PATH}")

find_package(OpenCV COMPONENTS core)
if(${OpenCV_VERSION_MINOR} LESS 9)
find_package(OpenCV REQUIRED core highgui imgproc ml)# gpufilters gpuarithm gpuimgproc ml)
else()
find_package(OpenCV REQUIRED core highgui imgproc ml features2d nonfree)
endif()

message("${OpenCV_VERSION}")

add_definitions(-DBOOST_ALL_NO_LIB)
add_definitions( -DBOOST_ALL_DYN_LINK )
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.54.0 COMPONENTS regex system filesystem)
if(Boost_FOUND)
message("${Boost_LIBRARIES}")
endif()

find_package(Eigen3 REQUIRED)
message("find eigen3" )
message("${EIGEN3_INCLUDE_DIR}" )

if("${OS}" MATCHES "MacOSX")
set(CMAKE_CXX_COMPILER "/opt/local/bin/g++")
else()
set(CMAKE_CXX_COMPILER "g++")
endif()
message("${CMAKE_CXX_COMPILER}")

set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -lm -lfftw3 -lfftw3_omp -fopenmp -ffast-math -msse2 -O0 -g3 -Wall -Wno-unused-variable -Wno-unused-but-set-variable -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE   "-std=c++11 -lm -lfftw3 -lfftw3_omp -fopenmp -ffast-math -msse2 -O3  -DNDEBUG -Wall -Wno-unused-variable -Wno-unused-but-set-variable -fPIC")

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo
MinSizeRel."
       FORCE )
endif()



add_subdirectory(modules)
#add_subdirectory(data)




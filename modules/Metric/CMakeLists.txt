#set(MTC_TENSOR_DIR ${MTC_MODULES_DIR}/TensorLib)
#set(MTC_METRIC_DIR ${MTC_MODULES_DIR}/Metric)
#set(MTC_UTIL_DIR ${MTC_MODULERS_DIR}/Utility)

include_directories(${MTC_TENSORLIB_DIR}/include ${MTC_UTILITY_DIR}/include ${MTC_METRIC_DIR}/include ${EIGEN3_PATH})

#link_directories(${MTC_TENSOR_DIR}/src)
set(common_link_dir  ${MTC_TENSORLIB_DIR} ${MTC_METRIC_DIR}/src ${MTC_UTILITY_DIR})
set(PROGRAMS Steerable fftw++ Steerable2 LRI algorithms IO Granulate Metric)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Debug/bin)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  link_directories(${CMAKE_SOURCE_DIR}/Debug/lib ${common_link_dir})
else()
  set(CMAKE_LIBRARY_:OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Release/lib)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  link_directories(${CMAKE_SOURCE_DIR}/Release/lib ${common_link_dir})
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Release/bin)
endif()

message("${ARCHIVE_OUTPUT_DIRECTORY}")


foreach(program ${PROGRAMS})
  #add_library(${program} src/${program}.cpp)
  #target_link_libraries(${program} ${OpenCV_LIBS} ${Boost_LIBRARIES} TensorLib util)
  set(all_src ${all_src} src/${program}.cpp)
  #install(TARGETS ${PROGRAMS}
#	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
#	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
#	ARCHIVE DESTINATION ${MTC_SOURCE_DIR}/lib
# )
endforeach(program)

add_library(MetricLib ${all_src}) #cannot be shared libary, boost complain
target_link_libraries(MetricLib ${OpenCV_LIBS} util ${Boost_LIBRARIES} TensorLib fftw3 fftw3_omp)





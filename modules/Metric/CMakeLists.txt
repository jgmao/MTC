set(MTC_TENSOR_DIR ${MTC_MODULES_DIR}/TensorLib)
set(MTC_METRIC_DIR ${MTC_MODULES_DIR}/Metric)
set(MTC_UTIL_DIR ${MTC_MODULERS_DIR}/Utility)

include_directories(${MTC_TENSOR_DIR}/include ${MTC_MODULES_DIR}/Utility/include ${MTC_METRIC_DIR}/include)

#link_directories(${MTC_TENSOR_DIR}/src)
set(common_link_dir  ${MTC_TENSOR_DIR} ${MTC_METRIC_DIR}/src ${MTC_UTIL_DIR})
set(PROGRAMS Steerable LRI algorithms Granulate)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Debug/bin)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  link_directories(${CMAKE_SOURCE_DIR}/Debug/lib ${common_link_dir})
else()
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Release/lib)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  link_directories(${CMAKE_SOURCE_DIR}/Release/lib ${common_link_dir})
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Release/bin)
endif()

message("${ARCHIVE_OUTPUT_DIRECTORY}")


foreach(program ${PROGRAMS})
  add_library(${program} src/${program}.cpp)
  target_link_libraries(${program} ${OpenCV_LIBS} TensorLite util)
  #install(TARGETS ${PROGRAMS}
#	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
#	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
#	ARCHIVE DESTINATION ${MTC_SOURCE_DIR}/lib
# )
endforeach(program)



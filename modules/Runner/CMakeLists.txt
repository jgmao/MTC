set(MTC_RUNNER_DIR ${MTC_MODULES_DIR}/Runner)
set(MTC_TENSOR_DIR ${MTC_MODULES_DIR}/TensorLib)
set(MTC_UTIL_DIR ${MTC_MODULES_DIR}/Utility)
set(MTC_METRIC_DIR ${MTC_MODULES_DIR}/Metric)

message("tensor dir is: ${MTC_TENSOR_DIR}")
include_directories(${MTC_RUNNER_DIR}/include ${MTC_TENSOR_DIR}/include ${MTC_UTIL_DIR}/include ${MTC_METRIC_DIR}/include)
set(com_link ${MTC_RUNNER_DIR}/src ${MTC_TENSOR_DIR} ${MTC_METRIC_DIR} ${MTC_UTIL_DIR})


set(PROGRAMS gpussim testTensor main)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Debug/bin)
  link_directories(${CMAKE_SOURCE_DIR}/Debug/lib ${com_link})
else()
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Release/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Release/bin)
        link_directories(${CMAKE_SOURCE_DIR}/Release/lib ${com_link})
endif()
 
add_library(Tester src/Tester.cpp)
target_link_libraries(Tester ${OpenCV_LIBS} TensorLite Size3 Cube Steerable algorithms util LRI Granulate)

foreach(program ${PROGRAMS})
  add_executable(${program} src/${program}.cpp)
  target_link_libraries(${program} ${OpenCV_LIBS} Cube Size3 Tester TensorLite Steerable LRI algorithms util LRI Granulate)
 # install(TARGETS ${PROGRAMS}
#	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
#	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
#  )
endforeach(program)



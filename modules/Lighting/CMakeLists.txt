#set(MTC_TENSOR_DIR ${MTC_MODULES_DIR}/TensorLib)
#set(MTC_UTIL_DIR ${MTC_MODULES_DIR}/Utility)
#set(MTC_LT_DIR ${MTC_MODULES_DIR}/Lighting)
#set(MTC_METRIC_DIR
include_directories(${MTC_TENSORLIB_DIR}/include ${MTC_UTILITY_DIR}/include ${MTC_LIGHTING_DIR}/include ${EIGEN3_INCLUDE_DIR} ${MTC_METRIC_DIR}/include)
#link_directories(${MTC_TENSOR_DIR}/src)
#message("${EIGEN3_INCLUDE_DIR}")
#set(common_link_dir  ${MTC_TENSOR_DIR} ${MTC_LT_DIR}/src ${MTC_UTIL_DIR})
set(PROGRAMS ThinPlateSpline PoissonSolver LightingVerify Lighting)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Debug/bin)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  link_directories(${CMAKE_SOURCE_DIR}/Debug/lib)
else()
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Release/lib)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  link_directories(${CMAKE_SOURCE_DIR}/Release/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Release/bin)
endif()

#message("${ARCHIVE_OUTPUT_DIRECTORY}")


foreach(program ${PROGRAMS})
  #add_library(${program} src/${program}.cpp)
  #target_link_libraries(${program} ${OpenCV_LIBS} ${Boost_LIBRARIES} TensorLib util)
  #install(TARGETS ${PROGRAMS}
  set(all_src ${all_src} src/${program}.cpp)
#	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
#	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
#	ARCHIVE DESTINATION ${MTC_SOURCE_DIR}/lib
# )
endforeach(program)
add_library(LightingLib SHARED ${all_src})
target_link_libraries(LightingLib ${OpenCV_LIBS} util ${Boost_LIBRARIES} TensorLib)

include_directories(${MTC_TENSORLIB_DIR}/include ${MTC_UTILITY_DIR}/include ${MTC_MTC_DIR}/include ${MTC_LIGHTING_DIR}/include ${MTC_METRIC_DIR}/include)

#link_directories(${MTC_TENSOR_DIR}/src)
set(common_link_dir  ${MTC_TENSORLIB_DIR} ${MTC_MTC_DIR}/src ${MTC_UTILITY_DIR})
set(PROGRAMS Statistics HuffTree MTC)

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
#  add_library(${program} src/${program}.cpp)
  set(all_src ${all_src} src/${program}.cpp)
#  target_link_libraries(${program} ${OpenCV_LIBS} ${Boost_LIBRARIES} TensorLib TensorExtend MetricLib LightingLib util)
  #install(TARGETS ${PROGRAMS}
#	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
#	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
#	ARCHIVE DESTINATION ${MTC_SOURCE_DIR}/lib
# )
endforeach(program)


add_library(MTC ${all_src}) #cannot be shared libary, boost complain
target_link_libraries(MTC ${OpenCV_LIBS} util ${Boost_LIBRARIES} TensorLib TensorExtend MetricLib LightingLib jpeg)

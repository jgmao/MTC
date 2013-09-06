#set(MTC_TENSOR_DIR ${MTC_MODULES_DIR}/TensorLib)
include_directories(${MTC_TENSORLIB_DIR}/include ${MTC_UTILITY_DIR}/include ${MTC_LIGHTING_DIR}/include ${MTC_METRIC_DIR}/include)
#link_directories(${MTC_TENSOR_DIR}/src)
set(PROGRAMS Size3 Cube TensorLite)
set(Other CubePlus QNode QTree QGrid)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Debug/bin)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  link_directories(${CMAKE_SOURCE_DIR}/Debug/lib ${MTC_TENSORLIB_DIR}/src)
else()
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Release/lib)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  link_directories(${CMAKE_SOURCE_DIR}/Release/lib ${MTC_TENSORLIB_DIR}/src)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Release/bin)
endif()


set (all_src "")
message("${ARCHIVE_OUTPUT_DIRECTORY}")


foreach(program ${PROGRAMS})
  add_library(${program} src/${program}.cpp)
  target_link_libraries(${program} ${OpenCV_LIBS})
  set(all_src ${all_src} src/${program}.cpp)
  #install(TARGETS ${PROGRAMS}
#	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
#	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
#	ARCHIVE DESTINATION ${MTC_SOURCE_DIR}/lib
# )
endforeach(program)


#build to a single one
add_library(TensorLib SHARED ${all_src})
target_link_libraries(TensorLib ${OpenCV_LIBS} util)

foreach(program ${Other})
  add_library(${program} src/${program}.cpp)
  target_link_libraries(${program} ${OpenCV_LIBS} Lighting)
  #set(all_src ${all_src} src/${program}.cpp)
  #install(TARGETS ${PROGRAMS}
#	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
#	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
#	ARCHIVE DESTINATION ${MTC_SOURCE_DIR}/lib
# )
endforeach(program)


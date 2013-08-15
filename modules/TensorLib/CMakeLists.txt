set(MTC_TENSOR_DIR ${MTC_MODULES_DIR}/TensorLib)
include_directories(${MTC_TENSOR_DIR}/include)
link_directories(${MTC_TENSOR_DIR}/src)
set(PROGRAMS Size3 Cube Tensor)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Debug/bin)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
else()
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Release/lib)
  set(ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Release/bin)
endif()

message("${ARCHIVE_OUTPUT_DIRECTORY}")
foreach(program ${PROGRAMS})
  add_library(${program} src/${program}.cpp)
  target_link_libraries(${program} ${OpenCV_LIBS})
  #install(TARGETS ${PROGRAMS}
#	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
#	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
#	ARCHIVE DESTINATION ${MTC_SOURCE_DIR}/lib
# )
endforeach(program)



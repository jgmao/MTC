set(MTC_RUNNER_DIR ${MTC_MODULES_DIR}/Runner)
include_directories(${MTC_RUNNER_DIR}/include)
link_directories(${MTC_RUNNER_DIR}/src )

set(PROGRAMS gpussim)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Debug/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Debug/bin)
else()
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Release/lib)
  set(LIB_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Release/lib)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/Release/bin)
endif()

foreach(program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${OpenCV_LIBS})
  install(TARGETS ${PROGRAMS}
	RUNTIME DESTINATION ${MTC_SOURCE_DIR}/bin
	LIBRARY DESTINATION ${MTC_SOURCE_DIR}/lib
  )
endforeach(program)



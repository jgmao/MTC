
set(MTC_MODULES_DIR ${MTC_SOURCE_DIR}/modules)
set(MTC_DATA_DIR ${MTC_DATA_DIR}/data)
message("${MTC_MODULES_DIR}")
#add modules
add_subdirectory(Runner)
add_subdirectory(TensorLib)
add_subdirectory(Utility)

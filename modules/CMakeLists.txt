
set(MTC_MODULES_DIR ${MTC_SOURCE_DIR}/modules)
set(MTC_DATA_DIR ${MTC_DATA_DIR}/data)
message("${MTC_MODULES_DIR}")

set(MTC_TENSORLIB_DIR ${MTC_MODULES_DIR}/TensorLib)
set(MTC_UTILITY_DIR ${MTC_MODULES_DIR}/Utility)
set(MTC_METRIC_DIR ${MTC_MODULES_DIR}/Metric)
set(MTC_RUNNER_DIR ${MTC_MODULES_DIR}/Runner)
set(MTC_LIGHTING_DIR ${MTC_MODULES_DIR}/Lighting)


#add modules


add_subdirectory(TensorLib)
add_subdirectory(Utility)
add_subdirectory(Metric)
add_subdirectory(Runner)
add_subdirectory(Lighting)

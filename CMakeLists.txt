cmake_minimum_required(VERSION 3.1)
project(pytetgen)

option(BUILD_EXECUTABLE "Build tetgen executable" ON)
option(BUILD_LIBRARY    "Build libtetgen library" ON)

if(BUILD_LIBRARY)
  add_library(tetgen STATIC tetgen.cxx predicates.cxx)
  target_compile_definitions(tetgen PRIVATE -DTETLIBRARY)
  # Generate position independent code
  set_target_properties(tetgen PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()


if(BUILD_EXECUTABLE)
  add_executable(tetgen_exec tetgen.cxx predicates.cxx)
  set_target_properties(tetgen_exec PROPERTIES OUTPUT_NAME tetgen)
endif()

# ================================== pytetgen

add_subdirectory(extern/pybind11)
pybind11_add_module(pytetgen src/pytetgen.cpp tetgen.cxx predicates.cxx)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(pytetgen
                            PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
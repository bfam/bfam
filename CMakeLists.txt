cmake_minimum_required(VERSION 2.8.7)

project(bfam C Fortran CXX)

include(InstallRequiredSystemLibraries)
include(FortranCInterface)
include(CMakePushCheckState)

find_package(Git)
find_package(MPI REQUIRED)

# Directory where customize cmake files reside
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/findopencl
  ${CMAKE_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}
)

find_package(CUDA)
find_package(OpenCL)
find_package(OpenMP)

##
## Version Stuff
##
if(GIT_FOUND)
  set(USE_GIT 1)
else(GIT_FOUND)
  set(USE_GIT 0)
endif(GIT_FOUND)
file(WRITE ${CMAKE_BINARY_DIR}/Version.cmake
"if(${USE_GIT} EQUAL 1)
  execute_process(
       COMMAND ${GIT_EXECUTABLE} describe --dirty ${GIT_DESCRIBE_FLAGS}
       WORKING_DIRECTORY \${WORKING_DIRECTORY}
       OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
       RESULT_VARIABLE GIT_DESCRIBE_RESULT
       ERROR_VARIABLE GIT_DESCRIBE_ERROR
       OUTPUT_STRIP_TRAILING_WHITESPACE
   )
else(${USE_GIT} EQUAL 1)
  set(GIT_DESCRIBE_RESULT -1)
endif(${USE_GIT} EQUAL 1)

if(GIT_DESCRIBE_RESULT EQUAL 0)
  file(WRITE ${CMAKE_BINARY_DIR}/source/VERSION \"\${GIT_DESCRIBE_VERSION}\")
else(GIT_DESCRIBE_RESULT EQUAL 0)
  file(READ  ${CMAKE_SOURCE_DIR}/VERSION GIT_DESCRIBE_VERSION)
endif(GIT_DESCRIBE_RESULT EQUAL 0)

string(REGEX REPLACE \"v([0-9]*)\\\\.([0-9]*)\\\\.(.*)\"
       \"\\\\1.\\\\2.\\\\3\" VERSION \"\${GIT_DESCRIBE_VERSION}\" )

string(REGEX REPLACE \"v([0-9]*)\\\\.([0-9]*)\\\\.(.*)\"
       \"\\\\1\" VERSION_MAJOR \"\${GIT_DESCRIBE_VERSION}\")
string(REGEX REPLACE \"v([0-9]*)\\\\.([0-9]*)\\\\.(.*)\"
       \"\\\\2\" VERSION_MINOR \"\${GIT_DESCRIBE_VERSION}\")
string(REGEX REPLACE \"v([0-9]*)\\\\.([0-9]*)\\\\.(.*)\"
       \"\\\\3\" VERSION_PATCH \"\${GIT_DESCRIBE_VERSION}\")
")

file(WRITE ${CMAKE_BINARY_DIR}/FillVersion.cmake
"INCLUDE(Version.cmake)
CONFIGURE_FILE(\${SRC} \${DST} @ONLY)
")


file(WRITE ${CMAKE_BINARY_DIR}/bfam_version_info.h.in
  "\#ifndef BFAM_VERSION_INFO_H\n"
  "\#define BFAM_VERSION_INFO_H\n"
  "\#define BFAM_VERSION_NUMBER \"@VERSION@\"\n"
  "\#endif\n"
)
add_custom_target(
    version
    ${CMAKE_COMMAND} -D SRC=${CMAKE_BINARY_DIR}/bfam_version_info.h.in
                     -D DST=${CMAKE_BINARY_DIR}/bfam_version_info.h
                     -D WORKING_DIRECTORY=${CMAKE_SOURCE_DIR}
                     -P ${CMAKE_BINARY_DIR}/FillVersion.cmake
)

set(WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
include(${CMAKE_BINARY_DIR}/Version.cmake)
set(bfam_VERSION "${VERSION}")
set(bfam_VERSION_MAJOR "${VERSION_MAJOR}")
set(bfam_VERSION_MINOR "${VERSION_MINOR}")
set(bfam_VERSION_PATCH "${VERSION_PATCH}")

if(UNIX)
  file(WRITE ${CMAKE_BINARY_DIR}/Dist.cmake
"message(STATUS \"Package Dist running...\")
include(CPackSourceConfig.cmake)
execute_process(COMMAND cpack -G TGZ --config CPackSourceConfig.cmake
  TIMEOUT 3600
  WORKING_DIRECTORY \${CMAKE_BINARY_DIR})
")
  #
  # Add custom target
  #
  add_custom_target(dist
    COMMAND ${CMAKE_COMMAND} .
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/Dist.cmake
    )
  add_dependencies(dist version)
endif(UNIX)

##
## Ctags
##
add_custom_target (tags
  COMMAND ctags -R --fields=+S src
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )

##
## Customize CMake defaults
##
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BUILD_DIR}/local" CACHE PATH
    "bfam install prefix" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

##
## Note debug mode
##
if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(BFAM_DEBUG TRUE)
else("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(BFAM_DEBUG FALSE)
endif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")

##
## OpenCL
##
if(OPENCL_FOUND)
  set(BFAM_HAVE_OPENCL ON)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif(OPENCL_FOUND)

##
## CUDA
##
if(CUDA_FOUND)
  set(BFAM_HAVE_CUDA ON)
  include_directories(${CUDA_INCLUDE_DIRS})
endif(CUDA_FOUND)


add_executable(bfam_header_gen src/bfam_header_gen.c)

function(FILE2HEADER FILE HEADER)
  add_custom_command (
    OUTPUT ${HEADER}
    COMMAND bfam_header_gen ${FILE} ${HEADER}
    DEPENDS bfam_header_gen ${FILE}
    )
endfunction()

##
## Set compiler flags
##

## Specify c99 compilers
if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endif("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
if("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")
endif("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")

## Enhance debug flags
if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common")
  set(CMAKE_C_FLAGS_DEBUG "-g3")
  set(CMAKE_C_FLAGS_RELEASE "-g0 -O2 -DNDEBUG")
endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common")
  set(CMAKE_C_FLAGS_DEBUG "-g3")
  set(CMAKE_C_FLAGS_RELEASE "-g0 -O2 -DNDEBUG")
endif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common")
  set(CMAKE_C_FLAGS_DEBUG "-g3")
  set(CMAKE_C_FLAGS_RELEASE "-g0 -O2 -DNDEBUG -ftz")
endif("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")

## Enhance stack trace message
if(LINUX)
  if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -rdynamic")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -rdynamic")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -rdynamic")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
endif(LINUX)


# Turn on vectorization flags
if("${CMAKE_BUILD_TYPE}" MATCHES "Release")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native -ftree-vectorize") #  -save-temps -funroll-all-loops -fsplit-ivs-in-unroller
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      # To be able to run with native instructions we need to use the Clang linker.
      # See http://stackoverflow.com/questions/9840207/how-to-use-avx-pclmulqdq-on-mac-os-x-lion
      # for more detials.
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wa,-q")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHost -restrict")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  if ("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Mvect")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
endif("${CMAKE_BUILD_TYPE}" MATCHES "Release")

option(WITH_WARNINGS "Compile with warnings" ON)
if(WITH_WARNINGS)
  # enable compiler warning flags
  if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-parameter -Wno-long-long -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -pedantic -Wextra -Wno-unknown-pragmas")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wextra -Wno-unused-parameter -Wno-long-long -Wno-unknown-pragmas")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -w2 -Wcheck")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  if ("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Minform=inform")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
endif(WITH_WARNINGS)

option(WITH_VECTOR_VERBOSE "Compile with info about vectorization" ON)
if(WITH_VECTOR_VERBOSE)
  # enable verbose vector flags
  if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorizer-verbose=1")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -vec-report3")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  if ("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Minfo=all")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "PGI")
endif(WITH_VECTOR_VERBOSE)

option(WITH_LONG_DOUBLE "Use long double for bfam_long_real_t" ON)
if(WITH_LONG_DOUBLE)
  set(BFAM_USE_LONG_DOUBLE TRUE)
endif(WITH_LONG_DOUBLE)

##
## Look for header files
##
INCLUDE (CheckIncludeFiles)
check_include_files(sysexits.h BFAM_HAVE_SYSEXITS_H)
check_include_files(sys/sysexits.h BFAM_HAVE_SYS_SYSEXITS_H)
check_include_files(time.h BFAM_HAVE_TIME_H)
check_include_files(sys/time.h BFAM_HAVE_SYS_TIME_H)
check_include_files(mach/mach_time.h BFAM_HAVE_MACH_MACH_TIME_H)

cmake_push_check_state()
if(${APPLE})
else(${APPLE})
  set(CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS} -D_XOPEN_SOURCE=600)
endif(${APPLE})
INCLUDE (CheckSymbolExists)
check_symbol_exists(CLOCK_MONOTONIC "time.h;sys/time.h" BFAM_HAVE_CLOCK_MONOTONIC)
INCLUDE (CheckFunctionExists)
check_function_exists(gettimeofday BFAM_HAVE_GETTIMEOFDAY)
cmake_pop_check_state()
if(UNIX AND NOT APPLE)
  set(TIME_LIBRARIES rt)
endif()

check_symbol_exists(mach_absolute_time "mach/mach_time.h" BFAM_HAVE_MACH_ABSOLUTE_TIME)

##
## Map Fortran SYMBOLS we need
##
FortranCInterface_HEADER(bfam_f2c.h MACRO_NAMESPACE "BFAM_" SYMBOL_NAMESPACE "BFAM_"
  SYMBOLS DSTEQR)

##
## Third-Party Libs
##
add_custom_target(build_bundled_libs)
add_custom_target(p4est_bundled_libs)

include(lua)
include_directories(${LUA_INCLUDE_DIR})

include(zlib)
include_directories(${ZLIB_INCLUDE_DIRS})

include(p4est)

option(WITH_OCCA "Compile with occa" ON)
if(WITH_OCCA)
  set(BFAM_HAVE_OCCA ON)
  include(occa2)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OCCA2_CPP_FLAGS} ${OCCA2_C_FLAGS}")
endif()

subdirs(third_party)

##
## Set build include path
##
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR})
include_directories(${MPI_C_INCLUDE_PATH})
include_directories(${MPI_Fortran_INCLUDE_PATH})

##
## Add MPI compile flags for the whole project
##
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAG}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAG}")

##
## BFAM Library
##
add_library(bfam
                 src/bfam_base.c
                 src/bfam_clock.c
                 src/bfam_communicator.c
                 src/bfam_critbit.c
                 src/bfam_dictionary.c
                 src/bfam_gopt.c
                 src/bfam_jacobi.c
                 src/bfam_log.c
                 src/bfam_opencl.c
                 src/bfam_version.c
                 src/bfam_subdomain.c
                 src/bfam_domain.c
                 src/bfam_vtk.c
                 src/bfam_timestep.c
                 src/bfam_timestep_lsrk.c
                 src/bfam_timestep_adams.c
                 src/bfam_timestep_local_adams.c
                 src/bfam_util.c
           )
add_dependencies(bfam build_bundled_libs version)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")

find_library(MATH_LIB m)
target_link_libraries(bfam ${MPI_C_LIBRARIES}
                           ${MPI_Fortran_LIBRARIES}
                           ${OCCA2_LIBRARIES}
                           ${P4EST_LIBRARIES}
                           ${LUA_LIBRARIES}
                           ${ZLIB_LIBRARIES}
                           ${TIME_LIBRARIES}
                           ${MATH_LIB})
if(OPENCL_FOUND)
  target_link_libraries(bfam ${OPENCL_LIBRARIES})
endif(OPENCL_FOUND)

if(CUDA_FOUND)
  target_link_libraries(bfam ${CUDA_LIBRARIES})
endif(CUDA_FOUND)

set(bfam_dgx_norders  1  2  3  4  5  6  7  8  9 10
                     11 12 13 14 15 16 17 18 19 20
                     21 22 23 24)

set(pxest_dimensions 2 3)

set(bfam_dimensions  0 1 2 3)

function(bfam_multiple_build define_name define_vals lib_name lib_src target)
  foreach(val ${define_vals})
    set(full_lib_name "${lib_name}_${val}")
    add_library(${full_lib_name} ${lib_src})

    set_target_properties(${full_lib_name} PROPERTIES COMPILE_FLAGS
      "-D${define_name}=${val}")

    target_link_libraries(${full_lib_name} bfam)
    target_link_libraries(${target} ${full_lib_name})
  endforeach()

  add_library(${lib_name} ${lib_src})
  target_link_libraries(${lib_name} bfam)
  target_link_libraries(${target} ${lib_name})
endfunction()

function(bfam_multiple_add_flag define_name define_val lib_vals lib_name)
  foreach(val ${lib_vals})
    set(full_lib_name "${lib_name}_${val}")

    set_property(TARGET ${full_lib_name} APPEND_STRING PROPERTY
      COMPILE_FLAGS " -D${define_name}=${define_val}")
  endforeach()

  set_property(TARGET ${lib_name} APPEND_STRING PROPERTY
    COMPILE_FLAGS " -D${define_name}=${define_val}")

endfunction()


##
## Examples
##
add_executable(bfam_gopt_usage examples/gopt/bfam_gopt_usage.c)
target_link_libraries(bfam_gopt_usage bfam)

add_executable(bfam_exam_logging examples/logging/bfam_exam_logging.c)
target_link_libraries(bfam_exam_logging bfam)

add_executable(bfam_exam_p4est examples/p4est/bfam_exam_p4est.c)
target_link_libraries(bfam_exam_p4est bfam)

add_executable(bfam_exam_occa examples/occa/bfam_exam_occa.c)
target_link_libraries(bfam_exam_occa bfam)
file(COPY examples/occa/addVectors.occa DESTINATION ${CMAKE_BINARY_DIR})

add_executable(bfam_exam_mesh_dgquad examples/mesh/bfam_exam_mesh_dgquad.c)
target_link_libraries(bfam_exam_mesh_dgquad bfam)

add_executable(bfam_exam_signal examples/signal/bfam_exam_signal.c)
target_link_libraries(bfam_exam_signal bfam)

add_executable(bfam_exam_vectorize examples/vectorize/bfam_exam_vectorize.c)
target_link_libraries(bfam_exam_vectorize bfam)

add_executable(bfam_exam_zpipe examples/zpipe/bfam_exam_zpipe.c)
target_link_libraries(bfam_exam_zpipe bfam)

add_executable(bfam_exam_clock examples/clock/bfam_exam_clock.c)
target_link_libraries(bfam_exam_clock bfam)

add_executable(beard_2d examples/beard/beard.c)
set_target_properties(beard_2d PROPERTIES COMPILE_FLAGS
  "-DBEARD_DGX_DIMENSION=2")

add_executable(beard_3d examples/beard/beard.c)
set_target_properties(beard_3d PROPERTIES COMPILE_FLAGS
  "-DBEARD_DGX_DIMENSION=3")

add_executable(blade_2d examples/blade/blade.c)
set_target_properties(blade_2d PROPERTIES COMPILE_FLAGS
  "-DBLADE_DGX_DIMENSION=2")

add_executable(blade_3d examples/blade/blade.c)
set_target_properties(blade_3d PROPERTIES COMPILE_FLAGS
  "-DBLADE_DGX_DIMENSION=3")

bfam_multiple_build("BFAM_PXEST_DIMENSION" "${pxest_dimensions}"
  bfam_domain_pxest
  src/bfam_domain_pxest.c
  bfam)

bfam_multiple_build("BFAM_DGX_DIMENSION" "${bfam_dimensions}"
  bfam_subdomain_dgx
  src/bfam_subdomain_dgx.c
  bfam)

bfam_multiple_build("NORDER" "${bfam_dgx_norders}"
  beard_dgx_rhs_2
  examples/beard/beard_dgx_rhs.c
  beard_2d)
bfam_multiple_add_flag("BEARD_DGX_DIMENSION" "2" "${bfam_dgx_norders}"
beard_dgx_rhs_2)


bfam_multiple_build("NORDER" "${bfam_dgx_norders}"
  beard_dgx_rhs_3
  examples/beard/beard_dgx_rhs.c
  beard_3d)
bfam_multiple_add_flag("BEARD_DGX_DIMENSION" "3" "${bfam_dgx_norders}"
beard_dgx_rhs_3)

bfam_multiple_build("NORDER" "${bfam_dgx_norders}"
  blade_dgx_rhs_2
  examples/blade/blade_dgx_rhs.c
  blade_2d)
bfam_multiple_add_flag("BLADE_DGX_DIMENSION" "2" "${bfam_dgx_norders}"
blade_dgx_rhs_2)


bfam_multiple_build("NORDER" "${bfam_dgx_norders}"
  blade_dgx_rhs_3
  examples/blade/blade_dgx_rhs.c
  blade_3d)
bfam_multiple_add_flag("BLADE_DGX_DIMENSION" "3" "${bfam_dgx_norders}"
blade_dgx_rhs_3)

##
## Tests
##
enable_testing()

add_executable(bfam_test_fail tests/bfam_test_fail.c)
target_link_libraries(bfam_test_fail bfam)

add_test(NAME bfam_test_fail_test
    COMMAND ${CMAKE_COMMAND}
    "-DTEST_PROG=${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:bfam_test_fail> ${MPIEXEC_POSTFLAGS}"
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_SOURCE_DIR}/cmake/runfailtest.cmake)

add_executable(bfam_test_domain tests/bfam_test_domain.c)
target_link_libraries(bfam_test_domain bfam)

add_test(bfam_test_domain_test bfam_test_domain)


add_executable(bfam_test_domain_fail_repeat tests/bfam_test_domain_fail_repeat.c)
target_link_libraries(bfam_test_domain_fail_repeat bfam)

add_test(NAME bfam_test_domain_fail_repeat
    COMMAND ${CMAKE_COMMAND}
    "-DTEST_PROG=${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:bfam_test_domain_fail_repeat> ${MPIEXEC_POSTFLAGS}"
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_SOURCE_DIR}/cmake/runfailtest.cmake)


add_executable(bfam_test_dictionary tests/bfam_test_dictionary.c)
target_link_libraries(bfam_test_dictionary bfam)

add_test(bfam_test_dictionary_test bfam_test_dictionary)

add_executable(bfam_test_critbit tests/bfam_test_critbit.c)
target_link_libraries(bfam_test_critbit bfam)

add_test(bfam_test_critbit_test bfam_test_critbit)

add_executable(bfam_test_bilinear tests/bfam_test_bilinear.c)
target_link_libraries(bfam_test_bilinear bfam)

add_test(bfam_test_bilinear_test bfam_test_bilinear)

add_executable(bfam_test_interpolation tests/bfam_test_interpolation.c)
target_link_libraries(bfam_test_interpolation bfam)

add_test(bfam_test_interpolation_test bfam_test_interpolation)

add_executable(bfam_test_malloc tests/bfam_test_malloc.c)
target_link_libraries(bfam_test_malloc bfam)

add_test(bfam_test_malloc_test bfam_test_malloc)

add_executable(bfam_test_aligned tests/bfam_test_aligned.c)
target_link_libraries(bfam_test_aligned bfam)

add_test(bfam_test_aligned_test bfam_test_aligned)

add_executable(bfam_test_communicator tests/bfam_test_communicator.c)
target_link_libraries(bfam_test_communicator bfam)

add_executable(bfam_test_mpi tests/bfam_test_mpi.c)
target_link_libraries(bfam_test_mpi bfam)

add_executable(bfam_test_jacobi tests/bfam_test_jacobi.c)
target_link_libraries(bfam_test_jacobi bfam)

add_test(bfam_test_jacobi_test bfam_test_jacobi)

add_executable(bfam_test_util tests/bfam_test_util.c)
target_link_libraries(bfam_test_util bfam)

add_test(bfam_test_util_test bfam_test_util)

add_executable(bfam_test_kron tests/bfam_test_kron.c)
target_link_libraries(bfam_test_kron bfam)

add_test(bfam_test_kron_test bfam_test_kron)

add_executable(bfam_test_mesh_dgx_geo tests/bfam_test_mesh_dgx_geo.c)
target_link_libraries(bfam_test_mesh_dgx_geo bfam)

add_test(bfam_test_mesh_dgx_geo bfam_test_mesh_dgx_geo)

add_executable(bfam_test_dgx tests/bfam_test_dgx.c)
target_link_libraries(bfam_test_dgx bfam)

add_test(bfam_test_dgx bfam_test_dgx)

add_executable(bfam_test_dgx_metric tests/bfam_test_dgx_metric.c)
target_link_libraries(bfam_test_dgx_metric bfam)

add_test(bfam_test_dgx_metric bfam_test_dgx_metric)

add_executable(bfam_test_mesh_dgx tests/bfam_test_mesh_dgx.c)
target_link_libraries(bfam_test_mesh_dgx bfam)

add_test(bfam_test_mesh_dgx_1
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_2
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_3
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx ${MPIEXEC_POSTFLAGS})

add_executable(bfam_test_mesh_dgx_3 tests/bfam_test_mesh_dgx_3.c)
target_link_libraries(bfam_test_mesh_dgx_3 bfam)

add_test(bfam_test_mesh_dgx_3_1
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_3 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_3_2
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_3 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_3_3
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_3 ${MPIEXEC_POSTFLAGS})

add_executable(bfam_test_mesh_dgx_adapt_2 tests/bfam_test_mesh_dgx_adapt_2.c)
target_link_libraries(bfam_test_mesh_dgx_adapt_2 bfam)

add_test(bfam_test_mesh_dgx_2_1
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_adapt_2 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_2_2
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_adapt_2 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_2_3
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_adapt_2 ${MPIEXEC_POSTFLAGS})

add_executable(bfam_test_mesh_dgx_adapt_3 tests/bfam_test_mesh_dgx_adapt_3.c)
target_link_libraries(bfam_test_mesh_dgx_adapt_3 bfam)

add_test(bfam_test_mesh_dgx_3_1
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_adapt_3 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_3_2
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_adapt_3 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_3_3
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_adapt_3 ${MPIEXEC_POSTFLAGS})


add_executable(bfam_test_mesh_dgx_glue_split_3
               tests/bfam_test_mesh_dgx_glue_split_3.c)
target_link_libraries(bfam_test_mesh_dgx_glue_split_3 bfam)

add_test(bfam_test_mesh_dgx_glue_split_3_1
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_glue_split_3 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_glue_split_3_2
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_glue_split_3 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_mesh_dgx_glue_split_3_3
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
  ${MPIEXEC_PREFLAGS} ./bfam_test_mesh_dgx_glue_split_3 ${MPIEXEC_POSTFLAGS})

add_test(bfam_test_communicator_test
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
  ${MPIEXEC_PREFLAGS} ./bfam_test_communicator ${MPIEXEC_POSTFLAGS})


add_test(bfam_test_mpi_test
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
  ${MPIEXEC_PREFLAGS} ./bfam_test_mpi ${MPIEXEC_POSTFLAGS})

file2header(${CMAKE_SOURCE_DIR}/tests/bfam_test_opencl_kernels.cl
            ${CMAKE_BINARY_DIR}/bfam_test_opencl_kernels.h)

add_executable(bfam_test_opencl tests/bfam_test_opencl.c
  ${CMAKE_BINARY_DIR}/bfam_test_opencl_kernels.h)
target_link_libraries(bfam_test_opencl bfam)
add_test(bfam_test_opencl_test bfam_test_opencl)

file2header(${CMAKE_SOURCE_DIR}/tests/bfam_test_lua_script.lua
            ${CMAKE_BINARY_DIR}/bfam_test_lua_script.h)

add_executable(bfam_test_lua tests/bfam_test_lua.c
  ${CMAKE_BINARY_DIR}/bfam_test_lua_script.h)
target_link_libraries(bfam_test_lua bfam)
add_test(bfam_test_lua_test bfam_test_lua)

##
## Print configuration
##
set(BFAM_COMPILE_LIST
  "             BFAM_DEBUG:                  ${BFAM_DEBUG}\\\\n"
  "             BFAM_HAVE_CUDA:              ${BFAM_HAVE_CUDA}\\\\n"
  "             BFAM_HAVE_OPENCL:            ${BFAM_HAVE_OPENCL}\\\\n"
  "             BFAM_HAVE_OCCA:              ${BFAM_HAVE_OCCA}\\\\n"
  "             BFAM_HAVE_SYSEXITS_H:        ${BFAM_HAVE_SYSEXITS_H}\\\\n"
  "             BFAM_HAVE_SYS_SYSEXITS_H:    ${BFAM_HAVE_SYS_SYSEXITS_H}\\\\n"
  "             CMAKE_C_COMPILER:            ${CMAKE_C_COMPILER}\\\\n"
  "             CMAKE_C_FLAGS:               ${CMAKE_C_FLAGS}\\\\n"
  "             CMAKE_C_FLAGS_DEBUG:         ${CMAKE_C_FLAGS_DEBUG}\\\\n"
  "             CMAKE_C_FLAGS_RELEASE:       ${CMAKE_C_FLAGS_RELEASE}\\\\n"
  "             CMAKE_CXX_COMPILER:          ${CMAKE_CXX_COMPILER}\\\\n"
  "             CMAKE_CXX_FLAGS:             ${CMAKE_CXX_FLAGS}\\\\n"
  "             CMAKE_CXX_FLAGS_DEBUG:       ${CMAKE_CXX_FLAGS_DEBUG}\\\\n"
  "             CMAKE_CXX_FLAGS_RELEASE:     ${CMAKE_CXX_FLAGS_RELEASE}\\\\n"
  "             CMAKE_Fortran_COMPILER:      ${CMAKE_Fortran_COMPILER}\\\\n"
  "             CMAKE_Fortran_FLAGS:         ${CMAKE_Fortran_FLAGS}\\\\n"
  "             CMAKE_Fortran_FLAGS_DEBUG:   ${CMAKE_Fortran_FLAGS_DEBUG}\\\\n"
  "             CMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}\\\\n"
  "             OCCA2_INCLUDE_DIRS:          ${OCCA2_INCLUDE_DIRS}\\\\n"
  "             OCCA2_LIBRARIES:             ${OCCA2_LIBRARIES}\\\\n"
  "             OCCA2_C_FLAGS:               ${OCCA2_C_FLAGS}\\\\n"
  "             OCCA2_CXX_FLAGS:             ${OCCA2_CXX_FLAGS}\\\\n"
  "             OCCA2_CPP_FLAGS:             ${OCCA2_CPP_FLAGS}\\\\n"
  "             LUA_LIBRARIES:               ${LUA_LIBRARIES}\\\\n"
  "             LUA_INCLUDE_DIR:             ${LUA_INCLUDE_DIR}\\\\n"
  "             OPENCL_FOUND:                ${OPENCL_FOUND}\\\\n"
  "             OPENCL_INCLUDE_DIRS:         ${OPENCL_INCLUDE_DIRS}\\\\n"
  "             OPENCL_LIBRARIES:            ${OPENCL_LIBRARIES}\\\\n"
  "             CUDA_FOUND:                  ${CUDA_FOUND}\\\\n"
  "             CUDA_VERSION:                ${CUDA_VERSION}\\\\n"
  "             CUDA_INCLUDE_DIRS:           ${CUDA_INCLUDE_DIRS}\\\\n"
  "             CUDA_CUDA_LIBRARY:           ${CUDA_CUDA_LIBRARY}\\\\n"
  "             CUDA_NVCC_EXECUTABLE:        ${CUDA_NVCC_EXECUTABLE}\\\\n"
  "             OPENMP_FOUND:                ${OPENMP_FOUND}\\\\n"
  "             OpenMP_C_FLAGS:              ${OpenMP_C_FLAGS}\\\\n"
  "             OpenMP_CXX_FLAGS:            ${OpenMP_CXX_FLAGS}\\\\n"
  "             MPI_C_FOUND:                 ${MPI_C_FOUND}\\\\n"
  "             MPI_C_COMPILER:              ${MPI_C_COMPILER}\\\\n"
  "             MPI_C_COMPILE_FLAGS:         ${MPI_C_COMPILE_FLAGS}\\\\n"
  "             MPI_C_INCLUDE_PATH:          ${MPI_C_INCLUDE_PATH}\\\\n"
  "             MPI_C_LINK_FLAGS:            ${MPI_C_LINK_FLAGS}\\\\n"
  "             MPI_C_LIBRARIES:             ${MPI_C_LIBRARIES}\\\\n"
  "             MPI_Fortran_FOUND:           ${MPI_Fortran_FOUND}\\\\n"
  "             MPI_Fortran_COMPILER:        ${MPI_Fortran_COMPILER}\\\\n"
  "             MPI_Fortran_COMPILE_FLAGS:   ${MPI_Fortran_COMPILE_FLAGS}\\\\n"
  "             MPI_Fortran_INCLUDE_PATH:    ${MPI_Fortran_INCLUDE_PATH}\\\\n"
  "             MPI_Fortran_LINK_FLAGS:      ${MPI_Fortran_LINK_FLAGS}\\\\n"
  "             MPI_Fortran_LIBRARIES:       ${MPI_Fortran_LIBRARIES}\\\\n"
  "             MPIEXEC:                     ${MPIEXEC}\\\\n"
  "             MPIEXEC_NUMPROC_FLAG:        ${MPIEXEC_NUMPROC_FLAG}\\\\n"
  "             MPIEXEC_PREFLAGS:            ${MPIEXEC_PREFLAGS}\\\\n"
  "             MPIEXEC_POSTFLAGS:           ${MPIEXEC_POSTFLAGS}")

function(JOIN VALUES GLUE OUTPUT)
  string (REGEX REPLACE "([^\\]|^);" "\\1${GLUE}" _TMP_STR "${VALUES}")
  string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction()

join("${BFAM_COMPILE_LIST}" "\" \"" BFAM_COMPILE_INFO)

set(slash "\\")
join("${bfam_dgx_norders}" ") X(" BFAM_DGX_NORDERS_MID)

set(slash "\\")
join("${bfam_dimensions}" ") X(" BFAM_DGX_DIMENSIONS_MID)

set(slash "\\")
join("${pxest_dimensions}" ") X(" BFAM_PXEST_DIMENSIONS_MID)

##
## Write BFAM config file
##
configure_file(${CMAKE_SOURCE_DIR}/src/bfam_config.h.in
               ${CMAKE_BINARY_DIR}/bfam_config.h)

# Set the CPack options
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "bfam")
set(CPACK_PACKAGE_VENDOR "")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${bfam_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_VERSION_MAJOR ${bfam_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${bfam_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${bfam_VERSION_PATCH})
set(CPACK_PACKAGE_INSTALL_DIRECTORY "bfam-${bfam_VERSION_FULL}")
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "${CMAKE_SOURCE_DIR};/;${CMAKE_BINARY_DIR}/source;/")
set(CPACK_SOURCE_IGNORE_FILES
"~$"
".sw[a-z]$"
".DS_Store"
"tags"
"cscope.*"
"^${PROJECT_SOURCE_DIR}/.git.*"
"^${PROJECT_SOURCE_DIR}/build.*"
)

# This should be the last line
include(CPack)

cmake_minimum_required(VERSION 2.8.8)

project(bfam C Fortran)

include(InstallRequiredSystemLibraries)
include(FortranCInterface)

find_package(Git)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)

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
  file(WRITE ${CMAKE_SOURCE_DIR}/VERSION \"\${GIT_DESCRIBE_VERSION}\")
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
## Set compiler flags
##
if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fno-common")
  set(CMAKE_C_FLAGS_DEBUG "-g3")
  set(CMAKE_C_FLAGS_RELEASE "-g0 -O2 -DNDEBUG")
endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")

if("${CMAKE_BUILD_TYPE}" MATCHES "Release")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native -ftree-vectorize -funroll-all-loops -fsplit-ivs-in-unroller")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
endif("${CMAKE_BUILD_TYPE}" MATCHES "Release")

option(WITH_WARNINGS "Compile with warnings" ON)
if(WITH_WARNINGS)
  # enable compiler warning flags
  if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-parameter -Wno-long-long -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -pedantic -Wextra")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wextra -Wno-unused-parameter -Wno-long-long")
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
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorizer-verbose=2")
  endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
endif(WITH_VECTOR_VERBOSE)

##
## Look for header files
##
check_include_files(sysexits.h BFAM_HAVE_SYSEXITS_H)
check_include_files(sys/sysexits.h BFAM_HAVE_SYS_SYSEXITS_H)

##
## Map Fortran SYMBOLS we need
##
FortranCInterface_HEADER(bfam_f2c.h MACRO_NAMESPACE "BFAM_" SYMBOL_NAMESPACE "BFAM_"
  SYMBOLS DSTEQR)

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
## Write BFAM config file
##
configure_file(${CMAKE_SOURCE_DIR}/src/bfam_config.h.in
               ${CMAKE_BINARY_DIR}/bfam_config.h)

##
## BFAM Library
##
add_library(bfam
                 src/bfam_base.c
                 src/bfam_log.c
                 src/bfam_version.c
           )
add_dependencies(bfam version)

set_property(TARGET bfam APPEND_STRING PROPERTY LINK_FLAGS ${MPI_C_LINK_FLAGS})
set_property(TARGET bfam APPEND_STRING PROPERTY LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
set_property(TARGET bfam APPEND_STRING PROPERTY LINK_FLAGS ${LAPACK_LINKER_FLAGS})
set_property(TARGET bfam APPEND_STRING PROPERTY LINK_FLAGS ${BLAS_LINKER_FLAGS})

target_link_libraries(bfam ${MPI_C_LIBRARIES}
                           ${MPI_Fortran_LIBRARIES}
                           ${LAPACK_LIBRARIES}
                           ${BLAS_LIBRARIES})

add_executable(bfam_exam_logging examples/logging/bfam_exam_logging.c)
target_link_libraries(bfam_exam_logging bfam)

##
## Tests
##
enable_testing()

add_executable(bfam_test_mpi tests/bfam_test_mpi.c)
target_link_libraries(bfam_test_mpi bfam)

add_test(bfam_test_mpi_test
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
  ${MPIEXEC_PREFLAGS} bfam_test_mpi ${MPIEXEC_POSTFLAGS})


##
## Print configuration
##
message("============== Configure information ==============")
message("BFAM_DEBUG:                ${BFAM_DEBUG}")
message("CMAKE_C_FLAGS:             ${CMAKE_C_FLAGS}")
message("CMAKE_C_FLAGS_DEBUG:       ${CMAKE_C_FLAGS_DEBUG}")
message("CMAKE_C_FLAGS_RELEASE:     ${CMAKE_C_FLAGS_RELEASE}")
message("BLAS_FOUND:                ${BLAS_FOUND}")
message("BLAS_LINKER_FLAGS:         ${BLAS_LINKER_FLAGS}")
message("BLAS_LIBRARIES:            ${BLAS_LIBRARIES}")
message("LAPACK_FOUND:              ${LAPACK_FOUND}")
message("LAPACK_LINKER_FLAGS:       ${LAPACK_LINKER_FLAGS}")
message("LAPACK_LIBRARIES:          ${LAPACK_LIBRARIES}")
message("MPI_C_FOUND:               ${MPI_C_FOUND}")
message("MPI_C_COMPILER:            ${MPI_C_COMPILER}")
message("MPI_C_COMPILE_FLAGS:       ${MPI_C_COMPILE_FLAGS}")
message("MPI_C_INCLUDE_PATH:        ${MPI_C_INCLUDE_PATH}")
message("MPI_C_LINK_FLAGS:          ${MPI_C_LINK_FLAGS}")
message("MPI_C_LIBRARIES:           ${MPI_C_LIBRARIES}")
message("MPI_Fortran_FOUND:         ${MPI_Fortran_FOUND}")
message("MPI_Fortran_COMPILER:      ${MPI_Fortran_COMPILER}")
message("MPI_Fortran_COMPILE_FLAGS: ${MPI_Fortran_COMPILE_FLAGS}")
message("MPI_Fortran_INCLUDE_PATH:  ${MPI_Fortran_INCLUDE_PATH}")
message("MPI_Fortran_LINK_FLAGS:    ${MPI_Fortran_LINK_FLAGS}")
message("MPI_Fortran_LIBRARIES:     ${MPI_Fortran_LIBRARIES}")
message("MPIEXEC:                   ${MPIEXEC}")
message("MPIEXEC_NUMPROC_FLAG:      ${MPIEXEC_NUMPROC_FLAG}")
message("MPIEXEC_PREFLAGS:          ${MPIEXEC_PREFLAGS}")
message("MPIEXEC_POSTFLAGS:         ${MPIEXEC_POSTFLAGS}")

# Set the CPack options
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "bfam")
set(CPACK_PACKAGE_VENDOR "")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${bfam_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_VERSION_MAJOR ${bfam_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${bfam_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${bfam_VERSION_PATCH})
set(CPACK_PACKAGE_INSTALL_DIRECTORY "bfam-${bfam_VERSION_FULL}")
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

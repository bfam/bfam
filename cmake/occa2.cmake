#
# Bundled occa2 paths.
#
set(OCCA2_BUNDLED_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(OCCA2_BUNDLED_LIBRARIES "${OCCA2_BUNDLED_PREFIX}/lib/libocca.so")

set(OCCA2_BUILD_PREFIX "${PROJECT_BINARY_DIR}/third_party/OCCA2")

set(OCCA2_C_FLAGS "")
set(OCCA2_CXX_FLAGS "")
set(OCCA2_CPP_FLAGS "-DLINUX_OS=1 -DOSX_OS=2 -DWINDOWS_OS=4 -DWINUX_OS=5")

macro(occa2_use_bundled)
  set(OCCA2_PREFIX "${OCCA2_BUNDLED_PREFIX}")
  set(OCCA2_INCLUDE_DIRS "${OCCA2_BUNDLED_PREFIX}/include")
  set(OCCA2_LIBRARIES "${OCCA2_BUNDLED_LIBRARIES}")

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_OS=OSX_OS")
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_OS=LINUX_OS")
  else()
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_OS=UNKNOWN_OS")
  endif()

  if (OPENMP_FOUND)
    set(OCCA2_C_FLAGS   "${OCCA2_C_FLAGS}   ${OpenMP_C_FLAGS}")
    set(OCCA2_CXX_FLAGS "${OCCA2_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_OPENMP_ENABLED=1")
  else()
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_OPENMP_ENABLED=0")
  endif()

  if (OPENCL_FOUND)
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_OPENCL_ENABLED=1")
  else()
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_OPENCL_ENABLED=0")
  endif()

  if (CUDA_FOUND)
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_CUDA_ENABLED=1")
  else()
    set(OCCA2_CPP_FLAGS "${OCCA2_CPP_FLAGS} -DOCCA_CUDA_ENABLED=0")
  endif()

  set(ENABLE_BUNDLED_OCCA2 True)
endmacro()

macro(occa2_try_system)
  find_path(OCCA2_OCCA2_INCLUDE_DIR occa2.h PATH_SUFFIXES include)
  find_library(OCCA2_OCCA2_LIB NAMES occa2 PATH_SUFFIXES lib)

  if(OCCA2_OCCA2_INCLUDE_DIRS AND OCCA2_OCCA2_LIB)
    message (STATUS "-- occa2 include: ${OCCA2_OCCA2_INCLUDE_DIR}, lib: ${OCCA2_OCCA2_LIB}")
    message (STATUS "Found a system-wide occa2.")
    set(OCCA2_INCLUDE_DIRS ${OCCA2_OCCA2_INCLUDE_DIR})
    set(OCCA2_LIBRARIES ${OCCA2_OCCA2_LIB})
    set(OCCA2_CPP_FLAGS "")
  else()
    message (FATAL_ERROR "Not found a system occa2")
    #occa2_use_bundled()
  endif()
endmacro()

#
# Check if there is a usable occa2 at the given prefix path.
#
macro (occa2_try_prefix)
  find_path(OCCA2_OCCA2_INCLUDE_DIR occa2.h ${OCCA2_PREFIX} NO_DEFAULT_PATH)
  find_library(OCCA2_OCCA2_LIB occa2 ${OCCA2_PREFIX} NO_DEFAULT_PATH)


  if(OCCA2_OCCA2_INCLUDE_DIRS AND OCCA2_OCCA2_LIB)
    set(OCCA2_INCLUDE_DIRS ${OCCA2_OCCA2_INCLUDE_DIR})
    set(OCCA2_LIBRARIES ${OCCA2_OCCA2_LIB})
    include_directories(${OCCA2_INCLUDE_DIRS})
    set(OCCA2_CPP_FLAGS "")
  else()
    message(FATAL_ERROR "Couldn't find occa2 in '${OCCA2_PREFIX}'")
  endif()
endmacro()

#
# occa2 options.
#
option(ENABLE_BUNDLED_OCCA2 "Enable building of the bundled occa2" ON)
option(OCCA2_PREFIX "Build with occa2 at the given path" "")

if(OCCA2_PREFIX AND ENABLE_BUNDLED_OCCA2)
  message(FATAL_ERROR "Options OCCA2_PREFIX and ENABLE_BUNDLED_OCCA2 "
    "are not compatible with each other.")
endif()

if (OCCA2_PREFIX)
  # trying to build with specified occa2.
  occa2_try_prefix()
elseif (NOT ENABLE_BUNDLED_OCCA2)
  # trying to build with system occa2, macro can turn on
  # building of occa2 bundled with the server source.
  occa2_try_system()
else()
  occa2_use_bundled()
endif()

include_directories(${OCCA2_INCLUDE_DIRS})

message(STATUS "Use occa2 includes: ${OCCA2_INCLUDE_DIRS}")
message(STATUS "Use occa2 library: ${OCCA2_LIBRARIES}")
message(STATUS "Use occa2 c flags: ${OCCA2_C_FLAGS}")
message(STATUS "Use occa2 cxx flags: ${OCCA2_CXX_FLAGS}")
message(STATUS "Use occa2 cpp flags: ${OCCA2_CPP_FLAGS}")

macro(occa2_build)
  set (occa2_buildoptions "OCCA_DIR='${OCCA2_BUILD_PREFIX}'")
  set (occa2_cxx "${CMAKE_CXX_COMPILER}")
  set (occa2_cxxflags "")

  if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set (occa2_cxxflags ${occa2_cxxflags} -O0 -g)
  else ()
    set (occa2_cxxflags ${occa2_cxxflags} -O3)
  endif()

  if (OPENMP_FOUND)
    set (occa2_buildoptions ${occa2_buildoptions} "OCCA_OPENMP_ENABLED='1'")
    set (occa2_cxxflags ${occa2_cxxflags} ${OpenMP_CXX_FLAGS})
  endif()

  set (occa2_buildoptions ${occa2_buildoptions} "CXX='${occa2_cxx}'")
  set (occa2_buildoptions ${occa2_buildoptions} "CXXFLAGS='${occa2_cxxflags}'")

  if ("${PROJECT_BINARY_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
    add_custom_command(OUTPUT ${OCCA2_BUNDLED_LIBRARIES}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/third_party/OCCA2
      COMMAND $(MAKE) ${occa2_buildoptions} clean
      COMMAND $(MAKE) ${occa2_buildoptions}
      COMMAND sh ${PROJECT_SOURCE_DIR}/cmake/install_occa.sh ${CMAKE_INSTALL_PREFIX}
      DEPENDS ${CMAKE_SOURCE_DIR}/CMakeCache.txt
              ${PROJECT_SOURCE_DIR}/cmake/install_occa.sh
      )
  else()
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/third_party/OCCA2
      COMMAND mkdir ${PROJECT_BINARY_DIR}/third_party/OCCA2
      )
    add_custom_command(OUTPUT ${OCCA2_BUNDLED_LIBRARIES}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/third_party/OCCA2
      COMMAND cp -r ${PROJECT_SOURCE_DIR}/third_party/OCCA2/* .
      COMMAND $(MAKE) ${occa2_buildoptions} clean
      COMMAND $(MAKE) ${occa2_buildoptions}
      COMMAND sh ${PROJECT_SOURCE_DIR}/cmake/install_occa.sh ${CMAKE_INSTALL_PREFIX}
      DEPENDS ${PROJECT_BINARY_DIR}/CMakeCache.txt ${PROJECT_BINARY_DIR}/third_party/OCCA2
              ${PROJECT_SOURCE_DIR}/cmake/install_occa.sh
      )
  endif()
  add_custom_target(libocca2 DEPENDS ${OCCA2_BUNDLED_LIBRARIES})
  add_dependencies(build_bundled_libs libocca2)
  unset (occa2_buildoptions)
endmacro()

if (ENABLE_BUNDLED_OCCA2)
  occa2_build()
endif()

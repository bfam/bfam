# Copyright (C) 2010, 2011, 2012 Tarantool/Box AUTHORS:
#
#   Aleksey Demakov, Aleksey Mashanov, Alexandre Kalendarev,
#   Damien Lefortier, Dmitry E. Oboukhov, Dmitry Simonenko,
#   Konstantin Osipov, Konstantin Shulgin, Mons Anderson,
#   Pavel Cherenkov, Roman Antipin, Roman Tokarev, Roman Tsisyk,
#   Teodor Sigaev, Timofey Khryukin, Yuriy Nevinitsin, Yuriy Vostrikov
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above
#    copyright notice, this list of conditions and the
#    following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials
#    provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# <COPYRIGHT HOLDER> OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
#
# LuaJIT configuration file.
#
# A copy of LuaJIT is maintained within BFAM
# source. It's located in third_party/luajit.
#
# Instead of this copy, BFAM can be compiled
# with a system-wide LuaJIT, or LuaJIT at a given
# prefix. This is used when compiling BFAM
# as part of a distribution, e.g. Debian.
#
# To explicitly request use of the bundled LuaJIT,
# add -DENABLE_BUNDLED_LUAJIT=True to CMake
# configuration flags.
# To explicitly request use of LuaJIT at a given
# prefix, use -DLUAJIT_PREFIX=/path/to/LuaJIT.
#
# These two options are incompatible with each other.
#
# If neither of the two options is given, this script
# first attempts to use the system-installed LuaJIT
# and, in case it is not present or can not be used,
# falls back to the bundled one.
#
# Adds CMake options: ENABLED_BUNDLED_LUAJIT, LUAJIT_PREFIX
# Exports CMake defines: LUAJIT_PREFIX, LUAJIT_INCLUDE, LUAJIT_LIB
# Modifies CMAKE_CFLAGS with -I${LUAJIT_INCLUDE}
#

#
# Bundled LuaJIT paths.
#
set (LUAJIT_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/third_party/luajit/src")
set (LUAJIT_BUNDLED_LIB "${LUAJIT_BUNDLED_PREFIX}/libluajit.a")

macro (luajit_use_bundled)
  set (LUAJIT_PREFIX "${LUAJIT_BUNDLED_PREFIX}")
  set (LUAJIT_INCLUDE "${PROJECT_SOURCE_DIR}/third_party/luajit/src")
  set (LUAJIT_LIB "${LUAJIT_BUNDLED_LIB}")
  set (ENABLE_BUNDLED_LUAJIT True)
endmacro()

#
# LuaJIT testing routine
# (see cmake/luatest.cpp for a description).
#
#
# Check if there is a system LuaJIT availaible and
# usable with the server (determined by a compiled test).
#
macro (luajit_try_system)
  find_path (LUAJIT_INCLUDE lj_obj.h PATH_SUFFIXES luajit-2.0 luajit)
  find_library (LUAJIT_LIB NAMES luajit luajit-5.1 PATH_SUFFIXES x86_64-linux-gnu)
  if (LUAJIT_INCLUDE AND LUAJIT_LIB)
    message (STATUS "include: ${LUAJIT_INCLUDE}, lib: ${LUAJIT_LIB}")
    message (STATUS "Found a system-wide LuaJIT.")
  else()
    message (FATAL_ERROR "Not found a system LuaJIT")
    #luajit_use_bundled()
  endif()
endmacro()

#
# Check if there is a usable LuaJIT at the given prefix path.
#
macro (luajit_try_prefix)
  find_path (LUAJIT_INCLUDE "lua.h" ${LUAJIT_PREFIX} NO_DEFAULT_PATH)
  find_library (LUAJIT_LIB "luajit" ${LUAJIT_PREFIX} NO_DEFAULT_PATH)
  if (LUAJIT_INCLUDE AND LUAJIT_LIB)
    include_directories("${LUAJIT_INCLUDE}")
  else()
    message (FATAL_ERROR "Couldn't find LuaJIT in '${LUAJIT_PREFIX}'")
  endif()
endmacro()

#
# LuaJIT options.
#
option(ENABLE_BUNDLED_LUAJIT "Enable building of the bundled LuaJIT" ON)
option(LUAJIT_PREFIX "Build with LuaJIT at the given path" "")

if (LUAJIT_PREFIX AND ENABLE_BUNDLED_LUAJIT)
  message (FATAL_ERROR "Options LUAJIT_PREFIX and ENABLE_BUNDLED_LUAJIT "
    "are not compatible with each other.")
endif()

if (LUAJIT_PREFIX)
  # trying to build with specified LuaJIT.
  luajit_try_prefix()
elseif (NOT ENABLE_BUNDLED_LUAJIT)
  # trying to build with system LuaJIT, macro can turn on
  # building of LuaJIT bundled with the server source.
  luajit_try_system()
else()
  luajit_use_bundled()
endif()

include_directories("${LUAJIT_INCLUDE}")

if(APPLE)
  set(LUAJIT_LINK_FLAGS "-pagezero_size 10000 -image_base 100000000")
endif(APPLE)

message (STATUS "Use LuaJIT includes: ${LUAJIT_INCLUDE}")
message (STATUS "Use LuaJIT library: ${LUAJIT_LIB}")
message (STATUS "Use LuaJIT link flags: ${LUAJIT_LINK_FLAGS}")

macro(luajit_build)
  set (luajit_buildoptions BUILDMODE=static)
  set (luajit_copt "")
  if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set (luajit_buildoptions ${luajit_buildoptions} CCDEBUG=${CC_DEBUG_OPT})
    set (luajit_copt ${luajit_copt} -O1)
    set (luajit_buildoptions ${luajit_buildoptions} XCFLAGS='-DLUA_USE_APICHECK -DLUA_USE_ASSERT')
  else ()
    set (luajit_copt ${luajit_copt} -O2)
  endif()
  set (luajit_copt ${luajit_copt})

  # XXX: Hack to avoid using ccc-analyzer
  if(${CMAKE_C_COMPILER} MATCHES "ccc-analyzer")
    if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
      set (luajit_cc gcc)
    endif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
      set (luajit_cc clang)
    endif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
      set (luajit_cc icc)
    endif("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
  else(${CMAKE_C_COMPILER} MATCHES "ccc-analyzer")
    set (luajit_cc ${CMAKE_C_COMPILER})
  endif(${CMAKE_C_COMPILER} MATCHES "ccc-analyzer")

  if (NOT luajit_cc)
    message (FATAL_ERROR "LuaJIT will not compile with default C compiler (cc)")
  endif()
  set(luajit_ldflags "${CMAKE_SHARED_LINKER_FLAGS}")
  separate_arguments(luajit_copt)
  separate_arguments(luajit_ldflags)
  set (luajit_buildoptions ${luajit_buildoptions} CC="${luajit_cc}" TARGET_CC="${luajit_cc}" CCOPT="${luajit_copt}")
  set (luajit_buildoptions ${luajit_buildoptions} Q='' LDFLAGS="${luajit_ldflags}")
  if ("${PROJECT_BINARY_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/third_party/luajit/src/libluajit.a
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/third_party/luajit
      COMMAND $(MAKE) clean
      COMMAND $(MAKE) -C src ${luajit_buildoptions}
      DEPENDS ${CMAKE_SOURCE_DIR}/CMakeCache.txt
      )
  else()
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/third_party/luajit
      COMMAND mkdir ${PROJECT_BINARY_DIR}/third_party/luajit
      )
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/third_party/luajit/src/libluajit.a
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/third_party/luajit
      COMMAND cp -r ${PROJECT_SOURCE_DIR}/third_party/luajit/* .
      COMMAND $(MAKE) clean
      COMMAND $(MAKE) -C src ${luajit_buildoptions}
      DEPENDS ${PROJECT_BINARY_DIR}/CMakeCache.txt ${PROJECT_BINARY_DIR}/third_party/luajit
      )
  endif()
  add_custom_target(libluajit
    DEPENDS ${PROJECT_BINARY_DIR}/third_party/luajit/src/libluajit.a
    )
  add_dependencies(build_bundled_libs libluajit)
  add_dependencies(p4est_bundled_libs libluajit)
  unset (luajit_buildoptions)
endmacro()



#
# Building shipped luajit only if there is no
# usable system one (see cmake/luajit.cmake) or by demand.
#
if (ENABLE_BUNDLED_LUAJIT)
  luajit_build()
endif()

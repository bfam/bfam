include (ExternalProject)

#
# Bundled lua paths.
#
set(LUA_BUNDLED_PREFIX "${PROJECT_BINARY_DIR}/third_party/lua/install")
set(LUA_BUNDLED_LIB    "${LUA_BUNDLED_PREFIX}/lib/liblua.a")

macro(lua_use_bundled)
  find_library ( DL_LIBRARY NAMES dl )
  if ( DL_LIBRARY )
    list ( APPEND LIBS ${DL_LIBRARY} )
  endif ( )
  set(LUA_PREFIX "${LUA_BUNDLED_PREFIX}")
  set(LUA_INCLUDE_DIR "${LUA_BUNDLED_PREFIX}/include")
  set(LUA_LIBRARIES "${LUA_BUNDLED_LIB}" "${DL_LIBRARY}")
  set(ENABLE_BUNDLED_LUA True)
endmacro()

macro(lua_try_system)
  find_package(Lua)
  if(LUA_FOUND)
    message(STATUS "include: ${LUA_INCLUDE_DIR}, lib: ${LUA_LIBRARIES}")
    message(STATUS "Found a system-wide lua.")
  else()
    message(FATAL_ERROR "Not found a system lua")
    #lua_use_bundled()
  endif()
endmacro()

#
# lua options.
#
option(ENABLE_BUNDLED_LUA "Enable building of the bundled lua" ON)

if (NOT ENABLE_BUNDLED_LUA)
  lua_try_system()
else()
  lua_use_bundled()
endif()

message(STATUS "Use lua include: ${LUA_INCLUDE_DIR}")
message(STATUS "Use lua library: ${LUA_LIBRARIES}")

macro(lua_build)
  ExternalProject_Add(lua
    PREFIX              ${CMAKE_BINARY_DIR}/third_party/lua
    URL                 ${CMAKE_SOURCE_DIR}/third_party/lua-5.2.3.tar.gz
    URL_MD5             fa17304d9d80870d19eb4461e962c2ee
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_BINARY_DIR}/third_party/lua/install
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
  )
  set_target_properties(lua PROPERTIES EXCLUDE_FROM_ALL ON)
  add_dependencies(build_bundled_libs lua)
  add_dependencies(p4est_bundled_libs lua)
endmacro()


if(ENABLE_BUNDLED_LUA)
  lua_build()
endif()

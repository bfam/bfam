//
// Example calling Lua from:
//
//   http://www.cs.usfca.edu/~galles/cs420/lecture/LuaAndC.html
//
#include <bfam.h>
#include <bfam_test_lua_script.h>

int main(int argc, char *argv[])
{
  lua_State *L;
  int status;

  // initialize Lua interpreter
  L = luaL_newstate();

  // load Lua base libraries (print / math / etc)
  luaL_openlibs(L);

  // run script
  status = luaL_dostring(L, bfam_test_lua_script);
  if (status)
  {
    fprintf(stderr, "Couldn't do script: %s\n", lua_tostring(L, -1));
    exit(1);
  }

  // Push the fib function on the top of the lua stack
  lua_getglobal(L, "fib");

  // Push the argument (the number 13) on the stack
  lua_pushnumber(L, 13);

  // call the function with 1 argument, returning a single result.  Note that
  // the function actually returns 2 results -- we just want one of them.  The
  // second result will *not* be pushed on the lua stack, so we don't need to
  // clean up after it
  lua_call(L, 1, 1);

  // Get the result from the lua stack
  int result = (int)lua_tointeger(L, -1);

  // Clean up.  If we don't do this last step, we'll leak stack memory.
  lua_pop(L, 1);

  // Cleanup:  Deallocate all space assocatated with the lua state
  lua_close(L);

  if (result == 233)
    return 0;
  else
    return 1;
}

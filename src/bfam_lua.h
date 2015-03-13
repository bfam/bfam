#ifndef BFAM_LUA_H
#define BFAM_LUA_H

#include <bfam_base.h>

/*
 * Som of the following functions were modified from the code at:
 *
 *    http://windrealm.org/tutorials/reading-a-lua-configuration-file-from-c.php
 */

/** Evaluates a Lua expression and returns the boolean result.
 *
 * If an error occurs or the result is not a boolean, def is returned.
 */
int bfam_lua_expr_boolean(lua_State *L, const char *expr, int def);

/** Evaluates a Lua expression and returns the integer result.
 *
 * If an error occurs or the result is not a integer, def is returned.
 */
lua_Integer bfam_lua_expr_integer(lua_State *L, const char *expr,
                                  lua_Integer def);

/** Evaluates a Lua expression and returns the number result.
 *
 * If an error occurs or the result is not a number, def is returned.
 */
lua_Number bfam_lua_expr_number(lua_State *L, const char *expr, lua_Number def);

/** Evaluates a Lua expression and returns the string result.
 *
 * If an error occurs or the result is not a string, a copy of def is returned.
 * The user of this function is responsible for cleaning up the memory.
 *
 */
char *bfam_lua_expr_string(lua_State *L, const char *expr, const char *def);

#endif

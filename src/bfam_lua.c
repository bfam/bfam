#include <bfam_base.h>
#include <bfam_lua.h>

#define BFAM_LUA_MAX_COMMAND_LEN 4096
#define BFAM_LUA_EVALEXP_VAR "XXX_bfam_evalexpr_XXX"

int bfam_lua_expr_boolean(lua_State *L, const char *expr, int def)
{
  int r = def;
  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
    if (lua_isboolean(L, -1))
      r = lua_toboolean(L, -1);
    /* remove lua_getglobal value */
    lua_pop(L, 1);
  }
  return r;
}

lua_Integer bfam_lua_expr_integer(lua_State *L, const char *expr,
                                  lua_Integer def)
{
  lua_Integer r = def;
  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
    if (lua_isnumber(L, -1))
      r = lua_tointeger(L, -1);
    /* remove lua_getglobal value */
    lua_pop(L, 1);
  }
  return r;
}

lua_Number bfam_lua_expr_number(lua_State *L, const char *expr, lua_Number def)
{
  lua_Number r = def;
  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
    if (lua_isnumber(L, -1))
      r = lua_tonumber(L, -1);
    /* remove lua_getglobal value */
    lua_pop(L, 1);
  }
  return r;
}

char *bfam_lua_expr_string(lua_State *L, const char *expr, const char *def)
{
  const char *r = def;
  size_t len = strlen(r);
  char buf[BFAM_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, BFAM_LUA_MAX_COMMAND_LEN, BFAM_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, BFAM_LUA_EVALEXP_VAR);
    if (lua_isstring(L, -1))
      r = lua_tolstring(L, -1, &len);

    /* remove lua_getglobal value */
    lua_pop(L, 1);
  }

  char *s = bfam_malloc((len + 1) * sizeof(char));
  strncpy(s, r, len + 1);
  s[len] = '\0';

  return s;
}

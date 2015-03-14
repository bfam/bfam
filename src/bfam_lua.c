#include <bfam_base.h>
#include <bfam_lua.h>
#include <bfam_log.h>

#define BFAM_LUA_MAX_COMMAND_LEN 4096
#define BFAM_LUA_EVALEXP_VAR "XXX_bfam_evalexpr_XXX"

int bfam_lua_global_function_call(lua_State *L, const char *name,
                                  const char *sig, ...)
{
  va_list vl;
  int num_arg = 0;
  int num_res = 0;

  va_start(vl, sig);

  lua_getglobal(L, name);

  if (!lua_isfunction(L, -1))
  {
    BFAM_ROOT_WARNING("function `%s' not found in lua file", name);
    lua_pop(L, 1);
    return 1;
  }

  for (num_arg = 0; sig[num_arg] && sig[num_arg] != '>'; num_arg++)
  {
    luaL_checkstack(L, 1, "too many arguments");

    switch (sig[num_arg])
    {
    case 'd':
      lua_pushnumber(L, va_arg(vl, double));
      break;
    case 'l':
      lua_pushnumber(L, (double)va_arg(vl, bfam_long_real_t));
      break;
    case 'r':
      lua_pushnumber(L, (double)va_arg(vl, bfam_real_t));
      break;
    case 'i':
      lua_pushinteger(L, va_arg(vl, int));
      break;
    case 's':
      lua_pushstring(L, va_arg(vl, char *));
      break;
    case '>':
      break;
    default:
      BFAM_ABORT("function '%s' invalid input argument (%c)", name,
                 sig[num_arg]);
    }
  }

  BFAM_ABORT_IF_NOT(sig[num_arg] == '>', "arguments for '%s' does not contain "
                                         " a '>' character",
                    name);

  num_res = (int)strlen(sig) - num_arg - 1;

  BFAM_ABORT_IF_NOT(lua_pcall(L, num_arg, num_res, 0) == 0,
                    "error running function %s: %s", name, lua_tostring(L, -1));

  for (int n = 0; n < num_res; n++)
  {
    switch (sig[num_arg + 1 + n])
    {
    case 'd':
      BFAM_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                        "for '%s' return %d expected number got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, double *) = (double)lua_tonumber(L, n - num_res);
      break;
    case 'l':
      BFAM_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                        "for '%s' return %d expected number got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, bfam_long_real_t *) =
          (bfam_long_real_t)lua_tonumber(L, n - num_res);
      break;
    case 'r':
      BFAM_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                        "for '%s' return %d expected number got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, bfam_real_t *) = (bfam_real_t)lua_tonumber(L, n - num_res);
      break;
    case 'i':
      BFAM_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                        "for '%s' return %d expected number got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, int *) = (int)lua_tointeger(L, n - num_res);
      break;
    case 's':
      BFAM_ABORT_IF_NOT(lua_isstring(L, n - num_res),
                        "for '%s' return %d expected string got '%s'", name, n,
                        lua_tostring(L, n - num_res));
      *va_arg(vl, const char **) = lua_tostring(L, n - num_res);
      break;
    default:
      BFAM_ABORT("function '%s' invalid output argument (%c)", name,
                 sig[num_arg]);
    }
  }

  lua_pop(L, num_res);

  va_end(vl);
  return 0;
}

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

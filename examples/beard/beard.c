#include <bfam.h>

#ifdef BEARD_DGX_DIMENSION

#if BEARD_DGX_DIMENSION == 2
#include <bfam_domain_pxest_2.h>
#include "beard_dgx_rhs_2.h"
#include <p4est_iterate.h>
#elif BEARD_DGX_DIMENSION == 3
#include <bfam_domain_pxest_3.h>
#include "beard_dgx_rhs_3.h"
#include <p8est_iterate.h>
#define P4EST_CONNECT_CORNER P8EST_CONNECT_CORNER
#else
#error "bad dimension"
#endif

#if BEARD_DGX_DIMENSION == 2

#define DIM 2
#define FDIM 1
#define BEARD_D3_AP(A1, A2) (A1)
#define BEARD_D3_OP(A) BFAM_NOOP()
#define bfam_domain_pxest_init_callback bfam_domain_pxest_init_callback_2
#define bfam_domain_pxest_new bfam_domain_pxest_new_2
#define bfam_domain_pxest_quad_to_glueid bfam_domain_pxest_quad_to_glueid_2
#define bfam_domain_pxest_split_dgx_subdomains                                 \
  bfam_domain_pxest_split_dgx_subdomains_2
#define bfam_domain_pxest_free bfam_domain_pxest_free_2
#define bfam_domain_pxest_t bfam_domain_pxest_t_2

#elif BEARD_DGX_DIMENSION == 3

#define DIM 3
#define FDIM 2
#define BEARD_D3_AP(A1, A2) (A1 A2)
#define BEARD_D3_OP(A) A
#define bfam_domain_pxest_init_callback bfam_domain_pxest_init_callback_3
#define bfam_domain_pxest_new bfam_domain_pxest_new_3
#define bfam_domain_pxest_quad_to_glueid bfam_domain_pxest_quad_to_glueid_3
#define bfam_domain_pxest_split_dgx_subdomains                                 \
  bfam_domain_pxest_split_dgx_subdomains_3
#define bfam_domain_pxest_free bfam_domain_pxest_free_3
#define bfam_domain_pxest_t bfam_domain_pxest_t_3

#else
#error "bad dimension"
#endif

const char *comm_args_face_scalars[] = {NULL};
const char *comm_args_scalars[] = {NULL};
const char *comm_args_vectors[] = {"v", NULL};
const char *comm_args_vector_components[] = {"v1", "v2", "v3", NULL};
const char *comm_args_tensors[] = {"T", NULL};
const char *comm_args_tensor_components[] = {"S11", "S12", "S13", "S22",
                                             "S23", "S33", NULL};

const char *elastic_fields[] = {"rho", "lam", "mu",  "v1",  "v2",  "v3", "S11",
                                "S22", "S33", "S12", "S13", "S23", NULL};

const char *sponge_fields[] = {"a_sponge", NULL};

const char **plastic_fields = NULL;
const char *duvaut_lions_plastic_fields[] = {"c",     "Tr",    "phi",   "nu",
                                             "pf",    "S11_0", "S12_0", "S13_0",
                                             "S22_0", "S23_0", "S33_0", NULL};
#define DUVAUT_LIONS_TAG "_Duvaut-Lions"

const char **friction_fields = NULL;
const char *slip_weakening_fields[] = {
    "Tp1_0", "Tp2_0", "Tp3_0",  "Tn_0",     "Tp1",  "Tp2",   "Tp3",
    "Tn",    "V",     "Vp1",    "Vp2",      "Vp3",  "S11_0", "S12_0",
    "S13_0", "S22_0", "S23_0",  "S33_0",    "pf_0", "Dc",    "fs",
    "fc",    "fd",    "c0",     "Dp",       "Dp1",  "Dp2",   "Dp3",
    "Dn",    "Trup",  "Tforce", "Tforce_0", NULL};
const char *rate_and_state_fields[] = {
    "Tp1_0", "Tp2_0", "Tp3_0", "Tn_0", "Tp1",   "Tp2",   "Tp3",   "Tn",
    "V",     "Vp1",   "Vp2",   "Vp3",  "S11_0", "S12_0", "S13_0", "S22_0",
    "S23_0", "S33_0", "pf_0",  "f0",   "V0",    "a",     "b",     "L",
    "psi",   "Trup",  "Dp",    "Dp1",  "Dp2",   "Dp3",   "Dn",    NULL};

const char **friction_rates = NULL;
const char *slip_weakening_rates[] = {"Dp", "Dp1", "Dp2", "Dp3", "Dn", NULL};
const char *rate_and_state_rates[] = {"Dp", "Dp1", "Dp2", "Dp3",
                                      "Dn", "psi", NULL};

const char *volume_vtk_tags[] = {"elastic", "plastic", "sponge", NULL};

#define BFAM_LOAD_FIELD_RESTRICT_ALIGNED(field, prefix, base, dictionary)      \
  bfam_real_t *restrict field;                                                 \
  {                                                                            \
    char bfam_load_field_name[BFAM_BUFSIZ];                                    \
    snprintf(bfam_load_field_name, BFAM_BUFSIZ, "%s%s", (prefix), (base));     \
    field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);   \
    BFAM_ASSERT(field != NULL);                                                \
  }                                                                            \
  BFAM_ASSUME_ALIGNED(field, 32);
#define BFAM_LOAD_FIELD_ALIGNED(field, prefix, base, dictionary)               \
  bfam_real_t *field;                                                          \
  {                                                                            \
    char bfam_load_field_name[BFAM_BUFSIZ];                                    \
    snprintf(bfam_load_field_name, BFAM_BUFSIZ, "%s%s", (prefix), (base));     \
    field = bfam_dictionary_get_value_ptr(dictionary, bfam_load_field_name);   \
    BFAM_ASSERT(field != NULL);                                                \
  }                                                                            \
  BFAM_ASSUME_ALIGNED(field, 32);

#if 0
/*
 * These algorithms are taken from
 * http://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/trilinear.html
 *
 * Accessed: 12/18/13
 */
static void
inverse_linear(bfam_real_t *r,
    const bfam_real_t x,
    const bfam_real_t x0, const bfam_real_t x1)
{
  r[0] = (2*x - (x0+x1)) / (x1 - x0);

#ifdef BFAM_DEBUG
  if(BFAM_REAL_ABS(r[0]) <= 1)
  {
    BFAM_LDEBUG("Found point "
        "(%+10.5"BFAM_REAL_PRIe ") "
        "with (r[0]) = "
        "(%+10.5"BFAM_REAL_PRIe ")",
        x, r[0]);
      BFAM_LDEBUG("in cell "
          "(%+10.5"BFAM_REAL_PRIe ") "
          "(%+10.5"BFAM_REAL_PRIe ")",
          x0,x1);
  }
#endif
}
#endif

// Function that we will expose to lua
static int L_jn (lua_State *L)
{
    int     n = luaL_checknumber(L, 1);
    double  d  = luaL_checknumber(L, 2);
    lua_pushnumber(L, jn(n,d)); // Push the result
    return 1;   // number of results
}


static void inverse_bilinear_normal(bfam_real_t *interp_coeff,
                                    const bfam_real_t x, const bfam_real_t y,
                                    const bfam_real_t z, const bfam_real_t nx,
                                    const bfam_real_t ny, const bfam_real_t nz,
                                    const bfam_real_t x0, const bfam_real_t x1,
                                    const bfam_real_t x2, const bfam_real_t x3,
                                    const bfam_real_t y0, const bfam_real_t y1,
                                    const bfam_real_t y2, const bfam_real_t y3,
                                    const bfam_real_t z0, const bfam_real_t z1,
                                    const bfam_real_t z2, const bfam_real_t z3)
{
  bfam_real_t r = 0;
  bfam_real_t s = 0;
  bfam_real_t t = 0;

  /* until we are curved this should converge in one iteration */
  {
    bfam_long_real_t f[] = {
        -x + t * nx +
            BFAM_REAL(0.25) * ((1 - r) * (1 - s) * x0 + (1 + r) * (1 - s) * x1 +
                               (1 - r) * (1 + s) * x2 + (1 + r) * (1 + s) * x3),
        -y + t * ny +
            BFAM_REAL(0.25) * ((1 - r) * (1 - s) * y0 + (1 + r) * (1 - s) * y1 +
                               (1 - r) * (1 + s) * y2 + (1 + r) * (1 + s) * y3),
        -z + t * nz +
            BFAM_REAL(0.25) *
                ((1 - r) * (1 - s) * z0 + (1 + r) * (1 - s) * z1 +
                 (1 - r) * (1 + s) * z2 + (1 + r) * (1 + s) * z3)};

    bfam_long_real_t J[] = {
        /* first column: df/dr */
        BFAM_REAL(0.25) *
            (-(1 - s) * x0 + (1 - s) * x1 + -(1 + s) * x2 + (1 + s) * x3),
        BFAM_REAL(0.25) *
            (-(1 - s) * y0 + (1 - s) * y1 + -(1 + s) * y2 + (1 + s) * y3),
        BFAM_REAL(0.25) *
            (-(1 - s) * z0 + (1 - s) * z1 + -(1 + s) * z2 + (1 + s) * z3),
        /* second column: df/ds */
        BFAM_REAL(0.25) *
            (-(1 - r) * x0 - (1 + r) * x1 + (1 - r) * x2 + (1 + r) * x3),
        BFAM_REAL(0.25) *
            (-(1 - r) * y0 - (1 + r) * y1 + (1 - r) * y2 + (1 + r) * y3),
        BFAM_REAL(0.25) *
            (-(1 - r) * z0 - (1 + r) * z1 + (1 - r) * z2 + (1 + r) * z3),
        /* third column: df/dt */
        nx, ny, nz,
    };
    bfam_long_real_t delta[3];
    bfam_util_backslash(3, 1, J, f, delta);

    r = (bfam_real_t)(r - delta[0]);
    s = (bfam_real_t)(s - delta[1]);
    t = (bfam_real_t)(t - delta[2]);
  }

  interp_coeff[0] = r;
  interp_coeff[1] = s;
  interp_coeff[2] = t;
}

#if DIM == 2
static void inverse_bilinear(bfam_real_t *r, const bfam_real_t x,
                             const bfam_real_t y, const bfam_real_t x0,
                             const bfam_real_t x1, const bfam_real_t x2,
                             const bfam_real_t x3, const bfam_real_t y0,
                             const bfam_real_t y1, const bfam_real_t y2,
                             const bfam_real_t y3)
{
  const bfam_real_t f0 = x3 + x2 + x1 + x0 - BFAM_REAL(4) * x;
  const bfam_real_t f1 = x3 - x2 + x1 - x0;
  const bfam_real_t f2 = x3 + x2 - x1 - x0;
  const bfam_real_t f3 = x3 - x2 - x1 + x0;

  const bfam_real_t g0 = y3 + y2 + y1 + y0 - BFAM_REAL(4) * y;
  const bfam_real_t g1 = y3 - y2 + y1 - y0;
  const bfam_real_t g2 = y3 + y2 - y1 - y0;
  const bfam_real_t g3 = y3 - y2 - y1 + y0;

  bfam_real_t a = 0;
  bfam_real_t b = 0;

  /* until we are curved this should converge in one iteration */
  {
    const bfam_real_t f = f0 + a * f1 + b * f2 + a * b * f3;
    const bfam_real_t g = g0 + a * g1 + b * g2 + a * b * g3;

    const bfam_real_t da =
        (f * (g2 + a * g3) - g * (f2 + a * f3)) /
        (-(f1 + b * f3) * (g2 + a * g3) + (g1 + b * g3) * (f2 + a * f3));
    const bfam_real_t db =
        (f * (g1 + b * g3) - g * (f1 + b * f3)) /
        (-(f2 + a * f3) * (g1 + b * g3) + (g2 + a * g3) * (f1 + b * f3));

    a += da;
    b += db;
  }

  r[0] = a;
  r[1] = b;

#ifdef BFAM_DEBUG
  if (BFAM_REAL_ABS(a) <= 1 && BFAM_REAL_ABS(b) <= 1)
  {
    BFAM_LDEBUG("Found point "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe ") "
                "with (r[0],r[1]) = "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe ")",
                x, y, r[0], r[1]);
    BFAM_LDEBUG("in cell "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe ")",
                x0, y0, x1, y1, x2, y2, x3, y3);
  }
#endif
}
#endif

#if DIM == 3
static void inverse_trilinear(
    bfam_real_t *r, const bfam_real_t x, const bfam_real_t y,
    const bfam_real_t z, const bfam_real_t x0, const bfam_real_t x1,
    const bfam_real_t x2, const bfam_real_t x3, const bfam_real_t x4,
    const bfam_real_t x5, const bfam_real_t x6, const bfam_real_t x7,
    const bfam_real_t y0, const bfam_real_t y1, const bfam_real_t y2,
    const bfam_real_t y3, const bfam_real_t y4, const bfam_real_t y5,
    const bfam_real_t y6, const bfam_real_t y7, const bfam_real_t z0,
    const bfam_real_t z1, const bfam_real_t z2, const bfam_real_t z3,
    const bfam_real_t z4, const bfam_real_t z5, const bfam_real_t z6,
    const bfam_real_t z7)
{
  const bfam_real_t f0 =
      x7 + x6 + x5 + x4 + x3 + x2 + x1 + x0 - BFAM_REAL(8) * x;
  const bfam_real_t f1 = x7 - x6 + x5 - x4 + x3 - x2 + x1 - x0;
  const bfam_real_t f2 = x7 + x6 - x5 - x4 + x3 + x2 - x1 - x0;
  const bfam_real_t f3 = x7 + x6 + x5 + x4 - x3 - x2 - x1 - x0;
  const bfam_real_t f4 = x7 - x6 - x5 + x4 + x3 - x2 - x1 + x0;
  const bfam_real_t f5 = x7 - x6 + x5 - x4 - x3 + x2 - x1 + x0;
  const bfam_real_t f6 = x7 + x6 - x5 - x4 - x3 - x2 + x1 + x0;
  const bfam_real_t f7 = x7 - x6 - x5 + x4 - x3 + x2 + x1 - x0;

  const bfam_real_t g0 =
      y7 + y6 + y5 + y4 + y3 + y2 + y1 + y0 - BFAM_REAL(8) * y;
  const bfam_real_t g1 = y7 - y6 + y5 - y4 + y3 - y2 + y1 - y0;
  const bfam_real_t g2 = y7 + y6 - y5 - y4 + y3 + y2 - y1 - y0;
  const bfam_real_t g3 = y7 + y6 + y5 + y4 - y3 - y2 - y1 - y0;
  const bfam_real_t g4 = y7 - y6 - y5 + y4 + y3 - y2 - y1 + y0;
  const bfam_real_t g5 = y7 - y6 + y5 - y4 - y3 + y2 - y1 + y0;
  const bfam_real_t g6 = y7 + y6 - y5 - y4 - y3 - y2 + y1 + y0;
  const bfam_real_t g7 = y7 - y6 - y5 + y4 - y3 + y2 + y1 - y0;

  const bfam_real_t h0 =
      z7 + z6 + z5 + z4 + z3 + z2 + z1 + z0 - BFAM_REAL(8) * z;
  const bfam_real_t h1 = z7 - z6 + z5 - z4 + z3 - z2 + z1 - z0;
  const bfam_real_t h2 = z7 + z6 - z5 - z4 + z3 + z2 - z1 - z0;
  const bfam_real_t h3 = z7 + z6 + z5 + z4 - z3 - z2 - z1 - z0;
  const bfam_real_t h4 = z7 - z6 - z5 + z4 + z3 - z2 - z1 + z0;
  const bfam_real_t h5 = z7 - z6 + z5 - z4 - z3 + z2 - z1 + z0;
  const bfam_real_t h6 = z7 + z6 - z5 - z4 - z3 - z2 + z1 + z0;
  const bfam_real_t h7 = z7 - z6 - z5 + z4 - z3 + z2 + z1 - z0;

  bfam_real_t a = 0;
  bfam_real_t b = 0;
  bfam_real_t c = 0;

  /* until we are curved this should converge in one iteration */
  const int MAX_ITERATIONS = 200;
  int diff = 0;
  for (int n = 0; n < MAX_ITERATIONS; n++)
  {
    const bfam_real_t f = f0 + a * f1 + b * f2 + c * f3 + a * b * f4 +
                          a * c * f5 + b * c * f6 + a * b * c * f7;
    const bfam_real_t g = g0 + a * g1 + b * g2 + c * g3 + a * b * g4 +
                          a * c * g5 + b * c * g6 + a * b * c * g7;
    const bfam_real_t h = h0 + a * h1 + b * h2 + c * h3 + a * b * h4 +
                          a * c * h5 + b * c * h6 + a * b * c * h7;
    diff = !BFAM_REAL_APPROX_EQ(BFAM_REAL(0), f, 10);
    diff += !BFAM_REAL_APPROX_EQ(BFAM_REAL(0), g, 10);
    diff += !BFAM_REAL_APPROX_EQ(BFAM_REAL(0), h, 10);
    if (!diff)
      break;

    const bfam_real_t a_f = f1 + b * f4 + c * f5 + b * c * f7;
    const bfam_real_t b_f = f2 + a * f4 + c * f6 + a * c * f7;
    const bfam_real_t c_f = f3 + a * f5 + b * f6 + a * b * f7;

    const bfam_real_t a_g = g1 + b * g4 + c * g5 + b * c * g7;
    const bfam_real_t b_g = g2 + a * g4 + c * g6 + a * c * g7;
    const bfam_real_t c_g = g3 + a * g5 + b * g6 + a * b * g7;

    const bfam_real_t a_h = h1 + b * h4 + c * h5 + b * c * h7;
    const bfam_real_t b_h = h2 + a * h4 + c * h6 + a * c * h7;
    const bfam_real_t c_h = h3 + a * h5 + b * h6 + a * b * h7;

    const bfam_real_t da =
        (g * (b_f * c_h - b_h * c_f)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f) -
        (f * (b_g * c_h - b_h * c_g)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f) -
        (h * (b_f * c_g - b_g * c_f)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f);
    const bfam_real_t db =
        (f * (a_g * c_h - a_h * c_g)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f) -
        (g * (a_f * c_h - a_h * c_f)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f) +
        (h * (a_f * c_g - a_g * c_f)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f);
    const bfam_real_t dc =
        (g * (a_f * b_h - a_h * b_f)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f) -
        (f * (a_g * b_h - a_h * b_g)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f) -
        (h * (a_f * b_g - a_g * b_f)) /
            (a_f * b_g * c_h - a_f * b_h * c_g - a_g * b_f * c_h +
             a_g * b_h * c_f + a_h * b_f * c_g - a_h * b_g * c_f);

    a += da;
    b += db;
    c += dc;
  }

  if (diff > 0)
  {
    a = bfam_real_nan("");
    b = bfam_real_nan("");
    c = bfam_real_nan("");
  }

  r[0] = a;
  r[1] = b;
  r[2] = c;

#ifdef BFAM_DEBUG
  if (BFAM_REAL_ABS(a) <= 1 && BFAM_REAL_ABS(b) <= 1 && BFAM_REAL_ABS(c) <= 1)
  {
    BFAM_LDEBUG("Found point "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ") "
                "with (r[0],r[1],r[2]) = "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ")",
                x, y, z, r[0], r[1], r[2]);
    BFAM_LDEBUG("in cell "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ") "
                "(%+10.5" BFAM_REAL_PRIe ",%+10.5" BFAM_REAL_PRIe
                ",%+10.5" BFAM_REAL_PRIe ")",
                x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5,
                y5, z5, x6, y6, z6, x7, y7, z7);
  }
#endif
}
#endif

typedef struct beard
{
  MPI_Comm mpicomm;
  int mpirank;
  int mpisize;

  p4est_connectivity_t *conn;
  bfam_domain_pxest_t *domain;

  bfam_ts_t *beard_ts;
  bfam_subdomain_comm_args_t *comm_args;

  bfam_dictionary_t *volume_stations;
  bfam_dictionary_t *fault_stations;
} beard_t;

struct lsrk_table
{
  const char *name;
  bfam_ts_lsrk_method_t lsrk_method;
} lsrk_table[] = {
      {"KC54", BFAM_TS_LSRK_KC54},
      {"FE", BFAM_TS_LSRK_FE},
      {"HEUN", BFAM_TS_LSRK_HEUN},
      {"W33", BFAM_TS_LSRK_W33},
      {NULL, BFAM_TS_LSRK_NOOP},
};

struct adams_table
{
  const char *name;
  bfam_ts_adams_method_t adams_method;
} adams_table[] = {
      {"Adams 1", BFAM_TS_ADAMS_1},
      {"Adams 2", BFAM_TS_ADAMS_2},
      {"Adams 3", BFAM_TS_ADAMS_3},
      {"Adams 4", BFAM_TS_ADAMS_4},
      {NULL, BFAM_TS_ADAMS_NOOP},
};

struct local_adams_table
{
  const char *name;
  bfam_ts_local_adams_method_t local_adams_method;
} local_adams_table[] = {
      {"Adams 1", BFAM_TS_LOCAL_ADAMS_1},
      {"Adams 2", BFAM_TS_LOCAL_ADAMS_2},
      {"Adams 3", BFAM_TS_LOCAL_ADAMS_3},
      {"Adams 4", BFAM_TS_LOCAL_ADAMS_4},
      {NULL, BFAM_TS_LOCAL_ADAMS_NOOP},
};

typedef struct brick_args
{
  bfam_locidx_t nx;
  bfam_locidx_t ny;
  bfam_locidx_t nz;
  bfam_locidx_t periodic_x;
  bfam_locidx_t periodic_y;
  bfam_locidx_t periodic_z;
  bfam_locidx_t p4est_brick;
  bfam_locidx_t *brick_to_tree;
  bfam_locidx_t *bc;
} brick_args_t;

typedef enum plastic_type
{
  NONE,
  DUVAUT_LIONS,
} plastic_type_t;

typedef struct prefs
{
  lua_State *L;

  int dimension;

  char output_prefix[BFAM_BUFSIZ];
  char data_directory[BFAM_BUFSIZ];
  int vtk_binary;
  int vtk_compress;
  int vtk_num_pnts;
  int output_file_fault;
  char default_boundary_tag[BFAM_BUFSIZ];

  brick_args_t *brick_args;
  char mesh_file[BFAM_BUFSIZ];

  bfam_ts_lsrk_method_t lsrk_method;
  char lsrk_name[BFAM_BUFSIZ];

  bfam_ts_adams_method_t adams_method;
  char adams_name[BFAM_BUFSIZ];

  bfam_ts_local_adams_method_t local_adams_method;
  char local_adams_name[BFAM_BUFSIZ];

  plastic_type_t plasticity;
} prefs_t;

/*
 * Lua helper functions
 * \param [in]  L              pointer to lua file
 * \param [in]  ret_val_cont   if true then continue on return type failure
 * \param [in]  name           function to execure in lua
 * \param [in]  sig            signature of the call
 * \param [in,out] var_args    inputs and pointers for outputs
 */
/*
 * Helper function for calling a lua function. Based on generic call function of
 * Listing 25.4-25.6 of
 * @book{Ierusalimschy2006Lua,
 *  author = {Ierusalimschy, Roberto},
 *  title = {Programming in Lua, Second Edition},
 *  year = {2006},
 *  isbn = {8590379825},
 *  publisher = {Lua.Org},
 * }
 */
static int lua_global_function_call(lua_State *L, const int ret_val_cont,
                                    const char *name, const char *sig, ...)
{
  va_list vl;
  int num_arg = 0;
  int num_res = 0;
  int ret_val = 0;

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

  num_res = strlen(sig) - num_arg - 1;

  BFAM_ABORT_IF_NOT(lua_pcall(L, num_arg, num_res, 0) == 0,
                    "error running function %s: %s", name, lua_tostring(L, -1));

  for (int n = 0; n < num_res; n++)
  {
    switch (sig[num_arg + 1 + n])
    {
    case 'r':
      if (lua_isnumber(L, n - num_res))
        *va_arg(vl, bfam_real_t *) = (bfam_real_t)lua_tonumber(L, n - num_res);
      else
        BFAM_ABORT_IF_NOT(ret_val_cont && ++ret_val,
                          "for '%s' return %d expected number got '%s'", name,
                          n, lua_tostring(L, n - num_res));
      break;
    case 'i':
      if (lua_isnumber(L, n - num_res))
        *va_arg(vl, int *) = lua_tointeger(L, n - num_res);
      else
        BFAM_ABORT_IF_NOT(ret_val_cont && ++ret_val,
                          "   for '%s' return %d expected number got '%s'",
                          name, n, lua_tostring(L, n - num_res));
      break;
    case 's':
      if (lua_isstring(L, n - num_res))
        *va_arg(vl, const char **) = lua_tostring(L, n - num_res);
      else
        BFAM_ABORT_IF_NOT(ret_val_cont && ++ret_val,
                          "for '%s' return %d expected string got '%s'", name,
                          n, lua_tostring(L, n - num_res));
      break;
    default:
      BFAM_ABORT("function '%s' invalid output argument (%c)", name,
                 sig[num_arg]);
    }
  }

  lua_pop(L, num_res);

  va_end(vl);
  return ret_val;
}
static int lua_get_global_int(lua_State *L, const char *name, int def)
{
  lua_getglobal(L, name);
  int result = def;
  if (!lua_isnumber(L, -1))
    BFAM_ROOT_WARNING("`%s' not found, using default %d", name, def);
  else
    result = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);
  return result;
}

static int lua_get_table_int(lua_State *L, const char *table, const char *name,
                             int def)
{
  int result = def;
  lua_getglobal(L, table);
  if (!lua_istable(L, -1))
    BFAM_ROOT_WARNING("table `%s' not found, using default %d", table, def);
  else
  {
    lua_pushstring(L, name);
    lua_gettable(L, -2);
    if (!lua_isnumber(L, -1))
      BFAM_ROOT_WARNING("table `%s' does not contain `%s', using default %d",
                        table, name, def);
    else
      result = (int)lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
  return result;
}

/*
 * Set up the beard preference file
 */
static prefs_t *new_prefs(const char *prefs_filename)
{
  prefs_t *prefs = bfam_malloc(sizeof(prefs_t));

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  if (luaL_loadfile(L, prefs_filename) || lua_pcall(L, 0, 0, 0))
    BFAM_LERROR("cannot run configuration file: `%s'", lua_tostring(L, -1));

  prefs->L = L;
  lua_pushcfunction(L, L_jn);
  lua_setglobal(L, "jn");


  BFAM_ASSERT(lua_gettop(L) == 0);

  prefs->dimension = DIM;

  BFAM_ABORT_IF_NOT(prefs->dimension == 2 || prefs->dimension == 3,
                    "beard not set up to handle dimension %d",
                    prefs->dimension);

  /* get output prefix */
  lua_getglobal(L, "output_prefix");
  strncpy(prefs->output_prefix, "solution", BFAM_BUFSIZ);
  if (lua_isstring(L, -1))
  {
    const char *output_prefix = lua_tostring(L, -1);
    strncpy(prefs->output_prefix, output_prefix, BFAM_BUFSIZ);
  }
  else
    BFAM_ROOT_WARNING("using default 'output_prefix': %s",
                      prefs->output_prefix);
  lua_pop(L, 1);

  /* get output prefix */
  lua_getglobal(L, "data_directory");
  strncpy(prefs->data_directory, ".", BFAM_BUFSIZ);
  if (lua_isstring(L, -1))
  {
    const char *data_directory = lua_tostring(L, -1);
    strncpy(prefs->data_directory, data_directory, BFAM_BUFSIZ);
  }
  else
    BFAM_ROOT_WARNING("using default 'data_directory': %s",
                      prefs->data_directory);
  lua_pop(L, 1);

  /* get vtk_binary */
  lua_getglobal(L, "vtk_binary");
  if (lua_isboolean(L, -1))
  {
    prefs->vtk_binary = lua_toboolean(L, -1);
  }
  else
  {
    prefs->vtk_binary = 1;
    BFAM_ROOT_WARNING("using default 'vtk_binary': %d", prefs->vtk_binary);
  }
  lua_pop(L, 1);

  /* get vtk_compress */
  lua_getglobal(L, "vtk_compress");
  if (lua_isboolean(L, -1))
  {
    prefs->vtk_compress = lua_toboolean(L, -1);
  }
  else
  {
    prefs->vtk_compress = 1;
    BFAM_ROOT_WARNING("using default 'vtk_compress': %d", prefs->vtk_compress);
  }
  lua_pop(L, 1);

  /* get vtk_num_pnts */
  lua_getglobal(L, "vtk_num_pnts");
  if (lua_isnumber(L, -1))
  {
    prefs->vtk_num_pnts = lua_tointeger(L, -1);
  }
  else
  {
    prefs->vtk_num_pnts = 0;
    BFAM_ROOT_WARNING("using default 'vtk_num_pnts': %d", prefs->vtk_num_pnts);
  }
  lua_pop(L, 1);

  /* get output_file_fault */
  lua_getglobal(L, "output_file_fault");
  if (lua_isnumber(L, -1))
  {
    prefs->output_file_fault = lua_tointeger(L, -1);
  }
  else
  {
    prefs->output_file_fault = 0;
    BFAM_ROOT_WARNING("using default 'output_file_fault': %d",
                      prefs->output_file_fault);
  }
  lua_pop(L, 1);

  /* get default boundary tag */
  lua_getglobal(L, "default_boundary_tag");
  strncpy(prefs->default_boundary_tag, "non-reflecting", BFAM_BUFSIZ);
  if (lua_isstring(L, -1))
  {
    const char *default_boundary_tag = lua_tostring(L, -1);
    strncpy(prefs->default_boundary_tag, default_boundary_tag, BFAM_BUFSIZ);
  }
  else
    BFAM_ROOT_WARNING("using default 'default_boundary_tag': %s",
                      prefs->default_boundary_tag);
  lua_pop(L, 1);

  /* get the time stepper type: we look for both lsrk_method and adams_method */
  prefs->lsrk_method = BFAM_TS_LSRK_NOOP;
  lua_getglobal(L, "lsrk_method");
  if (lua_isstring(L, -1))
  {
    prefs->lsrk_method = lsrk_table[0].lsrk_method;
    strncpy(prefs->lsrk_name, lsrk_table[0].name, BFAM_BUFSIZ);
    if (lua_isstring(L, -1))
    {
      int i;
      const char *lsrk_name = lua_tostring(L, -1);
      for (i = 0; lsrk_table[i].name != NULL; ++i)
      {
        if (strcmp(lsrk_name, lsrk_table[i].name) == 0)
          break;
      }

      if (lsrk_table[i].name == NULL)
        BFAM_ROOT_WARNING("invalid lsrk method name: `%s'; using default %s",
                          lsrk_name, prefs->lsrk_name);
      else
      {
        prefs->lsrk_method = lsrk_table[i].lsrk_method;
        strncpy(prefs->lsrk_name, lsrk_table[i].name, BFAM_BUFSIZ);
      }
    }
    else
    {
      BFAM_ROOT_WARNING("`lsrk method' not found, using default: %s",
                        prefs->lsrk_name);
    }
    BFAM_ASSERT(prefs->lsrk_method != BFAM_TS_LSRK_NOOP);
  }
  lua_pop(L, 1);

  prefs->adams_method = BFAM_TS_ADAMS_NOOP;
  lua_getglobal(L, "adams_method");
  if (lua_isstring(L, -1))
  {
    prefs->adams_method = adams_table[0].adams_method;
    strncpy(prefs->adams_name, adams_table[0].name, BFAM_BUFSIZ);
    if (lua_isstring(L, -1))
    {
      int i;
      const char *adams_name = lua_tostring(L, -1);
      for (i = 0; adams_table[i].name != NULL; ++i)
      {
        if (strcmp(adams_name, adams_table[i].name) == 0)
          break;
      }

      if (adams_table[i].name == NULL)
        BFAM_ROOT_WARNING("invalid adams method name: `%s'; using default %s",
                          adams_name, prefs->adams_name);
      else
      {
        prefs->adams_method = adams_table[i].adams_method;
        strncpy(prefs->adams_name, adams_table[i].name, BFAM_BUFSIZ);
      }
    }
    else
    {
      BFAM_ROOT_WARNING("`adams method' not found, using default: %s",
                        prefs->adams_name);
    }
    BFAM_ASSERT(prefs->adams_method != BFAM_TS_ADAMS_NOOP);
  }
  lua_pop(L, 1);

  prefs->local_adams_method = BFAM_TS_LOCAL_ADAMS_NOOP;
  lua_getglobal(L, "local_adams_method");
  if (lua_isstring(L, -1))
  {
    prefs->local_adams_method = local_adams_table[0].local_adams_method;
    strncpy(prefs->local_adams_name, local_adams_table[0].name, BFAM_BUFSIZ);
    if (lua_isstring(L, -1))
    {
      int i;
      const char *local_adams_name = lua_tostring(L, -1);
      for (i = 0; local_adams_table[i].name != NULL; ++i)
      {
        if (strcmp(local_adams_name, local_adams_table[i].name) == 0)
          break;
      }

      if (local_adams_table[i].name == NULL)
        BFAM_ROOT_WARNING("invalid local adams method name: `%s';"
                          " using default %s",
                          local_adams_name, prefs->local_adams_name);
      else
      {
        prefs->local_adams_method = local_adams_table[i].local_adams_method;
        strncpy(prefs->local_adams_name, local_adams_table[i].name,
                BFAM_BUFSIZ);
      }
    }
    else
    {
      BFAM_ROOT_WARNING("`local_adams method' not found, using default: %s",
                        prefs->local_adams_name);
    }
    BFAM_ASSERT(prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP);
  }
  lua_pop(L, 1);

  BFAM_ABORT_IF_NOT(
      (prefs->lsrk_method != BFAM_TS_LSRK_NOOP &&
       prefs->adams_method == BFAM_TS_ADAMS_NOOP &&
       prefs->local_adams_method == BFAM_TS_LOCAL_ADAMS_NOOP) ||
          (prefs->lsrk_method == BFAM_TS_LSRK_NOOP &&
           prefs->adams_method != BFAM_TS_ADAMS_NOOP &&
           prefs->local_adams_method == BFAM_TS_LOCAL_ADAMS_NOOP) ||
          (prefs->lsrk_method == BFAM_TS_LSRK_NOOP &&
           prefs->adams_method == BFAM_TS_ADAMS_NOOP &&
           prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP),
      "must have either LSRK, ADAMS, or LOCAL time stepper");

  /* get the connectivity type */
  lua_getglobal(L, "mesh_file");
  strncpy(prefs->mesh_file, "", BFAM_BUFSIZ);
  prefs->brick_args = NULL;
  if (lua_isstring(L, -1))
  {
    const char *mesh_file = lua_tostring(L, -1);
    strncpy(prefs->mesh_file, mesh_file, BFAM_BUFSIZ);
    BFAM_INFO("Using mesh file: %s", mesh_file);
    lua_pop(L, 1);
  }
  else
  {
    lua_pop(L, 1);
    prefs->brick_args = bfam_malloc(sizeof(brick_args_t));

    prefs->brick_args->nx = lua_get_table_int(prefs->L, "brick", "nx", 1);
    prefs->brick_args->ny = lua_get_table_int(prefs->L, "brick", "ny", 1);

    prefs->brick_args->periodic_x =
        lua_get_table_int(prefs->L, "brick", "periodic_x", 0);
    prefs->brick_args->periodic_y =
        lua_get_table_int(prefs->L, "brick", "periodic_y", 0);

    prefs->brick_args->p4est_brick =
        lua_get_table_int(prefs->L, "brick", "p4est_brick", 0);

    /* Get the boundary conditions */
    prefs->brick_args->bc = bfam_malloc(sizeof(bfam_locidx_t) * 6);
    prefs->brick_args->bc[0] = lua_get_table_int(prefs->L, "brick", "bc0", -1);
    prefs->brick_args->bc[1] = lua_get_table_int(prefs->L, "brick", "bc1", -1);
    prefs->brick_args->bc[2] = lua_get_table_int(prefs->L, "brick", "bc2", -1);
    prefs->brick_args->bc[3] = lua_get_table_int(prefs->L, "brick", "bc3", -1);
    prefs->brick_args->bc[4] = lua_get_table_int(prefs->L, "brick", "bc4", -1);
    prefs->brick_args->bc[5] = lua_get_table_int(prefs->L, "brick", "bc5", -1);

    prefs->brick_args->periodic_z = 0;

#if DIM == 2
#elif DIM == 3
    prefs->brick_args->nz = lua_get_table_int(prefs->L, "brick", "nz", 1);
    prefs->brick_args->periodic_z =
        lua_get_table_int(prefs->L, "brick", "periodic_z", 0);
#else
#error "Bad dimension"
#endif
  }

  /* Read (if they exist) the plasticity parameters */
  lua_getglobal(L, "plastic");
  prefs->plasticity = NONE;
  if (lua_istable(L, -1))
  {
    luaL_checktype(L, -1, LUA_TTABLE);

    lua_pushstring(L, "tag");
    lua_gettable(L, -2);
    if (lua_isstring(L, -1))
    {
      const char *tag = lua_tostring(L, -1);
      size_t match_len = lua_strlen(L, -1);
      if (match_len == 12 && 0 == strncmp(tag, "Duvaut-Lions", match_len))
      {
        prefs->plasticity = DUVAUT_LIONS;
        BFAM_ABORT_IF(plastic_fields, "Plastic fields already set!");
        plastic_fields = duvaut_lions_plastic_fields;
      }
      else
        BFAM_ABORT("invalid plasticity type '%s'", tag);
      lua_pop(L, 1);
    }
    else
      BFAM_ABORT("not 'tag' found in the 'plastic' table");
  }
  else
    BFAM_ROOT_WARNING("table `plastic' not found");
  lua_pop(L, 1);

  BFAM_ASSERT(lua_gettop(prefs->L) == 0);
  return prefs;
}

static void print_prefs(prefs_t *prefs)
{
  BFAM_ROOT_INFO("----------Preferences----------");
  BFAM_ROOT_INFO(" Beard Dimension          = %d", prefs->dimension);
  BFAM_ROOT_INFO(" Output prefix            = %s", prefs->output_prefix);
  BFAM_ROOT_INFO(" Default boundary tag     = %s", prefs->default_boundary_tag);
  if (prefs->lsrk_method != BFAM_TS_LSRK_NOOP)
    BFAM_ROOT_INFO(" Low Storage Time Stepper = %s", prefs->lsrk_name);
  if (prefs->adams_method != BFAM_TS_ADAMS_NOOP)
    BFAM_ROOT_INFO(" Adams-Bashforth Scheme   = %s", prefs->adams_name);
  if (prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
    BFAM_ROOT_INFO(" Local Adams-Bashforth Scheme   = %s",
                   prefs->local_adams_name);
  if (prefs->brick_args != NULL)
  {
    BFAM_ROOT_INFO(" brick arguments");
    BFAM_ROOT_INFO("  nx         = %d", prefs->brick_args->nx);
    BFAM_ROOT_INFO("  ny         = %d", prefs->brick_args->ny);
#if DIM == 2
#elif DIM == 3
    BFAM_ROOT_INFO("  nz         = %d", prefs->brick_args->nz);
#else
#error "Bad dimension"
#endif
    BFAM_ROOT_INFO("  periodic_x = %d", prefs->brick_args->periodic_x);
    BFAM_ROOT_INFO("  periodic_y = %d", prefs->brick_args->periodic_y);
#if DIM == 2
#elif DIM == 3
    BFAM_ROOT_INFO("  periodic_z = %d", prefs->brick_args->periodic_z);
#else
#error "Bad dimension"
#endif
    BFAM_ROOT_INFO("  p4est_brick = %d", prefs->brick_args->p4est_brick);
  }
  BFAM_ROOT_INFO("-------------------------------");
}

static void free_prefs(prefs_t *prefs)
{
  if (prefs->brick_args != NULL)
  {
    if (prefs->brick_args->bc != NULL)
      bfam_free(prefs->brick_args->bc);
    if (prefs->brick_args->brick_to_tree != NULL)
      bfam_free(prefs->brick_args->brick_to_tree);
    bfam_free(prefs->brick_args);
  }
  lua_close(prefs->L);
}

static void init_mpi(beard_t *beard, MPI_Comm mpicomm)
{
  beard->mpicomm = mpicomm;
  BFAM_MPI_CHECK(MPI_Comm_rank(mpicomm, &beard->mpirank));
  BFAM_MPI_CHECK(MPI_Comm_size(mpicomm, &beard->mpisize));
}

static int static_refine_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                            p4est_quadrant_t *quadrant)
{
  return 1;
}

static int refine_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t *quadrant)
{
  lua_State *L = (lua_State *)p4est->user_pointer;

  double vxyz[P4EST_CHILDREN * 3];

#if DIM == 2
  for (int iy = 0; iy < 2; iy++)
    for (int ix = 0; ix < 2; ix++)
    {
      int ox = ix * (1 << (P4EST_MAXLEVEL - quadrant->level));
      int oy = iy * (1 << (P4EST_MAXLEVEL - quadrant->level));
      p4est_qcoord_to_vertex(p4est->connectivity, which_tree, quadrant->x + ox,
                             quadrant->y + oy, &vxyz[3 * (ix + iy * 2)]);
    }

  int val = 0;
  int result = lua_global_function_call(
      L, 0, "refinement_function", "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ii"
                                   ">"
                                   "i",
      vxyz[0], vxyz[1], vxyz[2], vxyz[3], vxyz[4], vxyz[5], vxyz[6], vxyz[7],
      vxyz[8], vxyz[9], vxyz[10], vxyz[11], quadrant->level, which_tree, &val);
  BFAM_ABORT_IF_NOT(result == 0, "Expected 0 got %d", result);
#elif DIM == 3
  for (int iz = 0; iz < 2; iz++)
    for (int iy = 0; iy < 2; iy++)
      for (int ix = 0; ix < 2; ix++)
      {
        int ox = ix * (1 << (P4EST_MAXLEVEL - quadrant->level));
        int oy = iy * (1 << (P4EST_MAXLEVEL - quadrant->level));
        int oz = iz * (1 << (P4EST_MAXLEVEL - quadrant->level));
        p4est_qcoord_to_vertex(
            p4est->connectivity, which_tree, quadrant->x + ox, quadrant->y + oy,
            quadrant->z + oz, &vxyz[3 * (ix + 2 * (iy + 2 * iz))]);
      }

  int val = 0;
  int result = lua_global_function_call(
      L, 0, "refinement_function", "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ii"
                                   ">"
                                   "i",
      vxyz[0], vxyz[1], vxyz[2], vxyz[3], vxyz[4], vxyz[5], vxyz[6], vxyz[7],
      vxyz[8], vxyz[9], vxyz[10], vxyz[11], vxyz[12], vxyz[13], vxyz[14],
      vxyz[15], vxyz[16], vxyz[17], vxyz[18], vxyz[19], vxyz[20], vxyz[21],
      vxyz[22], vxyz[23], quadrant->level, which_tree, &val);
  BFAM_ABORT_IF_NOT(result == 0, "Expected 0 got %d", result);
#else
#error "Bad Dimension"
#endif

  return val;
}

typedef struct split_iter_data
{
  lua_State *L;
  bfam_locidx_t next_sub;
  bfam_locidx_t next_root;
  bfam_locidx_t loc;
  bfam_locidx_t *subdomain_id;
  bfam_dictionary_t *root_N;
  bfam_dictionary_t *root;
} split_iter_data_t;

static void get_element_order(p4est_iter_volume_info_t *info, void *arg)
{
  split_iter_data_t *data = (split_iter_data_t *)arg;

  double vxyz[P4EST_CHILDREN * 3];

  int N = 0;
  const char *root_str = "elastic";
#if DIM == 2
  for (int ix = 0; ix < 2; ix++)
    for (int iy = 0; iy < 2; iy++)
    {
      int ox = ix * (1 << (P4EST_MAXLEVEL - info->quad->level));
      int oy = iy * (1 << (P4EST_MAXLEVEL - info->quad->level));
      p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid,
                             info->quad->x + ox, info->quad->y + oy,
                             &vxyz[3 * (ix + 2 * iy)]);
    }

  int result = lua_global_function_call(
      data->L, 1, "element_order", "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ii"
                                   ">"
                                   "is",
      vxyz[0], vxyz[1], vxyz[2], vxyz[3], vxyz[4], vxyz[5], vxyz[6], vxyz[7],
      vxyz[8], vxyz[9], vxyz[10], vxyz[11], info->quad->level, info->treeid, &N,
      &root_str);
#elif DIM == 3
  for (int iz = 0; iz < 2; iz++)
    for (int iy = 0; iy < 2; iy++)
      for (int ix = 0; ix < 2; ix++)
      {
        int ox = ix * (1 << (P4EST_MAXLEVEL - info->quad->level));
        int oy = iy * (1 << (P4EST_MAXLEVEL - info->quad->level));
        int oz = iz * (1 << (P4EST_MAXLEVEL - info->quad->level));
        p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid,
                               info->quad->x + ox, info->quad->y + oy,
                               info->quad->z + oz,
                               &vxyz[3 * (ix + 2 * (iy + 2 * iz))]);
      }

  int result = lua_global_function_call(
      data->L, 0, "element_order", "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ddd"
                                   "ii"
                                   ">"
                                   "is",
      vxyz[0], vxyz[1], vxyz[2], vxyz[3], vxyz[4], vxyz[5], vxyz[6], vxyz[7],
      vxyz[8], vxyz[9], vxyz[10], vxyz[11], vxyz[12], vxyz[13], vxyz[14],
      vxyz[15], vxyz[16], vxyz[17], vxyz[18], vxyz[19], vxyz[20], vxyz[21],
      vxyz[22], vxyz[23], info->quad->level, info->treeid, &N, &root_str);
#else
#error "Bad Dimension"
#endif

  BFAM_ABORT_IF_NOT(result <= 1, "something went wrong with 'element_order'"
                                 " function");
  BFAM_ABORT_IF_NOT(N > 0, "order '%d' is not bigger than zero", N);

  bfam_locidx_t sub;
  bfam_locidx_t root;
  if (0 == bfam_dictionary_get_value_locidx(data->root, root_str, &root))
  {
    root = data->next_root;
    data->next_root++;
    bfam_dictionary_insert_locidx(data->root, root_str, root);
  }

  char order_str[BFAM_BUFSIZ];
  snprintf(order_str, BFAM_BUFSIZ, "%d_%d", root, N);
  if (0 == bfam_dictionary_get_value_locidx(data->root_N, order_str, &sub))
  {
    sub = data->next_sub;
    data->next_sub++;
    bfam_dictionary_insert_locidx(data->root_N, order_str, sub);
  }

  data->subdomain_id[data->loc] = sub;
  data->loc++;
}

/*
 * Read a line from a file. Obtained from:
 * http://stackoverflow.com/questions/314401/
 * how-to-read-a-line-from-the-console-in-c/314422#314422
 *
 * Using this avoids a dependence on IEEE Std 1003.1-2008 (``POSIX.1'') for the
 * getline function.
 */
static char *beard_getline_upper(FILE *stream)
{
  char *line = bfam_malloc(sizeof(char) * 1024), *linep = line;
  size_t lenmax = 1024, len = lenmax;
  int c;

  if (line == NULL)
    return NULL;

  for (;;)
  {
    c = fgetc(stream);
    c = toupper(c);
    if (c == EOF && linep == line)
    {
      bfam_free(linep);
      return NULL;
    }

    if (--len == 0)
    {
      len = lenmax;
      lenmax *= 2;
      char *linen = bfam_realloc(linep, sizeof(char) * lenmax);

      if (linen == NULL)
      {
        bfam_free(linep);
        return NULL;
      }
      line = linen + (line - linep);
      linep = linen;
    }
    if ((*line++ = (char)c) == '\n')
      break;
  }
  *line = '\0';
  return linep;
}

static int glueid_read_file(FILE *fid, bfam_locidx_t *tree_to_glueid,
                            const p4est_connectivity_t *conn)
{
  int lines_read = 0, lines_free = 0;
  char *line;
  int reading_side_set = 0;
  int glueid = -1;
  int faceid = -1;
  int abq_to_p4est[] = {4, 5, 2, 1, 3, 0};

  for (;;)
  {
    line = beard_getline_upper(fid);

    if (line == NULL)
    {
      break;
    }

    ++lines_read;
    /* check for control line */
    if (line[0] == '*')
    {
      reading_side_set = 0;
      glueid = -1;
      faceid = -1;

      if (strstr(line, "*ELSET"))
      {
        reading_side_set = 1;

        int retval =
            sscanf(line, "%*[^0-9]%d%*[^0-9]%d%*[^0-9]", &glueid, &faceid);
        BFAM_ABORT_IF_NOT(retval == 2, "wrong retval in glueid_read_file");

        /* convert to the p4est faceid */
        faceid = abq_to_p4est[faceid - 1];

        ++lines_free;
        bfam_free(line);
        continue;
      }
    }

    if (reading_side_set)
    {
      long long int treeid;
      int retval = 1;

      char *sv_ptr;
      char *toks = strtok_r(line, ", ", &sv_ptr);

      while (toks != NULL)
      {
        retval = sscanf(toks, "%lld", &treeid);
        if (retval <= 0)
          break;

        /* ABQ starts at 1 we start at zero */
        treeid--;
        tree_to_glueid[P4EST_FACES * treeid + faceid] = glueid;

        bfam_locidx_t ntreeid =
            conn->tree_to_tree[treeid * P4EST_FACES + faceid];

        /* tree_to_face contains the orientation as well hence mod */
        bfam_locidx_t nfaceid =
            conn->tree_to_face[treeid * P4EST_FACES + faceid] % P4EST_FACES;

        if (ntreeid != treeid)
        {
          tree_to_glueid[P4EST_FACES * ntreeid + nfaceid] = glueid;
        }

        toks = strtok_r(NULL, ", ", &sv_ptr);
      }
    }

    ++lines_free;
    bfam_free(line);
  }

  BFAM_LDEBUG("glueid_read_file: read %d / freed %d lines", lines_read,
              lines_free);
  return lines_read - lines_free;
}

static void init_tree_to_glueid(beard_t *beard, prefs_t *prefs,
                                bfam_locidx_t *tree_to_glueid)
{
  const p4est_connectivity_t *conn = beard->domain->pxest->connectivity;
  bfam_locidx_t num_tree_ids = P4EST_FACES * conn->num_trees;
  for (int k = 0; k < num_tree_ids; k++)
    tree_to_glueid[k] = -1;

  /* First we set the brick boundary conditions (if they exist) */
  brick_args_t *brick_args = prefs->brick_args;
  if (brick_args)
  {
    bfam_locidx_t faceid = 0;
    if (brick_args->bc[faceid] > 0)
    {
      BFAM_ABORT_IF_NOT(brick_args->brick_to_tree,
                        "brick to tree not defined"
                        " (this is likely due to 'p4est_brick' being set)");

      const bfam_locidx_t ix = 0;
#if DIM == 3
      for (bfam_locidx_t iz = 0; iz < brick_args->nz; iz++)
#endif
        for (bfam_locidx_t iy = 0; iy < brick_args->ny; iy++)
        {
          const bfam_locidx_t treeid = brick_args->brick_to_tree[BEARD_D3_AP(
              ix + brick_args->nx * iy, +brick_args->nx * brick_args->ny * iz)];
          tree_to_glueid[P4EST_FACES * treeid + faceid] =
              brick_args->bc[faceid];
        }
    }

    faceid = 1;
    if (brick_args->bc[faceid] > 0)
    {
      BFAM_ABORT_IF_NOT(brick_args->brick_to_tree,
                        "brick to tree not defined"
                        " (this is likely due to 'p4est_brick' being set)");

      const bfam_locidx_t ix = brick_args->nx - 1;
#if DIM == 3
      for (bfam_locidx_t iz = 0; iz < brick_args->nz; iz++)
#endif
        for (bfam_locidx_t iy = 0; iy < brick_args->ny; iy++)
        {
          const bfam_locidx_t treeid = brick_args->brick_to_tree[BEARD_D3_AP(
              ix + brick_args->nx * iy, +brick_args->nx * brick_args->ny * iz)];
          tree_to_glueid[P4EST_FACES * treeid + faceid] =
              brick_args->bc[faceid];
        }
    }

    faceid = 2;
    if (brick_args->bc[faceid] > 0)
    {
      BFAM_ABORT_IF_NOT(brick_args->brick_to_tree,
                        "brick to tree not defined"
                        " (this is likely due to 'p4est_brick' being set)");

      const bfam_locidx_t iy = 0;
#if DIM == 3
      for (bfam_locidx_t iz = 0; iz < brick_args->nz; iz++)
#endif
        for (bfam_locidx_t ix = 0; ix < brick_args->nx; ix++)
        {
          const bfam_locidx_t treeid = brick_args->brick_to_tree[BEARD_D3_AP(
              ix + brick_args->nx * iy, +brick_args->nx * brick_args->ny * iz)];
          tree_to_glueid[P4EST_FACES * treeid + faceid] =
              brick_args->bc[faceid];
        }
    }

    faceid = 3;
    if (brick_args->bc[faceid] > 0)
    {
      BFAM_ABORT_IF_NOT(brick_args->brick_to_tree,
                        "brick to tree not defined"
                        " (this is likely due to 'p4est_brick' being set)");

      const bfam_locidx_t iy = brick_args->ny - 1;
#if DIM == 3
      for (bfam_locidx_t iz = 0; iz < brick_args->nz; iz++)
#endif
        for (bfam_locidx_t ix = 0; ix < brick_args->nx; ix++)
        {
          const bfam_locidx_t treeid = brick_args->brick_to_tree[BEARD_D3_AP(
              ix + brick_args->nx * iy, +brick_args->nx * brick_args->ny * iz)];
          tree_to_glueid[P4EST_FACES * treeid + faceid] =
              brick_args->bc[faceid];
        }
    }

#if DIM == 3
    faceid = 4;
    if (brick_args->bc[faceid] > 0)
    {
      BFAM_ABORT_IF_NOT(brick_args->brick_to_tree,
                        "brick to tree not defined"
                        " (this is likely due to 'p4est_brick' being set)");

      const bfam_locidx_t iz = 0;
      for (bfam_locidx_t iy = 0; iy < brick_args->ny; iy++)
        for (bfam_locidx_t ix = 0; ix < brick_args->nx; ix++)
        {
          const bfam_locidx_t treeid = brick_args->brick_to_tree[BEARD_D3_AP(
              ix + brick_args->nx * iy, +brick_args->nx * brick_args->ny * iz)];
          tree_to_glueid[P4EST_FACES * treeid + faceid] =
              brick_args->bc[faceid];
        }
    }

    faceid = 5;
    if (brick_args->bc[faceid] > 0)
    {
      BFAM_ABORT_IF_NOT(brick_args->brick_to_tree,
                        "brick to tree not defined"
                        " (this is likely due to 'p4est_brick' being set)");

      const bfam_locidx_t iz = brick_args->nz - 1;
      for (bfam_locidx_t iy = 0; iy < brick_args->ny; iy++)
        for (bfam_locidx_t ix = 0; ix < brick_args->nx; ix++)
        {
          const bfam_locidx_t treeid = brick_args->brick_to_tree[BEARD_D3_AP(
              ix + brick_args->nx * iy, +brick_args->nx * brick_args->ny * iz)];
          tree_to_glueid[P4EST_FACES * treeid + faceid] =
              brick_args->bc[faceid];
        }
    }
#endif
  }

  /* Now we read the table from the Lua */
  lua_State *L = prefs->L;
#ifdef BFAM_DEBUG
  int top = lua_gettop(L);
#endif

  lua_getglobal(L, "glueid_treeid_faceid");

  if (lua_isstring(L, -1))
  {
    const char *tree_file = lua_tostring(L, -1);
    FILE *fid = fopen(tree_file, "r");
    BFAM_ABORT_IF(fid == NULL, "unable to read glueid file %s", tree_file);

    int retval = glueid_read_file(fid, tree_to_glueid, conn);
    BFAM_ABORT_IF(retval, "'glueid_read_file' returned %d on file %s", retval,
                  tree_file);

    fclose(fid);
  }
  else if (lua_istable(L, -1))
  {
    int n = lua_objlen(L, -1);
    BFAM_LDEBUG("glueid_treeid_faceid  #elem: %3d", n);

    BFAM_ABORT_IF_NOT(
        n % 3 == 0,
        "length of glueid_treeid_faceid should be a multiple of three");

    luaL_checktype(L, -1, LUA_TTABLE);

    for (int i = 1; i <= n; i += 3)
    {
      lua_rawgeti(L, -1, i + 0);
      int glueid = (int)lua_tointeger(L, -1);
      lua_pop(L, 1);

      lua_rawgeti(L, -1, i + 1);
      int treeid = (int)lua_tointeger(L, -1);
      lua_pop(L, 1);

      lua_rawgeti(L, -1, i + 2);
      int faceid = (int)lua_tointeger(L, -1);
      lua_pop(L, 1);

      BFAM_ABORT_IF(treeid < 0 || treeid > conn->num_trees,
                    "glueid_treeid_faceid: invalid tree id %d", treeid);

      BFAM_ABORT_IF(faceid < 0 || faceid > P4EST_FACES,
                    "glueid_treeid_faceid: invalid face id %d", faceid);
      if (prefs->brick_args && prefs->brick_args->brick_to_tree)
      {
        BFAM_LDEBUG("brick %04" BFAM_LOCIDX_PRId
                    " <--> tree %04" BFAM_LOCIDX_PRId,
                    treeid, prefs->brick_args->brick_to_tree[treeid]);
        treeid = prefs->brick_args->brick_to_tree[treeid];
      }

      if (glueid > 0)
      {
        tree_to_glueid[P4EST_FACES * treeid + faceid] = glueid;
      }
      else if (glueid < 0) /* negative glueid signifies to add neighbor */
      {

        bfam_locidx_t ntreeid =
            conn->tree_to_tree[treeid * P4EST_FACES + faceid];

        /* tree_to_face contains the orientation as well hence mod */
        bfam_locidx_t nfaceid =
            conn->tree_to_face[treeid * P4EST_FACES + faceid] % P4EST_FACES;

        tree_to_glueid[P4EST_FACES * treeid + faceid] = -glueid;
        if (ntreeid != treeid)
        {
          tree_to_glueid[P4EST_FACES * ntreeid + nfaceid] = -glueid;
        }
      }
    }
  }
  else
    BFAM_ROOT_WARNING("table or function `%s' not found",
                      "glueid_treeid_faceid");

  lua_pop(L, 1);
  BFAM_ASSERT(top == lua_gettop(L));
}

static int glue_order(const int N_m, const int N_p, const bfam_locidx_t uid,
                      void *args)
{
  prefs_t *prefs = (prefs_t *)args;
  lua_State *L = prefs->L;
  int N_g;
  lua_global_function_call(L, 0, "glue_order", "iii>i", N_m, N_p, (int)uid,
                           &N_g);
  return N_g;
}

static void transform_nodes(const bfam_locidx_t num_Vi,
                            const bfam_locidx_t num_pnts,
                            bfam_long_real_t **lxi, void *args)
{
  prefs_t *prefs = (prefs_t *)args;
  BFAM_ASSERT(num_Vi == prefs->dimension);

  lua_State *L = prefs->L;
  for (bfam_locidx_t k = 0; k < num_pnts; k++)
  {
    bfam_real_t x = (bfam_real_t)lxi[0][k];
    bfam_real_t y = (bfam_real_t)lxi[1][k];
#if DIM == 2
    lua_global_function_call(L, 0, "transform_nodes", "rr>rr", x, y, &x, &y);
#elif DIM == 3
    bfam_real_t z = (bfam_real_t)lxi[2][k];
    lua_global_function_call(L, 0, "transform_nodes", "rrr>rrr", x, y, z, &x,
                             &y, &z);
    lxi[2][k] = z;
#else
#error "Bad Dimension"
#endif

    lxi[0][k] = x;
    lxi[1][k] = y;
  }
}

typedef struct d2la_data
{
  bfam_locidx_t sz;
  bfam_locidx_t *A;
  bfam_locidx_t *B;
} d2la_data_t;

static int dict_to_locidx_array(const char *key_str, const char *val_str,
                                void *in_args)
{
  d2la_data_t *data = (d2la_data_t *)in_args;
  bfam_locidx_t V1;
  bfam_locidx_t V2;
  bfam_locidx_t k;

  int m = sscanf(key_str, "%d_%d", &V1, &V2);
  BFAM_ASSERT(m == 2);
  int n = sscanf(val_str, "%" BFAM_LOCIDX_SCNd, &k);
  BFAM_ASSERT(n == 1);

  BFAM_ASSERT(data->sz > k);
  data->A[k] = V1;
  data->B[k] = V2;

  return (n == 1) && (m == 2);
}

static int dict_to_volume_tag(const char *key_str, const char *val_str,
                              void *in_args)
{
  bfam_domain_t *dom = (bfam_domain_t *)in_args;
  bfam_locidx_t root;

  int n = sscanf(val_str, "%" BFAM_LOCIDX_SCNd, &root);
  BFAM_ASSERT(n == 1);
  char tag[BFAM_BUFSIZ];
  snprintf(tag, BFAM_BUFSIZ, "_volume_id_%jd", (intmax_t)root);
  const char *tags[] = {tag, NULL};
  bfam_domain_add_tag(dom, BFAM_DOMAIN_OR, tags, key_str);

  return (n == 1);
}

static void split_domain(beard_t *beard, prefs_t *prefs)
{
  bfam_domain_pxest_t *domain = beard->domain;
  if (domain->pxest->local_num_quadrants == 0)
    BFAM_ABORT("No quadrants");
  bfam_locidx_t *sub_ids =
      bfam_malloc(domain->pxest->local_num_quadrants * sizeof(bfam_locidx_t));

  bfam_dictionary_t root_N_dict;
  bfam_dictionary_t root_dict;
  bfam_dictionary_init(&root_N_dict);
  bfam_dictionary_init(&root_dict);
  split_iter_data_t data = {prefs->L, 0,            0,         0,
                            sub_ids,  &root_N_dict, &root_dict};
  p4est_iterate(domain->pxest, NULL, &data, get_element_order, NULL, NULL
#if DIM == 3
                ,
                NULL
#endif
                );

  bfam_locidx_t num_sub = root_N_dict.num_entries;
  bfam_locidx_t N[num_sub];
  bfam_locidx_t root[num_sub];

  d2la_data_t RN_data = {num_sub, root, N};
  bfam_dictionary_allprefixed(&root_N_dict, "", dict_to_locidx_array, &RN_data);

  bfam_locidx_t *tree_ids =
      bfam_malloc(P4EST_FACES * domain->pxest->connectivity->num_trees *
                  sizeof(bfam_locidx_t));

  init_tree_to_glueid(beard, prefs, tree_ids);

  bfam_locidx_t *glue_ids = bfam_malloc(
      P4EST_FACES * domain->pxest->local_num_quadrants * sizeof(bfam_locidx_t));

  bfam_domain_pxest_quad_to_glueid(domain->pxest, tree_ids, glue_ids);

  BFAM_ASSERT(lua_gettop(prefs->L) == 0);

  lua_getglobal(prefs->L, "transform_nodes");
  bfam_dgx_nodes_transform_t nodes_transform_fptr = NULL;
  void *nt_user_args = NULL;
  if (lua_isfunction(prefs->L, -1))
  {
    nodes_transform_fptr = &transform_nodes;
    nt_user_args = prefs;
  }
  else
    BFAM_ROOT_WARNING("function `%s' not found in lua file", "transform_nodes");
  lua_pop(prefs->L, 1);

  lua_getglobal(prefs->L, "glue_order");
  bfam_glue_order_t glue_order_fptr = NULL;
  void *go_user_args = NULL;
  if (lua_isfunction(prefs->L, -1))
  {
    glue_order_fptr = &glue_order;
    go_user_args = prefs;
  }
  else
    BFAM_ROOT_WARNING("function `%s' not found in lua file", "glue_order");
  lua_pop(prefs->L, 1);

  bfam_domain_pxest_split_dgx_subdomains(
      domain, num_sub, sub_ids, root, N, glue_ids, nodes_transform_fptr,
      nt_user_args, glue_order_fptr, go_user_args);
  BFAM_ASSERT(lua_gettop(prefs->L) == 0);

  bfam_dictionary_allprefixed(&root_dict, "", dict_to_volume_tag, domain);

  bfam_free(tree_ids);
  bfam_free(glue_ids);
  bfam_free(sub_ids);
  bfam_dictionary_clear(&root_N_dict);
  bfam_dictionary_clear(&root_dict);
}

typedef struct set_val_extended_args
{
  const char *fname;
  lua_State *L;
} set_val_extended_args_t;

static void field_set_val_extend(bfam_locidx_t npoints, const char *name,
                                 bfam_real_t t, bfam_real_t *restrict x,
                                 bfam_real_t *restrict y,
                                 bfam_real_t *restrict z,
                                 struct bfam_subdomain *s, void *arg,
                                 bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);
  bfam_real_t val = 0;

  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)s;
  BFAM_ASSERT(npoints == sub->Np * sub->K);

  set_val_extended_args_t *val_args = (set_val_extended_args_t *)arg;
  lua_State *L = val_args->L;
  lua_getglobal(L, val_args->fname);

  if (lua_isfunction(L, -1))
  {
    BFAM_ROOT_INFO("field_set_val: using '%s' as lua callback function",
                   val_args->fname);
    lua_pop(L, 1);
#if DIM == 2
    const bfam_real_t tmpz = 0;
#endif
    for (bfam_locidx_t elem = 0; elem < sub->K; ++elem)
    {
      const bfam_locidx_t off = elem * sub->Np;
      bfam_real_t xc = 0;
      bfam_real_t yc = 0;
      bfam_real_t zc = 0;
      switch (sub->dim)
      {
      case 1:
        for (bfam_locidx_t ix = 0; ix < 2; ix++)
        {
          const bfam_locidx_t ind = +ix * (sub->N);
          xc += x[off + ind];
          yc += y[off + ind];
#if DIM == 3
          zc += z[off + ind];
#endif
        }
        xc /= 2;
        yc /= 2;
        zc /= 2;
        break;
      case 2:
        for (bfam_locidx_t iy = 0; iy < 2; iy++)
          for (bfam_locidx_t ix = 0; ix < 2; ix++)
          {
            const bfam_locidx_t ind =
                +ix * (sub->N) + iy * (sub->N) * (sub->N + 1);

            xc += x[off + ind];
            yc += y[off + ind];
#if DIM == 3
            zc += z[off + ind];
#endif
          }
        xc /= 4;
        yc /= 4;
        zc /= 4;
        break;
      case 3:
        for (bfam_locidx_t iz = 0; iz < 2; iz++)
          for (bfam_locidx_t iy = 0; iy < 2; iy++)
            for (bfam_locidx_t ix = 0; ix < 2; ix++)
            {
              const bfam_locidx_t ind =
                  +ix * (sub->N) + iy * (sub->N) * (sub->N + 1) +
                  iz * (sub->N) * (sub->N + 1) * (sub->N + 1);
              xc += x[off + ind];
              yc += y[off + ind];
#if DIM == 3
              zc += z[off + ind];
#endif
            }
        xc /= 8;
        yc /= 8;
        zc /= 8;
        break;
      default:
        BFAM_ABORT("Unknown dimension: %" BFAM_LOCIDX_PRId,
                   (bfam_locidx_t)sub->dim);
      }

      for (bfam_locidx_t dof = 0; dof < sub->Np; dof++)
      {
        const bfam_locidx_t n = off + dof;
#if DIM == 2
        lua_global_function_call(L, 0, val_args->fname, "rrrrrrr>r", x[n], y[n],
                                 tmpz, t, xc, yc, zc, &field[n]);
#elif DIM == 3
        lua_global_function_call(L, 0, val_args->fname, "rrrrrrr>r", x[n], y[n],
                                 z[n], t, xc, yc, zc, &field[n]);
#else
#error "Bad Dimension"
#endif
      }
    }
    return;
  }
  else if (lua_isnumber(L, -1))
    val = (bfam_real_t)lua_tonumber(L, -1);
  else
    BFAM_ABORT("Did not find '%s' in lua as a function or number: "
               "lua message: '%s'",
               val_args->fname, lua_tostring(L, -1));
  lua_pop(L, 1);

  BFAM_ROOT_INFO("field_set_val: using '%s' with value %" BFAM_REAL_FMTe,
                 val_args->fname, val);

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = val;
}

static void field_set_val(bfam_locidx_t npoints, const char *name,
                          bfam_real_t t, bfam_real_t *restrict x,
                          bfam_real_t *restrict y, bfam_real_t *restrict z,
                          struct bfam_subdomain *s, void *arg,
                          bfam_real_t *restrict field)
{
  set_val_extended_args_t val_args;
  val_args.L = (lua_State *)arg;
  val_args.fname = name;
  field_set_val_extend(npoints, name, t, x, y, z, s, &val_args, field);
}

static void field_set_val_aux(bfam_locidx_t npoints, const char *name,
                              bfam_real_t time, bfam_real_t *restrict x,
                              bfam_real_t *restrict y, bfam_real_t *restrict z,
                              struct bfam_subdomain *s, void *arg,
                              bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(rho, "", "rho", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(lam, "", "lam", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(mu, "", "mu", &s->fields);

  if (strcmp(name, "rho_inv") == 0)
    for (bfam_locidx_t n = 0; n < npoints; ++n)
      field[n] = 1 / rho[n];
  else if (strcmp(name, "Zs") == 0)
    for (bfam_locidx_t n = 0; n < npoints; ++n)
      field[n] = BFAM_REAL_SQRT(rho[n] * mu[n]);
  else if (strcmp(name, "Zp") == 0)
    for (bfam_locidx_t n = 0; n < npoints; ++n)
      field[n] = BFAM_REAL_SQRT(rho[n] * (lam[n] + 2 * mu[n]));
  else
    BFAM_ABORT("Unknown auxilary field: '%s'", name);
}

typedef enum glue_info_type
{
  FRICTION,
  BOUNDARY,
  UNKNOWN,
} glue_info_type_t;

static void field_set_const(bfam_locidx_t npoints, const char *name,
                            bfam_real_t t, bfam_real_t *restrict x,
                            bfam_real_t *restrict y, bfam_real_t *restrict z,
                            struct bfam_subdomain *s, void *arg,
                            bfam_real_t *restrict field)
{
  bfam_real_t val = *((double *)arg);
  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = val;
}

static void field_set_friction_init_stress(
    bfam_locidx_t npoints, const char *name, bfam_real_t t,
    bfam_real_t *restrict x, bfam_real_t *restrict y, bfam_real_t *restrict z,
    struct bfam_subdomain *s, void *arg, bfam_real_t *restrict dummy_field)
{
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(nx1, "", "_grid_nx0", &s->glue_m->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(nx2, "", "_grid_nx1", &s->glue_m->fields);
#if DIM == 3
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(nx3, "", "_grid_nx2", &s->glue_m->fields);
#endif
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S11, "", "S11_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S12, "", "S12_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S13, "", "S13_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S22, "", "S22_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S23, "", "S23_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(S33, "", "S33_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tn, "", "Tn_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp1, "", "Tp1_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp2, "", "Tp2_0", &s->fields);
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(Tp3, "", "Tp3_0", &s->fields);

  for (bfam_locidx_t n = 0; n < npoints; ++n)
  {
    Tp1[n] = BEARD_D3_AP(S11[n] * nx1[n] + S12[n] * nx2[n], +S13[n] * nx3[n]);
    Tp2[n] = BEARD_D3_AP(S12[n] * nx1[n] + S22[n] * nx2[n], +S23[n] * nx3[n]);
    Tp3[n] = BEARD_D3_AP(S13[n] * nx1[n] + S23[n] * nx2[n], +S33[n] * nx3[n]);
    Tn[n] = BEARD_D3_AP(Tp1[n] * nx1[n] + Tp2[n] * nx2[n], +Tp3[n] * nx3[n]);
    Tp1[n] -= Tn[n] * nx1[n];
    Tp2[n] -= Tn[n] * nx2[n];
#if DIM == 3
    Tp3[n] -= Tn[n] * nx3[n];
#endif
  }
}

static void beard_grid_boundary(bfam_locidx_t npoints, const char *name,
                                bfam_real_t time, bfam_real_t *restrict x,
                                bfam_real_t *restrict y,
                                bfam_real_t *restrict z,
                                struct bfam_subdomain *s, void *arg,
                                bfam_real_t *restrict sJ)
{
  bfam_subdomain_dgx_t *sub_g = (bfam_subdomain_dgx_t *)s;
  bfam_subdomain_dgx_t *sub_m = (bfam_subdomain_dgx_t *)s->glue_m->sub_m;
  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)sub_g->base.glue_p;
  const bfam_locidx_t *restrict EToEm = glue_p->EToEm;
  const int8_t *restrict EToFm = glue_p->EToFm;

  /* get the fields we will need */
  bfam_dictionary_t *fields = &sub_g->base.glue_m->sub_m->fields;
  bfam_real_t *x_m = bfam_dictionary_get_value_ptr(fields, "_grid_x0");
  bfam_real_t *y_m = bfam_dictionary_get_value_ptr(fields, "_grid_x1");
#if DIM == 3
  bfam_real_t *z_m = bfam_dictionary_get_value_ptr(fields, "_grid_x2");
#endif

  BFAM_ASSERT(x_m != NULL);
  BFAM_ASSERT(y_m != NULL);
#if DIM == 3
  BFAM_ASSERT(z_m != NULL);
#endif

  const bfam_locidx_t K = sub_g->K;
  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    bfam_locidx_t *restrict fmask = sub_m->gmask[0][EToFm[k]];

    for (int n = 0; n < sub_g->Np; ++n)
    {
      x[n + k * sub_g->Np] = x_m[EToEm[k] * sub_m->Np + fmask[n]];
      y[n + k * sub_g->Np] = y_m[EToEm[k] * sub_m->Np + fmask[n]];
#if DIM == 3
      z[n + k * sub_g->Np] = z_m[EToEm[k] * sub_m->Np + fmask[n]];
#endif
    }
  }
}

static void beard_grid_glue(bfam_locidx_t npoints, const char *name,
                            bfam_real_t time, bfam_real_t *restrict x,
                            bfam_real_t *restrict y, bfam_real_t *restrict z,
                            struct bfam_subdomain *s, void *arg,
                            bfam_real_t *restrict tmp)
{
  bfam_subdomain_dgx_t *sub_g = (bfam_subdomain_dgx_t *)s;

  /* get the fields we will need */
  bfam_dictionary_t *fields_m = &sub_g->base.glue_m->fields;
  bfam_dictionary_t *fields_p = &sub_g->base.glue_p->fields;
  bfam_real_t *x_m = bfam_dictionary_get_value_ptr(fields_m, "_grid_x0");
  bfam_real_t *x_p = bfam_dictionary_get_value_ptr(fields_p, "_grid_x0");

  bfam_real_t *y_m = bfam_dictionary_get_value_ptr(fields_m, "_grid_x1");
  bfam_real_t *y_p = bfam_dictionary_get_value_ptr(fields_p, "_grid_x1");

#if DIM == 3
  bfam_real_t *z_m = bfam_dictionary_get_value_ptr(fields_m, "_grid_x2");
  bfam_real_t *z_p = bfam_dictionary_get_value_ptr(fields_p, "_grid_x2");
#endif

  BFAM_ASSERT(x_m != NULL);
  BFAM_ASSERT(x_p != NULL);
  BFAM_ASSERT(y_m != NULL);
  BFAM_ASSERT(y_p != NULL);

#if DIM == 3
  BFAM_ASSERT(z_m != NULL);
  BFAM_ASSERT(z_p != NULL);
#endif

  for (int n = 0; n < npoints; n++)
  {
    x[n] = BFAM_REAL(0.5) * (x_m[n] + x_p[n]);
    y[n] = BFAM_REAL(0.5) * (y_m[n] + y_p[n]);
#if DIM == 3
    z[n] = BFAM_REAL(0.5) * (z_m[n] + z_p[n]);
#endif
  }
}

static void domain_add_fields(beard_t *beard, prefs_t *prefs)
{

  const char *elastic[] = {"elastic*", "plastic*", "sponge*", NULL};
  bfam_domain_t *domain = (bfam_domain_t *)beard->domain;
  for (int f = 0; elastic_fields[f] != NULL; f++)
  {
    bfam_domain_add_field(domain, BFAM_DOMAIN_OR, elastic, elastic_fields[f]);
    bfam_domain_init_field(domain, BFAM_DOMAIN_OR, elastic, elastic_fields[f],
                           0, field_set_val, prefs->L);
  }
  const char *fields_aux[] = {"rho_inv", "Zs", "Zp", NULL};
  for (int f = 0; fields_aux[f] != NULL; f++)
  {
    bfam_domain_add_field(domain, BFAM_DOMAIN_OR, elastic, fields_aux[f]);
    bfam_domain_init_field(domain, BFAM_DOMAIN_OR, elastic, fields_aux[f], 0,
                           field_set_val_aux, NULL);
  }
  const char *sponge[] = {"sponge", NULL};
  for (int f = 0; sponge_fields[f] != NULL; f++)
  {
    bfam_domain_add_field(domain, BFAM_DOMAIN_OR, sponge, sponge_fields[f]);
    bfam_domain_init_field(domain, BFAM_DOMAIN_OR, sponge, sponge_fields[f], 0,
                           field_set_val, prefs->L);
  }

  /* Add plastic fields if we are supposed to */
  if (plastic_fields)
  {
    const char *plastic[] = {"plastic*", NULL};

    lua_State *L = prefs->L;
#ifdef BFAM_DEBUG
    int top = lua_gettop(L);
#endif
    lua_getglobal(L, "plastic");
    switch (prefs->plasticity)
    {
    case DUVAUT_LIONS:
      bfam_domain_add_tag((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                          plastic, DUVAUT_LIONS_TAG);
      break;
    default:
      BFAM_ABORT("Unknown plastic type for adding the tags");
    }

    for (int f = 0; plastic_fields[f] != NULL; f++)
    {
      bfam_domain_add_field(domain, BFAM_DOMAIN_OR, plastic, plastic_fields[f]);

      bfam_real_t value = 0;
      lua_pushstring(L, plastic_fields[f]);
      lua_gettable(L, -2);
      if (lua_isstring(L, -1) && !lua_isnumber(L, -1))
      {
        char fname[BFAM_BUFSIZ];
        strncpy(fname, lua_tostring(L, -1), BFAM_BUFSIZ);
        lua_pop(L, 1);

        set_val_extended_args_t args;
        args.fname = fname;
        args.L = L;
        bfam_domain_init_field(domain, BFAM_DOMAIN_OR, plastic,
                               plastic_fields[f], 0, field_set_val_extend,
                               &args);
      }
      else
      {
        if (lua_isnumber(L, -1))
        {
          value = (bfam_real_t)lua_tonumber(L, -1);
          BFAM_ROOT_INFO("plasticity: using '%s' with value %" BFAM_REAL_PRIe,
                         plastic_fields[f], value);
        }
        else
          BFAM_ROOT_WARNING(" plastic does not contain `%s',"
                            " using default %" BFAM_REAL_PRIe,
                            plastic_fields[f], value);
        bfam_domain_init_field(domain, BFAM_DOMAIN_OR, plastic,
                               plastic_fields[f], 0, field_set_const, &value);
        lua_pop(L, 1);
      }
    }
    lua_pop(L, -1);
#ifdef BFAM_DEBUG
    BFAM_ASSERT(top == lua_gettop(L));
#endif
  }

  /* add glue fields */
  const char *glue[] = {"_glue_parallel", "_glue_local", NULL};
  for (int f = 0; comm_args_scalars[f] != NULL; f++)
  {
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue,
                                comm_args_scalars[f]);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue,
                               comm_args_scalars[f]);
  }
  for (int f = 0; comm_args_vectors[f] != NULL; f++)
  {
    char name[BFAM_BUFSIZ];
    snprintf(name, BFAM_BUFSIZ, "%sn", comm_args_vectors[f]);
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, name);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, name);
    snprintf(name, BFAM_BUFSIZ, "%sp1", comm_args_vectors[f]);
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, name);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, name);
    snprintf(name, BFAM_BUFSIZ, "%sp2", comm_args_vectors[f]);
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, name);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, name);
    snprintf(name, BFAM_BUFSIZ, "%sp3", comm_args_vectors[f]);
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, name);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, name);
  }
  for (int f = 0; comm_args_tensors[f] != NULL; f++)
  {
    char name[BFAM_BUFSIZ];
    snprintf(name, BFAM_BUFSIZ, "%sn", comm_args_tensors[f]);
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, name);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, name);
    snprintf(name, BFAM_BUFSIZ, "%sp1", comm_args_tensors[f]);
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, name);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, name);
    snprintf(name, BFAM_BUFSIZ, "%sp2", comm_args_tensors[f]);
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, name);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, name);
    snprintf(name, BFAM_BUFSIZ, "%sp3", comm_args_tensors[f]);
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, name);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, name);
  }

  /* exchange material properties to glue */
  const char *glue_mat[] =
  {
    "Zs",
    "Zp",
    "_grid_x0",
    "_grid_x1",
#if DIM == 3
    "_grid_x2",
#endif
    NULL
  };
  for (bfam_locidx_t g = 0; glue_mat[g] != NULL; g++)
  {
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue, glue_mat[g]);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue, glue_mat[g]);
  }
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x0");
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x1");
#if DIM == 3
  bfam_domain_add_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x2");
#endif

  const char *glue_face_scalar[] =
  {
    "_grid_nx0",
    "_grid_nx1",
#if DIM == 3
    "_grid_nx2",
#endif
    NULL
  };
  for (bfam_locidx_t g = 0; glue_face_scalar[g] != NULL; g++)
  {
    bfam_domain_add_minus_field(domain, BFAM_DOMAIN_OR, glue,
                                glue_face_scalar[g]);
    bfam_domain_add_plus_field(domain, BFAM_DOMAIN_OR, glue,
                               glue_face_scalar[g]);
  }

  bfam_communicator_t material_comm;

  bfam_subdomain_comm_args_t mat_args;
  const char *mat_NULL[] = {NULL};
  mat_args.scalars_m = glue_mat;
  mat_args.vectors_m = mat_NULL;
  mat_args.vector_components_m = mat_NULL;
  mat_args.tensors_m = mat_NULL;
  mat_args.tensor_components_m = mat_NULL;
  mat_args.face_scalars_m = glue_face_scalar;

  mat_args.scalars_p = glue_mat;
  mat_args.vectors_p = mat_NULL;
  mat_args.vector_components_p = mat_NULL;
  mat_args.tensors_p = mat_NULL;
  mat_args.tensor_components_p = mat_NULL;
  mat_args.face_scalars_p = glue_face_scalar;

  mat_args.user_comm_info = NULL;
  mat_args.user_put_send_buffer = NULL;
  mat_args.user_get_recv_buffer = NULL;
  mat_args.user_data = NULL;

  mat_args.user_prefix_function = NULL;

  bfam_communicator_init(&material_comm, domain, BFAM_DOMAIN_OR, glue,
                         beard->mpicomm, 10, &mat_args);
  bfam_communicator_start(&material_comm);
  bfam_communicator_finish(&material_comm);
  bfam_communicator_free(&material_comm);

  /*
   * we can trick init fields into handling locations
   * to glue
   */
  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, glue, "_grid_x0", 0,
                         beard_grid_glue, NULL);

  const char *boundary[] = {"_glue_boundary", NULL};
  const char *boundary_fields[] =
  {
    "_grid_x0",
    "_grid_x1",
#if DIM == 3
    "_grid_x2",
#endif
    NULL
  };
  for (bfam_locidx_t g = 0; boundary_fields[g] != NULL; g++)
    bfam_domain_add_field(domain, BFAM_DOMAIN_OR, boundary, boundary_fields[g]);

  bfam_domain_init_field(domain, BFAM_DOMAIN_OR, boundary, "_grid_x1", 0,
                         beard_grid_boundary, NULL);

  /*
   * Loop through glue info
   */

  lua_State *L = prefs->L;
#ifdef BFAM_DEBUG
  int top = lua_gettop(L);
#endif

  lua_getglobal(L, "glue_info");
  if (!lua_istable(L, -1))
    BFAM_ROOT_WARNING("table `%s' not found", "glue_info");
  else
  {
    luaL_checktype(L, -1, LUA_TTABLE);

    int N_glueids = lua_objlen(L, -1);
    BFAM_LDEBUG("N_glueids: %3d", N_glueids);

    BFAM_ROOT_VERBOSE("Reading glue info from lua");
    for (int i = 1; i <= N_glueids; ++i)
    {
      BFAM_ROOT_VERBOSE("glue ID: %d", i);
      /*
       * get glue_id tag
       */
      const char *this_glue[2];
      char this_glue_tag[BFAM_BUFSIZ];
      this_glue[0] = this_glue_tag;
      this_glue[1] = NULL;
      snprintf(this_glue_tag, BFAM_BUFSIZ, "_glue_id_%jd", (intmax_t)i);

      lua_rawgeti(L, -1, i);
      luaL_checktype(L, -1, LUA_TTABLE);

      lua_pushstring(L, "type");
      lua_gettable(L, -2);
      luaL_checktype(L, -1, LUA_TSTRING);

      size_t match_len = BFAM_MIN(lua_strlen(L, -1), 8);

      glue_info_type_t type = UNKNOWN;
      if (0 == strncmp(lua_tostring(L, -1), "friction", match_len))
        type = FRICTION;
      else if (0 == strncmp(lua_tostring(L, -1), "boundary", match_len))
        type = BOUNDARY;

      BFAM_ROOT_VERBOSE("  type: %d", type);
      lua_pop(L, 1);

      lua_pushstring(L, "tag");
      lua_gettable(L, -2);
      if (lua_isstring(L, -1))
      {
        const char *tag = lua_tostring(L, -1);

        bfam_domain_add_tag((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                            this_glue, tag);
        if (type == FRICTION)
        {
          bfam_domain_add_tag((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                              this_glue, "friction");
          const char **tmp_fields = friction_fields;
          size_t match_len = lua_strlen(L, -1);
          if (0 == strncmp(lua_tostring(L, -1), "slip weakening", match_len))
          {
            friction_fields = slip_weakening_fields;
            friction_rates = slip_weakening_rates;
          }
          else if (0 == strncmp(lua_tostring(L, -1), "ageing law", match_len))
          {
            friction_fields = rate_and_state_fields;
            friction_rates = rate_and_state_rates;
          }
          else
            BFAM_ABORT("Uknown friction tag type: %s"
                       " (can be: 'slip weakening', 'ageing law'",
                       tag);
          BFAM_ABORT_IF_NOT(tmp_fields == NULL || friction_fields == tmp_fields,
                            "Not configured for multiple friction laws");
        }
      }
      else
        BFAM_ABORT("No tag for glue_id %d", i);
      lua_pop(L, 1);

      if (type == FRICTION)
      {
        for (int f = 0; friction_fields[f] != NULL; ++f)
        {
          bfam_domain_add_field(domain, BFAM_DOMAIN_OR, this_glue,
                                friction_fields[f]);

          bfam_real_t value = 0;
          if (strncmp(friction_fields[f], "Trup", BFAM_BUFSIZ) == 0)
            value = bfam_real_nan("");
          lua_pushstring(L, friction_fields[f]);
          lua_gettable(L, -2);
          if (lua_isstring(L, -1) && !lua_isnumber(L, -1))
          {
            char fname[BFAM_BUFSIZ];
            strncpy(fname, lua_tostring(L, -1), BFAM_BUFSIZ);
            lua_pop(L, 1);

            set_val_extended_args_t args;
            args.fname = fname;
            args.L = L;
            bfam_domain_init_field(domain, BFAM_DOMAIN_OR, this_glue,
                                   friction_fields[f], 0, field_set_val_extend,
                                   &args);
          }
          else
          {
            if (lua_isnumber(L, -1))
            {
              value = (bfam_real_t)lua_tonumber(L, -1);
            }
            else
            {
              BFAM_ROOT_WARNING(" glue %d does not contain `%s',"
                                " using default %" BFAM_REAL_PRIe,
                                i, friction_fields[f], value);
            }
            bfam_domain_init_field(domain, BFAM_DOMAIN_OR, this_glue,
                                   friction_fields[f], 0, field_set_const,
                                   &value);
            lua_pop(L, 1);
          }
        }

        bfam_domain_init_field(domain, BFAM_DOMAIN_OR, this_glue, "Tp1_0", 0,
                               field_set_friction_init_stress, NULL);
      }

      lua_pop(L, 1);
    }
  }

  /*
   * Set default boundary condition
   */
  const char *default_boundary_glue[] = {"_glue_boundary", "_glue_id_-1", NULL};
  bfam_domain_add_tag((bfam_domain_t *)beard->domain, BFAM_DOMAIN_AND,
                      default_boundary_glue, prefs->default_boundary_tag);

  lua_pop(L, -1);
  BFAM_ASSERT(top == lua_gettop(L));
}

static void init_domain(beard_t *beard, prefs_t *prefs)
{
  /* Set up the connectivity */
  brick_args_t *brick_args = prefs->brick_args;
  if (brick_args != NULL)
  {
    beard->conn = p4est_connectivity_new_brick(brick_args->nx, brick_args->ny,
#if DIM == 3
                                               brick_args->nz,
#endif
                                               brick_args->periodic_x,
                                               brick_args->periodic_y
#if DIM == 3
                                               ,
                                               brick_args->periodic_z
#endif
                                               );

    brick_args->brick_to_tree = NULL;
    if (!brick_args->p4est_brick)
    {
      const bfam_locidx_t num_trees = beard->conn->num_trees;
      brick_args->brick_to_tree =
          bfam_malloc(sizeof(bfam_locidx_t) * num_trees);
      const double *verts = beard->conn->vertices;
      const p4est_topidx_t *tree_to_verts = beard->conn->tree_to_vertex;
      for (bfam_locidx_t treeid = 0; treeid < num_trees; treeid++)
      {
#if DIM == 2
        bfam_locidx_t brick =
            (bfam_locidx_t)verts[3 * tree_to_verts[4 * treeid] + 0] +
            brick_args->nx *
                ((bfam_locidx_t)verts[3 * tree_to_verts[4 * treeid] + 1]);
#elif DIM == 3
        bfam_locidx_t brick =
            (bfam_locidx_t)verts[3 * tree_to_verts[8 * treeid] + 0] +
            brick_args->nx *
                ((bfam_locidx_t)verts[3 * tree_to_verts[8 * treeid] + 1] +
                 brick_args->ny *
                     ((bfam_locidx_t)verts[3 * tree_to_verts[8 * treeid] + 2]));
#endif
        brick_args->brick_to_tree[brick] = treeid;
        BFAM_LDEBUG("brick %04" BFAM_LOCIDX_PRId
                    " <--> tree %04" BFAM_LOCIDX_PRId,
                    brick, brick_args->brick_to_tree[brick]);
        BFAM_ASSERT(0 <= brick && brick < num_trees);
      }
    }
  }
  else if (strncmp(prefs->mesh_file, "", BFAM_BUFSIZ))
  {
    beard->conn = p4est_connectivity_read_inp(prefs->mesh_file);
  }
  else
    BFAM_ABORT("no connectivity");

  /* let the user modify the connectivity vertices */
  for (int i = 0; i < beard->conn->num_vertices; i++)
  {
    bfam_real_t x = beard->conn->vertices[i * 3 + 0];
    bfam_real_t y = beard->conn->vertices[i * 3 + 1];
    bfam_real_t z = beard->conn->vertices[i * 3 + 2];
    int result = lua_global_function_call(prefs->L, 0, "connectivity_vertices",
                                          "rrr>rrr", x, y, z, &x, &y, &z);
    if (result != 0)
      break;
    beard->conn->vertices[i * 3 + 0] = x;
    beard->conn->vertices[i * 3 + 1] = y;
    beard->conn->vertices[i * 3 + 2] = z;
  }
  /* create the domain */
  beard->domain = bfam_domain_pxest_new(beard->mpicomm, beard->conn);

  /* call user refinement function */
  lua_getglobal(prefs->L, "refinement_function");
  if (lua_isfunction(prefs->L, -1))
  {
    lua_pop(prefs->L, 1);
    void *current_user_pointer = beard->domain->pxest->user_pointer;
    beard->domain->pxest->user_pointer = prefs->L;
    p4est_refine(beard->domain->pxest, 1, refine_fn,
                 bfam_domain_pxest_init_callback);

    beard->domain->pxest->user_pointer = current_user_pointer;
  }
  else
    BFAM_ROOT_WARNING("function `%s' not found in lua file",
                      "refinement_function");

  p4est_balance(beard->domain->pxest, P4EST_CONNECT_CORNER, NULL);

  /*
   * This is to statically refine all cells of a balanced mesh, since it's
   * already balanced it will remain balanced
   */
  int stat_ref = lua_get_global_int(prefs->L, "static_refinement", 0);
  for (int i = 0; i < stat_ref; i++)
    p4est_refine(beard->domain->pxest, 0, static_refine_fn,
                 bfam_domain_pxest_init_callback);

  p4est_partition(beard->domain->pxest, 1, NULL);

  /* split the domain */
  split_domain(beard, prefs);

  /* add fields to the domnain */
  domain_add_fields(beard, prefs);
}

static void field_zero(bfam_locidx_t npoints, const char *name,
                       bfam_real_t time, bfam_real_t *restrict x,
                       bfam_real_t *restrict y, bfam_real_t *restrict z,
                       struct bfam_subdomain *s, void *arg,
                       bfam_real_t *restrict field)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(field, 32);

  for (bfam_locidx_t n = 0; n < npoints; ++n)
    field[n] = 0;
}

void aux_rates(bfam_subdomain_t *thisSubdomain, const char *prefix)
{
  if (bfam_subdomain_has_tag(thisSubdomain, "elastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "plastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "sponge*"))
  {
    const char *fields[] = {"v1",  "v2",  "v3",  "S11", "S22",
                            "S33", "S12", "S13", "S23", NULL};

    char field[BFAM_BUFSIZ];
    for (int f = 0; fields[f] != NULL; ++f)
    {
      snprintf(field, BFAM_BUFSIZ, "%s%s", prefix, fields[f]);
      thisSubdomain->field_add(thisSubdomain, field);
      bfam_subdomain_field_init(thisSubdomain, field, 0, field_zero, NULL);
    }
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "friction"))
  {
    char field[BFAM_BUFSIZ];
    for (int f = 0; friction_rates[f] != NULL; ++f)
    {
      snprintf(field, BFAM_BUFSIZ, "%s%s", prefix, friction_rates[f]);
      thisSubdomain->field_add(thisSubdomain, field);
      bfam_subdomain_field_init(thisSubdomain, field, 0, field_zero, NULL);
    }
  }
}

void glue_rates(bfam_subdomain_t *thisSubdomain, const char *prefix)
{
  if (bfam_subdomain_has_tag(thisSubdomain, "_glue_parallel") ||
      bfam_subdomain_has_tag(thisSubdomain, "_glue_local"))
  {
    BFAM_LDEBUG("Adding glue rate to %s", thisSubdomain->name);
    for (int f = 0; comm_args_scalars[f] != NULL; f++)
    {
      char field[BFAM_BUFSIZ];
      snprintf(field, BFAM_BUFSIZ, "%s%s", prefix, comm_args_scalars[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);
    }
    for (int f = 0; comm_args_vectors[f] != NULL; f++)
    {
      char field[BFAM_BUFSIZ];
      snprintf(field, BFAM_BUFSIZ, "%s%sn", prefix, comm_args_vectors[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);

      snprintf(field, BFAM_BUFSIZ, "%s%sp1", prefix, comm_args_vectors[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);

      snprintf(field, BFAM_BUFSIZ, "%s%sp2", prefix, comm_args_vectors[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);

      snprintf(field, BFAM_BUFSIZ, "%s%sp3", prefix, comm_args_vectors[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);
    }
    for (int f = 0; comm_args_tensors[f] != NULL; f++)
    {
      char field[BFAM_BUFSIZ];
      snprintf(field, BFAM_BUFSIZ, "%s%sn", prefix, comm_args_tensors[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);

      snprintf(field, BFAM_BUFSIZ, "%s%sp1", prefix, comm_args_tensors[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);

      snprintf(field, BFAM_BUFSIZ, "%s%sp2", prefix, comm_args_tensors[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);

      snprintf(field, BFAM_BUFSIZ, "%s%sp3", prefix, comm_args_tensors[f]);
      thisSubdomain->field_minus_add(thisSubdomain, field);
      thisSubdomain->field_plus_add(thisSubdomain, field);
    }
  }
}

void scale_rates_elastic(bfam_subdomain_dgx_t *sub, const char *rate_prefix,
                         const bfam_long_real_t a)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_scale_rates_elastic_2_##order(sub->N, sub, rate_prefix, a);      \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_scale_rates_elastic_3_##order(sub->N, sub, rate_prefix, a);      \
    break;
#else
#error "Bad Dimension"
#endif

  switch (sub->N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_scale_rates_elastic_2_(sub->N, sub, rate_prefix, a);
#elif DIM == 3
    beard_dgx_scale_rates_elastic_3_(sub->N, sub, rate_prefix, a);
#else
#error "Bad Dimension"
#endif

    break;
  }
#undef X
}

void scale_rates_friction(bfam_subdomain_dgx_t *sub, const char *rate_prefix,
                          const bfam_long_real_t a)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_scale_rates_interface_2_##order(sub->N, sub, rate_prefix, a,     \
                                              friction_rates);                 \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_scale_rates_interface_3_##order(sub->N, sub, rate_prefix, a,     \
                                              friction_rates);                 \
    break;
#else
#error "Bad Dimension"
#endif

  switch (sub->N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_scale_rates_interface_2_(sub->N, sub, rate_prefix, a,
                                       friction_rates);
#elif DIM == 3
    beard_dgx_scale_rates_interface_3_(sub->N, sub, rate_prefix, a,
                                       friction_rates);
#else
#error "Bad Dimension"
#endif

    break;
  }
#undef X
}

void scale_rates(bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
                 const bfam_long_real_t a)
{
  BFAM_ASSERT(rate_prefix);
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain, "_subdomain_dgx"));
  if (bfam_subdomain_has_tag(thisSubdomain, "elastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "plastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "sponge*"))
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;
    scale_rates_elastic(sub, rate_prefix, a);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "friction"))
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;
    scale_rates_friction(sub, rate_prefix, a);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "_glue_boundary"))
    ;
  else if (bfam_subdomain_has_tag(thisSubdomain, "_glue_parallel"))
    ;
  else if (bfam_subdomain_has_tag(thisSubdomain, "_glue_local"))
    ;
  else
    BFAM_ABORT("Unknown subdomain: %s", thisSubdomain->name);
}

static void rupture_time(int N, bfam_subdomain_dgx_t *sub,
                         const char *field_prefix, const bfam_long_real_t t,
                         const bfam_real_t Vrup)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_rupture_time_2_##order(N, sub, field_prefix, t, Vrup);           \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_rupture_time_3_##order(N, sub, field_prefix, t, Vrup);           \
    break;
#else
#error "Bad Dimension"
#endif

  switch (N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_rupture_time_2_(N, sub, field_prefix, t, Vrup);
#elif DIM == 3
    beard_dgx_rupture_time_3_(N, sub, field_prefix, t, Vrup);
#else
#error "Bad Dimension"
#endif
    break;
  }
#undef X
}

static void duvaut_lions_return_map(int N, bfam_subdomain_dgx_t *sub,
                                    const char *field_prefix,
                                    const bfam_long_real_t t,
                                    const bfam_long_real_t dt)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_duvaut_lions_return_map_2_##order(N, sub, field_prefix, t, dt);  \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_duvaut_lions_return_map_3_##order(N, sub, field_prefix, t, dt);  \
    break;
#else
#error "Bad Dimension"
#endif

  switch (N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_duvaut_lions_return_map_2_(N, sub, field_prefix, t, dt);
#elif DIM == 3
    beard_dgx_duvaut_lions_return_map_3_(N, sub, field_prefix, t, dt);
#else
#error "Bad Dimension"
#endif
    break;
  }
#undef X
}

static void intra_rhs_elastic(int N, bfam_subdomain_dgx_t *sub,
                              const char *rate_prefix, const char *field_prefix,
                              const bfam_long_real_t t)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_intra_rhs_elastic_2_##order(N, sub, rate_prefix, field_prefix,   \
                                          t);                                  \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_intra_rhs_elastic_3_##order(N, sub, rate_prefix, field_prefix,   \
                                          t);                                  \
    break;
#else
#error "Bad Dimension"
#endif

  switch (N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_intra_rhs_elastic_2_(N, sub, rate_prefix, field_prefix, t);
#elif DIM == 3
    beard_dgx_intra_rhs_elastic_3_(N, sub, rate_prefix, field_prefix, t);
#else
#error "Bad Dimension"
#endif
    break;
  }
#undef X
}

static void intra_rhs_sponge(int N, bfam_subdomain_dgx_t *sub,
                             const char *rate_prefix, const char *field_prefix,
                             const bfam_long_real_t t)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_intra_rhs_sponge_2_##order(N, sub, rate_prefix, field_prefix,    \
                                         t);                                   \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_intra_rhs_sponge_3_##order(N, sub, rate_prefix, field_prefix,    \
                                         t);                                   \
    break;
#else
#error "Bad Dimension"
#endif

  switch (N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_intra_rhs_sponge_2_(N, sub, rate_prefix, field_prefix, t);
#elif DIM == 3
    beard_dgx_intra_rhs_sponge_3_(N, sub, rate_prefix, field_prefix, t);
#else
#error "Bad Dimension"
#endif
    break;
  }
#undef X
}

void intra_rhs(bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
               const char *minus_rate_prefix, const char *field_prefix,
               const bfam_long_real_t t)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain, "_subdomain_dgx"));

  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;
  if (bfam_subdomain_has_tag(thisSubdomain, "elastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "plastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "sponge*"))
  {
    intra_rhs_elastic(sub->N, sub, rate_prefix, field_prefix, t);
    if (bfam_subdomain_has_tag(thisSubdomain, "sponge"))
      intra_rhs_sponge(sub->N, sub, rate_prefix, field_prefix, t);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "_glue_boundary") ||
           bfam_subdomain_has_tag(thisSubdomain, "_glue_parallel") ||
           bfam_subdomain_has_tag(thisSubdomain, "_glue_local"))
    ;
  else
    BFAM_ABORT("Unknown subdomain: %s", thisSubdomain->name);
}

void inter_rhs_boundary(int N, bfam_subdomain_dgx_t *sub,
                        const char *rate_prefix, const char *field_prefix,
                        const bfam_long_real_t t, const bfam_real_t R)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_inter_rhs_boundary_2_##order(N, sub, rate_prefix, field_prefix,  \
                                           t, R);                              \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_inter_rhs_boundary_3_##order(N, sub, rate_prefix, field_prefix,  \
                                           t, R);                              \
    break;
#else
#error "bad dimension"
#endif

  switch (N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_inter_rhs_boundary_2_(N, sub, rate_prefix, field_prefix, t, R);
#elif DIM == 3
    beard_dgx_inter_rhs_boundary_3_(N, sub, rate_prefix, field_prefix, t, R);
#else
#error "bad dimension"
#endif
    break;
  }
#undef X
}

void inter_rhs_slip_weakening_interface(int N, bfam_subdomain_dgx_t *sub,
                                        const char *rate_prefix,
                                        const char *minus_rate_prefix,
                                        const char *field_prefix,
                                        const bfam_long_real_t t)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_inter_rhs_slip_weakening_interface_2_##order(                    \
        N, sub, rate_prefix, minus_rate_prefix, field_prefix, t);              \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_inter_rhs_slip_weakening_interface_3_##order(                    \
        N, sub, rate_prefix, minus_rate_prefix, field_prefix, t);              \
    break;
#else
#error "bad dimension"
#endif

  switch (N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_inter_rhs_slip_weakening_interface_2_(
        N, sub, rate_prefix, minus_rate_prefix, field_prefix, t);
#elif DIM == 3
    beard_dgx_inter_rhs_slip_weakening_interface_3_(
        N, sub, rate_prefix, minus_rate_prefix, field_prefix, t);
#else
#error "bad dimension"
#endif
    break;
  }
#undef X
}

void inter_rhs_ageing_law_interface(int N, bfam_subdomain_dgx_t *sub,
                                    const char *rate_prefix,
                                    const char *minus_rate_prefix,
                                    const char *field_prefix,
                                    const bfam_long_real_t t)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_inter_rhs_ageing_law_interface_2_##order(                        \
        N, sub, rate_prefix, minus_rate_prefix, field_prefix, t);              \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_inter_rhs_ageing_law_interface_3_##order(                        \
        N, sub, rate_prefix, minus_rate_prefix, field_prefix, t);              \
    break;
#else
#error "bad dimension"
#endif

  switch (N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_inter_rhs_ageing_law_interface_2_(
        N, sub, rate_prefix, minus_rate_prefix, field_prefix, t);
#elif DIM == 3
    beard_dgx_inter_rhs_ageing_law_interface_3_(
        N, sub, rate_prefix, minus_rate_prefix, field_prefix, t);
#else
#error "bad dimension"
#endif
    break;
  }
#undef X
}

void inter_rhs_interface(int N, bfam_subdomain_dgx_t *sub,
                         const char *rate_prefix, const char *field_prefix,
                         const bfam_long_real_t t)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_inter_rhs_interface_2_##order(N, sub, rate_prefix, field_prefix, \
                                            t);                                \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_inter_rhs_interface_3_##order(N, sub, rate_prefix, field_prefix, \
                                            t);                                \
    break;
#else
#error "bad dimension"
#endif

  switch (N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_inter_rhs_interface_2_(N, sub, rate_prefix, field_prefix, t);
#elif DIM == 3
    beard_dgx_inter_rhs_interface_3_(N, sub, rate_prefix, field_prefix, t);
#else
#error "bad dimension"
#endif
    break;
  }
#undef X
}

void inter_rhs(bfam_subdomain_t *thisSubdomain, const char *rate_prefix,
               const char *minus_rate_prefix, const char *field_prefix,
               const bfam_long_real_t t)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain, "_subdomain_dgx"));

  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;
  if (bfam_subdomain_has_tag(thisSubdomain, "elastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "plastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "sponge*"))
    ;
  else if (bfam_subdomain_has_tag(thisSubdomain, "slip weakening"))
  {
    BFAM_ASSERT(rate_prefix);
    inter_rhs_slip_weakening_interface(
        ((bfam_subdomain_dgx_t *)sub->base.glue_m->sub_m)->N, sub, rate_prefix,
        minus_rate_prefix, field_prefix, t);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "ageing law"))
  {
    BFAM_ASSERT(rate_prefix);
    inter_rhs_ageing_law_interface(
        ((bfam_subdomain_dgx_t *)sub->base.glue_m->sub_m)->N, sub, rate_prefix,
        minus_rate_prefix, field_prefix, t);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "non-reflecting"))
  {
    BFAM_ASSERT(minus_rate_prefix);
    inter_rhs_boundary(((bfam_subdomain_dgx_t *)sub->base.glue_m->sub_m)->N,
                       sub, minus_rate_prefix, field_prefix, t, 0);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "free surface"))
  {
    BFAM_ASSERT(minus_rate_prefix);
    inter_rhs_boundary(((bfam_subdomain_dgx_t *)sub->base.glue_m->sub_m)->N,
                       sub, minus_rate_prefix, field_prefix, t, 1);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "rigid"))
  {
    BFAM_ASSERT(minus_rate_prefix);
    inter_rhs_boundary(((bfam_subdomain_dgx_t *)sub->base.glue_m->sub_m)->N,
                       sub, minus_rate_prefix, field_prefix, t, -1);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "_glue_parallel") ||
           bfam_subdomain_has_tag(thisSubdomain, "_glue_local"))
  {
    BFAM_ASSERT(minus_rate_prefix);
    inter_rhs_interface(((bfam_subdomain_dgx_t *)sub->base.glue_m->sub_m)->N,
                        sub, minus_rate_prefix, field_prefix, t);
  }
  else
    BFAM_ABORT("Unknown subdomain: %s", thisSubdomain->name);
}

void add_rates_friction(bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
                        const char *field_prefix_rhs, const char *rate_prefix,
                        const bfam_long_real_t a)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_add_rates_interface_2_##order(sub->N, sub, field_prefix_lhs,     \
                                            field_prefix_rhs, rate_prefix, a,  \
                                            friction_rates);                   \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_add_rates_interface_3_##order(sub->N, sub, field_prefix_lhs,     \
                                            field_prefix_rhs, rate_prefix, a,  \
                                            friction_rates);                   \
    break;
#else
#error "bad dimension"
#endif

  switch (sub->N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_add_rates_interface_2_(sub->N, sub, field_prefix_lhs,
                                     field_prefix_rhs, rate_prefix, a,
                                     friction_rates);
#elif DIM == 3
    beard_dgx_add_rates_interface_3_(sub->N, sub, field_prefix_lhs,
                                     field_prefix_rhs, rate_prefix, a,
                                     friction_rates);
#else
#error "bad dimension"
#endif
    break;
  }
#undef X
}

void add_rates_elastic(bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
                       const char *field_prefix_rhs, const char *rate_prefix,
                       const bfam_long_real_t a)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_add_rates_elastic_2_##order(sub->N, sub, field_prefix_lhs,       \
                                          field_prefix_rhs, rate_prefix, a);   \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_add_rates_elastic_3_##order(sub->N, sub, field_prefix_lhs,       \
                                          field_prefix_rhs, rate_prefix, a);   \
    break;
#else
#error "bad dimension"
#endif

  switch (sub->N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_add_rates_elastic_2_(sub->N, sub, field_prefix_lhs,
                                   field_prefix_rhs, rate_prefix, a);
#elif DIM == 3
    beard_dgx_add_rates_elastic_3_(sub->N, sub, field_prefix_lhs,
                                   field_prefix_rhs, rate_prefix, a);
#else
#error "bad dimension"
#endif
    break;
  }
#undef X
}

void add_rates_elastic_glue_p(bfam_subdomain_dgx_t *sub,
                              const char *field_prefix_lhs,
                              const char *field_prefix_rhs,
                              const char *rate_prefix, const bfam_long_real_t a)
{
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_add_rates_glue_p_2_##order(                                      \
        sub->N, sub, field_prefix_lhs, field_prefix_rhs, rate_prefix, a,       \
        comm_args_scalars, comm_args_vectors, comm_args_tensors);              \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_add_rates_glue_p_3_##order(                                      \
        sub->N, sub, field_prefix_lhs, field_prefix_rhs, rate_prefix, a,       \
        comm_args_scalars, comm_args_vectors, comm_args_tensors);              \
    break;
#else
#error "bad dimension"
#endif

  switch (sub->N)
  {
    BFAM_LIST_OF_DGX_NORDERS
  default:
#if DIM == 2
    beard_dgx_add_rates_glue_p_2_(
        sub->N, sub, field_prefix_lhs, field_prefix_rhs, rate_prefix, a,
        comm_args_scalars, comm_args_vectors, comm_args_tensors);
#elif DIM == 3
    beard_dgx_add_rates_glue_p_3_(
        sub->N, sub, field_prefix_lhs, field_prefix_rhs, rate_prefix, a,
        comm_args_scalars, comm_args_vectors, comm_args_tensors);
#else
#error "bad dimension"
#endif
    break;
  }
#undef X
}

void add_rates_glue_p(bfam_subdomain_t *thisSubdomain,
                      const char *field_prefix_lhs,
                      const char *field_prefix_rhs, const char *rate_prefix,
                      const bfam_long_real_t a)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain, "_subdomain_dgx"));
  if (bfam_subdomain_has_tag(thisSubdomain, "_glue_parallel") ||
      bfam_subdomain_has_tag(thisSubdomain, "_glue_local"))
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;
    add_rates_elastic_glue_p(sub, field_prefix_lhs, field_prefix_rhs,
                             rate_prefix, a);
  }
  else
    BFAM_ABORT("Unknown subdomain: %s", thisSubdomain->name);
}

void add_rates(bfam_subdomain_t *thisSubdomain, const char *field_prefix_lhs,
               const char *field_prefix_rhs, const char *rate_prefix,
               const bfam_long_real_t a)
{
  BFAM_ASSERT(bfam_subdomain_has_tag(thisSubdomain, "_subdomain_dgx"));
  if (bfam_subdomain_has_tag(thisSubdomain, "elastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "plastic*") ||
      bfam_subdomain_has_tag(thisSubdomain, "sponge*"))
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;
    add_rates_elastic(sub, field_prefix_lhs, field_prefix_rhs, rate_prefix, a);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "friction"))
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;
    add_rates_friction(sub, field_prefix_lhs, field_prefix_rhs, rate_prefix, a);
  }
  else if (bfam_subdomain_has_tag(thisSubdomain, "_glue_boundary") ||
           bfam_subdomain_has_tag(thisSubdomain, "_glue_parallel") ||
           bfam_subdomain_has_tag(thisSubdomain, "_glue_local"))
    ;
  else
    BFAM_ABORT("Unknown subdomain: %s", thisSubdomain->name);
}

static void init_time_stepper(beard_t *beard, prefs_t *prefs,
                              bfam_locidx_t num_lvl)
{
  beard->comm_args = bfam_malloc(sizeof(bfam_subdomain_comm_args_t));
  bfam_subdomain_comm_args_t *args = beard->comm_args;

  args->scalars_m = comm_args_scalars;
  args->vectors_m = comm_args_vectors;
  args->vector_components_m = comm_args_vector_components;
  args->tensors_m = comm_args_tensors;
  args->tensor_components_m = comm_args_tensor_components;
  args->face_scalars_m = comm_args_face_scalars;

  args->scalars_p = comm_args_scalars;
  args->vectors_p = comm_args_vectors;
  args->vector_components_p = comm_args_vector_components;
  args->tensors_p = comm_args_tensors;
  args->tensor_components_p = comm_args_tensor_components;
  args->face_scalars_p = comm_args_face_scalars;

  args->user_comm_info = NULL;
  args->user_put_send_buffer = NULL;
  args->user_get_recv_buffer = NULL;
  args->user_data = NULL;
  args->user_prefix_function = NULL;

  const char *timestep_tags[] = {"_volume", "_glue_parallel", "_glue_local",
                                 "_glue_boundary", NULL};
  const char *glue[] = {"_glue_parallel", "_glue_local", NULL};

  if (prefs->lsrk_method != BFAM_TS_LSRK_NOOP)
    beard->beard_ts = (bfam_ts_t *)bfam_ts_lsrk_new(
        (bfam_domain_t *)beard->domain, prefs->lsrk_method, BFAM_DOMAIN_OR,
        timestep_tags, BFAM_DOMAIN_OR, glue, beard->mpicomm, 10,
        beard->comm_args, &aux_rates, &scale_rates, &intra_rhs, &inter_rhs,
        &add_rates);
  else if (prefs->adams_method != BFAM_TS_ADAMS_NOOP)
    beard->beard_ts = (bfam_ts_t *)bfam_ts_adams_new(
        (bfam_domain_t *)beard->domain, prefs->adams_method, BFAM_DOMAIN_OR,
        timestep_tags, BFAM_DOMAIN_OR, glue, beard->mpicomm, 10,
        beard->comm_args, &aux_rates, &scale_rates, &intra_rhs, &inter_rhs,
        &add_rates, lua_get_global_int(prefs->L, "RK_init", 1));
  else if (prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
    beard->beard_ts = (bfam_ts_t *)bfam_ts_local_adams_new(
        (bfam_domain_t *)beard->domain, prefs->local_adams_method, num_lvl,
        BFAM_DOMAIN_OR, timestep_tags, BFAM_DOMAIN_OR, glue, beard->mpicomm, 10,
        beard->comm_args, &aux_rates, &glue_rates, &scale_rates, &intra_rhs,
        &inter_rhs, &add_rates, &add_rates_glue_p,
        lua_get_global_int(prefs->L, "RK_init", 1));
}

static void compute_subdomain_dt(bfam_locidx_t npoints, const char *name,
                                 bfam_real_t time, bfam_real_t *restrict x,
                                 bfam_real_t *restrict y,
                                 bfam_real_t *restrict z,
                                 struct bfam_subdomain *s, void *arg,
                                 bfam_real_t *restrict JI)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(JI, 32);
  bfam_real_t *dt = (bfam_real_t *)arg;

  bfam_real_t *restrict Jr0x0 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr0x0");
  BFAM_ASSUME_ALIGNED(Jr0x0, 32);

  bfam_real_t *restrict Jr0x1 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr0x1");
  BFAM_ASSUME_ALIGNED(Jr0x1, 32);

#if DIM == 3
  bfam_real_t *restrict Jr0x2 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr0x2");
  BFAM_ASSUME_ALIGNED(Jr0x2, 32);
#endif

  bfam_real_t *restrict Jr1x0 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr1x0");
  BFAM_ASSUME_ALIGNED(Jr1x0, 32);

  bfam_real_t *restrict Jr1x1 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr1x1");
  BFAM_ASSUME_ALIGNED(Jr1x1, 32);

#if DIM == 3
  bfam_real_t *restrict Jr1x2 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr1x2");
  BFAM_ASSUME_ALIGNED(Jr1x2, 32);
#endif

#if DIM == 3
  bfam_real_t *restrict Jr2x0 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr2x0");
  BFAM_ASSUME_ALIGNED(Jr2x0, 32);

  bfam_real_t *restrict Jr2x1 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr2x1");
  BFAM_ASSUME_ALIGNED(Jr2x1, 32);

  bfam_real_t *restrict Jr2x2 =
      bfam_dictionary_get_value_ptr(&s->fields, "_grid_Jr2x2");
  BFAM_ASSUME_ALIGNED(Jr2x2, 32);
#endif

  bfam_real_t *restrict mu = bfam_dictionary_get_value_ptr(&s->fields, "mu");
  BFAM_ASSUME_ALIGNED(mu, 32);

  bfam_real_t *restrict lam = bfam_dictionary_get_value_ptr(&s->fields, "lam");
  BFAM_ASSUME_ALIGNED(lam, 32);

  bfam_real_t *restrict rho = bfam_dictionary_get_value_ptr(&s->fields, "rho");
  BFAM_ASSUME_ALIGNED(rho, 32);

  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)s;
  bfam_real_t p2 = sub->N * sub->N;

  for (int n = 0; n < npoints; n++)
  {
    bfam_real_t cp = BFAM_REAL_SQRT((lam[n] + 2 * mu[n]) / rho[n]);
#if DIM == 2
    bfam_real_t hr1 =
        BFAM_REAL(1.0) / BFAM_REAL_SQRT(JI[n] * JI[n] * (Jr1x1[n] * Jr1x1[n] +
                                                         Jr1x0[n] * Jr1x0[n]));
    bfam_real_t hr0 =
        BFAM_REAL(1.0) / BFAM_REAL_SQRT(JI[n] * JI[n] * (Jr0x1[n] * Jr0x1[n] +
                                                         Jr0x0[n] * Jr0x0[n]));
    dt[0] = BFAM_MIN(dt[0], BFAM_MIN(hr1, hr0) / cp / p2);
#elif DIM == 3
    bfam_real_t hr2 =
        BFAM_REAL(1.0) / BFAM_REAL_SQRT(JI[n] * JI[n] * (Jr2x2[n] * Jr2x2[n] +
                                                         Jr2x1[n] * Jr2x1[n] +
                                                         Jr2x0[n] * Jr2x0[n]));
    bfam_real_t hr1 =
        BFAM_REAL(1.0) / BFAM_REAL_SQRT(JI[n] * JI[n] * (Jr1x2[n] * Jr1x2[n] +
                                                         Jr1x1[n] * Jr1x1[n] +
                                                         Jr1x0[n] * Jr1x0[n]));
    bfam_real_t hr0 =
        BFAM_REAL(1.0) / BFAM_REAL_SQRT(JI[n] * JI[n] * (Jr0x2[n] * Jr0x2[n] +
                                                         Jr0x1[n] * Jr0x1[n] +
                                                         Jr0x0[n] * Jr0x0[n]));
    dt[0] = BFAM_MIN(dt[0], BFAM_MIN(BFAM_MIN(hr2, hr1), hr0) / cp / p2);
#else
#error "bad dimension"
#endif
  }
}

static bfam_real_t compute_energy(beard_t *beard, prefs_t *prefs, bfam_real_t t,
                                  const char *prefix)
{
  const char *tags[] = {"_volume", NULL};
  bfam_subdomain_t *subs[beard->domain->base.numSubdomains];
  bfam_locidx_t num_subs = 0;
  bfam_domain_get_subdomains((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                             tags, beard->domain->base.numSubdomains, subs,
                             &num_subs);
  bfam_real_t energy = 0;
  bfam_real_t energy_local = 0;
  for (bfam_locidx_t s = 0; s < num_subs; s++)
  {
    bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subs[s];
#if DIM == 2
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_energy_2_##order(sub->N, &energy_local, sub, prefix);            \
    break;
#elif DIM == 3
#define X(order)                                                               \
  case order:                                                                  \
    beard_dgx_energy_3_##order(sub->N, &energy_local, sub, prefix);            \
    break;
#else
#error "Bad Dimension"
#endif

    switch (sub->N)
    {
      BFAM_LIST_OF_DGX_NORDERS
    default:
#if DIM == 2
      beard_dgx_energy_2_(sub->N, &energy_local, sub, prefix);
#elif DIM == 3
      beard_dgx_energy_3_(sub->N, &energy_local, sub, prefix);
#else
#error "Bad Dimension"
#endif
      break;
    }
#undef X
  }
  BFAM_MPI_CHECK(MPI_Reduce(&energy_local, &energy, 1, BFAM_REAL_MPI, MPI_SUM,
                            0, beard->mpicomm));
  if (beard->mpirank == 0)
    energy = BFAM_REAL_SQRT(energy);
  return energy;
}

typedef struct check_error_args
{
  char *field_prefix;
  lua_State *L;
} check_error_args_t;

static void check_error(bfam_locidx_t npoints, const char *name, bfam_real_t t,
                        bfam_real_t *restrict x, bfam_real_t *restrict y,
                        bfam_real_t *restrict z, struct bfam_subdomain *s,
                        void *arg, bfam_real_t *restrict err)
{
  BFAM_ASSUME_ALIGNED(x, 32);
  BFAM_ASSUME_ALIGNED(y, 32);
  BFAM_ASSUME_ALIGNED(z, 32);
  BFAM_ASSUME_ALIGNED(err, 32);

  check_error_args_t *err_args = (check_error_args_t *)arg;
  lua_State *L = err_args->L;
  lua_getglobal(L, name + 6);
  BFAM_ABORT_IF_NOT(lua_isfunction(L, -1),
                    "no callback function for initial condition and error: %s",
                    name + 6);
  lua_pop(L, 1);
  char fname[BFAM_BUFSIZ];
  snprintf(fname, BFAM_BUFSIZ, "%s%s", err_args->field_prefix, name + 6);
  bfam_real_t *restrict fld = bfam_dictionary_get_value_ptr(&s->fields, fname);
  BFAM_ABORT_IF(fld == NULL, "field '%s' not in fields for %s", fname, s->name);

#if DIM == 2
  bfam_real_t tmpz = 0;
#endif
  for (bfam_locidx_t n = 0; n < npoints; ++n)
  {
#if DIM == 2
    lua_global_function_call(L, 0, name + 6, "rrrr>r", x[n], y[n], tmpz, t,
                             &err[n]);
#elif DIM == 3
    lua_global_function_call(L, 0, name + 6, "rrrr>r", x[n], y[n], z[n], t,
                             &err[n]);
#else
#error "Bad Dimension"
#endif
    err[n] -= fld[n];
  }
}

typedef struct station_args
{
  prefs_t *prefs;
  const char **fields;
  bfam_real_t time;
} station_args_t;

static int beard_output_stations_flush(const char *key, void *val,
                                       void *in_args)
{
  bfam_subdomain_dgx_point_interp_t *point =
      (bfam_subdomain_dgx_point_interp_t *)val;

  fflush(point->file);

  return 0;
}

static int beard_output_stations(const char *key, void *val, void *in_args)
{
  bfam_subdomain_dgx_point_interp_t *point =
      (bfam_subdomain_dgx_point_interp_t *)val;

  station_args_t *args = (station_args_t *)in_args;
  BFAM_ASSERT(args);

  const char **fields = args->fields;

  if (point->num_interp == 1)
    bfam_subdomain_dgx_point_interp_fields_1(point, args->time, "", fields, 1);
  else if (point->num_interp == 2)
    bfam_subdomain_dgx_point_interp_fields_2(point, args->time, "", fields, 2);
  else if (point->num_interp == 3)
    bfam_subdomain_dgx_point_interp_fields_3(point, args->time, "", fields, 3);
  else
    BFAM_ABORT("Invalid number of interps %d (%s)", point->num_interp,
               point->filename);

  return 0;
}

static int beard_open_body_stations(const char *key, void *val, void *in_args)
{
  bfam_subdomain_dgx_point_interp_t *point =
      (bfam_subdomain_dgx_point_interp_t *)val;
  bfam_subdomain_dgx_point_interp_open_(point);

  station_args_t *args = (station_args_t *)in_args;
  BFAM_ASSERT(args);

  prefs_t *prefs = args->prefs;
  const char **fields = args->fields;
  FILE *file = point->file;
  BFAM_ABORT_IF_NOT(file, "problem with opening file %s", point->filename);
  fprintf(file, "# problem = %s\n", prefs->output_prefix);
  fprintf(file, "# date    = XXX\n");
  fprintf(file, "# code    = beard %dd\n", DIM);
  fprintf(file, "# version = %s\n", bfam_version_get());
  fprintf(file, "# body station = %s\n", key);
  fprintf(file, "# The line below lists the names of the data fields:\n");

  /* Print the interpolated point */
  fprintf(file, "# interpolated point:\n#    (");
  for (int k = 0; k < DIM; k++)
  {
    char fld_name[BFAM_BUFSIZ];
    snprintf(fld_name, BFAM_BUFSIZ, "_grid_x%d", k);
    fprintf(file, "%" BFAM_REAL_PRIe,
            BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_,
                               DIM)(point, "", fld_name, FDIM));
    if (k + 1 < DIM)
      fprintf(file, ", ");
  }
  fprintf(file, ")\n");

  fprintf(file, "t");
  for (int n = 0; fields[n]; n++)
    fprintf(file, " %s", fields[n]);
  fprintf(file, "\n");

  return beard_output_stations(key, val, args);
}

static int beard_open_fault_stations(const char *key, void *val, void *in_args)
{
  bfam_subdomain_dgx_point_interp_t *point =
      (bfam_subdomain_dgx_point_interp_t *)val;
  bfam_subdomain_dgx_point_interp_open_(point);

  station_args_t *args = (station_args_t *)in_args;
  BFAM_ASSERT(args);

  prefs_t *prefs = args->prefs;
  const char **fields = args->fields;
  FILE *file = point->file;
  BFAM_ABORT_IF_NOT(file, "problem with opening file %s", point->filename);
  fprintf(file, "# problem = %s\n", prefs->output_prefix);
  fprintf(file, "# date    = XXX\n");
  fprintf(file, "# code    = beard %dd\n", DIM);
  fprintf(file, "# version = %s\n", bfam_version_get());
  fprintf(file, "# fault station = %s\n", key);
  fprintf(file, "# The line below lists the names of the data fields:\n");

  /* Print the interpolated point */
  fprintf(file, "# interpolated point:\n#    (");
  for (int k = 0; k < DIM; k++)
  {
    char fld_name[BFAM_BUFSIZ];
    snprintf(fld_name, BFAM_BUFSIZ, "_grid_x%d", k);
    fprintf(file, "%" BFAM_REAL_PRIe,
            BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_,
                               FDIM)(point, "", fld_name, FDIM));
    if (k + 1 < DIM)
      fprintf(file, ", ");
  }
  fprintf(file, ")\n");

  /* Print the interpolated point */
  fprintf(file, "# interpolated normal:\n#    (");
  for (int k = 0; k < DIM; k++)
  {
    char fld_name[BFAM_BUFSIZ];
    snprintf(fld_name, BFAM_BUFSIZ, "_grid_nx%d", k);
    fprintf(file, "%" BFAM_REAL_PRIe,
            BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_m_,
                               FDIM)(point, "", fld_name, FDIM));
    if (k + 1 < DIM)
      fprintf(file, ", ");
  }
  fprintf(file, ")\n");

  fprintf(file, "t");
  for (int n = 0; fields[n]; n++)
    fprintf(file, " %s", fields[n]);
  fprintf(file, "\n");

  return beard_output_stations(key, val, args);
}

static void time_level_comm_info(bfam_subdomain_t *thisSubdomain,
                                 size_t *send_sz, size_t *recv_sz,
                                 void *comm_args)
{
  BFAM_ASSERT(comm_args);
  *send_sz += sizeof(int);
  *recv_sz += sizeof(int);
}

static void time_level_put_send_buffer(bfam_subdomain_t *thisSubdomain,
                                       void *buffer, size_t send_sz,
                                       void *comm_args)
{
  /* just some sanity check */
  BFAM_ASSERT(send_sz == sizeof(int));
  BFAM_ASSERT(comm_args);
  BFAM_ASSERT(thisSubdomain->glue_m);
  BFAM_ASSERT(thisSubdomain->glue_m->sub_m);

  /* determine the level of the minus side and add tag to minus side */
  bfam_subdomain_comm_args_t *args = (bfam_subdomain_comm_args_t *)comm_args;
  BFAM_ASSERT(args->user_data);
  int max_lvl = *(int *)args->user_data;

  int lvl = 0;
  char tag[BFAM_BUFSIZ];
  for (; lvl <= max_lvl; lvl++)
  {
    bfam_ts_local_adams_fill_level_tag(tag, BFAM_BUFSIZ, lvl);
    if (bfam_subdomain_has_tag(thisSubdomain->glue_m->sub_m, tag))
      break;
  }
  BFAM_ABORT_IF(lvl > max_lvl, "glue %s: "
                               "max number of levels searched in minus side %s "
                               "and no level tag found",
                thisSubdomain->name, thisSubdomain->glue_m->sub_m->name);

  bfam_subdomain_minus_add_tag(thisSubdomain, tag);

  /* store the level */
  *(int *)buffer = lvl;
}

static void time_level_get_recv_buffer(bfam_subdomain_t *thisSubdomain,
                                       void *buffer, size_t recv_sz,
                                       void *comm_args)
{
  BFAM_ASSERT(comm_args);
  BFAM_ASSERT(recv_sz == sizeof(int));
  int lvl = *(int *)buffer;
  char tag[BFAM_BUFSIZ];
  bfam_ts_local_adams_fill_level_tag(tag, BFAM_BUFSIZ, lvl);
  bfam_subdomain_plus_add_tag(thisSubdomain, tag);
}

static bfam_real_t compute_domain_dt(beard_t *beard, prefs_t *prefs,
                                     const char *volume[], const char *glue[],
                                     int *nsteps_ptr, int *ndisp_ptr,
                                     int *noutput_ptr, int *nfoutput_ptr,
                                     int *nstations_ptr, int *nerr_ptr)
{

  /* first we get all the volume subdomains */
  bfam_subdomain_t **subdomains =
      bfam_malloc(((bfam_domain_t *)beard->domain)->numSubdomains *
                  sizeof(bfam_subdomain_t **));

  bfam_locidx_t numSubdomains = 0;

  bfam_domain_get_subdomains(((bfam_domain_t *)beard->domain), BFAM_DOMAIN_OR,
                             volume,
                             ((bfam_domain_t *)beard->domain)->numSubdomains,
                             subdomains, &numSubdomains);

  bfam_real_t ldt[numSubdomains];
  bfam_real_t min_ldt = INFINITY;

  for (bfam_locidx_t s = 0; s < numSubdomains; ++s)
  {
    ldt[s] = INFINITY;
    bfam_subdomain_field_init(subdomains[s], "_grid_JI", 0,
                              compute_subdomain_dt, &ldt[s]);
    BFAM_INFO("local dt for %s is %" BFAM_REAL_FMTe, subdomains[s]->name,
              ldt[s]);
    min_ldt = BFAM_MIN(min_ldt, ldt[s]);
  }

  BFAM_INFO("min local dt is %" BFAM_REAL_FMTe, min_ldt);

  bfam_real_t min_global_dt = 0;
  BFAM_MPI_CHECK(MPI_Allreduce(&min_ldt, &min_global_dt, 1, BFAM_REAL_MPI,
                               MPI_MIN, beard->mpicomm));

  bfam_real_t dt = 0;
  int result = lua_global_function_call(
      prefs->L, 0, "time_step_parameters", "r>riiiii", min_global_dt, &dt,
      nsteps_ptr, ndisp_ptr, noutput_ptr, nfoutput_ptr, nstations_ptr);
  BFAM_ABORT_IF_NOT(result == 0,
                    "problem with lua call to 'time_step_parameters': "
                    "should be a function that takes dt "
                    "and returns dt, nsteps, ndisp, noutput,nfoutput");
  result = lua_global_function_call(prefs->L, 0, "nerr", "r>i", dt, nerr_ptr);
  if (*nerr_ptr > 0)
  {
    const char *err_flds[] = {
        "error_v1",  "error_v2",  "error_v3",  "error_S11", "error_S22",
        "error_S33", "error_S12", "error_S13", "error_S23", NULL};
    for (int f = 0; err_flds[f] != NULL; f++)
      bfam_domain_add_field((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                            volume, err_flds[f]);
  }

  bfam_locidx_t num_time_lvl = 0;
  if (prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
  {
    bfam_real_t max_ldt = dt;

    bfam_long_real_t dt_fudge =
        (bfam_long_real_t)dt / (bfam_long_real_t)min_global_dt;

    int max_time_level = lua_get_global_int(prefs->L, "max_time_level", 32);
    for (bfam_locidx_t s = 0; s < numSubdomains; ++s)
    {
      ldt[s] = (bfam_real_t)dt_fudge * ldt[s];

      /* keep double time_step until the level is big enough */
      bfam_locidx_t time_level = 0;
      if (ldt[s] != INFINITY)
        for (; (1 << (time_level + 1)) * dt < ldt[s]; time_level++)
          ;
      time_level = BFAM_MIN(time_level, max_time_level);
      BFAM_ASSERT(time_level >= 0);

      char tag[BFAM_BUFSIZ];
      bfam_ts_local_adams_fill_level_tag(tag, BFAM_BUFSIZ, time_level);
      bfam_subdomain_add_tag(subdomains[s], tag);

      /* set this guys real dt */
      ldt[s] = dt * (1 << time_level);

      BFAM_INFO("For %s is time level %d with local dt %" BFAM_REAL_FMTe,
                subdomains[s]->name, time_level, ldt[s]);

      max_ldt = BFAM_MAX(max_ldt, ldt[s]);

      num_time_lvl = BFAM_MAX(num_time_lvl, time_level + 1);
    }
    BFAM_INFO("local number of time levels %" BFAM_LOCIDX_PRId, num_time_lvl);

    BFAM_INFO("max local dt is %" BFAM_REAL_FMTe, max_ldt);

    BFAM_MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, &num_time_lvl, 1,
                                 BFAM_LOCIDX_MPI, MPI_MAX, beard->mpicomm));

    dt = (1 << (num_time_lvl - 1)) * dt;

    BFAM_INFO("number of time levels %" BFAM_LOCIDX_PRId
              " for global dt %" BFAM_REAL_FMTe,
              num_time_lvl, dt);

    bfam_communicator_t level_comm;

    bfam_subdomain_comm_args_t level_args;
    const char *level_NULL[] = {NULL};
    level_args.scalars_m = level_NULL;
    level_args.vectors_m = level_NULL;
    level_args.vector_components_m = level_NULL;
    level_args.tensors_m = level_NULL;
    level_args.tensor_components_m = level_NULL;
    level_args.face_scalars_m = level_NULL;

    level_args.scalars_p = level_NULL;
    level_args.vectors_p = level_NULL;
    level_args.vector_components_p = level_NULL;
    level_args.tensors_p = level_NULL;
    level_args.tensor_components_p = level_NULL;
    level_args.face_scalars_p = level_NULL;

    level_args.user_comm_info = time_level_comm_info;
    level_args.user_put_send_buffer = time_level_put_send_buffer;
    level_args.user_get_recv_buffer = time_level_get_recv_buffer;
    level_args.user_data = &num_time_lvl;

    level_args.user_prefix_function = NULL;

    bfam_communicator_init(&level_comm, (bfam_domain_t *)beard->domain,
                           BFAM_DOMAIN_OR, glue, beard->mpicomm, 10,
                           &level_args);
    bfam_communicator_start(&level_comm);
    bfam_communicator_finish(&level_comm);
    bfam_communicator_free(&level_comm);
  }

  bfam_free(subdomains);

  init_time_stepper(beard, prefs, num_time_lvl);

  return dt;
}

static void run_simulation(beard_t *beard, prefs_t *prefs)
{
  const char *volume_tags[] = {"_volume", NULL};
  const char *glue_tags[] = {"_glue_parallel", "_glue_local", NULL};
  const char *friction_tags[] = {"friction", NULL};

  int nsteps = 0;
  int ndisp = 0;
  int noutput = 0;
  int nfoutput = 0;
  int nstations = 0;
  int nerr = 0;

  bfam_real_t dt =
      compute_domain_dt(beard, prefs, volume_tags, glue_tags, &nsteps, &ndisp,
                        &noutput, &nfoutput, &nstations, &nerr);

  BFAM_ROOT_INFO("dt        = %" BFAM_REAL_FMTe, dt);
  BFAM_ROOT_INFO("nsteps    = %d", nsteps);
  BFAM_ROOT_INFO("ndisp     = %d", ndisp);
  BFAM_ROOT_INFO("noutput   = %d", noutput);
  BFAM_ROOT_INFO("nfoutput  = %d", nfoutput);
  BFAM_ROOT_INFO("nstations = %d", nstations);
  BFAM_ROOT_INFO("nerr      = %d", nerr);

  if (nstations >= 0)
  {
    const char *volume_station_fields[] = {"v1", "v2", "v3", NULL};
    station_args_t volume_args;
    volume_args.prefs = prefs;
    volume_args.fields = volume_station_fields;
    volume_args.time = 0;
    bfam_dictionary_allprefixed_ptr(beard->volume_stations, "",
                                    &beard_open_body_stations, &volume_args);

    const char *fault_station_fields[] = {"Tp1", "Tp2", "Tp3", "Tn", "V",
                                          "Vp1", "Vp2", "Vp3", "Dp", "Dp1",
                                          "Dp2", "Dp3", "Dn",  NULL};
    station_args_t fault_args;
    fault_args.prefs = prefs;
    fault_args.fields = fault_station_fields;
    fault_args.time = 0;
    bfam_dictionary_allprefixed_ptr(beard->fault_stations, "",
                                    &beard_open_fault_stations, &fault_args);
  }

  if (nfoutput >= 0)
  {
    char output[BFAM_BUFSIZ];
    snprintf(output, BFAM_BUFSIZ, "%s_fault_%05d", prefs->output_prefix, 0);
    bfam_vtk_write_file((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                        friction_tags, prefs->data_directory, output, (0) * dt,
                        friction_fields, NULL, NULL, prefs->vtk_binary,
                        prefs->vtk_compress, prefs->vtk_num_pnts);
  }

  /* compute the initial energy */
  bfam_real_t initial_energy = compute_energy(beard, prefs, 0, "");
  bfam_real_t energy = initial_energy;
  if (initial_energy < BFAM_REAL_EPS)
    initial_energy = -1;

  if (noutput >= 0)
  {
    char output[BFAM_BUFSIZ];

    /* dump the pxest mesh */
    snprintf(output, BFAM_BUFSIZ, "%s/%s_pxest_mesh", prefs->data_directory,
             prefs->output_prefix);
    p4est_vtk_write_all(beard->domain->pxest, NULL, 1, 1, 1, 1, 0, 0, 0,
                        output);

    snprintf(output, BFAM_BUFSIZ, "%s_%05d", prefs->output_prefix, 0);
    bfam_vtk_write_file((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                        volume_vtk_tags, prefs->data_directory, output,
                        (0) * dt, elastic_fields, NULL, NULL, prefs->vtk_binary,
                        prefs->vtk_compress, prefs->vtk_num_pnts);

    if (plastic_fields)
    {
      snprintf(output, BFAM_BUFSIZ, "%s_plasticity", prefs->output_prefix);
      bfam_vtk_write_file(
          (bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR, volume_vtk_tags,
          prefs->data_directory, output, (0) * dt, plastic_fields, NULL, NULL,
          prefs->vtk_binary, prefs->vtk_compress, prefs->vtk_num_pnts);
    }
  }

  BFAM_ASSERT(beard->beard_ts);
  const int NFLUSH = 10;
  int nflush = 0;
  const bfam_real_t first_step = bfam_clock();

  BFAM_ROOT_INFO("\x1B[%dm"
                 "time: %10.5" BFAM_REAL_PRIe " energy: %10.5" BFAM_REAL_PRIe
                 "\x1B[0m",
                 34, 0 * dt, energy);
  bfam_real_t time = 0;
  for (int s = 1; s <= nsteps; s++)
  {
    beard->beard_ts->step(beard->beard_ts, dt);
    time = s * dt;
    if (plastic_fields)
    {
      /* Use the return mapping algorithm to handle the plasticity */
      switch (prefs->plasticity)
      {
      case DUVAUT_LIONS:
      {
        /* get the subdomains being handled by DL plasticity */
        const char *DL_plastic[] = {DUVAUT_LIONS_TAG, NULL};
        bfam_subdomain_t *subs[beard->domain->base.numSubdomains];
        bfam_locidx_t num_subs = 0;
        bfam_domain_get_subdomains(
            (bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR, DL_plastic,
            beard->domain->base.numSubdomains, subs, &num_subs);

        BFAM_ABORT_IF_NOT(prefs->lsrk_method != BFAM_TS_LSRK_NOOP,
                          "right now plasticity requires LSRK time stepping");

        /* call the return map algorithm on these subdomains */
        for (bfam_locidx_t n = 0; n < num_subs; n++)
        {
          bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subs[n];
          duvaut_lions_return_map(sub->N, sub, "", time, dt);
        }
        break;
      }
      default:
        BFAM_ABORT("Unknown plastic type for adding the tags");
      }
    }

    /* Compute the rupture time */
    {
      bfam_subdomain_t *subs[beard->domain->base.numSubdomains];
      bfam_locidx_t num_subs = 0;
      bfam_domain_get_subdomains(
          (bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR, friction_tags,
          beard->domain->base.numSubdomains, subs, &num_subs);

      /* call the return map algorithm on these subdomains */
      for (bfam_locidx_t n = 0; n < num_subs; n++)
      {
        bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subs[n];

        /* Vrup is currently hard coded to 0.001*/
        rupture_time(sub->N, sub, "", time, 0.001);
      }
    }

    if (ndisp > 0 && s % ndisp == 0)
    {
      bfam_real_t new_energy = compute_energy(beard, prefs, time, "");
      if (initial_energy < 0)
      {
        BFAM_ROOT_INFO("\x1B[%dm"
                       "time: %10.5" BFAM_REAL_PRIe
                       " (step %8d of %8d with %10.5" BFAM_REAL_PRIe
                       " sec/step)"
                       "\x1B[0m",
                       34, time, s, nsteps, (bfam_clock() - first_step) / s);
      }
      else
      {
        int color = 32;
        if (new_energy > energy)
          color = 31;
        BFAM_ROOT_INFO(
            "\x1B[%dm"
            "time: %" BFAM_REAL_FMTe "\n"
            " init energy: %" BFAM_REAL_FMTe " energy: %" BFAM_REAL_FMTe
            " norm energy: %" BFAM_REAL_FMTe "\n"
            " current delta energy: %+" BFAM_REAL_FMTe
            " initial delta energy: %+" BFAM_REAL_FMTe "\n"
            "(%10.5" BFAM_REAL_PRIe "sec/step)\n"
            "\x1B[0m",
            color, time, initial_energy, new_energy,
            new_energy / initial_energy, (new_energy - energy) / initial_energy,
            energy / initial_energy - 1, (bfam_clock() - first_step) / s);
      }
      energy = new_energy;
    }
    if (nstations > 0 && s % nstations == 0)
    {
      const char *volume_station_fields[] = {"v1", "v2", "v3", NULL};
      station_args_t volume_args;
      volume_args.prefs = prefs;
      volume_args.fields = volume_station_fields;
      volume_args.time = time;
      bfam_dictionary_allprefixed_ptr(beard->volume_stations, "",
                                      &beard_output_stations, &volume_args);
      if (nflush % NFLUSH == 0)
        bfam_dictionary_allprefixed_ptr(beard->volume_stations, "",
                                        &beard_output_stations_flush,
                                        &volume_args);

      const char *fault_station_fields[] = {"Tp1", "Tp2", "Tp3", "Tn", "V",
                                            "Vp1", "Vp2", "Vp3", "Dp", "Dp1",
                                            "Dp2", "Dp3", "Dn",  NULL};
      station_args_t fault_args;
      fault_args.prefs = prefs;
      fault_args.fields = fault_station_fields;
      fault_args.time = time;
      bfam_dictionary_allprefixed_ptr(beard->fault_stations, "",
                                      &beard_output_stations, &fault_args);

      if (nflush % NFLUSH == 0)
        bfam_dictionary_allprefixed_ptr(beard->fault_stations, "",
                                        &beard_output_stations_flush,
                                        &fault_args);

      nflush++;
    }
    if (nfoutput > 0 && s % nfoutput == 0)
    {
      char output[BFAM_BUFSIZ];
      snprintf(output, BFAM_BUFSIZ, "%s_fault_%05d", prefs->output_prefix, s);
      bfam_vtk_write_file((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                          friction_tags, prefs->data_directory, output, time,
                          friction_fields, NULL, NULL, prefs->vtk_binary,
                          prefs->vtk_compress, prefs->vtk_num_pnts);
    }
    if (noutput > 0 && s % noutput == 0)
    {
      const char *fields[] = {"v1",  "v2",  "v3",  "S11", "S22",
                              "S33", "S12", "S13", "S23", NULL};
      char output[BFAM_BUFSIZ];
      snprintf(output, BFAM_BUFSIZ, "%s_%05d", prefs->output_prefix, s);
      bfam_vtk_write_file((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                          volume_vtk_tags, prefs->data_directory, output, time,
                          fields, NULL, NULL, prefs->vtk_binary,
                          prefs->vtk_compress, prefs->vtk_num_pnts);
    }
    if (nerr > 0 && s % nerr == 0)
    {
      check_error_args_t err_args;
      err_args.L = prefs->L;
      char prefix[] = "";
      err_args.field_prefix = prefix;
      const char *err_flds[] = {
          "error_v1",  "error_v2",  "error_v3",  "error_S11", "error_S22",
          "error_S33", "error_S12", "error_S13", "error_S23", NULL};
      for (int f = 0; err_flds[f] != NULL; f++)
        bfam_domain_init_field((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                               volume_vtk_tags, err_flds[f], time, check_error,
                               &err_args);
      bfam_real_t error = compute_energy(beard, prefs, time, "error_");
      bfam_real_t new_energy = compute_energy(beard, prefs, time, "");
      BFAM_ROOT_INFO("time: %" BFAM_REAL_FMTe " error: %" BFAM_REAL_FMTe
                     " d_energy: %" BFAM_REAL_FMTe,
                     time, error, (new_energy - initial_energy) / initial_energy);
      if (noutput > 0)
      {
        char err_output[BFAM_BUFSIZ];
        snprintf(err_output, BFAM_BUFSIZ, "%s_error_%05d", prefs->output_prefix,
                 s);
        bfam_vtk_write_file((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                            volume_vtk_tags, prefs->data_directory, err_output,
                            time, err_flds, NULL, NULL, prefs->vtk_binary,
                            prefs->vtk_compress, prefs->vtk_num_pnts);
        energy = new_energy;
      }
    }
  }
  if (prefs->output_file_fault)
  {
    char output[BFAM_BUFSIZ];
    snprintf(output, BFAM_BUFSIZ, "%s_final_fault", prefs->output_prefix);
    bfam_vtk_write_file((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                        friction_tags, prefs->data_directory, output, time,
                        friction_fields, NULL, NULL, prefs->vtk_binary,
                        prefs->vtk_compress, prefs->vtk_num_pnts);
  }
}

static int beard_free_stations(const char *key, void *val, void *arg)
{
  bfam_subdomain_dgx_point_interp_t *point =
      (bfam_subdomain_dgx_point_interp_t *)val;
  bfam_subdomain_dgx_point_interp_free_(point);
  bfam_free(point);
  return 0;
}

static void shave_beard(beard_t *beard, prefs_t *prefs)
{
  if (beard->volume_stations)
  {
    bfam_dictionary_allprefixed_ptr(beard->volume_stations, "",
                                    &beard_free_stations, NULL);
    bfam_dictionary_clear(beard->volume_stations);
    bfam_free(beard->volume_stations);
    beard->volume_stations = NULL;
  }
  if (beard->fault_stations)
  {
    bfam_dictionary_allprefixed_ptr(beard->fault_stations, "",
                                    &beard_free_stations, NULL);
    bfam_dictionary_clear(beard->fault_stations);
    bfam_free(beard->fault_stations);
    beard->fault_stations = NULL;
  }
  bfam_free(beard->comm_args);
  if (prefs->lsrk_method != BFAM_TS_LSRK_NOOP)
    bfam_ts_lsrk_free((bfam_ts_lsrk_t *)beard->beard_ts);
  if (prefs->adams_method != BFAM_TS_ADAMS_NOOP)
    bfam_ts_adams_free((bfam_ts_adams_t *)beard->beard_ts);
  if (prefs->local_adams_method != BFAM_TS_LOCAL_ADAMS_NOOP)
    bfam_ts_local_adams_free((bfam_ts_local_adams_t *)beard->beard_ts);
  bfam_free(beard->beard_ts);
  bfam_domain_pxest_free(beard->domain);
  bfam_free(beard->domain);
  p4est_connectivity_destroy(beard->conn);
}

static void init_fault_stations(beard_t *beard, prefs_t *prefs)
{
  lua_State *L = prefs->L;
#ifdef BFAM_DEBUG
  int top = lua_gettop(L);
#endif

  /* First we get the station info from lua */
  lua_getglobal(L, "fault_stations");

  BFAM_ASSERT(beard->fault_stations == NULL);
  beard->fault_stations = bfam_malloc(sizeof(bfam_dictionary_t));
  bfam_dictionary_init(beard->fault_stations);

  if (!lua_istable(L, -1))
    BFAM_ROOT_WARNING("table `%s' not found", "fault_stations");
  else
  {

    bfam_locidx_t nStationArgs = 1 + DIM + DIM + 1;

    const int length_stations = lua_objlen(L, -1);
    BFAM_LDEBUG("fault_stations  #elem: %3d", length_stations);

    BFAM_ABORT_IF_NOT(
        length_stations % nStationArgs == 0,
        "length of fault_stations should be a multiple of %" BFAM_LOCIDX_PRId,
        nStationArgs);

    const bfam_locidx_t num_stations = length_stations / nStationArgs;

    bfam_real_t xyz[DIM * num_stations];
    bfam_real_t nxf[DIM * num_stations];
    bfam_real_t tol[num_stations];

    char station_names[num_stations][BFAM_BUFSIZ];

    luaL_checktype(L, -1, LUA_TTABLE);

    for (bfam_locidx_t i = 0; i < num_stations; i++)
    {
      bfam_locidx_t offset = 1;

      /* Get the station name */
      lua_rawgeti(L, -1, nStationArgs * i + offset);
      BFAM_ABORT_IF_NOT(lua_isstring(L, -1),
                        "stations %d field 1 is not a string", i);
      strncpy(station_names[i], lua_tostring(L, -1), BFAM_BUFSIZ);
      lua_pop(L, 1);

      /* Get the station coordinates */
      for (bfam_locidx_t j = 0; j < DIM; j++)
      {
        offset++;
        lua_rawgeti(L, -1, nStationArgs * i + offset);
        BFAM_ABORT_IF_NOT(lua_isnumber(L, -1),
                          "stations %d field %d is not a number", i, offset);
        xyz[DIM * i + j] = (bfam_real_t)lua_tonumber(L, -1);
        lua_pop(L, 1);
      }

      /* Get the station normal projection direction */
      for (bfam_locidx_t j = 0; j < DIM; j++)
      {
        offset++;
        lua_rawgeti(L, -1, nStationArgs * i + offset);
        BFAM_ABORT_IF_NOT(lua_isnumber(L, -1),
                          "stations %d field %d is not a number", i, offset);
        nxf[DIM * i + j] = (bfam_real_t)lua_tonumber(L, -1);
        lua_pop(L, 1);
      }

      /* Get the station tolerance */
      offset++;
      lua_rawgeti(L, -1, nStationArgs * i + offset);
      BFAM_ABORT_IF_NOT(lua_isnumber(L, -1),
                        "stations %d field %d is not a number", i, offset);
      tol[i] = (bfam_real_t)lua_tonumber(L, -1);
      lua_pop(L, 1);

      BFAM_ASSERT(offset == nStationArgs);

#if DIM == 2
      BFAM_ROOT_LDEBUG("Station %04d (%s):"
                       " %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe
                       " %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe,
                       i, station_names[i], xyz[DIM * i + 0], xyz[DIM * i + 1],
                       nxf[DIM * i + 0], nxf[DIM * i + 1]);
#elif DIM == 3
      BFAM_ROOT_LDEBUG("Station %04d (%s):"
                       " %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe
                       " %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe
                       " %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe,
                       i, station_names[i], xyz[DIM * i + 0], xyz[DIM * i + 1],
                       xyz[DIM * i + 2], nxf[DIM * i + 0], nxf[DIM * i + 1],
                       nxf[DIM * i + 2]);
#endif
    }

    /* Now we handle finding where we interp */

    const char *fault[] = {"friction", NULL};
    bfam_subdomain_t *subs[beard->domain->base.numSubdomains];
    bfam_locidx_t num_subs = 0;
    bfam_domain_get_subdomains((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                               fault, beard->domain->base.numSubdomains, subs,
                               &num_subs);

    for (bfam_locidx_t n = 0; n < num_subs; n++)
    {
      bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t *)subs[n];
      const bfam_locidx_t num_cells = s->K;

      BFAM_LOAD_FIELD_RESTRICT_ALIGNED(x, "", "_grid_x0", &s->base.fields);
      BFAM_LOAD_FIELD_RESTRICT_ALIGNED(y, "", "_grid_x1", &s->base.fields);
#if DIM == 3
      BFAM_LOAD_FIELD_RESTRICT_ALIGNED(z, "", "_grid_x2", &s->base.fields);
#endif

      int **msk = s->gmask[s->numg - 1];

      for (bfam_locidx_t k = 0; k < num_cells; k++)
      {
        const bfam_real_t *x_e = x + k * s->Np;
        const bfam_real_t *y_e = y + k * s->Np;
        BEARD_D3_OP(const bfam_real_t *z_e = z + k * s->Np);
        for (bfam_locidx_t i = 0; i < num_stations; i++)
        {
#if DIM == 2
          bfam_real_t r[3];
          /* inverse_linear(r, xyz[i*DIM+0], x_e[msk[0][0]], x_e[msk[1][0]]);*/
          inverse_bilinear_normal(
              r, xyz[i * DIM + 0], xyz[i * DIM + 1], 0, nxf[i * DIM + 0],
              nxf[i * DIM + 1], 0, x_e[msk[0][0]], x_e[msk[1][0]],
              x_e[msk[0][0]], x_e[msk[1][0]], y_e[msk[0][0]], y_e[msk[1][0]],
              y_e[msk[0][0]], y_e[msk[1][0]], -1, -1, 1, 1);
#elif DIM == 3
          bfam_real_t r[3];

          inverse_bilinear_normal(
              r, xyz[i * DIM + 0], xyz[i * DIM + 1], xyz[i * DIM + 2],
              nxf[i * DIM + 0], nxf[i * DIM + 1], nxf[i * DIM + 2],
              x_e[msk[0][0]], x_e[msk[1][0]], x_e[msk[2][0]], x_e[msk[3][0]],
              y_e[msk[0][0]], y_e[msk[1][0]], y_e[msk[2][0]], y_e[msk[3][0]],
              z_e[msk[0][0]], z_e[msk[1][0]], z_e[msk[2][0]], z_e[msk[3][0]]);
#else
#error "Bad Dimension"
#endif

          /* if we are in the area, build the interpolant */
          if (BEARD_D3_AP(BFAM_REAL_ABS(r[0]) <= 1 + 100 * BFAM_REAL_EPS,
                          &&BFAM_REAL_ABS(r[1]) <= 1 + 100 * BFAM_REAL_EPS) &&
              BFAM_REAL_ABS(r[2]) <= tol[i])
          {
            /* Determine the station name */
            char filename[BFAM_BUFSIZ];
            snprintf(filename, BFAM_BUFSIZ, "%s/%s_%s_%s_%010d.dat",
                     prefs->data_directory, prefs->output_prefix,
                     station_names[i], s->base.name, k);

            /* Create the point and store it */
            bfam_subdomain_dgx_point_interp_t *point =
                BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_new_,
                                   FDIM)(s, k, r, filename, BFAM_BUFSIZ, FDIM);
            bfam_dictionary_insert_ptr(beard->fault_stations, filename, point);
          }
        }
      }
    }
  }

  lua_pop(L, 1);
  BFAM_ASSERT(top == lua_gettop(L));
}

static void init_volume_stations(beard_t *beard, prefs_t *prefs)
{
  lua_State *L = prefs->L;
#ifdef BFAM_DEBUG
  int top = lua_gettop(L);
#endif

  /* First we get the station info from lua */
  lua_getglobal(L, "volume_stations");

  BFAM_ASSERT(beard->volume_stations == NULL);
  beard->volume_stations = bfam_malloc(sizeof(bfam_dictionary_t));
  bfam_dictionary_init(beard->volume_stations);

  if (!lua_istable(L, -1))
    BFAM_ROOT_WARNING("table `%s' not found", "volume_stations");
  else
  {

    const int length_stations = lua_objlen(L, -1);
    BFAM_LDEBUG("volume_stations  #elem: %3d", length_stations);

    BFAM_ABORT_IF_NOT(length_stations % (DIM + 1) == 0,
                      "length of volume_stations should be a multiple of %d",
                      (DIM + 1));

    const int num_stations = length_stations / (DIM + 1);

    bfam_real_t xyz[DIM * num_stations];
    char station_names[num_stations][BFAM_BUFSIZ];

    luaL_checktype(L, -1, LUA_TTABLE);

    for (int i = 0; i < num_stations; i++)
    {
      lua_rawgeti(L, -1, (DIM + 1) * i + 1);
      BFAM_ABORT_IF_NOT(lua_isstring(L, -1),
                        "stations %d field 1 is not a string", i);
      strncpy(station_names[i], lua_tostring(L, -1), BFAM_BUFSIZ);
      lua_pop(L, 1);
      for (int j = 0; j < DIM; j++)
      {
        lua_rawgeti(L, -1, (DIM + 1) * i + j + 2);
        BFAM_ABORT_IF_NOT(lua_isnumber(L, -1),
                          "stations %d field %d is not a number", i, j + 1);
        xyz[DIM * i + j] = (bfam_real_t)lua_tonumber(L, -1);
        lua_pop(L, 1);
      }

#if DIM == 2
      BFAM_ROOT_LDEBUG("Station %04d (%s):"
                       " %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe,
                       i, station_names[i], xyz[DIM * i + 0], xyz[DIM * i + 1]);
#elif DIM == 3
      BFAM_ROOT_LDEBUG("Station %04d (%s):"
                       " %" BFAM_REAL_FMTe " %" BFAM_REAL_FMTe
                       " %" BFAM_REAL_FMTe,
                       i, station_names[i], xyz[DIM * i + 0], xyz[DIM * i + 1],
                       xyz[DIM * i + 2]);
#endif
    }

    /* Now we handle finding where we interp */

    const char *volume[] = {"_volume", NULL};
    bfam_subdomain_t *subs[beard->domain->base.numSubdomains];
    bfam_locidx_t num_subs = 0;
    bfam_domain_get_subdomains((bfam_domain_t *)beard->domain, BFAM_DOMAIN_OR,
                               volume, beard->domain->base.numSubdomains, subs,
                               &num_subs);

    for (int n = 0; n < num_subs; n++)
    {
      bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t *)subs[n];
      const int num_cells = s->K;

      BFAM_LOAD_FIELD_RESTRICT_ALIGNED(x, "", "_grid_x0", &s->base.fields);
      BFAM_LOAD_FIELD_RESTRICT_ALIGNED(y, "", "_grid_x1", &s->base.fields);
#if DIM == 3
      BFAM_LOAD_FIELD_RESTRICT_ALIGNED(z, "", "_grid_x2", &s->base.fields);
#endif

      int **msk = s->gmask[s->numg - 1];

      for (int k = 0; k < num_cells; k++)
      {
        const bfam_real_t *x_e = x + k * s->Np;
        const bfam_real_t *y_e = y + k * s->Np;
        BEARD_D3_OP(const bfam_real_t *z_e = z + k * s->Np);
        for (int i = 0; i < num_stations; i++)
        {
#if DIM == 2
          bfam_real_t r[2];
          inverse_bilinear(r, xyz[i * DIM + 0], xyz[i * DIM + 1],
                           x_e[msk[0][0]], x_e[msk[1][0]], x_e[msk[2][0]],
                           x_e[msk[3][0]], y_e[msk[0][0]], y_e[msk[1][0]],
                           y_e[msk[2][0]], y_e[msk[3][0]]);
#elif DIM == 3
          bfam_real_t r[3];
          inverse_trilinear(
              r, xyz[i * DIM + 0], xyz[i * DIM + 1], xyz[i * DIM + 2],
              x_e[msk[0][0]], x_e[msk[1][0]], x_e[msk[2][0]], x_e[msk[3][0]],
              x_e[msk[4][0]], x_e[msk[5][0]], x_e[msk[6][0]], x_e[msk[7][0]],
              y_e[msk[0][0]], y_e[msk[1][0]], y_e[msk[2][0]], y_e[msk[3][0]],
              y_e[msk[4][0]], y_e[msk[5][0]], y_e[msk[6][0]], y_e[msk[7][0]],
              z_e[msk[0][0]], z_e[msk[1][0]], z_e[msk[2][0]], z_e[msk[3][0]],
              z_e[msk[4][0]], z_e[msk[5][0]], z_e[msk[6][0]], z_e[msk[7][0]]);
#else
#error "Bad Dimension"
#endif

          /* if we are in the area, build the interpolant */
          if (BEARD_D3_AP(BFAM_REAL_ABS(r[0]) <= 1 + BFAM_REAL_EPS &&
                              BFAM_REAL_ABS(r[1]) <= 1 + BFAM_REAL_EPS,
                          &&BFAM_REAL_ABS(r[2]) <= 1 + BFAM_REAL_EPS))
          {
            char filename[BFAM_BUFSIZ];
            snprintf(filename, BFAM_BUFSIZ, "%s/%s_%s_%s_%010d.dat",
                     prefs->data_directory, prefs->output_prefix,
                     station_names[i], s->base.name, k);
            bfam_subdomain_dgx_point_interp_t *point =
                BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_new_,
                                   DIM)(s, k, r, filename, BFAM_BUFSIZ, DIM);

#ifdef BFAM_DEBUG
            bfam_real_t x =
                BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_,
                                   DIM)(point, "", "_grid_x0", DIM);
            BFAM_ABORT_IF_NOT(BFAM_REAL_APPROX_EQ(x, xyz[i * DIM + 0], 10),
                              "problem with the station %s interpolation in x."
                              "got: %" BFAM_REAL_FMTe
                              " expected: %" BFAM_REAL_FMTe,
                              station_names[i], x, xyz[i * DIM + 0]);

            bfam_real_t y =
                BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_,
                                   DIM)(point, "", "_grid_x1", DIM);
            BFAM_ABORT_IF_NOT(BFAM_REAL_APPROX_EQ(y, xyz[i * DIM + 1], 10),
                              "problem with the station %s interpolation in y."
                              " got: %" BFAM_REAL_FMTe
                              " expected: %" BFAM_REAL_FMTe,
                              station_names[i], y, xyz[i * DIM + 1]);

            BEARD_D3_OP(
                bfam_real_t z =
                    BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_,
                                       DIM)(point, "", "_grid_x2", DIM);
                BFAM_ABORT_IF_NOT(
                    BFAM_REAL_APPROX_EQ(z, xyz[i * DIM + 2], 100),
                    "problem with the station %s interpolation in z."
                    "got: %" BFAM_REAL_FMTe " expected: %" BFAM_REAL_FMTe,
                    station_names[i], z, xyz[i * DIM + 2]));
#endif

            bfam_dictionary_insert_ptr(beard->volume_stations, filename, point);
          }
        }
      }
    }
  }

  lua_pop(L, 1);
  BFAM_ASSERT(top == lua_gettop(L));
}

/*
 * run the beard
 */
static void run(MPI_Comm mpicomm, prefs_t *prefs)
{
  beard_t beard;
  beard.conn = NULL;
  beard.domain = NULL;
  beard.beard_ts = NULL;
  beard.comm_args = NULL;
  beard.volume_stations = NULL;
  beard.fault_stations = NULL;

  init_mpi(&beard, mpicomm);

  init_domain(&beard, prefs);

  init_volume_stations(&beard, prefs);

  init_fault_stations(&beard, prefs);

  run_simulation(&beard, prefs);

  shave_beard(&beard, prefs);
}

int main(int argc, char *argv[])
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;

  void *options = bfam_gopt_sort(
      &argc, (const char **)argv,
      bfam_gopt_start(bfam_gopt_option('h', 0, bfam_gopt_shorts('h', '?'),
                                       bfam_gopt_longs("help", "HELP")),
                      bfam_gopt_option('V', 0, bfam_gopt_shorts('V'),
                                       bfam_gopt_longs("version")),
                      bfam_gopt_option('v', BFAM_GOPT_REPEAT,
                                       bfam_gopt_shorts('v'),
                                       bfam_gopt_longs("verbose"))));

  const char *help_text =
      "  %s [options] prefs_file\n"
      "\n"
      "  there are four possible options to this program, some of which have\n"
      "  multiple names:\n"
      "\n"
      "    -h -? --help --HELP\n"
      "    -V --version\n"
      "    -v --verbose  (which may be repeated for more verbosity)\n"
      "\n";

  if (bfam_gopt(options, 'h'))
  {
    /*
     * if any of the help options was specified
     */
    BFAM_ROOT_INFO(help_text, argv[0]);
    exit(EXIT_SUCCESS);
  }

  if (bfam_gopt(options, 'V'))
  {
    BFAM_ROOT_INFO("BFAM Version: %s", bfam_version_get());
    BFAM_ROOT_INFO("BFAM Compile Info:\n" BFAM_COMPILE_INFO);
    exit(EXIT_SUCCESS);
  }

  int verbosity = bfam_gopt(options, 'v');

  BFAM_MPI_CHECK(MPI_Init(&argc, &argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  if (argc != 2)
  {
    BFAM_LERROR("Unexpected number of arguments.");
    BFAM_ROOT_INFO(help_text, argv[0]);
    exit(EXIT_FAILURE);
  }

  int logLevel = BFAM_MAX(BFAM_LL_INFO - verbosity, BFAM_LL_ALWAYS);

  bfam_log_init(rank, stdout, logLevel);
  bfam_signal_handler_set();

  int sc_log_priorities = BFAM_MAX(SC_LP_STATISTICS - verbosity, SC_LP_ALWAYS);
  sc_init(comm, 0, 0, NULL, sc_log_priorities);
  p4est_init(NULL, sc_log_priorities);

  prefs_t *prefs = new_prefs(argv[1]);
  print_prefs(prefs);

  run(comm, prefs);

  free_prefs(prefs);
  bfam_free(prefs);

  sc_finalize();
  BFAM_MPI_CHECK(MPI_Finalize());

  bfam_gopt_free(options);

  return EXIT_SUCCESS;
}

#endif

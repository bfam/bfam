#include <bfam_subdomain_dgx.h>
#include <bfam_jacobi.h>
#include <bfam_kron.h>
#include <bfam_log.h>
#include <bfam_util.h>
#include <bfam_vtk.h>

#ifndef BFAM_DGX_DIMENSION
#define BFAM_DGX_DIMENSION
#define USE_GENERIC_DGX_DIMENSION
#else
#define DIM (BFAM_DGX_DIMENSION)
#endif

#define bfam_subdomain_dgx_get_interpolator                                    \
  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_get_interpolator_, BFAM_DGX_DIMENSION)

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

int BFAM_APPEND_EXPAND(bfam_subdomain_dgx_clear_dgx_ops_dict_,
                       BFAM_DGX_DIMENSION)(const char *key, void *val,
                                           void *args)
{
  bfam_free_aligned(val);
  return 1;
}

int BFAM_APPEND_EXPAND(bfam_subdomain_dgx_clear_interpolation_dict_,
                       BFAM_DGX_DIMENSION)(const char *key, void *val,
                                           void *args)
{
  bfam_subdomain_dgx_interpolator_t *interp =
      (bfam_subdomain_dgx_interpolator_t *)val;
  for (bfam_locidx_t k = 0; k < interp->num_prj; k++)
    if (interp->prj[k])
      bfam_free_aligned(interp->prj[k]);
  bfam_free(interp->prj);
  for (bfam_locidx_t k = 0; k < interp->num_prj; k++)
    if (interp->mass_prj[k])
      bfam_free_aligned(interp->mass_prj[k]);
  bfam_free(interp->mass_prj);
  for (bfam_locidx_t k = 0; k < interp->num_prj; k++)
    if (interp->wi_mass_prj[k])
      bfam_free_aligned(interp->wi_mass_prj[k]);
  bfam_free(interp->wi_mass_prj);
  bfam_free(val);
  return 1;
}

static void init_interpolator(bfam_subdomain_dgx_interpolator_t *interp_a2b,
                              const int N_a, const int N_b)
{
  const bfam_locidx_t num_prj = 5;
  const bfam_locidx_t Np_a = N_a + 1;
  const bfam_locidx_t Np_b = N_b + 1;

  interp_a2b->N_src = N_a;
  interp_a2b->N_dst = N_b;
  interp_a2b->num_prj = num_prj;
  interp_a2b->prj = bfam_malloc(num_prj * sizeof(bfam_real_t *));

  /* Storage for the projection operators */
  for (bfam_locidx_t k = 0; k < num_prj; k++)
  {
    if (k == 0 && N_a == N_b)
      interp_a2b->prj[k] = NULL;
    else
      interp_a2b->prj[k] =
          bfam_malloc_aligned(Np_a * Np_b * sizeof(bfam_real_t));
  }

  interp_a2b->mass_prj = bfam_malloc(num_prj * sizeof(bfam_real_t *));

  interp_a2b->wi_mass_prj = bfam_malloc(num_prj * sizeof(bfam_real_t *));

  /* Storage for the projection operators */
  for (bfam_locidx_t k = 0; k < num_prj; k++)
    interp_a2b->mass_prj[k] =
        bfam_malloc_aligned(Np_a * Np_b * sizeof(bfam_real_t));
  for (bfam_locidx_t k = 0; k < num_prj; k++)
    interp_a2b->wi_mass_prj[k] =
        bfam_malloc_aligned(Np_a * Np_b * sizeof(bfam_real_t));
}

static void fill_grid_data(bfam_locidx_t N, bfam_long_real_t *lr,
                           bfam_long_real_t *lw, bfam_long_real_t *lV,
                           bfam_long_real_t *M)
{
  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);
  bfam_jacobi_p_vandermonde(0, 0, N, N + 1, lr, lV);
  bfam_jacobi_p_mass(0, 0, N, lV, M);
}

/* Build interpolation and projection between space a and g. Space a is the
 * lower order space and space g is the higher order space.
 */
static void fill_interp_proj_data(bfam_locidx_t N_a, bfam_locidx_t N_g,
                                  bfam_long_real_t *V_a, bfam_long_real_t *M_a,
                                  bfam_long_real_t *M_g, bfam_long_real_t *lr_g,
                                  bfam_long_real_t *I_a2g,
                                  bfam_long_real_t *P_g2a)
{
  BFAM_ASSERT(N_g >= N_a);
  bfam_locidx_t Np_a = N_a + 1;
  bfam_locidx_t Np_g = N_g + 1;

  bfam_jacobi_p_interpolation(0, 0, N_a, Np_g, lr_g, V_a, I_a2g);

  bfam_long_real_t MP_g2a[Np_g * Np_a];
  bfam_long_real_t Vt_MP_g2a[Np_g * Np_a];
  for (bfam_locidx_t n = 0; n < Np_g * Np_a; n++)
  {
    MP_g2a[n] = 0;
    Vt_MP_g2a[n] = 0;
    P_g2a[n] = 0;
  }

  /* MP_g2a = M_a * P_g2a = I_a2g^{T} * M_g */
  /* [Np_a X Np_g] = [Np_a X Np_g] [Np_g X Np_g] */
  bfam_util_mTmmult(Np_a, Np_g, Np_g, I_a2g, Np_g, M_g, Np_g, MP_g2a, Np_a);

  /* Vt_MP_g2a = V_a^{-1} * P_g2a = V_a^{T} MP_g2a */
  /* [Np_a X Np_g] = [Np_a X Np_a] [Np_a X Np_g] */
  bfam_util_mTmmult(Np_a, Np_g, Np_a, V_a, Np_a, MP_g2a, Np_a, Vt_MP_g2a, Np_a);

  /* P_g2a = V_a * V_a^{T} I_a2g^{T} * M_g */
  /* [Np_a X Np_g] = [Np_a X Np_a] [Np_a X Np_g] */
  bfam_util_mmmult(Np_a, Np_g, Np_a, V_a, Np_a, Vt_MP_g2a, Np_a, P_g2a, Np_a);

  /*
  if (N_g != N_a)
  {
    BFAM_INFO("%d -> %d", (int)N_g, (int)N_a);
    for (int i = 0; i < Np_a; i++)
      for (int j = 0; j < Np_g; j++)
        BFAM_INFO("P_g2a[%d][%d] = %+" BFAM_REAL_PRIe, i, j,
                  (bfam_real_t)P_g2a[i * Np_g + j]);
  }
  */
}

static void fill_hanging_data(bfam_locidx_t N, bfam_locidx_t Np,
                              bfam_long_real_t *lr_f, bfam_long_real_t *M,
                              bfam_long_real_t *lV, bfam_long_real_t *P_b2f,
                              bfam_long_real_t *P_t2f, bfam_long_real_t *I_f2b,
                              bfam_long_real_t *I_f2t)
{
  const bfam_long_real_t HALF = BFAM_LONG_REAL(0.5);

  /* Grid for the top and bottom */
  bfam_long_real_t lr_t[Np];
  bfam_long_real_t lr_b[Np];

  for (bfam_locidx_t k = 0; k < Np; k++)
  {
    lr_t[k] = HALF * (lr_f[k] + 1);
    lr_b[k] = HALF * (lr_f[k] - 1);
  }

  /* Mass matrix for the top and bottom */
  bfam_long_real_t Mh[Np * Np];
  for (bfam_locidx_t k = 0; k < Np * Np; k++)
    Mh[k] = HALF * M[k];

  fill_interp_proj_data(N, N, lV, M, Mh, lr_t, I_f2t, P_t2f);
  fill_interp_proj_data(N, N, lV, M, Mh, lr_b, I_f2b, P_b2f);
}

static void multiply_projections(const int N_b, const int N_a, const int N_g,
                                 bfam_long_real_t *P_g2b,
                                 bfam_long_real_t *P_g2g,
                                 bfam_long_real_t *P_a2g, bfam_real_t *P_a2b,
                                 bfam_long_real_t *M_b, bfam_real_t *MP_a2b,
                                 bfam_long_real_t *w_b, bfam_real_t *wiMP_a2b)
{
  const int Np_b = N_b + 1;
  const int Np_g = N_g + 1;
  const int Np_a = N_a + 1;

  BFAM_ASSERT(MP_a2b);

  /* In the case the target is NULL return */
  if (!P_a2b)
  {
    BFAM_ASSERT(Np_a == Np_b);
    for (bfam_locidx_t n = 0; n < Np_a * Np_b; n++)
      MP_a2b[n] = (bfam_real_t)(M_b[n]);
    for (bfam_locidx_t j = 0; j < Np_a; j++)
      for (bfam_locidx_t i = 0; i < Np_b; i++)
        wiMP_a2b[i + j * Np_b] = (bfam_real_t)(M_b[i + j * Np_b] / w_b[i]);
    return;
  }

  /* First set up the long storage for multiplication */

  bfam_long_real_t tmp[Np_b * Np_a];
  bfam_long_real_t *l_P = tmp;
  bfam_long_real_t l_MP[Np_b * Np_a];

  for (bfam_locidx_t n = 0; n < Np_a * Np_b; n++)
  {
    tmp[n] = 0;
    l_MP[n] = 0;
  }
  /* No glue to b (thus b is glue) */
  if (!P_g2b && P_g2g && P_a2g)
    bfam_util_mmmult(Np_b, Np_a, Np_g, P_g2g, Np_b, P_a2g, Np_g, l_P, Np_b);
  else if (P_g2b && P_g2g && !P_a2g)
    bfam_util_mmmult(Np_b, Np_a, Np_g, P_g2b, Np_b, P_g2g, Np_g, l_P, Np_b);
  else if (!P_g2b && !P_g2g && P_a2g)
    l_P = P_a2g;
  else if (!P_g2b && P_g2g && !P_a2g)
    l_P = P_g2g;
  else if (P_g2b && !P_g2g && !P_a2g)
    l_P = P_g2b;
  else
    BFAM_ABORT("Case of all NULL or all not NULL is not handled");

  /* MP_a2b = M_b * P_a2b */
  /* [Np_b X Np_a] = [Np_b X Np_b] [Np_b X Np_a]  */
  bfam_util_mmmult(Np_b, Np_a, Np_b, M_b, Np_b, l_P, Np_b, l_MP, Np_b);
  for (bfam_locidx_t n = 0; n < Np_a * Np_b; n++)
  {
    P_a2b[n] = (bfam_real_t)l_P[n];
    MP_a2b[n] = (bfam_real_t)(l_MP[n]);
  }
  for (bfam_locidx_t j = 0; j < Np_a; j++)
    for (bfam_locidx_t i = 0; i < Np_b; i++)
      wiMP_a2b[i + j * Np_b] = (bfam_real_t)(l_MP[i + j * Np_b] / w_b[i]);
}

static void create_interpolators(bfam_subdomain_dgx_interpolator_t *interp_a2b,
                                 bfam_subdomain_dgx_interpolator_t *interp_b2a,
                                 const int N_a, const int N_b)
{
  /* Figure out the glue space and the number of points */
  const int N_g = BFAM_MAX(N_a, N_b);
  const int Np_g = N_g + 1;
  const bfam_locidx_t Np_a = N_a + 1;
  const bfam_locidx_t Np_b = N_b + 1;

  /* initialize the interpolators */
  init_interpolator(interp_a2b, N_a, N_b);
  if (interp_b2a)
    init_interpolator(interp_b2a, N_b, N_a);

  /* Set up the reference grids for the side a */
  bfam_long_real_t lr_a[Np_a];
  bfam_long_real_t lw_a[Np_a];
  bfam_long_real_t V_a[Np_a * Np_a];
  bfam_long_real_t M_a[Np_a * Np_a];
  fill_grid_data(N_a, lr_a, lw_a, V_a, M_a);

  /* Set up the reference grids for the side b */
  bfam_long_real_t lr_b[Np_b];
  bfam_long_real_t lw_b[Np_b];
  bfam_long_real_t V_b[Np_b * Np_b];
  bfam_long_real_t M_b[Np_b * Np_b];
  fill_grid_data(N_b, lr_b, lw_b, V_b, M_b);

  /* Set up the reference grids for the glue space */
  bfam_long_real_t lr_g[Np_g];
  bfam_long_real_t lw_g[Np_g];
  bfam_long_real_t V_g[Np_g * Np_g];
  bfam_long_real_t M_g[Np_g * Np_g];
  fill_grid_data(N_g, lr_g, lw_g, V_g, M_g);

  /* Interpolate to the hanging faces */
  bfam_long_real_t *prj_g[5];
  prj_g[0] = NULL;
  for (bfam_locidx_t k = 1; k < 5; k++)
    prj_g[k] = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_g);

  fill_hanging_data(N_g, Np_g, lr_g, M_g, V_g, prj_g[1], prj_g[2], prj_g[3],
                    prj_g[4]);

  /* Interpolate to the intermediate space and projection back */
  bfam_long_real_t *I_a2g = NULL;
  bfam_long_real_t *P_g2a = NULL;
  bfam_long_real_t *I_b2g = NULL;
  bfam_long_real_t *P_g2b = NULL;
  if (N_a < N_g && N_b == N_g)
  {
    I_a2g = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_a);
    P_g2a = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_a);
    fill_interp_proj_data(N_a, N_g, V_a, M_a, M_g, lr_g, I_a2g, P_g2a);
  }
  if (N_a == N_g && N_b < N_g)
  {
    I_b2g = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_b);
    P_g2b = bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_g * Np_b);
    fill_interp_proj_data(N_b, N_g, V_b, M_b, M_g, lr_g, I_b2g, P_g2b);
  }

  for (bfam_locidx_t k = 0; k < 5; k++)
    multiply_projections(N_b, N_a, N_g, P_g2b, prj_g[k], I_a2g,
                         interp_a2b->prj[k], M_b, interp_a2b->mass_prj[k], lw_b,
                         interp_a2b->wi_mass_prj[k]);
  if (interp_b2a)
    for (bfam_locidx_t k = 0; k < 5; k++)
      multiply_projections(N_a, N_b, N_g, P_g2a, prj_g[k], I_b2g,
                           interp_b2a->prj[k], M_a, interp_b2a->mass_prj[k],
                           lw_a, interp_b2a->wi_mass_prj[k]);

  if (I_a2g)
    bfam_free_aligned(I_a2g);
  if (P_g2a)
    bfam_free_aligned(P_g2a);
  if (I_b2g)
    bfam_free_aligned(I_b2g);
  if (P_g2b)
    bfam_free_aligned(P_g2b);
  bfam_free_aligned(prj_g[1]);
  bfam_free_aligned(prj_g[2]);
  bfam_free_aligned(prj_g[3]);
  bfam_free_aligned(prj_g[4]);
}

bfam_subdomain_dgx_interpolator_t *
    BFAM_APPEND_EXPAND(bfam_subdomain_dgx_get_interpolator_,
                       BFAM_DGX_DIMENSION)(bfam_dictionary_t *N2N,
                                           const int N_src, const int N_dst,
                                           const int inDIM)
{
  char str[BFAM_BUFSIZ];
  snprintf(str, BFAM_BUFSIZ, "%d_to_%d", N_src, N_dst);

  bfam_subdomain_dgx_interpolator_t *interp =
      (bfam_subdomain_dgx_interpolator_t *)bfam_dictionary_get_value_ptr(N2N,
                                                                         str);
  if (!interp)
  {
    interp = bfam_malloc(sizeof(bfam_subdomain_dgx_interpolator_t));
    bfam_subdomain_dgx_interpolator_t *interp2 = NULL;
    if (N_src != N_dst)
      interp2 = bfam_malloc(sizeof(bfam_subdomain_dgx_interpolator_t));
    create_interpolators(interp, interp2, N_src, N_dst);
    BFAM_VERBOSE(">>>>>> Interpolator `%s' created", str);
    int rval = bfam_dictionary_insert_ptr(N2N, str, interp);
    BFAM_ASSERT(rval != 1);

    if (interp2)
    {
      char str2[BFAM_BUFSIZ];
      snprintf(str2, BFAM_BUFSIZ, "%d_to_%d", N_dst, N_src);
      BFAM_VERBOSE(">>>>>> Interpolator `%s' created", str2);
      rval = bfam_dictionary_insert_ptr(N2N, str2, interp2);
      BFAM_ASSERT(rval != 1);
    }
  }

  return interp;
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_interpolate_data_,
                        BFAM_DGX_DIMENSION)(const bfam_real_t *src,
                                            int8_t c_src, int N_src,
                                            bfam_real_t *dst, int8_t c_dst,
                                            int N_dst, uint8_t flags,
                                            bfam_dictionary_t *N2N,
                                            const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_interpolate_data");
  const int DIM = inDIM;
#endif

  BFAM_VERBOSE("flag:%2d", (int)flags);
  BFAM_VERBOSE("  src: c:%2d N:%02d p:%p", (int)c_src, (int)N_src, src);
  BFAM_VERBOSE("  dst: c:%2d N:%02d p:%p", (int)c_dst, (int)N_dst, dst);

  bfam_subdomain_dgx_interpolator_t *interp =
      bfam_subdomain_dgx_get_interpolator(N2N, N_src, N_dst, inDIM);
  BFAM_ASSERT(interp);

  bfam_real_t *Px = NULL;
  bfam_real_t *Py = NULL;
  bfam_real_t *Pz = NULL;

  switch (flags)
  {
  case BFAM_FLAG_COARSEN:
    BFAM_VERBOSE("  coarsen in h");
    Px = interp->prj[(c_src / 1) % 2 + 1];
    if (DIM > 1)
      Py = interp->prj[(c_src / 2) % 2 + 1];
    if (DIM > 2)
      Pz = interp->prj[(c_src / 4) % 2 + 1];
    break;
  case BFAM_FLAG_REFINE:
    BFAM_VERBOSE("  refine in h");
    Px = interp->prj[(c_dst / 1) % 2 + 3];
    if (DIM > 1)
      Py = interp->prj[(c_dst / 2) % 2 + 3];
    if (DIM > 2)
      Pz = interp->prj[(c_dst / 4) % 2 + 3];
    break;
  case BFAM_FLAG_SAME:
    BFAM_VERBOSE("  no change in h");
    Px = interp->prj[0];
    if (DIM > 1)
      Py = interp->prj[0];
    if (DIM > 2)
      Pz = interp->prj[0];
    break;
  default:
    BFAM_ABORT("These are not the droids you are looking for");
  }

#ifdef BFAM_DEBUG
  if (DIM == 2)
    BFAM_ASSERT((Px && Py) || !(Px || Py));
  if (DIM == 3)
    BFAM_ASSERT((Px && Py && Pz) || !(Px || Py || Pz));
#endif

  /* In this case we actually need to the do multiplication */
  if (Px)
  {
    const bfam_locidx_t N1_src = N_src + 1;
    const bfam_locidx_t N1_dst = N_dst + 1;
    if (DIM == 1)
      for (bfam_locidx_t jx = 0; jx < N1_src; jx++)
        for (bfam_locidx_t ix = 0; ix < N1_dst; ix++)
          dst[ix] += Px[jx * N1_dst + ix] * src[jx];
    else if (DIM == 2)
      for (bfam_locidx_t jy = 0; jy < N1_src; jy++)
        for (bfam_locidx_t jx = 0; jx < N1_src; jx++)
          for (bfam_locidx_t iy = 0; iy < N1_dst; iy++)
            for (bfam_locidx_t ix = 0; ix < N1_dst; ix++)
              dst[ix + iy * N1_dst] += Py[jy * N1_dst + iy] *
                                       Px[jx * N1_dst + ix] *
                                       src[jx + jy * N1_src];
    else if (DIM == 3)
      for (bfam_locidx_t jz = 0; jz < N1_src; jz++)
        for (bfam_locidx_t jy = 0; jy < N1_src; jy++)
          for (bfam_locidx_t jx = 0; jx < N1_src; jx++)
            for (bfam_locidx_t iz = 0; iz < N1_dst; iz++)
              for (bfam_locidx_t iy = 0; iy < N1_dst; iy++)
                for (bfam_locidx_t ix = 0; ix < N1_dst; ix++)
                  dst[ix + iy * N1_dst + iz * N1_dst * N1_dst] +=
                      Pz[jz * N1_dst + iz] * Py[jy * N1_dst + iy] *
                      Px[jx * N1_dst + ix] *
                      src[jx + jy * N1_src + jz * N1_src * N1_src];
  }
  /* In this case we just do a memcopy */
  else
  {
    BFAM_ASSERT(N_src == N_dst);
    if (DIM == 1)
      memcpy(dst, src, sizeof(bfam_real_t) * (N_src + 1));
    else if (DIM == 2)
      memcpy(dst, src, sizeof(bfam_real_t) * (N_src + 1) * (N_src + 1));
    else if (DIM == 3)
      memcpy(dst, src,
             sizeof(bfam_real_t) * (N_src + 1) * (N_src + 1) * (N_src + 1));
  }
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_free_,
                        BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_point_interp_t *point)
{
  point->sub = NULL;
  point->elem = -1;

  if (point->file)
    fclose(point->file);
  point->file = NULL;

  strncpy(point->filename, "", BFAM_BUFSIZ);

  for (int n = 0; n < point->num_interp; n++)
    bfam_free_aligned(point->interp[n]);
  bfam_free_aligned(point->interp);
  point->interp = NULL;

  point->num_interp = 0;
}

static bfam_real_t BFAM_APPEND_EXPAND(point_interp_, BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_point_interp_t *point, bfam_real_t *source,
    const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_point_interp_field_m");
  const int DIM = inDIM;
#endif
  const bfam_locidx_t elem = point->elem;

  bfam_real_t d = 0;

  bfam_subdomain_dgx_t *sub = point->sub;

  const bfam_real_t *s = source + elem * sub->Np;
  bfam_real_t **interp = point->interp;

  int Nrp = point->sub->N + 1;
  if (DIM == 1)
    for (int l = 0; l < Nrp; l++)
      d += interp[0][l] * s[l];
  else if (DIM == 2)
    for (int m = 0; m < Nrp; m++)
      for (int l = 0; l < Nrp; l++)
        d += interp[0][l] * interp[1][m] * s[m * Nrp + l];
  else if (DIM == 3)
    for (int n = 0; n < Nrp; n++)
      for (int m = 0; m < Nrp; m++)
        for (int l = 0; l < Nrp; l++)
          d += interp[0][l] * interp[1][m] * interp[2][n] *
               s[(n * Nrp + m) * Nrp + l];
  else
    BFAM_ABORT("Cannot handle dim = %d", DIM);

  return d;
}

bfam_real_t BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_m_,
                               BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_point_interp_t *point, const char *prefix,
    const char *field, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_point_interp_field_m");
  const int DIM = inDIM;
#endif

  bfam_subdomain_dgx_t *sub = point->sub;
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(source, prefix, field,
                                   &sub->base.glue_m->fields);

  return BFAM_APPEND_EXPAND(point_interp_, BFAM_DGX_DIMENSION)(point, source,
                                                               DIM);
}

bfam_real_t BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_,
                               BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_point_interp_t *point, const char *prefix,
    const char *field, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_point_interp_field");
  const int DIM = inDIM;
#endif

  bfam_subdomain_dgx_t *sub = point->sub;
  BFAM_LOAD_FIELD_RESTRICT_ALIGNED(source, prefix, field, &sub->base.fields);
  return BFAM_APPEND_EXPAND(point_interp_, BFAM_DGX_DIMENSION)(point, source,
                                                               DIM);
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_fields_,
                        BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_point_interp_t *point, bfam_real_t t, const char *prefix,
    const char **fields, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_point_interp_fields");
// const int DIM = inDIM;
#endif
  BFAM_ASSERT(point->file);

  fprintf(point->file, "%24.16" BFAM_REAL_PRIe, t);
  for (int n = 0; fields[n]; n++)
  {
    const bfam_real_t v =
        BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_field_,
                           BFAM_DGX_DIMENSION)(point, prefix, fields[n], inDIM);
    fprintf(point->file, " %24.16" BFAM_REAL_PRIe, v);
  }
  fprintf(point->file, "\n");
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_open_,
                        BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_point_interp_t *point)
{
  point->file = fopen(point->filename, "w");
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_init_,
                        BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_point_interp_t *point, bfam_subdomain_dgx_t *sub,
    const bfam_locidx_t elem, const bfam_real_t *r, const char *filename,
    const int filename_size, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_point_interp_init");
  const int DIM = inDIM;
#endif

  point->sub = sub;
  point->elem = elem;
  BFAM_ASSERT(point->elem < sub->K);
  point->num_interp = DIM;
  point->interp = bfam_malloc_aligned(DIM * sizeof(bfam_real_t *));

  for (int n = 0; n < DIM; n++)
    point->interp[n] = bfam_malloc_aligned((sub->N + 1) * sizeof(bfam_real_t));

  strncpy(point->filename, filename, BFAM_BUFSIZ);
  point->file = NULL;

  bfam_long_real_t *cal_interp =
      bfam_malloc_aligned(sizeof(bfam_long_real_t) * (sub->N + 1));
  bfam_long_real_t lr[1] = {0};

  for (int n = 0; n < DIM; n++)
  {
    lr[0] = r[n];
    bfam_jacobi_p_interpolation(0, 0, sub->N, 1, lr, sub->lV, cal_interp);
    for (int m = 0; m < (sub->N + 1); m++)
      point->interp[n][m] = (bfam_real_t)cal_interp[m];
  }

  bfam_free_aligned(cal_interp);
}

bfam_subdomain_dgx_point_interp_t *BFAM_APPEND_EXPAND(
    bfam_subdomain_dgx_point_interp_new_,
    BFAM_DGX_DIMENSION)(bfam_subdomain_dgx_t *sub, const bfam_locidx_t elem,
                        const bfam_real_t *r, const char *filename,
                        const int filename_size, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_point_interp_new");
  const int DIM = inDIM;
#endif

  bfam_subdomain_dgx_point_interp_t *point =
      bfam_malloc(sizeof(bfam_subdomain_dgx_point_interp_t));

  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_point_interp_init_, BFAM_DGX_DIMENSION)(
      point, sub, elem, r, filename, filename_size, DIM);

  return point;
}

typedef struct bfam_subdomain_dgx_get_put_data
{
  bfam_subdomain_dgx_t *sub;
  bfam_real_t *buffer;
  size_t size;
  size_t field;
} bfam_subdomain_dgx_get_put_data_t;

static int bfam_subdomain_dgx_put_scalar_fields_p(const char *key, void *val,
                                                  void *arg)
{
  bfam_subdomain_dgx_get_put_data_t *data =
      (bfam_subdomain_dgx_get_put_data_t *)arg;

  const bfam_locidx_t K = data->sub->K;
  const int Np = data->sub->Np;

  const size_t buffer_offset = data->field * Np * K;
  const bfam_real_t *restrict recv_field = data->buffer + buffer_offset;

  BFAM_ASSERT((data->field + 1) * Np * K * sizeof(bfam_real_t) <= data->size);

  bfam_real_t *restrict glue_field = val;

  BFAM_ASSUME_ALIGNED(glue_field, 32);

  memcpy(glue_field, recv_field, K * Np * sizeof(bfam_real_t));

  ++data->field;
  return 0;
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_add_rates_glue_p_,
                        BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_t *sub, const char *field_prefix_lhs,
    const char *field_prefix_rhs, const char *rate_prefix,
    const bfam_long_real_t a_l, const char **scalars, const char **vectors,
    const char **tensors)
{

  const bfam_real_t a = (bfam_real_t)a_l;
  bfam_locidx_t num_pts = sub->K * sub->Np;

  for (int s = 0; scalars[s] != NULL; s++)
  {
    char key[BFAM_BUFSIZ];

    snprintf(key, BFAM_BUFSIZ, "%s%s", field_prefix_lhs, scalars[s]);
    bfam_real_t *lhs = (bfam_real_t *)bfam_dictionary_get_value_ptr(
        &sub->base.glue_p->fields, key);

    snprintf(key, BFAM_BUFSIZ, "%s%s", field_prefix_rhs, scalars[s]);
    bfam_real_t *rhs = (bfam_real_t *)bfam_dictionary_get_value_ptr(
        &sub->base.glue_p->fields, key);

    snprintf(key, BFAM_BUFSIZ, "%s%s", rate_prefix, scalars[s]);
    bfam_real_t *rate = (bfam_real_t *)bfam_dictionary_get_value_ptr(
        &sub->base.glue_p->fields, key);
    BFAM_ASSERT(lhs && rhs && rate);

    for (int n = 0; n < num_pts; n++)
      lhs[n] = rhs[n] + a * rate[n];
  }
  const char *vt_postfix[] = {"n", "p1", "p2", "p3", NULL};
  for (int v = 0; vectors[v] != NULL; v++)
  {
    for (int c = 0; vt_postfix[c] != NULL; c++)
    {
      char key[BFAM_BUFSIZ];
      snprintf(key, BFAM_BUFSIZ, "%s%s%s", field_prefix_lhs, vectors[v],
               vt_postfix[c]);
      bfam_real_t *lhs = (bfam_real_t *)bfam_dictionary_get_value_ptr(
          &sub->base.glue_p->fields, key);

      snprintf(key, BFAM_BUFSIZ, "%s%s%s", field_prefix_rhs, vectors[v],
               vt_postfix[c]);
      bfam_real_t *rhs = (bfam_real_t *)bfam_dictionary_get_value_ptr(
          &sub->base.glue_p->fields, key);

      snprintf(key, BFAM_BUFSIZ, "%s%s%s", rate_prefix, vectors[v],
               vt_postfix[c]);
      bfam_real_t *rate = (bfam_real_t *)bfam_dictionary_get_value_ptr(
          &sub->base.glue_p->fields, key);

      BFAM_ASSERT(lhs && rhs && rate);
      for (int n = 0; n < num_pts; n++)
        lhs[n] = rhs[n] + a * rate[n];
    }
  }
  for (int t = 0; tensors[t] != NULL; t++)
  {
    for (int c = 0; vt_postfix[c] != NULL; c++)
    {
      char key[BFAM_BUFSIZ];
      snprintf(key, BFAM_BUFSIZ, "%s%s%s", field_prefix_lhs, tensors[t],
               vt_postfix[c]);
      bfam_real_t *lhs = (bfam_real_t *)bfam_dictionary_get_value_ptr(
          &sub->base.glue_p->fields, key);

      snprintf(key, BFAM_BUFSIZ, "%s%s%s", field_prefix_rhs, tensors[t],
               vt_postfix[c]);
      bfam_real_t *rhs = (bfam_real_t *)bfam_dictionary_get_value_ptr(
          &sub->base.glue_p->fields, key);

      snprintf(key, BFAM_BUFSIZ, "%s%s%s", rate_prefix, tensors[t],
               vt_postfix[c]);
      bfam_real_t *rate = (bfam_real_t *)bfam_dictionary_get_value_ptr(
          &sub->base.glue_p->fields, key);

      BFAM_ASSERT(lhs && rhs && rate);
      for (int n = 0; n < num_pts; n++)
        lhs[n] = rhs[n] + a * rate[n];
    }
  }
}

static void bfam_subdomain_dgx_get_recv_buffer(bfam_subdomain_t *thisSubdomain,
                                               void *buffer, size_t recv_sz,
                                               void *comm_args)
{
  bfam_subdomain_dgx_get_put_data_t data;

  data.sub = (bfam_subdomain_dgx_t *)thisSubdomain;
  data.buffer = (bfam_real_t *)buffer;
  data.size = recv_sz;
  data.field = 0;

  if (comm_args == NULL)
  {
    BFAM_ASSERT(recv_sz ==
                data.sub->base.glue_p->fields.num_entries * data.sub->K *
                    data.sub->Np * sizeof(bfam_real_t));

    /*
     * Fill fields_p from the recv buffer.
     */
    bfam_dictionary_allprefixed_ptr(&data.sub->base.glue_p->fields, "",
                                    &bfam_subdomain_dgx_put_scalar_fields_p,
                                    &data);

    BFAM_ASSERT(recv_sz ==
                data.sub->base.glue_p->fields.num_entries * data.sub->K *
                    data.sub->Np * sizeof(bfam_real_t));
  }
  else
  {
    char prefix[BFAM_BUFSIZ];
    prefix[0] = '\0';

    bfam_subdomain_comm_args_t *args = (bfam_subdomain_comm_args_t *)comm_args;
    if (args->user_prefix_function)
      args->user_prefix_function(thisSubdomain, prefix, BFAM_BUFSIZ,
                                 args->user_data);

    for (int s = 0; args->scalars_p[s] != NULL; s++)
    {
      char key[BFAM_BUFSIZ];
      snprintf(key, BFAM_BUFSIZ, "%s%s", prefix, args->scalars_p[s]);
      void *field =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_p->fields, key);
      bfam_subdomain_dgx_put_scalar_fields_p(key, field, &data);
    }
    const char *vt_postfix[] = {"n", "p1", "p2", "p3", NULL};
    for (int v = 0; args->vectors_p[v] != NULL; v++)
    {
      for (int c = 0; vt_postfix[c] != NULL; c++)
      {
        char key[BFAM_BUFSIZ];
        snprintf(key, BFAM_BUFSIZ, "%s%s%s", prefix, args->vectors_p[v],
                 vt_postfix[c]);
        void *field =
            bfam_dictionary_get_value_ptr(&data.sub->base.glue_p->fields, key);
        bfam_subdomain_dgx_put_scalar_fields_p(key, field, &data);
      }
    }
    for (int t = 0; args->tensors_p[t] != NULL; t++)
    {
      for (int c = 0; vt_postfix[c] != NULL; c++)
      {
        char key[BFAM_BUFSIZ];
        snprintf(key, BFAM_BUFSIZ, "%s%s%s", prefix, args->tensors_p[t],
                 vt_postfix[c]);
        void *field =
            bfam_dictionary_get_value_ptr(&data.sub->base.glue_p->fields, key);
        bfam_subdomain_dgx_put_scalar_fields_p(key, field, &data);
      }
    }
    for (int s = 0; args->face_scalars_p[s] != NULL; s++)
    {
      const char *key = args->face_scalars_p[s];
      void *field =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_p->fields, key);
      bfam_subdomain_dgx_put_scalar_fields_p(key, field, &data);
    }
    if (args->user_get_recv_buffer)
    {
      const size_t buffer_offset = data.field * data.sub->Np * data.sub->K;
      args->user_get_recv_buffer(thisSubdomain, data.buffer + buffer_offset,
                                 recv_sz - buffer_offset * sizeof(bfam_real_t),
                                 comm_args);
    }
  }
}

static void bfam_subdomain_dgx_fill_scalar_fields_m(const char *key, void *val,
                                                    void *arg)
{
  bfam_subdomain_dgx_get_put_data_t *data =
      (bfam_subdomain_dgx_get_put_data_t *)arg;

  bfam_subdomain_dgx_t *sub = data->sub;

#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_fill_scalar_fields_m");
  const int DIM = sub->dim;
#endif

  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_p;

  bfam_subdomain_dgx_glue_data_t *glue_m =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_m;

  bfam_subdomain_dgx_t *sub_m = (bfam_subdomain_dgx_t *)glue_m->base.sub_m;

  const bfam_locidx_t K = sub->K;
  const int Np = sub->Np;

  const bfam_locidx_t *restrict EToEp = glue_p->EToEp;
  const bfam_locidx_t *restrict EToEm = glue_p->EToEm;
  const int8_t *restrict EToFm = glue_p->EToFm;
  const int8_t *restrict EToHm = glue_p->EToHm;
  const int8_t *restrict EToOm = glue_p->EToOm;
  BFAM_ASSUME_ALIGNED(EToEp, 32);
  BFAM_ASSUME_ALIGNED(EToEm, 32);
  BFAM_ASSUME_ALIGNED(EToFm, 32);
  BFAM_ASSUME_ALIGNED(EToHm, 32);
  BFAM_ASSUME_ALIGNED(EToOm, 32);

  const int sub_m_Np = sub_m->Np;

  const bfam_real_t *restrict sub_m_field =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, key);

  bfam_real_t *restrict glue_field = val;

  BFAM_ASSERT(sub_m_field != NULL);
  BFAM_ASSERT(glue_field != NULL);

  BFAM_ASSUME_ALIGNED(sub_m_field, 32);
  BFAM_ASSUME_ALIGNED(glue_field, 32);

  BFAM_ASSUME_ALIGNED(EToEp, 32);
  BFAM_ASSUME_ALIGNED(EToEm, 32);
  BFAM_ASSUME_ALIGNED(EToFm, 32);
  BFAM_ASSUME_ALIGNED(EToHm, 32);
  BFAM_ASSUME_ALIGNED(EToOm, 32);

  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    BFAM_ASSERT(EToEp[k] < sub->K);
    BFAM_ASSERT(EToEm[k] < sub_m->K);
    BFAM_ASSERT(EToFm[k] < sub_m->Ng[0]);
    /* BFAM_ASSERT(EToHm[k] < sub_m->Nh); */
    /* BFAM_ASSERT(EToOm[k] < sub_m->No); */

    bfam_locidx_t *restrict fmask = sub_m->gmask[0][EToFm[k]];

    const bfam_real_t *restrict sub_m_elem = sub_m_field + EToEm[k] * sub_m_Np;

    bfam_real_t *restrict glue_elem = glue_field + k * Np;

    /*
     * Interpolate.
     */
    if (EToHm[k] || glue_m->interpolation[0])
    {
      for (int n = 0; n < Np; ++n)
        glue_elem[n] = 0;
      if (DIM == 1)
      {
        /*
         * Decide which interpolation operation to use.
         */
        const bfam_real_t *restrict interpolation =
            glue_m->interpolation[EToHm[k]];
        BFAM_ASSUME_ALIGNED(interpolation, 32);

        /*
         * XXX: Replace with something faster;
         */
        for (int j = 0; j < sub_m->Ngp[0]; ++j)
          for (int i = 0; i < Np; ++i)
            glue_elem[i] += interpolation[j * Np + i] * sub_m_elem[fmask[j]];
      }
      else if (DIM == 2)
      {
        int I1 = (EToHm[k] == 0) ? 0 : (EToHm[k] - 1) / 2 + 1;
        int I2 = (EToHm[k] == 0) ? 0 : (EToHm[k] - 1) % 2 + 1;
        const bfam_real_t *interp1 = glue_m->interpolation[I1];
        const bfam_real_t *interp2 = glue_m->interpolation[I2];
        for (int m = 0; m < sub_m->N + 1; m++)
          for (int l = 0; l < sub_m->N + 1; l++)
            for (int j = 0; j < sub->N + 1; j++)
              for (int i = 0; i < sub->N + 1; i++)
                glue_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] *
                    sub_m_elem[fmask[m * (sub_m->N + 1) + l]];
      }
      else
        BFAM_ABORT("Cannot handle dim = %d", DIM);
    }
    else
    {
      for (int i = 0; i < Np; ++i)
        glue_elem[i] = sub_m_elem[fmask[i]];
    }
  }
}

static void bfam_subdomain_dgx_fill_glue_m(void *val, void *arg)
{
  bfam_subdomain_dgx_get_put_data_t *data =
      (bfam_subdomain_dgx_get_put_data_t *)arg;

  /* if the buffer is NULL we don't do anything */
  if (!data->buffer)
    return;

  bfam_subdomain_dgx_t *sub = data->sub;

#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_fill_glue_m");
#endif

  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_p;

  const bfam_locidx_t K = sub->K;
  const int Np = sub->Np;

  const bfam_locidx_t *restrict EToEp = glue_p->EToEp;
  const int8_t *restrict EToOm = glue_p->EToOm;
  BFAM_ASSUME_ALIGNED(EToEp, 32);
  BFAM_ASSUME_ALIGNED(EToOm, 32);

  const size_t buffer_offset = data->field * Np * K;

  bfam_real_t *restrict send_field = data->buffer + buffer_offset;

  BFAM_ASSERT((data->field + 1) * Np * K * sizeof(bfam_real_t) <= data->size);

  bfam_real_t *restrict glue_field = val;

  BFAM_ASSERT(send_field != NULL);
  BFAM_ASSERT(glue_field != NULL);

  BFAM_ASSUME_ALIGNED(glue_field, 32);

  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    BFAM_ASSERT(EToEp[k] < sub->K);

    bfam_real_t *restrict send_elem = send_field + EToEp[k] * Np;
    bfam_real_t *restrict glue_elem = glue_field + k * Np;

    /*
     * Copy data to send buffer based on orientation.
     */
    BFAM_ASSERT(EToOm[k] >= 0 && EToOm[k] < glue_p->num_orient);
    if (EToOm[k])
      for (int n = 0; n < Np; ++n)
        send_elem[n] = glue_elem[glue_p->mapOm[EToOm[k]][n]];
    else
      memcpy(send_elem, glue_elem, Np * sizeof(bfam_real_t));
  }

  ++data->field;
}

static int bfam_subdomain_dgx_get_scalar_fields_m(const char *key, void *val,
                                                  void *arg)
{
  bfam_subdomain_dgx_fill_scalar_fields_m(key, val, arg);
  bfam_subdomain_dgx_fill_glue_m(val, arg);
  return 0;
}

static void bfam_subdomain_dgx_fill_vector_fields_m(const char *prefix,
                                                    const char **comp, void *vn,
                                                    void *vp1, void *vp2,
                                                    void *vp3, void *arg)
{
  BFAM_ASSUME_ALIGNED(vn, 32);
  BFAM_ASSUME_ALIGNED(vp1, 32);
  BFAM_ASSUME_ALIGNED(vp2, 32);
  BFAM_ASSUME_ALIGNED(vp3, 32);

  bfam_subdomain_dgx_get_put_data_t *data =
      (bfam_subdomain_dgx_get_put_data_t *)arg;

  bfam_subdomain_dgx_t *sub = data->sub;

#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_fill_vector_fields_m");
  const int DIM = sub->dim;
#endif

  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_p;

  bfam_subdomain_dgx_glue_data_t *glue_m =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_m;

  bfam_subdomain_dgx_t *sub_m = (bfam_subdomain_dgx_t *)glue_m->base.sub_m;

  const bfam_locidx_t K = sub->K;
  const int Np = sub->Np;

  const int sub_m_Nfp = sub_m->Ngp[0];
  const int sub_m_Nf = sub_m->Ng[0];
  const int sub_m_Nrp = sub_m->N + 1;

  const bfam_locidx_t *restrict EToEp = glue_p->EToEp;
  const bfam_locidx_t *restrict EToEm = glue_p->EToEm;
  const int8_t *restrict EToFm = glue_p->EToFm;
  const int8_t *restrict EToHm = glue_p->EToHm;
  const int8_t *restrict EToOm = glue_p->EToOm;
  BFAM_ASSUME_ALIGNED(EToEp, 32);
  BFAM_ASSUME_ALIGNED(EToEm, 32);
  BFAM_ASSUME_ALIGNED(EToFm, 32);
  BFAM_ASSUME_ALIGNED(EToHm, 32);
  BFAM_ASSUME_ALIGNED(EToOm, 32);

  const int sub_m_Np = sub_m->Np;

  char str[BFAM_BUFSIZ];
  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[0]);
  const bfam_real_t *v1 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(v1 != NULL);
  BFAM_ASSUME_ALIGNED(v1, 32);

  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[1]);
  const bfam_real_t *v2 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(v2 != NULL);
  BFAM_ASSUME_ALIGNED(v2, 32);

  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[2]);
  const bfam_real_t *v3 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(v3 != NULL);
  BFAM_ASSUME_ALIGNED(v3, 32);

  bfam_real_t *restrict n1 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields_face, "_grid_nx0");
  BFAM_ASSERT(n1 != NULL);
  BFAM_ASSUME_ALIGNED(n1, 32);

  bfam_real_t *restrict n2 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields_face, "_grid_nx1");
  BFAM_ASSERT(n2 != NULL);
  BFAM_ASSUME_ALIGNED(n2, 32);

  bfam_real_t *restrict n3 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields_face, "_grid_nx2");
  if (DIM > 2)
    BFAM_ASSERT(n3 != NULL);
  BFAM_ASSUME_ALIGNED(n3, 32);

  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    BFAM_ASSERT(EToEp[k] < sub->K);
    BFAM_ASSERT(EToEm[k] < sub_m->K);
    BFAM_ASSERT(EToFm[k] < sub_m->Ng[0]);
    /* BFAM_ASSERT(EToHm[k] < sub_m->Nh); */
    /* BFAM_ASSERT(EToOm[k] < sub_m->No); */

    int8_t face = EToFm[k];
    bfam_locidx_t *restrict fmask = sub_m->gmask[0][face];

    const bfam_real_t *v1_m_elem = v1 + EToEm[k] * sub_m_Np;
    const bfam_real_t *v2_m_elem = v2 + EToEm[k] * sub_m_Np;
    const bfam_real_t *v3_m_elem = v3 + EToEm[k] * sub_m_Np;

    bfam_real_t *restrict vn_g_elem = (bfam_real_t *)vn + k * Np;
    bfam_real_t *restrict vp1_g_elem = (bfam_real_t *)vp1 + k * Np;
    bfam_real_t *restrict vp2_g_elem = (bfam_real_t *)vp2 + k * Np;
    bfam_real_t *restrict vp3_g_elem = (bfam_real_t *)vp3 + k * Np;

    /*
     * Interpolate.
     */
    if (EToHm[k] || glue_m->interpolation[0])
    {
      for (int n = 0; n < Np; ++n)
      {
        vn_g_elem[n] = 0;
        vp1_g_elem[n] = 0;
        vp2_g_elem[n] = 0;
        vp3_g_elem[n] = 0;
      }
      if (DIM == 1)
      {
        /*
         * Decide which interpolation operation to use.
         */
        const bfam_real_t *restrict interpolation =
            glue_m->interpolation[EToHm[k]];
        BFAM_ASSUME_ALIGNED(interpolation, 32);
        for (int j = 0; j < sub_m_Nfp; ++j)
        {
          const bfam_locidx_t f = j + sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);

          const bfam_real_t vn_e =
              (n1[f] * v1_m_elem[fmask[j]] + n2[f] * v2_m_elem[fmask[j]]);
          const bfam_real_t vp1_e = v1_m_elem[fmask[j]] - vn_e * n1[f];
          const bfam_real_t vp2_e = v2_m_elem[fmask[j]] - vn_e * n2[f];
          const bfam_real_t vp3_e = v3_m_elem[fmask[j]];

          for (int i = 0; i < Np; ++i)
          {
            vn_g_elem[i] += interpolation[j * Np + i] * vn_e;
            vp1_g_elem[i] += interpolation[j * Np + i] * vp1_e;
            vp2_g_elem[i] += interpolation[j * Np + i] * vp2_e;
            vp3_g_elem[i] += interpolation[j * Np + i] * vp3_e;
          }
        }
      }
      else if (DIM == 2)
      {
        int I1 = (EToHm[k] == 0) ? 0 : (EToHm[k] - 1) / 2 + 1;
        int I2 = (EToHm[k] == 0) ? 0 : (EToHm[k] - 1) % 2 + 1;
        const bfam_real_t *interp1 = glue_m->interpolation[I1];
        const bfam_real_t *interp2 = glue_m->interpolation[I2];
        for (int m = 0; m < sub_m_Nrp; m++)
          for (int l = 0; l < sub_m_Nrp; l++)
          {
            const bfam_locidx_t f =
                l + m * (sub_m_Nrp)+sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);
            const int n = m * (sub_m->N + 1) + l;
            const bfam_real_t vn_e =
                (n1[f] * v1_m_elem[fmask[n]] + n2[f] * v2_m_elem[fmask[n]] +
                 n3[f] * v3_m_elem[fmask[n]]);
            const bfam_real_t vp1_e = v1_m_elem[fmask[n]] - vn_e * n1[f];
            const bfam_real_t vp2_e = v2_m_elem[fmask[n]] - vn_e * n2[f];
            const bfam_real_t vp3_e = v3_m_elem[fmask[n]] - vn_e * n3[f];

            for (int j = 0; j < sub->N + 1; j++)
              for (int i = 0; i < sub->N + 1; i++)
              {
                vn_g_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] * vn_e;
                vp1_g_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] * vp1_e;
                vp2_g_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] * vp2_e;
                vp3_g_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] * vp3_e;
              }
          }
      }
      else
        BFAM_ABORT("Cannot handle dim = %d", DIM);
    }
    else if (DIM == 1)
    {
      for (int j = 0; j < sub_m_Nfp; ++j)
      {
        const bfam_locidx_t f = j + sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);

        vn_g_elem[j] =
            (n1[f] * v1_m_elem[fmask[j]] + n2[f] * v2_m_elem[fmask[j]]);
        vp1_g_elem[j] = v1_m_elem[fmask[j]] - vn_g_elem[j] * n1[f];
        vp2_g_elem[j] = v2_m_elem[fmask[j]] - vn_g_elem[j] * n2[f];
        vp3_g_elem[j] = v3_m_elem[fmask[j]];
      }
    }
    else if (DIM == 2)
    {
      for (int m = 0; m < sub_m_Nrp; m++)
        for (int l = 0; l < sub_m_Nrp; l++)
        {
          const bfam_locidx_t f =
              l + m * (sub_m_Nrp)+sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);

          const int n = m * (sub_m->N + 1) + l;
          vn_g_elem[n] =
              (n1[f] * v1_m_elem[fmask[n]] + n2[f] * v2_m_elem[fmask[n]] +
               n3[f] * v3_m_elem[fmask[n]]);
          vp1_g_elem[n] = v1_m_elem[fmask[n]] - vn_g_elem[n] * n1[f];
          vp2_g_elem[n] = v2_m_elem[fmask[n]] - vn_g_elem[n] * n2[f];
          vp3_g_elem[n] = v3_m_elem[fmask[n]] - vn_g_elem[n] * n3[f];
        }
    }
  }
}

static int bfam_subdomain_dgx_get_vector_fields_m(const char *prefix,
                                                  const char **comp, void *vn,
                                                  void *vp1, void *vp2,
                                                  void *vp3, void *arg)
{
  bfam_subdomain_dgx_fill_vector_fields_m(prefix, comp, vn, vp1, vp2, vp3, arg);
  bfam_subdomain_dgx_fill_glue_m(vn, arg);
  bfam_subdomain_dgx_fill_glue_m(vp1, arg);
  bfam_subdomain_dgx_fill_glue_m(vp2, arg);
  bfam_subdomain_dgx_fill_glue_m(vp3, arg);
  return 0;
}

static void bfam_subdomain_dgx_fill_tensor_fields_m(const char *prefix,
                                                    const char **comp, void *Tn,
                                                    void *Tp1, void *Tp2,
                                                    void *Tp3, void *arg)
{
  BFAM_ASSUME_ALIGNED(Tn, 32);
  BFAM_ASSUME_ALIGNED(Tp1, 32);
  BFAM_ASSUME_ALIGNED(Tp2, 32);
  BFAM_ASSUME_ALIGNED(Tp3, 32);

  bfam_subdomain_dgx_get_put_data_t *data =
      (bfam_subdomain_dgx_get_put_data_t *)arg;

  bfam_subdomain_dgx_t *sub = data->sub;

#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_fill_tensor_fields_m");
  const int DIM = sub->dim;
#endif

  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_p;

  bfam_subdomain_dgx_glue_data_t *glue_m =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_m;

  bfam_subdomain_dgx_t *sub_m = (bfam_subdomain_dgx_t *)glue_m->base.sub_m;

  const bfam_locidx_t K = sub->K;
  const int Np = sub->Np;

  const int sub_m_Nfp = sub_m->Ngp[0];
  const int sub_m_Nf = sub_m->Ng[0];
  const int sub_m_Nrp = sub_m->N + 1;

  const bfam_locidx_t *restrict EToEp = glue_p->EToEp;
  const bfam_locidx_t *restrict EToEm = glue_p->EToEm;
  const int8_t *restrict EToFm = glue_p->EToFm;
  const int8_t *restrict EToHm = glue_p->EToHm;
  const int8_t *restrict EToOm = glue_p->EToOm;
  BFAM_ASSUME_ALIGNED(EToEp, 32);
  BFAM_ASSUME_ALIGNED(EToEm, 32);
  BFAM_ASSUME_ALIGNED(EToFm, 32);
  BFAM_ASSUME_ALIGNED(EToHm, 32);
  BFAM_ASSUME_ALIGNED(EToOm, 32);

  const int sub_m_Np = sub_m->Np;

  char str[BFAM_BUFSIZ];
  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[0]);
  const bfam_real_t *S11 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(S11 != NULL);
  BFAM_ASSUME_ALIGNED(S11, 32);

  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[1]);
  const bfam_real_t *S12 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(S12 != NULL);
  BFAM_ASSUME_ALIGNED(S12, 32);

  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[2]);
  const bfam_real_t *S13 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(S13 != NULL);
  BFAM_ASSUME_ALIGNED(S13, 32);

  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[3]);
  const bfam_real_t *S22 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(S22 != NULL);
  BFAM_ASSUME_ALIGNED(S22, 32);

  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[4]);
  const bfam_real_t *S23 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(S23 != NULL);
  BFAM_ASSUME_ALIGNED(S23, 32);

  snprintf(str, BFAM_BUFSIZ, "%s%s", prefix, comp[5]);
  const bfam_real_t *S33 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields, str);
  BFAM_ASSERT(S33 != NULL);
  BFAM_ASSUME_ALIGNED(S33, 32);

  bfam_real_t *restrict n1 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields_face, "_grid_nx0");
  BFAM_ASSERT(n1 != NULL);
  BFAM_ASSUME_ALIGNED(n1, 32);

  bfam_real_t *restrict n2 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields_face, "_grid_nx1");
  BFAM_ASSERT(n2 != NULL);
  BFAM_ASSUME_ALIGNED(n2, 32);

  bfam_real_t *restrict n3 =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields_face, "_grid_nx2");
  if (DIM > 2)
    BFAM_ASSERT(n3 != NULL);
  BFAM_ASSUME_ALIGNED(n3, 32);

  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    BFAM_ASSERT(EToEp[k] < sub->K);
    BFAM_ASSERT(EToEm[k] < sub_m->K);
    BFAM_ASSERT(EToFm[k] < sub_m->Ng[0]);
    /* BFAM_ASSERT(EToHm[k] < sub_m->Nh); */
    /* BFAM_ASSERT(EToOm[k] < sub_m->No); */

    int8_t face = EToFm[k];
    bfam_locidx_t *restrict fmask = sub_m->gmask[0][face];

    const bfam_real_t *S11_m_elem = S11 + EToEm[k] * sub_m_Np;
    const bfam_real_t *S12_m_elem = S12 + EToEm[k] * sub_m_Np;
    const bfam_real_t *S13_m_elem = S13 + EToEm[k] * sub_m_Np;
    const bfam_real_t *S22_m_elem = S22 + EToEm[k] * sub_m_Np;
    const bfam_real_t *S23_m_elem = S23 + EToEm[k] * sub_m_Np;
    const bfam_real_t *S33_m_elem = S33 + EToEm[k] * sub_m_Np;

    bfam_real_t *restrict Tn_g_elem = (bfam_real_t *)Tn + k * Np;
    bfam_real_t *restrict Tp1_g_elem = (bfam_real_t *)Tp1 + k * Np;
    bfam_real_t *restrict Tp2_g_elem = (bfam_real_t *)Tp2 + k * Np;
    bfam_real_t *restrict Tp3_g_elem = (bfam_real_t *)Tp3 + k * Np;

    /*
     * Interpolate.
     */
    if (EToHm[k] || glue_m->interpolation[0])
    {
      for (int n = 0; n < Np; ++n)
      {
        Tn_g_elem[n] = 0;
        Tp1_g_elem[n] = 0;
        Tp2_g_elem[n] = 0;
        Tp3_g_elem[n] = 0;
      }
      if (DIM == 1)
      {
        /*
         * Decide which interpolation operation to use.
         */
        const bfam_real_t *restrict interpolation =
            glue_m->interpolation[EToHm[k]];
        BFAM_ASSUME_ALIGNED(interpolation, 32);
        for (int j = 0; j < sub_m_Nfp; ++j)
        {
          const bfam_locidx_t f = j + sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);
          const bfam_locidx_t fm = fmask[j];

          bfam_real_t Tp1_e = S11_m_elem[fm] * n1[f] + S12_m_elem[fm] * n2[f];
          bfam_real_t Tp2_e = S12_m_elem[fm] * n1[f] + S22_m_elem[fm] * n2[f];
          bfam_real_t Tp3_e = S13_m_elem[fm] * n1[f] + S23_m_elem[fm] * n2[f];
          bfam_real_t Tn_e = n1[f] * Tp1_e + n2[f] * Tp2_e;

          Tp1_e -= n1[f] * Tn_e;
          Tp2_e -= n2[f] * Tn_e;

          for (int i = 0; i < Np; ++i)
          {
            Tn_g_elem[i] += interpolation[j * Np + i] * Tn_e;
            Tp1_g_elem[i] += interpolation[j * Np + i] * Tp1_e;
            Tp2_g_elem[i] += interpolation[j * Np + i] * Tp2_e;
            Tp3_g_elem[i] += interpolation[j * Np + i] * Tp3_e;
          }
        }
      }
      else if (DIM == 2)
      {
        int I1 = (EToHm[k] == 0) ? 0 : (EToHm[k] - 1) / 2 + 1;
        int I2 = (EToHm[k] == 0) ? 0 : (EToHm[k] - 1) % 2 + 1;
        const bfam_real_t *interp1 = glue_m->interpolation[I1];
        const bfam_real_t *interp2 = glue_m->interpolation[I2];
        for (int m = 0; m < sub_m_Nrp; m++)
          for (int l = 0; l < sub_m_Nrp; l++)
          {
            const bfam_locidx_t f =
                l + m * (sub_m_Nrp)+sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);
            const int n = m * (sub_m->N + 1) + l;
            const bfam_locidx_t fm = fmask[n];

            bfam_real_t Tp1_e = S11_m_elem[fm] * n1[f] +
                                S12_m_elem[fm] * n2[f] + S13_m_elem[fm] * n3[f];
            bfam_real_t Tp2_e = S12_m_elem[fm] * n1[f] +
                                S22_m_elem[fm] * n2[f] + S23_m_elem[fm] * n3[f];
            bfam_real_t Tp3_e = S13_m_elem[fm] * n1[f] +
                                S23_m_elem[fm] * n2[f] + S33_m_elem[fm] * n3[f];
            bfam_real_t Tn_e = n1[f] * Tp1_e + n2[f] * Tp2_e + n3[f] * Tp3_e;

            Tp1_e -= n1[f] * Tn_e;
            Tp2_e -= n2[f] * Tn_e;
            Tp3_e -= n3[f] * Tn_e;

            for (int j = 0; j < sub->N + 1; j++)
              for (int i = 0; i < sub->N + 1; i++)
              {
                Tn_g_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] * Tn_e;
                Tp1_g_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] * Tp1_e;
                Tp2_g_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] * Tp2_e;
                Tp3_g_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] * Tp3_e;
              }
          }
      }
      else
        BFAM_ABORT("Cannot handle dim = %d", DIM);
    }
    else if (DIM == 1)
    {
      for (int j = 0; j < sub_m_Nfp; ++j)
      {
        const bfam_locidx_t f = j + sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);
        const bfam_locidx_t fm = fmask[j];

        Tp1_g_elem[j] = S11_m_elem[fm] * n1[f] + S12_m_elem[fm] * n2[f];
        Tp2_g_elem[j] = S12_m_elem[fm] * n1[f] + S22_m_elem[fm] * n2[f];
        Tp3_g_elem[j] = S13_m_elem[fm] * n1[f] + S23_m_elem[fm] * n2[f];
        Tn_g_elem[j] = n1[f] * Tp1_g_elem[j] + n2[f] * Tp2_g_elem[j];

        Tp1_g_elem[j] -= n1[f] * Tn_g_elem[j];
        Tp2_g_elem[j] -= n2[f] * Tn_g_elem[j];
      }
    }
    else if (DIM == 2)
    {
      for (int m = 0; m < sub_m_Nrp; m++)
        for (int l = 0; l < sub_m_Nrp; l++)
        {
          const bfam_locidx_t f =
              l + m * (sub_m_Nrp)+sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);
          const int n = m * (sub_m->N + 1) + l;
          const bfam_locidx_t fm = fmask[n];

          Tp1_g_elem[n] = S11_m_elem[fm] * n1[f] + S12_m_elem[fm] * n2[f] +
                          S13_m_elem[fm] * n3[f];
          Tp2_g_elem[n] = S12_m_elem[fm] * n1[f] + S22_m_elem[fm] * n2[f] +
                          S23_m_elem[fm] * n3[f];
          Tp3_g_elem[n] = S13_m_elem[fm] * n1[f] + S23_m_elem[fm] * n2[f] +
                          S33_m_elem[fm] * n3[f];
          Tn_g_elem[n] = n1[f] * Tp1_g_elem[n] + n2[f] * Tp2_g_elem[n] +
                         n3[f] * Tp3_g_elem[n];

          Tp1_g_elem[n] -= n1[f] * Tn_g_elem[n];
          Tp2_g_elem[n] -= n2[f] * Tn_g_elem[n];
          Tp3_g_elem[n] -= n3[f] * Tn_g_elem[n];
        }
    }
  }
}

static int bfam_subdomain_dgx_get_tensor_fields_m(const char *prefix,
                                                  const char **comp, void *Tn,
                                                  void *Tp1, void *Tp2,
                                                  void *Tp3, void *arg)
{
  bfam_subdomain_dgx_fill_tensor_fields_m(prefix, comp, Tn, Tp1, Tp2, Tp3, arg);
  bfam_subdomain_dgx_fill_glue_m(Tn, arg);
  bfam_subdomain_dgx_fill_glue_m(Tp1, arg);
  bfam_subdomain_dgx_fill_glue_m(Tp2, arg);
  bfam_subdomain_dgx_fill_glue_m(Tp3, arg);
  return 0;
}

static void bfam_subdomain_dgx_fill_face_scalar_fields_m(const char *key,
                                                         void *val, void *arg)
{
  bfam_subdomain_dgx_get_put_data_t *data =
      (bfam_subdomain_dgx_get_put_data_t *)arg;

  bfam_subdomain_dgx_t *sub = data->sub;

#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_fill_face_scalar_fields_m");
  const int DIM = sub->dim;
#endif

  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_p;

  bfam_subdomain_dgx_glue_data_t *glue_m =
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_m;

  bfam_subdomain_dgx_t *sub_m = (bfam_subdomain_dgx_t *)glue_m->base.sub_m;

  const bfam_locidx_t K = sub->K;
  const int Np = sub->Np;
  const int sub_m_Nf = sub_m->Ng[0];
  const int sub_m_Nfp = sub_m->Ngp[0];

  const bfam_locidx_t *restrict EToEp = glue_p->EToEp;
  const bfam_locidx_t *restrict EToEm = glue_p->EToEm;
  const int8_t *restrict EToFm = glue_p->EToFm;
  const int8_t *restrict EToHm = glue_p->EToHm;
  const int8_t *restrict EToOm = glue_p->EToOm;
  BFAM_ASSUME_ALIGNED(EToEp, 32);
  BFAM_ASSUME_ALIGNED(EToEm, 32);
  BFAM_ASSUME_ALIGNED(EToFm, 32);
  BFAM_ASSUME_ALIGNED(EToHm, 32);
  BFAM_ASSUME_ALIGNED(EToOm, 32);

  const bfam_real_t *restrict sub_m_face_field =
      bfam_dictionary_get_value_ptr(&sub_m->base.fields_face, key);

  bfam_real_t *restrict glue_field = val;

  BFAM_ASSERT(sub_m_face_field != NULL);
  BFAM_ASSERT(glue_field != NULL);

  BFAM_ASSUME_ALIGNED(sub_m_face_field, 32);
  BFAM_ASSUME_ALIGNED(glue_field, 32);

  BFAM_ASSUME_ALIGNED(EToEp, 32);
  BFAM_ASSUME_ALIGNED(EToEm, 32);
  BFAM_ASSUME_ALIGNED(EToFm, 32);
  BFAM_ASSUME_ALIGNED(EToHm, 32);
  BFAM_ASSUME_ALIGNED(EToOm, 32);

  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    BFAM_ASSERT(EToEp[k] < sub->K);
    BFAM_ASSERT(EToEm[k] < sub_m->K);
    BFAM_ASSERT(EToFm[k] < sub_m->Ng[0]);
    /* BFAM_ASSERT(EToHm[k] < sub_m->Nh); */
    /* BFAM_ASSERT(EToOm[k] < sub_m->No); */

    int8_t face = EToFm[k];
    const bfam_real_t *restrict sub_m_face_elem =
        sub_m_face_field + sub_m_Nfp * (face + sub_m_Nf * EToEm[k]);

    bfam_real_t *restrict glue_elem = glue_field + k * Np;

    /*
     * Interpolate.
     */
    if (EToHm[k] || glue_m->interpolation[0])
    {
      for (int n = 0; n < Np; ++n)
        glue_elem[n] = 0;
      if (DIM == 1)
      {
        /*
         * Decide which interpolation operation to use.
         */
        const bfam_real_t *restrict interpolation =
            glue_m->interpolation[EToHm[k]];
        BFAM_ASSUME_ALIGNED(interpolation, 32);

        /*
         * XXX: Replace with something faster;
         */
        for (int j = 0; j < sub_m->Ngp[0]; ++j)
          for (int i = 0; i < Np; ++i)
            glue_elem[i] += interpolation[j * Np + i] * sub_m_face_elem[j];
      }
      else if (DIM == 2)
      {
        int I1 = (EToHm[k] == 0) ? 0 : (EToHm[k] - 1) / 2 + 1;
        int I2 = (EToHm[k] == 0) ? 0 : (EToHm[k] - 1) % 2 + 1;
        const bfam_real_t *interp1 = glue_m->interpolation[I1];
        const bfam_real_t *interp2 = glue_m->interpolation[I2];
        for (int m = 0; m < sub_m->N + 1; m++)
          for (int l = 0; l < sub_m->N + 1; l++)
            for (int j = 0; j < sub->N + 1; j++)
              for (int i = 0; i < sub->N + 1; i++)
                glue_elem[j * (sub->N + 1) + i] +=
                    interp1[(sub->N + 1) * m + j] *
                    interp2[(sub->N + 1) * l + i] *
                    sub_m_face_elem[m * (sub_m->N + 1) + l];
      }
      else
        BFAM_ABORT("Cannot handle dim = %d", DIM);
    }
    else
    {
      for (int i = 0; i < Np; ++i)
        glue_elem[i] = sub_m_face_elem[i];
    }
  }
}

static int bfam_subdomain_dgx_get_face_scalar_fields_m(const char *key,
                                                       void *val, void *arg)
{
  bfam_subdomain_dgx_fill_face_scalar_fields_m(key, val, arg);
  bfam_subdomain_dgx_fill_glue_m(val, arg);
  return 0;
}

static void bfam_subdomain_dgx_put_send_buffer(bfam_subdomain_t *thisSubdomain,
                                               void *buffer, size_t send_sz,
                                               void *comm_args)
{
  bfam_subdomain_dgx_get_put_data_t data;

  data.sub = (bfam_subdomain_dgx_t *)thisSubdomain;
  data.buffer = (bfam_real_t *)buffer;
  data.size = send_sz;
  data.field = 0;

  if (comm_args == NULL)
  {
    BFAM_ASSERT(!buffer ||
                send_sz ==
                    data.sub->base.glue_m->fields.num_entries * data.sub->K *
                        data.sub->Np * sizeof(bfam_real_t));

    /*
     * Fill fields_m and the send buffer from sub_m.
     */
    bfam_dictionary_allprefixed_ptr(&data.sub->base.glue_m->fields, "",
                                    &bfam_subdomain_dgx_get_scalar_fields_m,
                                    &data);
  }
  else
  {
    char prefix[BFAM_BUFSIZ];
    prefix[0] = '\0';

    bfam_subdomain_comm_args_t *args = (bfam_subdomain_comm_args_t *)comm_args;
    if (args->user_prefix_function)
      args->user_prefix_function(thisSubdomain, prefix, BFAM_BUFSIZ,
                                 args->user_data);

    for (int s = 0; args->scalars_m[s] != NULL; s++)
    {
      char key[BFAM_BUFSIZ];
      snprintf(key, BFAM_BUFSIZ, "%s%s", prefix, args->scalars_m[s]);
      void *field =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, key);
      bfam_subdomain_dgx_get_scalar_fields_m(key, field, &data);
    }
    for (int v = 0; args->vectors_m[v] != NULL; v++)
    {
      const char *vec_prefix = args->vectors_m[v];
      char str[BFAM_BUFSIZ];
      snprintf(str, BFAM_BUFSIZ, "%s%sn", prefix, vec_prefix);
      void *vn =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, str);
      snprintf(str, BFAM_BUFSIZ, "%s%sp1", prefix, vec_prefix);
      void *vp1 =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, str);
      snprintf(str, BFAM_BUFSIZ, "%s%sp2", prefix, vec_prefix);
      void *vp2 =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, str);
      snprintf(str, BFAM_BUFSIZ, "%s%sp3", prefix, vec_prefix);
      void *vp3 =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, str);
      BFAM_ASSERT(vn != NULL && vp1 != NULL && vp2 != NULL && vp3 != NULL);
      const char **comps = args->vector_components_m + 3 * v;
      bfam_subdomain_dgx_get_vector_fields_m(prefix, comps, vn, vp1, vp2, vp3,
                                             &data);
    }
    for (int t = 0; args->tensors_m[t] != NULL; t++)
    {
      const char *ten_prefix = args->tensors_m[t];
      char str[BFAM_BUFSIZ];
      snprintf(str, BFAM_BUFSIZ, "%s%sn", prefix, ten_prefix);
      void *Tn =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, str);
      snprintf(str, BFAM_BUFSIZ, "%s%sp1", prefix, ten_prefix);
      void *Tp1 =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, str);
      snprintf(str, BFAM_BUFSIZ, "%s%sp2", prefix, ten_prefix);
      void *Tp2 =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, str);
      snprintf(str, BFAM_BUFSIZ, "%s%sp3", prefix, ten_prefix);
      void *Tp3 =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, str);
      BFAM_ASSERT(Tn != NULL && Tp1 != NULL && Tp2 != NULL && Tp3 != NULL);
      const char **comps = args->tensor_components_m + 3 * t;
      bfam_subdomain_dgx_get_tensor_fields_m(prefix, comps, Tn, Tp1, Tp2, Tp3,
                                             &data);
    }
    for (int s = 0; args->face_scalars_m[s] != NULL; s++)
    {
      const char *key = args->face_scalars_m[s];
      void *field =
          bfam_dictionary_get_value_ptr(&data.sub->base.glue_m->fields, key);
      bfam_subdomain_dgx_get_face_scalar_fields_m(key, field, &data);
    }

    if (args->user_put_send_buffer)
    {
      const size_t buffer_offset = data.field * data.sub->Np * data.sub->K;
      args->user_put_send_buffer(thisSubdomain, data.buffer + buffer_offset,
                                 send_sz - buffer_offset * sizeof(bfam_real_t),
                                 comm_args);
    }
  }
}

static void bfam_subdomain_dgx_comm_info(bfam_subdomain_t *thisSubdomain,
                                         int *rank, bfam_locidx_t *sort,
                                         int num_sort, size_t *send_sz,
                                         size_t *recv_sz, void *comm_args)
{
  BFAM_ASSERT(num_sort > 2);
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;

#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_comm_info");
/* const int DIM = sub->dim; */
#endif

  BFAM_ASSERT(sub->base.glue_m && sub->base.glue_p);

  *rank = sub->base.glue_p->rank;
  sort[0] = sub->base.glue_p->id; /* neighbor ID */
  sort[1] = sub->base.glue_m->id; /* my ID */
  sort[2] = sub->base.uid;        /* my user ID */

  size_t send_num = sub->base.glue_m->fields.num_entries * sub->K * sub->Np;
  size_t recv_num = sub->base.glue_p->fields.num_entries * sub->K * sub->Np;
  BFAM_ASSERT(send_num == recv_num);

  if (comm_args != NULL)
  {
    bfam_subdomain_comm_args_t *args = (bfam_subdomain_comm_args_t *)comm_args;

    int count = 0;
    for (int i = 0; args->scalars_m[i] != NULL; i++)
      count++;
    for (int i = 0; args->vectors_m[i] != NULL; i++)
      count += 4;
    for (int i = 0; args->tensors_m[i] != NULL; i++)
      count += 4;
    for (int i = 0; args->face_scalars_m[i] != NULL; i++)
      count++;
    send_num = count * sub->K * sub->Np;

    count = 0;
    for (int i = 0; args->scalars_p[i] != NULL; i++)
      count++;
    for (int i = 0; args->vectors_p[i] != NULL; i++)
      count += 4;
    for (int i = 0; args->tensors_p[i] != NULL; i++)
      count += 4;
    for (int i = 0; args->face_scalars_p[i] != NULL; i++)
      count++;
    recv_num = count * sub->K * sub->Np;
  }

  *send_sz = send_num * sizeof(bfam_real_t);
  *recv_sz = recv_num * sizeof(bfam_real_t);

  if (comm_args != NULL)
  {
    bfam_subdomain_comm_args_t *args = (bfam_subdomain_comm_args_t *)comm_args;

    if (args->user_comm_info)
      args->user_comm_info(thisSubdomain, send_sz, recv_sz, comm_args);
  }

  BFAM_LDEBUG(
      " rank %3d   ns %3jd   ms %3jd   uid %3jd   send_sz %3zd   recv_sz %3zd",
      *rank, (intmax_t)sort[0], (intmax_t)sort[1], (intmax_t)sort[2], *send_sz,
      *recv_sz);
}

static void
bfam_subdomain_dgx_field_init(bfam_subdomain_t *subdomain, const char *name,
                              bfam_real_t time,
                              bfam_subdomain_init_field_t init_field, void *arg)
{
  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t *)subdomain;

  bfam_real_t *field = bfam_dictionary_get_value_ptr(&s->base.fields, name);

  BFAM_ABORT_IF(field == NULL, "Init: Field %s not found in subdomain %s", name,
                subdomain->name);

  bfam_locidx_t fieldLength = s->Np * s->K;

  bfam_real_t *restrict x0 =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x0");
  bfam_real_t *restrict x1 =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x1");
  bfam_real_t *restrict x2 =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x2");
  /* this is just to check is we need to update the init function */
  BFAM_ASSERT(NULL ==
              bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x3"));

  init_field(fieldLength, name, time, x0, x1, x2, subdomain, arg, field);
}

static void bfam_subdomain_dgx_vtk_interp(bfam_locidx_t K, int N_d,
                                          bfam_real_t *restrict d, int N_s,
                                          const bfam_real_t *restrict s,
                                          const bfam_real_t *restrict interp,
                                          int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_gmask_set");
  const int DIM = inDIM;
#endif

  BFAM_ASSUME_ALIGNED(d, 32);
  BFAM_ASSUME_ALIGNED(s, 32);
  BFAM_ASSUME_ALIGNED(interp, 32);
  for (bfam_locidx_t elem = 0; elem < K; ++elem)
  {
    const int Np_d = bfam_ipow(N_d + 1, DIM);
    const int Np_s = bfam_ipow(N_s + 1, DIM);

    const bfam_locidx_t o_d = elem * Np_d;
    const bfam_locidx_t o_s = elem * Np_s;

    for (bfam_locidx_t n = 0; n < Np_d; n++)
      d[o_d + n] = 0;
    if (DIM == 1)
      for (int l = 0; l < N_s + 1; l++)
        for (int i = 0; i < N_d + 1; i++)
          d[o_d + i] += interp[(N_d + 1) * l + i] * s[o_s + l];
    else if (DIM == 2)
      for (int m = 0; m < N_s + 1; m++)
        for (int l = 0; l < N_s + 1; l++)
          for (int j = 0; j < N_d + 1; j++)
            for (int i = 0; i < N_d + 1; i++)
              d[o_d + j * (N_d + 1) + i] += interp[(N_d + 1) * m + j] *
                                            interp[(N_d + 1) * l + i] *
                                            s[o_s + m * (N_s + 1) + l];
    else if (DIM == 3)
      for (int n = 0; n < N_s + 1; n++)
        for (int m = 0; m < N_s + 1; m++)
          for (int l = 0; l < N_s + 1; l++)
            for (int k = 0; k < N_d + 1; k++)
              for (int j = 0; j < N_d + 1; j++)
                for (int i = 0; i < N_d + 1; i++)
                  d[o_d + k * (N_d + 1) * (N_d + 1) + j * (N_d + 1) + i] +=
                      interp[(N_d + 1) * l + i] * interp[(N_d + 1) * m + j] *
                      interp[(N_d + 1) * n + k] *
                      s[o_s + n * (N_s + 1) * (N_s + 1) + m * (N_s + 1) + l];
    else
      BFAM_ABORT("Cannot handle dim = %d", DIM);
  }
}

static int bfam_subdomain_dgx_vtk_write_vtu_piece(
    bfam_subdomain_t *subdomain, FILE *file, bfam_real_t time,
    const char **scalars, const char **vectors, const char **components,
    int writeBinary, int writeCompressed, int rank, bfam_locidx_t id,
    int Np_write)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)subdomain;

  BFAM_LDEBUG("Handling vtk for subdomain %s", subdomain->name);

  const char *format;

  if (writeBinary)
    format = "binary";
  else
    format = "ascii";

  const bfam_locidx_t K = sub->K;
  int N_vtk = sub->N;
  int Np_vtk = sub->Np;
  bfam_real_t *interp = NULL;

  bfam_real_t *restrict stor1 = NULL;
  bfam_real_t *restrict stor2 = NULL;
  bfam_real_t *restrict stor3 = NULL;

  if (Np_write > 0)
  {
    BFAM_ABORT_IF_NOT(Np_write > 1, "Np_write = %d is not valid", Np_write);

    N_vtk = Np_write - 1;
    Np_vtk = bfam_ipow(Np_write, sub->dim);

    interp =
        bfam_malloc_aligned(sizeof(bfam_real_t) * (sub->N + 1) * (N_vtk + 1));

    bfam_long_real_t *cal_interp = bfam_malloc_aligned(
        sizeof(bfam_long_real_t) * (sub->N + 1) * (N_vtk + 1));
    bfam_long_real_t *lr =
        bfam_malloc_aligned(sizeof(bfam_long_real_t) * Np_write);

    for (int r = 0; r < Np_write; r++)
      lr[r] = -1 + 2 * (bfam_long_real_t)r / (Np_write - 1);

    bfam_jacobi_p_interpolation(0, 0, sub->N, Np_write, lr, sub->lV,
                                cal_interp);

    for (int n = 0; n < (sub->N + 1) * (N_vtk + 1); n++)
      interp[n] = (bfam_real_t)cal_interp[n];

    stor1 = bfam_malloc_aligned(sizeof(bfam_real_t) * Np_vtk * K);
    stor2 = bfam_malloc_aligned(sizeof(bfam_real_t) * Np_vtk * K);
    stor3 = bfam_malloc_aligned(sizeof(bfam_real_t) * Np_vtk * K);

    bfam_free_aligned(lr);
    bfam_free_aligned(cal_interp);
  }

  const int Ncorners = sub->Ng[sub->numg - 1];

  const bfam_locidx_t Ncells = K * bfam_ipow(N_vtk, sub->dim);
  const bfam_locidx_t Ntotal = K * Np_vtk;

  bfam_real_t *restrict x =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x0");
  bfam_real_t *restrict y =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x1");
  bfam_real_t *restrict z =
      bfam_dictionary_get_value_ptr(&subdomain->fields, "_grid_x2");

  if (interp == NULL)
  {
    stor1 = x;
    stor2 = y;
    stor3 = z;
  }
  else
  {
    bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor1, sub->N, x, interp, sub->dim);
    bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor2, sub->N, y, interp, sub->dim);
    if (z != NULL)
      bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor3, sub->N, z, interp,
                                    sub->dim);
    else
      for (bfam_locidx_t k = 0; k < Np_vtk * K; k++)
        stor3[k] = 0;
  }

  fprintf(file, "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
          (intmax_t)Ntotal, (intmax_t)Ncells);

  /*
   * Points
   */
  fprintf(file, "      <Points>\n");

  bfam_vtk_write_real_vector_data_array(file, "Position", writeBinary,
                                        writeCompressed, Ntotal, stor1, stor2,
                                        stor3);

  fprintf(file, "      </Points>\n");

  /*
   * Cells
   */
  fprintf(file, "      <Cells>\n");

  /*
   * Connectivity
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"connectivity\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  if (writeBinary)
  {
    size_t cellsSize = Ncells * Ncorners * sizeof(bfam_locidx_t);
    bfam_locidx_t *cells = bfam_malloc_aligned(cellsSize);

    if (sub->dim == 1)
      for (bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for (int n = 0; n < N_vtk; ++n)
        {
          cells[i++] = Np_vtk * k + (n + 0);
          cells[i++] = Np_vtk * k + (n + 1);
        }
      }
    else if (sub->dim == 2)
      for (bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for (int m = 0; m < N_vtk; ++m)
        {
          for (int n = 0; n < N_vtk; ++n)
          {
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 0);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 1);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 0);
            cells[i++] = Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 1);
          }
        }
      }
    else if (sub->dim == 3)
      for (bfam_locidx_t k = 0, i = 0; k < K; ++k)
      {
        for (int l = 0; l < N_vtk; ++l)
        {
          for (int m = 0; m < N_vtk; ++m)
          {
            for (int n = 0; n < N_vtk; ++n)
            {
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                           (N_vtk + 1) * (m + 0) + (n + 0);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                           (N_vtk + 1) * (m + 0) + (n + 1);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                           (N_vtk + 1) * (m + 1) + (n + 0);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                           (N_vtk + 1) * (m + 1) + (n + 1);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                           (N_vtk + 1) * (m + 0) + (n + 0);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                           (N_vtk + 1) * (m + 0) + (n + 1);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                           (N_vtk + 1) * (m + 1) + (n + 0);
              cells[i++] = Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                           (N_vtk + 1) * (m + 1) + (n + 1);
            }
          }
        }
      }
    else
      BFAM_ABORT("not implemented for dim = %d", sub->dim);

    fprintf(file, "          ");
    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)cells,
                                          cellsSize);
    fprintf(file, "\n");
    if (rval)
      BFAM_WARNING("Error encoding cells");

    bfam_free_aligned(cells);
  }
  else
  {
    if (sub->dim == 1)
      for (bfam_locidx_t k = 0; k < K; ++k)
        for (int n = 0; n < N_vtk; ++n)
          fprintf(file, "          %8jd %8jd\n", (intmax_t)Np_vtk * k + (n + 0),
                  (intmax_t)Np_vtk * k + (n + 1));
    else if (sub->dim == 2)
      for (bfam_locidx_t k = 0; k < K; ++k)
        for (int m = 0; m < N_vtk; ++m)
          for (int n = 0; n < N_vtk; ++n)
            fprintf(file, "          %8jd %8jd %8jd %8jd\n",
                    (intmax_t)Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 0),
                    (intmax_t)Np_vtk * k + (N_vtk + 1) * (m + 0) + (n + 1),
                    (intmax_t)Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 0),
                    (intmax_t)Np_vtk * k + (N_vtk + 1) * (m + 1) + (n + 1));
    else if (sub->dim == 3)
      for (bfam_locidx_t k = 0; k < K; ++k)
        for (int l = 0; l < N_vtk; ++l)
          for (int m = 0; m < N_vtk; ++m)
            for (int n = 0; n < N_vtk; ++n)
              fprintf(
                  file, "          %8jd %8jd %8jd %8jd %8jd %8jd %8jd %8jd\n",
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                      (N_vtk + 1) * (m + 0) + (n + 0),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                      (N_vtk + 1) * (m + 0) + (n + 1),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                      (N_vtk + 1) * (m + 1) + (n + 0),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 0) +
                      (N_vtk + 1) * (m + 1) + (n + 1),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                      (N_vtk + 1) * (m + 0) + (n + 0),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                      (N_vtk + 1) * (m + 0) + (n + 1),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                      (N_vtk + 1) * (m + 1) + (n + 0),
                  (intmax_t)Np_vtk * k + (N_vtk + 1) * (N_vtk + 1) * (l + 1) +
                      (N_vtk + 1) * (m + 1) + (n + 1));
    else
      BFAM_ABORT("not implemented for dim = %d %d", sub->dim);
  }
  fprintf(file, "        </DataArray>\n");

  /*
   * Offsets
   */
  fprintf(file, "        <DataArray type=\"%s\" Name=\"offsets\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t offsetsSize = Ncells * sizeof(bfam_locidx_t);
    bfam_locidx_t *offsets = bfam_malloc_aligned(offsetsSize);

    for (bfam_locidx_t i = 1; i <= Ncells; ++i)
      offsets[i - 1] = Ncorners * i;

    int rval = bfam_vtk_write_binary_data(writeCompressed, file,
                                          (char *)offsets, offsetsSize);
    if (rval)
      BFAM_WARNING("Error encoding offsets");

    bfam_free_aligned(offsets);
  }
  else
  {
    for (bfam_locidx_t i = 1, sk = 1; i <= Ncells; ++i, ++sk)
    {
      fprintf(file, " %8jd", (intmax_t)(Ncorners * i));
      if (!(sk % 20) && i != Ncells)
        fprintf(file, "\n          ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  /*
   * Types
   */
  fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\""
                " format=\"%s\">\n",
          format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t typesSize = Ncells * sizeof(uint8_t);
    uint8_t *types = bfam_malloc_aligned(typesSize);

    if (sub->dim == 1)
      for (bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 3; /* VTK_LINE */
    else if (sub->dim == 2)
      for (bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 8; /* VTK_PIXEL */
    else if (sub->dim == 3)
      for (bfam_locidx_t i = 0; i < Ncells; ++i)
        types[i] = 11; /* VTK_VOXEL */
    else
      BFAM_ABORT("cannot handle dim = %d", sub->dim);

    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)types,
                                          typesSize);
    if (rval)
      BFAM_WARNING("Error encoding types");

    bfam_free_aligned(types);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      if (sub->dim == 1)
        fprintf(file, " 3"); /* VTK_LINE */
      else if (sub->dim == 2)
        fprintf(file, " 8"); /* VTK_PIXEL */
      else if (sub->dim == 3)
        fprintf(file, " 11"); /* VTK_VOXEL */
      else
        BFAM_ABORT("cannot handle dim = %d", sub->dim);
      if (!(sk % 20) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "      </Cells>\n");

  /*
   * Cell Data
   */
  fprintf(file, "      <CellData Scalars=\"time,mpirank,subdomain_id\">\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"time\""
                " format=\"%s\">\n",
          BFAM_REAL_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t timesize = Ncells * sizeof(bfam_real_t);
    bfam_real_t *times = bfam_malloc_aligned(timesize);

    for (bfam_locidx_t i = 0; i < Ncells; ++i)
      times[i] = time;

    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)times,
                                          timesize);
    if (rval)
      BFAM_WARNING("Error encoding times");

    bfam_free_aligned(times);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %" BFAM_REAL_FMTe, time);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"mpirank\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t ranksSize = Ncells * sizeof(bfam_locidx_t);
    bfam_locidx_t *ranks = bfam_malloc_aligned(ranksSize);

    for (bfam_locidx_t i = 0; i < Ncells; ++i)
      ranks[i] = rank;

    int rval = bfam_vtk_write_binary_data(writeCompressed, file, (char *)ranks,
                                          ranksSize);
    if (rval)
      BFAM_WARNING("Error encoding ranks");

    bfam_free_aligned(ranks);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)rank);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  fprintf(file, "        <DataArray type=\"%s\" Name=\"subdomain_id\""
                " format=\"%s\">\n",
          BFAM_LOCIDX_VTK, format);
  fprintf(file, "          ");
  if (writeBinary)
  {
    size_t idsSize = Ncells * sizeof(bfam_locidx_t);
    bfam_locidx_t *ids = bfam_malloc_aligned(idsSize);

    for (bfam_locidx_t i = 0; i < Ncells; ++i)
      ids[i] = subdomain->id;

    int rval =
        bfam_vtk_write_binary_data(writeCompressed, file, (char *)ids, idsSize);
    if (rval)
      BFAM_WARNING("Error encoding ids");

    bfam_free_aligned(ids);
  }
  else
  {
    for (bfam_locidx_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)subdomain->id);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");

  fprintf(file, "      </CellData>\n");

  char pointscalars[BFAM_BUFSIZ];
  bfam_util_strcsl(pointscalars, scalars);

  char pointvectors[BFAM_BUFSIZ];
  bfam_util_strcsl(pointvectors, vectors);

  fprintf(file, "      <PointData Scalars=\"%s\" Vectors=\"%s\">\n",
          pointscalars, pointvectors);

  if (scalars)
  {
    for (size_t s = 0; scalars[s]; ++s)
    {
      bfam_real_t *sdata =
          bfam_dictionary_get_value_ptr(&subdomain->fields, scalars[s]);
      BFAM_ABORT_IF(sdata == NULL, "VTK: Field %s not in subdomain %s",
                    scalars[s], subdomain->name);
      if (interp == NULL)
      {
        stor1 = sdata;
      }
      else
      {
        bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor1, sub->N, sdata, interp,
                                      sub->dim);
      }

      bfam_vtk_write_real_scalar_data_array(file, scalars[s], writeBinary,
                                            writeCompressed, Ntotal, stor1);
    }
  }

  if (vectors)
  {
    for (size_t v = 0; vectors[v]; ++v)
    {

      bfam_real_t *v1 = bfam_dictionary_get_value_ptr(&subdomain->fields,
                                                      components[3 * v + 0]);
      bfam_real_t *v2 = bfam_dictionary_get_value_ptr(&subdomain->fields,
                                                      components[3 * v + 1]);
      bfam_real_t *v3 = bfam_dictionary_get_value_ptr(&subdomain->fields,
                                                      components[3 * v + 2]);

      BFAM_ABORT_IF(v1 == NULL, "VTK: Field %s not in subdomain %s",
                    components[3 * v + 0], subdomain->name);
      BFAM_ABORT_IF(v2 == NULL, "VTK: Field %s not in subdomain %s",
                    components[3 * v + 1], subdomain->name);
      BFAM_ABORT_IF(v3 == NULL, "VTK: Field %s not in subdomain %s",
                    components[3 * v + 2], subdomain->name);
      if (interp == NULL)
      {
        stor1 = v1;
        stor2 = v2;
        stor3 = v3;
      }
      else
      {
        bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor1, sub->N, v1, interp,
                                      sub->dim);
        bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor2, sub->N, v2, interp,
                                      sub->dim);
        bfam_subdomain_dgx_vtk_interp(K, N_vtk, stor3, sub->N, v3, interp,
                                      sub->dim);
      }

      bfam_vtk_write_real_vector_data_array(file, vectors[v], writeBinary,
                                            writeCompressed, Ntotal, stor1,
                                            stor2, stor3);
    }
  }

  fprintf(file, "      </PointData>\n");
  fprintf(file, "    </Piece>\n");

  if (interp != NULL)
  {
    bfam_free_aligned(interp);
    bfam_free_aligned(stor1);
    bfam_free_aligned(stor2);
    bfam_free_aligned(stor3);
  }
  return 1;
}

static int bfam_subdomain_dgx_field_face_add(bfam_subdomain_t *subdomain,
                                             const char *name)
{
  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t *)subdomain;

  if (bfam_dictionary_get_value_ptr(&s->base.fields_face, name))
    return 1;

  size_t fieldSize = s->Ng[0] * s->Ngp[0] * s->K * sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);

  int rval = bfam_dictionary_insert_ptr(&s->base.fields_face, name, field);

  BFAM_ASSERT(rval != 1);

  if (rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static int bfam_subdomain_dgx_field_add(bfam_subdomain_t *subdomain,
                                        const char *name)
{
  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t *)subdomain;

  if (bfam_dictionary_get_value_ptr(&s->base.fields, name))
    return 1;

  size_t fieldSize = s->Np * s->K * sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);
#ifdef BFAM_DEBUG
  for (int i = 0; i < s->Np * s->K; i++)
    field[i] = bfam_real_nan("");
#endif

  int rval = bfam_dictionary_insert_ptr(&s->base.fields, name, field);

  BFAM_ASSERT(rval != 1);

  if (rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static inline int ***bfam_subdomain_dgx_gmask_set(const int numg, const int N,
                                                  int *Np, int *Ng, int *Ngp,
                                                  int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_gmask_set");
  const int DIM = inDIM;
#endif
  BFAM_ABORT_IF(DIM > 3 || DIM < 0,
                "bfam_subdomain_dgx_gmask_set cannot handle dim = %d", DIM);

  if (DIM == 0)
  {
    *Np = 1;
    return NULL;
  }

  /* this could probably be made generic for arbitrary dimensions, but until
   * that's needed... */
  switch (DIM)
  {
  case 1:
    *Np = N + 1;

    /* just corners */
    Ng[0] = 2;
    Ngp[0] = 1;
    break;

  case 2:
    *Np = (N + 1) * (N + 1);

    /* edges */
    Ng[0] = 4;
    Ngp[0] = N + 1;

    /* corners */
    Ng[1] = 4;
    Ngp[1] = 1;
    break;

  case 3:
    *Np = (N + 1) * (N + 1) * (N + 1);

    /* faces */
    Ng[0] = 6;
    Ngp[0] = (N + 1) * (N + 1);

    /* edges */
    Ng[1] = 12;
    Ngp[1] = N + 1;

    /* corners */
    Ng[2] = 8;
    Ngp[2] = 1;
    break;

  default:
    BFAM_ABORT("cannot handle dim = %d", DIM);
  }

  int ***gmask = bfam_malloc_aligned(numg * sizeof(int **));
  for (int g = 0; g < numg; g++)
  {
    gmask[g] = bfam_malloc_aligned(Ng[g] * sizeof(int *));
    for (int i = 0; i < Ng[g]; i++)
      gmask[g][i] = bfam_malloc_aligned(Ngp[g] * sizeof(int));
  }

  switch (DIM)
  {
  case 1:
    gmask[0][0][0] = 0;
    gmask[0][1][0] = N;
    break;

  case 2:
    /* edges */
    for (int i = 0; i < N + 1; ++i)
      gmask[0][0][i] = i * (N + 1);
    for (int i = 0; i < N + 1; ++i)
      gmask[0][1][i] = (i + 1) * (N + 1) - 1;
    for (int i = 0; i < N + 1; ++i)
      gmask[0][2][i] = i;
    for (int i = 0; i < N + 1; ++i)
      gmask[0][3][i] = (N + 1) * N + i;

    /* corners */
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 2; ++i)
        gmask[1][i + j * 2][0] = i * N + j * (N + 1) * N;
    break;

  case 3:
    /* This could all probably be cleaned up... */

    /* faces */
    {
      int n, i, j, k, f = -1;

      n = 0;
      i = 0;
      f++;
      for (k = 0; k < N + 1; k++)
        for (j = 0; j < N + 1; j++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      i = N;
      f++;
      for (k = 0; k < N + 1; k++)
        for (j = 0; j < N + 1; j++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      j = 0;
      f++;
      for (k = 0; k < N + 1; k++)
        for (i = 0; i < N + 1; i++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      j = N;
      f++;
      for (k = 0; k < N + 1; k++)
        for (i = 0; i < N + 1; i++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      k = 0;
      f++;
      for (j = 0; j < N + 1; j++)
        for (i = 0; i < N + 1; i++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);

      n = 0;
      k = N;
      f++;
      for (j = 0; j < N + 1; j++)
        for (i = 0; i < N + 1; i++)
          gmask[0][f][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);
    }

    /* edges */
    {
      int n, i, j, k, e = 0;

      for (k = 0; k < N + 1; k += N)
        for (j = 0; j < N + 1; j += N)
        {
          n = 0;
          for (i = 0; i < N + 1; i++)
            gmask[1][e][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);
          e++;
        }
      for (k = 0; k < N + 1; k += N)
        for (i = 0; i < N + 1; i += N)
        {
          n = 0;
          for (j = 0; j < N + 1; j++)
            gmask[1][e][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);
          e++;
        }
      for (j = 0; j < N + 1; j += N)
        for (i = 0; i < N + 1; i += N)
        {
          n = 0;
          for (k = 0; k < N + 1; k++)
            gmask[1][e][n++] = i + j * (N + 1) + k * (N + 1) * (N + 1);
          e++;
        }
    }

    /* corners */
    for (int k = 0, c = 0; k < N + 1; k += N)
      for (int j = 0; j < N + 1; j += N)
        for (int i = 0; i < N + 1; i += N)
          gmask[2][c++][0] = i + j * (N + 1) + k * (N + 1) * (N + 1);

    break;

  default:
    BFAM_ABORT("cannot handle dim = %d", DIM);
  }

  return gmask;
}

static void bfam_subdomain_dgx_buildmaps(
    int N, bfam_locidx_t K, int Np, int Nfp, int Nfaces,
    const bfam_locidx_t *EToE, const int8_t *EToF, int ***gmask,
    bfam_locidx_t *restrict vmapP, bfam_locidx_t *restrict vmapM, int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_buildmaps");
  const int DIM = inDIM;
#endif

  BFAM_ASSERT(DIM == inDIM);

  for (bfam_locidx_t k1 = 0, sk = 0; k1 < K; ++k1)
  {
    for (int8_t f1 = 0; f1 < Nfaces; ++f1)
    {
      bfam_locidx_t k2 = EToE[Nfaces * k1 + f1];
      int8_t f2 = (int8_t)(EToF[Nfaces * k1 + f1] % Nfaces);
      int8_t o = (int8_t)(EToF[Nfaces * k1 + f1] / Nfaces);

      for (int n = 0; n < Nfp; ++n)
      {
        vmapM[sk + n] = Np * k1 + gmask[0][f1][n];
      }

      switch (DIM)
      {
      case 1:
        /* Orientation does not matter in 1D */
        for (int n = 0; n < Nfp; ++n)
        {
          vmapP[sk + n] = Np * k2 + gmask[0][f2][n];
        }
        break;
      case 2:
        for (int n = 0; n < Nfp; ++n)
        {
          if (o)
            vmapP[sk + n] = Np * k2 + gmask[0][f2][Nfp - 1 - n];
          else
            vmapP[sk + n] = Np * k2 + gmask[0][f2][n];
        }
        break;
      case 3:
      {
        const int Nrp = N + 1;
        BFAM_ASSERT(Nfp == Nrp * Nrp);

        int oidx = -1;

        const int8_t nr = bfam_p8est_FToF_code[f1][f2];
        const int8_t ns =
            (k1 == k2 && f1 == f2) ? 0 : bfam_p8est_code_to_perm[nr][o];

        for (int j = 0, n = 0; j < Nrp; ++j)
        {
          for (int i = 0; i < Nrp; ++i, ++n)
          {
            int ir = Nrp - (i + 1);
            int jr = Nrp - (j + 1);
            switch (ns)
            {
            case 0:
              oidx = i + j * Nrp;
              break;
            case 1:
              oidx = j + i * Nrp;
              break;
            case 2:
              oidx = ir + j * Nrp;
              break;
            case 3:
              oidx = jr + i * Nrp;
              break;
            case 4:
              oidx = j + ir * Nrp;
              break;
            case 5:
              oidx = i + jr * Nrp;
              break;
            case 6:
              oidx = jr + ir * Nrp;
              break;
            case 7:
              oidx = ir + jr * Nrp;
              break;
            default:
              BFAM_ABORT("invalid orientation %d", o);
            }
            vmapP[sk + n] = Np * k2 + gmask[0][f2][oidx];
          }
        }
      }
      break;
      default:
        BFAM_ABORT("cannot handle dim = %d", DIM);
      }

      sk += Nfp;
    }
  }
}

/* set all subdomain values to something logical */
static void bfam_subdomain_dgx_null_all_values(bfam_subdomain_dgx_t *sub)
{
  sub->K = 0;
  sub->N = -1;
  sub->Np = 0;
  sub->Ngp = NULL;
  sub->numg = 0;
  sub->Ng = NULL;
  sub->r = NULL;
  sub->w = NULL;
  sub->wi = NULL;
  sub->Dr = NULL;
  sub->lDr = NULL;
  sub->lr = NULL;
  sub->lw = NULL;
  sub->lV = NULL;
  sub->K = 0;
  sub->vmapM = NULL;
  sub->vmapP = NULL;
  sub->gmask = NULL;
  sub->EToQ = NULL;
}

static void
bfam_subdomain_dgx_generic_init(bfam_subdomain_dgx_t *subdomain,
                                const bfam_locidx_t id, const bfam_locidx_t uid,
                                const char *name, const int N,
                                const bfam_locidx_t K, bfam_dictionary_t *N2N,
                                bfam_dictionary_t *dgx_ops, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_init");
  const int DIM = inDIM;
#endif

  BFAM_ASSERT(DIM == inDIM);
  BFAM_ABORT_IF(DIM < 0, "dimension %d is not possible in bfam", DIM);
  BFAM_ABORT_IF(DIM == 0 && N != 0,
                "if DIM < 1 then N must be zero (i.e., constant");

  bfam_subdomain_init(&subdomain->base, id, uid, name);
  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");
  char dim_str[BFAM_BUFSIZ];
  snprintf(dim_str, BFAM_BUFSIZ, "_dimension_%d", DIM);
  bfam_subdomain_add_tag(&subdomain->base, dim_str);

  bfam_subdomain_dgx_null_all_values(subdomain);

  subdomain->dim = DIM;

  subdomain->base.free =
      BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_, BFAM_DGX_DIMENSION);
  subdomain->base.vtk_write_vtu_piece = bfam_subdomain_dgx_vtk_write_vtu_piece;
  subdomain->base.field_add = bfam_subdomain_dgx_field_add;
  subdomain->base.field_face_add = bfam_subdomain_dgx_field_face_add;
  subdomain->base.field_init = bfam_subdomain_dgx_field_init;
  subdomain->base.glue_comm_info = bfam_subdomain_dgx_comm_info;
  subdomain->base.glue_put_send_buffer = bfam_subdomain_dgx_put_send_buffer;
  subdomain->base.glue_get_recv_buffer = bfam_subdomain_dgx_get_recv_buffer;

  subdomain->numg = DIM;
  const int numg = subdomain->numg;

  int *Ng = NULL;
  if (numg > 0)
    Ng = bfam_malloc_aligned(sizeof(int) * numg);
  subdomain->Ng = Ng;

  int *Ngp = NULL;
  if (numg > 0)
    Ngp = bfam_malloc_aligned(sizeof(int) * numg);
  subdomain->Ngp = Ngp;

  subdomain->gmask =
      bfam_subdomain_dgx_gmask_set(numg, N, &subdomain->Np, Ng, Ngp, DIM);

  subdomain->hadapt = NULL;
  subdomain->padapt = NULL;
  subdomain->q_id = NULL;
  subdomain->lvl = NULL;

  if (DIM > 0)
  {
    subdomain->K = K;
    subdomain->N = N;
    subdomain->hadapt = bfam_malloc_aligned(K * sizeof(uint8_t));
    subdomain->padapt = bfam_malloc_aligned(K * sizeof(int8_t));
    subdomain->q_id = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
    subdomain->lvl = bfam_malloc_aligned(K * sizeof(int8_t));
    for (bfam_locidx_t k = 0; k < K; ++k)
    {
      subdomain->hadapt[k] = BFAM_FLAG_SAME;
      subdomain->padapt[k] = (int8_t)N;
      subdomain->q_id[k] = -1;
      subdomain->lvl[k] = -1;
    }

    /*
     * TODO: query dictionary to see if we need to create this, otherwise used
     * stored values
     */
    BFAM_ASSERT(dgx_ops);

    char name[BFAM_BUFSIZ];
    snprintf(name, BFAM_BUFSIZ, "lr_%d", N);
    if (!bfam_dictionary_contains(dgx_ops, name))
    {

      const int Nrp = N + 1;
      bfam_long_real_t *lr =
          bfam_malloc_aligned(Nrp * sizeof(bfam_long_real_t));
      bfam_long_real_t *lw =
          bfam_malloc_aligned(Nrp * sizeof(bfam_long_real_t));
      bfam_long_real_t *lV =
          bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_long_real_t));
      bfam_long_real_t *lDr =
          bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_long_real_t));
      bfam_real_t *Dr = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_real_t));
      bfam_real_t *r = bfam_malloc_aligned(Nrp * sizeof(bfam_real_t));
      bfam_real_t *w = bfam_malloc_aligned(Nrp * sizeof(bfam_real_t));
      bfam_real_t *wi = bfam_malloc_aligned(Nrp * sizeof(bfam_real_t));

      bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);
      bfam_jacobi_p_vandermonde(0, 0, N, Nrp, lr, lV);
      bfam_jacobi_p_differentiation(0, 0, N, Nrp, lr, lV, lDr);

      /* store the volume stuff */
      for (int n = 0; n < Nrp; ++n)
      {
        r[n] = (bfam_real_t)lr[n];
        w[n] = (bfam_real_t)lw[n];
        wi[n] = (bfam_real_t)(1.0l / lw[n]);
      }
      for (int n = 0; n < Nrp * Nrp; ++n)
      {
        Dr[n] = (bfam_real_t)lDr[n];
      }

      int rval = 1;

      snprintf(name, BFAM_BUFSIZ, "lr_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, lr);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "lw_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, lw);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "lV_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, lV);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "lDr_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, lDr);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "Dr_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, Dr);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "r_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, r);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "w_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, w);
      BFAM_ASSERT(rval != 1);

      snprintf(name, BFAM_BUFSIZ, "wi_%d", N);
      rval = bfam_dictionary_insert_ptr(dgx_ops, name, wi);
      BFAM_ASSERT(rval != 1);
    }

    snprintf(name, BFAM_BUFSIZ, "lr_%d", N);
    subdomain->lr = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->lr != NULL);

    snprintf(name, BFAM_BUFSIZ, "lw_%d", N);
    subdomain->lw = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->lw != NULL);

    snprintf(name, BFAM_BUFSIZ, "lV_%d", N);
    subdomain->lV = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->lV != NULL);

    snprintf(name, BFAM_BUFSIZ, "lDr_%d", N);
    subdomain->lDr = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->lDr != NULL);

    snprintf(name, BFAM_BUFSIZ, "Dr_%d", N);
    subdomain->Dr = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->Dr != NULL);

    snprintf(name, BFAM_BUFSIZ, "r_%d", N);
    subdomain->r = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->r != NULL);

    snprintf(name, BFAM_BUFSIZ, "w_%d", N);
    subdomain->w = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->w != NULL);

    snprintf(name, BFAM_BUFSIZ, "wi_%d", N);
    subdomain->wi = bfam_dictionary_get_value_ptr(dgx_ops, name);
    BFAM_ASSERT(subdomain->wi != NULL);
  }
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_, BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_t *subdomain, const bfam_locidx_t id,
    const bfam_locidx_t uid, const char *name, const int N,
    const bfam_locidx_t K, const bfam_locidx_t *EToQ, const bfam_locidx_t *EToE,
    const int8_t *EToF, bfam_dictionary_t *N2N, bfam_dictionary_t *dgx_ops,
    const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_init");
  const int DIM = inDIM;
#endif
  bfam_subdomain_dgx_generic_init(subdomain, id, uid, name, N, K, N2N, dgx_ops,
                                  inDIM);

  const int *Ng = subdomain->Ng;
  const int *Ngp = subdomain->Ngp;
  const int Np = subdomain->Np;

  if (EToQ)
  {
    subdomain->EToQ = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
    for (bfam_locidx_t k = 0; k < K; k++)
      subdomain->EToQ[k] = EToQ[k];
  }

  if (DIM > 0)
  {
    /* store the face stuff */
    subdomain->vmapP =
        bfam_malloc_aligned(K * Ngp[0] * Ng[0] * sizeof(bfam_locidx_t));
    subdomain->vmapM =
        bfam_malloc_aligned(K * Ngp[0] * Ng[0] * sizeof(bfam_locidx_t));

    bfam_subdomain_dgx_buildmaps(N, K, Np, Ngp[0], Ng[0], EToE, EToF,
                                 subdomain->gmask, subdomain->vmapP,
                                 subdomain->vmapM, DIM);
  }
}

bfam_subdomain_dgx_t *BFAM_APPEND_EXPAND(bfam_subdomain_dgx_new_,
                                         BFAM_DGX_DIMENSION)(
    const bfam_locidx_t id, const bfam_locidx_t uid, const char *name,
    const int N, const bfam_locidx_t K, const bfam_locidx_t *EToQ,
    const bfam_locidx_t *EToE, const int8_t *EToF, bfam_dictionary_t *N2N,
    bfam_dictionary_t *dgx_ops, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_new");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);

  bfam_subdomain_dgx_t *newSubdomain =
      bfam_malloc(sizeof(bfam_subdomain_dgx_t));

  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_, BFAM_DGX_DIMENSION)(
      newSubdomain, id, uid, name, N, K, EToQ, EToE, EToF, N2N, dgx_ops, DIM);
  return newSubdomain;
}

static int bfam_subdomain_dgx_free_fields(const char *key, void *val, void *arg)
{
  bfam_free_aligned(val);

  return 1;
}

static void bfam_subdomain_dgx_free_glue(bfam_subdomain_dgx_glue_data_t *glue)
{
  if (glue)
  {
    if (glue->interpolation)
      bfam_free_aligned(glue->interpolation);
    if (glue->massprojection)
      bfam_free_aligned(glue->massprojection);
    if (glue->projection)
      bfam_free_aligned(glue->projection);
    if (glue->exact_mass)
      bfam_free_aligned(glue->exact_mass);

    if (glue->EToEp)
      bfam_free_aligned(glue->EToEp);
    if (glue->EToHp)
      bfam_free_aligned(glue->EToHp);
    if (glue->EToEm)
      bfam_free_aligned(glue->EToEm);
    if (glue->EToFm)
      bfam_free_aligned(glue->EToFm);
    if (glue->EToHm)
      bfam_free_aligned(glue->EToHm);
    if (glue->EToOm)
      bfam_free_aligned(glue->EToOm);
    if (glue->mapOm)
    {
      for (int n = 0; n < glue->num_orient; n++)
        bfam_free_aligned(glue->mapOm[n]);
      bfam_free_aligned(glue->mapOm);
    }
    bfam_critbit0_clear(&glue->base.tags);
    bfam_free(glue);
  }
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_free_,
                        BFAM_DGX_DIMENSION)(bfam_subdomain_t *thisSubdomain)
{
  bfam_subdomain_dgx_t *sub = (bfam_subdomain_dgx_t *)thisSubdomain;

  bfam_dictionary_allprefixed_ptr(&sub->base.fields, "",
                                  &bfam_subdomain_dgx_free_fields, NULL);
  if (sub->base.glue_p)
    bfam_dictionary_allprefixed_ptr(&sub->base.glue_p->fields, "",
                                    &bfam_subdomain_dgx_free_fields, NULL);
  if (sub->base.glue_m)
    bfam_dictionary_allprefixed_ptr(&sub->base.glue_m->fields, "",
                                    &bfam_subdomain_dgx_free_fields, NULL);
  bfam_dictionary_allprefixed_ptr(&sub->base.fields_face, "",
                                  &bfam_subdomain_dgx_free_fields, NULL);

  bfam_subdomain_free(thisSubdomain);

  if (sub->gmask)
  {
    for (int g = 0; g < sub->numg; g++)
    {
      for (int i = 0; i < sub->Ng[g]; i++)
        bfam_free_aligned(sub->gmask[g][i]);
      bfam_free_aligned(sub->gmask[g]);
    }
    bfam_free_aligned(sub->gmask);
    sub->gmask = NULL;
  }
  if (sub->Ng)
    bfam_free_aligned(sub->Ng);
  sub->Ng = NULL;
  if (sub->Ngp)
    bfam_free_aligned(sub->Ngp);
  sub->Ngp = NULL;

  if (sub->EToQ)
    bfam_free_aligned(sub->EToQ);
  sub->EToQ = NULL;
  if (sub->vmapP)
    bfam_free_aligned(sub->vmapP);
  sub->vmapP = NULL;
  if (sub->vmapM)
    bfam_free_aligned(sub->vmapM);
  sub->vmapM = NULL;

  if (sub->hadapt)
    bfam_free_aligned(sub->hadapt);
  sub->hadapt = NULL;
  if (sub->padapt)
    bfam_free_aligned(sub->padapt);
  sub->padapt = NULL;

  if (sub->q_id)
    bfam_free_aligned(sub->q_id);
  sub->q_id = NULL;

  if (sub->lvl)
    bfam_free_aligned(sub->lvl);
  sub->lvl = NULL;

  bfam_subdomain_dgx_null_all_values(sub);

  bfam_subdomain_dgx_free_glue(
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_m);
  bfam_subdomain_dgx_free_glue(
      (bfam_subdomain_dgx_glue_data_t *)sub->base.glue_p);
}

bfam_subdomain_dgx_t *BFAM_APPEND_EXPAND(bfam_subdomain_dgx_glue_new_,
                                         BFAM_DGX_DIMENSION)(
    const bfam_locidx_t id, const bfam_locidx_t uid, const char *name,
    const int N_m, const int N_p, const int N_g, const bfam_locidx_t rank_m,
    const bfam_locidx_t rank_p, const bfam_locidx_t id_m,
    const bfam_locidx_t id_p, bfam_subdomain_dgx_t *sub_m,
    bfam_locidx_t *ktok_m, const bfam_locidx_t K,
    bfam_subdomain_face_map_entry_t *mapping, bfam_dictionary_t *N2N,
    bfam_dictionary_t *dgx_ops, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_glue_new");
#ifdef BFAM_DEBUG
  const int DIM = inDIM;
#endif
#endif
  BFAM_ASSERT(DIM == inDIM);

  bfam_subdomain_dgx_t *newSubdomain =
      bfam_malloc(sizeof(bfam_subdomain_dgx_t));

  BFAM_APPEND_EXPAND(bfam_subdomain_dgx_glue_init_, BFAM_DGX_DIMENSION)(
      newSubdomain, id, uid, name, N_m, N_p, N_g, rank_m, rank_p, id_m, id_p,
      sub_m, ktok_m, K, mapping, N2N, dgx_ops, inDIM);
  return newSubdomain;
}

static int bfam_subdomain_dgx_glue_field_minus_add(bfam_subdomain_t *subdomain,
                                                   const char *name)
{
  BFAM_ASSERT(subdomain->glue_m);

  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t *)subdomain;

  if (bfam_dictionary_get_value_ptr(&s->base.glue_m->fields, name))
    return 1;

  size_t fieldSize = s->Np * s->K * sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);
#ifdef BFAM_DEBUG
  for (int i = 0; i < s->Np * s->K; i++)
    field[i] = bfam_real_nan("");
#endif

  int rval = bfam_dictionary_insert_ptr(&s->base.glue_m->fields, name, field);

  BFAM_ASSERT(rval != 1);

  if (rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static int bfam_subdomain_dgx_glue_field_plus_add(bfam_subdomain_t *subdomain,
                                                  const char *name)
{
  BFAM_ASSERT(subdomain->glue_p);

  bfam_subdomain_dgx_t *s = (bfam_subdomain_dgx_t *)subdomain;

  if (bfam_dictionary_get_value_ptr(&s->base.glue_p->fields, name))
    return 1;

  size_t fieldSize = s->Np * s->K * sizeof(bfam_real_t);
  bfam_real_t *field = bfam_malloc_aligned(fieldSize);
#ifdef BFAM_DEBUG
  for (int i = 0; i < s->Np * s->K; i++)
    field[i] = bfam_real_nan("");
#endif

  int rval = bfam_dictionary_insert_ptr(&s->base.glue_p->fields, name, field);

  BFAM_ASSERT(rval != 1);

  if (rval == 0)
    bfam_free_aligned(field);

  return rval;
}

static void bfam_subdomain_dgx_glue_generic_init(
    bfam_subdomain_dgx_glue_data_t *glue, const bfam_locidx_t rank,
    const bfam_locidx_t id, const bfam_locidx_t id_s,
    bfam_subdomain_dgx_t *sub_m, int inDIM)
{
  bfam_subdomain_glue_init(&glue->base, rank, id, id_s,
                           (bfam_subdomain_t *)sub_m);
  glue->EToEp = NULL;
  glue->EToHp = NULL;
  glue->EToEm = NULL;
  glue->EToFm = NULL;
  glue->EToHm = NULL;
  glue->EToOm = NULL;
  glue->mapOm = NULL;
  glue->num_orient = 0;
  glue->num_interp = 0;
  glue->interpolation = NULL;
  glue->projection = NULL;
  glue->massprojection = NULL;
  glue->exact_mass = NULL;
  glue->base.tags.root = NULL;
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_glue_init_, BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_t *subdomain, const bfam_locidx_t id,
    const bfam_locidx_t uid, const char *name, const int N_m, const int N_p,
    const int N_g, const bfam_locidx_t rank_m, const bfam_locidx_t rank_p,
    const bfam_locidx_t id_m, const bfam_locidx_t id_p,
    bfam_subdomain_dgx_t *sub_m, bfam_locidx_t *ktok_m, const bfam_locidx_t K,
    bfam_subdomain_face_map_entry_t *mapping, bfam_dictionary_t *N2N,
    bfam_dictionary_t *dgx_ops, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_glue_init");
  const int DIM = inDIM;
#endif
  BFAM_ASSERT(DIM == inDIM);
  BFAM_ASSERT(DIM > 0);

  bfam_subdomain_dgx_generic_init(subdomain, id, uid, name, N_g, K, N2N,
                                  dgx_ops, inDIM);

#ifdef BFAM_DEBUG
  {
    BFAM_LDEBUG("glue mapping: %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s",
                "np", "nk", "nf", "nh", "gi", "i", "k", "f", "h", "o", "id");
    for (bfam_locidx_t k = 0; k < K; ++k)
    {
      BFAM_LDEBUG("glue mapping: %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd %5jd "
                  "%5jd %5jd",
                  (intmax_t)mapping[k].np, (intmax_t)mapping[k].nk,
                  (intmax_t)mapping[k].nf, (intmax_t)mapping[k].nh,
                  (intmax_t)mapping[k].gi, (intmax_t)mapping[k].i,
                  (intmax_t)mapping[k].k, (intmax_t)mapping[k].f,
                  (intmax_t)mapping[k].h, (intmax_t)mapping[k].o,
                  (intmax_t)mapping[k].id);
    }
  }
#endif

  BFAM_ABORT_IF(N_g < N_m || N_g < N_p,
                "glue space must currently be a higher order space"
                " N_g = %d, N_m = %d, N_p = %d",
                N_g, N_m, N_p);

  subdomain->base.field_minus_add = bfam_subdomain_dgx_glue_field_minus_add;
  subdomain->base.field_plus_add = bfam_subdomain_dgx_glue_field_plus_add;

  bfam_subdomain_add_tag(&subdomain->base, "_subdomain_dgx");

  BFAM_ASSERT(subdomain->base.glue_m == NULL);

  subdomain->base.glue_m = bfam_malloc(sizeof(bfam_subdomain_dgx_glue_data_t));
  bfam_subdomain_dgx_glue_data_t *glue_m =
      (bfam_subdomain_dgx_glue_data_t *)subdomain->base.glue_m;
  bfam_subdomain_dgx_glue_generic_init(
      glue_m, rank_m, id_m, (bfam_locidx_t)(imaxabs(id_m) - 1), sub_m, DIM);

  subdomain->base.glue_p = bfam_malloc(sizeof(bfam_subdomain_dgx_glue_data_t));
  bfam_subdomain_dgx_glue_data_t *glue_p =
      (bfam_subdomain_dgx_glue_data_t *)subdomain->base.glue_p;
  bfam_subdomain_dgx_glue_generic_init(
      glue_p, rank_p, id_p, (bfam_locidx_t)(imaxabs(id_p) - 1), NULL, DIM);

  const int num_interp = 3;
  glue_m->num_interp = num_interp;

  const int N = subdomain->N;
  const int Nrp = N + 1;

  bfam_subdomain_dgx_interpolator_t *proj_m2g =
      bfam_subdomain_dgx_get_interpolator(N2N, N_m, N, DIM);
  bfam_subdomain_dgx_interpolator_t *proj_g2m =
      bfam_subdomain_dgx_get_interpolator(N2N, N, N_m, DIM);

  BFAM_ASSERT(proj_m2g);
  BFAM_ASSERT(proj_g2m);

  glue_m->interpolation =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));
  glue_m->interpolation[0] = proj_m2g->prj[0];
  glue_m->interpolation[1] = proj_m2g->prj[3];
  glue_m->interpolation[2] = proj_m2g->prj[4];

  glue_m->projection =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));
  glue_m->projection[0] = proj_g2m->prj[0];
  glue_m->projection[1] = proj_g2m->prj[1];
  glue_m->projection[2] = proj_g2m->prj[2];

  glue_m->massprojection =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));
  glue_m->massprojection[0] = proj_g2m->mass_prj[0];
  glue_m->massprojection[1] = proj_g2m->mass_prj[1];
  glue_m->massprojection[2] = proj_g2m->mass_prj[2];
  const int sub_m_Nrp = N_m + 1;
  for (int i = 0; i < 3; i++)
  {
    if (glue_m->massprojection[i])
      for (int n = 0; n < Nrp * sub_m_Nrp; ++n)
        BFAM_INFO("MP: %d :: (%d, %d) :: %d :: %e", (int)i, (int)Nrp,
                  (int)sub_m_Nrp, (int)n, (double)glue_m->massprojection[i][n]);
    if (glue_m->projection[i])
      for (int n = 0; n < Nrp * sub_m_Nrp; ++n)
        BFAM_INFO(" P: %d :: (%d, %d) :: %d :: %e", (int)i, (int)Nrp,
                  (int)sub_m_Nrp, (int)n, (double)glue_m->projection[i][n]);
  }

#if 0


  const bfam_long_real_t projection_scale[3] = {
      BFAM_LONG_REAL(1.0), BFAM_LONG_REAL(0.5), BFAM_LONG_REAL(0.5)};

  bfam_long_real_t **lr;
  lr = bfam_malloc_aligned(num_interp * sizeof(bfam_long_real_t *));

  for (int i = 0; i < num_interp; ++i)
    lr[i] = bfam_malloc_aligned(Nrp * sizeof(bfam_long_real_t));

  bfam_long_real_t *restrict lw;
  lw = bfam_malloc_aligned(Nrp * sizeof(bfam_long_real_t));

  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr[0], lw);

  const bfam_long_real_t half = BFAM_LONG_REAL(0.5);
  for (int n = 0; n < Nrp; ++n)
  {
    lr[1][n] = -half + half * lr[0][n];
    lr[2][n] = half + half * lr[0][n];
  }

  bfam_long_real_t *restrict sub_m_lr, *restrict sub_m_lw;
  sub_m_lr = bfam_malloc_aligned(sub_m_Nrp * sizeof(bfam_long_real_t));
  sub_m_lw = bfam_malloc_aligned(sub_m_Nrp * sizeof(bfam_long_real_t));

  bfam_jacobi_gauss_lobatto_quadrature(0, 0, N_m, sub_m_lr, sub_m_lw);

  bfam_long_real_t *restrict sub_m_V;
  sub_m_V =
      bfam_malloc_aligned(sub_m_Nrp * sub_m_Nrp * sizeof(bfam_long_real_t));

  bfam_jacobi_p_vandermonde(0, 0, N_m, N_m + 1, sub_m_lr, sub_m_V);

  bfam_long_real_t **interpolation =
      bfam_malloc_aligned(num_interp * sizeof(bfam_long_real_t *));

  for (int i = 0; i < num_interp; ++i)
  {
    interpolation[i] =
        bfam_malloc_aligned(Nrp * sub_m_Nrp * sizeof(bfam_long_real_t));

    bfam_jacobi_p_interpolation(0, 0, N_m, Nrp, lr[i], sub_m_V,
                                interpolation[i]);
  }

  bfam_long_real_t **massprojection =
      bfam_malloc_aligned(num_interp * sizeof(bfam_long_real_t *));

  bfam_long_real_t *restrict lV;
  lV = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_long_real_t));
  bfam_jacobi_p_vandermonde(0, 0, N, N + 1, lr[0], lV);

  bfam_long_real_t *restrict mass =
      bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_long_real_t));
  bfam_jacobi_p_mass(0, 0, N, lV, mass);

  for (int i = 0; i < num_interp; ++i)
  {
    massprojection[i] =
        bfam_malloc_aligned(Nrp * sub_m_Nrp * sizeof(bfam_long_real_t));
    for (int n = 0; n < Nrp * sub_m_Nrp; n++)
      massprojection[i][n] = 0;

    bfam_util_mTmmult(sub_m_Nrp, Nrp, Nrp, interpolation[i], Nrp, mass, Nrp,
                      massprojection[i], sub_m_Nrp);
    for (int n = 0; n < Nrp * sub_m_Nrp; n++)
      massprojection[i][n] *= projection_scale[i];
  }

  bfam_long_real_t **projection =
      bfam_malloc_aligned(num_interp * sizeof(bfam_long_real_t *));

  bfam_long_real_t *Vt_MP =
      bfam_malloc_aligned(Nrp * sub_m_Nrp * sizeof(bfam_long_real_t));

  for (int i = 0; i < num_interp; ++i)
  {
    projection[i] =
        bfam_malloc_aligned(Nrp * sub_m_Nrp * sizeof(bfam_long_real_t));
    for (int n = 0; n < Nrp * sub_m_Nrp; n++)
      projection[i][n] = 0;
    for (int n = 0; n < Nrp * sub_m_Nrp; n++)
      Vt_MP[n] = 0;

    bfam_util_mTmmult(sub_m_Nrp, Nrp, sub_m_Nrp, sub_m_V, sub_m_Nrp,
                      massprojection[i], sub_m_Nrp, Vt_MP, sub_m_Nrp);
    bfam_util_mmmult(sub_m_Nrp, Nrp, sub_m_Nrp, sub_m_V, sub_m_Nrp, Vt_MP,
                     sub_m_Nrp, projection[i], sub_m_Nrp);
  }

  glue_m->interpolation =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));

  glue_m->massprojection =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));

  glue_m->projection =
      bfam_malloc_aligned(glue_m->num_interp * sizeof(bfam_real_t *));

  glue_m->exact_mass = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_real_t));
  for (int n = 0; n < Nrp * Nrp; ++n)
    glue_m->exact_mass[n] = (bfam_real_t)mass[n];

  for (int i = 0; i < glue_m->num_interp; ++i)
  {
    if (i == 0 && N_m == N)
    {
      /*
       * Identity interpolation and projection operator
       */
      glue_m->interpolation[i] = NULL;
    }
    else
    {
      glue_m->interpolation[i] =
          bfam_malloc_aligned(Nrp * sub_m_Nrp * sizeof(bfam_real_t));
      for (int n = 0; n < Nrp * sub_m_Nrp; ++n)
        glue_m->interpolation[i][n] = (bfam_real_t)interpolation[i][n];
    }

    glue_m->massprojection[i] =
        bfam_malloc_aligned(Nrp * sub_m_Nrp * sizeof(bfam_real_t));

    for (int n = 0; n < Nrp * sub_m_Nrp; ++n)
      glue_m->massprojection[i][n] = (bfam_real_t)massprojection[i][n];

    glue_m->projection[i] =
        bfam_malloc_aligned(Nrp * sub_m_Nrp * sizeof(bfam_real_t));

    for (int n = 0; n < Nrp * sub_m_Nrp; ++n)
      glue_m->projection[i][n] = (bfam_real_t)projection[i][n];
  }
#endif

  /* STOP */

  glue_p->same_order = (N_p == N) && (N_m == N);
  glue_p->EToEp = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
  glue_p->EToHp = bfam_malloc_aligned(K * sizeof(int8_t));
  glue_p->EToEm = bfam_malloc_aligned(K * sizeof(bfam_locidx_t));
  glue_p->EToFm = bfam_malloc_aligned(K * sizeof(int8_t));
  glue_p->EToHm = bfam_malloc_aligned(K * sizeof(int8_t));
  glue_p->EToOm = bfam_malloc_aligned(K * sizeof(int8_t));
  if (DIM == 1)
    glue_p->num_orient = 2;
  else if (DIM == 2)
    glue_p->num_orient = 8;
  else
    BFAM_ABORT("Cannot handle dim = %d", DIM);

  glue_p->mapOm =
      bfam_malloc_aligned(glue_p->num_orient * sizeof(bfam_locidx_t *));
  if (DIM == 1)
  {
    for (int n = 0; n < glue_p->num_orient; n++)
      glue_p->mapOm[n] = bfam_malloc_aligned(Nrp * sizeof(bfam_locidx_t));
    for (int n = 0; n < Nrp; n++)
    {
      glue_p->mapOm[0][n] = n;
      glue_p->mapOm[1][n] = Nrp - (n + 1);
    }
  }
  else if (DIM == 2)
  {
    for (int n = 0; n < glue_p->num_orient; n++)
      glue_p->mapOm[n] = bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_locidx_t));
    for (int j = 0; j < Nrp; j++)
      for (int i = 0; i < Nrp; i++)
      {
        int ir = Nrp - (i + 1);
        int jr = Nrp - (j + 1);

        glue_p->mapOm[0][i + j * Nrp] = i + j * Nrp;
        glue_p->mapOm[1][i + j * Nrp] = j + i * Nrp;
        glue_p->mapOm[2][i + j * Nrp] = ir + j * Nrp;
        glue_p->mapOm[3][i + j * Nrp] = j + ir * Nrp;
        glue_p->mapOm[4][i + j * Nrp] = jr + i * Nrp;
        glue_p->mapOm[5][i + j * Nrp] = i + jr * Nrp;
        glue_p->mapOm[6][i + j * Nrp] = jr + ir * Nrp;
        glue_p->mapOm[7][i + j * Nrp] = ir + jr * Nrp;
      }
  }
  else
    BFAM_ABORT("Cannot handle dim = %d", DIM);

  qsort(mapping, K, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_send_cmp);

  for (bfam_locidx_t k = 0; k < K; ++k)
    mapping[k].i = k;

  qsort(mapping, K, sizeof(bfam_subdomain_face_map_entry_t),
        bfam_subdomain_face_recv_cmp);

  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    glue_p->EToEm[k] = ktok_m[mapping[k].k];
    glue_p->EToFm[k] = mapping[k].f;
    glue_p->EToHm[k] = mapping[k].h;
    glue_p->EToOm[k] = mapping[k].o;

    glue_p->EToEp[k] = mapping[k].i;
    glue_p->EToHp[k] = mapping[k].nh;
  }

#ifdef BFAM_DEBUG
  for (bfam_locidx_t k = 0; k < K; ++k)
    BFAM_ASSERT(mapping[k].s == glue_m->base.id_s &&
                mapping[k].ns == glue_p->base.id_s);
  BFAM_LDEBUG("    k EToEp EToHp EToEm EToFm EToHm EToOm");
  for (bfam_locidx_t k = 0; k < K; ++k)
  {
    BFAM_LDEBUG("%5d %5d %5d %5d %5d %5d %5d", k, glue_p->EToEp[k],
                glue_p->EToHp[k], glue_p->EToEm[k], glue_p->EToFm[k],
                glue_p->EToHm[k], glue_p->EToOm[k]);
  }
#endif

  for (int i = 0; i < num_interp; ++i)
  {
    // bfam_free_aligned(lr[i]);
    // bfam_free_aligned(interpolation[i]);
    // bfam_free_aligned(massprojection[i]);
    // bfam_free_aligned(projection[i]);
  }

  // bfam_free_aligned(lV);
  // bfam_free_aligned(mass);

  // bfam_free_aligned(lr);
  // bfam_free_aligned(lw);

  // bfam_free_aligned(interpolation);
  // bfam_free_aligned(massprojection);
  // bfam_free_aligned(projection);

  // bfam_free_aligned(sub_m_lr);
  // bfam_free_aligned(sub_m_lw);
  // bfam_free_aligned(sub_m_V);
  // bfam_free_aligned(Vt_MP);
}

static void bfam_subdomain_dgx_geo(
    int N, bfam_locidx_t K, int Np, int ***gmask, const int *restrict Ng,
    const int *restrict Ngp, int num_Vi,
    /*const*/ bfam_long_real_t **restrict xi,
    /*const*/ bfam_long_real_t *restrict Dr, bfam_long_real_t **restrict Jrx,
    bfam_long_real_t *restrict J, bfam_long_real_t **restrict ni,
    bfam_long_real_t *restrict sJ, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_geo");
  const int DIM = inDIM;
#endif
  BFAM_ABORT_IF(num_Vi != DIM && !(J == NULL && ni == NULL && sJ == NULL),
                "[J,ni,sJ] != NULL not implemented for this case");
  BFAM_ABORT_IF_NOT((J == NULL) == (ni == NULL) && (J == NULL) == (sJ == NULL),
                    "[J,ni,sJ] must all be NULL or not NULL");

  const int Nfaces = Ng[0];
  const int Nfp = Ngp[0];

  /*
  BFAM_ASSUME_ALIGNED( x, 32);
  BFAM_ASSUME_ALIGNED( y, 32);
  BFAM_ASSUME_ALIGNED(Dr, 32);
  BFAM_ASSUME_ALIGNED(Jrx, 32);
  BFAM_ASSUME_ALIGNED(Jsx, 32);
  BFAM_ASSUME_ALIGNED(Jry, 32);
  BFAM_ASSUME_ALIGNED(Jsy, 32);
  BFAM_ASSUME_ALIGNED( J, 32);
  */

  BFAM_ASSERT(DIM > 0);
  for (bfam_locidx_t k = 0, vsk = 0, fsk = 0; k < K; ++k)
  {
    if (DIM == 1)
      for (int v = 0; v < num_Vi; v++)
      {
        for (int i = 0; i < N + 1; i++)
          Jrx[v][vsk + i] = 0;
        bfam_util_mvmult(N + 1, N + 1, Dr, N + 1, xi[v] + vsk, Jrx[v] + vsk);
      }
    else if (DIM == 2)
      for (int v = 0; v < num_Vi; v++)
      {
        BFAM_KRON_IXA(N + 1, Dr, xi[v] + vsk, Jrx[0 + 2 * v] + vsk); /* xr */
        BFAM_KRON_AXI(N + 1, Dr, xi[v] + vsk, Jrx[1 + 2 * v] + vsk); /* xs */
      }
    else if (DIM == 3)
      for (int v = 0; v < num_Vi; v++)
      {
        BFAM_KRON_IXIXA(N + 1, Dr, xi[v] + vsk,
                        Jrx[0 + 3 * v] + vsk); /* xi_r */
        BFAM_KRON_IXAXI(N + 1, Dr, xi[v] + vsk,
                        Jrx[1 + 3 * v] + vsk); /* xi_s */
        BFAM_KRON_AXIXI(N + 1, Dr, xi[v] + vsk,
                        Jrx[2 + 3 * v] + vsk); /* xi_t */
      }
    else
      BFAM_ABORT("Cannot handle dim = %d", DIM);

    if (J)
    {
      if (DIM == 1)
      {
        for (int n = 0; n < Np; ++n)
          J[n + vsk] = BFAM_LONG_REAL_ABS(Jrx[0][n + vsk]);
        for (int n = 0; n < Nfp; ++n)
        {
          const bfam_locidx_t fidx0 = fsk + 0 * Nfp + n;
          const bfam_locidx_t fidx1 = fsk + 1 * Nfp + n;

          const bfam_locidx_t vidx0 = vsk + gmask[0][0][n];
          const bfam_locidx_t vidx1 = vsk + gmask[0][1][n];

          /* face 0 */
          ni[0][fidx0] = -Jrx[0][vidx0]; /* -sy */

          /* face 1 */
          ni[0][fidx1] = Jrx[0][vidx1]; /*  sy */

          sJ[fidx0] = BFAM_LONG_REAL_ABS(ni[0][fidx0]);
          sJ[fidx1] = BFAM_LONG_REAL_ABS(ni[0][fidx1]);

          ni[0][fidx0] /= sJ[fidx0];
          ni[0][fidx1] /= sJ[fidx1];
        }
      }
      else if (DIM == 2)
      {
        for (int n = 0; n < Np; ++n)
        {
          bfam_locidx_t idx = n + vsk;

          const bfam_long_real_t xr = Jrx[0][idx];
          const bfam_long_real_t xs = Jrx[1][idx];
          const bfam_long_real_t yr = Jrx[2][idx];
          const bfam_long_real_t ys = Jrx[3][idx];

          /* xr*ys - xs*yr */
          J[idx] = xr * ys - xs * yr;

          /* J*rx = ys */
          Jrx[0][idx] = ys;

          /* J*ry = -xs */
          Jrx[1][idx] = -xs;

          /* J*sx = -yr */
          Jrx[2][idx] = -yr;

          /* J*sy = xr */
          Jrx[3][idx] = xr;
        }
        for (int n = 0; n < Nfp; ++n)
        {
          const bfam_locidx_t fidx0 = fsk + 0 * Nfp + n;
          const bfam_locidx_t fidx1 = fsk + 1 * Nfp + n;
          const bfam_locidx_t fidx2 = fsk + 2 * Nfp + n;
          const bfam_locidx_t fidx3 = fsk + 3 * Nfp + n;

          const bfam_locidx_t vidx0 = vsk + gmask[0][0][n];
          const bfam_locidx_t vidx1 = vsk + gmask[0][1][n];
          const bfam_locidx_t vidx2 = vsk + gmask[0][2][n];
          const bfam_locidx_t vidx3 = vsk + gmask[0][3][n];

          /* rx = 0; ry = 1; sx = 2; sy = 3 */

          /* face 0 */
          ni[0][fidx0] = -Jrx[0][vidx0]; /* -Jrx/sJ */
          ni[1][fidx0] = -Jrx[1][vidx0]; /* -Jry/sJ */

          /* face 1 */
          ni[0][fidx1] = Jrx[0][vidx1]; /*  Jrx/sJ */
          ni[1][fidx1] = Jrx[1][vidx1]; /*  Jry/sJ */

          /* face 2 */
          ni[0][fidx2] = -Jrx[2][vidx2]; /* -Jsx/sJ */
          ni[1][fidx2] = -Jrx[3][vidx2]; /* -Jsy/sJ */

          /* face 3 */
          ni[0][fidx3] = Jrx[2][vidx3]; /*  Jsx/sJ */
          ni[1][fidx3] = Jrx[3][vidx3]; /*  Jsy/sJ */

          sJ[fidx0] = BFAM_LONG_REAL_HYPOT(ni[0][fidx0], ni[1][fidx0]);
          sJ[fidx1] = BFAM_LONG_REAL_HYPOT(ni[0][fidx1], ni[1][fidx1]);
          sJ[fidx2] = BFAM_LONG_REAL_HYPOT(ni[0][fidx2], ni[1][fidx2]);
          sJ[fidx3] = BFAM_LONG_REAL_HYPOT(ni[0][fidx3], ni[1][fidx3]);

          ni[0][fidx0] /= sJ[fidx0];
          ni[1][fidx0] /= sJ[fidx0];
          ni[0][fidx1] /= sJ[fidx1];
          ni[1][fidx1] /= sJ[fidx1];
          ni[0][fidx2] /= sJ[fidx2];
          ni[1][fidx2] /= sJ[fidx2];
          ni[0][fidx3] /= sJ[fidx3];
          ni[1][fidx3] /= sJ[fidx3];
        }
      }
      else if (DIM == 3)
      {
        bfam_long_real_t u1[Np];
        bfam_long_real_t u2[Np];
        bfam_long_real_t u3[Np];

        bfam_long_real_t v1[Np];
        bfam_long_real_t v2[Np];
        bfam_long_real_t v3[Np];

        bfam_long_real_t w1[Np];
        bfam_long_real_t w2[Np];
        bfam_long_real_t w3[Np];

        for (int n = 0; n < Np; ++n)
        {
          bfam_locidx_t idx = n + vsk;

          const bfam_long_real_t xr = Jrx[0][idx];
          const bfam_long_real_t xs = Jrx[1][idx];
          const bfam_long_real_t xt = Jrx[2][idx];
          const bfam_long_real_t yr = Jrx[3][idx];
          const bfam_long_real_t ys = Jrx[4][idx];
          const bfam_long_real_t yt = Jrx[5][idx];
          const bfam_long_real_t zr = Jrx[6][idx];
          const bfam_long_real_t zs = Jrx[7][idx];
          const bfam_long_real_t zt = Jrx[8][idx];

          /*       xr*(ys*zt-yt*zs) - xs*(yr*zt-yt*zr) + xt*(yr*zs-ys*zr) */
          J[idx] = xr * (ys * zt - yt * zs) - xs * (yr * zt - yt * zr) +
                   xt * (yr * zs - ys * zr);

          const bfam_long_real_t x = xi[0][idx];
          const bfam_long_real_t y = xi[1][idx];
          const bfam_long_real_t z = xi[2][idx];

          /* u = z \nabla y - y \nabla z */
          u1[n] = z * yr - y * zr;
          u2[n] = z * ys - y * zs;
          u3[n] = z * yt - y * zt;

          /* v = x \nabla z - z \nabla x */
          v1[n] = x * zr - z * xr;
          v2[n] = x * zs - z * xs;
          v3[n] = x * zt - z * xt;

          /* w = y \nabla x - x \nabla y */
          w1[n] = y * xr - x * yr;
          w2[n] = y * xs - x * ys;
          w3[n] = y * xt - x * yt;

#if 0
          /* J*rx     =  ( ys * zt - yt * zs ) */
          Jrx[0][idx] = (ys * zt - yt * zs);

          /* J*ry     = -( xs * zt - xt * zs ) */
          Jrx[1][idx] = -(xs * zt - xt * zs);

          /* J*rz     =  ( xs * yt - xt * ys ) */
          Jrx[2][idx] = (xs * yt - xt * ys);

          /* J*sx     = -( yr * zt - yt * zr ) */
          Jrx[3][idx] = -(yr * zt - yt * zr);

          /* J*sy     =  ( xr * zt - xt * zr ) */
          Jrx[4][idx] = (xr * zt - xt * zr);

          /* J*sz     = -( xr * yt - xt * yr ) */
          Jrx[5][idx] = -(xr * yt - xt * yr);

          /* J*tx     =  ( yr * zs - ys * zr ) */
          Jrx[6][idx] = (yr * zs - ys * zr);

          /* J*ty     = -( xr * zs - xs * zr ) */
          Jrx[7][idx] = -(xr * zs - xs * zr);

          /* J*tz     =  ( xr * ys - xs * yr ) */
          Jrx[8][idx] = (xr * ys - xs * yr);
#endif
        }

        bfam_long_real_t tmp1[Np];
        bfam_long_real_t tmp2[Np];

        // DR: BFAM_KRON_IXIXA(N + 1, Dr, x, dx);
        // DS: BFAM_KRON_IXAXI(N + 1, Dr, x, dx);
        // DT: BFAM_KRON_AXIXI(N + 1, Dr, x, dx);

        /* <Jrx,Jsx,Jtx> = (1/2) \nabla \times \vec{u} */
        BFAM_KRON_AXIXI(N + 1, Dr, u2, tmp1);
        BFAM_KRON_IXAXI(N + 1, Dr, u3, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[0][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        BFAM_KRON_IXIXA(N + 1, Dr, u3, tmp1);
        BFAM_KRON_AXIXI(N + 1, Dr, u1, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[3][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        BFAM_KRON_IXAXI(N + 1, Dr, u1, tmp1);
        BFAM_KRON_IXIXA(N + 1, Dr, u2, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[6][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        /* <Jry,Jsy,Jty> = (1/2) \nabla \times \vec{v} */
        BFAM_KRON_AXIXI(N + 1, Dr, v2, tmp1);
        BFAM_KRON_IXAXI(N + 1, Dr, v3, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[1][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        BFAM_KRON_IXIXA(N + 1, Dr, v3, tmp1);
        BFAM_KRON_AXIXI(N + 1, Dr, v1, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[4][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        BFAM_KRON_IXAXI(N + 1, Dr, v1, tmp1);
        BFAM_KRON_IXIXA(N + 1, Dr, v2, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[7][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        /* <Jrz,Jsz,Jtz> = (1/2) \nabla \times \vec{w} */
        BFAM_KRON_AXIXI(N + 1, Dr, w2, tmp1);
        BFAM_KRON_IXAXI(N + 1, Dr, w3, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[2][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        BFAM_KRON_IXIXA(N + 1, Dr, w3, tmp1);
        BFAM_KRON_AXIXI(N + 1, Dr, w1, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[5][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        BFAM_KRON_IXAXI(N + 1, Dr, w1, tmp1);
        BFAM_KRON_IXIXA(N + 1, Dr, w2, tmp2);
        for (int n = 0; n < Np; n++)
          Jrx[8][n + vsk] = 0.5 * (tmp1[n] - tmp2[n]);

        for (int n = 0; n < Nfp; ++n)
        {
          const bfam_locidx_t fidx0 = fsk + 0 * Nfp + n;
          const bfam_locidx_t fidx1 = fsk + 1 * Nfp + n;
          const bfam_locidx_t fidx2 = fsk + 2 * Nfp + n;
          const bfam_locidx_t fidx3 = fsk + 3 * Nfp + n;
          const bfam_locidx_t fidx4 = fsk + 4 * Nfp + n;
          const bfam_locidx_t fidx5 = fsk + 5 * Nfp + n;

          const bfam_locidx_t vidx0 = vsk + gmask[0][0][n];
          const bfam_locidx_t vidx1 = vsk + gmask[0][1][n];
          const bfam_locidx_t vidx2 = vsk + gmask[0][2][n];
          const bfam_locidx_t vidx3 = vsk + gmask[0][3][n];
          const bfam_locidx_t vidx4 = vsk + gmask[0][4][n];
          const bfam_locidx_t vidx5 = vsk + gmask[0][5][n];

          /* face 0 */
          ni[0][fidx0] = -Jrx[0][vidx0]; /* -Jrx/sJ */
          ni[1][fidx0] = -Jrx[1][vidx0]; /* -Jry/sJ */
          ni[2][fidx0] = -Jrx[2][vidx0]; /* -Jrz/sJ */

          /* face 1 */
          ni[0][fidx1] = Jrx[0][vidx1]; /*  Jrx/sJ */
          ni[1][fidx1] = Jrx[1][vidx1]; /*  Jry/sJ */
          ni[2][fidx1] = Jrx[2][vidx1]; /*  Jrz/sJ */

          /* face 2 */
          ni[0][fidx2] = -Jrx[3][vidx2]; /* -Jsx/sJ */
          ni[1][fidx2] = -Jrx[4][vidx2]; /* -Jsy/sJ */
          ni[2][fidx2] = -Jrx[5][vidx2]; /* -Jsz/sJ */

          /* face 3 */
          ni[0][fidx3] = Jrx[3][vidx3]; /*  Jsx/sJ */
          ni[1][fidx3] = Jrx[4][vidx3]; /*  Jsy/sJ */
          ni[2][fidx3] = Jrx[5][vidx3]; /*  Jsz/sJ */

          /* face 4 */
          ni[0][fidx4] = -Jrx[6][vidx4]; /* -Jtx/sJ */
          ni[1][fidx4] = -Jrx[7][vidx4]; /* -Jty/sJ */
          ni[2][fidx4] = -Jrx[8][vidx4]; /* -Jtz/sJ */

          /* face 5 */
          ni[0][fidx5] = Jrx[6][vidx5]; /*  Jtx/sJ */
          ni[1][fidx5] = Jrx[7][vidx5]; /*  Jty/sJ */
          ni[2][fidx5] = Jrx[8][vidx5]; /*  Jtz/sJ */

          sJ[fidx0] = BFAM_LONG_REAL_HYPOT(
              ni[0][fidx0], BFAM_LONG_REAL_HYPOT(ni[1][fidx0], ni[2][fidx0]));
          sJ[fidx1] = BFAM_LONG_REAL_HYPOT(
              ni[0][fidx1], BFAM_LONG_REAL_HYPOT(ni[1][fidx1], ni[2][fidx1]));
          sJ[fidx2] = BFAM_LONG_REAL_HYPOT(
              ni[0][fidx2], BFAM_LONG_REAL_HYPOT(ni[1][fidx2], ni[2][fidx2]));
          sJ[fidx3] = BFAM_LONG_REAL_HYPOT(
              ni[0][fidx3], BFAM_LONG_REAL_HYPOT(ni[1][fidx3], ni[2][fidx3]));
          sJ[fidx4] = BFAM_LONG_REAL_HYPOT(
              ni[0][fidx4], BFAM_LONG_REAL_HYPOT(ni[1][fidx4], ni[2][fidx4]));
          sJ[fidx5] = BFAM_LONG_REAL_HYPOT(
              ni[0][fidx5], BFAM_LONG_REAL_HYPOT(ni[1][fidx5], ni[2][fidx5]));

          ni[0][fidx0] /= sJ[fidx0];
          ni[1][fidx0] /= sJ[fidx0];
          ni[2][fidx0] /= sJ[fidx0];
          ni[0][fidx1] /= sJ[fidx1];
          ni[1][fidx1] /= sJ[fidx1];
          ni[2][fidx1] /= sJ[fidx1];
          ni[0][fidx2] /= sJ[fidx2];
          ni[1][fidx2] /= sJ[fidx2];
          ni[2][fidx2] /= sJ[fidx2];
          ni[0][fidx3] /= sJ[fidx3];
          ni[1][fidx3] /= sJ[fidx3];
          ni[2][fidx3] /= sJ[fidx3];
          ni[0][fidx4] /= sJ[fidx4];
          ni[1][fidx4] /= sJ[fidx4];
          ni[2][fidx4] /= sJ[fidx4];
          ni[0][fidx5] /= sJ[fidx5];
          ni[1][fidx5] /= sJ[fidx5];
          ni[2][fidx5] /= sJ[fidx5];
        }
      }
    }

    vsk += Np;
    fsk += Nfaces * Nfp;
  }
}

void BFAM_APPEND_EXPAND(bfam_subdomain_dgx_init_grid_, BFAM_DGX_DIMENSION)(
    bfam_subdomain_dgx_t *subdomain, const int num_Vi,
    const bfam_long_real_t **Vi, const bfam_locidx_t *EToV,
    void (*nodes_transform)(const bfam_locidx_t num_Vi,
                            const bfam_locidx_t num_pnts,
                            bfam_long_real_t **lxi, void *user_args),
    void *user_args, const int inDIM)
{
#ifdef USE_GENERIC_DGX_DIMENSION
  BFAM_WARNING("Using generic bfam_subdomain_dgx_init_grid");
  const int DIM = inDIM;
#endif

  const int *Ng = subdomain->Ng;
  const int *Ngp = subdomain->Ngp;
  const int numg = subdomain->numg;
  const int Np = subdomain->Np;
  const int N = subdomain->N;
  const bfam_locidx_t K = subdomain->K;

  if (DIM > 0)
  {
    const int Nrp = N + 1;
    bfam_long_real_t *lr, *lw;
    lr = bfam_malloc_aligned(Nrp * sizeof(bfam_long_real_t));
    lw = bfam_malloc_aligned(Nrp * sizeof(bfam_long_real_t));

    bfam_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);

    bfam_long_real_t **lxi =
        bfam_malloc_aligned(num_Vi * sizeof(bfam_long_real_t *));

    for (int i = 0; i < num_Vi; i++)
      lxi[i] = bfam_malloc_aligned(K * Np * sizeof(bfam_long_real_t));

    /* Loop over all the elements and set up the grid*/
    int Ncorners = Ng[numg - 1];
    for (bfam_locidx_t k = 0; k < K; ++k)
    {
      const bfam_locidx_t *v = EToV + Ncorners * k;
      bfam_long_real_t w[Ncorners];

      if (DIM == 1)
        for (int n = 0; n < Nrp; ++n)
        {
          int offset = n;
          w[0] = 1 - lr[n];
          w[1] = 1 + lr[n];
          for (int i = 0; i < num_Vi; i++)
          {
            lxi[i][Np * k + offset] = 0;
            for (int c = 0; c < Ncorners; c++)
              lxi[i][Np * k + offset] += w[c] * Vi[i][v[c]];
            lxi[i][Np * k + offset] /= Ncorners;
          }
        }
      else if (DIM == 2)
        for (int n = 0; n < Nrp; ++n)
          for (int m = 0; m < Nrp; ++m)
          {
            int offset = n * Nrp + m;
            w[0] = (1 - lr[m]) * (1 - lr[n]);
            w[1] = (1 + lr[m]) * (1 - lr[n]);
            w[2] = (1 - lr[m]) * (1 + lr[n]);
            w[3] = (1 + lr[m]) * (1 + lr[n]);
            for (int i = 0; i < num_Vi; i++)
            {
              lxi[i][Np * k + offset] = 0;
              for (int c = 0; c < Ncorners; c++)
                lxi[i][Np * k + offset] += w[c] * Vi[i][v[c]];
              lxi[i][Np * k + offset] /= Ncorners;
            }
          }
      else if (DIM == 3)
        for (int n = 0; n < Nrp; ++n)
          for (int m = 0; m < Nrp; ++m)
            for (int l = 0; l < Nrp; ++l)
            {
              int offset = n * Nrp * Nrp + m * Nrp + l;
              w[0] = (1 - lr[l]) * (1 - lr[m]) * (1 - lr[n]);
              w[1] = (1 + lr[l]) * (1 - lr[m]) * (1 - lr[n]);
              w[2] = (1 - lr[l]) * (1 + lr[m]) * (1 - lr[n]);
              w[3] = (1 + lr[l]) * (1 + lr[m]) * (1 - lr[n]);
              w[4] = (1 - lr[l]) * (1 - lr[m]) * (1 + lr[n]);
              w[5] = (1 + lr[l]) * (1 - lr[m]) * (1 + lr[n]);
              w[6] = (1 - lr[l]) * (1 + lr[m]) * (1 + lr[n]);
              w[7] = (1 + lr[l]) * (1 + lr[m]) * (1 + lr[n]);
              for (int i = 0; i < num_Vi; i++)
              {
                lxi[i][Np * k + offset] = 0;
                for (int c = 0; c < Ncorners; c++)
                  lxi[i][Np * k + offset] += w[c] * Vi[i][v[c]];
                lxi[i][Np * k + offset] /= Ncorners;
              }
            }
      else
        BFAM_ABORT("not setup of dim = %d", DIM);
    }

    if (nodes_transform)
      nodes_transform(num_Vi, K * Np, lxi, user_args);

    bfam_long_real_t *restrict lV = subdomain->lV;

    bfam_long_real_t *restrict D =
        bfam_malloc_aligned(Nrp * Nrp * sizeof(bfam_long_real_t));

    bfam_jacobi_p_differentiation(0, 0, N, Nrp, lr, lV, D);

    bfam_long_real_t **lJrx =
        bfam_malloc_aligned(num_Vi * DIM * sizeof(bfam_long_real_t *));
    for (int n = 0; n < num_Vi * DIM; n++)
      lJrx[n] = bfam_malloc_aligned(K * Np * sizeof(bfam_long_real_t));

    bfam_long_real_t *lJ = NULL;
    bfam_long_real_t **lni = NULL;
    bfam_long_real_t *lsJ = NULL;
    if (DIM == num_Vi)
    {
      lJ = bfam_malloc_aligned(K * Np * sizeof(bfam_long_real_t));

      /* Ng[0] = number of faces, Ngp[0] = number of face points */
      lni = bfam_malloc_aligned(num_Vi * sizeof(bfam_long_real_t *));
      for (int n = 0; n < num_Vi; n++)
        lni[n] =
            bfam_malloc_aligned(K * Ng[0] * Ngp[0] * sizeof(bfam_long_real_t));

      lsJ = bfam_malloc_aligned(K * Ng[0] * Ngp[0] * sizeof(bfam_long_real_t));
    }

    bfam_subdomain_dgx_geo(N, K, Np, subdomain->gmask, Ng, Ngp, num_Vi, lxi, D,
                           lJrx, lJ, lni, lsJ, DIM);

    /* store the grid */
    for (int i = 0; i < num_Vi; i++)
    {
      char name[BFAM_BUFSIZ];
      snprintf(name, BFAM_BUFSIZ, "_grid_x%d", i);
      int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
      BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s", name);
      bfam_real_t *restrict xi =
          bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
      for (int n = 0; n < K * Np; ++n)
      {
        xi[n] = (bfam_real_t)lxi[i][n];
      }
    }

    /* store the metric stuff */
    if (lJ)
    {
      for (int d = 0; d < DIM; d++)
      {
        for (int v = 0; v < num_Vi; v++)
        {
          char name[BFAM_BUFSIZ];
          snprintf(name, BFAM_BUFSIZ, "_grid_Jr%dx%d", d, v);
          int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
          BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s", name);
          bfam_real_t *restrict Jrx =
              bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
          for (int n = 0; n < K * Np; ++n)
            Jrx[n] = (bfam_real_t)lJrx[v + d * num_Vi][n];
        }
      }
      {
        char name[] = "_grid_J";
        int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
        BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s", name);
        bfam_real_t *restrict J =
            bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
        for (int n = 0; n < K * Np; ++n)
          J[n] = (bfam_real_t)lJ[n];
      }
      {
        char name[] = "_grid_JI";
        int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
        BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s", name);
        bfam_real_t *restrict JI =
            bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
        for (int n = 0; n < K * Np; ++n)
          JI[n] = (bfam_real_t)(BFAM_LONG_REAL(1.0) / lJ[n]);
      }
      BFAM_ASSERT(lni != NULL && lsJ != NULL);
      for (int v = 0; v < num_Vi; v++)
      {
        char name[BFAM_BUFSIZ];
        snprintf(name, BFAM_BUFSIZ, "_grid_nx%d", v);
        int rval = bfam_subdomain_dgx_field_face_add(&subdomain->base, name);
        BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s", name);
        bfam_real_t *restrict ni =
            bfam_dictionary_get_value_ptr(&subdomain->base.fields_face, name);
        for (int n = 0; n < K * Ng[0] * Ngp[0]; ++n)
          ni[n] = (bfam_real_t)lni[v][n];
      }
      {
        char name[] = "_grid_sJ";
        int rval = bfam_subdomain_dgx_field_face_add(&subdomain->base, name);
        BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s", name);
        bfam_real_t *restrict sJ =
            bfam_dictionary_get_value_ptr(&subdomain->base.fields_face, name);

        for (int n = 0; n < K * Ng[0] * Ngp[0]; ++n)
          sJ[n] = (bfam_real_t)lsJ[n];
      }
    }
    else
    {
      /* In this Jrx really has dx/dr */
      for (int v = 0; v < num_Vi; v++)
      {
        for (int d = 0; d < DIM; d++)
        {
          char name[BFAM_BUFSIZ];
          snprintf(name, BFAM_BUFSIZ, "_grid_x%dr%d", v, d);
          int rval = bfam_subdomain_dgx_field_add(&subdomain->base, name);
          BFAM_ABORT_IF_NOT(rval == 2, "Error adding %s", name);
          bfam_real_t *restrict xr =
              bfam_dictionary_get_value_ptr(&subdomain->base.fields, name);
          for (int n = 0; n < K * Np; ++n)
            xr[n] = (bfam_real_t)lJrx[d + v * DIM][n];
        }
      }
      BFAM_ASSERT(lni == NULL && lsJ == NULL);
    }

    /* free stuff */
    if (lsJ)
      bfam_free_aligned(lsJ);
    if (lni)
    {
      for (int n = 0; n < num_Vi; n++)
        bfam_free_aligned(lni[n]);
      bfam_free_aligned(lni);
    }

    if (lJ)
      bfam_free_aligned(lJ);

    for (int n = 0; n < num_Vi * DIM; n++)
      bfam_free_aligned(lJrx[n]);
    bfam_free_aligned(lJrx);

    bfam_free_aligned(D);

    bfam_free_aligned(lr);
    bfam_free_aligned(lw);

    for (int i = 0; i < num_Vi; i++)
      bfam_free_aligned(lxi[i]);
    bfam_free_aligned(lxi);
  }
}

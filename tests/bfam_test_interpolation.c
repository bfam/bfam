#include <bfam.h>

void
bilinear_test()
{
  bfam_gloidx_t N[3]  = {44,24,28};
  bfam_long_real_t x[(N[0]+1)*(N[1]+1)*(N[2]+1)];
  bfam_long_real_t y[(N[0]+1)*(N[1]+1)*(N[2]+1)];
  bfam_long_real_t z[(N[0]+1)*(N[1]+1)*(N[2]+1)];


  bfam_long_real_t xc[8] = {-0.5, 0.5,-0.3, 0.3,-1.0, 0.0, 1.0, 0.0};
  bfam_long_real_t yc[8] = {-0.1,-0.2, 0.1, 0.2, 0.0,-1.0, 0.0, 1.0};
  bfam_long_real_t zc[8] = {-1.0,-1.0,-1.1,-0.9, 1.0, 1.2, 1.1, 0.7};

  bfam_util_linear_blend(x,NULL,NULL,1,N,NULL,NULL,xc,yc,zc);
  BFAM_ABORT_IF_NOT(fabs(x[N[0]/2]) < 1e-15,"1d interp fail");

  bfam_util_linear_blend(x,y,NULL,2,N,NULL,NULL,xc,yc,zc);
  BFAM_ABORT_IF_NOT(fabs(x[N[0]/2+N[1]/2*(N[0]+1)]) < 1e-15 &&
      fabs(y[N[0]/2+N[1]/2*(N[0]+1)]) < 1e-15, "2-D interp fail");

  bfam_util_linear_blend(x,y,z,3,N,NULL,NULL,xc,yc,zc);
  BFAM_ABORT_IF_NOT(
      fabs(x[N[0]/2+(N[1]/2+N[2]/2*(N[1]+1))*(N[0]+1)]) < 1e-15 &&
      fabs(y[N[0]/2+(N[1]/2+N[2]/2*(N[1]+1))*(N[0]+1)]) < 1e-15 &&
      fabs(z[N[0]/2+(N[1]/2+N[2]/2*(N[1]+1))*(N[0]+1)]) < 1e-15,
      "3-D interp fail");

  bfam_locidx_t Nl[3]  = { 6, 2,10};
  bfam_gloidx_t gx[3]  = {19,11, 9};
  bfam_long_real_t xl[(Nl[0]+1)*(Nl[1]+1)*(Nl[2]+1)];
  bfam_long_real_t yl[(Nl[0]+1)*(Nl[1]+1)*(Nl[2]+1)];
  bfam_long_real_t zl[(Nl[0]+1)*(Nl[1]+1)*(Nl[2]+1)];

  bfam_util_linear_blend(xl,NULL,NULL,1,N,Nl,gx,xc,yc,zc);
  BFAM_ABORT_IF_NOT(fabs(xl[Nl[0]/2]) < 1e-15,"1d interp fail");

  bfam_util_linear_blend(xl,yl,NULL,2,N,Nl,gx,xc,yc,zc);
  BFAM_ABORT_IF_NOT(fabs(xl[Nl[0]/2+Nl[1]/2*(Nl[0]+1)]) < 1e-15 &&
      fabs(yl[Nl[0]/2+Nl[1]/2*(Nl[0]+1)]) < 1e-15, "2-D interp fail");


  bfam_util_linear_blend(xl,yl,zl,2,N,Nl,gx,xc,yc,zc);
  BFAM_ABORT_IF_NOT(
      fabs(xl[Nl[0]/2+(Nl[1]/2)*(Nl[0]+1)]) < 1e-15 &&
      fabs(yl[Nl[0]/2+(Nl[1]/2)*(Nl[0]+1)]) < 1e-15 &&
      fabs(1+zl[Nl[0]/2+(Nl[1]/2)*(Nl[0]+1)]) < 1e-15,
      "3-D interp fail");

  bfam_util_linear_blend(xl,yl,zl,3,N,Nl,gx,xc,yc,zc);
  BFAM_ABORT_IF_NOT(
      fabs(xl[Nl[0]/2+(Nl[1]/2+Nl[2]/2*(Nl[1]+1))*(Nl[0]+1)]) < 1e-15 &&
      fabs(yl[Nl[0]/2+(Nl[1]/2+Nl[2]/2*(Nl[1]+1))*(Nl[0]+1)]) < 1e-15 &&
      fabs(zl[Nl[0]/2+(Nl[1]/2+Nl[2]/2*(Nl[1]+1))*(Nl[0]+1)]) < 1e-15,
      "3-D interp fail");
}

int
main (int argc, char *argv[])
{
  bilinear_test();
  return EXIT_SUCCESS;
}

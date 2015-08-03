/*
 * This is adapted from the OCCA2 add vectors example.
 */

#include <bfam.h>
#include <occa_c.h>

int main(int argc, char **argv)
{
  int rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  BFAM_MPI_CHECK(MPI_Init(&argc,&argv));
  BFAM_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  bfam_log_init(rank, stdout, BFAM_LL_DEFAULT);
  bfam_signal_handler_set();

  occaPrintAvailableDevices();

  int entries = 5;
  int i;

  float *a  = (float*) malloc(entries*sizeof(float));
  float *b  = (float*) malloc(entries*sizeof(float));
  float *ab = (float*) malloc(entries*sizeof(float));

  for(i = 0; i < entries; ++i){
    a[i]  = (float)i;
    b[i]  = (float)(1 - i);
    ab[i] = 0;
  }

  occaDevice device;
  occaKernel addVectors;
  occaMemory o_a, o_b, o_ab;

  //---[ Device setup with string flags ]-------------------
  // const char *deviceInfo   = "mode = Serial";
  const char *deviceInfo   = "mode = OpenMP  , schedule = compact, chunk = 10";
  // const char *deviceInfo   = "mode = OpenCL  , platformID = 0, deviceID = 0";
  // const char *deviceInfo   = "mode = CUDA    , deviceID = 0";
  // const char *deviceInfo   = "mode = COI     , deviceID = 0";

  device = occaCreateDevice(deviceInfo);

  o_a  = occaDeviceMalloc(device, entries*sizeof(float), NULL);
  o_b  = occaDeviceMalloc(device, entries*sizeof(float), NULL);
  o_ab = occaDeviceMalloc(device, entries*sizeof(float), NULL);

  occaKernelInfo info = occaCreateKernelInfo();
  occaKernelInfoAddDefine(info, "DIMENSION", occaInt(10));

  addVectors = occaDeviceBuildKernel(device,
                                     "addVectors.occa", "addVectors",
                                     info);

  int dims = 1;
  occaDim itemsPerGroup, groups;

  itemsPerGroup.x = 2;
  groups.x        = (entries + itemsPerGroup.x - 1)/itemsPerGroup.x;

  occaKernelSetWorkingDims(addVectors,
                           dims, itemsPerGroup, groups);

  occaCopyPtrToMem(o_a, a, entries*sizeof(float), 0);
  occaCopyPtrToMem(o_b, b, occaAutoSize, occaNoOffset);

  occaKernelRun(addVectors,
               occaInt(entries), o_a, o_b, o_ab);

  occaCopyMemToPtr(ab, o_ab, occaAutoSize, occaNoOffset);

  for(i = 0; i < 5; ++i)
    printf("%d = %f\n", i, ab[i]);

  for(i = 0; i < entries; ++i){
    if(ab[i] != (a[i] + b[i]))
      exit(1);
  }

  free(a);
  free(b);
  free(ab);

  occaKernelFree(addVectors);
  occaMemoryFree(o_a);
  occaMemoryFree(o_b);
  occaMemoryFree(o_ab);
  occaDeviceFree(device);

  BFAM_MPI_CHECK(MPI_Finalize());

  return EXIT_SUCCESS;
}


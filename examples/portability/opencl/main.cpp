#include <iostream>

#include "occa.hpp"

int main(int argc, char **argv){
  int entries = 5;

  //---[ Init OpenCL ]------------------
  cl_int error;

  cl_platform_id clPlatformID = occa::cl::platformID(0);
  cl_device_id clDeviceID     = occa::cl::deviceID(0,0);

  cl_context clContext = clCreateContext(NULL,
                                         1, &clDeviceID,
                                         NULL, NULL, &error);
  OCCA_CL_CHECK("Device: Creating Context", error);

  cl_command_queue clStream = clCreateCommandQueue(clContext,
                                                   clDeviceID,
                                                   CL_QUEUE_PROFILING_ENABLE, &error);
  OCCA_CL_CHECK("Device: createStream", error);

  cl_mem cl_a = clCreateBuffer(clContext,
                               CL_MEM_READ_WRITE,
                               entries*sizeof(float), NULL, &error);

  cl_mem cl_b = clCreateBuffer(clContext,
                               CL_MEM_READ_WRITE,
                               entries*sizeof(float), NULL, &error);

  cl_mem cl_ab = clCreateBuffer(clContext,
                                CL_MEM_READ_WRITE,
                                entries*sizeof(float), NULL, &error);
  //====================================

  float *a  = new float[entries];
  float *b  = new float[entries];
  float *ab = new float[entries];

  occa::device device = occa::cl::wrapDevice(clPlatformID,
                                             clDeviceID,
                                             clContext);

  occa::stream stream = device.wrapStream(&clStream);
  device.setStream(stream);

  occa::kernel addVectors;
  occa::memory o_a, o_b, o_ab;

  for(int i = 0; i < entries; ++i){
    a[i]  = i;
    b[i]  = 1 - i;
    ab[i] = 0;
  }

  o_a  = device.wrapMemory(&cl_a , entries*sizeof(float));
  o_b  = device.wrapMemory(&cl_b , entries*sizeof(float));
  o_ab = device.wrapMemory(&cl_ab, entries*sizeof(float));

  addVectors = device.buildKernelFromSource("addVectors.occa",
                                            "addVectors");

  int dims = 1;
  int itemsPerGroup(2);
  int groups((entries + itemsPerGroup - 1)/itemsPerGroup);

  addVectors.setWorkingDims(dims, itemsPerGroup, groups);

  o_a.copyFrom(a);
  o_b.copyFrom(b);

  occa::initTimer(device);

  occa::tic("addVectors");

  addVectors(entries, o_a, o_b, o_ab);

  double elapsedTime = occa::toc("addVectors", addVectors);

  o_ab.copyTo(ab);

  std::cout << "Elapsed time = " << elapsedTime << " s" << std::endl;

  occa::printTimer();

  for(int i = 0; i < 5; ++i)
    std::cout << i << ": " << ab[i] << '\n';

  addVectors.free();
  o_a.free();
  o_b.free();
  o_ab.free();

  delete [] a;
  delete [] b;
  delete [] ab;

  return 0;
}

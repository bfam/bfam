#include <iostream>

#include "occa.hpp"

int main(int argc, char **argv){
  int entries = 8;

  float *a  = new float[entries];
  float *b  = new float[entries];
  float *ab = new float[entries];

  for(int i = 0; i < entries; ++i){
    a[i]  = i;
    b[i]  = 1 - i;
    ab[i] = 0;
  }

  std::string mode = "CUDA";
  int platformID = 0;
  int deviceID   = 0;

  occa::device device;
  occa::kernel addVectors;
  occa::memory o_a, o_b, o_ab;

  occa::stream streamA, streamB;

  device.setup(mode, platformID, deviceID);

  streamA = device.getStream();
  streamB = device.createStream();

  o_a  = device.malloc(entries*sizeof(float));
  o_b  = device.malloc(entries*sizeof(float));
  o_ab = device.malloc(entries*sizeof(float));

  addVectors = device.buildKernelFromSource("addVectors.occa",
                                            "addVectors");

  int dims = 1;
  int itemsPerGroup(2);
  int groups(4);

  addVectors.setWorkingDims(dims, itemsPerGroup, groups);

  o_a.copyFrom(a);
  o_b.copyFrom(b);

  device.setStream(streamA);
  addVectors(entries, 0, o_a, o_b, o_ab);

  device.setStream(streamB);
  addVectors(entries, 4, o_a, o_b, o_ab);

  o_ab.copyTo(ab);

  for(int i = 0; i < entries; ++i)
    std::cout << i << ": " << ab[i] << '\n';

  delete [] a;
  delete [] b;
  delete [] ab;

  addVectors.free();
  o_a.free();
  o_b.free();
  o_ab.free();
  device.free();
}

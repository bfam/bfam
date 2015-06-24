#include <iostream>

#include "occa.hpp"

int main(int argc, char **argv){
  int entries = 10000; // Not divisible
  int p_Nred = 256;
  int reducedEntries = (entries + p_Nred - 1)/p_Nred;

  float *a    = new float[entries];
  float *aRed = new float[reducedEntries];

  float trueRed = 0;

  for(int i = 0; i < entries; ++i){
    a[i]     = 1;
    trueRed += a[i];
  }

  for(int i = 0; i < reducedEntries; ++i)
    aRed[i] = 0;

  occa::device device;
  occa::kernel reduction;
  occa::memory o_a, o_aRed;

  occa::kernelInfo reductionInfo;

  device.setup("mode = Serial");

  o_a  = device.malloc(entries*sizeof(float));
  o_aRed = device.malloc(reducedEntries*sizeof(float));

  reductionInfo.addDefine("p_Nred", p_Nred);

#if 1
  reduction = device.buildKernelFromSource("reduction.okl",
                                           "reduction",
                                           reductionInfo);
#else
  reduction = device.buildKernelFromSource("reduction.cu",
                                           "reduction",
                                           reductionInfo);

  size_t dims     = 1;
  occa::dim inner(p_Nred);
  occa::dim outer((entries + p_Nred - 1) / p_Nred);

  reduction.setWorkingDims(dims, inner, outer);
#endif

  o_a.copyFrom(a);

  occa::initTimer(device);

  occa::tic("reduction");

  reduction(entries, o_a, o_aRed);

  double elapsedTime = occa::toc("reduction", reduction);

  o_aRed.copyTo(aRed);

  std::cout << "Elapsed time = " << elapsedTime << " s" << std::endl;

  occa::printTimer();

  for(int i = 1; i < reducedEntries; ++i)
    aRed[0] += aRed[i];

  if(aRed[0] != trueRed){
    std::cout << "aRed[0] = " << aRed[0] << '\n'
              << "trueRed = " << trueRed << '\n';

    std::cout << "Reduction failed\n";
    throw 1;
  }
  else
    std::cout << "Reduction(a) = " << aRed[0] << '\n';

  delete [] a;
  delete [] aRed;

  reduction.free();
  o_a.free();
  o_aRed.free();

  device.free();

  return 0;
}

#ifndef TENSORHELPER_H
#define TENSORHELPER_H
#include <TensorLite.h>
namespace tensor{

  EXPORTLIB void CombineImage(const vector<string>& infilenames, string& outfilename, int step=5);
//  struct GPUBuffer                                  // Optimized GPU versions
//  {   // Data allocations are very expensive on GPU. Use a buffer to solve: allocate once reuse later.
//      gpu::GpuMat gI1, gI2, gs, t1,t2;
//  };

}
#endif // TENSORHELPER_H

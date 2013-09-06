#ifndef TENSORHELPER_H
#define TENSORHELPER_H
#include <TensorLite.h>
namespace tensor{

  EXPORTLIB void CombineImage(const vector<string>& infilenames, string& outfilename, int step=5);


}
#endif // TENSORHELPER_H

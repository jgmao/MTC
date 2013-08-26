#include <TensorLite.h>
#include <MetricData.h>
#include <dirent.h>
#ifndef GRANULATE_H
#define GRANULATE_H

namespace metric {

#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif


class Granulate
{
public:
  Granulate(){}
  ~Granulate(){}
  bool readFiles(string path, string searchExt);
  bool generateGrid();
  vector<string> filenames;
  vector<Tensor<uchar,1> > images;

static const vector<Size3> blkSizes;
static const vector<int> sz_idx ;

};

}


#endif // GRANULATE_H

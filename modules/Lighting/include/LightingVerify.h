#include <string>

#include <sys/stat.h>

#include <Metric.h>
#include <fstream>
#include <PoissonSolver.h>
#include <Lighting.h>


#ifndef LIGHTINGVERIFY_H
#define LIGHTINGVERIFY_H
using namespace std;
using namespace metric;
using namespace tensor;
namespace lighting{
class LightingVerify
{
public:
  LightingVerify(void);
  ~LightingVerify(void);
  void readfile(string filename, string ext="tiff");
  vector<double>& doPLC(int blksize, int bdsize, int fsize);//blksize size of block <= imagesize, bdsize is the size of boundary, fsize is size to find food average
  vector<double>& doAdaptivePLC(int blksize, int bdsize, int fsize, double thrd);
  vector<Tensor<double,1>>& estimateLight(vector<Tensor<double,1>>& images,vector<Tensor<double,1>>& rst, int blksize, int bdsize, int fsize);
  void estimateLight(Tensor<double,1>& im,Tensor<double,1>& rst, int blksize, int bdsize, int fsize);
  void saveResults(string folder="rst/");
  string fpath;
  int blksize,bdsize,fsize;
  vector<FootType> feet;
protected:
   vector<string> fnames;
   vector<Tensor<double,1>> orignalWithBorder;
   vector<Tensor<double,1>> originalImages;
   vector<Tensor<double,1>> correctedImages;
   vector<Tensor<double,1>> lightingOrg;
   vector<Tensor<double,1>> lightingRst;
public:
   vector<double> MSEs; 
   Lighting lt;
};
}
#endif

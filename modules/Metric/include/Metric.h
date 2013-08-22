#include <string>
#include <math.h>
#include <TensorLite.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>
#include <regex>
#ifndef DEC_FLAG
#define DEC_FLAG
#define DEC_FLAG_HARD 0
#define DEC_FLAG_SOFT 1
#endif
using namespace std;
using namespace boost::math;
static const double pi = 3.14159265358979323846;

class Metric
{
public:
  Metric();
public:
  bool readFiles(string path, string searchPattern, string searchExt);
  double computeMetric(const Tensor<double,1>& im1, const Tensor<double,1>& im2);
  double computeMetricLSE(const Tensor<double,1>& im1, const Tensor<double,1>& im2);
  int computeStats(string path, string searchPattern, string searchExt);
  int computeFeatures(string searchPattern, string searchExt);
  int computeFeatures(const vector<vector<Tensor<double,2>>>& stats);
  int computeStats(const Tensor<double,1>& im1, const Tensor<double,1>& im2);
  Mat getFeature(const Mat& f, int subband, int feature);
  Mat& estimateML(const Mat& f, int subband, int feature,vector<vector<Mat_<double>>>& gamma, vector<vector<double>>& lambda);
  Mat& computeLLR(const Mat& f, Mat& llr, int subband, int feature, vector<vector<Mat_<double>>>& gamma0, vector<vector<Mat_<double>>>& gamma1, vector<vector<double>>& lambda0, vector<vector<double>>& lambda1);
  Mat& makeDec(const Mat& LLR, Mat& dec, double thred = 0, int flag=0); 
  void trainMetric(string path, string searchPattern, string searchExt);
  void trainMetirc(string path,string scorefilepath="");
  void loadParams(void);
  float predict(const Mat& f);
  
  Size3 subwinSize,subwinStep;
  string searchPattern, searchExt,searchPath;
  int featureNum;
  int Nor,Nsc;
  bool subsample;
  bool changeWin;
protected:
  vector<vector<Mat_<double> > > gamma1,gamma0;
  vector<vector<double> > lambda1,lambda0;
  vector<Tensor<double,2>>  mu; //each entry corresponds to a subband, the order of subband is high to low, 0 to pi and then low pass high pass
  vector<Tensor<double,2>> sigma;
  vector<Tensor<double,2>> r01;
  vector<Tensor<double,2>> r10;
  vector<string> filenames;
  map<string,int> clusternames;
  int clustercount;
  int coeffcount0, coeffcount1;
  int sameCount,diffCount;
  int pairNum; //number of pairs
  int coeffNum; //number of total coeffs in a feature, that is #subbands * #feature_in_band * featureNum for all subbands
  vector<vector<Tensor<double,2> > > stats;
  vector<int> coeff_in_band;
  vector<vector<Tensor<double,2>>> sameClusterFeature;
  vector<vector<Tensor<double,2>>> diffClusterFeature;
  //f1 is the same cluster feature
  //f0 is diff cluster feature
  //vector<Mat> f1,f0;
  Mat f1,f0;
  Mat LLR0,LLR1,dec0,dec1,llr,dec;
  Mat label;
  Mat f;
  SVM s;
  SVMParams params;
  Mat weights;
  void saveGamma(string filename,const vector<vector<Mat_<double>>>& gamma);
  void loadGamma(string filename, vector<vector<Mat_<double>>>& gamma);
  void saveLambda(string filename, const vector<vector<double>>& lambda);
  void loadLambda(string filename, vector<vector<double>>& lambda);
  bool loadClassifier(string filename);
};

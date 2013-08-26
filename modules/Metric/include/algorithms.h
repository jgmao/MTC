#include <Steerable.h>
#include <TensorLite.h>
#include <LRI.h>
#include <MetricData.h>
#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
using namespace tensor;
#ifndef ALGORITHMS_H
#define ALGORITHMS_H
namespace metric{
  
  EXPORTLIB double ComputeMSE(const Mat& tsA, const Mat& tsB) ; // if lambda = 0, then it is normal MSE, otherwise it is MSE with constraint
  EXPORTLIB bool ComputeAIM(const Mat& tsA, const Mat& tsB, Scalar thrd); // adaptive interpolation metric	for PQI
  EXPORTLIB Tensor<uchar,1> CompareElement(const Mat& A, const Mat& B, int flag = CMP_LE);//@
  EXPORTLIB Tensor<uchar,1> CompareElement(const Mat& A, Scalar thrd, int flag=CMP_LE);//@// the return type is definited
  EXPORTLIB double ComputeSAD(const Mat& tsA, const Mat& tsB) ;
  EXPORTLIB double ComputePSNR(const Mat& tsA, const Mat& tsB);
  EXPORTLIB double ComputeLRI(const Mat& tsA, const Mat& tsB);
  EXPORTLIB double ComputeSTSIM2(const Mat& tsA, const Mat& tsB, 
                             const Size3& subWinSize = Size3(8,8,1),
                             const Size3& subWinStep = Size3(1,1,1),
                             int nLevel=3, int nDir=4,
                             bool downsample = false,
                             FilterBoundary boundary_cut = FilterBoundary::FILTER_BOUND_HALF,
                             FeaturePoolType stsim2_pool_type = FeaturePoolType::FEATURE_POOL_MIN,
                             MetricModifier stsim2_modifer = MetricModifier::STSIM2_BASELINE,
                             bool debug=false);
  EXPORTLIB double ComputeSTSIM3_LSE(const Mat& tsA, const Mat& tsB, 
                             const Size3& subWinSize = Size3(8,8,1),
                             const Size3& subWinStep = Size3(1,1,1),
                             int nLevel=3, int nDir=4,
                             bool downsample = false,
                             FilterBoundary boundary_cut = FilterBoundary::FILTER_BOUND_HALF,
                             bool debug=false);
  EXPORTLIB vector<Tensor<double,2> > ComputeStatistics(const Mat& ts,
                             const Size3& subWinSize = Size3(8,8,1),
                             const Size3& subWinStep = Size3(1,1,1),
                             bool subsample = false,
                             bool changeWin=false,
                             int nLevel = 3, int nDir=4,
                             FilterBoundary boundary_cut=FilterBoundary::FILTER_BOUND_FULL,
                             bool compute00=false
                             ) ;
  // this need SVM Metric implementation EXPORTLIB double ComputeSVMMetric(const Mat& tsA, const Mat& tsB, const Size3& subWinSize, const Size3& subWinStep);
/*! legacy
  EXPORTLIB double Compare(const Mat& ts, int criteria,
                         double param1=3, double param2=4, int param3 = 1,
                         int param4 = (int)Stsim2PoolType::STSIM2_POOL_MIN,
                         int param5 = (int)MetricModifier::STSIM2_BASELINE, cv::Mat& param6 = cv::Mat());
*/
  EXPORTLIB double Compare(const Mat& tsA, const Mat& tsB, CompareCriteria criteria, const Size3& subWinSize, const Size3& subWinStep,
                         Printable_t param1=Printable_t(0),
                         Printable_t param2=Printable_t(0), Printable_t param3=Printable_t(0),
                         Printable_t param4=Printable_t(0), Printable_t param5=Printable_t(0),
                         Printable_t param6=Printable_t(0),
                         bool debug=false);

  EXPORTLIB cv::Mat EstimateVarForMahalanobis(const Mat& ts, Size3 wsize, Size3 stepsize);
  EXPORTLIB double ComputeMahalanobis(const Mat& tsA, const Mat& tsB, Size3 subWinSize, Size3 subWinStep, const Mat& iMcovar, FilterBoundary boundary_cut = FilterBoundary::FILTER_BOUND_HALF);
  EXPORTLIB Tensor<double,2> ComputeRho(const Mat& im11, const Mat& im12, const Size3& subWinSize, const Size3& subWinStep);
  EXPORTLIB Tensor<double,1> ComputeCrossTerm(const Mat& im11, const Mat& im12, const Mat& im21, const Mat& im22, const Size3& subWinSize, const Size3& subWinStep);
}
#endif

#include <Steerable.h>
#include <TensorLite.h>

#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
using namespace tensor;
namespace metric{

  EXPORTLIB double ComputeMSE(const Mat& tsA, const Mat& tsB) ; // if lambda = 0, then it is normal MSE, otherwise it is MSE with constraint
        EXPORTLIB bool ComputeAIM(const Mat& tsA, const Mat& tsB, Scalar thrd); // adaptive interpolation metric	for PQI
        EXPORTLIB Tensor<uchar,1> CompareElement(const Mat& A, const Mat& B, int flag = CMP_LE);//@
        EXPORTLIB Tensor<uchar,1> CompareElement(const Mat& A, Scalar thrd, int flag=CMP_LE);//@// the return type is definited
//	EXPORTLIB double ComputeSAD(const Tensor<T,cn> ts) const;
//	EXPORTLIB double ComputePSNR(const Tensor<T,cn> ts) const;
//	EXPORTLIB double ComputeLRI(const Tensor<T,cn> ts) const;
//	EXPORTLIB double ComputeSSIM(const Tensor<T,cn> ts, const Size3& subWinSize = Size3(8,8,1), const Size3& subWinStep = Size3(1,1,1),int nLevel=3, int nDir=4, int boundary_cut = FILTER_BOUND_HALF, int ststi2_pool_type = STSIM2_POOL_MIN, int stsim2_modifer = STSIM2_BASELINE) const;
//	EXPORTLIB vector<Tensor<double,2> > ComputeStatistics(const Size3& subWinSize = Size3(8,8,1), const Size3& subWinStep = Size3(1,1,1), bool subsample = false, int nLevel = 3, int nDir=4, bool changeWin=false) const;
// 	EXPORTLIB double ComputeSVMMetric(const Tensor<T,cn> ts, const Size3& subWinSize, const Size3& subWinStep) const;
//	EXPORTLIB double Compare(const Tensor<T,cn> ts, int criteria, double param1=3, double param2=4, int param3 = 1, int param4 = STSIM2_POOL_MIN, int param5 = STSIM2_BASELINE, cv::Mat& param6 = cv::Mat()) const;
// 	EXPORTLIB double Compare(const Tensor<T,cn> ts, int criteria, Printable_t param1=Printable_t(0), Printable_t param2=Printable_t(0), Printable_t param3=Printable_t(0), Printable_t param4=Printable_t(0), Printable_t param5=Printable_t(0), Printable_t param6=Printable_t(0)) const;
// 	EXPORTLIB cv::Mat EstimateVarForMahalanobis(Size3 wsize, Size3 stepsize);
// 	EXPORTLIB double ComputeMahalanobis(const Tensor<T,cn> ts, Size3 subWinSize, Size3 subWinStep, const Mat& iMcovar) const;

}

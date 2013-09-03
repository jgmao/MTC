#ifndef __POISSON_SOLVER_H__
#define __POISSON_SOLVER_H__
#include <map>
#include <opencv2/opencv.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
#include <Eigen/Core>
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>
#include "MyLib.h"
#include "CommonDef.h"
//#include <Eigen/UmfPackSupport>
//#include <Eigen/SuperLUSupport>
#include <vector>
#include <Eigen/SparseCholesky>


namespace lighting{
using namespace std;
using namespace tensor;
class FootType
{
public:
  int x; // x-row
  int y; // y-col
  double value;
  FootType(double v, int x, int y){this->value=v; this->x = x; this->y =y;}
  FootType(double v){this->value=v;}
  FootType(const pair<Point3i,FootItem>& f){this->value=f.second.value;this->x = f.first.x; this->y = f.first.y;}
  FootType(){};
  ~FootType(){}

};
class PoissonSolver
{
  cv::Mat _cand, _tar, _mask, _ref;
  double foot,upfoot,leftfoot;
  vector<FootType> feet;
  cv::Rect mask_roi1;
  cv::Mat mask1;
  cv::Mat dst1;
  cv::Mat tar1;
  cv::Mat cand1;
  cv::Mat laplacian_tar,laplacian_cand; //
  cv::Mat NI;
  cv::Mat ref1;
  int ch;
  double footThred;
  std::map<int,int> mp;
 
  template <typename T>
  bool buildLCMatrix(Eigen::SparseMatrix<T> &A, Eigen::Matrix<T, Eigen::Dynamic, 1> &b, 
                   Eigen::Matrix<T, Eigen::Dynamic, 1> &u);//in PLC only left and up boundary conditions are known
  template <typename T>
  bool buildPostLCMatrix(Eigen::SparseMatrix<T> &A, Eigen::Matrix<T, Eigen::Dynamic, 1> &b, 
                   Eigen::Matrix<T, Eigen::Dynamic, 1> &u); // in post PLC, all the boundary conditions are known
  template <typename T>
  bool buildStichingMatrix(Eigen::SparseMatrix<T> &A, Eigen::Matrix<T, Eigen::Dynamic, 1> &b, 
                   Eigen::Matrix<T, Eigen::Dynamic, 1> &u);
  template <typename T> 
  bool solve(const Eigen::SparseMatrix<T> &A, const Eigen::Matrix<T, Eigen::Dynamic, 1> &b,
             Eigen::Matrix<T, Eigen::Dynamic, 1> &u);
  template <typename T>
  bool copyResult(Eigen::Matrix<T, Eigen::Dynamic, 1> &u);

  void computeLaplacian(cv::Mat& _dst, int offx=-1, int offy=-1);
  vector<double> interpFeet(vector<FootType>& feet, Point spos, Point epos);//spos start point, epos end point, the interpolate line is between two points
public:
  EXPORTLIB PoissonSolver();
  EXPORTLIB PoissonSolver(const cv::Mat &cand, const cv::Mat &target, const cv::Mat &mask, const cv::Mat& ref,  const double& foot, const double& upfoot=-1, const double& leftfoot = -1 );
  EXPORTLIB PoissonSolver(const cv::Mat &cand, const cv::Mat &target, const cv::Mat &mask, const cv::Mat& ref,  const vector<FootType>& feet);
  EXPORTLIB PoissonSolver(const cv::Mat &cand, const cv::Mat &target, const cv::Mat &mask);
  EXPORTLIB ~PoissonSolver() {};
  EXPORTLIB bool setImages(const cv::Mat &src, const cv::Mat &target, const cv::Mat &mask);
  EXPORTLIB void copyTo(PoissonSolver &b) const;
  EXPORTLIB PoissonSolver clone() const;
  EXPORTLIB bool poissonLightCorrection(cv::Mat &_dst, int offx=-1, int offy=-1);
  EXPORTLIB bool PostPLC(cv::Mat &_dst);
  EXPORTLIB void gradientStiching(cv::Mat & _dst, int offx=-1, int offy=-1);//, cv::Mat & target, const cv::Mat & mask);
  // bool smoothComplete(cv::Mat &dst); // implemented easily with zero boundary conditions.
  EXPORTLIB void test(void);
  EXPORTLIB void poissonLightChange(cv::Mat &_dst);
  EXPORTLIB void addFoot(const FootType& f);
};
}
#endif /*  __POISSON_SOLVER_H__ */

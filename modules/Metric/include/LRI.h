#include <opencv2/opencv.hpp>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace cv;
using namespace std;
#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
namespace metric{
class LRI
{
public:
  EXPORTLIB LRI(void);
  EXPORTLIB ~LRI(void);
private:
  double var(const Mat& im);
  int* lri8k(double image[][9], double th);
  double computeLRI(const Mat& im1, const Mat& im2);
  int lbp81(double image[][3]);
  double computeLBP(const Mat& im1, const Mat& im2);
  double* grad(const Mat& im);
  double sedest(const Mat& im1, const Mat& im2);
  double lp(const Mat& im1, const Mat& im2);
  int* a;
  double* gr;
public:
  EXPORTLIB double computeNewMetric(const Mat& im1, const Mat& im2);
};
}

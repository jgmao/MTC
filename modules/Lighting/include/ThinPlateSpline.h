#pragma once
#include <opencv2/opencv.hpp>
#ifndef THINPLATESPLINE_H
#define THINPLATESPLINE_H
using namespace cv;
using namespace std;
namespace lighting{
class ThinPlateSpline
{
protected:
  Mat data;
  int m,n;
  double p;
  int datatype;

  Mat g;
  Mat& buildRegulation(double p);
public:
    Mat B;
  Mat& solve(double p);
  Mat& solve();
  void load(Mat& data);
public:
  ThinPlateSpline(Mat& data);
  ThinPlateSpline(void);
  ~ThinPlateSpline(void);
};
}
#endif

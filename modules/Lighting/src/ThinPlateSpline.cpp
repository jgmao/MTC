#include "ThinPlateSpline.h"
//#include "mylib.h"
namespace lighting {

ThinPlateSpline::ThinPlateSpline(void)
{
  this->m=0;
  this->n=0;

}


ThinPlateSpline::~ThinPlateSpline(void)
{



}


ThinPlateSpline::ThinPlateSpline(Mat& data)
{
  load(data);
}

void ThinPlateSpline::load(Mat& data)
{
  datatype = data.type();
  if (data.type()!=CV_64F)
    data.convertTo(this->data,CV_64F);
  else
    this->data = data;
  this->n = data.size().height;
  this->m = data.size().width;
}

Mat& ThinPlateSpline::buildRegulation(double p)
{
  this->B = Mat::zeros(this->m,this->n,CV_64F);
  double temp = 0;
  for (int i=0; i<this->m; i++)
    for (int j=0; j<this->n; j++)
    {
      temp = (2-cos(3.1416*double(i)/double(m))-cos(3.1416*double(j)/double(n)));
      temp = ( p + 4*(1-p)*temp*temp);
      B.at<double>(i,j) =1.0/temp;
    }
  return this->B;
}

Mat& ThinPlateSpline::solve(double p)
{
  this->p = p;
  return solve();
}

Mat& ThinPlateSpline::solve()
{
  buildRegulation(this->p);
  //mylib::DisplayMat(B,"1overB",true);
  Mat temp1,temp2;
  cv::dct(this->data,temp1);
  //mylib::DisplayMat(temp1,"dct",true);
  temp1 = temp1*p;
  //mylib::DisplayMat(temp1,"dctxp",true);
  cv::multiply(temp1,this->B,temp2);
  //mylib::DisplayMat(temp2,"dctxpoverB",true);
  cv::idct(temp2,this->g,DFT_SCALE);
  //mylib::DisplayMat(this->g, "idct",true);
  if (datatype!=CV_64F)
    this->g.convertTo(g,datatype);
  return this->g;
}
}

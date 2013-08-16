#include <Tensor.h>
using namespace tensor;

int main(void)
{
  Tensor<double,1> ts("/home/guoxin/Projects/MTC/data/baboon.pgm");
  ts.Display();
  Tensor<double,2> cts = ts.ToComplex();
  Tensor<double,2> f = cts.DFT();
  f.Print();
  f = f.Abs();
  Mat temp = f.GetFrame(0);
  vector<Mat> tv;
  cv::split(temp,tv);
  Tensor<double,1>(tv[0]).Display(); 
  return 0;
}

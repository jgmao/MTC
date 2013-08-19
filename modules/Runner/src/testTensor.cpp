#include <Tensor.h>
using namespace tensor;

int main(void)
{
  Tensor<double,1> ts("/home/guoxin/Projects/MTC/data/texture1.png");
  //ts.Display();
  Tensor<double,2> cts = ts.ToComplex();
  gpu::Stream stream;
  BufferGPU gbuf;
  cts.stream = &stream;
  cts.gbuf = &gbuf;
  double time = (double)getTickCount();
  Tensor<double,2> f;
  for (int i=0; i<1000; i++)
  {
    f = cts.DFT();
    cout<<i<<endl;
  }
  f.stream = &stream;
  f.gbuf = &gbuf;
  f = f.IDFT();

  time = 1000*((double)getTickCount() - time)/getTickFrequency();
  cout<<"time = " <<time<<"ms"<<endl;
  Mat temp = f.GetFrame(0);
  vector<Mat> tv;
  cv::split(temp,tv);
  Tensor<double,1>(tv[0]).Display(); 
  return 0;
}




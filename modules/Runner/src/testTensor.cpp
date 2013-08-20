#include <TensorLite.h>
#include <Steerable.h>
#include <algorithms.h>
using namespace tensor;
using namespace metric;
int main(void)
{
  Tensor<double,1> ts("/home/guoxin/Projects/MTC/data/texture1.png");
  //ts.Display();
  Tensor<double,2> cts = ts.ToComplex();
  double time = (double)getTickCount();
  Tensor<double,2> f;
  for (int i=0; i<1000; i++)
  {
    f = cts.DFT();
    cout<<i<<endl;
  }
  f = f.IDFT();

  time = 1000*((double)getTickCount() - time)/getTickFrequency();
  cout<<"time = " <<time<<"ms"<<endl;
  Mat temp = f.GetFrame(0);
  vector<Mat> tv;
  cv::split(temp,tv);
  Tensor<double,1>(tv[0]).Display(); 
  Cube c(0,0,0,10,10,1);
  cout<<c.area()<<endl;
  Steerable sp;
  sp.buildSCFpyr(f,3,4,1,false);
  vector<Tensor<double,2> >& pyr = sp.getSpaceDomainPyr();
  cv::split(pyr[12],tv);
  Tensor<double,1>(tv[0]).Display(); 
  Tensor<double,1> ts2("/home/guoxin/Projects/MTC/data/texture2.png");
  cout<<ComputeMSE(ts,ts2)<<endl;
  return 0;
}




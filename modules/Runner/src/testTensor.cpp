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
  ts2.Print("t2",true);

  //test LocalMean
  Mat flatker = Mat::ones(Size(16,16),CV_64F)/256;
  //verify local mean
  BufferGPU gbuf;
  Size3 subWin(16,16,1);
  auto mu1 = ts.LocalMean(flatker,subWin);
  auto mu2 = ts.LocalMeanGPU(flatker,gbuf,subWin);
  (mu1-mu2).Print();
  mu1.Print();
  auto var1 = ts.LocalVariance(mu1,flatker,subWin);
  auto var2 = ts.LocalVarianceGPU(mu2,flatker,gbuf,subWin);
  var1.Print();
  var2.Print();
  (var1-var2).Print();
//  mylib::DisplayMat(flatker);
  int TIMES=10000;
     time = (double)getTickCount();
  for (int i=0; i< TIMES; i++)
    {
  Tensor<double,1> mu=ts.LocalMean(flatker,Size3(16,16,1));
  ts.LocalVariance(mu,flatker,Size3(16,16,1));
    }
  time = 1000*((double)getTickCount() - time)/getTickFrequency();
     time /= TIMES;
 cout << "Time of CPU (averaged for " << TIMES << " runs): " << time << " milliseconds."<<endl;


time = (double)getTickCount();
for (int i=0; i< TIMES; i++)
{
Tensor<double,1> mu=ts.LocalMeanGPU(flatker,gbuf,Size3(16,16,1));
ts.LocalVarianceGPU(mu,flatker,gbuf,Size3(16,16,1));
}
time = 1000*((double)getTickCount() - time)/getTickFrequency();
 time /= TIMES;
cout << "Time of GPU (averaged for " << TIMES << " runs): " << time << " milliseconds."<<endl;



  return 0;



}




#include <algorithms.h>
namespace metric
{


template <class T, size_t cn> double ComputeMSE(const Tensor<T,cn>& tsA, const Tensor<T,cn>& tsB) 
{
	double distance=0;
	CV_DbgAssert(tsA.channels() == cn);
	CV_DbgAssert(tsA.size()==tsB.size());
	Scalar temp; 
	Mat v;
	Mat v1,v2;
  	Mat v1_dft, v2_dft, v_dft;
  //int pnum = 4;
  //vector<thread> threads;
  //vector<Vec<T,cn>> sums(pnum);
	for (int z = 0; z < tsA.size().depth; z++)
	{
		//cv::absdiff(GetFrame(z),ts[z],v);
   /* for (int p=0; p<pnum; p++)
    {
      threads.push_back(thread([&](int p, int z, int H, int W){
        Vec<T,cn> tt=0;
        for (int i=p; i<H; i+=pnum)
          for (int j=p; j<W; j+=pnum)
          {
            tt =  this->operator()(i,j,z)-ts(i,j,z);
            sums[p] += mylib::VecMul<T,cn>(tt,tt);
          }
      },p,z,ts.size().height,ts.size().width));
    }*/
    tsA.GetFrameRef(z,v2);
    tsB.GetFrameRef(z,v2);
    v = v1-v2;
    cv::pow(v,2.0,v);
    //cv::multiply(v1,v2,v);
    //cv::dft(v1,v1_dft);
    //cv::dft(v2,v2_dft);
   /* for (int p=0; p<pnum;p++)
    {
      threads[p].join();
      temp+=sums[p];
    }*/
    //temp[0]+=cv::norm(v1,v2,NORM_L1);
		temp +=cv::sum(v);
	}
	for (int i =0; i < cn; i++)
		distance+= temp[i];
	return distance;///tsSize.volumn()/cn;
}


template  ComputeMSE<double,1>(const Tensor<double,1>& tsA, const Tensor<double,1>& tsB);

}

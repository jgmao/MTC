#include <algorithms.h>
namespace metric
{


	double ComputeMSE(const Mat& tsA, const Mat& tsB)
	{
	  double distance=0;
	  CV_Assert(tsA.channels() == tsB.channels());
	  CV_Assert(tsA.size()==tsB.size());
	  Scalar temp=0;
	  Mat v;
	  v = tsA - tsB;
	  cv::pow(v,2.0,v);
	  temp +=cv::sum(v);
	  for (int i =0; i < tsA.channels(); i++)
	    distance+= temp[i];
	  return distance/double(tsA.size().area())/double(tsA.channels());///tsSize.volumn()/cn;
	}



}

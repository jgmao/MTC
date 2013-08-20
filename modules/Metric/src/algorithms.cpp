#include <algorithms.h>
namespace metric
{

  //! compute Mean Squre error between two Mats
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


  //! compare Mat with another Mat
  Tensor<uchar,1>  CompareElement(const Mat& A, const Mat& B, int flag) {
    CV_Assert(A.size()==B.size());
    CV_Assert(A.channels()==B.channels());
    cv::Mat tempMat;
                //mylib::DisplayMat(A);
    cv::compare(A,B,tempMat,flag);
                //mylib::DisplayMat(tempMat);
    return tempMat;
  }

  //! compare matrix with a scalar
  Tensor<uchar,1> CompareElement(const Mat& A, Scalar thrd, int flag) {
    Mat B(A.size(),A.type(),thrd);
    return CompareElement(A,B,flag);
  }

  //! compute Adaptive Interpolate Metric of D. Neuhoff
  bool ComputeAIM(const Mat& tsA, const Mat& tsB, Scalar thrd)
  {
    Mat temp;
    cv::absdiff(tsA,tsB,temp);
    Tensor<uchar,1> tempMap = CompareElement(temp, thrd);
    Tensor<uchar,1> extMap = tempMap.ExtendBoundary(Size3(1,1,0),255);
    //extMap.Print();
    //extMap.Print(true);
    for (int z = 0 ; z < tempMap.size().depth; z++)
      {
        for (int x = 1; x < tempMap.size().height+1; x++)
          for (int y=1; y < tempMap.size().width+1; y++)
            {
              if (extMap(x,y,z)[0]== 0) //test 1 failed
                {
                  if (extMap(x+1,y,z)[0]/255 + extMap(x,y+1,z)[0]/255
                      + extMap(x-1,y,z)[0]/255 + extMap(x,y-1,z)[0]/255 < 3)
                    return false;
                }
            }
      }
    return true;
  }


}

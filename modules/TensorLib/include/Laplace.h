#ifndef LAPLACE_H
#define LAPLACE_H
#include <opencv2/opencv.hpp>
#include <boost/math/special_functions/bessel.hpp>
class Laplace
{
public:
  Laplace()
  {
    this->mu =Mat();
    this->gamma  =Mat();
    this->history = Mat();
    lambda = 0;
  }

  ~Laplace(){}

  Mat mu;
  Mat gamma;
  Mat igamma;
  double lambda;
  void addData(Mat d)
  {
    //make sure d is 1x2
    history.push_back(d.clone());
  }
  void erase(void)
  {
    history = history.rowRange(1,history.size().height-1);
    //history.erase(history.begin());
  }

  uint size(void)
  {
    return history.size().height;
  }
  void mle1D(void)
  {
    Mat history1D = history.reshape(0,size()*2);
    Mat var;
    cv::meanStdDev(history1D,mu,var);
    Mat ret;
    cv::reduce(cv::abs(history1D-mu.at<double>(0,0)),ret,0,CV_REDUCE_SUM);
    //cout<<ret<<endl;
    lambda = ret.at<double>(0,0)/history1D.size().height;
  }

  void mle(void)
  {
    Mat covar;
   // cout<<size()<<endl;
    cv::calcCovarMatrix(history,covar,mu,cv::COVAR_NORMAL|cv::COVAR_ROWS,history.type());
    lambda = std::sqrt(determinant(covar));
    gamma = covar/lambda;
    igamma = gamma.inv();
    //cout<<gamma<<endl;
    //cout<<igamma<<endl;
    //for (int i=0; i< N; i++)
    //{
    //  mu += history[i];
    //}
    //mu/=double(N);
    //for (int i=0; i< N; i++)
    //{
    //  b += std::abs(history[i]-mu);
    //}
    //b/=double(N);
  }

  void progressEst(double d)
  {
 //   mu*=double(N);
 //   mu+=d;
 //   N++;
 //   mu/=double(N);
 //   b=0;
 //   history.push_back(d);
 //   for(double& x : history)
 //   {
 //     b += std::abs(x-mu);
 //   }
 //   b/=double(N);
 }
  void clear(void)
  {
    this->mu =Mat();
    this->gamma  =Mat();
    this->igamma = Mat();
    this->history = Mat();
    lambda = 0;
  }

  double pdf(Mat d)
  {
    //cout<<mu<<endl;
    //cout<<gamma<<endl;
    //cout<<igamma<<endl;
      double q =cv::Mahalanobis(d,mu,igamma);
      //cout<<q<<endl;
      double K = boost::math::cyl_bessel_k<double,double>(0,q);
      //cout<<K<<endl;
      return K/lambda/3.1415926;
   //return 0.5*exp(-std::abs(d-mu)/b)/b;
  }
  double pdf1D(double d)
  {
    //d must 1x1 matrix
    return std::exp(-std::abs(d-mu.at<double>(0.0))/lambda)/2/lambda;
  }

  double cdfinv(double p)
  {
      return 0;
//    p = p-0.5;
//    double s = (p<0)?-1:1;
 //   return mu - b*s*std::log(1 - 2*abs(p));
  }

protected:

  //int memory;
  Mat history;
};

#endif // LAPLACE_H

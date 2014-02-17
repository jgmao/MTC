#ifndef GAUSSIAN_H
#define GAUSSIAN_H
class Gaussian
{
public:
  Gaussian()
  {
    this->mu =Mat();
    this->sigma  =Mat();
    this->history = Mat();
  }

  ~Gaussian(){}

  Mat mu;
  Mat sigma;
  Mat isigma;
  double detSigma;
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
    cv::meanStdDev(history1D,mu,sigma);
    detSigma = sigma.at<double>(0,0);
    //not it is std, not var
    //var = cv::pow(var,2.0);
  }

  void mle(void)
  {
    cv::calcCovarMatrix(history,sigma,mu,cv::COVAR_NORMAL|cv::COVAR_ROWS,history.type());
    isigma = sigma.inv();
    detSigma = 2*3.1415926*sqrt(cv::determinant(sigma));
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
    this->sigma  =Mat();
    this->isigma = Mat();
    this->history = Mat();
  }

  double pdf(Mat d)
  {
      double q =cv::Mahalanobis(d,mu,isigma);
      return std::exp(-0.5*q)/detSigma;
  }
  double pdf1D(double d)
  {
      return std::exp(-(d-mu.at<double>(0,0))*(d-mu.at<double>(0,0))/2/sigma.at<double>(0,0)/sigma.at<double>(0,0))/sigma.at<double>(0,0)/std::sqrt(2*3.1415926);
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
#endif // GAUSSIAN_H

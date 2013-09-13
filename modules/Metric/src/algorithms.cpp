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

  double ComputeSAD(const Mat& tsA, const Mat& tsB) {
    double distance=0;
    CV_DbgAssert(tsA.channels() == tsB.channels());
    CV_DbgAssert(tsA.size()==tsB.size());
    Scalar temp=0;
    Mat v;
    cv::absdiff(tsA,tsB,v);
    temp +=cv::sum(v);
    for (int i =0; i < tsA.channels(); i++)
      distance+= temp[i];
    return distance/double(tsA.size().area())/double(tsA.channels());
  }

  double ComputePSNR(const Mat& tsA, const Mat& tsB) {
    double distance = ComputeMSE(tsA,tsB);
    //distance/= double(ts.size().area());//do it since ComputeMSE does not normalize
    return 10*log10(std::pow(255.0,2.0)/distance);
  }
  //!gj01072013_2 add LRI based metric to computeLRI()
  double ComputeLRI(const Mat& tsA, const Mat& tsB)
  {
    LRI lri;
    double rst=0;
    rst=lri.computeNewMetric(tsA,tsB);
    return rst;
  }


  vector<Tensor<double,2> > ComputeStatistics(const Mat& ts, const Size3& subWinSize, const Size3& subWinStep, bool subsample, bool changeWin,int nLevel, int nDir,FilterBoundary boundary_cut, FeaturePoolType pooltype, bool compute00 )
  {
    typedef vector<Tensor<double,2> > vT;
    Steerable sp;
    Tensor<double,2> tc = Tensor<double,1>(ts).ToComplex();
    if (boundary_cut == FilterBoundary::FILTER_BOUND_EXTEND) //zero padding and then cut center
    {
      tc=tc.ExtendBoundary(tc.size()/2);
    }
    sp.buildSCFpyr(tc,nLevel,nDir,1,subsample);
    //! instead of use Space domain, use freq domain
    //vT& pyr = sp.getSpaceDomainPyr();
    vT& pyr = sp.getPyr();
    cout<<"subband size : "<<pyr.size()<<endl;
    //cout<<"size is: "<<ts.size()<<endl;
    //cout<<"print extended size: ";
    //tc.size().Print();
    if ( boundary_cut == FilterBoundary::FILTER_BOUND_HALF||boundary_cut == FilterBoundary::FILTER_BOUND_EXTEND)//modify Dec 27 2011, not cut half but cut half + boundary
      {
        for (unsigned int i = 0; i< pyr.size(); i++)
          {
            int hh = pyr[i].size().height/2;
            int ww = pyr[i].size().width/2;
            int dd = pyr[i].size().depth;
            Cube roi(pyr[i].size().height/4, pyr[i].size().width/4, 0, hh,ww,dd);
            pyr[i] = pyr[i](roi);//no size changed
          }
      }
    else if ( boundary_cut == FilterBoundary::FILTER_BOUND_VALID)
      {
        for (unsigned int i = 0; i< pyr.size(); i++)
          {
            int hh = pyr[i].size().height/2 + pyr[i].size().height/4;
            int ww = pyr[i].size().width/2  + pyr[i].size().width/4;
            int dd = pyr[i].size().depth/2  + pyr[i].size().depth/4;
            if (dd<1)
              dd=1;
            Cube roi(pyr[i].size().height/8, pyr[i].size().width/8, 0, hh,ww,dd); // temporaly use this /8
            pyr[i] = pyr[i](roi);//no size changed
          }
      }


    int bands = pyr.size();
    //    Size3 initSize = pyr[0].size();
    Size3 sz = pyr[0].size();
    double weight = subWinSize.height*subWinSize.width;
    double scalar = weight/(weight-1);
    //          cv::Mat temp;
    int index=0;
    Size3 subWinStepLv = subWinStep;
    Size3 subWinSizeLv = subWinSize;
    int lvl=0;
    vT mu(bands);
    vT sigma2(bands);
    vT rho01(bands);
    vT rho10(bands);
    int crossbandNum = nDir*(nLevel-1) + nLevel* (mylib::combination(nDir,2));
    vT rho00 (crossbandNum);
    vT statistics(4*bands);
    for (index = 0; index < (int)pyr.size(); index++)
      {
       // cout<<"===========doing idx "<<index<<endl;
        //pyr[index].Print("pry.txt",true); //checked, the pyr are correct
        sz = pyr[index].size();
        //cout<<"subwinsizeLv ";
        //subWinSizeLv.Print();
        if (subWinSizeLv.height>pyr[index].size().height||subWinSizeLv.width>pyr[index].size().width)
          {
            lvl = lvl+1;
            subWinSizeLv = pyr[index].size();
            subWinStepLv = pyr[index].size();
          }
        else if (index== int(pyr.size()-1)) //highpass
          {
            lvl = 0;
            subWinSizeLv = subWinSize;
            subWinStepLv = subWinStep;
          }
        else
          lvl = index/nDir;

        if (changeWin) //change window size and not lowpass band
          {
            subWinSizeLv = Size3(subWinSize.height/(1<<lvl),subWinSize.width/(1<<lvl),subWinSize.depth);
            subWinStepLv = Size3(subWinStep.height/(1<<lvl),subWinStep.width/(1<<lvl),subWinStep.depth);
          }
        weight = subWinSizeLv.height*subWinSizeLv.width;
        scalar = weight/(weight-1);
        //subWinSizeLv.Print();
        cv::Mat flatKel1=cv::Mat(subWinSizeLv,pyr[index][0].type()-((pyr[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1.0/weight));
        cv::Mat flatKel2 = cv::Mat(subWinSizeLv,pyr[index][0].type()-((pyr[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1.0/(weight-1)));
        //mylib::DisplayMat(flatKel1);
        // pyr[index].Print("pyr",true);
        mu[index] = pyr[index].LocalMean(flatKel1,subWinStepLv);
        //check mu
        // mu[index].Print();
        sigma2[index] = pyr[index].LocalVariance(mu[index]*scalar,flatKel2,subWinStepLv).Sqrt();//take the sqrt to get std
        //is rho01 rho10 a good choice?????? YES 20130722
        rho01[index] = ComputeRho(pyr[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth)),
            pyr[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth)),
            subWinSizeLv-Size3(0,1,0), subWinStepLv);
        rho10[index] = ComputeRho(pyr[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth)),
            pyr[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth)),
            subWinSizeLv-Size3(1,0,0), subWinStepLv);
        statistics[index]=mu[index].Clone();
        statistics[index+bands]=sigma2[index].Clone();
        statistics[index+2*bands] = rho01[index].Clone();
        statistics[index+3*bands] = rho10[index].Clone();

      }
    if (compute00){
        index = 0;
        for (int dr = 0; dr < nDir; dr++)
          for (int lv = 0; lv < nLevel-1; lv++)
            {
              rho00[index] = ComputeRho( pyr[ lv*nDir + dr], pyr[ (lv+1)*nDir + dr],
                  subWinSizeLv,subWinStepLv);
              //cout<<index<<endl;
              index++;
            }


        for (int lv=0; lv< nLevel; lv++)
          for (int dr = 0; dr < nDir-1; dr++)
            for ( int p = dr+1; p < nDir; p++)
              {
                rho00[index] = ComputeRho(pyr[lv*nDir + dr],pyr[lv*nDir + p],
                    subWinSizeLv, subWinStepLv);
               // cout<<index<<endl;
                index++;
              }
        for (index=0; index<crossbandNum; index++)
          {
            statistics.push_back(rho00[index].Clone());
          }
    }

    cout<<"size of feature for window size: ";
    subWinSize.Print();
    cout<<"is "<<statistics[0].size()<<endl;
    for (unsigned int i=0; i<statistics.size(); i++)
    {
        //statistics[i].Print();
        if (pooltype==FeaturePoolType::FEATURE_POOL_AVE)
          statistics[i]=Tensor<double,2>(Size3(1,1,1),statistics[i].Mean());
        else if (pooltype==FeaturePoolType::FEATURE_POOL_MIN)
          statistics[i]=Tensor<double,2>(Size3(1,1,1),statistics[i].Abs().Min());

        //statistics[i].Print();
    }
    return statistics;
  }


  Tensor<double,2> ComputeRho(const Mat& im11, const Mat& im12, const Size3& subWinSize, const Size3& subWinStep)
  {

    if (im11.size()!= im12.size())
      {
#if (CV_MINOR_VERSION >5)
        CV_Error(Error::StsUnmatchedSizes,"im11,im12 must in same size");
#else
        CV_Error(CV_StsUnmatchedSizes,"im11,im12 must in same size");
#endif
      }

    const double C = 0.0001;
    //im11.Print();
    //im12.Print();
    //im21.Print();
    //im22.Print();
    Tensor<double,2> A(im11);
    Tensor<double,2> B(im12);
    Tensor<double,2> mu11,mu12;
    Tensor<double,2> sigma_cross;
    Tensor<double,2> rho;
    Tensor<double,2> sigma11,sigma12;
    //          Mat flatKel = Mat(Size(subWinSize.width,subWinSize.height),CV_64F,Scalar(1/double(subWinSize.width*subWinSize.height)));

    mu11 = A.LocalMean(subWinSize,subWinStep);      //   mu11.Print();
    mu12 = B.LocalMean(subWinSize,subWinStep);//  mu12.Print();
    sigma11 = A.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
    sigma12 = B.LocalVariance(mu12,subWinSize,subWinStep);// sigma12.Print();
    sigma_cross = ((A * B.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
    rho = (sigma_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
    return rho;
  }

  double ComputeSTSIM3_LSE(const Mat& tsA, const Mat& tsB, const Size3& subWinSize, const Size3& subWinStep, int nLevel, int nDir,bool downsample, FilterBoundary boundary_cut, bool debug)
  {
    double rst=0;

    fstream debugfile;
    if (debug)
      {
        //    string PID = boost::lexical_cast<string>(_getpid());

        debugfile.open("/home/guoxin/Projects/MTC/temp/ssim_terms.txt",ios::app);
        debugfile.precision(3);
        //debugfile.width(10);
        debugfile<<"=========inter-band========="<<endl;
        debugfile<<left<<"Band"<<"\t"<<"L"<<"\t"<<"C"<<"\t"<<"C01"<<"\t"<<"C10"<<"\t"<<"Pool"<<endl;
      }

    //Metric mc;
    int coeffNum=0;
    vector<int> coeff_in_band;
    vector<vector<Tensor<double,2>>> stats = vector<vector<Tensor<double,2>>>(2);
    stats[0] = ComputeStatistics(tsA, subWinSize,subWinStep,downsample,false,nLevel,nDir,boundary_cut);
    stats[1] = ComputeStatistics(tsB, subWinSize,subWinStep,downsample,false,nLevel,nDir,boundary_cut);
    for (auto t : stats[0])
      {
        coeffNum += t.size().area();
        coeff_in_band.push_back( t.size().area() );
      }
    Mat f = Mat(coeff_in_band[0],(nLevel*nDir+2)*4*2,CV_64F);
    Rect roi;
    for (unsigned int k=0; k<stats[0].size(); k++)
      {
        roi = Rect(k*2,0,2,coeff_in_band[k]);
        Mat temp=   (stats[0][k] - stats[1][k]).GetFrameRef(0).clone();
        if (k/(nLevel*nDir+2)==0)//means max 255
          temp = temp/255;
        else if (k/(nLevel+nDir+2)==1)//variance max is given by max((mu-min_pix)*(max-pix - mu)) about 128^2
          temp = temp/16384;
        vector<Mat> tt;
        cv::split(temp,tt);
        //abs
        tt[0]=cv::abs(tt[0]);
        tt[1]=cv::abs(tt[1]);
        //satruation
        //cv::threshold(tt[0],tt[0],1,1,cv::THRESH_TRUNC);
        //cv::threshold(tt[1],tt[1],1,1,cv::THRESH_TRUNC);
        cv::merge(tt,temp);
        //mylib::DisplayMat(temp);
        //inverse
        if (k/(nLevel*nDir+2)==0||k/(nLevel*nDir+2)==1)
          temp = Scalar(1,1) - temp;
        //mylib::DisplayMat(temp);
        //save to feature vector
        temp.reshape(1,stats[0][k].size().area())
            .copyTo(f(roi));
        /*  (stats[0][k] - stats[1][k]).Abs()
          .GetFrameRef(0).clone()
          .reshape(1,stats[0][k].size().area())
          .copyTo(f(roi));*/
        //mylib::DisplayMat(f(roi));
      }
    Mat A(1,coeffNum*2,CV_64F);
    Mat temp(coeffNum,2,CV_64F);
    int copypos =0;
    for (unsigned int j=0; j<coeff_in_band.size();j++)
      {
        Rect roi(j*2,0,2,coeff_in_band[j]);
        Rect toroi(0,copypos,2,coeff_in_band[j]);
        f(roi).copyTo(temp(toroi));
        copypos+=coeff_in_band[j];
      }
    //mylib::DisplayMat(temp);
    Mat stsim2_lse_weight = Mat(A.cols,A.rows,CV_64F,(void*)&STSIM2_LSE_WEIGHT);
    //mylib::DisplayMat(stsim2_lse_weight,"weight",true);
    //mylib::DisplayMat(A,"A",true);
    temp.reshape(1,1).copyTo(A.row(0));
    Mat scorerst = A*stsim2_lse_weight;
    // mylib::DisplayMat(scorerst);
    rst = scorerst.at<double>(0,0)/10;
    return rst;
  }


  double ComputeSTSIM2(const Mat& tsA, const Mat& tsB,
                       const Size3& subWinSize,
                       const Size3& subWinStep,
                       int nLevel, int nDir,
                       bool downsample,
                       FilterBoundary boundary_cut ,
                       FeaturePoolType stsim2_pool_type ,
                       MetricModifier stsim2_modifer, bool debug)
  {
    CV_Assert(tsA.size()==tsB.size());
    CV_Assert(tsA.channels()==tsB.channels());
    fstream debugfile;
    if (debug)
      {
        //    string PID = boost::lexical_cast<string>(_getpid());

        debugfile.open("/home/guoxin/Projects/MTC/temp/ssim_terms.txt",ios::app|ios::out);
        debugfile.precision(3);
        //debugfile.width(10);
        debugfile<<"=========inter-band========="<<endl;
        debugfile<<left<<"Band"<<"\t"<<"L"<<"\t"<<"C"<<"\t"<<"C01"<<"\t"<<"C10"<<"\t"<<"Pool"<<endl;
      }

    double rst=0;
    const double C0 = 0.001;
    const double C1 = 0.001;
    //const double C2 = 0.001;a
    //const double C3 = 0.001;
    Steerable spA, spB;
    Tensor<double,2> A;
    Tensor<double,2> B;
    if (tsA.channels()==1)
      {
        A = Tensor<double,1>(tsA).ToComplex();
        B = Tensor<double,1>(tsB).ToComplex();
      }
    else
      {
        CV_Assert(tsA.channels()==2);
        A = Tensor<double,2>(tsA);
        B = Tensor<double,2>(tsB);
      }
    //A.Print("A");
    //B.Print("B");
    spA.buildSCFpyr(A,nLevel,nDir,1,downsample);//A = this is orgExt
    spB.buildSCFpyr(B,nLevel,nDir,1,downsample);//B = ts is candExt
    vector<Tensor<double,2>>& pyrA = spA.getSpaceDomainPyr();
    vector<Tensor<double,2>>& pyrB = spB.getSpaceDomainPyr();
    
    if ( boundary_cut == FilterBoundary::FILTER_BOUND_HALF)//modify Dec 27 2011, not cut half but cut half + boundary
      {
        for (unsigned int i = 0; i< pyrA.size(); i++)
          {
            //pyrA[i].Print();
            //dec 27 2011, extend the block by size of block size + bounday (below and right)
            //int hh = subWinSize.height*(int)ceil(float(pyrA[i].size().height)/2.0/float(subWinSize.height));
            //int ww = subWinSize.width*(int)ceil(float(pyrA[i].size().width)/2.0/float(subWinSize.width));
            //int dd = subWinSize.depth*(int)ceil(float(pyrA[i].size().depth)/2.0/float(subWinSize.depth));
            //just half, dec 28 2012
            int hh = pyrA[i].size().height/2;
            int ww = pyrA[i].size().width/2;
            int dd = pyrA[i].size().depth;
            Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, hh,ww,dd);
            //Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, pyrA[i].size().height/2, pyrA[i].size().width/2, pyrA[i].size().depth);
            pyrA[i] = pyrA[i](roi);//no size changed
            //pyrA[i].Print();
            pyrB[i] = pyrB[i](roi);
            //pyrB[i].Print();
          }
      }
    else if ( boundary_cut == FilterBoundary::FILTER_BOUND_VALID)
      {
        for (unsigned int i = 0; i< pyrA.size(); i++)
          {
            int hh = pyrA[i].size().height/2 + pyrA[i].size().height/4;
            int ww = pyrA[i].size().width/2  + pyrA[i].size().width/4;
            int dd = pyrA[i].size().depth/2  + pyrA[i].size().depth/4;
            if (dd<1)
              dd=1;
            //just half, dec 28 2012
            //int hh = pyrA[i].size().height/2;
            //int ww = pyrA[i].size().width/2;
            //int dd = pyrA[i].size().depth;
            Cube roi(pyrA[i].size().height/8, pyrA[i].size().width/8, 0, hh,ww,dd); // temporaly use this /8
            //Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, pyrA[i].size().height/2, pyrA[i].size().width/2, pyrA[i].size().depth);
            pyrA[i] = pyrA[i](roi);//no size changed
            //pyrA[i].Print();
            pyrB[i] = pyrB[i](roi);
            //pyrB[i].Print();
          }
      }
    //pyrA[0].size().Print();
    //pyrB[0].size().Print();
    vector<Tensor<double,2>> mu_A(pyrA.size());
    vector<Tensor<double,2>> mu_B(pyrB.size());
    vector<Tensor<double,2>> sigma2_A(pyrA.size());
    vector<Tensor<double,2>> sigma2_B(pyrA.size());
    vector<Tensor<double,1>> L(pyrA.size());
    vector<Tensor<double,1>> C(pyrB.size());
    vector<Tensor<double,1>> C01(pyrA.size());
    vector<Tensor<double,1>> C10(pyrB.size());
    int crossbandNum = nDir*(nLevel-1) + nLevel* (mylib::combination(nDir,2));
    vector<Tensor<double,1>> C00 (crossbandNum);
    int index = 0;
    Size3 initSize = pyrA[0].size();
    Size3 sz = pyrA[0].size();
    //double weight;// = sz.height*sz.width;
    cv::Mat temp;
    Size3 subWinStepLv = subWinStep;//include PLC boundary
    Size3 subWinSizeLv = subWinSize;//include PLC boundary

  //  cout<<subWinStepLv<<endl;
    //if (boundary_cut == FILTER_BOUND_VALID)
    //{
    //  subWinSizeLv = subWinSizeLv + subWinSizeLv/2; //2*overlap + size
    //  subWinStepLv = subWinStepLv + subWinStepLv/2;
    //}
    //Tensor<double,1> rstMat(pyrA[0].size()/subWinStep - subWinSize+Size3(1,1,1));
    Size3 tempsz = Size3_<int>(Size3_<double>(pyrA[0].size()-subWinSizeLv)/Size3_<double>(subWinStepLv)) + Size3(1,1,1);
    Tensor<double,1> rstMat(tempsz,0);//gj20130120 set inital value to 1 for productive pooling
    //set 0 for addictive pooling
    //rstMat.Print();
//    cout<<subWinSizeLv<<endl;
//    cout<<subWinStepLv<<endl;
    for (index = 0; index < (int)pyrA.size(); index++)
      {
        int lvl=0;
        if (downsample)
          {
            if (index ==(int)pyrA.size()-1) //HP part
              {
                lvl = 0;
              }
            else
              {
                lvl = index/nDir;
              }
            subWinSizeLv.height = max(subWinSize.height>>lvl,1);
            subWinSizeLv.width = max(subWinSize.width>>lvl,1);
            subWinStepLv.height =max( subWinStep.height>>lvl,1);
            subWinStepLv.width = max(subWinStep.width>>lvl,1);

            //if (boundary_cut == FILTER_BOUND_VALID)
            //{
            //  subWinSizeLv = subWinSizeLv + subWinSizeLv/2; //2*overlap + size
            //  subWinStepLv = subWinStepLv + subWinStepLv/2;
            //}
            sz.height = max(initSize.height>>lvl,1);
            sz.width = max(initSize.width>>lvl,1);
            //mylib::DisplayMat(gaussKernel);
            //weight = subWinSizeLv.volumn();
          }
        //cv::Mat flatKel=cv::Mat(subWinSizeLv,pyrA[index][0].type()-((pyrA[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1.0/weight));
        cv::Mat gaussKernel = mylib::GenGaussKer(subWinSizeLv.height,double(subWinSizeLv.height)/6.0,CV_64F);

        //pyrA[index].Print("pyrA");
        // Tensor<double,1>(gaussKernel).Print("gauss");
        //cout<<"subWinStepLv: ";
        //subWinStepLv.Print();

        //pyrA[index].Print("pa");
        //pyrB[index].Print("pb");
        //! 20130910 before localMean
        //char tempc;
        //cin>>tempc;
        mu_A[index] = pyrA[index].LocalMean(gaussKernel,subWinStepLv);
        if (debug)
                mu_A[index].Print("mua");
        mu_B[index] = pyrB[index].LocalMean(gaussKernel,subWinStepLv);
        if (debug)
                mu_B[index].Print("mub");

        //mu_A[index] = pyrA[index].LocalMean(flatKel,subWinStepLv);    //Dec 30 2012, use flat (usuall) average
        //mu_B[index] = pyrB[index].LocalMean(flatKel,subWinStepLv);    //Dec 30 2012, use flat (usuall) average
        /* if (sz.height==32)
    {
      Tensor<double,2>(mu_A[index]).Print();
      Tensor<double,2>(mu_B[index]).Print();
    }*/
        if (stsim2_modifer == MetricModifier::STSIM2_BASELINE||stsim2_modifer==MetricModifier::STSIM2_TUNE)
        {
          //L[index] = (mu_A[index].Abs() * mu_B[index].Abs() *2 + C0).Real() / (mu_A[index]*mu_A[index].Conjugate() + mu_B[index]*mu_B[index].Conjugate() + C0).Real();
          L[index] = ((mu_A[index] * mu_B[index]).Abs() *2 + C0).Real() / (mu_A[index].Square() + mu_B[index].Square() + C0).Real();
        }
        else if (stsim2_modifer == MetricModifier::STSIM2_NEW_L1)
          {

            double clipThrd = (mu_A[index].Abs().Sum()/double(mu_A[index].size().volumn()))[0]*0.1;
            if (clipThrd<4)
              clipThrd=4;
            clipThrd = clipThrd*clipThrd;
            L[index] = ((mu_A[index]-mu_B[index])*((mu_A[index]-mu_B[index]).Conjugate())).Real().MaxClip(clipThrd);
            L[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>(clipThrd)) - L[index];
            L[index] = L[index]/clipThrd;
            if(debug&&index==4)
              {
                mu_A[index].Print("mu_A");
                mu_B[index].Print("mu_B");
                (mu_A[index]-mu_B[index]).Abs().Real().Print("muA-muB");
                L[index].Print("L");
              }
          }
        else if (stsim2_modifer ==MetricModifier:: STSIM2_NEW_L2)
          {
            //L[index] = (One - (mu_A[index]-mu_B[index]).Abs()/510).Real();
            //L[index] = Tensor<double,1>(mu_A[index].size(),1)- (mu_A[index]-mu_B[index]).Abs().Real()/510;
            //L[index].Print();
            if (index == 12)//DC band only
              {
                double clipThrd = (mu_A[index].Abs().Sum()/double(mu_A[index].size().volumn()))[0]*0.1;
                if (clipThrd<4)
                  clipThrd=4;
                //mu_A[index].Print("mu A");
                //mu_B[index].Print("mu B");
                clipThrd = clipThrd*clipThrd;
                L[index] = ((mu_A[index]-mu_B[index])*((mu_A[index]-mu_B[index]).Conjugate())).Real().MaxClip(clipThrd);
                //L[index] = mu_A[index].AbsDiff(mu_B[index]).Real().MaxClip(clipThrd);
                //L[index].Print();
                L[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>(clipThrd)) - L[index];
                //L[index].Print();
                L[index] = L[index]/clipThrd;
                //if (debug)
                //	L[index].Print();
              }
            else
              L[index] = (mu_A[index].Abs() * mu_B[index].Abs() *2 + C0).Real() / (mu_A[index]*mu_A[index].Conjugate() + mu_B[index]*mu_B[index].Conjugate() + C0).Real();
          }
        else if (stsim2_modifer ==MetricModifier:: STSIM2_NEW_L3)
          {
            if (index >=8 && index <= 12)//DC band only and lower band
              {
                double clipThrd =(mu_A[index].Abs().Sum()/double(mu_A[index].size().volumn()))[0]*0.1;
                if (clipThrd<4)
                  clipThrd=4;
                //mu_A[index].Print("mu A");
                //mu_B[index].Print("mu B");
                clipThrd = clipThrd*clipThrd;
                L[index] = ((mu_A[index]-mu_B[index])*((mu_A[index]-mu_B[index]).Conjugate())).Real().MaxClip(clipThrd);
                //L[index] = mu_A[index].AbsDiff(mu_B[index]).Real().MaxClip(clipThrd);
                //L[index].Print();
                L[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>(clipThrd)) - L[index];
                //L[index].Print();
                L[index] = L[index]/clipThrd;
                //if (debug)
                //	L[index].Print();
              }
            else
              L[index] = (mu_A[index].Abs() * mu_B[index].Abs() *2 + C0).Real() / (mu_A[index]*mu_A[index].Conjugate() + mu_B[index]*mu_B[index].Conjugate() + C0).Real();

          }
        else if (stsim2_modifer == MetricModifier::STSIM2_NEW_L4)
          {
            //because the brightness perception is based on the ration, not the differences, so try to use another metric
            double clipMin = mu_A[index].Abs()(0,0,0)[0]*0.3;
            double clipMax = 4;
            //clipThrd>4?clipThrd=clipThrd:clipThrd=4;
            //clipThrd = clipThrd*clipThrd;
            L[index] = ((mu_A[index]-mu_B[index])/mu_B[index]).Abs().Real().MaxClip(clipMax).MinClip(clipMin);

            //L[index].Print();
            //L[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>(clipMax)) - L[index];
            //L[index].Print();
            //L[index] = L[index]/(clipMax-clipMin);
            double tempclip = log(clipMax-clipMin);
            for (int ii=0; ii<L[index].size().height;ii++)
              for (int jj=0; jj<L[index].size().width;jj++)
                L[index](ii,jj,0)[0]=log(L[index](ii,jj,0)[0])/tempclip;
            //L[index].Print();
          }

        //mu_A[index].Abs().Print();
        //mu_B[index].Abs().Print();
        //(mu_A[index]-mu_B[index]).Print();
        //(mu_A[index]-mu_B[index]).Conjugate().Print();
        //((mu_A[index]-mu_B[index])*((mu_A[index]-mu_B[index]).Conjugate())).Print();
        //pyrA[index].Abs().Real().Print("band A");
        //mu_A[index].Abs().Real().Print("mu_A abs");
        //(mu_A[index].Abs() * mu_B[index].Abs() *2).Print("num");
        //(mu_A[index]*mu_A[index].Conjugate() + mu_B[index]*mu_B[index].Conjugate() + C0).Real().Print("dem");
        //L[index].Print("L");
        //.Max(Vec<double,2>(clipThrd,0));
        //L[index] = Tensor<double,clipThrd - L[index];
        //max((mu_A[index]-mu_b[inxedx])*(mu_A[index]-mu_b[inxedx])/

        //L[index].Print();
        //here computed the abs(pyrA.* pyrA) that is (actuall) not equal to the true variance
        //for the complex number, it will be a covariance mx instead of var
        //which one is better should be studied.
        //use gaussian kernel (traditional)
        if (debug&&index==2)
          {
            //pyrA[index].Print("A");
            //pyrB[index].Print("B");
            //mu_A[index].Print("muA");
            //mu_B[index].Print("muB");
          }
        sigma2_A[index] = pyrA[index].LocalVariance(mu_A[index],gaussKernel,subWinStepLv);
        sigma2_B[index] = pyrB[index].LocalVariance(mu_B[index],gaussKernel,subWinStepLv);
//        if (debug&&index==2)
//          {
//            sigma2_A[index].Print();
//            sigma2_B[index].Print();
////            auto temp2 = (sigma2_A[index].Real()*sigma2_B[index].Real()).Sqrt() + C1;
////            for (int ii=0; ii< rstMat.size().height; ii++)
////              for (int jj=0; jj<rstMat.size().width; jj++)
////                {
////                  //double m = std::numeric_limits<double>::min();
////                  //double temp = (sigma2_A[index](ii,jj,0)[0]*sigma2_B[index](ii,jj,0)[0]+C1);
////                  cout<<left<<index<<"\t"<<temp2(ii,jj,0)[0]<<endl;
////                }
//          }
        //use flat kernel Dec 30 2012
        //sigma2_A[index] = pyrA[index].LocalVariance(mu_A[index],flatKel,subWinStepLv);
        //sigma2_B[index] = pyrB[index].LocalVariance(mu_B[index],flatKel,subWinStepLv);
        //sigma2_B[index].Print();
        //Try ignore the C term in order to work on smooth/smooth or smooth/texture region, gjin aug 24, 2011
        //C[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>::all(1));
        C[index] = ((sigma2_A[index].Real()*sigma2_B[index].Real()).Sqrt()*2 + C1)/(sigma2_A[index]+sigma2_B[index] +C1).Real();


//        if(debug&&index==2)
//          {
//            ((sigma2_A[index].Real()*sigma2_B[index].Real()).Sqrt()*2 + C1).Real().Print("num");
//            (sigma2_A[index]+sigma2_B[index] +C1).Real().Print("denum");
//            (((sigma2_A[index].Real()*sigma2_B[index].Real()).Sqrt()*2 + C1)/(sigma2_A[index]+sigma2_B[index] +C1).Real()).Print();
//            //C[index].Print();
//          }
        C01[index] = ComputeCrossTerm(pyrA[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth)),
            pyrA[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth)),
            pyrB[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth)),
            pyrB[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth)),
            subWinSizeLv-Size3(0,1,0), subWinStepLv);
        //C01[index].Print();
        //cout<<"pyrA[idx].size ";pyrA[index].size().Print();
        //pyrB[index].size().Print();
        //cout<<"sz ";sz.Print();
        C10[index] = ComputeCrossTerm(pyrA[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth)),
            pyrA[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth)),
            pyrB[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth)),
            pyrB[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth)),
            subWinSizeLv-Size3(1,0,0), subWinStepLv);
        //C10[index].Print();
//        cout<<"size "<<tsA.size()<<endl;
//        cout<<"sub "<<subWinSizeLv<<endl;
//        cout<<"step "<<subWinStepLv<<endl;
//        cout<<"L "<<L[index].size()<<endl;
//        cout<<"C "<<C[index].size()<<endl;
//        cout<<"C01 "<<C01[index].size()<<endl;
//        cout<<"C10 "<<C10[index].size()<<endl;
 
        Tensor<double,1> tempRstMat;
        if (stsim2_modifer == MetricModifier::STSIM2_TUNE)
          {
            tempRstMat = L[index]*TUNE_WEIGHT[index*4]+C[index]*TUNE_WEIGHT[index*4+1]+C01[index]*TUNE_WEIGHT[index*4+2]+C10[index]*TUNE_WEIGHT[index*4+3];
          }
        else
          tempRstMat=(L[index]*C[index]*C01[index]*C10[index]).Pow(0.25);//change weight for C gj20130120
        //tempRstMat.Print();
        //if (stsim2_modifer == STSIM2_BASELINE)
        //  tempRstMat = tempRstMat.Pow(0.25);// ((L[index]*C[index]*C01[index]*C10[index]).Pow(0.25));
        //else
        //{//    rstMat += ((L[index]*C[index]*C01[index]*C10[index]).Pow(0.25))*()
        //   if (index == 12)
        //	   tempRstMat = tempRstMat.Pow(0.25)*(14.0/3.0);//+= ((L[index]*C[index]*C01[index]*C10[index]).Pow(0.25))*(14.0/3);//(3.5)
        //   else
        //     tempRstMat = tempRstMat.Pow(0.25)*(28.0/39.0);// += ((L[index]*C[index]*C01[index]*C10[index]).Pow(0.25))*(28.0/39.0);//;*(42.0/52.0);
        //}

        rstMat=rstMat+tempRstMat; //gj20130120 to emphasize every term, use product instead of sum
        /*if (index==12) //lp
        rstMat += tempRstMat*(7);
      else if (index ==13 ) //hp
        rstMat += tempRstMat*(0.875);
      else if (index >=8) //level 1
        rstMat += tempRstMat*(0.875);
      else if (index >=4) //level 2
        rstMat += tempRstMat*(0.4375);
      else //level 3
        rstMat += tempRstMat*(0.2188);
*/
        //}
        //else if (stsim2_modifer == STSIM2_NEW_L2)
        // {
        //
        // }
        //else if (stsim2_modifer == STSIM2_NEW_L3)
        //{
        //
        // }
        //if (debug)
        //	rstMat.Print();
        //print paramters
        if (debug)
          {
            for (int ii=0; ii< rstMat.size().height; ii++)
              for (int jj=0; jj<rstMat.size().width; jj++)
                debugfile<<left<<index<<"\t"<<L[index](ii,jj,0)[0]<<"\t"<<C[index](ii,jj,0)[0]<<"\t"<<C01[index](ii,jj,0)[0]<<"\t"<<C10[index](ii,jj,0)[0]<<"\t"<<tempRstMat(ii,jj,0)[0]<<endl;
          }
      }

    Tensor<double,1> rstMatBackup;
    if (stsim2_modifer==MetricModifier::STSIM2_TUNE)
      rstMatBackup = rstMat;
    else
      {
        //rstMatBackup = (rstMat/double(pyrA.size())).Clone();///;.Pow(1.0/double(index)) gj20130120 use product, so do geometric normalization
        rstMatBackup = rstMat.Clone();
        rstMat = Tensor<double,1>(Size3(Size3_<double>(pyrA[0].size()-subWinSizeLv)/Size3_<double>(subWinStepLv)) + Size3(1,1,1),0);
      }
    if (debug)
      rstMatBackup.Print();
    if (debug)
      {
        debugfile<<"STSIM-1 score"<<"\t\t\t\t\t"<<rstMatBackup.Min()[0]/index/*diabled when using product pooling/index*/<<endl;
        debugfile<<"====== cross band ======="<<endl;
        debugfile<<"Band A"<<"\t"<<"Band B"<<"\t"<<"C00"<<endl;
      }

    if (stsim2_modifer!=MetricModifier::STSIM2_TUNE)
      {
        index = 0;
        for (int dr = 0; dr < nDir; dr++)
          for (int lv = 0; lv < nLevel-1; lv++)
            {
              C00[index] = ComputeCrossTerm( pyrA[ lv*nDir + dr].Abs(), pyrA[ (lv+1)*nDir + dr].Abs(),
                  pyrB[ lv*nDir + dr].Abs(), pyrB[ (lv+1)*nDir + dr].Abs(),
                  subWinSizeLv,subWinStepLv);
              //rst+= C00[index].Mean();
              rstMat+= C00[index];
              //  C00[index].Print();
              if(debug)
                {
                  for (int ii=0; ii< rstMat.size().height; ii++)
                    for (int jj=0; jj<rstMat.size().width; jj++)
                      debugfile<<lv*nDir+dr<<"\t"<<(lv+1)*nDir+dr<<"\t"<<C00[index](ii,jj,00)[0]<<endl;
                }
              index++;
            }

        for (int lv=0; lv< nLevel; lv++)
          for (int dr = 0; dr < nDir-1; dr++)
            for ( int p = dr+1; p < nDir; p++)
              {
                C00[index] = ComputeCrossTerm(pyrA[lv*nDir + dr].Abs(),pyrA[lv*nDir + p].Abs(),
                    pyrB[lv*nDir + dr].Abs(),pyrB[lv*nDir + p].Abs(),
                    subWinSizeLv, subWinStepLv);
                //rst+=C00[index].Mean();
                rstMat+= C00[index];
                //C00[index].Print();
                if (debug)
                  {
                    for (int ii=0; ii< rstMat.size().height; ii++)
                      for (int jj=0; jj<rstMat.size().width; jj++)
                        debugfile<<lv*nDir+dr<<"\t"<<(lv)*nDir+p<<"\t"<<C00[index](ii,jj,00)[0]<<endl;
                  }
                index++;
              }
        //rstMat = rstMat/double(crossbandNum);
        if (debug)
          rstMat.Print();
        double totalbands = double(pyrA.size()+crossbandNum);

        //gj20120120 use product in STSIM-1 but still use summation in C00, so the pooling is different//
        //rstMat = rstMatBackup*double(pyrA.size())/totalbands + rstMat*double(crossbandNum)/totalbands;

        rstMat = (rstMatBackup + rstMat)/totalbands;

        // turn off cross band, so set rstMat
        //rstMat = rstMatBackup; // Dec28 2012, turn off cross band

        //  rstMat = rstMat/double(pyrA.size());
        //pooling STSIM terms
        //if (rstMat.size().area() == 4)
        // rstMat.Print();
      }
    if (stsim2_pool_type == FeaturePoolType::FEATURE_POOL_AVE)
      {
        rst = rstMat.Mean()[0];
      }
    else
      {
        rst = rstMat.Min()[0];
      }
    if(debug)
      {
        debugfile<<"STSIM-2 score"<<"\t\t\t\t\t"<<rst<<endl;
        //return rst/(crossbandNum+pyrA.size())/sz.depth;
        debugfile.close();
        rstMat.Print();
      }
    return rst;
  }


  Mat EstimateVarForMahalanobis(const Mat& ts, Size3 wsize, Size3 stepsize)
  {
    Mat covar;
    Mat mean;
    //compute 2D only
    CV_Assert(wsize.depth==1);
    CV_Assert(ts.channels()==1);
    int numberofparam = 4*14;
    int numberofvector = (ts.size().height)/stepsize.height*(ts.size().width)/stepsize.width;
    //CV_MAKE_TYPE(CV_64F,cn);
    //CV_MAKETYPE(
    cv::Mat samples(numberofparam,numberofvector,CV_64F);
    Cube roi(Point3i(0,0,0),wsize);
    Tensor<double,1> block;
    int nLevel = 3;
    int nDir = 4;
    int nband = nLevel*nDir +2;
    vector<Tensor<double,2>> mu_A(nband);
    vector<Tensor<double,2>> sigma2_A(nband);
    vector<Tensor<double,2>> rho01(nband), rho10(nband);
    int crossbandNum = nDir*(nLevel-1) + nLevel* (mylib::combination(nDir,2));
    vector<Tensor<double,1>> C00 (crossbandNum);
    //int index = 0;
    Size3 subWinSize = wsize;
    Size3 subWinStep = stepsize;
    cv::Mat flatKel;//=cv::Mat(sz,CV_64F2[0].type()-((pyrA[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1/weight));
    cv::Mat gaussKernel = mylib::GenGaussKer(wsize.height,double(wsize.height)/6.0,CV_64F);
    cv::Mat temp;
    Tensor<double,1> A(ts);
    double C = 0.001;
    int col =0;
    for (int i=0; i<ts.size().height; i+=wsize.height)
      {
        for (int j=0; j<ts.size().width; j+=wsize.width)
          {
            block = A.GetBlockRef(roi);
            //Tensor<double,1>(block).Display();
            Steerable spA;
            spA.buildSCFpyr(Tensor<double,1>(block).ToComplex(),3,4,1,false);
            vector<Tensor<double,2>>& pyrA = spA.getSpaceDomainPyr();
            for (unsigned int index=0; index<pyrA.size(); index++)
              {
                // pyrA[index].Print();
                mu_A[index] = pyrA[index].LocalMean(gaussKernel,subWinStep);
                sigma2_A[index] = pyrA[index].LocalVariance(mu_A[index],gaussKernel,subWinStep);
                Size3 sz = pyrA[index].size();
                Tensor<double,2> im11 = pyrA[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth));
                Tensor<double,2> im12 = pyrA[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth));
                flatKel = Mat(Size(wsize.width-1,wsize.height),CV_64F,Scalar(1/double((wsize.width-1)*wsize.height)));
                Tensor<double,2> mu11 = im11.LocalMean(subWinSize,subWinStep);         //mu11.Print();
                Tensor<double,2> mu12 = im12.LocalMean(subWinSize,subWinStep);
                Tensor<double,2> sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
                Tensor<double,2> sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep);
                Tensor<double,2> sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
                rho01[index] = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
                im11 = pyrA[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth));
                im12 = pyrA[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth));
                flatKel = Mat(Size(subWinSize.width,subWinSize.height-1),CV_64F,Scalar(1/double(subWinSize.width*(subWinSize.height-1))));
                sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
                sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep);
                sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
                rho10[index] = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
                //cout<< mu_A[index].Abs().Real()(0,0,0)[0]<<endl<<sigma2_A[index].Abs().Real()(0,0,0)[0]<<endl<< rho01[index].Abs().Real()(0,0,0)[0] <<endl<< rho10[index].Abs().Real()(0,0,0)[0]<<endl;
                samples.at<double>(index*4,col) = mu_A[index].Abs().Real()(0,0,0)[0];
                samples.at<double>(index*4+1, col) = sqrt(sigma2_A[index].Abs().Real()(0,0,0)[0]);
                samples.at<double>(index*4+2,col) = rho01[index].Abs().Real()(0,0,0)[0];
                samples.at<double>(index*4+3,col) = rho10[index].Abs().Real()(0,0,0)[0];
              }
            cout<<roi.x<<"," <<roi.y<<endl;
            roi.y +=wsize.width;
            col++;
          }
        roi.x += wsize.height;
        roi.y =0;
        //col++;
      }
    //mylib::DisplayMat(samples, "samples",true);
    cv::calcCovarMatrix(samples,covar,mean,CV_COVAR_NORMAL|CV_COVAR_COLS);
    // mylib::DisplayMat(covar,"covar",true);
    //mylib::DisplayMat(mean,"mean",true);
    cv::Mat icovar;
    cv::invert(covar,icovar,DECOMP_SVD);
    //Tensor<double,1>(icovar).Print("icovar",true);
    return icovar;
  }


  double ComputeMahalanobis(const Mat& tsA,const Mat& tsB, Size3 subWinSize, Size3 subWinStep, const Mat& iMcovar, FilterBoundary boundary_cut)
  {
    const double C = 0.001;
    Steerable spA, spB;
    int nDir = 4; int nLevel = 3;
    int crossbandNum = nDir*(nLevel-1) + nLevel* (mylib::combination(nDir,2));
    Tensor<double,1> A(tsA),B(tsB);
    spA.buildSCFpyr(A.ToComplex(),nLevel,nDir,1,false);
    spB.buildSCFpyr(B.ToComplex(),nLevel,nDir,1,false);
    vector<Tensor<double,2>>& pyrA = spA.getSpaceDomainPyr();
    vector<Tensor<double,2>>& pyrB = spB.getSpaceDomainPyr();
    //int boundary_cut = FILTER_BOUND_HALF;
    if ( boundary_cut == FilterBoundary::FILTER_BOUND_HALF)//modify Dec 27 2011, not cut half but cut half + boundary
      {
        for (unsigned int i = 0; i< pyrA.size(); i++)
          {
            //pyrA[i].Print();
            //dec 27 2011, extend the block by size of block size + bounday (below and right)
            int hh = subWinSize.height*(int)ceil(float(pyrA[i].size().height)/2.0/float(subWinSize.height));
            int ww = subWinSize.width*(int)ceil(float(pyrA[i].size().width)/2.0/float(subWinSize.width));
            int dd = subWinSize.depth*(int)ceil(float(pyrA[i].size().depth)/2.0/float(subWinSize.depth));
            Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, hh,ww,dd);
            //Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, pyrA[i].size().height/2, pyrA[i].size().width/2, pyrA[i].size().depth);
            pyrA[i] = pyrA[i](roi);//no size changed
            //pyrA[i].Print();
            pyrB[i] = pyrB[i](roi);
            //pyrB[i].Print();
          }
      }
    Size3 sz = pyrA[0].size();
    //Tensor<double,1> rstMat(pyrA[0].size()/subWinStep - subWinSize+Size3(1,1,1));
    Tensor<double,1> rstMat(Size3(Size3_<double>(pyrA[0].size()-subWinSize)/Size3_<double>(subWinStep)) + Size3(1,1,1));
    //rstMat.Print();
    vector<Tensor<double,2>> mu_A(pyrA.size());
    vector<Tensor<double,2>> mu_B(pyrB.size());
    Tensor<double,2> mu11,mu12,mu21,mu22;
    vector<Tensor<double,2>> sigma2_A(pyrA.size());
    vector<Tensor<double,2>> sigma2_B(pyrA.size());
    vector<Tensor<double,2>> rho01_A(pyrA.size());
    vector<Tensor<double,2>> rho01_B(pyrA.size());
    vector<Tensor<double,2>> rho10_A(pyrA.size());
    vector<Tensor<double,2>> rho10_B(pyrA.size());
    Mat featureA(pyrA.size()*4,1,CV_64F);
    Mat featureB(pyrA.size()*4,1,CV_64F);
    Tensor<double,2> sigma1_cross,sigma2_cross;
    Tensor<double,2> sigma11,sigma12,sigma21,sigma22;

    vector<Tensor<double,1>> C00 (crossbandNum);
    int index = 0;

    double weight = sz.height*sz.width;
    cv::Mat flatKel=cv::Mat(sz,pyrA[index][0].type()-((pyrA[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1/weight));
    cv::Mat gaussKernel = mylib::GenGaussKer(subWinSize.height,double(subWinSize.height)/6.0,CV_64F);
    cv::Mat temp;
    for (index = 0; index < (int)pyrA.size(); index++)
      {
        mu_A[index] = pyrA[index].LocalMean(gaussKernel,subWinStep);
        mu_B[index] = pyrB[index].LocalMean(gaussKernel,subWinStep);
        sigma2_A[index] = pyrA[index].LocalVariance(mu_A[index],gaussKernel,subWinStep);
        sigma2_B[index] = pyrB[index].LocalVariance(mu_B[index],gaussKernel,subWinStep);
        Tensor<double,2> im11 = pyrA[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth));
        Tensor<double,2> im12 = pyrA[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth));
        Tensor<double,2> im21 = pyrB[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth));
        Tensor<double,2> im22 = pyrB[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth)),
            flatKel = Mat(Size(subWinSize.width-1,subWinSize.height),CV_64F,Scalar(1/double((subWinSize.width-1)*subWinSize.height)));
        mu11 = im11.LocalMean(subWinSize,subWinStep);         //mu11.Print();
        mu12 = im12.LocalMean(subWinSize,subWinStep);
        mu21 = im21.LocalMean(subWinSize,subWinStep);
        mu22 = im22.LocalMean(subWinSize,subWinStep);
        sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
        sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep);
        sigma21 = im21.LocalVariance(mu21,subWinSize,subWinStep);
        sigma22 = im22.LocalVariance(mu22,subWinSize,subWinStep);
        sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
        sigma2_cross = ((im21 * im22.Conjugate()).LocalMean(subWinSize,subWinStep) - mu21 * mu22.Conjugate());// sigma2_cross.Print();
        rho01_A[index] = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
        rho01_B[index] = (sigma2_cross + Vec2d(C,0))/((sigma21*sigma22).Sqrt()+Vec2d(C,0));// rho2.Print();
        im11 = pyrA[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth));
        im12 = pyrA[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth));
        im21 = pyrB[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth));
        im22 = pyrB[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth));
        flatKel = Mat(Size(subWinSize.width,subWinSize.height-1),CV_64F,Scalar(1/double(subWinSize.width*(subWinSize.height-1))));
        mu11 = im11.LocalMean(subWinSize,subWinStep);         //mu11.Print();
        mu12 = im12.LocalMean(subWinSize,subWinStep);
        mu21 = im21.LocalMean(subWinSize,subWinStep);
        mu22 = im22.LocalMean(subWinSize,subWinStep);
        sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
        sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep);
        sigma21 = im21.LocalVariance(mu21,subWinSize,subWinStep);
        sigma22 = im22.LocalVariance(mu22,subWinSize,subWinStep);
        //    sigma11.Print("sigma11");
        //   sigma12.Print("sigma12");
        //   sigma1_cross.Print("sigma1_cross");
        sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
        sigma2_cross = ((im21 * im22.Conjugate()).LocalMean(subWinSize,subWinStep) - mu21 * mu22.Conjugate());// sigma2_cross.Print();
        rho10_A[index] = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
        rho10_B[index] = (sigma2_cross + Vec2d(C,0))/((sigma21*sigma22).Sqrt()+Vec2d(C,0));// rho2.Print();

        // rho10_A[index].Print();
        // rho10_A[index].Abs().Real().Print();
        // cout<<rho10_A[index].Abs().Real()(0,0,0)[0]<<endl;;
        featureA.at<double>(index*4,0) = mu_A[index].Abs().Real()(0,0,0)[0];
        featureA.at<double>(index*4+1,0) = sqrt(sigma2_A[index].Abs().Real()(0,0,0)[0]);
        featureA.at<double>(index*4+2,0) = rho01_A[index].Abs().Real()(0,0,0)[0];
        featureA.at<double>(index*4+3,0) = rho10_A[index].Abs().Real()(0,0,0)[0];
        featureB.at<double>(index*4,0) = mu_B[index].Abs().Real()(0,0,0)[0];
        featureB.at<double>(index*4+1,0) = sqrt(sigma2_B[index].Abs().Real()(0,0,0)[0]);
        featureB.at<double>(index*4+2,0) = rho01_B[index].Abs().Real()(0,0,0)[0];
        featureB.at<double>(index*4+3,0) = rho10_B[index].Abs().Real()(0,0,0)[0];
      }
    //return rst/(crossbandNum+pyrA.size())/sz.depth;
    // mylib::DisplayMat(featureA,"fa");
    // mylib::DisplayMat(featureB,"fb");
    // mylib::DisplayMat((featureA-featureB).t());
    // mylib::DisplayMat(iMcovar,"imcovar");
    /*  Mat temp2(1,iMcovar.size().width,CV_64F,Scalar(0));
 // mylib::DisplayMat(temp2,"temp2");
  cv::gemm(featureA-featureB,iMcovar,1,featureA.t(),0,temp2,GEMM_1_T);
  Mat F = (featureA - featureB);
  mylib::DisplayMat(F,"F");
  rst = 0;
  for (int jj=0; jj<iMcovar.size().width; jj++)
  {
    temp2.at<double>(0,jj)=0;
    for (int ii=0; ii<iMcovar.size().height; ii++)
    {
      temp2.at<double>(0,jj) += F.at<double>(ii,0)*iMcovar.at<double>(ii,jj);
      cout<<"+"<<F.at<double>(ii,0)<<" x "<<iMcovar.at<double>(ii,jj)<<" = "<< temp2.at<double>(0,jj)<<endl;

    }
  }
  mylib::DisplayMat(temp2,"temp2");
  mylib::DisplayMat((featureA-featureB).t()*iMcovar);
  temp2 = (featureA - featureB).t()*iMcovar*(featureA-featureB);
  rst = temp2.at<double>(0,0);*/
    return cv::Mahalanobis(featureA,featureB,iMcovar);
  }

  //! this is the metric parser
  double Compare(const Mat& tsA, const Mat& tsB, CompareCriteria criteria, const Size3& subWinSize,const Size3& subWinStep, Printable_t param1, Printable_t param2, Printable_t param3, Printable_t param4, Printable_t param5, Printable_t param6, bool debug)
  {
    /*va_list ap;
  va_start(ap,criteria);
  vector<string> argname;
  vector<Printable_t> argdata;
  Printable_t printable;
  int count=0;
  for (count=0;count<paramNum;count++)
  {
    char* s = va_arg(ap,char*);
    printable.s = s;// va_arg(ap,char*);
    argname.push_back(printable.s);
      //cout<<Printable.s<<endl;;
    if (criteria==COMPARE_CRITERIA_SSIM) //(!strcmp(printable.s,"nDir")||!strcmp(printable.s,"nLevel")||!strcmp(printable.s,"pool_type")||!strcmp(printable.s,"boundary_cut"))
        printable.i = va_arg(ap,int);
    else if (criteria==COMPARE_CRITERIA_MAHALANOBIS)
        printable.m = va_arg(ap,cv::Mat*);
    else
        printable.d = va_arg(ap,double);
    argdata.push_back(printable);
      //cout<<Printable.i<<endl;
  }
  va_end(ap);
  */
    CV_Assert(tsA.size()==tsB.size());
    CV_Assert(tsA.channels()==tsB.channels());
    CV_Assert(tsA.channels()==1);
    Tensor<double,1> A(tsA);
    Tensor<double,1> B(tsB);
    if (criteria == CompareCriteria::COMPARE_CRITERIA_MSE)
      {
        int modifier = param1.i;//get the modifer
        int bsize = param2.i;
        int osize = param3.i;
        if (modifier == (int) MetricModifier::SE_MSE)
          {
            Tensor<double,1> A_low,A_right,B_low,B_right;
            A_low = A.GetBlock(Cube(osize+bsize,osize,0,osize,bsize,1));
            A_right = A.GetBlock(Cube(osize,osize+bsize,0,bsize,osize,1));
            B_low = B.GetBlock(Cube(osize+bsize,osize,0,osize,bsize,1));
            B_right=B.GetBlock(Cube(osize,osize+bsize,0,bsize,osize,1));
            return ComputeMSE(A_low,B_low)+ComputeMSE(A_right,B_right);
          }
        return ComputeMSE(tsA,tsB);
      }
    else if (criteria == CompareCriteria::COMPARE_CRITERIA_INTERP)
      {
        //CV_DbgAssert(paramNum==1);

        bool rst = ComputeAIM(tsA,tsB,param1.d);

        if (rst)
          return 1.0;
        else
          return 0;
      }
    else if (criteria == CompareCriteria::COMPARE_CRITERIA_SSIM)//|STSIM2_BASELINE|SSIM2_DCT|STSIM2_ADT_NATIVE|STSIM2_NEW_L1|STSIM2_NEW_L2|STSIM2_NEW_L3|STSIM2_NEW_L4|STSIM2_SE_MSE)
      {
        //CV_DbgAssert(paramNum==5);
        if (subWinSize == Size3(0,0,0))
          CV_Error(CV_StsBadSize,"subWinSize is un-set");
        int nDir, nLevel, bd_cut, pool_type, modifier;
        bool downsample;
        nDir = param2.i;
        nLevel = param1.i;
        bd_cut = param3.i;
        pool_type = param4.i;
        modifier = param5.i;
        downsample = param6.b;
        /*
    count=0;
    for (string name : argname)
    {
      if (!name.compare("nDir"))
      {
        nDir = argdata[count].i;
      }
      else if (!name.compare("nLevel"))
      {
        nLevel = argdata[count].i;
      }
      else if (!name.compare("boundary_cut"))
      {
        bd_cut=argdata[count].i;
      }
      else if (!name.compare("pool_type"))
      {
        pool_type = argdata[count].i;
      }
      else if (!name.compare("metric_modifier"))
      {
        modifier = argdata[count].i;
      }
      else
      {
        CV_DbgAssert(0);
      }
      count++;
    }
    */
        Tensor<double,2> T1 = A.ToComplex();
        Size3 tempSubSize,tempSubStep;
        if (T1.size().height<=subWinSize.height&&bd_cut==(int)FilterBoundary::FILTER_BOUND_HALF)
          {
            tempSubSize = T1.size()/2+Size3(0,0,1);
            tempSubStep = T1.size()/2+Size3(0,0,1);
          }
        else
          {
            tempSubSize = subWinSize;
            tempSubStep = tempSubSize;
            //!20130912 do not exam boundary  tempSubStep = Size3(subWinSize.height/4,subWinSize.width/4,1);/*subWinStep;*///20130521 use 1/4 of subwinSize as the step size in order to do blk + LU boundary
          }



        if (modifier==(int)MetricModifier::STSIM3_LSE) //gj01142013
        {
            return ComputeSTSIM3_LSE(tsA,tsB,tempSubSize,tempSubStep,nLevel,nDir,downsample,FilterBoundary(bd_cut),debug);
          }
        else
          return ComputeSTSIM2(tsA,tsB,tempSubSize,tempSubStep,nLevel,nDir,downsample,/*FILTER_BOUND_HALF*/FilterBoundary(bd_cut), FeaturePoolType(pool_type), MetricModifier(modifier),debug);
        //else
        //  return T1.ComputeSSIM(T2,subWinSize,subWinStep,(int)param1,(int)param2, FILTER_BOUND_VALID, param4, param5);
      }
    else if (criteria == CompareCriteria::COMPARE_CRITERIA_SUBJECTIVE)
      {
        //TBD!
        return 1.0;
      }
    else if (criteria== CompareCriteria::COMPARE_CRITERIA_SAD)
      {
        return ComputeSAD(tsA,tsB);
      }
    else if (criteria== CompareCriteria::COMPARE_CRITERIA_MAHALANOBIS)
      {
        //CV_DbgAssert(count==2);
        return ComputeMahalanobis(tsA,tsB,subWinSize,subWinStep,*param1.m);
      }
    else if (criteria== CompareCriteria::COMPARE_CRITERIA_LRI)
      {
        return ComputeLRI(tsA,tsB);
      }
    //! 20130820  temp un-available
    //else if (criteria == COMPARE_CRITERIA_SVM)
    //{
    //  return ComputeSVMMetric(ts,subWinSize,subWinStep);
    // }
    else
      return INT_MAX;
    return INT_MAX;
  }

  Tensor<double,1> ComputeCrossTerm(const Mat& im11, const Mat& im12, const Mat& im21, const Mat& im22, const Size3& subWinSize, const Size3& subWinStep)
  {

    if (im11.size()!= im12.size() ||
        im11.size()!= im21.size() ||
        im11.size()!= im22.size())
      {
        CV_Error(CV_StsUnmatchedSizes,"im11,im12,im21,im22 must in same size");
      }
    CV_Assert(im11.channels()==im12.channels());
    CV_Assert(im21.channels()==im22.channels());
    CV_Assert(im11.channels()==im22.channels());
    CV_Assert(im11.channels()==2);
    //im11.Print();
    //im12.Print();
    //im21.Print();
    //im22.Print();
    const double C = 0.001;
    Tensor<double,2> mu11,mu12,mu21,mu22;
    Tensor<double,2> sigma1_cross,sigma2_cross;
    Tensor<double,2> rho1,rho2;
    Tensor<double,2> sigma11,sigma12,sigma21,sigma22;
    Tensor<double,1> rst;
    Mat flatKel;
    Tensor<double,2> t11(im11), t12(im12), t21(im21), t22(im22);
    //if (im11.size().height> im11.size().width)
    //	flatKel = Mat(Size(subWinSize.width-1,subWinSize.height),CV_64F,Scalar(1/double((subWinSize.width-1)*subWinSize.height)));
    //else if (im11.size().height < im11.size().width)
    //	flatKel = Mat(Size(subWinSize.width,subWinSize.height-1),CV_64F,Scalar(1/double(subWinSize.width*(subWinSize.height-1))));
    //else
    flatKel = Mat(Size(subWinSize.width,subWinSize.height),CV_64F,Scalar(1/double(subWinSize.width*subWinSize.height)));
    //mu11 = im11.LocalMean(flatKel,subWinStep); mu11.Print();
    mu11 = t11.LocalMean(subWinSize,subWinStep);      //   mu11.Print();
    //mu12 = im12.LocalMean(flatKel,subWinStep); //mu12.Print();
    //mu21 = im21.LocalMean(flatKel,subWinStep); //mu21.Print();
    //mu22 = im22.LocalMean(flatKel,subWinStep); //mu22.Print();
    mu12 = t12.LocalMean(subWinSize,subWinStep);//  mu12.Print();
    mu21 = t21.LocalMean(subWinSize,subWinStep);//  mu21.Print();
    mu22 = t22.LocalMean(subWinSize,subWinStep);//  mu22.Print();
    sigma11 = t11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
    //sigma11 = im11.LocalVariance(mu11,flatKel,subWinStep); sigma11.Print();
    //sigma12 = im12.LocalVariance(mu12,flatKel,subWinStep); //sigma12.Print();
    //sigma21 = im21.LocalVariance(mu21,flatKel,subWinStep); //sigma21.Print();
    //sigma22 = im22.LocalVariance(mu22,flatKel,subWinStep); //sigma22.Print();
    sigma12 = t12.LocalVariance(mu12,subWinSize,subWinStep); //sigma12.Print();
    sigma21 = t21.LocalVariance(mu21,subWinSize,subWinStep); //sigma21.Print();
    sigma22 = t22.LocalVariance(mu22,subWinSize,subWinStep); //sigma22.Print();

    //sigma1_cross = ((im11 * im12.Conjugate()).Filter2D(flatKel).SubSample(subWinStep) - mu11 * mu12.Conjugate()); sigma1_cross.Print();
    sigma1_cross = ((t11 * t12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
    sigma2_cross = ((t21 * t22.Conjugate()).LocalMean(subWinSize,subWinStep) - mu21 * mu22.Conjugate()); //sigma2_cross.Print();
    //((sigma11*sigma12).Sqrt()+Vec2d(C,0)).Print();
    //((sigma21*sigma22).Sqrt()).Print();

    rho1 = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print("rho1");
    //(sigma11*sigma12).Sqrt().Print();
    rho2 = (sigma2_cross + Vec2d(C,0))/((sigma21*sigma22).Sqrt()+Vec2d(C,0)); //rho2.Print("rho2");
    rst = Tensor<double,1>(rho1.size(),Vec<double,1>::all(1)) - (rho1-rho2).Abs().Real()/2;
    //rst.Print();
    return rst;
  }




}//end of tensor

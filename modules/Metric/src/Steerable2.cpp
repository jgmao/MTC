#include <Steerable2.h>
#include <fftw++.h>
#include <Imlib.h>
#include <Array.h>
#include <output.h>
namespace metric{
Steerable2::Steerable2()
{

}
Steerable2::Steerable2(c_real_ref im)
{
  this->im = im;
  h = im.size().height;
  w = im.size().width;
  align = sizeof(FComplex);
  K = 1;
  maxscale = 3;
  int h1 = h;
  int w1 = w;

  fLP0 = new  Ad(h,w,align);
  fHP0 = new Ad(h,w,align);
  fHPtmp = new Ad(h,w,align);
  fLPtmp = new Ad(h,w,align);
  L0 = vector<Ac*>(maxscale);
  L1 = vector<Ac*>(maxscale+1);
  finput = vector<Ac*>(maxscale);
  fim = vector<Ac*>(maxscale);
  conv = vector<Ac*>(maxscale);
  ftmp = vector<Ac*>(maxscale);
  Hp = new Ac(h,w,align);
  Lp = new Ac(h,w,align);
  fB = vector<Ad*>(K);
  for (int k=0; k<K;k++)
    fB[k] = new Ad(h,w,align);
  for (int i=0; i<maxscale;i++)
  {

      w1 = w;//pow(2.0,i);
      h1 = h;//pow(2.0,i);
      finput[i]=new Ac(h1,w1,align);
      fim[i] = new Ac(h1,w1,align);
      conv[i] = new Ac(h1,w1,align);
      ftmp[i] = new Ac(h1,w1,align);
      L0[i] = new Ac(h1,w1,align);
      L1[i] = new Ac(h1,w1,align);

  }
  for (int k=0; k<maxscale*K+2;k++)
    pyr_space.push_back(Tensor<double,2>(h,w));
  //w1 = w/pow(2.0,maxscale);
  //h1 = h/pow(2.0,maxscale);

  L1[maxscale] = new Ac(h,w,align);
  B = vector<Ac*>(K*maxscale+2);
  A = vector<Ac*>(K);

  for (int i = 0; i <maxscale; i++)
  {
      h1 = h; // pow(2.0, i);
      w1 = w; // pow(2.0, i);

      for (int j = 0; j < K; j++)
      {

        B[K*i+j] = new Ac(h1,w1,align);
        if (i==0)
          A[j] = new Ac(h,w,align);
      }

  }
  B[K*maxscale] = new Ac(h,w,align);
  B[K*maxscale+1] = new Ac(h1,w1,align);

  for (int i=0; i< h; i++)
    for (int j=0; j<w;j++)
    {
      (*finput[0])(i,j).re=im.at<double>(i,j);
      (*finput[0])(i,j).im=0.0;

      (*L1[0])(i,j).re = im.at<double>(i,j);
      (*L1[0])(i,j).im = 0.0;
    }
  buildFilters(maxscale,K,h,w);
}

Steerable2::~Steerable2()
{
        //free(finput);
  for (int i=0; i<maxscale;i++)
  {

      delete finput[i];
      delete fim[i];
      delete conv[i];
      delete ftmp[i];
      delete L0[i];
      delete L1[i];
      for (int j=0; j<K;j++)
      {
        if (i==0)
          {
            delete fB[j];
          delete A[j];
          }
      }
  }
  for (int i=0; i<(int) B.size();i++)
    delete B[i];
  delete L1[maxscale];
}
int Steerable2::buildFilters(int maxscale, int K, int h, int w)
{
  CV_Assert(h==w);
  //note there is no downsampling
  double normfactor = sqrt(2.0 * K - 1) * exp(lgamma(K)) / sqrt(K * (exp(lgamma(2.0 * K))));
  fLP = vector<Ad*>(maxscale+1);
  fHP = vector<Ad*>(maxscale+1);
  int a = w;
  for (int i=0; i<maxscale+1;i++)
  {
      a /=2;
      fLP[i] = new Ad(h,w,align);
      fHP[i] = new Ad(h,w,align);
      genLPfilter(*fLP[i],w,h,a/2,a);
      genHPfilter(*fHP[i],w,h,a/2,a);
  }
  for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
      {
        int k, l, m, n, t;
        double theta;
        if (i <  w/2) {k =     i;}
        if (i >= w/2) {k = i - w;}
        if (j <  h/2) {l =     j;}
        if (j >= h/2) {l = j - h;}

        if (i <  w/2) {m = i + w/2;}
        if (i >= w/2) {m = i - w/2;}
        if (j <  h/2) {n = j + h/2;}
        if (j >= h/2) {n = j - h/2;}

        for (t = 0; t < K; t ++)
          {
            theta = std::atan2(l, k) - t * M_PI / K;
            (*fB[t])(j,i) = normfactor * pow(2.0 * cos(theta), K-1);
          }
      }
  return 1;
}

int Steerable2::deleteFilters(void)
{
  for (int i=0; i<(int)fLP.size();i++)
    delete fLP[i];
  for (int j=0; j<(int)fHP.size();j++)
    delete fHP[j];
  for (int k=0; k<K; k++)
    delete fB[k];
  return 1;
}

int Steerable2::decompose()
  {
          int i, j, k, t, w1, h1, w2, h2, scale;
          double a, normfactor;
          fftwpp::fftw::maxthreads = get_max_threads();
          fftwpp::fft2d Forward(h,w,-1,*finput[0],*fim[0]);
#if STEER_DBG
          cout<<*finput[0]<<endl;
#endif
          Forward.fft(*finput[0],*fim[0]); //forward

#if STEER_DBG
          cout<<*fim[0]<<endl;
          Tensor<double,1> temp(h,w);
          for (j=0; j<h;j++)
            for (i=0;i<w;i++)
              temp(j,i) = 255*(*fLP[0])(j,i);
          temp.SaveBlock("fLP0.tif");
          for (j=0; j<h;j++)
            for (i=0;i<w;i++)
              temp(j,i) = 255*(*fHP[0])(j,i);
          temp.SaveBlock("fHP0.tif");

          vector<Ad*> tmpar(K);
          for (int k=0; k<K+1;k++)
          {
            tmpar[k] = new Ad(h,w,align);
          }
          for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
              {
                int k, l, m,n,t;
                if (i <  w/2) {k =     i;}
                if (i >= w/2) {k = i - w;}
                if (j <  h/2) {l =     j;}
                if (j >= h/2) {l = j - h;}

                if (i <  w/2) {m = i + w/2;}
                if (i >= w/2) {m = i - w/2;}
                if (j <  h/2) {n = j + h/2;}
                if (j >= h/2) {n = j - h/2;}
                (*tmpar[K])(n,m) = (*fHP[0])(j,i);
              }
          saveBand(*tmpar[K],"fB0Hp.tif");
          for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
            {
               int k, l, m, n, t ;
               if (i <  w/2) {k =     i;}
               if (i >= w/2) {k = i - w;}
               if (j <  h/2) {l =     j;}
               if (j >= h/2) {l = j - h;}

               if (i <  w/2) {m = i + w/2;}
               if (i >= w/2) {m = i - w/2;}
               if (j <  h/2) {n = j + h/2;}
               if (j >= h/2) {n = j - h/2;}

               for (t = 0; t < K; t ++)
               {
                  (*tmpar[t])(n,m) = (*fB[t])(j,i) * (*fLP[0])(j,i) * (*fHP[1])(j,i);
               }
            }
          for (int t=0; t<K;t++)
            saveBand(*tmpar[t],"fB0"+std::to_string(t)+".tif");
#endif
         //get Hp
          fourier2spatialband1(w,h,*fHP[0],*B[K*maxscale],*conv[0],*fim[0],*ftmp[0]);
          for (t = 0; t < K; t ++)
          {
            fourier2spatialband3(K, w, h, *fLP[0], *fHP[1], *fB[t], *B[t], *conv[0], *fim[0], *ftmp[0]);
#if STEER_DBG
            cout<<"conv"<<t<<endl;
            cout<<*conv[0]<<endl;
            cout<<"B "<<t<<endl;
            cout<<*B[t]<<endl;
#endif
          }
          fourier2spatialband1(w, h, *fLP[1], *L1[0], *conv[0], *fim[0], *ftmp[0]);
#if STEER_DBG
          cout<<"conv_next"<<endl;
          cout<<*conv[0]<<endl;
#endif
//          //subsampling for the next scale
//          //==============================

          w2 = w; // 2;
          h2 = h; // 2;

//          //iterative part of the decomposition
//          //===================================
            w1 = w;
            h1 = h;
            for (scale = 1; scale < maxscale; scale ++)
            {
               //don't do downsample
               //down2Freq(*conv[scale-1],*fim[scale],h1,w1);
               // instead use the same
               for (j=0; j< h1;j++)
                 for (i=0;i<w1;i++)
                   (*fim[scale])(j,i) = (*conv[scale-1])(j,i);
               //w1 /=2 ;//current level size
               //h1 /=2 ;
               #if STEER_DBG
               for (j=0; j<h;j++)
                 for (i=0;i<w;i++)
                   temp(j,i) = 255*(*fLP[scale])(j,i);
               temp.SaveBlock("fLPtmp"+std::to_string(scale)+".tif");
               for (j=0; j<h;j++)
                 for (i=0;i<w;i++)
                   temp(j,i) = 255*(*fHP[scale])(j,i);
               temp.SaveBlock("fHPtmp"+std::to_string(scale)+".tif");


               for (j = 0; j < h1; j++)
                   for (i = 0; i < w1; i++)
                   {
                        int k, l, m,n,t;
                        //double theta;
                        if (i <  w1/2) {k =     i;}
                        if (i >= w1/2) {k = i - w1;}
                        if (j <  h1/2) {l =     j;}
                        if (j >= h1/2) {l = j - h1;}

                        if (i <  w/2) {m = i + w/2;}
                        if (i >= w/2) {m = i - w/2;}
                        if (j <  h/2) {n = j + h/2;}
                        if (j >= h/2) {n = j - h/2;}
                        for (t = 0; t < K; t ++)
                        {
                            //theta = atan2(l, k) - t * M_PI / K;
                            //(*fB[t])(j,i) = normfactor * pow(2.0 * cos(theta), K-1);
                            (*tmpar[t])(n,m) = (*fB[t])(j,i) * (*fLP[scale])(j,i) * (*fHP[scale+1])(j,i);
                        }
                    }
               for (t=0;t<K;t++)
                 saveBand(*tmpar[0],"fB"+std::to_string(scale)+std::to_string(t)+".tif");
#endif
               for (t = 0; t < K; t ++)
               {
                   fourier2spatialband2(w1, h1, *fHP[scale+1], *fB[t], *B[t+scale*K], *conv[scale], *fim[scale], *ftmp[scale]);
               }
               if (scale==maxscale-1)//get LP
               {
                 fourier2spatialband1(w1, h1, *fLP[maxscale], *B[maxscale*K+1], *conv[scale], *fim[scale], *ftmp[scale]);
#if STEER_DBG
                 for (j = 0; j < h1; j++)
                   for (i = 0; i < w1; i++)
                     {
                       int k, l, m,n,t;
                       //double theta;
                       if (i <  w1/2) {k =     i;}
                       if (i >= w1/2) {k = i - w1;}
                       if (j <  h1/2) {l =     j;}
                       if (j >= h1/2) {l = j - h1;}

                       if (i <  w/2) {m = i + w/2;}
                       if (i >= w/2) {m = i - w/2;}
                       if (j <  h/2) {n = j + h/2;}
                       if (j >= h/2) {n = j - h/2;}
                       (*tmpar[K])(n,m) = (*fLP[scale+1])(j,i);
                     }
                     saveBand(*tmpar[K],"fB3Lp.tif");
#endif
               }
               else
                 fourier2spatialband1(w1, h1, *fLP[scale+1], *L1[scale], *conv[scale], *fim[scale], *ftmp[scale]);

          }

         // convert to Mat
          int hb,wb;
          for (k=0;k<(int)B.size();k++)
          {
            hb = B[k]->Nx();
            wb = B[k]->Ny();
            //pyr_space.push_back(Tensor<double,2>(h,w));
            for (int j=0; j<hb;j++)
              for (int i=0; i<wb; i++)
              {
                  pyr_space[k](j,i)=Vec<double,2>((*B[k])(j,i).re,(*B[k])(j,i).im);
              }
            pyr_space[k]=pyr_space[k].Crop(Point3i(h/4,w/4,0),Size3(h/2,w/2,1));
#if STEER_DBG
            Tensor<double,1> absband = pyr_space[k].Abs().Real();
            double mintmp,maxtmp;
            cv::minMaxLoc(absband,&mintmp,&maxtmp);
            absband = absband - mintmp;
            absband =absband/(maxtmp-mintmp);
            absband = absband*255;

            absband.SaveBlock("band"+std::to_string(k)+".tif");
#endif
          }
          #if STEER_DBG
          cout<<"input"<<endl;
          cout<<im<<endl;
          pyr_space[13].Print("LP");
#endif
          return 1;
  }

int Steerable2::down2Freq(Array::array2<FComplex>& L0, Array::array2<FComplex>&L1, int h0, int w0)
{
  //L0 is prior level
  //w0 and h0 are width and height of L0
  //L1 is next level
  //this only useful when using fftw++ for represent half (left w/2+1) real to complex fft
  int h1 = h0/2;
  int w1 = w0/2;
  int i,j;
  for (j=0; j<h1/2+1;j++)
    for (i=0; i< w1/2+1; i++)
      L1(j,i) = L0(j,i);

  for (j=h1/2+1; j<h1;j++)
    for (i=0; i<w1/2+1; i++)
      L1(j,i) = L0(j+h1,i);

  for (j=0; j<h1/2+1;j++)
    for (i=w1/2+1; i<w1; i++)
      L1(j,i) = L0(j,i+w1);

  for (j=h1/2+1; j<h1;j++)
    for (i=w1/2+1;i<w1;i++)
      L1(j,i) = L0(j+h1,i+w1);
  return 1;
}

int Steerable2::genLPfilter(Array::array2<double>& LP0, int w, int h, double x1, double x2)
  {
     double r, x;
     int i, j, k, l;

     for (j = 0; j < h; j++)
       {
          for (i = 0; i < w; i++)
            {
              if (i <= w/2) {k = i;}
              if (i >  w/2) {k = i - w;}

	      if (j <= h/2) {l = j;}
	      if (j >  h/2) {l = j - h;}

	      r = sqrt(k*k + l*l);

	      if (r < x1) {LP0(j,i) = 1.0;}
	      if (r > x2) {LP0(j,i) = 0.0;}

	      if ((r >= x1) && (r <= x2))
		{
		  x = M_PI / 4.0 * (1.0 + (r - x1) / (x2 - x1));
		  LP0(j,i) = cos(M_PI / 2.0 * (log(x * 4.0 / M_PI)/ log(2.0)));
		}
	    }
      }
    return 1;
  }

  int Steerable2::genHPfilter(Array::array2<double>& HP, int w, int h, double x1, double x2)
  {
     double r, x;
     int i, j, k, l;

     for (j = 0; j < h; j++)
       {
          for (i = 0; i < w; i++)
            {

	      if (i <= w/2) {k = i;}
	      if (i >  w/2) {k = i - w;}

	      if (j <= h/2) {l = j;}
	      if (j >  h/2) {l = j - h;}

	      r = sqrt(k*k + l*l);

	      if (r < x1) {HP(j,i) = 0.0;}
	      if (r > x2) {HP(j,i) = 1.0;}
	      if ((r >= x1) && (r <= x2))
		{
		  x = M_PI / 4.0 * (1.0 + (r - x1) / (x2 - x1));
		  HP(j,i) = cos(M_PI / 2.0 * (log (x * 2.0 / M_PI) / log(2.0)));
		}
	  }
      }
    return 1;
  }



  int Steerable2::fourier2spatialband3(int K, int w, int h,
                                         Array::array2<double>& otf1,
                                         Array::array2<double>& otf2,
                                         Array::array2<double>& otf3,
                                         Array::array2<FComplex>& BP,
                                         Array::array2<FComplex>& conv,
                                         Array::array2<FComplex>& fim,
                                         Array::array2<FComplex>& ftmp)
     {
        int i, j;
        div_t tmp = div(K, 2);
//        //int wp = w/2+1;
//        cout<<"otf1"<<endl;
//        cout<<otf1<<endl;
//        cout<<"otf2"<<endl;
//        cout<<otf2<<endl;
//        cout<<"otf3"<<endl;
//        cout<<otf3<<endl;
//        cout<<"fim"<<endl;
//        cout<<fim<<endl;
        for (j = 0; j < h; j++)
        {
            for (i = 0; i < w; i++)
              {
                int x = h-j-1;
                int y = w-i-1;
                conv(j,i).re = fim(j,i).re * otf1(x,y) * otf2(x,y) * otf3(x,y);
                conv(j,i).im = fim(j,i).im * otf1(x,y) * otf2(x,y) * otf3(x,y);
              }
        }
       // cout<<"conv"<<endl;
       // cout<<conv<<endl;
        fftwpp::fft2d Backward(h,w,1,conv,BP);
        Backward.fftNormalized(conv,BP);
      //  cout<<"BP"<<endl;
     //   cout<<BP<<endl;
//        if (tmp.rem == 1)
//          {
//            for (j = 0; j < h; j++)
//              {
//                for (i = 0; i < w; i++)
//                {
//                  BP[i + w * j] = ftmp[i + w * j][0];
//                }
//              }
//          }
//        if (tmp.rem == 0)
//          {
//            for (j = 0; j < h; j++)
//            {
//              for (i = 0; i < w; i++)
//                {
//                  BP[i + w * j] = ftmp[i + w * j][1];
//                }
//            }
//          }
        return 1;
     }

  int Steerable2::fourier2spatialband2(int w, int h, Array::array2<double>& otf1,
                                         Array::array2<double>& otf2, Array::array2<FComplex>& BP,
                                         Array::array2<FComplex>& conv, Array::array2<FComplex>& fim,
                                         Array::array2<FComplex>& ftmp
                                         )
  {
        int i, j;
        //int wp = w/2+1;
//        cout<<"otf1"<<endl;
//        cout<<otf1<<endl;
//        cout<<"otf2"<<endl;
//        cout<<otf2<<endl;
//        cout<<"fim"<<endl;
//        cout<<fim<<endl;
        for (j = 0; j < h; j++)
        {
            for (i = 0; i < w; i++)
              {
                int x = h-j-1;
                int y = w-i-1;
                conv(j,i).re = fim(j,i).re * otf1(x,y) * otf2(x,y);
                conv(j,i).im = fim(j,i).im * otf1(x,y) * otf2(x,y);
              }
        }
        //cout<<"conv"<<endl;
       // cout<<conv<<endl;
        fftwpp::fft2d Backward(h,w,1,conv,BP);
        Backward.fftNormalized(conv,BP);
       // cout<<"BP"<<endl;
       // cout<<BP<<endl;
        return 1;
  }

  int Steerable2::fourier2spatialband1(int w, int h, Array::array2<double>& otf1,
                                         Array::array2<FComplex>& BP,
                                         Array::array2<FComplex>& conv, Array::array2<FComplex>& fim,
                                        Array::array2<FComplex>& ftmp
                                         )
  {
        int i, j;
//        int wp = w/2+1;
        for (j = 0; j < h; j++)
        {
            for (i = 0; i < w; i++)
              {
                int x = h-j-1;
                int y = w-i-1;
                conv(j,i).re = fim(j,i).re * otf1(x,y);
                conv(j,i).im = fim(j,i).im * otf1(x,y);
              }
        }
        fftwpp::fft2d Backward(h,w,1,conv,BP);
        Backward.fftNormalized(conv,BP);
        return 1;
  }

  void Steerable2::saveBand(Ad &band,string name)
  {
    Tensor<double,1> temp(band.Nx(),band.Ny());
    int h = band.Nx();
    int w = band.Ny();
    for (int j=0; j<h;j++)
      for (int i=0;i<w;i++)
        temp(j,i) = band(j,i);
    double mintmp,maxtmp;
    minMaxLoc(temp,&mintmp,&maxtmp);
    temp = temp - mintmp;
    temp = temp / (maxtmp-mintmp);
    temp = temp*255;
    temp.SaveBlock(name);
  }

}

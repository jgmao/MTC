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
  K = 4;
  maxscale = 3;
  int h1 = h;
  int w1 = w;
  L0 = vector<Ac*>(maxscale);
  L1 = vector<Ac*>(maxscale+1);
  finput = vector<Ac*>(maxscale);
  fim = vector<Ac*>(maxscale);
  conv = vector<Ac*>(maxscale);
  ftmp = vector<Ac*>(maxscale);

  for (int i=0; i<maxscale;i++)
  {

      w1 = w/pow(2.0,i);
      h1 = h/pow(2.0,i);
      finput[i]=new Ac(h1,w1,align);
      fim[i] = new Ac(h1,w1,align);
      conv[i] = new Ac(h1,w1,align);
      ftmp[i] = new Ac(h1,w1,align);
      L0[i] = new Ac(h1,w1,align);
      L1[i] = new Ac(h1,w1,align);

  }

  w1 = w/pow(2.0,maxscale);
  h1 = h/pow(2.0,maxscale);

  L1[maxscale] = new Ac(h,w,align);
  B = vector<Ac*>(maxscale);
  B0 = vector<Ac*>(maxscale);
  B1 = vector<Ac*>(maxscale);
  B2 = vector<Ac*>(maxscale);
  A = vector<Ac*>(K);

  for (int i = 0; i <maxscale; i++)
  {
      h1 = h / pow(2.0, i);
      w1 = w / pow(2.0, i);
    //  B[i] = vector<Ac*>(K);

      for (int j = 0; j < K; j++)
      {

        B[K*i+j] = new Ac(h1,w1,align);
        //B.push_back(Ac(h1,w1,align));
        if (i==0)
          A[j] = new Ac(h,w,align);
        if (i==0)
          B0[j] = new Ac(h1,w1,align);
        if (i==1)
          B1[j] = new Ac(h1,w1,align);
        if (i==2)
          B2[j] = new Ac(h1,w1,align);
      }

  }

  for (int i=0; i< h; i++)
    for (int j=0; j<w;j++)
    {
      (*finput[0])(i,j).re=im.at<double>(i,j);
      (*finput[0])(i,j).im=0.0;

      (*L1[0])(i,j).re = im.at<double>(i,j);
      (*L1[0])(i,j).im = 0.0;
    }
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
        delete B[i][j];

        if (i==0)
          delete A[j];
      }
  }
  for (int j=0; j<K;j++)
  {
    delete B0[j];
    delete B1[j];
    delete B2[j];
  }
  delete L1[maxscale];
}

int Steerable2::decompose()
  {
          int i, j, k, t, w1, h1, w2, h2, scale;
          double a, normfactor;
          normfactor = sqrt(2.0 * K - 1) * exp(lgamma(K)) / sqrt(K * (exp(lgamma(2.0 * K))));
          fftwpp::fftw::maxthreads = get_max_threads();

          fftwpp::fft2d Forward(h,w,-1,*finput[0],*fim[0]);
//          fftwpp::crfft2d ifft(h,w,fim,finput);
          cout<<*finput[0]<<endl;

          Forward.fft(*finput[0],*fim[0]); //forward
          cout<<*fim[0]<<endl;
          //ifft.fftNormalized(fim,finput);
          //cout<<finput<<endl;


          //fftw_complex *fim, *conv, *ftmp, *finput;


//          int i, j, k, t, w1, h1, w2, h2, scale;
//          Mat fLP0, fHP0, fLPtmp, fHPtmp;
//          Mat fB, tmpar, tmparL, tmparH, tmparall;
//          double a, normfactor;

//          //exp(lgamma(x)) = (x-1)!
//          //=======================

//          normfactor = sqrt(2.0 * K - 1) * exp(lgamma(K)) / sqrt(K * (exp(lgamma(2.0 * K))));

//          f
//          fim    = fftw_malloc(sizeof(fftw_complex) * w  * h);
//          conv   = fftw_malloc(sizeof(fftw_complex) * w  * h);
//          ftmp   = fftw_malloc(sizeof(fftw_complex) * w  * h);
//          finput = fftw_malloc(sizeof(fftw_complex) * w  * h);

//          //initialize input
//          //================

//          for (j = 0; j < h; j++)
//            {
//              for (i = 0; i < w; i ++)
//                {
//                  L1[0][i + w * j] = im(i,j);
//                  finput[i + w * j][0] = im(i,j);
//                  finput[i + w * j][1] = 0.0;
//                }
//            }

//          //scale = 0
//          //=========

          a = w / 2;

//          p = fftw_plan_dft_2d(h, w, finput, fim, -1, FFTW_ESTIMATE); //1st fft
//          fftw_execute(p);
          Array::array2<double> fLP0(h,w,align);
          Array::array2<double> fHP0(h,w,align);
          Array::array2<double> fHPtmp(h,w,align);
          Array::array2<double> fLPtmp(h,w,align);
          //Array::array2<double> tmparall(h,w,align);

          vector<Array::array2<double>> fB(K,Array::array2<double>(h,w,align));
          //vector<Array::array2<double>> tmpar(K,Array::array2<double>(h,w,align));

          Tensor<double,1> temp(h,w);
          genLPfilter(fLP0,   w, h, a/2, a);
          genHPfilter(fHP0,   w, h, a/2, a);
          genHPfilter(fHPtmp, w, h, a/4, a/2);
          genLPfilter(fLPtmp, w, h, a/4, a/2);
          //cout<<fLP0<<endl<<fHP0<<endl<<fHPtmp<<endl<<fLPtmp<<endl;

          for (j=0; j<h;j++)
            for (i=0;i<w;i++)
              temp(j,i) = 255*fLP0(j,i);
          temp.SaveBlock("fLP0.tif");
          for (j=0; j<h;j++)
            for (i=0;i<w;i++)
              temp(j,i) = 255*fHP0(j,i);
          temp.SaveBlock("fHP0.tif");
          for (j=0; j<h;j++)
            for (i=0;i<w;i++)
              temp(j,i) = 255*fLPtmp(j,i);
          temp.SaveBlock("fLPtmp0.tif");
          for (j=0; j<h;j++)
            for (i=0;i<w;i++)
              temp(j,i) = 255*fHPtmp(j,i);
          temp.SaveBlock("fHPtmp0.tif");

          for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
            {
               int k, l, m, n, t, index2, index3;
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
                  fB[t](j,i) = normfactor * pow(2.0 * cos(theta), K-1);
                  //tmpar[t](n,m) = fB[t](j,i) * fLP0(j,i) * fHPtmp(j,i);
               }
            }
//          for (j = 0; j < h; j++)
//            for (i = 0; i < w; i++)
//            {
//              double ttemp;
//              ttemp = 0.0;
//              for (t = 0; t < K; t ++)
//              {
//                ttemp += pow(fB[t](j,i), 2.0);
//              }
//              tmparall(j,i) =   pow(fHP0(j,i), 2.0)
//                                  + (ttemp * pow(fHPtmp(j,i), 2.0)
//                                  + pow(fLPtmp(j,i), 2.0)) * pow(fLP0(j,i), 2.0);

//            }
//          cout<<"tmpar[0]--------------\n";
//          cout<<tmpar[0]<<endl;
//          for (j=0; j<h;j++)
//            for (i=0;i<w;i++)
//              temp(j,i) = tmparall(j,i);
//          double mintmp;// = temp.Min()[0];
//          double maxtmp;// = cv::max(temp);

//         minMaxLoc( temp, &mintmp, &maxtmp );
//         for (j=0; j<h;j++)
//            for (i=0;i<w;i++)
//              temp(j,i) = 255*(temp(j,i)[0]-mintmp)/(maxtmp-mintmp);
//          temp.SaveBlock("fBPcheck0.tif");
//          for (j=0; j<h;j++)
//            for (i=0;i<w;i++)
//              temp(j,i) = 255*tmpar[0](j,i);
//          temp.SaveBlock("fB1new0.tif");
//          for (j=0; j<h;j++)
//            for (i=0;i<w;i++)
//              temp(j,i) = 255*tmpar[1](j,i);
//          temp.SaveBlock("fB2new0.tif");
//          for (j=0; j<h;j++)
//            for (i=0;i<w;i++)
//              temp(j,i) = 255*tmpar[2](j,i);
//          temp.SaveBlock("fB3new0.tif");
//          for (j=0; j<h;j++)
//            for (i=0;i<w;i++)
//              temp(j,i) = 255*tmpar[3](j,i);
//          temp.SaveBlock("fB4new0.tif");
//          //...

//          output (tmparL, w, h, id, im, "fLnew.ppm", 1);
//          output (tmparH, w, h, id, im, "fHnew.ppm", 1);

//          free(tmpar);  	free(tmparL);	free(tmparH);

          for (t = 0; t < K; t ++)
          {

            //fourier2spatialband2BP(K, w, h, fHP0, fB[t],  *A[t], *conv[0], *fim[0], *ftmp[0]);
            fourier2spatialband3(K, w, h, fLP0, fHPtmp, fB[t], *B[t], *conv[0], *fim[0], *ftmp[0]);
            cout<<"B"<<t<<endl;
            cout<<B[t]<<endl;
          }
          fourier2spatialband2(w, h, fLP0, fLPtmp, *L1[0], *conv[0], *fim[0], *ftmp[0]);
          cout<<"conv0"<<endl;
          cout<<*conv[0]<<endl;
//          //subsampling for the next scale
//          //==============================

          w2 = w / 2;
          h2 = h / 2;

//          for (j = 0; j < h2; j ++)
//          {
//            for (i = 0; i < w2; i ++)
//            {
//                L1[1]->operator()(j,i) = L1[0]->operator()(2*j,2*i);
//            }
//          }


//          fftw_free(fim);	fftw_free(conv);	fftw_free(ftmp);	fftw_free(finput);
//          free(fLP0);	free(fHP0);	free(fLPtmp);	free(fHPtmp);	free(fB);

//          //iterative part of the decomposition
//          //===================================
            w1 = w;
            h1 = h;
            for (scale = 1; scale < maxscale; scale ++)
            {
               down2Freq(*conv[scale-1],*fim[scale],h1,w1);
               w1 /=2 ;//current level size
               h1 /=2 ;

               a = w1 / 2;

               //fim    = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
               //conv   = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
               //ftmp   = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
               //finput = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));

//               for (j = 0; j < h1; j++)
//                 for (i = 0; i < w1; i++)
//                 {
//                   finput[scale]->operator()(j,i) = L1[scale-1]->operator()(j,i);
//                 }
//               cout<<"down2 in spatial"<<endl;
//               cout<<*finput[scale]<<endl;

            //   Array::array2<FComplex> temp(h1,w1,align);
            //   fftwpp::fft2d fft2(h1,w1,-1,*finput[scale],temp);
            //   fft2.fft(*finput[scale],temp);
            //   cout<<"use fft-ifft-fft"<<endl;
            //   cout<<temp<<endl;

               cout<<"d2 in freq"<<endl;
               cout<<*fim[scale]<<endl;
//               p = fftw_plan_dft_2d(h1, w1, finput, fim, -1, FFTW_ESTIMATE); fftw_execute(p);
               Array::array2<double> fHPtmp(h1,w1,align);
               Array::array2<double> fLPtmp(h1,w1,align);

               vector<Array::array2<double>> fB(K,Array::array2<double>(h1,w1,align));

//               fLPtmp = malloc(sizeof(double) * w1 * h1);
//               fHPtmp = malloc(sizeof(double) * w1 * h1);

//               fB = malloc(sizeof(double *) * K);
//               for (k = 0; k < K; k ++)
//                 {
//                   fB[k] = malloc(sizeof(double) * w1 * h1);
//                 }

               genHPfilter(fHPtmp, w1, h1, a/4, a/2);
               genLPfilter(fLPtmp, w1, h1, a/4, a/2);

               for (j = 0; j < h1; j++)
                   for (i = 0; i < w1; i++)
                   {
                        int k, l, t;
                        double theta;
                        if (i <  w1/2) {k =     i;}
                        if (i >= w1/2) {k = i - w1;}
                        if (j <  h1/2) {l =     j;}
                        if (j >= h1/2) {l = j - h1;}

                        for (t = 0; t < K; t ++)
                        {
                            theta = atan2(l, k) - t * M_PI / K;
                            fB[t](j,i) = normfactor * pow(2.0 * cos(theta), K-1);
                        }
                    }
               for (t = 0; t < K; t ++)
               {
                   Array::array2<FComplex> tempArr(h1,w1,align);
                   fourier2spatialband2(w1, h1, fHPtmp, fB[t], *B[t+scale*K], *conv[scale], *fim[scale], *ftmp[scale]);
               }
               fourier2spatialband1(w1, h1, fLPtmp, *L1[scale], *conv[scale], *fim[scale], *ftmp[scale]);

//               //subsampling for the next scale
//               //==============================

//               w2 = w / pow(2, scale+1);
//               h2 = h / pow(2, scale+1);

//               for (j = 0; j < h2; j ++)
//                 for (i = 0; i < w2; i ++)
//                 {
//                   L1[scale + 1]->operator()(j,i) = L1[scale]->operator()(2*j,2*i);
//                 }

//               fftw_free(fim);	fftw_free(conv);	fftw_free(ftmp);	fftw_free(finput);
//               free(fLPtmp);	free(fHPtmp);		free(fB);
          }

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
        //int wp = w/2+1;
        cout<<"otf1"<<endl;
        cout<<otf1<<endl;
        cout<<"otf2"<<endl;
        cout<<otf2<<endl;
        cout<<"otf3"<<endl;
        cout<<otf3<<endl;
        cout<<"fim"<<endl;
        cout<<fim<<endl;
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
        cout<<"conv"<<endl;
        cout<<conv<<endl;
        fftwpp::fft2d Backward(h,w,1,conv,BP);
        Backward.fftNormalized(conv,BP);
        cout<<"BP"<<endl;
        cout<<BP<<endl;
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
        cout<<"otf1"<<endl;
        cout<<otf1<<endl;
        cout<<"otf2"<<endl;
        cout<<otf2<<endl;
        cout<<"fim"<<endl;
        cout<<fim<<endl;
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
        cout<<"conv"<<endl;
        cout<<conv<<endl;
        fftwpp::fft2d Backward(h,w,1,conv,BP);
        Backward.fftNormalized(conv,BP);
        cout<<"BP"<<endl;
        cout<<BP<<endl;
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

}

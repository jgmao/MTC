#include <Steerable2.h>
#include <fftw++.h>
#include <Array.h>
namespace metric{
Steerable2::Steerable2()
{

}
Steerable2::Steerable2(c_real_ref im)
{
        this->im = im;
}

Steerable2::~Steerable2()
{

}

int Steerable2::decompose(int maxscale, int K, Mat im, Mat* L0, Mat* L1, Mat* B, Mat* A)
//  int decompose(int maxscale, int K, double *kanaal, int w, int h, double **L0, double **L1, double *H0,
//                  double ***B, double  **A, ImlibData  *id, ImlibImage *im)
  {
          fftwpp::fftw::maxthreads = get_max_threads();
          unsigned int h = im.size().height;
          int w = im.size().width;
          size_t align = sizeof(FComplex);
          Array::array2<double> finput(h,w,align);
          Array::array2<FComplex> fim(h,w/2+1,align);
          fftwpp::rcfft2d plan(h,w,finput,fim);
          fftwpp::crfft2d ifft(h,w,fim,finput);
          for (unsigned int i=0; i< h; i++)
            for (unsigned int j=0; j<w;j++)
              finput(i,j)=im.at<double>(i,j);
          cout<<finput<<endl;

          plan.fft(finput,fim);
          cout<<fim<<endl;
          ifft.fftNormalized(fim,finput);
          cout<<finput<<endl;
          unsigned int nx=h, ny=w;
          unsigned int nyp=ny/2+1;

          Array::array2<double> f(nx,ny,align);
          Array::array2<FComplex> g(nx,nyp,align);

          fftwpp::rcfft2d Forward(nx,ny,f,g);
          fftwpp::crfft2d Backward(nx,ny,g,f);

          for(unsigned int i=0; i < nx; i++)
            for(unsigned int j=0; j < ny; j++)
              f(i,j)=im.at<double>(i,j);//i+j;

          cout << f << endl;

          Forward.fft(f,g);

          cout << g << endl;

          Backward.fftNormalized(g,f);

          cout << f << endl;
//          fftw_plan p;
//          fftw_complex *fim, *conv, *ftmp, *finput;



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

//          a = w / 2;

//          p = fftw_plan_dft_2d(h, w, finput, fim, -1, FFTW_ESTIMATE); //1st fft
//          fftw_execute(p);

//          fLP0     = malloc(sizeof(double) * w * h);
//          fHP0     = malloc(sizeof(double) * w * h);
//          fLPtmp   = malloc(sizeof(double) * w * h);
//          fHPtmp   = malloc(sizeof(double) * w * h);
//          tmparL   = malloc(sizeof(double) * w * h);
//          tmparH   = malloc(sizeof(double) * w * h);
//          tmparall = malloc(sizeof(double) * w * h);

//          fB    = malloc(sizeof(double *) * K);
//          tmpar = malloc(sizeof(double *) * K);


//          for (k = 0; k < K; k ++)
//            {
//              fB[k]    = malloc(sizeof(double) * w * h);
//              tmpar[k] = malloc(sizeof(double) * w * h);
//            }

//          tmparL = malloc(sizeof(double) * w * h);
//          tmparH = malloc(sizeof(double) * w * h);


//          genLPfilter(fLP0,   w, h, a/2, a);
//          genHPfilter(fHP0,   w, h, a/2, a);
//          genHPfilter(fHPtmp, w, h, a/4, a/2);
//          genLPfilter(fLPtmp, w, h, a/4, a/2);

//          for (j = 0; j < h; j++)
//            {
//              for (i = 0; i < w; i++)
//                {
//                  int index = i + w * j;
//                  int k, l, m, n, t, index2, index3;
//                  double theta;
//                  if (i <  w/2) {k =     i;}
//                  if (i >= w/2) {k = i - w;}
//                  if (j <  h/2) {l =     j;}
//                  if (j >= h/2) {l = j - h;}

//                  if (i <  w/2) {m = i + w/2;}
//                  if (i >= w/2) {m = i - w/2;}
//                  if (j <  h/2) {n = j + h/2;}
//                  if (j >= h/2) {n = j - h/2;}

//                  index2 = k + w * l;
//                  index3 = m + w * n;

//                  for (t = 0; t < K; t ++)
//                    {
//                      theta = atan2(l, k) - t * M_PI / K;
//                      fB[t][index] = normfactor * pow(2.0 * cos(theta), K-1);
//                      tmpar[t][index3] = fB[t][index] * fLP0[index] * fHPtmp[index];
//                    }
//               }
//            }
//          for (j = 0; j < h; j++)
//            {
//              for (i = 0; i < w; i++)
//                {
//                  int index = i + w * j;
//                  double ttemp;
//                  ttemp = 0.0;
//                  for (t = 0; t < K; t ++)
//                    {
//                      ttemp += pow(fB[t][index], 2.0);
//                    }
//                  tmparall[index] =   pow(fHP0[index], 2.0)
//                                    + (ttemp * pow(fHPtmp[index], 2.0)
//                                    + pow(fLPtmp[index], 2.0)) * pow(fLP0[index], 2.0);

//                }
//            }
//          output (tmparall, w, h, id, im, "fBPcheck.ppm", 1);


//                      output (tmpar[0], w, h, id, im, "fB1new.ppm", 1);
//          if (K > 1) {output (tmpar[1], w, h, id, im, "fB2new.ppm", 1);}
//          if (K > 2) {output (tmpar[2], w, h, id, im, "fB3new.ppm", 1);}
//          if (K > 3) {output (tmpar[3], w, h, id, im, "fB4new.ppm", 1);}
//          if (K > 4) {output (tmpar[4], w, h, id, im, "fB5new.ppm", 1);}
//          if (K > 5) {output (tmpar[5], w, h, id, im, "fB6new.ppm", 1);}

//          //...

//          output (tmparL, w, h, id, im, "fLnew.ppm", 1);
//          output (tmparH, w, h, id, im, "fHnew.ppm", 1);

//          free(tmpar);  	free(tmparL);	free(tmparH);

//          for (t = 0; t < K; t ++)
//            {
//               fourier2spatialband2BP(K, w, h, fHP0, fB[t],  A[t], conv, fim, ftmp);
//               fourier2spatialband3BP(K, w, h, fLP0, fHPtmp, fB[t], B[t][0], conv, fim, ftmp);
//            }

//          fourier2spatialband2LP(K, w, h, fLP0, fLPtmp, L1[0], conv, fim, ftmp);

//          //subsampling for the next scale
//          //==============================

//          w2 = w / pow(2, 1);
//          h2 = h / pow(2, 1);

//          for (j = 0; j < h2; j ++)
//            {
//              for (i = 0; i < w2; i ++)
//                {
//                   L1[1][i + w2 * j] = L1[0][2 * (i + 2 * w2 * j)];
//                }
//            }

//          fftw_free(fim);	fftw_free(conv);	fftw_free(ftmp);	fftw_free(finput);
//          free(fLP0);	free(fHP0);	free(fLPtmp);	free(fHPtmp);	free(fB);

//          //iterative part of the decomposition
//          //===================================

//          for (scale = 1; scale < maxscale; scale ++)
//            {
//               w1 = w / pow(2, scale);
//               h1 = h / pow(2, scale);

//               printf("Scale = %d\n", scale);

//               a = w1 / 2;

//               fim    = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
//               conv   = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
//               ftmp   = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
//               finput = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));

//               for (j = 0; j < h1; j++)
//                 {
//                   for (i = 0; i < w1; i++)
//                     {
//                        int index = i + w1 * j;
//                        finput[index][0] = L1[scale][index];
//                        finput[index][1] = 0.0;
//                     }
//                 }

//               p = fftw_plan_dft_2d(h1, w1, finput, fim, -1, FFTW_ESTIMATE); fftw_execute(p);

//               fLPtmp = malloc(sizeof(double) * w1 * h1);
//               fHPtmp = malloc(sizeof(double) * w1 * h1);

//               fB = malloc(sizeof(double *) * K);
//               for (k = 0; k < K; k ++)
//                 {
//                   fB[k] = malloc(sizeof(double) * w1 * h1);
//                 }

//               genHPfilter(fHPtmp, w1, h1, a/4, a/2);
//               genLPfilter(fLPtmp, w1, h1, a/4, a/2);

//               for (j = 0; j < h1; j++)
//                 {
//                   for (i = 0; i < w1; i++)
//                     {
//                        int index = i + w1 * j;
//                        int k, l, t;
//                        double theta;
//                        if (i <  w1/2) {k =     i;}
//                        if (i >= w1/2) {k = i - w1;}
//                        if (j <  h1/2) {l =     j;}
//                        if (j >= h1/2) {l = j - h1;}

//                        for (t = 0; t < K; t ++)
//                          {
//                            theta = atan2(l, k) - t * M_PI / K;
//                            fB[t][index] = normfactor * pow(2.0 * cos(theta), K-1);
//                          }
//                    }
//                 }
//               for (t = 0; t < K; t ++)
//                 {
//                   fourier2spatialband2BP(K, w1, h1, fHPtmp, fB[t], B[t][scale], conv, fim, ftmp);
//                 }
//               fourier2spatialband1LP(K, w1, h1, fLPtmp, L1[scale], conv, fim, ftmp);

//               //subsampling for the next scale
//               //==============================

//               w2 = w / pow(2, scale+1);
//               h2 = h / pow(2, scale+1);

//               for (j = 0; j < h2; j ++)
//                 {
//                   for (i = 0; i < w2; i ++)
//                     {
//                        L1[scale + 1][i + w2 * j] = L1[scale][2 * (i + 2 * w2 * j)];
//                     }
//                 }

//               fftw_free(fim);	fftw_free(conv);	fftw_free(ftmp);	fftw_free(finput);
//               free(fLPtmp);	free(fHPtmp);		free(fB);
//            }
          return 1;
   }



}

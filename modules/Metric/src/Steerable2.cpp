#include <Steerable2.h>

namespace metric{

vector<Steerable2::data_type>& Steerable2::buildSCFpyr(Steerable2::c_data_ref im, int nLevel, int nDir, int twidth, bool subsample)
{
  int w, h, w1a, h1a, rr, i, j, verif, maxscale, K;
  double drempel;
  double *imin, **L0, **L1, *H0, *recim, *input2;
  double **Ar, **Ac;

  double ***Br, ***Bc;
  double **mag, **ori;

  //ImlibData  *id, *idout;
  //ImlibImage *im, *imout;
  unsigned char *alfa, *largeim;

  rr = im.size().height/2; //padding;
  maxscale = nLevel;
  K = nDir;
  w = im.size().width;
  h = im.size().height;
  w1a = w + 2 * rr;
  h1a = h + 2 * rr;

  imin   = malloc(sizeof(double) * w1a * h1a);
  H0     = malloc(sizeof(double) * w1a * h1a);
  input2 = malloc(sizeof(double) * w1a * h1a);

  Ar = malloc(sizeof(double *) * K);
  Ac = malloc(sizeof(double *) * K);

  Br = malloc(sizeof(double **) * K);
  Bc = malloc(sizeof(double **) * K);
  for (i = 0; i < K; i++)
  {
    Br[i] = malloc(sizeof(double *) * maxscale);
    Bc[i] = malloc(sizeof(double *) * maxscale);
  }

  for (i = 0; i < K; i++)
  {
    for (j = 0; j < maxscale; j++)
    {
      int w1, h1;
      w1 = w1a / pow(2.0, j);
      h1 = h1a / pow(2.0, j);
      Br[i][j] = malloc(sizeof(double) * w1 * h1);
      Bc[i][j] = malloc(sizeof(double) * w1 * h1);
    }
    Ar[i] = malloc(sizeof(double) * w1a * h1a);
    Ac[i] = malloc(sizeof(double) * w1a * h1a);
  }
  recim  = malloc(sizeof(double) * w1a * h1a);
  largeim= malloc(sizeof(int) * w1a * h1a);

  mag   = malloc(sizeof(double *) *  maxscale);
  ori   = malloc(sizeof(double *) *  maxscale);
  L0    = malloc(sizeof(double *) *  maxscale);
  L1    = malloc(sizeof(double *) * (maxscale + 1));

  for (i = 0; i < maxscale; i ++)
  {
    int w1, h1;
    w1 = w1a / pow(2.0, i);
    h1 = h1a / pow(2.0, i);
    mag[i]  = malloc(sizeof(double) * w1 * h1);
    ori[i]  = malloc(sizeof(double) * w1 * h1);
    L0[i]   = malloc(sizeof(double) * w1 * h1);
  }
  for (i = 0; i < maxscale + 1; i ++)
  {
    int w1, h1;
    w1 = w1a / pow(2.0, i);
    h1 = h1a / pow(2.0, i);
    L1[i]    = malloc(sizeof(double) * w1 * h1);
  }
  kanaal(imin, w, h, im.data);
  expand (imin, input2, w, h, rr);
  verif = decompose(maxscale, K, input2, w1a, h1a, L0, L1, H0,
                              Br, Bc, Ar, Ac, idout, imout);

  free(H0);	free(imin);	free(recim);
  free(Ar);	free(Ac);	free(Br);	free(Bc);
  free(L0);	free(L1);


}


double Steerable2::steera(double theta, double *B1, double *B2, double *B3, double *B4, int w, int h, int x, int y)
 {
    double k1a, k1b, k1c, k1d, outsteera;
    int index;
    k1a = (2.0 * cos((theta - 0.0 * M_PI / 4.0)) + 2.0 * cos(3.0 * (theta - 0.0 * M_PI / 4.0))) / 4.0;
    k1b = (2.0 * cos((theta - 1.0 * M_PI / 4.0)) + 2.0 * cos(3.0 * (theta - 1.0 * M_PI / 4.0))) / 4.0;
    k1c = (2.0 * cos((theta - 2.0 * M_PI / 4.0)) + 2.0 * cos(3.0 * (theta - 2.0 * M_PI / 4.0))) / 4.0;
    k1d = (2.0 * cos((theta - 3.0 * M_PI / 4.0)) + 2.0 * cos(3.0 * (theta - 3.0 * M_PI / 4.0))) / 4.0;
    index = x + w * y;
    outsteera = k1a * B1[index] + k1b * B2[index] + k1c * B3[index] + k1d * B4[index];
    return outsteera;
 }

double Steerable2::steerb(double theta, double *hB1, double *hB2, double *hB3, double *hB4, double *hB5,
       int w, int h, int x, int y)
 {
    double k1a, k1b, k1c, k1d, k1e, outsteerb;
    int index;
    k1a = (1.0 + 2.0 * cos(2.0 * (theta - 0.0 * M_PI / 5.0)) + 2.0 * cos(4.0 * (theta - 0.0 * M_PI / 5.0))) / 5.0;
    k1b = (1.0 + 2.0 * cos(2.0 * (theta - 1.0 * M_PI / 5.0)) + 2.0 * cos(4.0 * (theta - 1.0 * M_PI / 5.0))) / 5.0;
    k1c = (1.0 + 2.0 * cos(2.0 * (theta - 2.0 * M_PI / 5.0)) + 2.0 * cos(4.0 * (theta - 2.0 * M_PI / 5.0))) / 5.0;
    k1d = (1.0 + 2.0 * cos(2.0 * (theta - 3.0 * M_PI / 5.0)) + 2.0 * cos(4.0 * (theta - 3.0 * M_PI / 5.0))) / 5.0;
    k1e = (1.0 + 2.0 * cos(2.0 * (theta - 4.0 * M_PI / 5.0)) + 2.0 * cos(4.0 * (theta - 4.0 * M_PI / 5.0))) / 5.0;
    index = x + w * y;
    outsteerb = k1a * hB1[index] + k1b * hB2[index] + k1c * hB3[index] + k1d * hB4[index] + k1e * hB5[index];
    return outsteerb;
 }

//Generate bandpass image from FT with one transfer function mask, real output
//============================================================================

int Steerable2::fourier2spatialband1(int w, int h, double *otf1, double *BP, fftw_complex *conv, fftw_complex *fim, fftw_complex *ftmp)
  {
     fftw_plan p;
     int i, j;
     for (j = 0; j < h; j++)
       {
         for (i = 0; i < w; i++)
           {
             int index = i + w * j;
             conv[index][0] = fim[index][0] / (w * h) * otf1[index];
             conv[index][1] = fim[index][1] / (w * h) * otf1[index];
           }
       }
     p = fftw_plan_dft_2d(h, w, conv, ftmp, 1, FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);
     for (j = 0; j < h; j++) {for (i = 0; i < w; i++) {BP[i + w * j] = ftmp[i + w * j][0];}}
     return 1;
  }

//Generate bandpass image from FT with two transfer function masks, complex output
//================================================================================

int Steerable2::fourier2spatialband2(int w, int h, double *otf1, double *otf2, double *BPr, double *BPc, fftw_complex *conv, fftw_complex *fim, fftw_complex *ftmp)
  {
     fftw_plan p;
     int i, j;
     for (j = 0; j < h; j++)
       {
         for (i = 0; i < w; i++)
           {
             int index = i + w * j;
             conv[index][0] = fim[index][0] / (w * h) * otf1[index] * otf2[index];
             conv[index][1] = fim[index][1] / (w * h) * otf1[index] * otf2[index];
           }
       }
     p = fftw_plan_dft_2d(h, w, conv, ftmp, 1, FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);
     for (j = 0; j < h; j++) {for (i = 0; i < w; i++) {BPr[i + w * j] = ftmp[i + w * j][0];}}
     for (j = 0; j < h; j++) {for (i = 0; i < w; i++) {BPc[i + w * j] = ftmp[i + w * j][1];}}
     return 1;
  }
//Generate bandpass image from FT with two transfer function masks, real output
//=============================================================================

int Steerable2::fourier2spatialband2a(int w, int h, double *otf1, double *otf2, double *BP, fftw_complex *conv, fftw_complex *fim, fftw_complex *ftmp)
  {
     fftw_plan p;
     int i, j;
     double maxr, minr, maxc, minc;
     maxr = 0.0; minr = 1000; maxc = 0.0; minc = 1000.0;
     for (j = 0; j < h; j++)
       {
         for (i = 0; i < w; i++)
           {
             int index = i + w * j;
             conv[index][0] = fim[index][0] / (w * h) * otf1[index] * otf2[index];
             conv[index][1] = fim[index][1] / (w * h) * otf1[index] * otf2[index];
           }
       }
     p = fftw_plan_dft_2d(h, w, conv, ftmp, 1, FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);
     for (j = 0; j < h; j++) {for (i = 0; i < w; i++) {BP[i + w * j] = ftmp[i + w * j][0];}}
     return 1;
  }

//Generate bandpass image from FT with three transfer function masks, complex output
//==================================================================================

int Steerable2::fourier2spatialband3(int w, int h, double *otf1, double *otf2, double *otf3, double *BPr, double *BPc, fftw_complex *conv, fftw_complex *fim, fftw_complex *ftmp)
  {
     fftw_plan p;
     int i, j;
     for (j = 0; j < h; j++)
       {
         for (i = 0; i < w; i++)
           {
             int index = i + w * j;
             conv[index][0] = fim[index][0] / (w * h) * otf1[index] * otf2[index] * otf3[index];
             conv[index][1] = fim[index][1] / (w * h) * otf1[index] * otf2[index] * otf3[index];
           }
       }
     p = fftw_plan_dft_2d(h, w, conv, ftmp, 1, FFTW_ESTIMATE);
     fftw_execute(p);
     fftw_destroy_plan(p);
     for (j = 0; j < h; j++) {for (i = 0; i < w; i++) {BPr[i + w * j] = ftmp[i + w * j][0];}}
     for (j = 0; j < h; j++) {for (i = 0; i < w; i++) {BPc[i + w * j] = ftmp[i + w * j][1];}}
     return 1;
  }

//Reconstruction substep, input real and complex part of oriented
//subbands for fourier transform for adding the different subbands
//================================================================

int Steerable2::reconststep(fftw_complex *fourtmp, fftw_complex *fourBP, double *BPr, double *BPc, int w, int h)
 {
   fftw_plan p;
   int i, j;
   for (j = 0; j < h; j ++)
     {
       for (i = 0; i < w; i++)
         {
           fourtmp[i + w * j][0] = BPr[i + w * j];
           fourtmp[i + w * j][1] = BPc[i + w * j];
         }
     }
   p = fftw_plan_dft_2d(h, w, fourtmp, fourBP, -1, FFTW_ESTIMATE);
   fftw_execute(p);
   fftw_destroy_plan(p);
   return 1;
 }

//Reconstruction substep, input with only real part subbands (LP
//and HP) subbands for fourier transform for adding the different subbands
//========================================================================

int Steerable2::reconststepa(fftw_complex *fourtmp, fftw_complex *fourBP, double *BP, int w, int h)
 {
   fftw_plan p;
   int i, j;
   for (j = 0; j < h; j ++)
     {
       for (i = 0; i < w; i++)
         {
           fourtmp[i + w * j][0] = BP[i + w * j];
           fourtmp[i + w * j][1] = 0.0;
         }
     }
   p = fftw_plan_dft_2d(h, w, fourtmp, fourBP, -1, FFTW_ESTIMATE);
   fftw_execute(p);
   fftw_destroy_plan(p);
   return 1;
 }

int Steerable2::genHPfilter(double *HP, int w, int h, double x1, double x2)
 {
    double r, x;
    int i, j, k, l;

    for (j = 0; j < h; j++)
      {
         for (i = 0; i < w; i++)
           {
             int index = i + w * j;

             if (i <= w/2) {k = i;}
             if (i >  w/2) {k = i - w;}

             if (j <= h/2) {l = j;}
             if (j >  h/2) {l = j - h;}

             r = sqrt(k*k + l*l);

             if (r < x1) {HP[index] = 0.0;}
             if (r > x2) {HP[index] = 1.0;}
             if ((r >= x1) && (r <= x2))
               {
                 x = M_PI / 4.0 * (1.0 + (r - x1) / (x2 - x1));
                 HP[index] = cos(M_PI / 2.0 * (log (x * 2.0 / M_PI) / log(2.0)));
               }
         }
     }
   return 1;
 }


int Steerable2::genLPfilter(double *LP0, int w, int h, double x1, double x2)
 {
    double r, x;
    int i, j, k, l;

    for (j = 0; j < h; j++)
      {
         for (i = 0; i < w; i++)
           {
             int index = i + w * j;

             if (i <= w/2) {k = i;}
             if (i >  w/2) {k = i - w;}

             if (j <= h/2) {l = j;}
             if (j >  h/2) {l = j - h;}

             r = sqrt(k*k + l*l);

             if (r < x1) {LP0[index] = 1.0;}
             if (r > x2) {LP0[index] = 0.0;}

             if ((r >= x1) && (r <= x2))
               {
                 x = M_PI / 4.0 * (1.0 + (r - x1) / (x2 - x1));
                 LP0[index] = cos(M_PI / 2.0 * (log(x * 4.0 / M_PI)/ log(2.0)));
               }
           }
     }
   return 1;
 }
int Steerable2::decompose(int maxscale, int K, double *kanaal, int w, int h, double **L0, double **L1, double *H0,
               double ***Br, double ***Bc, double  **Ar, double  **Ac, ImlibData  *id, ImlibImage *im)

 {
       fftw_plan p;
       fftw_complex *fim, *conv, *ftmp, *finput;

       int i, j, k, t, w1, h1, w2, h2, scale;
       double *fLP0, *fHP0, *fLPtmp, *fHPtmp;
       double **fB, *fBP, **tmpar, *tmparL, *tmparH;
       double a, kfact, kfact2, normfactor;
       double normC;
       normC = sqrt(2.0 * K -1);

       fim    = fftw_malloc(sizeof(fftw_complex) * w  * h);
       conv   = fftw_malloc(sizeof(fftw_complex) * w  * h);
       ftmp   = fftw_malloc(sizeof(fftw_complex) * w  * h);
       finput = fftw_malloc(sizeof(fftw_complex) * w  * h);

       //initialize input
       //================

       for (j = 0; j < h; j++)
         {
           for (i = 0; i < w; i ++)
             {
               L1[0][i + w * j] = kanaal[i + w * j];
               finput[i + w * j][0] = kanaal[i + w * j];
               finput[i + w * j][1] = 0.0;
             }
         }

       //scale = 0
       //=========

       a = w / 2;

       p = fftw_plan_dft_2d(h, w, finput, fim, -1, FFTW_ESTIMATE);
       fftw_execute(p);

       fLP0   = malloc(sizeof(double) * w * h);
       fHP0   = malloc(sizeof(double) * w * h);
       fLPtmp = malloc(sizeof(double) * w * h);
       fHPtmp = malloc(sizeof(double) * w * h);
       tmparL = malloc(sizeof(double) * w * h);
       tmparH = malloc(sizeof(double) * w * h);

       fB    = malloc(sizeof(double *) * K);
       tmpar = malloc(sizeof(double *) * K);


       for (k = 0; k < K; k ++)
         {
           fB[k]    = malloc(sizeof(double) * w * h);
           tmpar[k] = malloc(sizeof(double) * w * h);
         }

       tmparL = malloc(sizeof(double) * w * h);
       tmparH = malloc(sizeof(double) * w * h);

       //exp(lgamma(x)) = (x-1)!
       //=======================

       kfact  = exp(lgamma(K));
       kfact2 = exp(lgamma(2.0 * K));
       normfactor = normC * kfact / sqrt(K * (kfact2));

       genLPfilter(fLP0,   w, h, a/2, a);
       genHPfilter(fHP0,   w, h, a/2, a);
       genHPfilter(fHPtmp, w, h, a/4, a/2);
       genLPfilter(fLPtmp, w, h, a/4, a/2);

       for (j = 0; j < h; j++)
         {
           for (i = 0; i < w; i++)
             {
               int index = i + w * j;
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

               index2 = k + w * l;
               index3 = m + w * n;

               for (t = 0; t < K; t ++)
                 {
                   theta = atan2(l, k) - t * M_PI / K;
                   if ((fabs(theta) < M_PI / 2.0) || (fabs(theta) > 3.0 * M_PI / 2.0))
                     {
                       fB[t][index] = normfactor * pow(2.0 * cos(theta), K-1);
                     }
                   else
                     {
                       fB[t][index] = 0.0;
                     }
                   tmpar[t][index3] = fB[t][index] * fLP0[index] * fHPtmp[index];
                 }
            }
         }

                   output (tmpar[0], w, h, id, im, "fB1new.ppm", 1);
       if (K > 1) {output (tmpar[1], w, h, id, im, "fB2new.ppm", 1);}
       if (K > 2) {output (tmpar[2], w, h, id, im, "fB3new.ppm", 1);}
       if (K > 3) {output (tmpar[3], w, h, id, im, "fB4new.ppm", 1);}
       output (tmparL, w, h, id, im, "fLnew.ppm", 1);
       output (tmparH, w, h, id, im, "fHnew.ppm", 1);

       free(tmpar);  	free(tmparL);	free(tmparH);

       for (t = 0; t < K; t ++)
         {
            fourier2spatialband2(w, h, fHP0, fB[t],  Ar[t], Ac[t], conv, fim, ftmp);
            fourier2spatialband3(w, h, fLP0, fHPtmp, fB[t], Br[t][0], Bc[t][0], conv, fim, ftmp);
         }

       fourier2spatialband2a(w, h, fLP0, fLPtmp, L1[0], conv, fim, ftmp);

       //subsampling for the next scale
       //==============================

       w2 = w / pow(2, 1);
       h2 = h / pow(2, 1);

       for (j = 0; j < h2; j ++)
         {
           for (i = 0; i < w2; i ++)
             {
                L1[1][i + w2 * j] = L1[0][2 * (i + 2 * w2 * j)];
             }
         }

       fftw_free(fim);	fftw_free(conv);	fftw_free(ftmp);	fftw_free(finput);
       free(fLP0);	free(fHP0);	free(fLPtmp);	free(fHPtmp);	free(fB);

       //iterative part of the decomposition
       //===================================

       for (scale = 1; scale < maxscale; scale ++)
         {
            w1 = w / pow(2, scale);
            h1 = h / pow(2, scale);

            printf("Scale = %d\n", scale);

            a = w1 / 2;

            fim    = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
            conv   = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
            ftmp   = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
            finput = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));

            for (j = 0; j < h1; j++)
              {
                for (i = 0; i < w1; i++)
                  {
                     int index = i + w1 * j;
                     finput[index][0] = L1[scale][index];
                     finput[index][1] = 0.0;
                  }
              }

            p = fftw_plan_dft_2d(h1, w1, finput, fim, -1, FFTW_ESTIMATE); fftw_execute(p);

            fLPtmp = malloc(sizeof(double) * w1 * h1);
            fHPtmp = malloc(sizeof(double) * w1 * h1);

            fB = malloc(sizeof(double *) * K);
            for (k = 0; k < K; k ++)
              {
                fB[k] = malloc(sizeof(double) * w1 * h1);
              }

	    kfact  = exp(lgamma(K));
	    kfact2 = exp(lgamma(2.0 * K));
	    normfactor = normC * kfact / sqrt(K * (kfact2));

	    genHPfilter(fHPtmp, w1, h1, a/4, a/2);
	    genLPfilter(fLPtmp, w1, h1, a/4, a/2);

            for (j = 0; j < h1; j++)
              {
                for (i = 0; i < w1; i++)
                  {
                     int index = i + w1 * j;
                     int k, l, t;
                     double theta;
                     if (i <  w1/2) {k =     i;}
                     if (i >= w1/2) {k = i - w1;}
                     if (j <  h1/2) {l =     j;}
                     if (j >= h1/2) {l = j - h1;}

                     for (t = 0; t < K; t ++)
                       {
                         theta = atan2(l, k) - t * M_PI / K;
                         if ((fabs(theta) < M_PI / 2.0) || (fabs(theta) > 3.0 * M_PI / 2.0))
                           {
                             fB[t][index] = normfactor * pow(2.0 * cos(theta), K-1);
                           }
                         else
                           {
                             fB[t][index] = 0.0;
                           }
                       }
                 }
              }
            for (t = 0; t < K; t ++)
              {
                fourier2spatialband2(w1, h1, fHPtmp, fB[t],  Br[t][scale], Bc[t][scale], conv, fim, ftmp);
              }
            fourier2spatialband1(w1, h1, fLPtmp, L1[scale], conv, fim, ftmp);

            //subsampling for the next scale
            //==============================

            w2 = w / pow(2, scale+1);
            h2 = h / pow(2, scale+1);

            for (j = 0; j < h2; j ++)
              {
                for (i = 0; i < w2; i ++)
                  {
                     L1[scale + 1][i + w2 * j] = L1[scale][2 * (i + 2 * w2 * j)];
                  }
              }

            fftw_free(fim);	fftw_free(conv);	fftw_free(ftmp);	fftw_free(finput);
            free(fLPtmp);	free(fHPtmp);		free(fB);
         }
       return 1;
}
int Steerable2::reconstruct(int maxscale, int K, double *kanaal, int w, int h, double **L0, double **L1, double *H0,
	       double ***Br, double ***Bc, double  **Ar, double  **Ac,
	       ImlibData  *id, ImlibImage *im)
 {
       fftw_plan  p, pinv;
       fftw_complex *fourL0, *fourL1, *fourH0, *fourreconst;
       fftw_complex *fourtmp;
       fftw_complex **fourB, **fourA;

       double *fLP0, *fHP0, *fLPtmp, *fHPtmp, **fB;
       int i, j, k, l, t, w1, h1, w2, h2, scale;
       double a, kfact, kfact2, normfactor;
       double normC;
       normC = sqrt(2.0 * K -1);

       //recursive part of the reconstruction
       //====================================

       for (scale = maxscale - 1; scale >= 1; scale --)
         {
           w1 = w / pow(2, scale);
           h1 = h / pow(2, scale);

           w2 = w / pow(2, scale + 1);
           h2 = h / pow(2, scale + 1);

           a = w1 / 2;
           kfact  = exp(lgamma(K));
           kfact2 = exp(lgamma(2.0 * K));
           normfactor = normC * kfact / sqrt(K * (kfact2));

           fB = malloc(sizeof(double *) * K);
           for (t = 0; t < K; t ++)
             {
               fB[t] = malloc(sizeof(double) * w1 * h1);
             }

           fLPtmp = malloc(sizeof(double) * w1 * h1);
           fHPtmp = malloc(sizeof(double) * w1 * h1);

	   genHPfilter(fHPtmp, w1, h1, a/4, a/2);
	   genLPfilter(fLPtmp, w1, h1, a/4, a/2);

	   // upsampling, and multiplying the energy with 4
	   // to compensate for energy loss by replacing
	   // three quart of the pixels with 0
	   //==============================================

	   for (j = 0; j < h1; j++)
	     {
	       for (i = 0; i < w1; i++)
		 {
		   L1[scale][i + w1 * j] = 0.0;
		 }
	     }
	   for (j = 0; j < h2; j++)
	     {
	       for (i = 0; i < w2; i++)
		 {
		   L1[scale][2 * (i + 2 * w2 * j)] = 4.0 * L1[scale + 1][i + w2 * j];
		 }
	     }

	   for (j = 0; j < h1; j++)
	     {
	       for (i = 0; i < w1; i++)
		 {
		    int index = i + w1 * j;
		    double theta;
		    if (i <  w1/2) {k =      i;}
		    if (i >= w1/2) {k = i - w1;}
		    if (j <  h1/2) {l =      j;}
		    if (j >= h1/2) {l = j - h1;}
		    for (t = 0; t < K; t ++)
		      {
			theta = atan2(l, k) - t * M_PI / K;
			if ((fabs(theta) < M_PI / 2.0) || (fabs(theta) > 3.0 * M_PI / 2.0))
			  {
			    fB[t][index] = normfactor * pow(2.0 * cos(theta), K-1);
			  }
			else
			  {
			    fB[t][index] = 0.0;
			  }
		      }

               }
             }
           fourB = malloc(sizeof(fftw_complex *) * K);
           for (t = 0; t < K; t ++)
             {
               fourB[t] = malloc(sizeof(fftw_complex) * w1 * h1);
             }

	   fourL0      = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
	   fourL1      = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
	   fourtmp     = fftw_malloc(sizeof(fftw_complex) * w1 * h1);

	   fourreconst = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
	   for (t = 0; t < K; t ++)
	     {
	       reconststep(fourtmp, fourB[t], Br[t][scale], Bc[t][scale], w1, h1);
	     }
	   reconststepa(fourtmp, fourL1, L1[scale], w1, h1);

	   for (j = 0; j < h1; j++)
	     {
	       for (i = 0; i < w1; i++)
		 {
		   int index1 = i + w1 * j;
		   int index2 = index1;//(w1-i-1) + w1 * (h1-j-1);

		   double rectmp_re, rectmp_im;

		   rectmp_re = 0.0;
		   rectmp_im = 0.0;

		   for (t = 0; t < K; t++)
		     {
		       rectmp_re += fB[t][index2] * fourB[t][index1][0];
		       rectmp_im += fB[t][index2] * fourB[t][index1][1];
		     }

		   fourreconst[index1][0] = (    rectmp_re * 2.0 * fHPtmp[index2]
					       + fLPtmp[index2]  * fourL1[index1][0]  )  / (w1 * h1);

		   fourreconst[index1][1] = (    rectmp_im * 2.0 * fHPtmp[index2]
					       + fLPtmp[index2]  * fourL1[index1][1]  )  / (w1 * h1);
		 }
	     }

	   pinv = fftw_plan_dft_2d(h1, w1, fourreconst, fourtmp, 1, FFTW_ESTIMATE);
	   fftw_execute(pinv);
	   fftw_destroy_plan(pinv);

	   for (j = 0; j < h1; j ++)
	     {
	       for (i = 0; i < w1; i++)
		 {
		   L1[scale][i + w1 * j] = fourtmp[i + w1 * j][0];
		 }
	     }

	   free(fLPtmp);	free(fHPtmp);
	   free(fB);

           fftw_free(fourB);
           fftw_free(fourL0);	fftw_free(fourL1);
           fftw_free(fourreconst);	fftw_free(fourtmp);
        }

   //reconstruction of scale 0
   //=========================

   scale = 0;

   w1 = w / pow(2, scale);
   h1 = h / pow(2, scale);

   w2 = w / pow(2, scale + 1);
   h2 = h / pow(2, scale + 1);

   a = w1 / 2;
   kfact  = exp(lgamma(K));
   kfact2 = exp(lgamma(2.0 * K));
   normfactor = normC * kfact / sqrt(K * (kfact2));
   fB = malloc(sizeof(double *) * K);
   for (t = 0; t < K; t ++)
     {
       fB[t] = malloc(sizeof(double) * w1 * h1);
     }

   fLPtmp = malloc(sizeof(double) * w1 * h1);
   fHPtmp = malloc(sizeof(double) * w1 * h1);
   fLP0   = malloc(sizeof(double) * w1 * h1);
   fHP0   = malloc(sizeof(double) * w1 * h1);

   genHPfilter(fHPtmp, w1, h1, a/4, a/2);
   genLPfilter(fLPtmp, w1, h1, a/4, a/2);
   genHPfilter(fHP0,   w1, h1, a/2, a);
   genLPfilter(fLP0,   w1, h1, a/2, a);

   //upsampling, and multiplying the energy with 4
   //to compensate for energy loss by replacing
   //three quart of the pixels with 0
   //===============================================

   for (j = 0; j < h1; j++)
     {
       for (i = 0; i < w1; i++)
         {
           L1[scale][i + w1 * j] = 0.0;
         }
     }
   for (j = 0; j < h2; j++)
     {
       for (i = 0; i < w2; i++)
         {
           L1[scale][2 * (i + 2 * w2 * j)] = 4.0 * L1[scale + 1][i + w2 * j];
         }
     }

   for (j = 0; j < h1; j++)
     {
       for (i = 0; i < w1; i++)
         {
            int index = i + w1 * j;
            double theta;
            if (i <  w1/2) {k =      i;}
            if (i >= w1/2) {k = i - w1;}
            if (j <  h1/2) {l =      j;}
            if (j >= h1/2) {l = j - h1;}
            for (t = 0; t < K; t ++)
              {
                theta = atan2(l, k) - t * M_PI / K;
                if ((fabs(theta) < M_PI / 2.0) || (fabs(theta) > 3.0 * M_PI / 2.0))
                  {
                    fB[t][index] = normfactor * pow(2.0 * cos(theta), K-1);
                  }
                else
                  {
                    fB[t][index] = 0.0;
                  }
              }
        }
     }

   fourB = malloc(sizeof(fftw_complex *) * K);
   fourA = malloc(sizeof(fftw_complex *) * K);

   for (t = 0; t < K; t ++)
     {
       fourB[t] = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
       fourA[t] = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
     }

   fourL0      = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
   fourH0      = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
   fourL1      = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
   fourreconst = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
   fourtmp     = fftw_malloc(sizeof(fftw_complex) * w1 * h1);
   for (t = 0; t < K; t ++)
     {
       reconststep(fourtmp, fourB[t], Br[t][scale], Bc[t][scale], w1, h1);
       reconststep(fourtmp, fourA[t], Ar[t], Ac[t], w1, h1);
     }

   for (j = 0; j < h1; j ++)
     {
       for (i = 0; i < w1; i++)
         {
           fourtmp[i + w1 * j][0] = L1[scale][i + w1 * j]; fourtmp[i + w1 * j][1] = 0.0;
         }
     }
   p = fftw_plan_dft_2d(h1, w1, fourtmp, fourL1, -1, FFTW_ESTIMATE);
   fftw_execute(p);

   for (j = 0; j < h1; j++)
     {
       for (i = 0; i < w1; i++)
         {
           int index1 = i + w1 * j;
           int index2 = index1;//(w1-i-1) + w1 * (h1-j-1);
           double Arectmp_re, Arectmp_im;
           double Brectmp_re, Brectmp_im;

           Arectmp_re = 0.0;
           Arectmp_im = 0.0;
           Brectmp_re = 0.0;
           Brectmp_im = 0.0;

           for (t = 0; t < K; t++)
             {

               Brectmp_re += fB[t][index2] * fourB[t][index1][0];
               Brectmp_im += fB[t][index2] * fourB[t][index1][1];
               Arectmp_re += fB[t][index2] * fourA[t][index1][0];
               Arectmp_im += fB[t][index2] * fourA[t][index1][1];
             }

           fourreconst[index1][0] = (  (   Brectmp_re * 2.0 * fHPtmp[index2]
                                         + fLPtmp[index2]  * fourL1[index1][0]) * fLP0[index2]
                                      + Arectmp_re * 2.0 * fHP0[index2]
                                    )   / (w1 * h1);

           fourreconst[index1][1] = (  (  (Brectmp_im) * 2.0 * fHPtmp[index2]
                                         + fLPtmp[index2]  * fourL1[index1][1]) * fLP0[index2]
                                      + (Arectmp_im) * 2.0 * fHP0[index2]
                                    ) / (w1 * h1);

         }
     }

   pinv = fftw_plan_dft_2d(h1, w1, fourreconst, fourtmp, 1, FFTW_ESTIMATE);
   fftw_execute(pinv);
   fftw_destroy_plan(pinv);

   for (j = 0; j < h1; j ++)
     {
       for (i = 0; i < w1; i++)
         {
           kanaal[i + w1 * j] = fourtmp[i + w1 * j][0];
         }
     }

   free(fLP0);		free(fHP0);	free(fLPtmp);	free(fHPtmp);
   free(fB);

   fftw_free(fourB);
   fftw_free(fourA);
   fftw_free(fourL0);	fftw_free(fourH0);  	fftw_free(fourL1);  	fftw_free(fourtmp);   	fftw_free(fourreconst);

   return 1;
}

int Steerable2::kanaal(double *out, int w, int h, ImlibImage *im)
  {
        int i, j, index;

	for(i = 0; i < w; i++)
	  {
	      for (j = 0; j < h; j++)
		{
		   index = (i + w * j) * 3;
		   out[i + w * j] = (double) (im->rgb_data[index]);
		}
	}
	return 1;
  }

int Steerable2::decompose(int maxscale, int K, double *kanaal, int w, int h, double **L0, double **L1, double *H0,
                double ***Br, double ***Bc, double  **Ar, double  **Ac, ImlibData  *id, ImlibImage *im)

  {
        fftw_plan p;
        fftw_complex *fim, *conv, *ftmp, *finput;

	int i, j, k, t, w1, h1, w2, h2, scale;
	double *fLP0, *fHP0, *fLPtmp, *fHPtmp;
	double **fB, *fBP, **tmpar, *tmparL, *tmparH;
	double a, kfact, kfact2, normfactor;
	double normC;
	normC = sqrt(2.0 * K -1);

	fim    = fftw_malloc(sizeof(fftw_complex) * w  * h);
	conv   = fftw_malloc(sizeof(fftw_complex) * w  * h);
	ftmp   = fftw_malloc(sizeof(fftw_complex) * w  * h);
	finput = fftw_malloc(sizeof(fftw_complex) * w  * h);

	//initialize input
	//================

	for (j = 0; j < h; j++)
	  {
	    for (i = 0; i < w; i ++)
	      {
		L1[0][i + w * j] = kanaal[i + w * j];
		finput[i + w * j][0] = kanaal[i + w * j];
		finput[i + w * j][1] = 0.0;
	      }
	  }

	//scale = 0
	//=========

	a = w / 2;

	p = fftw_plan_dft_2d(h, w, finput, fim, -1, FFTW_ESTIMATE);
	fftw_execute(p);

	fLP0   = malloc(sizeof(double) * w * h);
	fHP0   = malloc(sizeof(double) * w * h);
	fLPtmp = malloc(sizeof(double) * w * h);
	fHPtmp = malloc(sizeof(double) * w * h);
	tmparL = malloc(sizeof(double) * w * h);
	tmparH = malloc(sizeof(double) * w * h);

	fB    = malloc(sizeof(double *) * K);
	tmpar = malloc(sizeof(double *) * K);


	for (k = 0; k < K; k ++)
	  {
	    fB[k]    = malloc(sizeof(double) * w * h);
	    tmpar[k] = malloc(sizeof(double) * w * h);
	  }

	tmparL = malloc(sizeof(double) * w * h);
	tmparH = malloc(sizeof(double) * w * h);

	//exp(lgamma(x)) = (x-1)!
	//=======================

	kfact  = exp(lgamma(K));
	kfact2 = exp(lgamma(2.0 * K));
	normfactor = normC * kfact / sqrt(K * (kfact2));

	genLPfilter(fLP0,   w, h, a/2, a);
	genHPfilter(fHP0,   w, h, a/2, a);
	genHPfilter(fHPtmp, w, h, a/4, a/2);
	genLPfilter(fLPtmp, w, h, a/4, a/2);

	for (j = 0; j < h; j++)
	  {
	    for (i = 0; i < w; i++)
	      {
		int index = i + w * j;
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

		index2 = k + w * l;
		index3 = m + w * n;

		for (t = 0; t < K; t ++)
		  {
		    theta = atan2(l, k) - t * M_PI / K;
		    if ((fabs(theta) < M_PI / 2.0) || (fabs(theta) > 3.0 * M_PI / 2.0))
		      {
			fB[t][index] = normfactor * pow(2.0 * cos(theta), K-1);
		      }
		    else
		      {
			fB[t][index] = 0.0;
		      }
		    tmpar[t][index3] = fB[t][index] * fLP0[index] * fHPtmp[index];
		  }
	     }
	  }

		    output (tmpar[0], w, h, id, im, "fB1new.ppm", 1);
	if (K > 1) {output (tmpar[1], w, h, id, im, "fB2new.ppm", 1);}
	if (K > 2) {output (tmpar[2], w, h, id, im, "fB3new.ppm", 1);}
	if (K > 3) {output (tmpar[3], w, h, id, im, "fB4new.ppm", 1);}
	output (tmparL, w, h, id, im, "fLnew.ppm", 1);
	output (tmparH, w, h, id, im, "fHnew.ppm", 1);

	free(tmpar);  	free(tmparL);	free(tmparH);

	for (t = 0; t < K; t ++)
	  {
	     fourier2spatialband2(w, h, fHP0, fB[t],  Ar[t], Ac[t], conv, fim, ftmp);
	     fourier2spatialband3(w, h, fLP0, fHPtmp, fB[t], Br[t][0], Bc[t][0], conv, fim, ftmp);
	  }

	fourier2spatialband2a(w, h, fLP0, fLPtmp, L1[0], conv, fim, ftmp);

	//subsampling for the next scale
	//==============================

	w2 = w / pow(2, 1);
	h2 = h / pow(2, 1);

	for (j = 0; j < h2; j ++)
	  {
	    for (i = 0; i < w2; i ++)
	      {
		 L1[1][i + w2 * j] = L1[0][2 * (i + 2 * w2 * j)];
	      }
	  }

	fftw_free(fim);	fftw_free(conv);	fftw_free(ftmp);	fftw_free(finput);
	free(fLP0);	free(fHP0);	free(fLPtmp);	free(fHPtmp);	free(fB);

	//iterative part of the decomposition
	//===================================

	for (scale = 1; scale < maxscale; scale ++)
	  {
	     w1 = w / pow(2, scale);
	     h1 = h / pow(2, scale);

	     printf("Scale = %d\n", scale);

	     a = w1 / 2;

	     fim    = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
	     conv   = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
	     ftmp   = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));
	     finput = fftw_malloc(sizeof(fftw_complex) * (w1  * h1));

	     for (j = 0; j < h1; j++)
	       {
		 for (i = 0; i < w1; i++)
		   {
		      int index = i + w1 * j;
		      finput[index][0] = L1[scale][index];
		      finput[index][1] = 0.0;
		   }
	       }

	     p = fftw_plan_dft_2d(h1, w1, finput, fim, -1, FFTW_ESTIMATE); fftw_execute(p);

	     fLPtmp = malloc(sizeof(double) * w1 * h1);
	     fHPtmp = malloc(sizeof(double) * w1 * h1);

	     fB = malloc(sizeof(double *) * K);
	     for (k = 0; k < K; k ++)
	       {
		 fB[k] = malloc(sizeof(double) * w1 * h1);
	       }

	     kfact  = exp(lgamma(K));
	     kfact2 = exp(lgamma(2.0 * K));
	     normfactor = normC * kfact / sqrt(K * (kfact2));

	     genHPfilter(fHPtmp, w1, h1, a/4, a/2);
	     genLPfilter(fLPtmp, w1, h1, a/4, a/2);

	     for (j = 0; j < h1; j++)
	       {
		 for (i = 0; i < w1; i++)
		   {
		      int index = i + w1 * j;
		      int k, l, t;
		      double theta;
		      if (i <  w1/2) {k =     i;}
		      if (i >= w1/2) {k = i - w1;}
		      if (j <  h1/2) {l =     j;}
		      if (j >= h1/2) {l = j - h1;}

		      for (t = 0; t < K; t ++)
			{
			  theta = atan2(l, k) - t * M_PI / K;
			  if ((fabs(theta) < M_PI / 2.0) || (fabs(theta) > 3.0 * M_PI / 2.0))
			    {
			      fB[t][index] = normfactor * pow(2.0 * cos(theta), K-1);
			    }
			  else
			    {
			      fB[t][index] = 0.0;
			    }
			}
		  }
	       }
	     for (t = 0; t < K; t ++)
	       {
		 fourier2spatialband2(w1, h1, fHPtmp, fB[t],  Br[t][scale], Bc[t][scale], conv, fim, ftmp);
	       }
	     fourier2spatialband1(w1, h1, fLPtmp, L1[scale], conv, fim, ftmp);

	     //subsampling for the next scale
	     //==============================

	     w2 = w / pow(2, scale+1);
	     h2 = h / pow(2, scale+1);

	     for (j = 0; j < h2; j ++)
	       {
		 for (i = 0; i < w2; i ++)
		   {
		      L1[scale + 1][i + w2 * j] = L1[scale][2 * (i + 2 * w2 * j)];
		   }
	       }

	     fftw_free(fim);	fftw_free(conv);	fftw_free(ftmp);	fftw_free(finput);
	     free(fLPtmp);	free(fHPtmp);		free(fB);
	  }
	return 1;
 }

void Steerable2::expand (double *inputar, double *outputar, int w, int h, int rr)
  {
    int i, j, ii, jj, w1a, h1a;

    w1a = w + 2 * rr;
    h1a = h + 2 * rr;
    for (jj = 0; jj < h1a; jj++)
      {
        for (ii = 0; ii < w1a; ii++)
          {
            i = ii - rr;
            if (i < 0)
              {
                i = rr - (ii + 1);
              }
            if (i > w - 1)
              {
                i = rr + 2 * w - 1 - ii;
              }
            j = jj - rr;
            if (j < 0)
              {
                j = rr - (jj + 1);
              }
            if (j > h - 1)
              {
                j = rr + 2 * h - 1 - jj;
              }
            outputar[ii + w1a * jj] = inputar[i + w * j];
          }
      }
  }









}

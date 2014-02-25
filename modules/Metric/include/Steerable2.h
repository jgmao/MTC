
#pragma once
#define  _CRT_SECURE_NO_WARNINGS
#ifndef STEERABLE2_H
#define STEERABLE2_H

#include "TensorLite.h"
#include "MyLib.h"
#include <vector>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>

#include <cmath>
#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
using namespace tensor;
namespace metric{
class Steerable2
{
public:
	typedef double value_type;
	typedef Tensor<value_type,1> real_data_type;
	typedef Tensor<value_type,2> data_type;
	typedef const data_type&  c_data_ref;
	typedef data_type& data_ref;
	typedef const real_data_type& c_real_ref;
public:
	EXPORTLIB Steerable2(void);
	EXPORTLIB ~Steerable2(void);
	Steerable2(c_real_ref im); //input must be a 1-channel image
protected:
	int max_ht;
	int order;
	int twidth;
	int nbands;
	Size dims;
	cv::Point ctr;
	data_type yramp,xramp,angle,log_rad;
	real_data_type im;
	data_type xrcos,yrcos,yicos;
	data_type lo0mask;
	data_type imdft,lo0dft,hi0dft;
	data_type hi0;
	vector<data_type> pyr_freq;
	vector<data_type> pyr_space;
public:
	cv::Mat raisedCosine(int maskSize, double passPosition = 0.5, double  stopPosition = 1, int octave = 1, int channels = 1);
	cv::Mat angleFilter(int maskSize, int nBands, int direction);
	EXPORTLIB vector<data_type>& buildSCFpyr(c_data_ref im, int nLevel, int nDir, int twidth, bool subsample = true);
	void buildSCFpyrLevs(data_ref loDft, vector<data_type>& rst, int nLevel, int nDir, int twidth, bool subsample = true);
	cv::Mat getChannel(cv::Mat& im, int n);
	cv::Mat toComplex(cv::Mat& A, cv::Mat& B);
	cv::Mat toComplex(cv::Mat& A);
	//STBlock toComplex(STBlock& A);
	EXPORTLIB cv::Mat getMagnitude(cv::Mat& A);
	//STBlock getMagnitude(STBlock & A);
	cv::Mat DFT(cv::Mat &A);
	//STBlock DFT(STBlock& A);
	cv::Mat IDFT(cv::Mat &A);
	//STBlock IDFT(STBlock & A);
	cv::Mat DFTShift(cv::Mat &A);	//only works when size is even
	//STBlock DFTShift(STBlock & A);	//only works when size is even
	double steera(double theta, double *B1, double *B2, double *B3, double *B4, int w, int h, int x, int y);
	double steerb(double theta, double *hB1, double *hB2, double *hB3, double *hB4, double *hB5, int w, int h, int x, int y);

	int fourier2spatialband1(int w, int h, double *otf1, double *BP, fftw_complex *conv, fftw_complex *fim, fftw_complex *ftmp);
	int fourier2spatialband2(int w, int h, double *otf1, double *otf2, double *BPr, double *BPc, fftw_complex *conv, fftw_complex *fim, fftw_complex *ftmp);
	int fourier2spatialband2a(int w, int h, double *otf1, double *otf2, double *BP, fftw_complex *conv, fftw_complex *fim, fftw_complex *ftmp);
	int fourier2spatialband3(int w, int h, double *otf1, double *otf2, double *otf3, double *BPr, double *BPc, fftw_complex *conv, fftw_complex *fim, fftw_complex *ftmp);
	int  reconststep(fftw_complex *fourtmp, fftw_complex *fourBP, double *BPr, double *BPc, int w, int h);
	int reconststepa(fftw_complex *fourtmp, fftw_complex *fourBP, double *BP, int w, int h);
	void expand (double *inputar, double *outputar, int w, int h, int rr);
	//int kanaal(double *out, int w, int h, ImlibImage *im) ;
	int genHPfilter(double *HP,  int w, int h, double x1, double x2);
	int genLPfilter(double *LP0, int w, int h, double x1, double x2);
	int decompose(int maxscale, int K, Mat im, Mat* L0, Mat* L1, Mat* B, Mat* A);
	//int decompose(int maxscale, int K, double *kanaal, int w, int h, double **L0, double **L1, double *H0,
//			double ***Br, double ***Bc, double  **Ar, double  **Ac, ImlibData  *id, ImlibImage *im);

//	int reconstruct(int maxscale, int K, double *kanaal, int w, int h, double **L0, double **L1, double *H0,
		//	double ***Br, double ***Bc, double  **Ar, double  **Ac,
			//ImlibData  *id, ImlibImage *im);

	//int variance(double *arin, double *var, int w, int h);


	//STBlock DisplayLevel(int n, vector<STBlock> & pyr);
	////STBlock DisplayAll(bool freqDomain = 1);
	//STBlock DisplayLowpass(vector<STBlock> & pyr);
	//STBlock DisplayHighpass(vector<STBlock> & pyr);
	cv::Mat normalize(cv::Mat& A);
	//STBlock normalize(STBlock& A);
	cv::Mat convertTo(cv::Mat &A, int type = CV_8U);
	//STBlock convertTo(STBlock &A, int type = CV_8U);
	EXPORTLIB vector<data_type>& getSpaceDomainPyr(void); //get space domain
	EXPORTLIB vector<data_type>& getPyr(void);
};
}
#endif


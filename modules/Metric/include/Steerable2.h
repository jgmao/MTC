
#pragma once
#define  _CRT_SECURE_NO_WARNINGS
#ifndef STEERABLE2_H
#define STEERABLE2_H
#define STEER_DBG 1
#include "TensorLite.h"
#include "MyLib.h"
#include <vector>
#include <fftw++.h>
#include <math.h>
#include <stdlib.h>
#include <Array.h>
#include <cmath>
#include <FilterBank.h>
#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
using namespace tensor;
namespace metric{
class Steerable2:public FilterBank
{
public:

	typedef vector<Array::array2<double>> Vd;
	typedef vector<Array::array2<FComplex>> Vc;
	typedef Array::array2<double> Ad;
	typedef Array::array2<FComplex> Ac;

public:
	EXPORTLIB Steerable2(void);
	EXPORTLIB ~Steerable2(void);
	Steerable2(c_real_ref im); //input must be a 1-channel image
protected:

	int maxscale;
	int K;//orientation
	int rr;
	vector<Ac*> finput;
	vector<Ac*> L0,L1;
	vector<Ac*> A;
	vector<Ac*> B;
	//vector<Ac*> B0,B1,B2,LH;
	vector<Ac> test;
	//Ac* Lp, *Hp;
	vector<Ad*> fLP;//set of LP filters
	vector<Ad*> fHP;//set of HP filters
	//Ad *fLP0,*fHP0,*fHPtmp,*fLPtmp;
	vector<Ad*> fB;
	int h,w;
	//for each scale there is a set of temp Ac s
	vector<Ac*> fim;
	vector<Ac*> conv;
	vector<Ac*> ftmp;
	Mat im;
	size_t align;
	vector<data_type> pyr_freq;
	vector<data_type> pyr_space;
	int startLevel = 0;
	//fftwpp::fft2d* Forward;
	// the start level 0 corresponds to the filter of largest posiible block size
	// to deal with variate block size
	//for examle, 64x64 is level 0, 32x32 is level 1, 16x16 is level 2...
	// only one set of filters are necessary
protected:
	int deleteFilters(void);
	void saveBand(Ad& band, string name);
	int down2Freq(Array::array2<FComplex>& L0, Array::array2<FComplex>&L1, int h0, int w0);
	int fourier2spatialband2(int w, int h, Array::array2<double>& otf1,
						 Array::array2<double>& otf2, Array::array2<FComplex>& BP,
						 Array::array2<FComplex>& conv, Array::array2<FComplex>& fim,
						 Array::array2<FComplex>& ftmp);
	int fourier2spatialband3(int K, int w, int h,
				    Array::array2<double>& otf1,
				    Array::array2<double>& otf2,
				    Array::array2<double>& otf3,
				    Array::array2<FComplex>& BP,
				    Array::array2<FComplex>& conv,
				    Array::array2<FComplex>& fim,
				    Array::array2<FComplex>& ftmp);

	int fourier2spatialband1(int w, int h, Array::array2<double>& otf1,
					       Array::array2<FComplex>& BP,
					       Array::array2<FComplex>& conv, Array::array2<FComplex>& fim,
					       Array::array2<FComplex>& ftmp);

	int genHPfilter(Array::array2<double>& HP,  int w, int h, double x1, double x2);
	int genLPfilter(Array::array2<double>& LP0, int w, int h, double x1, double x2);
public:
	int setSize(int h, int w);
	int buildFilters(int maxscale, int K);
	int expand(int rh,int rw, int border, value_type value=0);
	int decompose(void);
	int updateData(c_data_ref data);
	int getDir();
	int getLevel();
	EXPORTLIB vector<data_type>& getSpaceDomainPyr(void); //get space domain
	//EXPORTLIB vector<data_type>& getPyr(void);
};
}
#endif


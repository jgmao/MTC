////////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                     License Agreement for TensorLib
//
// Copyright (C) 2010-2011, Guoxin Jin, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use only in binary forms with or without modification,
// are permitted provided that the following conditions are met:
//
//	 * The user is a member of research group under the direction of Prof. 
//	   Thrasos Pappas in Northwestern University, IL or Prof. David Neuhoff in 
//	   University of Michigan.

//	 * The redistribution is not used for commercial applications.

//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the copyright holders or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
// By technical support and documentaion is provided by the copyright holders as will 
// in the Wiki site: https://sites.google.com/site/tensorlib/
//
//
////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#define  _CRT_SECURE_NO_WARNINGS
#ifndef STEERABLE_H
#define STEERABLE_H

//#include "Size3.h"
#include "TensorLite.h"
#include "MyLib.h"
#include <vector>
#include <cmath>
//2D STSIM2
//for one input image
//#define DEBUG
#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
using namespace tensor;
namespace metric{
class Steerable
{
public:
	typedef double value_type;
	typedef Tensor<value_type,1> real_data_type;
	typedef Tensor<value_type,2> data_type;
	typedef const data_type&  c_data_ref;
	typedef data_type& data_ref;
	typedef const real_data_type& c_real_ref;
public:
	EXPORTLIB Steerable(void);
	EXPORTLIB ~Steerable(void);
	Steerable(c_real_ref im); //input must be a 1-channel image
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

	//STBlock DisplayLevel(int n, vector<STBlock> & pyr);
	////STBlock DisplayAll(bool freqDomain = 1);
	//STBlock DisplayLowpass(vector<STBlock> & pyr);
	//STBlock DisplayHighpass(vector<STBlock> & pyr);
	cv::Mat normalize(cv::Mat& A);
	//STBlock normalize(STBlock& A);
	cv::Mat convertTo(cv::Mat &A, int type = CV_8U);
	//STBlock convertTo(STBlock &A, int type = CV_8U);
	EXPORTLIB vector<data_type>& getSpaceDomainPyr(void); //get space domain
};
} 
#endif

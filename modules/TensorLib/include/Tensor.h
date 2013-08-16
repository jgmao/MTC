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

//critial version notes:
//This is second version of Tensor, which inherited directly from cv::Mat
//because Mat will support N-dim element. 
//From this version, the storage structure of 3D block is changed:
//1. operator= and operator() will return a reference not copy! (consistent to OpenCV)
//2. in 3D block, 1st dim use to represent the time-axis.
//3. Clone() method is changed, more good for deep copy

/// 20130813 change to GCC version
/// use GCC4.8 and c++11
/// the purpose of Tensor<T,cn>is majorly force on support both real and complex number data
/// sepearte the simple operation to complex operation (such as STSIM)
/// the complex operation will use GPU
/// and also try to seperate the CommonDef def's for tensor (data structure) and coder(MTC)
/// keep the operation done using Mat 
#pragma once
#ifndef TENSOR_H
#define TENSOR_H
#include "CommonDef.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
//#include <cv.h>
//#include "highgui.h"
#include "Size3.h"
#include <iostream>
#include <fstream>
#include "MyLib.h"
#include "Cube.h"
#include <iomanip>
#include <limits>
#include <cstdarg>
#include <mutex>
#include <stdlib.h>
#include <thread>
//20130815#include "ThinPlateSpline.h"
//#include "Metric.h"
//#include "MetricBase.h"
//#include <boost\format.hpp>
using namespace cv;
using namespace std;

#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif

namespace tensor{
template<class T, size_t cn> class Tensor: public Mat
{
  //Metric mc;
public:
  //20130815 static Mat stsim2_lse_weight;
public:
  //static Metric* mc;
	typedef Vec<T,cn> value_type;
	typedef const Vec<T,cn> c_value_type;
	typedef value_type& ref_type;
	typedef c_value_type& c_ref_type;
public:
/////////////2.1 constructor, copy  assignment and type conversion ////////////////////////

	EXPORTLIB Tensor<T,cn>(void);
	EXPORTLIB ~Tensor<T,cn>(void);
	EXPORTLIB Tensor<T,cn>(int height, int width, int depth, c_ref_type val = value_type::all(0));	
	EXPORTLIB Tensor<T,cn>(const Size3& size, c_ref_type val) ;
        EXPORTLIB Tensor<T,cn>(const Size3& size);
	EXPORTLIB Tensor<T,cn>(const Tensor<T,cn>& ts);
	EXPORTLIB Tensor<T,cn>(const Mat& mt);
	EXPORTLIB Tensor<T,cn>(const string cFileName);
	EXPORTLIB Tensor<T,cn>& operator= (const Tensor<T,cn>& ts);
	EXPORTLIB Tensor<T,2> ToComplex(void) const;
////////////2.2 public utilities	/////////////////////////////////////////////////////
	EXPORTLIB void AssertRange(const Point3i& pos, const Size3& sz) const;
	template<size_t cn2> Tensor<T,cn2> CvtColor(int flag = COLOR_RGB2GRAY);//convert color space
	EXPORTLIB bool IsInside(const Point3i& pos) const;
	EXPORTLIB void Display(int flag =1) const;
	EXPORTLIB void Display(int sec, int flag) const;
	//move this to utility module inside Tensor 
  	//EXPORTLIB void DisplayAll(vector<Tensor<T,cn>>& vts, int row, int col,bool save=false, string savename="./rst.tiff") const;
	EXPORTLIB void Load(string cFileName);
	EXPORTLIB void SetFrame(int i, const Mat& frm);
	EXPORTLIB void SetFileName(const string& name);
	EXPORTLIB void SetOffset(const Point3i& pos);
	EXPORTLIB Tensor<T,cn>& SetBlock(const Point3i& pos, const Tensor<T,cn>& ts);
	EXPORTLIB Tensor<T,cn>& SetBlock(const Tensor<T,cn>& ts);
	EXPORTLIB void SaveBlock(const string& cFilename, bool isGray = false);
	EXPORTLIB void Print(const string& fname = "print_log", bool tofileonly=false) const;
	/* 20130815 move to Metric Module 
	EXPORTLIB void SetSubWinSize(const Size3& sz);
	EXPORTLIB void SetSubWinStep(const Size3& sz);
	*/
///////////2.3 data accessors //////////////////////////////////////////////////////////
	EXPORTLIB Size3 size(void) const;
	EXPORTLIB ref_type operator() (int x, int y, int z);
	EXPORTLIB c_ref_type operator() (int x, int y, int z) const;
	EXPORTLIB Tensor<T,cn> operator()(const Cube& roi);
	EXPORTLIB Mat GetFrame(int i);
  	EXPORTLIB Mat& GetFrameRef(int i) ;
  	EXPORTLIB void GetFrameRef(int i, Mat& rst) const;
	EXPORTLIB const Mat GetFrame(int i) const;
	EXPORTLIB const Mat operator[](int i) const;
	EXPORTLIB ref_type operator[](const Point3i& pos);  
	EXPORTLIB c_ref_type operator[](const Point3i& pos) const;
	EXPORTLIB Mat operator[](int i);
	EXPORTLIB Tensor<T,cn> GetBlock(const Cube& roi);
  	//EXPORTLIB Tensor<T,cn>& GetBlockRef(const Cube& roi);
	EXPORTLIB Tensor<T,cn> GetBlock(const Cube& roi) const;
	EXPORTLIB void Ref(const Cube& roi, Tensor<T,cn>& dst) const;
	EXPORTLIB Point3i offset(void) const;
	EXPORTLIB Tensor<T,cn> Clone() const;
	EXPORTLIB Tensor<T,cn> Crop(const Point3i& pos, const Size3& sz) const; //crop is same as GetBlock but return a copy of data
	EXPORTLIB Tensor<T,cn> Row(int i); //this is untested
	//20130815  move to Metric module
	/*
	EXPORTLIB Size3 GetSubWinSize(void) const;
	EXPORTLIB Size3 GetSubWinStep(void) const;
	*/
	//20130815 move to lighting module
	/*
	EXPORTLIB vector<Vec<T,cn>> GetTagLighting(void) const;
	EXPORTLIB vector<Vec<T,cn>> GetCanLighting(void) const;
	*/
////////2.4 arithmatic operations //////////////////////////////////////////
	EXPORTLIB Tensor<T,cn> operator*(const Tensor<T,cn>& ts ) const;
	EXPORTLIB Tensor<T,cn> operator/(const Tensor<T,cn>& ts ) const;
	EXPORTLIB Tensor<T,cn> operator+(const Tensor<T,cn>& ts ) const;
	EXPORTLIB Tensor<T,cn> operator-(const Tensor<T,cn>& ts ) const;

	EXPORTLIB Tensor<T,cn> operator+(c_ref_type s ) const;
	EXPORTLIB Tensor<T,cn> operator-(c_ref_type s ) const;	
	EXPORTLIB Tensor<T,cn> operator*(c_ref_type s ) const;
	EXPORTLIB Tensor<T,cn> operator/(c_ref_type s ) const;

	EXPORTLIB Tensor<T,cn> operator+(T s ) const;//@//
	EXPORTLIB Tensor<T,cn> operator-(T s ) const;//@//	
	EXPORTLIB Tensor<T,cn> operator*(T s ) const;//@//
	EXPORTLIB Tensor<T,cn> operator/(T s ) const;//@//

	EXPORTLIB value_type All(T val) const;//@
	EXPORTLIB Tensor<T,cn>Abs(void) const; // if data type is not float/double this will throw type exception
	EXPORTLIB Tensor<T,cn>AbsDiff(c_ref_type s) const;//@
	EXPORTLIB Tensor<T,cn>AbsDiff(const Tensor<T,cn>& ts) const;//@

	EXPORTLIB Tensor<T,cn>Conjugate(void) const;//@
	EXPORTLIB Tensor<T,cn>ExtendBoundary(Size3 extSz = Size3(1,1,0), value_type val = value_type::all(0)) const;
	EXPORTLIB Tensor<T,cn>ExtendHalfBoundary(Size3 extSz = Size3(1,1,0), value_type val = value_type::all(0), bool which_side = 0) const; // which side 0: left and upper 1: right and below
	//20130815 use GPU algorithm here to compute 
	EXPORTLIB Vec<T,cn> Mean(void) const;
	EXPORTLIB Vec<T,cn> Var(void) const;
	EXPORTLIB Vec<T,cn> LRMean(void) const; //mean of lower right triangle
	EXPORTLIB Vec<T,cn> URMean(void) const; //mean of lower right triangle
	EXPORTLIB Vec<T,cn> PercentageMean(double low=0, double high=1) const; //mean taken from (low - high) percentage
	EXPORTLIB Tensor<T,cn>LocalMean(const Mat& ker, const Size3& subWinStep = Size3(1,1,1)) const;//only support depth =1 kernel
	EXPORTLIB Tensor<T,cn>LocalMean(const Size3& subWinSize,const Size3& subWinStep = Size3(1,1,1)) const;
	EXPORTLIB Tensor<T,cn>LocalVariance( const Tensor<T,cn>& mu, const Mat& ker, const Size3& subWinStep = Size3(1,1,1)) const;
	EXPORTLIB Tensor<T,cn>LocalVariance( const Tensor<T,cn>& mu, const Size3& subWinSize, const Size3& subWinStep = Size3(1,1,1)) const;
	EXPORTLIB Tensor<T,cn>Pow(double p) const;//@
  	EXPORTLIB Tensor<T,cn>Log(void) const;//take the ln
 	EXPORTLIB Tensor<T,cn>Exp(void) const;//take the e^()
	EXPORTLIB Tensor<T,1> Real(void) const; //a copy of real plane
	EXPORTLIB Tensor<T,1> Imag(void) const; //a copy of imag plane
	EXPORTLIB Tensor<T,cn> Sqrt(void) const; //will throw exception for T = uchar
	EXPORTLIB Vec<T,cn> Sum(void) const;
	EXPORTLIB Tensor<T,cn> Transpose(void) const;//@
	EXPORTLIB Tensor<T,cn> SubSample(const Size3& subrate);
	EXPORTLIB Vec<T,cn> Min() const;
	EXPORTLIB Tensor<T,cn> MaxClip(const Vec<T,cn>& value) const;
  	EXPORTLIB Tensor<T,cn> MinClip(const Vec<T,cn>& value) const;
	EXPORTLIB Tensor<T,cn> Flip(bool byX, bool byY, bool byZ=false) const;
////////2.5 comparison and metric //////////////////////////////////////////////
///// 20130815 change this out side the Tensor class, use another module
/*
	EXPORTLIB bool ComputeAIM(const Tensor<T,cn> ts, double thrd) const; // adaptive interpolation metric	for PQI
	EXPORTLIB Mat CompareElement(const Mat& inMat, double thrd, int flag = CMP_LE) const;//@
	EXPORTLIB Tensor<uchar,1> CompareElement(double thrd, int flag=CMP_LT) const;//@// the return type is definited
	EXPORTLIB double ComputeMSE(const Tensor<T,cn> ts) const; // if lambda = 0, then it is normal MSE, otherwise it is MSE with constraint
	EXPORTLIB double ComputeSAD(const Tensor<T,cn> ts) const;
	EXPORTLIB double ComputePSNR(const Tensor<T,cn> ts) const;
  	EXPORTLIB double ComputeLRI(const Tensor<T,cn> ts) const;
	EXPORTLIB double ComputeSSIM(const Tensor<T,cn> ts, const Size3& subWinSize = Size3(8,8,1), const Size3& subWinStep = Size3(1,1,1),int nLevel=3, int nDir=4, int boundary_cut = FILTER_BOUND_HALF, int ststi2_pool_type = STSIM2_POOL_MIN, int stsim2_modifer = STSIM2_BASELINE) const;
  	EXPORTLIB vector<Tensor<double,2> > ComputeStatistics(const Size3& subWinSize = Size3(8,8,1), const Size3& subWinStep = Size3(1,1,1), bool subsample = false, int nLevel = 3, int nDir=4, bool changeWin=false) const;
  	EXPORTLIB double ComputeSVMMetric(const Tensor<T,cn> ts, const Size3& subWinSize, const Size3& subWinStep) const;
	//EXPORTLIB double Compare(const Tensor<T,cn> ts, int criteria, double param1=3, double param2=4, int param3 = 1, int param4 = STSIM2_POOL_MIN, int param5 = STSIM2_BASELINE, cv::Mat& param6 = cv::Mat()) const;
  	EXPORTLIB double Compare(const Tensor<T,cn> ts, int criteria, Printable_t param1=Printable_t(0), Printable_t param2=Printable_t(0), Printable_t param3=Printable_t(0), Printable_t param4=Printable_t(0), Printable_t param5=Printable_t(0), Printable_t param6=Printable_t(0)) const;
  	EXPORTLIB cv::Mat EstimateVarForMahalanobis(Size3 wsize, Size3 stepsize);
  	EXPORTLIB double ComputeMahalanobis(const Tensor<T,cn> ts, Size3 subWinSize, Size3 subWinStep, const Mat& iMcovar) const;
*/
/////////2.6 DSP ////////////////////////////////////////////////////////////
	//20130815 use GPU here
	EXPORTLIB Tensor<T,cn> FreqComplexFilter(cv::Mat& kernel, bool conj=false) const;
	EXPORTLIB Tensor<T,cn> DFT(void) const;
	EXPORTLIB Tensor<T,cn> IDFT(void) const; //@//only 2D IDFT on each frame
	EXPORTLIB Mat DFTShift(const Mat& A) const; //@// a none template member
	EXPORTLIB Tensor<T,cn> DFTShift(void) const; //@// only 2D Shift on each frame
	EXPORTLIB Tensor<T,cn> Laplacian(void) const;
	EXPORTLIB Tensor<T,cn> Filter2D(const Mat& ker, int boundary = (int)FilterBoundary::FILTER_BOUND_VALID) const;
	EXPORTLIB Tensor<T,cn> Normalize(void) const;
	//20130815 move lighting correction block out side Tensor
	/*
	EXPORTLIB Tensor<T,cn> LSFitting(int order = 2) const;
	EXPORTLIB Tensor<T,cn> BuildLightingPlane(const Tensor<T,cn>& param, int order = 2) const;
	EXPORTLIB Tensor<T,cn> LightingCorrection(const Tensor<T,cn>& candid,bool saveCodeLength = true);
	EXPORTLIB Tensor<T,cn> LightingCorrection(const Tensor<T,cn>& candid, const Tensor<T,cn>& VQCodebook);
	EXPORTLIB int GetLightingCodeLength() const;
	EXPORTLIB void SetLightingCodeLength(int l);
	EXPORTLIB Tensor<T,cn> ComputeSAT(void) const;//compute sumed-area-table
	EXPORTLIB Tensor<T,cn>& ComputeSAT(Tensor<T,cn>& SAT, Point3i& sPos, Point3i& ePos) const;
  	EXPORTLIB Tensor<T,cn> ComputeTPSS(double p) const;//thin plate spline smoothing
	EXPORTLIB void RecodeLighting(void);
	*/
/////////////////// 2.7 geometric distortions 
// 20130815 move distortion outside Tensor
/*
	EXPORTLIB Tensor<T,cn> MicroShift(Size3 win, double maxdev);
  	EXPORTLIB Tensor<T,cn> MicroRotate(Size3 win, double maxdev);
  	EXPORTLIB value_type Bilinear(int i,int j, int offset_i, int offset_j, int ww, int hei1, int wid1, int hei2, int wid2, double phi);
*/

/// 2.8 other //////////////
public:
	bool debugtrigger;
private:
///// 20130815 move to Metirc module
/*
	Tensor<T,cn> SearchCodeword(const Tensor<T,cn>& val, const Tensor<T,cn>& VQCodeBook);
	double ComputeCrossTerm(const cv::Mat &im11,const cv::Mat &im12,const cv::Mat &im21,const cv::Mat &im22) const;
	Tensor<double,1> ComputeRho(const Tensor<double,2>& im11, const Tensor<double,2>& im12, const Size3& subWinSize, const Size3& subWinStep) const;
	Tensor<double,1> ComputeCrossTerm(const Tensor<double,2>& im11, const Tensor<double,2>& im12,Tensor<double,2>& im21, const Tensor<double,2>& im22, const Size3& subWinSize, const Size3& subWinStep) const;
*/
///// 2.9 Protected Data
protected:
	Size3 tsSize;
	Point3i tsOffset;
	string cFileName;
	Mat mxFrame; // a temp (header) to store the frame info
	const Mat cmxFrame;
	//20130815 move to Metric module and Lighthing Module
	/*
	Size3 subWinSize; //window size of STSIM, default (8,8,1)
	Size3 subWinStep; //step of sliding window, default (1,1,1)
	vector<Vec<T,cn>> lightTag; //original
	vector<Vec<T,cn>> lightCan; //candidate
	vector<double> lightingDCTCoeffStat;
	int codeLength;
	*/	
  	mutable std::mutex tsmutex;
  
};

////////////////////////
template<class T, size_t cn> template<size_t cn2> Tensor<T,cn2> Tensor<T,cn>::CvtColor(int flag)
{
	
		Tensor<float,cn2> rst(size());
		Tensor<float,cn> floatIn = Tensor<float,cn>(*this);
		for (int i=0; i< size().depth; i++)
		{
			cv::cvtColor(floatIn[i],rst[i],flag);
		}
		return Tensor<T,cn2>(rst);
}
}
#endif

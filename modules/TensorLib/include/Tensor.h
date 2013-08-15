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

#pragma once
#ifndef TENSOR_H
#define TENSOR_H
#include "CommonDef.h"
#include <opencv2/opencv.hpp>
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
#include "ThinPlateSpline.h"
//#include "Metric.h"
//#include "MetricBase.h"
//#include <boost\format.hpp>
using namespace cv;
using namespace std;
template<class T, size_t cn> class Tensor : public Mat
{
  //Metric mc;
public:
  static Mat stsim2_lse_weight;
public:
  //static Metric* mc;
	typedef Vec<T,cn> value_type;
	typedef const Vec<T,cn> c_value_type;
	typedef value_type& ref_type;
	typedef c_value_type& c_ref_type;
public:
/////////////2.1 constructor, copy  assignment and type conversion ////////////////////////

	__declspec(dllexport) Tensor(void);
	__declspec(dllexport) ~Tensor(void);
	__declspec(dllexport) Tensor(int height, int width, int depth, c_ref_type val = value_type::all(0));	
	__declspec(dllexport) Tensor(const Size3& size, c_ref_type val) ;
  __declspec(dllexport) Tensor(const Size3& size);
	__declspec(dllexport) Tensor(const Tensor& ts);
	__declspec(dllexport) Tensor(const Mat& mt);
	__declspec(dllexport) Tensor(const string cFileName);
	__declspec(dllexport) Tensor& operator= (const Tensor& ts);
	__declspec(dllexport) Tensor<T,2> ToComplex(void) const;
////////////2.2 public utilities	/////////////////////////////////////////////////////
	__declspec(dllexport) void AssertRange(const Point3i& pos, const Size3& sz) const;
	template<size_t cn2> Tensor<T,cn2> CvtColor(int flag = CV_RGB2GRAY);//convert color space
	__declspec(dllexport) bool IsInside(const Point3i& pos) const;
	__declspec(dllexport) void Display(int flag =1) const;
	__declspec(dllexport) void Display(int sec, int flag) const;
  __declspec(dllexport) void DisplayAll(vector<Tensor<T,cn>>& vts, int row, int col,bool save=false, string savename="./rst.tiff") const;
	__declspec(dllexport) void Load(string cFileName);
	__declspec(dllexport) void SetFrame(int i, const Mat& frm);
	__declspec(dllexport) void SetFileName(const string& name);
	__declspec(dllexport) void SetOffset(const Point3i& pos);
	__declspec(dllexport) Tensor& SetBlock(const Point3i& pos, const Tensor<T,cn>& ts);
	__declspec(dllexport) Tensor& SetBlock(const Tensor<T,cn>& ts);
	__declspec(dllexport) void SaveBlock(const string& cFilename, bool isGray = false);
	__declspec(dllexport) void Print(const string& fname = "print_log", bool tofileonly=false) const;
	__declspec(dllexport) void SetSubWinSize(const Size3& sz);
	__declspec(dllexport) void SetSubWinStep(const Size3& sz);
///////////2.3 data accessors //////////////////////////////////////////////////////////
	__declspec(dllexport) Size3 size(void) const;
	__declspec(dllexport) ref_type operator() (int x, int y, int z);
	__declspec(dllexport) c_ref_type operator() (int x, int y, int z) const;
	__declspec(dllexport) Tensor operator()(const Cube& roi);
	__declspec(dllexport) Mat GetFrame(int i);
  __declspec(dllexport)  Mat& GetFrameRef(int i) ;
  __declspec(dllexport) void GetFrameRef(int i, Mat& rst) const;
	__declspec(dllexport) const Mat GetFrame(int i) const;
	__declspec(dllexport) const Mat operator[](int i) const;
	__declspec(dllexport) ref_type operator[](const Point3i& pos);  
	__declspec(dllexport) c_ref_type operator[](const Point3i& pos) const;
	__declspec(dllexport) Mat operator[](int i);
	__declspec(dllexport) Tensor GetBlock(const Cube& roi);
  //__declspec(dllexport) Tensor<T,cn>& GetBlockRef(const Cube& roi);
	__declspec(dllexport) Tensor GetBlock(const Cube& roi) const;
	__declspec(dllexport) void Ref(const Cube& roi, Tensor& dst) const;
	__declspec(dllexport) Point3i offset(void) const;
	__declspec(dllexport) Tensor Clone() const;
	__declspec(dllexport) Tensor Crop(const Point3i& pos, const Size3& sz) const; //crop is same as GetBlock but return a copy of data
	__declspec(dllexport) Tensor Row(int i); //this is untested
	__declspec(dllexport) Size3 GetSubWinSize(void) const;
	__declspec(dllexport) Size3 GetSubWinStep(void) const;
	__declspec(dllexport) vector<Vec<T,cn>> GetTagLighting(void) const;
	__declspec(dllexport) vector<Vec<T,cn>> GetCanLighting(void) const;
////////2.4 arithmatic operations //////////////////////////////////////////
	__declspec(dllexport) Tensor operator*(const Tensor& ts ) const;
	__declspec(dllexport) Tensor operator/(const Tensor& ts ) const;
	__declspec(dllexport) Tensor operator+(const Tensor& ts ) const;
	__declspec(dllexport) Tensor operator-(const Tensor& ts ) const;

	__declspec(dllexport) Tensor operator+(c_ref_type s ) const;
	__declspec(dllexport) Tensor operator-(c_ref_type s ) const;	
	__declspec(dllexport) Tensor operator*(c_ref_type s ) const;
	__declspec(dllexport) Tensor operator/(c_ref_type s ) const;

	__declspec(dllexport) Tensor operator+(T s ) const;//@//
	__declspec(dllexport) Tensor operator-(T s ) const;//@//	
	__declspec(dllexport) Tensor operator*(T s ) const;//@//
	__declspec(dllexport) Tensor operator/(T s ) const;//@//

	__declspec(dllexport) value_type All(T val) const;//@
	__declspec(dllexport) Tensor Abs(void) const; // if data type is not float/double this will throw type exception
	__declspec(dllexport) Tensor AbsDiff(c_ref_type s) const;//@
	__declspec(dllexport) Tensor AbsDiff(const Tensor& ts) const;//@

	__declspec(dllexport) Tensor Conjugate(void) const;//@
	__declspec(dllexport) Tensor ExtendBoundary(Size3 extSz = Size3(1,1,0), value_type val = value_type::all(0)) const;
	__declspec(dllexport) Tensor ExtendHalfBoundary(Size3 extSz = Size3(1,1,0), value_type val = value_type::all(0), bool which_side = 0) const; // which side 0: left and upper 1: right and below
	__declspec(dllexport) Vec<T,cn> Mean(void) const;
	__declspec(dllexport) Vec<T,cn> Var(void) const;
	__declspec(dllexport) Vec<T,cn> LRMean(void) const; //mean of lower right triangle
	__declspec(dllexport) Vec<T,cn> URMean(void) const; //mean of lower right triangle
	__declspec(dllexport) Vec<T,cn> PercentageMean(double low=0, double high=1) const; //mean taken from (low - high) percentage
	__declspec(dllexport) Tensor LocalMean(const Mat& ker, const Size3& subWinStep = Size3(1,1,1)) const;//only support depth =1 kernel
	__declspec(dllexport) Tensor LocalMean(const Size3& subWinSize,const Size3& subWinStep = Size3(1,1,1)) const;
	__declspec(dllexport) Tensor LocalVariance( const Tensor<T,cn>& mu, const Mat& ker, const Size3& subWinStep = Size3(1,1,1)) const;
	__declspec(dllexport) Tensor LocalVariance( const Tensor<T,cn>& mu, const Size3& subWinSize, const Size3& subWinStep = Size3(1,1,1)) const;
	__declspec(dllexport) Tensor Pow(double p) const;//@
  __declspec(dllexport) Tensor Log(void) const;//take the ln
  __declspec(dllexport) Tensor Exp(void) const;//take the e^()
	__declspec(dllexport) Tensor<T,1> Real(void) const; //a copy of real plane
	__declspec(dllexport) Tensor<T,1> Imag(void) const; //a copy of imag plane
	__declspec(dllexport) Tensor<T,cn> Sqrt(void) const; //will throw exception for T = uchar
	__declspec(dllexport) Vec<T,cn> Sum(void) const;
	__declspec(dllexport) Tensor Transpose(void) const;//@
	__declspec(dllexport) Tensor SubSample(const Size3& subrate);
	__declspec(dllexport) Vec<T,cn> Min() const;
	__declspec(dllexport) Tensor<T,cn> MaxClip(const Vec<T,cn>& value) const;
  __declspec(dllexport) Tensor<T,cn> MinClip(const Vec<T,cn>& value) const;
	__declspec(dllexport) Tensor<T,cn> Flip(bool byX, bool byY, bool byZ=false) const;
////////2.5 comparison and metric //////////////////////////////////////////////
	__declspec(dllexport) bool ComputeAIM(const Tensor& ts, double thrd) const; // adaptive interpolation metric	for PQI
	__declspec(dllexport) Mat CompareElement(const Mat& inMat, double thrd, int flag = CMP_LE) const;//@
	__declspec(dllexport) Tensor<uchar,1> CompareElement(double thrd, int flag=CMP_LT) const;//@// the return type is definited
	__declspec(dllexport) double ComputeMSE(const Tensor& ts) const; // if lambda = 0, then it is normal MSE, otherwise it is MSE with constraint
	__declspec(dllexport) double ComputeSAD(const Tensor& ts) const;
	__declspec(dllexport) double ComputePSNR(const Tensor& ts) const;
  __declspec(dllexport) double ComputeLRI(const Tensor& ts) const;
	__declspec(dllexport) double ComputeSSIM(const Tensor& ts, const Size3& subWinSize = Size3(8,8,1), const Size3& subWinStep = Size3(1,1,1),int nLevel=3, int nDir=4, int boundary_cut = FILTER_BOUND_HALF, int ststi2_pool_type = STSIM2_POOL_MIN, int stsim2_modifer = STSIM2_BASELINE) const;
  __declspec(dllexport) vector<Tensor<double,2> > ComputeStatistics(const Size3& subWinSize = Size3(8,8,1), const Size3& subWinStep = Size3(1,1,1), bool subsample = false, int nLevel = 3, int nDir=4, bool changeWin=false) const;
  __declspec(dllexport) double ComputeSVMMetric(const Tensor& ts, const Size3& subWinSize, const Size3& subWinStep) const;
	//__declspec(dllexport) double Compare(const Tensor& ts, int criteria, double param1=3, double param2=4, int param3 = 1, int param4 = STSIM2_POOL_MIN, int param5 = STSIM2_BASELINE, cv::Mat& param6 = cv::Mat()) const;
  __declspec(dllexport) double Compare(const Tensor& ts, int criteria, Printable_t param1=Printable_t(0), Printable_t param2=Printable_t(0), Printable_t param3=Printable_t(0), Printable_t param4=Printable_t(0), Printable_t param5=Printable_t(0), Printable_t param6=Printable_t(0)) const;
  __declspec(dllexport) cv::Mat EstimateVarForMahalanobis(Size3 wsize, Size3 stepsize);
  __declspec(dllexport) double ComputeMahalanobis(const Tensor& ts, Size3 subWinSize, Size3 subWinStep, const Mat& iMcovar) const;

/////////2.6 DSP ////////////////////////////////////////////////////////////
	__declspec(dllexport) Tensor<T,cn> FreqComplexFilter(cv::Mat& kernel, bool conj=false) const;
	__declspec(dllexport) Tensor<T,cn> DFT(void) const;
	__declspec(dllexport) Tensor<T,cn> IDFT(void) const; //@//only 2D IDFT on each frame
	__declspec(dllexport) Mat DFTShift(const Mat& A) const; //@// a none template member
	__declspec(dllexport) Tensor<T,cn> DFTShift(void) const; //@// only 2D Shift on each frame
	__declspec(dllexport) Tensor<T,cn> Laplacian(void) const;
	__declspec(dllexport) Tensor<T,cn> Filter2D(const Mat& ker, int boundary = FILTER_BOUND_VALID) const;
	__declspec(dllexport) Tensor<T,cn> Normalize(void) const;
	__declspec(dllexport) Tensor<T,cn> LSFitting(int order = 2) const;
	__declspec(dllexport) Tensor<T,cn> BuildLightingPlane(const Tensor<T,cn>& param, int order = 2) const;
	__declspec(dllexport) Tensor<T,cn> LightingCorrection(const Tensor<T,cn>& candid,bool saveCodeLength = true);
	__declspec(dllexport) Tensor<T,cn> LightingCorrection(const Tensor<T,cn>& candid, const Tensor<T,cn>& VQCodebook);
	__declspec(dllexport) int GetLightingCodeLength() const;
	__declspec(dllexport) void SetLightingCodeLength(int l);
	__declspec(dllexport) Tensor<T,cn> ComputeSAT(void) const;//compute sumed-area-table
	__declspec(dllexport) Tensor<T,cn>& ComputeSAT(Tensor<T,cn>& SAT, Point3i& sPos, Point3i& ePos) const;
  __declspec(dllexport) Tensor<T,cn> ComputeTPSS(double p) const;//thin plate spline smoothing
/////////////////// 2.7 geometric distortions 
  __declspec(dllexport) Tensor<T,cn> MicroShift(Size3 win, double maxdev);
  __declspec(dllexport) Tensor<T,cn> MicroRotate(Size3 win, double maxdev);
  __declspec(dllexport) value_type Bilinear(int i,int j, int offset_i, int offset_j, int ww, int hei1, int wid1, int hei2, int wid2, double phi);
// other
	__declspec(dllexport) void RecodeLighting(void);
public:
	bool debugtrigger;
  //static Metric mc;
private:
	Tensor<T,cn> SearchCodeword(const Tensor<T,cn>& val, const Tensor<T,cn>& VQCodeBook);
	double ComputeCrossTerm(const cv::Mat &im11,const cv::Mat &im12,const cv::Mat &im21,const cv::Mat &im22) const;
	Tensor<double,1> ComputeRho(const Tensor<double,2>& im11, const Tensor<double,2>& im12, const Size3& subWinSize, const Size3& subWinStep) const;
	Tensor<double,1> ComputeCrossTerm(const Tensor<double,2>& im11, const Tensor<double,2>& im12,Tensor<double,2>& im21, const Tensor<double,2>& im22, const Size3& subWinSize, const Size3& subWinStep) const;

protected:
	Size3 tsSize;
	Point3i tsOffset;
	string cFileName;
	Mat mxFrame; // a temp (header) to store the frame info
  const Mat cmxFrame;
	Size3 subWinSize; //window size of STSIM, default (8,8,1)
	Size3 subWinStep; //step of sliding window, default (1,1,1)
	vector<Vec<T,cn>> lightTag; //original
	vector<Vec<T,cn>> lightCan; //candidate
	vector<double> lightingDCTCoeffStat;
	int codeLength;
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

#endif

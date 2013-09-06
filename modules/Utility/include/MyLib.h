#pragma once
//#pragma warning (disable : 4267)
#ifndef MYLIB_H
#define MYLIB_H

#include <vector>
//#include <cv.h>
#include <opencv2/opencv.hpp>
//#include <boost/lexical_cast.hpp>
#include <list>
#include <string>
#include <map>
#include <cstdarg>
//#include "Tensor.h"

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

//#include <highgui.h>
#include <iostream>
//#include "Size3.h"
#include <fstream>
using namespace std;
using namespace cv;
#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
namespace mylib
{
////test 
/*20130813 temp disable to seperate from tensorlib, move this to tensor lib
template <typename T> Size3_<T> Ceil(const Size3_<T>& sz)
{
	return Size3_<T>(std::ceil(sz.height),std::ceil(sz.width),std::ceil(sz.depth));
}
*/
template <typename T> cv::Mat VecToMat(vector<vector<T>> mx, int elementType)
{
	//be very careful about the element type, use CV_32F in most case

	typename vector<vector<T>>::iterator iter1 = mx.begin();
	typename vector<T>::iterator iter2 = iter1->begin();
	cv::Mat rst(mx.size(),iter1->size(),CV_64F);
	int i=0;
	int j=0;
	for (iter1=mx.begin();iter1!= mx.end(); iter1++)
	{
		//iter2= iter1->begin();
		for (iter2=iter1->begin();iter2!=iter1->end(); iter2++)
		{
			rst.at<double>(i,j)=mx[i][j];
			j++;
		}
		i++;
		j=0;
	}
	return rst;
}
//
//
//template <typename T> vector<vector<vector<T>>> conv(vector<vector<vector<T>>> src, vector<vector<vector<T>>> kernel, int flag);
//template <typename T> vector<vector<T>> conv(vector<vector<T>> src, vector<vector<T>> kernel, int flag);
template <typename T> vector<T> conv(vector<T> src, vector<T> kernel, int flag)
{
	int lengthSrc = src.size();
	int lengthKel = kernel.size();
	vector<T> rst(lengthSrc+lengthKel-1);
	for (int t=0;t<lengthSrc+lengthKel-1;t++)
	{
		for (int i=0;i<=t; i++)
		{
			if (t-i < lengthKel && i < lengthSrc)
				rst[t]+=src[i]*kernel[t-i];
		}
	}
	return rst;
}
template <typename T> vector<vector<T>> crossproduct(vector<T>& left, vector<T>& right)
{
	int lengthLeft = left.size();
	int lengthRight = right.size();
	vector<vector<T>> rst(lengthLeft,vector<T>(lengthRight));
	for (int i =0; i<lengthLeft; i++)
		for (int j=0; j< lengthRight;j++)
		{
			rst[i][j]=left[i]*right[j];
		}
	return rst;
}
EXPORTLIB void DisplayMat(const cv::Mat &m, const string & name="undef", bool tofileonly=false);
EXPORTLIB cv::Mat Filter(cv::Mat src, cv::Mat kernel, int flag);
//STBlock Filter(STBlock src, STBlock kernel, int flag);
//
EXPORTLIB cv::Mat BinomialKernel(int size, int nDataType = CV_32F);

EXPORTLIB cv::Mat CrossProduct(const cv::Mat &colVec, const cv::Mat &rowVec);
EXPORTLIB double WeightedSum(const cv::Mat& A, const cv::Mat& W);
EXPORTLIB std::complex<double> WeightedSumComplex(const cv::Mat& A, const cv::Mat& W);


EXPORTLIB cv::Mat GenGaussKer(int size, double sigma, int type);

template <typename T> Vec<T, 3> Vec3Diff(Vec<T,3> &A, Vec<T, 3> &B)
{
	Vec<T, 3> C;
	C[0] = A[0]-B[0];
	C[1] = A[1]-B[1];
	C[2] = A[2]-B[2];
	return C;
}

template <typename T, size_t cn> inline Vec<T,cn> VecDiv(const Vec<T,cn>& A, const T& b)
{
	Vec<T,cn> C;
	for (int i=0; i< cn; i++) C[i] = saturate_cast<T>(A[i]/b);
	return C;
}
template <typename T, size_t cn> inline Vec<T,cn> VecMul(const Vec<T,cn>& A, const Vec<T,cn>& B) 
{
	Vec<T,cn> C;
	for (int i=0; i< cn; i++) C[i] = saturate_cast<T>(A[i]*B[i]);
	return C;
}

template <typename T, size_t cn> inline Vec<T,cn> VecDiv(const Vec<T,cn>& A, const Vec<T,cn>& B)
{
	Vec<T,cn> C;
	for (int i=0; i< cn; i++) C[i] = saturate_cast<T>(A[i]/B[i]);
	return C;
}

template <typename T, size_t cn> inline Vec<T,cn> VecPow(const Vec<T,cn>& A, double p)
{
	Vec<T,cn> C;
	for (int i=0; i< cn; i++)
		C[i] = saturate_cast<T>(pow(double(A[i]),p));
	return C;

}


template <typename T> Vec<T, 3> Vec3Div(Vec<T,3> &A, T b)
{
	Vec<T, 3> C;
	C[0] = saturate_cast<T>(A[0]/b);
	C[1] = saturate_cast<T>(A[1]/b);
	C[2] = saturate_cast<T>(A[2]/b);
	return C;
}

template <typename T> Vec<T,3> Vec3Div(Vec<T,3> &A, int b)
{
	return Vec3Div<T>(A,T(b));
}

template <typename T, int n> bool VecCmp(Vec<T,n> &A, T b, bool lessThan = true)
{
	bool rst = true;
	if (lessThan)
	{
		for (int i=0; i< n; i++)
			rst = rst && (A[i] < b);
		if (rst)
			return true;
		else
			return false;
	}
	else
	{
		for (int i=0; i< n; i++)
			rst = rst && (A[i] > b);
		if (rst)
			return true;
		else
			return false;
	}
}

template <typename T, int n> Vec<T,n> VecAbs(Vec<T,n> &A)
{
	Vec<T,n> C;
	for (int i=0; i< n; i++)
		C[i] = abs(A[i]);
	return C;
}

EXPORTLIB double AverageScalar3(cv::Scalar s);

EXPORTLIB unsigned long factorial(unsigned long n);
EXPORTLIB double factorial(int n);
//STBlock FreqComplexFilter(STBlock& src, cv::Mat& kernel);
EXPORTLIB cv::Mat FreqComplexFilter(const cv::Mat& src, cv::Mat& kernel, bool conj=false);

//template<typename T, size_t cn> Tensor_<T,cn> FreqComplexFilter(const Tensor_<T,cn>& src, cv::Mat& kernel)
//{
//	Tensor_<T,cn> rst(src.size());
//	if ( src.channel() != kernel.channels())
//		cout<<"error, the src and kerenl are not with same channels\n";
//	else if (src.channel()== 2 || src.channel()==6)
//	{
//		cv::Mat temp;
//		int i=0;
//		Tensor_<T,cn>::c_iter_type citer;
//		for (citer = src.begin(); citer!= src.end(); citer++)
//		{
//			cv::mulSpectrums(kernel,*citer,temp,CV_DXT_ROWS);
//			rst.SetFrame(i, temp);
//			i++;
//		}
//	}
//	else
//		cout<<"error, the input are not complex image\n";
//
//	return rst;
//
//}
template<typename _T> cv::Scalar complexToScalar(std::complex<_T> val)
{
	cv::Scalar_<_T> rst;
	rst[0] = val.real();
	rst[1] = val.imag();
	return rst;
}

template<typename _T> cv::Complex<_T> vecToComplex(const cv::Vec<_T,2>& val)
{
	return cv::Complex<_T>(val[0],val[1]);
}

template<typename _T> cv::Vec<_T,2> complexToVec(const cv::Complex<_T>& val)
{
	return cv::Vec<_T,2>(val.re,val.im);
}
/*
EXPORTLIB Printable_t ParseMyArgs(va_list& ap)
{
  int count=0;
  while (*ap)
  {
     if (count%2==0)
    {
      Printable.s = va_arg(ap,char*);
      cout<<Printable.s<<endl;;
    }
    else
    {
      Printable.i= va_arg(ap,int);
      cout<<Printable.i<<endl;
    }
    count++;
  }}
  */
EXPORTLIB cv::Mat readMatFromTxt(string filename, int H=0, int W=0);
EXPORTLIB cv::Mat complexMul(const cv::Mat& A, const cv::Mat& B);
EXPORTLIB cv::Mat toComplex(const cv::Mat& A, const cv::Mat& B);
EXPORTLIB cv::Mat toComplex(const cv::Mat& A);
EXPORTLIB cv::Mat getChannel(const cv::Mat& im, int n);
EXPORTLIB int combination(int fromNum, int selectNum);

EXPORTLIB void CombineImage(const vector<string>& infilenames, string& outfilename, int step = 5);

//#include "Tensor.h"
//20130816 enum class FilterBoundary :int {FILTER_BOUND_EXTEND,FILTER_BOUND_SAME, FILTER_BOUND_VALID};

//template<typename T, size_t cn> Tensor_<T,cn> Filter2D(const Tensor_<T,cn> & ts, const Tensor_<T,cn> & ker)
//{
//	Tensor_<T,cn> rst;
//	return rst;
//}

//
//template<typename T, size_t cn> Tensor_<T,cn> Filter2D(const Tensor_<T,cn> & ts, const Tensor_<T,cn>& ker, int boundary)
//{
//	if (ker.size().depth!=1 && ker.size()().depth!= ts.size().depth)
//		CV_Errorr(CV_StsNotImplemented, "unsupport 3D filtering");
//	else
//	{ 
//		Tensor_<T,cn> rst(ts.size());
//		for (int z=0; z< ts.size().depth; z++)
//		{
//			if (ker.size().depth ==1)
//				cv::filter2D(ts[z],rst[z],-1,ker[0]);
//			else
//				cv::filter2D(ts[z],rst[z],-1,ker[z]);
//		}
//		if (boundary == FILTER_BOUND_EXTEND)
//		{
//			CV_Error(CV_StsNotImplemented,"unsupport extend boundary version");
//			return rst;
//		}
//		else if (boundary == FILTER_BOUND_VALID)
//		{
//			Size3 validSize = ts.size() - ker.size() + Size3(1,1,1);
//			Point3i validOffset = ker.size()/2;
//			return rst.Crop(validOffset,validSize);
//		}
//		else
//			return rst;
//	}
//}

	
    
}



#endif

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
#define _CRT_SECURE_NO_WARNINGS
#ifndef SIZE3_H
#define SIZE3_H
//#include <cv.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <iostream>
using namespace std;
using namespace cv;

#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
namespace tensor{
template <typename Tp> class EXPORTLIB Size2_ : public Size_<Tp>
{
public:
	//the new 2D block size, the different to cv::Size_ is
	//the order ot initialization, is (height, width)
	//that is more consist with the order of 2D point (x, y);
	Size2_<Tp>();//defalut constructor
	Size2_<Tp>(Tp height, Tp width);
	Size2_<Tp>(const Size2_<Tp>& sz);//copy
	Size2_<Tp>(const Size_<Tp>& sz);//convert from cv::Size
	Size2_<Tp>(const Point_<Tp>& pt);//convert form cv Point
	Size2_<Tp>& operator= (const Size2_<Tp>& sz);
	Size2_<Tp>& operator= (const Size_<Tp>& sz);
	//inherite var Tp width and Tp heigth from cv::Size_
	//conversion form another data type
	template<typename Tp2> operator Size2_<Tp2>() const;
	//calc area
	Tp area(void) const;
	//arithmatic operators
	Size2_<Tp> operator* (Tp b);
	Size2_<Tp> operator* (const Size2_<Tp>& b);
	Size2_<Tp> operator/ (Tp b);
	Size2_<Tp> operator/ (const Size2_<Tp>& b);
	Size2_<Tp> operator+ (const Size2_<Tp>& b);
	Size2_<Tp> operator- (const Size2_<Tp>& b);
	bool operator== (const Size2_<Tp>& b);
	Size2_<Tp> Floor() const;
	Size2_<Tp> Ceil() const;
};

/////////////////////////////////////Size3_<Tp>/////////////////////////////
template <typename Tp> class EXPORTLIB Size3_ : public Size2_<Tp>
{
public:
	//3D block size is inherited from Size2_
	//the order of init. is (height, width, depth)
	 Size3_<Tp>();
	Size3_<Tp>(Tp height, Tp width, Tp depth);
	Size3_<Tp>(const Size3_<Tp>& sz);
	//conver from 3D point
	Size3_<Tp>(const Point3_<Tp>& pt);
	Size3_<Tp>& operator = (const Size3_<Tp>& sz);
	//can trivally convert from Size_2
	Size3_<Tp>(const Size2_<Tp>& sz2);
	//convert from cv::Size
	Size3_<Tp>(const Size_<Tp> & cvsz);
	//conversion from another data type
	template<typename Tp2> operator Size3_<Tp2>() const;
	//calc volumn
	Tp volumn() const;
	Size3_<Tp> operator* (Tp b);
	Size3_<Tp> operator* (const Size3_<Tp>& b);
	Size3_<Tp> operator/ (Tp b);
	Size3_<Tp> operator/ (const Size3_<Tp>& b);
	Size3_<Tp> operator+ (const Size3_<Tp>& b);
	Size3_<Tp> operator- (const Size3_<Tp>& b) const;
	bool operator== (const Size3_<Tp>& b);
	//Point3_<Tp> operator- (const Point3_<Tp>& p);
	Point3_<Tp> Point3() const;
	//members
	Tp depth;
	void Print(void) const;
	Size3_<Tp> Floor() const;
	Size3_<Tp> Ceil() const;
};	

//this two type conversion are very special, i cannot dllexport them, don't know the reason
//will leave for future research. 

template<typename Tp> template<typename Tp2> inline Size2_<Tp>::operator Size2_<Tp2>() const
{
	return Size2_<Tp2>(saturate_cast<Tp2>(this->height),saturate_cast<Tp2>(this->width));
}

template<typename Tp> template<typename Tp2> inline Size3_<Tp>::operator Size3_<Tp2>() const
{
	return Size3_<Tp2>(saturate_cast<Tp2>(this->height),saturate_cast<Tp2>(this->width),saturate_cast<Tp2>(this->depth));
}
typedef	Size3_<int> Size3;
typedef Size2_<int> Size2;
}
#endif


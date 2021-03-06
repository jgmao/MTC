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
#define _CRT_SECURE_NO_DEPRECATE
#ifndef CUBE_H
#define CUBE_H
//#include <cv.h>
#include <opencv2/opencv.hpp>
#include "Size3.h"

#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
using namespace std;

namespace tensor{
template <typename Tp> class Cube_ : public Point3_<Tp>, public Size3_<Tp>
{
public:
	EXPORTLIB Cube_(void);
	EXPORTLIB Cube_(Tp x, Tp y, Tp z, Tp height, Tp width, Tp depth);
	EXPORTLIB Cube_(const Point3_<Tp>& pos, const Size3_<Tp>& sz);
//	EXPORTLIB Cube_(const Point_<Tp>& pos, const Size_<Tp>&sz);
	EXPORTLIB Cube_(const Rect_<Tp>& r);
	EXPORTLIB Cube_(const Cube_<Tp>& c);
	EXPORTLIB Cube_& operator = (const Cube_<Tp>& c);
	EXPORTLIB Point3_<Tp> offset(void) const;
	EXPORTLIB Size3_<Tp> size(void) const;
	EXPORTLIB Rect_<Tp> toRect(void) const;
	EXPORTLIB ~Cube_(void);
	EXPORTLIB void Print(void) const;
};

typedef Cube_<int> Cube;

///////////////////////////////////////
}
#endif

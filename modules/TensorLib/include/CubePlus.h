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

#ifndef CUBEPLUS_H
#define CUBEPLUS_H

#include "Cube.h"
#include "TensorLite.h"
namespace tensor{
template <typename T, size_t cn> class CubePlus_ : public Cube
{
protected:
	Tensor<T,cn> content;
  Tensor<T,cn> extraContent;
public:
	EXPORTLIB CubePlus_(void);
	EXPORTLIB CubePlus_(const Cube& c);
	EXPORTLIB CubePlus_(const Cube& c, const Tensor<T,cn>& src);
	EXPORTLIB CubePlus_(const Tensor<T,cn>& t);
  EXPORTLIB CubePlus_(const Tensor<T,cn>& t, const Tensor<T,cn>& extra);
        EXPORTLIB CubePlus_& operator= (const CubePlus_& c);
        EXPORTLIB CubePlus_& operator= (const Cube& c);
        EXPORTLIB void SetContent(const Tensor<T,cn>& t);
  EXPORTLIB void SetExtraContent(const Tensor<T,cn>& t);
        EXPORTLIB const Tensor<T,cn>& GetContent(void) const;
  EXPORTLIB const Tensor<T,cn>& GetExtraContent(void) const;

};


typedef CubePlus_<double,1> CubePlus;
}
#endif

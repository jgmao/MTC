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
#ifndef COMMON_DEF_H
#define COMMON_DEF_H
#include <opencv2/opencv.hpp>
#include <opencv2/gpu/gpu.hpp>
#include <vector>
namespace tensor{
//Marcos
  //debug
  //typedef uchar uchar;
  #ifndef DEBUG_SYMBOL
  #define DEBUG_SYMBOL
  #define DEBUG_X 848//27*32
  #define DEBUG_Y 176//
  #define DEBUG_DISP_X -128
  #define DEBUG_DISP_Y 0
  #define DEBUG_SIZE 16
  /*extern int DEBUG_X ;
  extern int DEBUG_Y;
  extern int DEBUG_DISP_X;
  extern int DEBUG_DISP_Y;
  extern int DEBUG_SIZE;*/
  #endif
#define DEBUG_20131226	0
#ifndef USE_GPU
//#define USE_GPU
#endif

#define PARALLEL_MATCHING   0

#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif



#ifndef BUFFERGPU
#define BUFFERGPU
class BufferGPU// genearal structure
{
public:
  cv::gpu::GpuMat gI1,gI2; //input
  cv::gpu::GpuMat t1,t2; //intermediate
  cv::gpu::GpuMat gs; //rst
  std::vector<cv::gpu::GpuMat> v;
  cv::gpu::Stream stream;
};

#endif
// 20130816 only define the Marcors and constants used in TensorLib
#ifndef FILTER_BOUNDARY
#define FILTER_BOUNDARY
enum class FilterBoundary :int {FILTER_BOUND_EXTEND,FILTER_BOUND_SAME, FILTER_BOUND_VALID, FILTER_BOUND_FULL, FILTER_BOUND_HALF};
#endif

#ifndef COMPARE_CRITERIA
#define COMPARE_CRITERIA

enum class CompareCriteria :int {COMPARE_CRITERIA_MSE,
			    COMPARE_CRITERIA_SSIM_DCT,
			    COMPARE_CRITERIA_SSIM ,
			    COMPARE_CRITERIA_INTERP ,
			    COMPARE_CRITERIA_SUBJECTIVE ,
			    COMPARE_CRITERIA_SAD ,
			    COMPARE_CRITERIA_MAHALANOBIS ,
			    COMPARE_CRITERIA_LRI ,
			    COMPARE_CRITERIA_SVM,
			    COMPARE_CRITERIA_OTHER};
#endif

#ifndef PRINTABLE
#define PRINTABLE


  union Printable_t {
         int     i;
         float   f;
         double  d;
         char    c;
         char   *s;
         bool    b;
         cv::Mat *m;
         Printable_t(int a_):i(a_){};
  };

#endif

#ifndef QUALITY_TEST
#define QUALITY_TEST
enum class MetricModifier :int {STSIM2_BASELINE,
                                SSIM2_DCT,
                                PQI_METRIC,
                                STSIM2_ADT_NAIVE,
                                STSIM2_NEW_L1, STSIM2_NEW_L2, STSIM2_NEW_L3, STSIM2_NEW_L4,
                                MSE_BASELINE,
                                SAD_BASELINE,
                                MAHALANOBIS_DIST,
                                LRI_METRIC,
                                SE_MSE,
                                STSIM2_SE_MSE,
                                STSIM2_TUNE,
                                SVM_METRIC,
                                STSIM3_LSE,
                                STSIM2_PART,//partial STSIM2 for fast
                                Q_UNDEF};

#ifndef STSIM2_TYPE
#define STSIM2_TYPE

enum class FeaturePoolType :int {FEATURE_POOL_AVE,FEATURE_POOL_MIN,FEATURE_POOL_ALL};

#endif


#ifndef BOUND_DIR
#define BOUND_DIR
enum class BoundDir :int {UP,DOWN,LEFT,RIGHT,FRONT,BACK,NONE};//may be extern?
#endif

#endif
#ifndef FOOT_ITEM
#define FOOT_ITEM
struct FootItem
{
    double value;
    int bits;
    std::string bitstring;
};
class ComparePoint3i //class used for map<Point3i, FootItem> ordering
{
public:
  EXPORTLIB bool operator() (const cv::Point3i& a, const cv::Point3i& b) const
  {
    //if (x.z < y.z)
    //  return true;

    if (a.x < b.x)
      return true;
    else if (a.x == b.x)
    {
      if (a.y < b.y)
        return true;
      else
        return false;
    }
    else
      return false;

  };

};
#endif

#ifndef DIRECTION
#define DIRECTION

enum class Directions	 : int {DIRECTION_VERTICAL,
									DIRECTION_HORIZONTAL,
									DIRECTION_CENTER,
									DIRECTION_OTHER};

enum class CornerPos: int {CORNER_NW,
							  CORNER_SE};
#endif

#ifndef SOURCE_CODING_METHOD
#define SOURCE_CODING_METHOD
enum class SrcCodingMethod : int { JPEG_HUFF, HUFFMAN_CODE, UNARY_CODE};

#endif


#ifndef CODING_METHOD
#define CODING_METHOD
enum  class CodingMethodNames : int {CODING_JPEG, CODING_MTC, CODING_PQI, CODING_JPEG_DEGRADE, CODING_OTHER, NO_CODING};
#endif


#ifndef BLENDING_METHOD
#define BLENDING_METHOD
enum class BlendingMethod : int {SHORTEST_PATH_BLENDING, GRADIENT_BLENDING, GRAPHCUT_BLENDING,NO_BLENDING};
enum class BlendingLocation : int { FORWARD_BLENDING, POST_BLENDING_RIGHT, POST_BLENDING_LOW, CUSTOM_BLENDING};
#endif

#ifndef SIDE_MATCHING_METHOD
#define SIDE_MATCHING_METHOD
enum class MatchingMethod : int { MATCHING_MSE, MATCHING_SAT, MATCHING_VAR, MATCHING_SAD, MATCHING_MSE_CONSTRAINT,MATCHING_HIERARCHY, MATCHING_DIRECT, MATCHING_OPENCV, MATCHING_STSIM, MATCHING_STSIM_PART};
#endif
//!20131010 disable this since STSIM2 do not have power to detect lighint erro
#define SUBSTRACT_BOUND_MEAN 0 //substract mean when do matching

}
#endif


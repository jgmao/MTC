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
#ifndef QNODE_H
#define QNODE_H


//#define INT_RST
//#define INT_QUAN
#include "CommonDef.h"
#include "TensorLite.h"
#include "LinkQueue.h"
#include "Lighting.h"
#include "PostBlendSet.h"
#include "PoissonSolver.h"
using namespace lighting;
namespace tensor{
template< class T, size_t cn> class QNode :	public Tensor<T,cn>
{
public:
	typedef typename Tensor<T,cn>::value_type value_type;
	typedef vector<Tensor<T,cn>> bound_type;
	typedef bound_type& bound_ref;
  typedef Vec<T,cn> data_type;
	typedef vector<vector<vector<double>>> Matrix3Df;
	typedef vector<vector<vector<int>>> Matrix3Di;

    EXPORTLIB QNode<T,cn>(void);
    EXPORTLIB ~QNode<T,cn>(void);
    EXPORTLIB QNode<T, cn>(const cv::Mat &mt);
    EXPORTLIB QNode<T,cn>(const QNode& nd);
    //EXPORTLIB QNode(const QTree& qt);
    EXPORTLIB QNode<T,cn>(const Size3& sz, const Vec<T,cn>& val = data_type::all(0));
    EXPORTLIB QNode<T,cn>(const Size3& sz, const Size3& overlapSize);
    EXPORTLIB QNode<T,cn>(const Size3& sz, double overlapRatio);
    //EXPORTLIB QNode(const QTree& qt, const Size3& overlapSize);
	//ref init from a tensor
  EXPORTLIB QNode<T,cn> (const Tensor<T,cn>& ts, const Size3& sz, const Point3i& pos, const Size3& boundarySize = Size3(1,1,0));
    EXPORTLIB QNode<T,cn> & operator= (const QNode<T,cn>& nd);
    EXPORTLIB void SetOverlapSize(const Size3& overlapSize);
    EXPORTLIB QNode<T,cn>(string cFileName);
    EXPORTLIB Size3 overlap(void) const;
    EXPORTLIB void SetFoot(const Vec<T,cn>& ft);
    EXPORTLIB void SetFoot(const Tensor<T,cn>& src);
  EXPORTLIB void AddFoot(const pair<Point3i,FootItem>& f);
  EXPORTLIB void clearFoot(void);
  EXPORTLIB vector<pair<Point3i,FootItem> > getFeetCopy(void) const;
  EXPORTLIB int getFeetNumber(void) const;
  EXPORTLIB void SetUpFoot(const value_type& up);
  EXPORTLIB void SetLeftFoot(const value_type& left);
    EXPORTLIB Vec<T,cn> GetFoot() const;
  EXPORTLIB pair<Point3i,FootItem> GetFoot(int i) const;
  EXPORTLIB void SetFoot(const pair<Point3i,FootItem>& f, int i);
  EXPORTLIB Vec<T,cn> GetApproxFoot() const;
public:
    EXPORTLIB vector<vector<int>>& Quilting(QNode<T,cn>& ts);
    EXPORTLIB vector<vector<int>>& Quilting(vector<Tensor<T,cn>>& boundA, vector<Tensor<T,cn>>& boundB, BoundDir dir=BoundDir::NONE);
    EXPORTLIB QNode<T,cn>& LinearInterp(double qsize,bool computeFoot=true, Directions direction = Directions::DIRECTION_OTHER);
    EXPORTLIB QNode<T,cn>& LinearInterpGeneric(double qsize);
    //EXPORTLIB void SetFootRegion(int region=0);
    //EXPORTLIB void SetFootComputeMethod(int method=0);
    EXPORTLIB int BitLength(void) const;
    EXPORTLIB string GetBits(void) const;
    EXPORTLIB Tensor<T,cn>& GetBoundaryUp(void);
    EXPORTLIB const Tensor<T,cn>& GetBoundaryUp(void) const;
    EXPORTLIB Vec<T,cn>& GetBoundaryUp(const Point3i& pos);
    EXPORTLIB Tensor<T,cn>& GetBoundaryLeft(void);
    EXPORTLIB const Tensor<T,cn>& GetBoundaryLeft(void) const;
    EXPORTLIB Vec<T,cn>& GetBoundaryLeft(const Point3i& pos);
    EXPORTLIB Vec<double,cn> Quantize(const Vec<double,cn>& val, double qsize, SrcCodingMethod coding_method = SrcCodingMethod::JPEG_HUFF);
    EXPORTLIB Vec<int, cn> UniformQuantize(const value_type& val, double qsize);
    EXPORTLIB Tensor<T,cn>& GetBoundary(int i);
    EXPORTLIB const Tensor<T,cn>& GetBoundary(int i) const;
    EXPORTLIB int GetBoundarySize(void) const;
    EXPORTLIB Point3i GetFootPos(void) const;
  EXPORTLIB Point3i GetFootPos(void);
    EXPORTLIB void Ref(const Tensor<T,cn>& ts, const Cube& roi, const Size3& boundarySize);
    EXPORTLIB Tensor<T,cn> LightingCorrection(const Tensor<T,cn>& changeTo, bool saveCodeLength=true);
    EXPORTLIB QNode<T,cn>& LightingCorrection2(const QNode<T,cn>& changeTo);
    EXPORTLIB Tensor<T,cn> LightingCorrection(const Tensor<T,cn>& changeTo, const Tensor<T,cn>& VQCodebook);
    EXPORTLIB QNode<T,cn>& PoissonLightingCorrection(const QNode<T,cn>& changeTo,const Tensor<T,cn>& ref, const Tensor<T,cn>& rec, Size3 boundSize = Size3(0,0,0), int footRegion=0, int footMethod=0, double qsize=0);
    //EXPORTLIB QNode<T,cn>& PostPLC(const QNode<T,cn>& changeTo);


    EXPORTLIB void GradientStitching(const QNode<T,cn> & changeTo, BlendingLocation blendPos, Size3 boundSize,const  Tensor<T,1>& tempmask);


    EXPORTLIB void LightingCorrection2(const Tensor<T,cn>& fromLight, const Tensor<T,cn>& toLight);
    EXPORTLIB QNode<T,cn> Clone() const;
    EXPORTLIB void ComputeFoot(double qsize, Directions direction= Directions::DIRECTION_OTHER, SrcCodingMethod coding_method = SrcCodingMethod::JPEG_HUFF);//direction_other means take both horizontal and vertical
    EXPORTLIB Vec<T,cn> ComputePLCFoot(const Tensor<T,cn>& ref, const Point3i& pos, const Size3& sz, int caseNo=1, int method=1);
    EXPORTLIB void InitBound(double qsize); //init. node boundary when reach the image bound (upper and left)
    EXPORTLIB Tensor<T,cn> GetExtendTensor(bool up=1, bool left = 1, bool down= 1, bool right =1) const;
  EXPORTLIB Tensor<T,cn> GetExtendTensor(Vec<T,cn> padding, bool up=1, bool left = 1, bool down= 1, bool right =1) const;

private:
	void MakeVeryFirstNode(double qsize);
	void InterpZeroRow(void);
	void InterpZeroCol(void);
	bool AssertExtend(const QNode& extNode, int n =1) const;
	void InterpLastRow(void);
	void InterpLastCol(void);
	void InterpCenter(void);

	int VariableLengthCoding(int index);

	int UnaryCoding(int index);

	//////////////////////////
	vector<vector<int>> Blending(Matrix3Di path, QNode<T,cn>& ts, int criteria =0); // ts is the other patch which will blend to 
	vector<vector<int>> Blending(Matrix3Di path, vector<Tensor<T,cn>>& boundA,vector<Tensor<T,cn>>& boundB, BoundDir dir);
	Matrix3Di FindPath(const QNode<T,cn>& ts);
	Matrix3Di FindPath(const vector<Tensor<T,cn>>& boundA, vector<Tensor<T,cn>>& boundB,BoundDir dir=BoundDir::NONE);
	double ShortestPathDP(int i, int j, const cv::Mat &emx, cv::Mat &dpmx, vector<vector<int>> &pathPlane, Directions direction= Directions::DIRECTION_VERTICAL);
	void TreatCrossRegion(Matrix3Di & path, CornerPos corner = CornerPos::CORNER_NW);
	void MakePath(vector<vector<int>>& pathPlane, cv::Mat &dpmx, Directions direction);
protected:
	Size3 overlapSize;
	Vec<double,cn> footError;
	Vec<double,cn> nodeFoot;
	Vec<double,cn> approxFoot;
  Vec<double,cn> leftFoot;
  Vec<double,cn> upFoot;
  Vec<double,cn> antiFoot;

	Point3i footPos;
	string quanBit;
	int bits;
	CodingMethodNames predicted_method; //0, not predicted, 1, quilted, 2, PQI
	double qsize;//quantize size
	bound_type boundary;//for 2D still image, only two bound exist
	//int footComputeRegion;
	//int footComputeMethod;
  value_type upNodeFoot;
	value_type leftNodeFoot;
  vector<pair<Point3i,FootItem>> feet;
public:
 Lighting lt;
public: ///change it for debug use,!!!!! change it back for future
	Tensor<T,cn> leftBound, upBound, lowBound, rightBound;
	vector<Tensor<T,cn>*> bounds;
public:
	vector<vector<int>> seam; // the position of the ridge
   Tensor<T,cn> eNode;//extended node
    
};
}
#endif

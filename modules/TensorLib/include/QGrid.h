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
#ifndef QGRID_H
#define QGRID_H

#include "QTree.h"
#include <thread>
#include <memory>
#include "algorithms.h"
#include "Laplace.h"

namespace tensor{
template<class T, size_t cn>
class QGrid
{
public:
	typedef vector<vector<vector<QTree<T,cn>>>> Grid;
	EXPORTLIB QGrid(void);
	EXPORTLIB QGrid& operator=(const QGrid& g);
	EXPORTLIB QGrid(const QGrid& g);
	//QGrid(const Size3& blockSize);
	EXPORTLIB QGrid(const string& cFileName, const Size3& blockSize, const Size3& overlapSize);
	EXPORTLIB QGrid(const string& cFileName, const Size3& blockSize, const Size3_<double>& overlapRatio);
	EXPORTLIB QGrid(const string& cFileName);
	EXPORTLIB void Display(int flag=1);
	EXPORTLIB void DisplayCasualMap(int flag=1);
	EXPORTLIB void SetCausalMap(const QNode<T,cn>& nd);
	EXPORTLIB bool IsInsideCausalRegion(const Point3i& curPos, const QNode<T,cn>& nd, const int multiBound=3) const;
	EXPORTLIB ~QGrid(void);
	EXPORTLIB QTree<T,cn>& GetNode(const Point3i& pos) ;
	EXPORTLIB void SetNode(const Point3i& gridPos, const Tensor<T,cn>& ts, const Cube& roi, const Size3& overlapSize);
	EXPORTLIB vector<pair<Point3i, double> > BoundaryMatching(QNode<T,cn>& qNode,MatchingMethod matching_method = MatchingMethod::MATCHING_MSE, double matching_thrd = 0/*accept all*/, Size3 subWinSize=Size3(16,16,1));
	EXPORTLIB void ReInitGrid(void);
	EXPORTLIB Size3 GetGridSize(void) const;
	EXPORTLIB void SetVarThrd(double t1 =1300);
	EXPORTLIB void UpdateRstSqr(const Point3i& sPos,const Point3i& ePos);
	EXPORTLIB void UpdateRstSqr(const Cube& cb);
	//EXPORTLIB Size3 GetSubWinSize(void) const;
	//EXPORTLIB Size3 GetSubWinStep(void) const;
	//EXPORTLIB void SetSubWinSize(const Size3& sz);
	//EXPORTLIB void SetSubWinStep(const Size3& sz);
  EXPORTLIB void ComputeBoundMean(const Tensor<T,cn>& tarUp, const Tensor<T,cn>& tarLeft, const Point3i& offset, const Size3& size, const Size3& overlap,Tensor<T,cn>& upMean, Tensor<T,cn>& leftMean);
private:
	void InitGrid(void);

protected:
	Laplace L1Model[2];//left H0, up H0, left H1, up H1
	int lenH0,lenH1;
	std::list<bool> queueLen;
	const uint L1_train_max =100;
	int L1_train_len;
	bool train;
	//Size3 subWinSize;
	//Size3 subWinStep;
	string cFileName;
	Size3 blockSize; //init block size
	Size3 overlapSize; // int overlap size
	//Point3i curPos;
	Size3 gridSize;
	Grid grid;
	Tensor<T,cn> ensemble;
	Tensor<T,cn> rst;
  Tensor<T,cn> ensembleExt;
  //Tensor<T,cn> rstExt;
	Tensor<T,1> causalMap;
	int candidNum;
	Size3 searchStep; //matching step size
	Tensor<double,cn> localVarMap;
	double varThrd1;
	Tensor<T,cn> rstSqr;
  Tensor<T,cn> rstBorderSqr;
  std::mutex mymux;

};
}

#endif

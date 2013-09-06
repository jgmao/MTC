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
#ifndef QTREE_H
#define QTREE_H

#include "QNode.h"
namespace tensor{
template<class T, size_t cn> class QTree : public QNode<T,cn>
{
public:
	typedef typename std::vector<QTree<T,cn> >::iterator qiter;
	typedef typename std::vector<QTree<T,cn> >::const_iterator citer;
  typedef Vec<T,cn> data_type;
	EXPORTLIB QTree(void);
	EXPORTLIB ~QTree(void);
	EXPORTLIB QTree(const cv::Mat& mt);
	EXPORTLIB QTree(const Tensor<T,cn>& ts);
	EXPORTLIB QTree(const QTree& qt);
	EXPORTLIB QTree(string cFileName);
	EXPORTLIB QTree(const Size3& sz, const Vec<T,cn>& val=data_type::all(0));
	EXPORTLIB QTree& operator= (const QTree& qt);
  EXPORTLIB QTree(const QNode<T,cn>& nd);
	//template<class T2> operator QTree<T2,cn>() const;
	/////////
	EXPORTLIB void InitEntries(void);
	EXPORTLIB void InitEntries(const QTree& qt);
public:
	EXPORTLIB qiter GetSubTree(int pos);
	//qiter GetSubTree(int lv, int pos)const;
	EXPORTLIB qiter GetSubTreeBeginning(void) ;
	EXPORTLIB qiter GetSubTreeEnding(void) ;
	EXPORTLIB int GetSubTreeSize(void) const;
	EXPORTLIB int GetPeerTreeSize(void) const;
	EXPORTLIB int GetTreeLevel(void) const;
	EXPORTLIB int GetTreePos(void) const;
	EXPORTLIB void SetTreeLevel(int lv);
	EXPORTLIB void SetTreePos(int pos);
	EXPORTLIB void Split(void);
	EXPORTLIB void Split(bool direction);//0 vertical 1 horizontal
	EXPORTLIB QTree * NextLeaf(void);
	//EXPORTLIB QTree<T,cn> Clone() const;
	//////////////////////
public:
	bool splitMark;
protected:
	int branchNum;
	int treeLevel;
	int treePos;
	vector<QTree> subTrees;
	qiter subTreeIter;
	QTree * parentTree;
	QTree * nextTree;
public:
  class BitCount
  {
  protected:
    CodingMethodNames coding_method;
    //int coding_bit_len;
    string coding_bit;
    string foot_bit;
    string dec_bit;
    //int tree_bit;
    string bitstring;
  public:
    BitCount(){foot_bit="";coding_method=CodingMethodNames::NO_CODING;coding_bit="";dec_bit="";bitstring="";}
    ~BitCount(){}
    void SetFootBit(string str){foot_bit=str;AppendString(str);}
    void SetCodeMethod(CodingMethodNames method){coding_method=method;}
    void SetCodeBit(string str){coding_bit = str;AppendString(str);}
    void SetDecBit(string str){dec_bit=str; AppendString(str);}
    void AddCodeBit(string str){coding_bit+=str;AppendString(str);}
    int GetDecBitN(void)const{return dec_bit.length();}
    int GetFootBitN(void)const{return foot_bit.length();}
    int GetCodingBitN(void)const{return coding_bit.length();}
    CodingMethodNames GetCodingMethod(void)const{return coding_method;}
  private:
    void AppendString(string str){bitstring+=str;}
  }bitcount;
  EXPORTLIB int CollectBits(void) const;
  EXPORTLIB int CollectTreeBits(void) const;
  EXPORTLIB int CollectModeBits(void) const;
  EXPORTLIB vector<int> CollectMTCBits(int level=0) const;
  EXPORTLIB vector<int> CollectJPEGBits(void) const;
  EXPORTLIB int CollectPQIBits(void) const;
  EXPORTLIB vector<int> GetMTCN(int level=0) const;
  EXPORTLIB vector<int> GetJPEGN() const;
  EXPORTLIB vector<int> GetPQIN(int level=0) const;
};

//template<class T, size_t cn> template<class T2>
//QTree<T,cn>::operator QTree<T2, cn>() const
//{
//	Tensor<T2,cn> rst(Tensor<T2,cn>());
//	return rst;
//}
}
#endif

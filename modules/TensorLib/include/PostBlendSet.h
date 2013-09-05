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

#ifndef POSTBLENDSET_H
#define POSTBLENDSET_H


#include <map>
#include "CubePlus.h"
#include "Size3.h"

namespace tensor{

using namespace cv;
class PBRecord
{
public:
	Point3i offset;
	BoundDir direction;
	bool operator< (const PBRecord& rhs) const
	{
		if (offset.x < rhs.offset.x)
			return true;
		else if ((offset.x == rhs.offset.x)&&(offset.y < rhs.offset.y))
			return true;
		else if ((offset == rhs.offset)&&(direction<rhs.direction))
			return true;
		else
			return false;
	}
	PBRecord():offset(Point3i()),direction(BoundDir::UP){}
	PBRecord(const Point3i& pos, BoundDir dir):offset(pos),direction(dir){}
};

class PBSet:public map<PBRecord,CubePlus>
{
	typedef pair<PBRecord,CubePlus> value_type;
public:
	void insert ( const value_type& x )
	{
		iterator it = find(x.first);
		if (it!=end())
		{
			if(it->second.size()!=x.second.size())
				erase(x.first,x.second);
		}
		map<PBRecord,CubePlus>::insert(x);
	}
	//bool isB2sA(const key_type& key, const Cube& val)
	//{
	//	iterator it = find(key);
	//	if (val.height<it->second.height && it->first.direction==LEFT)
	//		return true;
	//	else if (val.width < it->second.width && it->first.direction==UP)
	//		return true;
	//	else
	//		return false;
	//}
	//
	void erase(const key_type& key, const CubePlus& val)
	{
		iterator it = find(key);
		if (it!=end())
		{
      //int extSize = max(it->second.GetExtraContent().size().height,it->second.GetExtraContent().size().width);
      //int osize = min(it->second.GetExtraContent().size().height,it->second.GetExtraContent().size().width);
			//big to small case A 
			if (val.height< it->second.height && it->first.direction==BoundDir::LEFT) //Big to small case A - vertical
				this->operator[](PBRecord(it->first.offset+Point3i(val.height,0,0),it->first.direction)) = CubePlus(it->second.GetContent().Crop(Point3i(val.height,0,0),it->second.size()/2+Size3(0,0,1)),it->second.GetExtraContent()/*.Crop(Point3i(extSize/2-1,0,0),Size3(extSize/2+1,osize/2+1,1))*/);//-Size3(val.height,0,0));
			if (val.width< it->second.width && it->first.direction == BoundDir::UP)//Big to small case A - horizontal
				this->operator[](PBRecord(it->first.offset+Point3i(0,val.width,0),it->first.direction))=CubePlus(it->second.GetContent().Crop(Point3i(0,val.width,0),it->second.size()/2+Size3(0,0,1)),it->second.GetExtraContent()/*.Crop(Point3i(0,extSize/2-1,0),Size3(osize/2+1,extSize/2+1,1))*/); //also update extra boundary
				//(val.offset()+Point3i(0,val.width,0),it->second.size()/2+Size3(0,0,1));//-Size3(0,val.width,0));
			map<key_type,CubePlus>::erase(it); //same size and small to big case A are nothing special
		}
		else 
		{

			//big to small case B
			it = find(key_type(key.offset-Point3i(val.height,0,0),key.direction)); //vertical
			if (it!=end() && it->second.height > val.height)
				it->second = CubePlus(it->second.GetContent().Crop(Point3i(/*val.height*/0,0,0),it->second.size() - Size3(val.height,val.width,0)),it->second.GetExtraContent()/*.Crop(Point3i(it->second.GetExtraContent().size().height/2-1,0,0),it->second.GetExtraContent().size()/2 + Size3(1,1,1))*/);
				//Cube(it->second.offset()+Point3i(val.height,0,0),it->second.size() - Size3(val.height,val.width,0));
			it = find(key_type(key.offset-Point3i(0,val.width,0),key.direction)); //horizontal
			if (it!=end() && it->second.width > val.width)
				it->second = CubePlus(it->second.GetContent().Crop(Point3i(0,/*val.width*/0,0),it->second.size() - Size3(val.height,val.width,0)),it->second.GetExtraContent()/*.Crop(Point3i(0,it->second.GetExtraContent().size().width/2-1,0),it->second.GetExtraContent().size()/2 + Size3(1,1,1))*/);
				//Cube(it->second.offset()+Point3i(0,val.width,0),it->second.size() - Size3(val.height,val.width,0));
			//small to big case B
			int i=val.height/2;
			while(i>0)
			{
				it = find(key_type(key.offset+Point3i(i,0,0),key.direction));//vertical
				if (it!=end()&&it->second.height < val.height)
				{
					map<key_type,CubePlus>::erase(it);
					break;
				}
				i=i/2;
			}

			i=val.width/2;
			while(i>0)
			{
				it = find(key_type(key.offset+Point3i(0,i,0),key.direction));//horizontal
				if (it!=end()&&it->second.width < val.width)
				{
					map<key_type,CubePlus>::erase(it);
					break;
				}
				i=i/2;
			}
		}

	}

};
}
#endif

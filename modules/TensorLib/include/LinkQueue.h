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

#include <iostream>
#include <vector>
//#include <cv.h>
#include <opencv2/opencv.hpp>
#include <mutex>
using namespace std;
struct CandidateRecord
{
  double data;
  cv::Point3i addr;
};
class LinkStruct
{
public:
   std::mutex mymutex;
  LinkStruct(){};
  ~LinkStruct(){};
  virtual int getLength()=0;
  virtual void push(double d,const cv::Point3i pos)=0;
  virtual void pop(void)=0;
  virtual double getData(unsigned int i)=0;
  virtual  void set(unsigned int i, double d, Point3i addr)=0;
  virtual void insert(double d, const cv::Point3i pos, int insertPos)=0;
  virtual 	int compareInsert(double d, cv::Point3i pos)=0;
  virtual vector<cv::Point3i> GetAddress(void)=0;
  virtual cv::Point3i GetAddress(unsigned int i)=0;
  virtual CandidateRecord getMin(void)=0;
  virtual int getMinIdx(void)=0;
};
class LinkArray:public LinkStruct
{
private:
  list<CandidateRecord> rec;
public:
  LinkArray(){};
  ~LinkArray(){};
  int getLength()
  {
    return rec.size();
  }
  void push(double d,const cv::Point3i pos)
	{
    CandidateRecord temp;
    temp.addr = pos;
    temp.data = d;
    rec.push_back(temp);
	}
  void pop(void)
	{
    std::lock_guard<std::mutex> guard(mymutex);
    rec.pop_back();
	}
  CandidateRecord getMin(void)
  {
    return *(--rec.end());
  }
   int getMinIdx(void)
  {
    return rec.size()-1;

   }
  double getData(unsigned int i)
	{
    std::lock_guard<std::mutex> guard(mymutex);
 		CV_DbgAssert(i<rec.size());
    list<CandidateRecord>::iterator iter=rec.begin();
    for (unsigned int ii=0; ii< i; ii++)
      iter++;
		return iter->data; 
	}

	void set(unsigned int i, double d, Point3i addr)
	{
		CV_DbgAssert(i<rec.size());
    list<CandidateRecord>::iterator iter=rec.begin();
    for (unsigned int ii=0; ii< i; ii++)
      iter++;
		iter->data=d;
		iter->addr=addr;
	}

	void insert(double d, const cv::Point3i pos, int insertPos)
	{
    CV_DbgAssert(insertPos<=rec.size());//allow for the last position
    list<CandidateRecord>::iterator iter=rec.begin();
    for ( int ii=0; ii< insertPos; ii++)
      iter++;
    CandidateRecord temp;
    temp.addr = pos;
    temp.data = d;
    rec.insert(iter,temp);
	}
  
	int compareInsert(double d, cv::Point3i pos)
	{
    std::lock_guard<std::mutex> guard(mymutex);
    list<CandidateRecord>::iterator iter=rec.begin();
    CandidateRecord temp;
    temp.addr = pos;
    temp.data = d;
		int count=0;
		for (iter=rec.begin();iter!=rec.end();iter++)
		{
			if (iter->data<d)
			{
        rec.insert(iter,temp);
				return count;
			}
			count++;
		}
		rec.push_back(temp);
		return count;
	}

	vector<cv::Point3i> GetAddress(void)
	{
    vector<cv::Point3i> address;
    list<CandidateRecord>::iterator iter=rec.begin();
    for (iter=rec.begin();iter!=rec.end();iter++)
		{
      address.push_back(iter->addr);
    }
		return address;
	}
  cv::Point3i GetAddress(unsigned int i)
	{
   	CV_DbgAssert(i<rec.size());
    list<CandidateRecord>::iterator iter=rec.begin();
    for (unsigned int ii=0; ii< i; ii++)
      iter++;
		return iter->addr; 
	}

};
class LinkQueue:public LinkStruct
{
private:
  //vector<CandidateRecord> rec;
	vector<double> data;
	vector<cv::Point3i> address; 
	unsigned int length;
//	struct node
//	{
//		double data;
//		int x,y,z; //address
//		node *link;
//	}*p;
//
//
public:
//
	LinkQueue(void)
	{
//		p=NULL;
		length = 0;
	}
//
	LinkQueue(double d)
	{
//		p=new node;
//		p->data=d;
//		p->x=0;
//		p->y=0;
//		p->z=0;
//		p->link = NULL;
		data.push_back(d);
		address.push_back(cv::Point3i(-1,-1,-1));
		length = 1;
	}
	LinkQueue(double d,const cv::Point3i pos)
	{
		data.push_back(d);
		address.push_back(pos);
		length =1;
	}
//	LinkQueue(double d, int a, int b, int c)
//	{
//		p=new node;
//		p->data=d;
//		p->x=a;
//		p->y=b;
//		p->z=c;
//		p->link = NULL;
//	}

	LinkQueue(int num)
	{
		for (int i=0;i<num;i++)
		{
			data.push_back(DBL_MAX);
			address.push_back(cv::Point3i(-1,-1,-1));
		}
		length = num;
	}
//	~LinkQueue(void)
//	{
//	}
//
//
	void setData(unsigned int i, double d)
	{
		data[i]=d;
	}

	void SetMaxLength(int num)
	{
		length = num;
	}

	void push(double d,const cv::Point3i pos)
	{

		data.push_back(d);
		address.push_back(pos);
		if (data.size()>length) //queue full, note length is fixed
		{
			data.erase(data.begin());
			address.erase(address.begin());
		}
	}
	void pop(void)
	{
    std::lock_guard<std::mutex> guard(mymutex);
		length--;
		data.erase(data.begin());
		address.erase(address.begin());
	}
	int getLength(void)
	{
		return length;
	}
	double getData(unsigned int i)
	{
    std::lock_guard<std::mutex> guard(mymutex);
		return data[i]; 
	}
	void set(unsigned int i, double d, Point3i addr)
	{
		CV_DbgAssert(i<length);
		data[i]=d;
		address[i]=addr;
	}

	void insert(double d, const cv::Point3i pos, int insertPos)
	{

		data.insert(data.begin()+insertPos,d);
		address.insert(address.begin()+insertPos,pos);
		if (data.size()>length) //queue full
		{
			data.erase(data.begin());
			address.erase(address.begin());
		}
	}

	int compareInsert(double d, cv::Point3i pos)
	{
    std::lock_guard<std::mutex> guard(mymutex);
		vector<double>::iterator iter=data.begin();
		int count=0;
		for (iter=data.begin();iter!=data.end();iter++)
		{
			if ((*iter)<d)
			{
				if (count==0)
				{
					return count;//d greater than any one in the queue, pass
				}
				data.insert(iter,d);
				address.insert(address.begin()+count,pos);
				if (data.size()>length) //queue full
				{
					data.erase(data.begin());
					address.erase(address.begin());
				}
				return count;
			}
			count++;
		}
		push(d,pos);
		return count;
	}

	vector<cv::Point3i> GetAddress(void)
	{
		return address;
	}

  cv::Point3i GetAddress(unsigned int i)
  {
   return address[i];
  }

  CandidateRecord getMin(void)
  {
    CandidateRecord rec;
    rec.addr = address[length-1];
    rec.data = data[length-1];
    return rec;
  }
  int getMinIdx(void)
  {
    return length-1;

   }
//	void append(double d, int a, int b, int c)
//	{
//		node *q,*t;
//		if( p == NULL )
//		{
//			p = new node;
//			p->data = d;
//			p->x=a;
//			p->y=b;
//			p->z=c;
//			p->link = NULL;
//		}
//		else
//		{
//			q = p;
//			while( q->link != NULL )
//				q = q->link;
//			t = new node;
//			t->data = d;
//			t->x =a;
//			t->y =b;
//			t->z =c;
//			t->link = NULL;
//			q->link = t;
//		}
// 
//	}
//	#
//	void add_as_first(double d)
//	{
//     node *q;
//	 q = new node;
//	 q->data = d;
//	 q->link = p;
//	 p = q;
//	}
};

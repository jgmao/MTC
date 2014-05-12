#include "TensorLite.h"
//20130815 decoule Tensor to high level operations
//#include "Steerable.h"
//#include "LRI.h"
//#include "Metric.h"
namespace tensor{
//20130815 template<class T, size_t cn> Mat Tensor<T,cn>::stsim2_lse_weight = Mat(); 
/////////////2.1 constructor, copy  assignment and type conversion ////////////////////////
//BufferGPU gbuf;
template<class T, size_t cn> Tensor<T,cn>::Tensor(void):Mat()
{
	this->tsOffset = Point3i();
	this->cFileName = "unknown";
	this->tsSize = Size3();
	this->mxFrame = Mat();
	//subWinSize = Size3();
	//subWinStep = Size3();
	this->debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(int height, int width, int depth, typename Tensor<T,cn>::c_ref_type val): Mat()
{
    *this = Mat(height,width,CV_MAKETYPE(DataType<T>::depth,cn),val);
	this->tsSize = Size3(height,width,depth);
	this->cFileName = "unknown";
	this->tsOffset = Point3i();
	this->mxFrame = Mat(height,width,CV_MAKETYPE(DataType<T>::depth,cn));
	this->debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(const Size3& size, const Vec<T,cn>& val): Mat()
{
	*this = Mat(size, CV_MAKETYPE(DataType<T>::depth,cn),val);
	this->tsSize = size;
	this->cFileName = "unknown";
	this->tsOffset = Point3i();
	this->mxFrame = Mat(size.height,size.width,CV_MAKETYPE(DataType<T>::depth,cn));
	this->debugtrigger=false;
}
//template<class T, size_t cn> Tensor<T,cn>::Tensor(const Size3& size):Mat()
//{
//	*this = Mat(size, CV_MAKETYPE(DataType<T>::depth,cn));
//	this->tsSize = size;
//	this->cFileName = "unknown";
//	this->tsOffset = Point3i();
//	this->mxFrame = Mat(size.height,size.width,CV_MAKETYPE(DataType<T>::depth,cn));
//	//subWinSize = Size3();
//	//subWinStep = Size3();
//	this->debugtrigger=false;
//}
template<class T, size_t cn> Tensor<T,cn>::Tensor(const Tensor<T,cn>& ts): Mat(ts)
{
	this->tsSize = ts.tsSize;
	if ( ts.type() != CV_MAKETYPE(DataType<T>::depth,cn))
	{
		Tensor<T,cn> rst(tsSize);
		this->convertTo(rst,DataType<T>::depth);
		*this = rst;
	}
	this->cFileName = ts.cFileName;
	this->tsOffset = ts.tsOffset;
	this->mxFrame = ts.mxFrame;
	this->debugtrigger=ts.debugtrigger;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(const Mat& mt): Mat(mt)
{
	//include the type conversion function
	CV_Assert(mt.dims==2);
  this->tsSize = Size3(mt.rows,mt.cols,1);
        this->SetFrame(0,mt);
	this->cFileName = "unknown";
	this->tsOffset = Point3i();
	if (mt.type() != CV_MAKETYPE(DataType<T>::depth,cn))
	{
		Tensor<T,cn> rst(this->tsSize);
		this->convertTo(rst,DataType<T>::depth);
		*this = rst;
	}
	this->debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(const string cFileName)//:subTensor(Tensor_<T,cn>())
{

	this->Load(cFileName);
	this->tsOffset = Point3i();

}
template<class T, size_t cn> Tensor<T,cn>::~Tensor(void)
{
}

template<class T, size_t cn> Tensor<T,2> Tensor<T,cn>::ToComplex(void) const
{
	CV_Assert( 2 == 2*cn || 2 ==cn);
	if(cn == 2)
		return *this;
	Tensor<T,2> rst(size());
	rst.SetFileName(this->cFileName);
	for (int i=0; i< size().depth; i++)
	{
		rst.SetFrame(i,mylib::toComplex(this->operator[](i)));
	}
	rst.debugtrigger = this->debugtrigger;
	return rst;
}

template<class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::operator= (const Tensor<T,cn>& ts)
{
	this->Mat::operator=(ts);
	this->tsSize = ts.tsSize;
	this->cFileName = ts.cFileName;
	this->tsOffset = ts.tsOffset;
	this->debugtrigger=ts.debugtrigger;
	return *this;
}

////////////2.2 public utilities	/////////////////////////////////////////////////////

template<class T, size_t cn> void Tensor<T,cn>::Display(int flag) const
{
	string tempName;
	CV_Assert(cn == this->channels());
	Tensor<uchar,cn> tempTs;
	if (this->depth() != DataType<uchar>::depth)
		tempTs = Tensor<uchar,cn>(*this); //directly type conversion
	else
		tempTs = *this;
      	
	if (cFileName.empty())
		tempName="Unknown";
#if (CV_MINOR_VERSION < 5)
	cv::namedWindow(tempName,flag|CV_WINDOW_FREERATIO/*|CV_GUI_EXPANDED*/);
#else
	cv::namedWindow(tempName,flag|WINDOW_FREERATIO/*|CV_GUI_EXPANDED*/);//if Qt enabled, you can uncomment this
#endif
	cv::moveWindow(tempName.c_str(),100,100);
	
	if (tsSize.depth==1)
	{
		cv::imshow(tempName,tempTs.GetFrame(0));
		cv::waitKey();
		cv::destroyWindow(tempName.c_str());
	}
	else
	{
		for (int i=0; i< tsSize.depth; i++)
		{
			cv::imshow(tempName, tempTs.GetFrame(i));
			if(i == tsSize.depth -1)
				i = 0;
			if (cv::waitKey(30)>=0) break;
		}
		cv::destroyWindow(tempName.c_str());
	}

}

template<class T, size_t cn> void Tensor<T,cn>::Display(int sec, int flag) const
{
	string tempName;
	CV_Assert(cn == this->channels());

	Tensor<uchar,cn> tempTs;
	if (this->depth() != DataType<uchar>::depth)
		tempTs = Tensor<uchar,cn>(*this); //directly type conversion
	else
		tempTs = *this;
	

	if (cFileName.empty())
		tempName="Unknown";
#if (CV_MINOR_VERSION > 5)
	cv::namedWindow(tempName,flag|WINDOW_KEEPRATIO/*|CV_GUI_EXPANDED*/);//if Qt enabled, you can uncomment this
#else
	cv::namedWindow(tempName,flag|CV_WINDOW_KEEPRATIO/*|CV_GUI_EXPANDED*/);//if Qt enabled, you can uncomment this
#endif
	cv::moveWindow(tempName.c_str(),100,100);

	if (tsSize.depth==1)
	{
		cv::imshow(tempName,tempTs.GetFrame(0));
		cv::waitKey(sec);
		cv::destroyWindow(tempName.c_str());
	}
	else
	{
		for (int i=0; i< tsSize.depth; i++)
		{
			cv::imshow(tempName, tempTs.GetFrame(i));
			if(i == tsSize.depth -1)
				i = 0;
			if (cv::waitKey(sec)>=0) break;
		}
		cv::destroyWindow(tempName.c_str());
	}

}
//!20130815 move this to utility module inside Tensor (not mylib)

template<typename T, size_t cn> void Tensor<T,cn>::DisplayAll(std::vector<Tensor<T,cn>>& vts, int row, int col, bool save, string savename) const
{
  Size3 padsz(0,0,0);
  CV_Assert(vts.size()<=row*col);
  for (Tensor<T,cn>& ts :vts)
  {
    if (padsz.height<ts.size().height)
      padsz.height=ts.size().height;
    if (padsz.width<ts.size().width)
      padsz.width=ts.size().width;
  }
  padsz.height+=padsz.height/4;
  padsz.width+=padsz.width/4;
  Tensor<T,cn> all(row*padsz.height,col*padsz.width,1,255);
  int count=0;
  int i=0,j=0;
  for (Tensor<T,cn>& ts :vts)
  {
    i=count/col; j=count%col;
    all.SetBlock(Point3i(i*padsz.height,j*padsz.width,0),ts);
    count++;
  }
  if (save)
    all.SaveBlock(savename);
  else
    all.Display();
}


template<class T, size_t cn> void Tensor<T,cn>::Load(string cFileName) 
{
	this->cFileName = cFileName;
	cv::VideoCapture cap;
	cap.open(cFileName);
	cap>>this->mxFrame;//read next frame
	cv::Mat grayone;
	if (!this->mxFrame.data) //if cannot be read as a stream
	{
		cout<<"cannot be read by stream, use imread()"<<endl;
		if (cn == 3)
			this->mxFrame = cv::imread(cFileName,1);
		else if (cn == 1)
			this->mxFrame = cv::imread(cFileName,0);
		else if (cn==2) //complex
		{
			this->mxFrame = cv::imread(cFileName,0);
			Mat azMat = Mat::zeros(this->mxFrame.size(),this->mxFrame.type());
			vector<Mat> vcm;
			vcm.push_back(this->mxFrame);
			vcm.push_back(azMat);
			cv::merge(vcm,this->mxFrame);
		}
		else
			this->mxFrame = cv::imread(cFileName,-1);
		CV_Assert(!this->data && this->mxFrame.channels() == cn ); 
		const int sz[] = {this->mxFrame.rows,this->mxFrame.cols};
		this->mxFrame.convertTo(this->mxFrame,DataType<T>::depth);//new added
		this->create(2,sz,this->mxFrame.type());
		this->tsSize = Size3(sz[0],sz[1],1);
		SetFrame(0,this->mxFrame.clone());
		cout<<"size: "; 
		this->tsSize.Print();
		cout<<"type: "<<this->type();
		cout<<", channel: "<<this->channels()<<endl;
	}
	else
	{
		std::cout<<"read as stream\n";
		//cout<<mxFrame.channels()<<","<<mxFrame.type()<<endl;
		if (cn==1&&this->mxFrame.channels()==3)
		{
			cv::cvtColor(this->mxFrame,grayone,COLOR_RGB2GRAY);
			this->mxFrame = grayone;
		//	cout<<mxFrame.type();
		}
		if (cn==2) //complex
		{
			cv::cvtColor(this->mxFrame,grayone,COLOR_RGB2GRAY);
			Mat azMat = Mat::zeros(grayone.size(),grayone.type());
			vector<Mat> vcm;
			vcm.push_back(grayone);
			vcm.push_back(azMat);
			cv::merge(vcm,this->mxFrame);
		}
		int sz[] = {this->mxFrame.rows,this->mxFrame.cols};
		cout<<this->mxFrame.size()<<endl;
    this->mxFrame.convertTo(this->mxFrame,DataType<T>::depth);//new added
    		//cout<<mxFrame.channels()<<","<<mxFrame.type()<<endl;
                this->create(2, sz,this->mxFrame.type());
                this->tsSize = Size3(sz[0],sz[1],1);
		this->SetFrame(0,this->mxFrame.clone());
	  //never treat as a video in TensorLite
		//cap>>this->mxFrame;
		//while(mxFrame.data)
		//{	
		//	if (cn==1&&mxFrame.channels()>1)
		//	{
		//		cv::cvtColor(mxFrame,grayone,COLOR_RGB2GRAY);
		//		mxFrame = grayone;	
		//	}
  //    Mat tempFrame(3,sz,mxFrame.type(),mxFrame.data);
		//	this->push_back(tempFrame.clone()); //untested
		//	cap>>mxFrame;
		//	tsSize.depth++; //measure the length of video
		//}
		*this = Tensor<T,cn>(*this);//new added
		cout<<"read a fream"<<endl;
	}
	cap.release();
	this->debugtrigger=false;
	//this->Display();
}

//////////////////////// 2.3 accessors
template<class T, size_t cn> Size3 Tensor<T,cn>::size(void) const 
{
	return this->tsSize;
}

//! data accesor will return the reference to the pixel
template<class T, size_t cn> typename Tensor<T,cn>::ref_type Tensor<T,cn>::operator() (int x, int y, int z)
{
	CV_Assert( CV_MAKETYPE(DataType<T>::depth,cn) == this->type());
	return this->at<value_type>(x,y);

}

template<class T, size_t cn> typename Tensor<T,cn>::c_ref_type Tensor<T,cn>::operator() (int x, int y, int z) const
{
	//int t = this->type();
	//int t2 = CV_MAKETYPE(DataType<T>::depth,cn) ;
	CV_Assert( CV_MAKETYPE(DataType<T>::depth,cn) == this->type());

	return this->at<value_type>(x,y);

}

template<class T, size_t cn> typename Tensor<T,cn>::ref_type Tensor<T,cn>::operator[] (const Point3i& pos)
{
	return this->operator()(pos.x,pos.y,pos.z);
}
template<class T, size_t cn> typename Tensor<T,cn>::c_ref_type Tensor<T,cn>::operator[](const Point3i& pos) const
{
	return this->operator()(pos.x,pos.y,pos.z);
}

//! 20130820 I don't know if GetFrame should return a reference or a clone
/*! just set it as a Clone now, and GetFrameRef will give a reference
	  operator[] will return a reference too
*/
template< class T, size_t cn> Mat Tensor<T,cn>::GetFrame(int i)
{

	if (dims <=2 )
		return this->clone();//note this is a clone
	else
	{
		CV_Assert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
    
		//mxFrame = this->row(i);
    		this->mxFrame =  Mat(*this, Range(i, i+1), Range::all());
    //mylib::DisplayMat(mxFrame);
		Mat rst = Mat(2,this->mxFrame.size.p+1,this->mxFrame.type(),this->mxFrame.data,this->mxFrame.step.p+1);
    //if (rst.size().height<10)
    //mylib::DisplayMat(rst);
		return rst;
	}
}

template< class T, size_t cn> Mat Tensor<T,cn>::GetFrameRef(int i)
{
	if (dims <=2 )
		return *this;
	else
	{
		CV_Assert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
		//mxFrame = this->row(i);
     //Mat* arrays[] = {this};
    		Mat* arrays[] = {(Mat*)this};
    		Mat planes[1];
    		NAryMatIterator it((const Mat**)arrays, planes,1);
    		this->mxFrame = it.planes[i].reshape(cn,this->size().height);
    //if (mxFrame.size().height<0)
    //mylib::DisplayMat(mxFrame);
    		return this->mxFrame;
	}
}

template< class T, size_t cn> const Mat Tensor<T,cn>::GetFrameRef(int i) const
{
	if (dims <=2 )
		return *this;
	else
	{

		CV_Assert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
		//mxFrame = this->row(i);
     //Mat* arrays[] = {this};
                Mat* arrays[] = {(Mat*)this};
                Mat planes[1];
                NAryMatIterator it((const Mat**)arrays, planes,1);
                Mat rst = it.planes[i].reshape(cn,this->size().height);

    //if (mxFrame.size().height<0)
    //mylib::DisplayMat(mxFrame);
                return rst;
        }
}

template< class T, size_t cn> void Tensor<T,cn>::GetFrameRef(int i, Mat& rst) const 
{
	if (dims <=2 )
		rst= *this;
	else
	{
		CV_Assert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
		//mxFrame = this->row(i);
     //Mat* arrays[] = {this};
    		Mat* arrays[] = {(Mat*)this};
   	 	Mat planes[1];
    		NAryMatIterator it((const Mat**)arrays, planes,1);
    		rst = it.planes[i].reshape(cn,this->size().height);
	}
}

template< class T, size_t cn> const Mat Tensor<T,cn>::GetFrame(int i) const
{
	if (dims <=2 )
		return this->clone();
	else
	{
		CV_Assert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
    		Mat tempMat = Mat(*this, Range(i, i+1), Range::all());
	  //Mat tempMat = this->row(i);
		Mat rst = Mat(2, tempMat.size.p+1,tempMat.type(),tempMat.data,tempMat.step.p+1);
		return rst;
	}
}

template< class T, size_t cn> Mat Tensor<T,cn>::operator[](int i)
{
	return this->GetFrameRef(i);
}

template< class T, size_t cn> const Mat Tensor<T,cn>::operator[](int i) const
{
	return this->GetFrameRef(i);
}

template< class T, size_t cn> Point3i Tensor<T,cn>::offset() const
{
	return this->tsOffset;
}

template< class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Clone(void) const
{
  	std::lock_guard<std::mutex> lock(this->tsmutex);
  	Tensor<T,cn> rst(this->Mat::clone());
                rst.SetFileName(this->cFileName);
                rst.SetOffset(this->tsOffset);

	return rst;
}
//! SetFrame is a copyTo opertaion, if you don't need a copy, don't use this.
template< class T, size_t cn> void Tensor<T,cn>::SetFrame(int i, const Mat& frm)
{
	//use to modify(or assignment by copy) one frame
	if (dims <=2 )
		frm.copyTo(*this);
	else
	{
		CV_Assert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
		this->mxFrame = GetFrame(i); //get a reference to the frame
		frm.copyTo(this->mxFrame);
		//make sure frm is nonempty
	//CV_Assert(frm.total()!=0&&frm.isContinuous() && mxFrame.isContinuous());
        //size_t planeSize = frm.elemSize() * frm.rows * frm.cols;
		//memcpy(mxFrame.data,frm.data,planeSize);
	}
}


template<class T, size_t cn> void Tensor<T,cn>::Print(const string& fname, bool tofileonly) const
{

	fstream logfile;
	logfile.open("./temp/"+fname+".txt",ios::out);
	//bool tofileonly=false;
	if (!tofileonly)
	{
		cout<<"File: "<<fname<<" Display "<<this-> size().height <<" "<<this-> size().width<<" Data type: "<<this-> type()<<endl;
		cout<<"cn="<<cn<<endl;
	}
	//logfile<<"File: "<<cFileName<<"Display "<<size().height <<" "<<size().width<<" Data type: "<<type()<<endl;
	
	for (int d=0; d< this-> size().depth; d++) 
	{
		if (!tofileonly)	
			cout<<"frame "<<d<<" :"<<endl;
		//logfile<<"frame "<<d<<" :"<<endl;
		for (int c=0; c< cn;c++)
		{
			if (!tofileonly)
				cout<<"channel "<<c<<" :"<<endl;
			//logfile<<"channel "<<c<<" :"<<endl;
			for (int x = 0; x < this-> size().height; x++)
			{
				for (int y = 0; y < this-> size().width; y++)
				{
					if (!tofileonly)
						cout<< double(this->operator()(x,y,d)[c]) <<"\t";
					logfile<< double(this->operator()(x,y,d)[c])<<"\t";
				}
				if (!tofileonly)
					cout<<";"<<endl;
				//logfile<<";"<<endl;
        			logfile<<endl;
			}
			if (!tofileonly)
				cout<<"------------end of channel "<<c<<" -----------------"<<endl;
			//logfile<<";------------end of channel "<<c<<" -----------------"<<endl;

		}
		if (!tofileonly)	
			cout<<"------------end of frame "<<d<<" -----------------"<<endl;
		//logfile<<";------------end of frame "<<d<<" -----------------"<<endl;
	}
	logfile.close();
}



template <class T, size_t cn> void Tensor<T,cn>::SetFileName(const string& name)
{
	this->cFileName = name;
}

template <class T, size_t cn> void Tensor<T,cn>::SetOffset(const Point3i& pos)
{
	this->tsOffset = pos;
}

template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::GetBlockRef(const Cube& roi)
{
	AssertRange(roi.offset(),roi.size());
	Rect r = roi.toRect();
	Tensor<T,cn> rst = this->Mat::operator()(r);
	rst.SetFileName(this->cFileName);
	rst.SetOffset(roi.offset()+this->offset());
  return rst;

}

template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::GetBlockRef(const Cube& roi) const
{
	AssertRange(roi.offset(),roi.size());
	Rect r = roi.toRect();
	Tensor<T,cn> rst = this->Mat::operator()(r);
	rst.SetFileName(this->cFileName);
	rst.SetOffset(roi.offset()+this->offset());
  return rst;
}



//! 20130820 I need to test this
//! get block get a copy of subblock
//! in openCV Rect::x is the ----> axis
//! and Rect::y is | axis
//!                v
//! the structure is Rect::Rect(--> , | , width, height)
//!                                   v
//! so in TensorLib, the convert from Cube to Rect is
//! Rect(Cube::y, Cube::x, Cube::width, Cube::height)
//! I made a converter Cube::toRect() to do it.
template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::GetBlock(const Cube& roi)
{
	Rect r = roi.toRect();
	Tensor<T,cn> rst(roi.size());
	rst.SetFrame(0,this->Mat::operator()(r));
	rst.SetFileName(this->cFileName);
	rst.SetOffset(roi.offset());
	return rst;
}

template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::GetBlock(const Cube& roi) const
{
	Rect r = roi.toRect();
	Tensor<T,cn> rst(roi.size());
	rst.SetFrame(0,this->Mat::operator()(r));
	rst.SetFileName(this->cFileName);
	rst.SetOffset(roi.offset()+this->offset());
	return rst;
}

template<class T, size_t cn> void Tensor<T,cn>::Ref(const Cube& roi, Tensor<T,cn>& dst) const
{
  Rect r = roi.toRect();
	dst.Mat::operator=(this->Mat::operator()(r));
	dst.tsSize = roi.size();
	dst.tsOffset = roi.offset();
	dst.SetFileName(this->cFileName);
	dst.SetOffset(roi.offset()+this->offset());
}

template <class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::SetBlock(const Point3i& pos, const Tensor<T,cn>& ts)
{
	CV_Assert(ts.type() == type());
	AssertRange(pos,ts.size());
	Cube roi(pos,ts.size());
	cv::Mat tempMat = GetBlockRef(roi);
	ts.copyTo(tempMat);
	return *this;
}
template <class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::SetBlock(const Tensor<T,cn>& ts)
{
	return SetBlock(Point3i(0,0,0),ts);
}
template<class T, size_t cn> void Tensor<T,cn>::AssertRange(const Point3i& pos, const Size3& sz) const
{
	CV_Assert( pos.x + sz.height <= this-> tsSize.height &&
				  pos.y + sz.width  <= this-> tsSize.width &&
				  pos.z + sz.depth  <= this-> tsSize.depth);
}

template<class T, size_t cn> void Tensor<T,cn>::SaveBlock(const string& cFileName, bool isGray)
{
	if (size().depth>1)
	{
#if (CV_MINOR_VERSION > 5)
		cv::VideoWriter vd(cFileName+".avi",cv::VideoWriter::fourcc('D','I','B',' '),30, cv::Size(size().width,this-> size().height),true);
#else
	    cv::VideoWriter vd(cFileName+".avi", CV_FOURCC('D','I','B',' '),30, cv::Size(size().width,this-> size().height),true);
#endif
		for (int i=0;i<this-> size().depth;i++)
		{

			if (isGray)
			{
				Mat grayone = this-> GetFrame(i);
				cv::cvtColor(grayone,grayone,COLOR_RGB2GRAY);
				vd<<grayone;
			}
			else
				vd<<this-> GetFrame(i);

		}
	}
	else
	{
		if (isGray)
		{
			Mat grayone = this-> GetFrameRef(0);
			cv::cvtColor(grayone,grayone,COLOR_RGB2GRAY);
			cv::imwrite(cFileName,grayone);
		}
		else
			cv::imwrite(cFileName,this-> GetFrameRef(0));
	}
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Crop(const Point3i& pos, const Size3& sz) const
{
  //cout<<"pos"<<pos<<"sz"<<sz<<endl;
  return this-> GetBlockRef(Cube(pos,sz)).Clone();
}
/* 20130815
template<class T, size_t cn> void Tensor<T,cn>::SetSubWinSize(const Size3& sz)
{
	subWinSize = sz;
	//make the value at least 1
	if (subWinSize.height == 0)
		subWinSize.height = 1;
	if (subWinSize.width == 0)
		subWinSize.width = 1;
	if (subWinSize.depth == 0)
		subWinSize.depth = 1;

}

template<class T, size_t cn> void Tensor<T,cn>::SetSubWinStep(const Size3& sz)
{
	subWinStep  = sz;
	//make the value at least 1
	if (subWinStep.height == 0)
		subWinStep.height = 1;
	if (subWinStep.width == 0)
		subWinStep.width = 1;
	if (subWinStep.depth == 0)
		subWinStep.depth = 1;
}
template<class T, size_t cn> vector<Vec<T,cn>> Tensor<T,cn>::GetTagLighting(void) const
{
		return lightTag;
}

template<class T, size_t cn> vector<Vec<T,cn>> Tensor<T,cn>::GetCanLighting(void) const
{
		return lightCan;
}
template<class T, size_t cn> Size3 Tensor<T,cn>::GetSubWinSize(void) const
{
		return subWinSize;
}
template<class T, size_t cn> Size3 Tensor<T,cn>::GetSubWinStep(void)const
{
		return subWinStep;
}
*/
template<class T, size_t cn> bool Tensor<T,cn>::IsInside(const Point3i& pos) const
{
	if ( pos.x >= tsOffset.x && pos.x < tsOffset.x + this-> tsSize.height &&
		pos.y >= tsOffset.y && pos.y < tsOffset.y + this-> tsSize.width  &&
		pos.z >= tsOffset.z && pos.z < tsOffset.z + this-> tsSize.depth )
		return true;
	else
		return 0;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Row(int i)
{
	Tensor<T,cn> rst(Size3(1,this-> tsSize.width,this-> tsSize.depth));
	if (this-> tsSize.depth == 1)
		rst = this->row(i);
	else
	{
		Range r[] = {Range::all(),Range(i,i+1),Range::all()};
		rst = this->Mat::operator()(r);
	}
	//this->Ref(Cube(i,0,0,1,tsSize.width,tsSize.depth), *rst);
	return rst;
}

//! This and GetBlock is different to GetFrame, since it returns a Tensor not Mat
//! so when I talk about return a reference, it means the Mat data is reference, but
//! the value of tsSize, offset are all copy
//! operator(Cube) return a reference
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator()(const Cube& roi)
{	
	return GetBlockRef(roi);
}

////////2.4 arithmatic operations //////////////////////////////////////////

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator+(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) + Tensor<T,cn>(this->size(),s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator-(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) - Tensor<T,cn>(this->size(),s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator*(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) * Tensor<T,cn>(size().height,size().width, size().depth,s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator/(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) / Tensor<T,cn>(size().height,size().width, size().depth,s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator*(const Tensor<T,cn>& ts) const
{

	Tensor<T,cn> rst = Clone();
	//rst.Print("clone");
	for (int z =0; z < size().depth; z++)
	{
		if (cn==2 || cn == 6)
			cv::mulSpectrums(rst[z],ts[z],rst[z],DFT_ROWS,false); // use CCS strucrue for complex mul
		else
    		{

                        rst.SetFrame(z,GetFrameRef(z).mul(ts[z]));
    		}
	}
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator/(const Tensor<T,cn>& ts) const
{
	Tensor<T,cn> rst = Clone();
	if (cn==2|| cn==6) //complex div
	{
		for (int i=0; i<tsSize.depth; i++)
		{
			Mat tempMat;
			//mylib::DisplayMat(rst.GetFrame(i));
			cv::mulSpectrums(ts.GetFrameRef(i),ts.GetFrameRef(i),tempMat,DFT_ROWS,true);
			//mylib::DisplayMat(tempMat);
			vector<Mat> temp;
			cv::split(tempMat, temp);
			temp[1] = temp[0];
			merge(temp,tempMat);
			//mylib::DisplayMat(tempMat);
			cv::mulSpectrums(rst.GetFrameRef(i),ts.GetFrameRef(i),rst.GetFrameRef(i),DFT_ROWS,true);
			//mylib::DisplayMat(rst.GetFrame(i));
			//Vec<T,cn> a= rst(0,0,0);
			//Vec<T,cn> b = tempMat.at<Vec<T,cn>> (0,0);
			rst.GetFrameRef(i) /= tempMat;
			//mylib::DisplayMat(rst.GetFrame(i));
		}
	}
	else
	{
		for (int i = 0; i< tsSize.depth; i++)
			rst.GetFrameRef(i) /= ts.GetFrameRef(i);

	}
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator+(const Tensor<T,cn>& ts) const
{
	Tensor<T,cn> rst = Clone();
	rst += ts;
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator-(const Tensor<T,cn>& ts) const
{
	Tensor<T,cn> rst =Clone(); // copy all parameters
	rst -= ts;
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator+(T s) const
{
	return (*this)+Tensor<T,cn>(this->size(),s);
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator-(T s) const
{
	return (*this)-Tensor<T,cn>(this->size(),s);
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator*(T s) const
{
	return (*this)*Tensor<T,cn>(this->size(),s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator/(T s) const
{
	return (*this)/Tensor<T,cn>(this->size(),s);
}

template<class T, size_t cn> typename Tensor<T,cn>::value_type Tensor<T,cn>::All(T val) const
{
	return value_type::all(val);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Abs(void) const
{
	
	Tensor<T,cn> rst(this->size());
	if (cn == 2) // if data type is not float/double this will throw type exception
	{
		for (int i=0; i< size().depth; i++)
		{
			cv::mulSpectrums(this->GetFrameRef(i),this->GetFrameRef(i),rst.GetFrameRef(i),DFT_ROWS,true);
			cv::pow(rst.GetFrameRef(i),0.5,rst.GetFrameRef(i));
		}
	}
	else
	{
		for (int i=0; i< size().depth; i++)
			rst.SetFrame(i,cv::abs(this->GetFrameRef(i)));
	}
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Square(void) const
{

  Tensor<T,cn> rst(this->size());
  if (cn == 2) // if data type is not float/double this will throw type exception
  {
     cv::mulSpectrums(this->GetFrameRef(0),this->GetFrameRef(0),rst.GetFrameRef(0),DFT_ROWS,true);
  }
  else
  {
     Mat temp =  this->GetFrameRef(0).mul(this->GetFrameRef(0));
     rst = temp;
  }
  return rst;
}

//template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::SquareGPU(BufferGPU& gbuf) const
//{

//  Tensor<T,cn> rst(this->size());
//  gbuf.gI1.upload(*this);
//  if (cn == 2) // if data type is not float/double this will throw type exception
//  {
//     gpu::mulSpectrums(gbuf.gI1,gbuf.gI1,gbuf.gs,DFT_ROWS,true,gbuf.stream);
//  }
//  else
//  {
//     gpu::multiply(gbuf.gI1,gbuf.gI1,gbuf.gs,1,-1,gbuf.stream);
//  }
//  gbuf.gs.download(rst,gbuf.stream);
//  return rst;
//}


template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::AbsDiff(c_ref_type s) const
{
//	if (DataType<T>::depth != CV_32F && DataType<T>::depth != CV_64F)
//		CV_Error_(CV_StsInplaceNotSupported, ("depth is %d, only support float:%d or doule:%d",DataType<T>::depth,DataType<float>::depth,DataType<double>::depth));
	Tensor<T,cn> rst(this->size());
	Scalar v(s);
	if (cn == 2) // complex
	{
		rst = (*this-s).Abs();
	}
	else
	{
		for (int i=0; i< size().depth; i++)
			cv::absdiff(this->GetFrameRef(i),v,rst.GetFrameRef(i));
	}
	return rst;
}


template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::AbsDiff(const Tensor& ts) const
{
//	if (DataType<T>::depth != CV_32F && DataType<T>::depth != CV_64F)
//		CV_Error_(CV_StsInplaceNotSupported, ("depth is %d, only support float:%d or doule:%d",DataType<T>::depth,DataType<float>::depth,DataType<double>::depth));
//
	Tensor<T,cn> rst(this->size());

	if (cn==2) //complex
	{
		rst = (*this - ts).Abs();
	}
	else
	{
		for (int i=0; i< size().depth; i++)
			cv::absdiff(this->GetFrameRef(i),ts.GetFrameRef(i),rst.GetFrameRef(i));
	}
	return rst;
}

template<class T, size_t cn> Tensor<T,1> Tensor<T,cn>::Real(void) const
{
	if (cn==2)
	{
		Tensor<T,1> rst(tsSize);
		vector<Mat> temp;
		cv::split(*this,temp);
		rst = temp[0];
		return rst;
	}
	else
		return this->Clone();
}
template<class T, size_t cn> Tensor<T,1> Tensor<T,cn>::Imag(void) const
{
	if (cn==2)
	{
		Tensor<T,1> rst(tsSize);
		vector<Mat> temp;
		cv::split(*this,temp);
		rst = temp[1];
		return rst;
	}
	else
		return Tensor<T,1>(this->size());
}


template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Conjugate(void) const
{
	Tensor<T,cn> rst(this->size());
	if (cn == 2) //complex
	{
		vector<Mat> temp;
		cv::split(*this,temp);
		temp[1] = -temp[1];
		cv::merge(temp,rst);
	}
	else
		rst = this->Clone();
	return rst;
}
//template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::HorDev(void) const
//{
//    Tensor<T,cn> rst = this->operator(Cube(0,0,0,this->tsSize.height,this->tsSize.width-1,1))-
//            this->operator(Cube(0,1,0,this->tsSize.height,this->tsSize.width-1,1));
//    return rst.Clone();
//}

//template <class T, size_t cn> Tensor<T<cn> Tensor<T,cn>::VerDev(void) const
//{
//Tensor<T,cn> rst = this->operator(Cube(0,0,0,this->tsSize.height-1,this->tsSize.width,1))-
//        this->operator(Cube(0,1,0,this->tsSize.height-1,this->tsSize.width-,1));
//return rst.Clone();

//}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Pow(double p) const
{
	Tensor<T,cn> rst(this->size());
	if (cn == 2)
	{
		int ipower = cvRound(p);
		if( fabs(ipower - p) < DBL_EPSILON )
		{
			if( ipower < 0 )
			{
				for(int i=0; i< tsSize.depth; i++)
					divide( 1., this->GetFrameRef(i), rst.GetFrameRef(i));
				if( ipower == -1 )
					return rst;
				ipower = -ipower;
			}
			else
				rst = this->Clone();

			switch( ipower )
			{
			case 0:
				rst = Tensor<T,2>(this-> tsSize.height,this-> tsSize.width,this-> tsSize.depth,Vec<T,2>(1,0));
				return rst;
			case 1:
				return rst;
			case 2:
				return rst*rst;
			default:
				CV_Error_(CV_StsNotImplemented,("power is %d, but only -2,-1,0,1,2 power of complex number are implemented",p));
			}
        }
		else
		{
			CV_Error_(CV_StsNotImplemented,("power is %d, but only -2,-1,0,1,2 power of complex number are implemented",p));
		}
   	}
	else
	{
		cv::pow(*this,p,rst);
	}
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Log(void) const
{

  Tensor<T,cn> rst(this->size());
  if (cn==2)
  {
   vector<Mat> temp;
   cv::split(*this,temp);
   cv::log(this->Abs(),rst);
   cv::phase(temp[1],temp[0],temp[1]);
   temp[0] = cv::Mat::zeros(temp[0].size(),temp[0].type()); 
   Mat temp2;
   cv::merge(temp,temp2);
   rst = rst+Tensor<T,cn>(temp2);
  }
  else
    cv::log(*this,rst);
  return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Exp(void) const
{
  Tensor<T,cn> rst(this->size());
/*temp unavailable 20130816 
 if (cn==2)
  {
    vector<Mat> temp;
    cv::split(*this, temp);
    Mat temp2;
    cv::exp(temp[0],temp2);
    temp[0] = cv::cos(temp[1]);
    temp[1] = cv::sin(temp[1]);
    temp[0] = temp[0].mul(temp2);
    temp[1] = temp[1].mul(temp2);
    cv::mege(temp,temp2);
    rst = Tensor<T,cn>(temp2);
  }
  else
*/
  cv::exp(*this,rst);
  return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Sqrt(void) const
{
	Tensor<T,cn> rst(this->size());
	cv::sqrt(*this,rst);
	return rst;
}

template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::Sum(void) const
{
	CV_Assert(cn <=4);
	cv::Scalar s = cv::sum(*this);
	Vec<T,cn> rst;
	//T rst=0;
	for (int i=0; i<cn; i++)
		rst[i] = saturate_cast<T>(s[i]);
	return rst;
}

template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::Mean(void) const
{
	Vec<T,cn> temp = this->Sum();
	T a = T(this->tsSize.area());
	return mylib::VecDiv<T,cn>(temp,a);
}

template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::Var(void) const
{
	Vec<T,cn> mu, var;
	//this->Print();
	cv::meanStdDev(*this,mu,var);
	return mylib::VecMul<T,cn>(var,var);
}

template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::LRMean(void) const
{
	Vec<T,cn> tempmean=0;
	//double tempmean=0;
	//CV_Assert(cn==1);
	int count=0;
		for (int ii=0;ii<this->size().height; ii++)
			for (int jj=size().height-ii-1;jj<size().width;jj++)
				for (int kk=0; kk<cn; kk++)
				{
				//cout<<this->operator()(ii,jj,0)[0];
				tempmean[kk]+=(this->operator()(ii,jj,0)[kk]);
				count++;
				}
		tempmean/=double(count);
	return tempmean;
}
template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::URMean(void) const
{
	Vec<T,cn> tempmean=0;
	//double tempmean=0;
	//CV_Assert(cn==1);
	int count=0;
		for (int ii=0;ii<this->size().height; ii++)
			for (int jj=0; jj <= ii; jj++)
				for (int kk=0; kk<cn; kk++)
				{
				//cout<<this->operator()(ii,jj,0)[0];
				tempmean[kk]+=(this->operator()(ii,jj,0)[kk]);
				count++;
				}
		tempmean/=double(count);
	return tempmean;
}

template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::PercentageMean(double low, double high) const
{
	cv::Vec<T,cn> mean=0;
	int total = this->size().volumn();
	cv::Mat vecMat(1, total,CV_MAKETYPE(DataType<T>::depth,cn));
//	T* pv = vecMat->ptr<T>();
	int count=0;
	//cv::Vec<T,cn> temp;
	for (int i=0; i< this->size().height; i++)
		for (int j=0; j< this->size().width;j++)
			for (int k=0; k<this->size().depth;k++)
			{
				vecMat.at<Vec<T,cn>>(0,count)= this->operator()(i,j,k);
				count++;
			}
	cv::Mat sortMat;
	cv::sort(vecMat,sortMat,CV_SORT_ASCENDING);
	//mylib::DisplayMat(vecMat);
	//mylib::DisplayMat(sortMat);
	int start = int(total*low);
	int end = int(total*high);
	if (low == high)//percentile
	{
		mean = sortMat.at<Vec<T,cn>>(0,start);
	}
	else
	{
		for (int i=start; i<end; i++)
		{
			//for (int j=0;j<cn;j++)
			mean += sortMat.at<Vec<T,cn>>(0,i);
		}
		mean/=(end-start);
	}
	return mean;

}
template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Transpose(void) const
{
	Tensor<T,cn> rst;
	if (tsSize.height==tsSize.width)
		rst = this->Clone();
	else
	{
		rst = Tensor<T,cn>(tsSize.width,tsSize.height,tsSize.depth);
		rst.tsOffset = tsOffset;
		//rst.subWinSize=subWinSize;
		//rst.subWinStep=subWinStep;
	}
	//	(Size3(this->size().width,this->size().height,this->size().depth));
	for (int i=0; i< size().depth; i++)
		cv::transpose(GetFrameRef(i),rst[i]);
	return rst;
}


template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::ExtendBoundary(Size3 extSz, typename Tensor<T,cn>::value_type val) const
{
	Tensor<T,cn> rst(size()+(extSz*2),val);
	rst.tsOffset = Point3i(extSz.height,extSz.width,extSz.depth);
  //cout<<"20130909 extsz :";extSz.Print();
  //cout<<"20130909 rstsz :";rst.size().Print();
	rst.SetBlock(rst.tsOffset,*this);
	//rst.SetSubWinSize(subWinSize);
	//rst.SetSubWinStep(subWinStep);
	return rst;
}

template <class T, size_t cn>
Tensor<T,cn> Tensor<T,cn>::ExtendHalfBoundary(Size3 extSz, typename Tensor<T,cn>::value_type val, bool which_side) const
{
	Tensor<T,cn> rst(size()+(extSz),val);
	if (which_side)
		rst.tsOffset = Point3i();
	else
		rst.tsOffset = Point3i(extSz.height,extSz.width,extSz.depth);
	rst.SetBlock(rst.tsOffset,*this);
	//rst.SetSubWinSize(subWinSize);
	//rst.SetSubWinStep(subWinStep);
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalMean(const Mat& ker,const Size3& subWinStep) const
{
//! 20130910 implement sliding window with steps
  Size3 sz(this->tsSize.height / subWinStep.height, this->tsSize.width/subWinStep.width, this->tsSize.depth/subWinStep.depth);
  Tensor<T,cn> rst(sz);
  //Cube roi(0,0,0,subWinStep.height, subWinStep.width, subWinStep.depth);
  Rect roi(0,0,ker.size().width,ker.size().height);
  vector<Mat> mats;
  cv::split(this->GetFrameRef(0),mats);
  for (unsigned int k=0; k< cn; k++)
  {
  for (unsigned int i=0;roi.y+roi.height<=this->tsSize.height;  ++i, roi.y+=subWinStep.height)
  {
    for (unsigned int j=0; roi.x+roi.width<=this->tsSize.width; ++j, roi.x+=subWinStep.width)
    {
        //mylib::DisplayMat(mats[k](roi));
        rst(i,j,0)[k]=(T)ker.dot(mats[k](roi));
    }
    roi.x=0;//reset to left
  }
  roi.y=0;
  }
  //auto rst2 = Filter2D(ker).SubSample(subWinStep);
  //cout<<(rst2-rst).Abs().Sum()[0]<<endl;//verify
  return rst;
}
//template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalMeanGPU(const Mat& ker,BufferGPU& gbuf,const Size3& subWinStep) const
//{
////! 20130910 implement sliding window with steps
//  Size3 sz(this->tsSize.height / subWinStep.height, this->tsSize.width/subWinStep.width, this->tsSize.depth/subWinStep.depth);
//  Tensor<T,cn> rst(sz);
//  gbuf.gI1.upload(*this);
//  gbuf.gI2.upload(ker);
//  gpu::Stream stream;
//  gpu::split(gbuf.gI1,gbuf.v,stream);
//  //Cube roi(0,0,0,subWinStep.height, subWinStep.width, subWinStep.depth);
//  Rect roi(0,0,ker.size().width,ker.size().height);
//  for (unsigned  int k=0; k<cn; k++)
//    {
//      for (unsigned int i=0;roi.y+roi.height<=this->tsSize.height;  ++i, roi.y+=subWinStep.height)
//        {
//          for (unsigned int j=0; roi.x+roi.width<=this->tsSize.width; ++j, roi.x+=subWinStep.width)
//            {
//              gpu::multiply(gbuf.v[k](roi),gbuf.gI2,gbuf.gs,1,-1,stream);
//              Scalar temp = gpu::sum(gbuf.gs);
//              rst(i,j,0)[k] = (T)temp[0];//since split channels, all data are 1 channel.
//            }
//          roi.x=0;//reset to left
//        }
//      roi.y=0;
//    }
//  return rst;
//}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalMean(const Size3& subWinSize,const Size3& subWinStep) const
{
	//tsSize.operator-(subWinSize);
	Size3_<double> temp1 = tsSize-subWinSize;
        Size3_<double> temp2 = subWinStep;
	//(Size3_<double>(Size3(tsSize-subWinSize))/Size3_<double>(Size3(subWinStep))).Ceil()+Size3(1,1,1);
	Size3 temp3 = (temp1/temp2).Ceil();
	Size3 rstSize = temp3 + Size3(1,1,1);
	Tensor<T,cn> rst(rstSize);
	for (int x=0; x< rstSize.height; x++)
		for (int y=0; y< rstSize.width; y++)
			for (int z=0; z<rstSize.depth; z++)
			{
				rst(x,y,z)= GetBlock(Cube(x*subWinStep.height,y*subWinStep.width,z*subWinStep.depth,
					x*subWinStep.height+subWinSize.height<tsSize.height?subWinSize.height:tsSize.height-x*subWinStep.height,
					y*subWinStep.width +subWinSize.width <tsSize.width ?subWinSize.width :tsSize.width -y*subWinStep.width,
					z*subWinStep.depth +subWinSize.depth <tsSize.depth ?subWinSize.depth :tsSize.depth -z*subWinStep.width)
					).Mean();
			}
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalVariance( const Tensor<T,cn>& mu, const Mat& ker, const Size3& subWinStep) const
{
    //! 20130911 implement sliding window with steps
    Size3 sz(this->tsSize.height / subWinStep.height, this->tsSize.width/subWinStep.width, this->tsSize.depth/subWinStep.depth);
    Tensor<T,cn> rst(sz);
    // rst.Print("rst init");
    //Cube roi(0,0,0,subWinStep.height, subWinStep.width, subWinStep.depth);
    Rect roi(0,0,ker.size().width,ker.size().height);
    //Tensor<T,cn> conj =(*this) * this->Conjugate();
    Tensor<T,cn> conj = this->Square();
    //conj.Print("conj");
    vector<Mat> mats;
    cv::split(conj,mats);

    for (unsigned int i=0;roi.y+roi.height<=this->tsSize.height;  ++i, roi.y+=subWinStep.height)
    {
      for (unsigned int j=0; roi.x+roi.width<=this->tsSize.width; ++j, roi.x+=subWinStep.width)
      {
          rst(i,j,0)[0]=(T)ker.dot(mats[0](roi));
      }
      roi.x=0;//reset to left
    }
   //auto temp2 =  ((*this)*(this->Conjugate())).Filter2D(ker).SubSample(subWinStep)- (mu*mu.Conjugate());
   auto temp = rst - mu.Square();
   //mu.Square().Print("mu2");
   //temp.Print("temp");
   //cout<<(temp-temp2).Abs().Sum()[0]<<endl;//verify
  return temp.Abs();
}
//template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalVarianceGPU( const Tensor<T,cn>& mu, const Mat& ker, BufferGPU& gbuf, const Size3& subWinStep) const
//{
//    //! 20130911 implement sliding window with steps
//    Size3 sz(this->tsSize.height / subWinStep.height, this->tsSize.width/subWinStep.width, this->tsSize.depth/subWinStep.depth);
//    Tensor<T,cn> rst(sz);

//    this->SquareGPU(gbuf);
//    gpu::Stream stream;
//    gpu::split(gbuf.gs,gbuf.v,stream);
//    gbuf.gI1.upload(ker);
//    Rect roi(0,0,ker.size().width,ker.size().height);
//    for (unsigned int i=0;roi.y+roi.height<=this->tsSize.height;  ++i, roi.y+=subWinStep.height)
//    {
//      for (unsigned int j=0; roi.x+roi.width<=this->tsSize.width; ++j, roi.x+=subWinStep.width)
//      {
//          gpu::multiply(gbuf.v[0](roi),gbuf.gI1,gbuf.gs,1,-1,stream);
//          Scalar temp = gpu::sum(gbuf.gs,gbuf.t1);
//          rst(i,j,0)[0]=(T)temp[0];
//          //rst(i,j,0)[1]=0;
//      }
//      roi.x=0;//reset to left
//    }
//  rst = rst - mu.SquareGPU(gbuf);
//  return rst.Abs();
//}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalVariance( const Tensor<T,cn>& mu, const Size3& subWinSize, const Size3& subWinStep) const
{
	return ((*this)*(this->Conjugate())).LocalMean(subWinSize,subWinStep)- (mu*mu.Conjugate());
}



////////2.5 comparison and metric //////////////////////////////////////////////
/// move to Metric module

/////////2.6 DSP ////////////////////////////////////////////////////////////
//
/*
 Various border types, image boundaries are denoted with '|'

 * BORDER_REPLICATE:     aaaaaa|abcdefgh|hhhhhhh
 * BORDER_REFLECT:       fedcba|abcdefgh|hgfedcb
 * BORDER_REFLECT_101:   gfedcb|abcdefgh|gfedcba
 * BORDER_WRAP:          cdefgh|abcdefgh|abcdefg
 * BORDER_CONSTANT:      iiiiii|abcdefgh|iiiiiii  with some specified 'i'
 */

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Laplacian(void) const
{
	Tensor<T,cn> rst(tsSize);
	for (int z=0; z<tsSize.depth; z++)
		cv::Laplacian(this->GetFrameRef(z),rst.GetFrameRef(z),this->depth(),1,1,0,BORDER_REFLECT_101);
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Filter2D(const Mat& ker, int boundary) const
{
		Tensor<T,cn> rst(tsSize);
		if(DataType<T>::type != ker.type())
		{
			CV_Error(CV_StsUnsupportedFormat, "bad type");
			//return rst;
		}

		for (int z=0; z< tsSize.depth; z++)
		{
			cv::filter2D(this->GetFrameRef(z),rst.GetFrameRef(z),-1,ker);
		}
		if (boundary == (int)FilterBoundary::FILTER_BOUND_EXTEND)
		{
			CV_Error(CV_StsNotImplemented,"unsupport extend boundary version");
			return rst;
		}
		else if (boundary == (int)FilterBoundary::FILTER_BOUND_VALID)
		{
		  //tsSize.Print();
		  //cout<<ker.size()<<endl;
			Size3 validSize(tsSize.height - ker.size().height +1, tsSize.width - ker.size().width+1, 1);
			Point3i validOffset(ker.size().height/2, ker.size().width/2,0);
			//validSize.Print();
			//cout<<validOffset<<endl;
			return rst.Crop(validOffset,validSize);
		}
		else
			return rst;
	
}
//this is slow
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::SubSample(const Size3& subrate)
{
	Size3_<double> temp1 = tsSize;
	Size3_<double> temp2 = subrate;
	Size3 rstSize = (temp1/temp2).Ceil();
	//temp1.Print();
	//temp2.Print();
	//rstSize.Print();
	if (rstSize.height==0 || rstSize.width==0 || rstSize.depth==0)
		CV_Error(CV_StsBadSize, "sample rate too high");
	Tensor<T,cn> rst(rstSize);
	for (int k=0; k*subrate.depth< tsSize.depth; k++)
		for (int i =0; i*subrate.height< tsSize.height; i++)
			for (int j=0; j*subrate.width< tsSize.width; j++)
			{
				rst(i,j,k) = this->operator()(i*subrate.height,j*subrate.width,k*subrate.depth);
			}
	return rst;
}
template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::Min() const
{
	if (cn != 1)
		CV_Error(CV_StsNotImplemented,"do not work for multi-channel yet");
	double p=DBL_MAX;
	double temp=0;
	for (int z = 0; z<this->size().depth; z++)
	{
		cv::minMaxLoc(this->GetFrameRef(z),&temp);
		if (temp < p)
			p= temp;
	}
	Vec<T,cn> rst;
	rst[0] = cv::saturate_cast<T>(p);
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::MaxClip(const Vec<T,cn>& value) const
{
	Tensor<T,cn> rst = this->Clone();
	for (int i=0; i< this->size().height; i++)
		for (int j=0; j<this->size().width; j++)
			for (int k=0; k<this->size().depth;k++)
			{
				//Vec<T,cn> temp = rst(i,j,k);
				//temp = (T)cv::norm(rst(i,j,k));
				//temp = (T)cv::norm(value);
				if (cv::norm(rst(i,j,k)) > cv::norm(value))
					rst(i,j,k)=value;
			}
	//rst.Print();
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::MinClip(const Vec<T,cn>& value) const
{
	Tensor<T,cn> rst = this->Clone();
	for (int i=0; i< this->size().height; i++)
		for (int j=0; j<this->size().width; j++)
			for (int k=0; k<this->size().depth;k++)
			{
				//Vec<T,cn> temp = rst(i,j,k);
				//temp = (T)cv::norm(rst(i,j,k));
				//temp = (T)cv::norm(value);
				if (cv::norm(rst(i,j,k)) < cv::norm(value))
					rst(i,j,k)=value;
			}
	//rst.Print();
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Flip(bool byX, bool byY, bool byZ) const
{
	Tensor<T,cn> rst = this->Clone();
	CV_Assert(!byZ);
	if (byX)
	{
	for (int y=0; y< this->size().width; y++)
		for (int x=0; x< this->size().height/2; x++)
		{	
			rst(x,y,0) = this->operator()(tsSize.height-x-1,y,0);
			rst(tsSize.height-x-1,y,0)=this->operator()(x,y,0);
		}
	}
	if (byY)
	{
	for (int x=0; x< this->size().height; x++)
		for (int y=0; y< this->size().width/2; y++)
		{	
			rst(x,y,0) = this->operator()(x,tsSize.width-y-1,0);
			rst(x,tsSize.width-y-1,0)=this->operator()(x,y,0);
		}
	}
	
	return rst;
}	
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::FreqComplexFilter(cv::Mat& kernel, bool conj) const
{
	Tensor<T,cn> rst(size());
	CV_Assert(this->channels() == kernel.channels() && this->channels()== 2); 
	for (int i=0; i< size().depth; i++)
		rst.SetFrame(i,mylib::FreqComplexFilter(GetFrameRef(i),kernel,conj));
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::DFT(void) const
{

	Tensor<T,cn> rst(size());
	CV_Assert(cn ==2 || cn== 6); 
	for (int i=0; i< size().depth; i++)
	{
		#ifdef USE_GPU
		  gbuf.gI1.upload(GetFrameRef(i));
		  gbuf.gI1.convertTo(gbuf.t1, CV_32FC2, gbuf.stream);
		  gpu::dft(gbuf.t1, gbuf.t2, this->tsSize,0,gbuf.stream);
		  gbuf.t2.convertTo(gbuf.gI2, CV_64FC2,gbuf.stream);//make this adaptive
		  gbuf.gI2.download(rst[i]);
		#else
		cv::dft(GetFrameRef(i),rst[i]);
		#endif

	}
	//gbuf.release();
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::IDFT(void) const

{
	Tensor<T,cn> rst(size());
	CV_Assert(cn ==2 || cn== 6); 
	for (int i=0; i< size().depth; i++)
	{

		#ifdef USE_GPU
		gbuf.gI1.upload(GetFrameRef(i));
		gbuf.gI1.convertTo(gbuf.t1, CV_32FC2, gbuf.stream);
		gpu::dft(gbuf.t1, gbuf.t2, this->tsSize, DFT_INVERSE, gbuf.stream);
		gbuf.t2.convertTo(gbuf.gI2, CV_64FC2,gbuf.stream);
		gpu::split(gbuf.gI2, gbuf.v, gbuf.stream);
		gpu::divide(gbuf.v[0], double(rst[i].size().height*rst[i].size().width), gbuf.v[0],1, -1, gbuf.stream);
		gpu::divide(gbuf.v[1], double(rst[i].size().height*rst[i].size().width), gbuf.v[1],1, -1, gbuf.stream);
		gpu::merge(gbuf.v, gbuf.gI2, gbuf.stream);
		gbuf.gI2.download(rst[i]);
		//vector<Mat> temp;
                //cv::split(rst[i],temp);
                //Tensor<double,1>(temp[0]/temp[0].size().height/temp[0].size().width).Display(); 
                //mylib::DisplayMat(temp[0]);
                //cv::divide(this->tsSize.area(), rst[i], rst[i], -1);
		#else	
		cv::dft(GetFrameRef(i),rst[i],DFT_INVERSE|DFT_SCALE);
		#endif
	}

	return rst;
}

template<class T, size_t cn> Mat Tensor<T,cn>::DFTShift(const cv::Mat& A) const
{
	cv::Mat rst = cv::Mat::zeros(A.size(),A.type());
	//only works when size is even
	cv::Rect roi1(0,0,A.size().width/2, A.size().height/2);
	cv::Rect roi2(A.size().height/2,0,A.size().width/2, A.size().height/2);
	cv::Rect roi3(0,A.size().width/2,A.size().width/2, A.size().height/2);
	cv::Rect roi4(A.size().height/2, A.size().width/2,A.size().width/2, A.size().height/2);
	A(roi4).copyTo(rst(roi1));
	A(roi3).copyTo(rst(roi2));
	A(roi2).copyTo(rst(roi3));
	A(roi1).copyTo(rst(roi4));
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::DFTShift(void) const
{
	Tensor<T,cn> rst(size());
	for (int i=0; i< size().depth; i++)
	{
		rst.SetFrame(i,DFTShift(GetFrameRef(i)));
	}
	return rst;
}


template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Normalize(void) const
{
	CV_Assert(cn == 1);
	Tensor<T,cn> rst(size());
	for (int i=0; i< size().depth; i++)
	{
		cv::normalize(GetFrameRef(i),rst[i],0,255,NORM_MINMAX);
	}
	return rst;
}
///////////////
template class Tensor<uchar,1>;
template class Tensor<uchar,2>;
template class Tensor<uchar,3>;
template class Tensor<float,1>;
template class Tensor<float,2>;
template class Tensor<float,3>;
template class Tensor<double,1>;
template class Tensor<double,2>;
template class Tensor<double,3>;


/*
#ifdef USE_GPU
template<> typename Tensor<uchar,1>::BufferGPU Tensor<uchar,1>::gbuf;//declare static buffer
template<> typename Tensor<uchar,2>::BufferGPU Tensor<uchar,2>::gbuf;//declare static buffer
template<> typename Tensor<uchar,3>::BufferGPU Tensor<uchar,3>::gbuf;//declare static buffer
template<> typename Tensor<float,1>::BufferGPU Tensor<float,1>::gbuf;//declare static buffer
template<> typename Tensor<float,2>::BufferGPU Tensor<float,2>::gbuf;//declare static buffer
template<> typename Tensor<float,3>::BufferGPU Tensor<float,3>::gbuf;//declare static buffer
template<> typename Tensor<double,1>::BufferGPU Tensor<double,1>::gbuf;//declare static buffer
template<> typename Tensor<double,2>::BufferGPU Tensor<double,2>::gbuf;//declare static buffer
template<> typename Tensor<double,3>::BufferGPU Tensor<double,3>::gbuf;//declare static buffer
#endif
*/
}

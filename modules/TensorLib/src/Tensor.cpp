#include "Tensor.h"
#include "Steerable.h"
#include "LRI.h"
#include "Metric.h"

template<class T, size_t cn> Mat Tensor<T,cn>::stsim2_lse_weight = Mat(); 
/////////////2.1 constructor, copy  assignment and type conversion ////////////////////////
template<class T, size_t cn> Tensor<T,cn>::Tensor(void):Mat()
{
	tsOffset = Point3i();
	cFileName = "unknown";
	tsSize = Size3();
	mxFrame = Mat();
	subWinSize = Size3();
	subWinStep = Size3();
	debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(int height, int width, int depth, typename Tensor<T,cn>::c_ref_type val): Mat()
{
	const int sz[] = {depth, height, width}; //important , order changed
	*this = Mat(3, sz,  CV_MAKETYPE(DataType<T>::depth,cn), Scalar(val));
	tsSize = Size3(height,width,depth);
	cFileName = "unknown";
	tsOffset = Point3i();
	mxFrame = Mat(height,width,CV_MAKETYPE(DataType<T>::depth,cn));
	subWinSize = Size3();
	subWinStep = Size3();
		debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(const Size3& size, const Vec<T,cn>& val): Mat()
{
	const int sz[] = {size.depth, size.height,size.width};
	*this = Mat(3, sz, CV_MAKETYPE(DataType<T>::depth,cn),Scalar(val));
	this->tsSize = size;
	cFileName = "unknown";
	tsOffset = Point3i();
	mxFrame = Mat(size.height,size.width,CV_MAKETYPE(DataType<T>::depth,cn));
	subWinSize = Size3();
	subWinStep = Size3();
		debugtrigger=false;
}
template<class T, size_t cn> Tensor<T,cn>::Tensor(const Size3& size):Mat()
{
 	const int sz[] = {size.depth, size.height,size.width};
	*this = Mat(3, sz, CV_MAKETYPE(DataType<T>::depth,cn));
	this->tsSize = size;
	cFileName = "unknown";
	tsOffset = Point3i();
	mxFrame = Mat(size.height,size.width,CV_MAKETYPE(DataType<T>::depth,cn));
	subWinSize = Size3();
	subWinStep = Size3();
	debugtrigger=false;
}
template<class T, size_t cn> Tensor<T,cn>::Tensor(const Tensor& ts): Mat(ts)
{
	this->tsSize = ts.tsSize;
	if ( ts.type() != CV_MAKETYPE(DataType<T>::depth,cn))
	{
		Tensor<T,cn> rst(tsSize);
		this->convertTo(rst,DataType<T>::depth);
		*this = rst;
	}
	this->cFileName = ts.cFileName;
	tsOffset = ts.tsOffset;
	mxFrame = ts.mxFrame;
	SetSubWinSize(ts.subWinSize);
	SetSubWinStep(ts.subWinStep); 
	lightTag = ts.lightTag;
	lightCan = ts.lightCan;
		debugtrigger=ts.debugtrigger;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(const Mat& mt): Mat(mt)
{
	//include the type conversion function
	if (mt.dims == 2)
	{
		this->tsSize = Size3(mt.rows,mt.cols,1);
		int sz[] = {1,mt.rows,mt.cols};
		this->create(3,sz,mt.type());
		SetFrame(0,mt);
		
	}
	else
	{
		this->tsSize = Size3(mt.size[1],mt.size[2],mt.size[0]);
	}
	cFileName = "unknown";
	tsOffset = Point3i();
	if (mt.type() != CV_MAKETYPE(DataType<T>::depth,cn))
	{
		Tensor<T,cn> rst(size());
		this->convertTo(rst,DataType<T>::depth);
		*this = rst;
	}
	subWinSize = Size3();
	subWinStep = Size3();
	//mxFrame = Mat(mt.size().height,mt.size().width,CV_MAKETYPE(DataType<T>::depth,cn));
		debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(const string cFileName)//:subTensor(Tensor_<T,cn>())
{

	Load(cFileName);
	tsOffset = Point3i();

}
template<class T, size_t cn> Tensor<T,cn>::~Tensor(void)
{
}

template<class T, size_t cn> Tensor<T,2> Tensor<T,cn>::ToComplex(void) const
{
	CV_DbgAssert( 2 == 2*cn || 2 ==cn);
	if(cn == 2)
		return *this;
	Tensor<T,2> rst(size());
	rst.SetFileName(this->cFileName);
	for (int i=0; i< size().depth; i++)
	{
		rst.SetFrame(i,mylib::toComplex(GetFrame(i)));
	}
	rst.debugtrigger = this->debugtrigger;
	return rst;
}

template<class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::operator= (const Tensor<T,cn>& ts)
{
//	CV_DbgAssert(ts.type() == CV_MAKETYPE(this->depth(), this->channels()));
	this->Mat::operator=(ts);
	tsSize = ts.tsSize;
	cFileName = ts.cFileName;
	tsOffset = ts.tsOffset;
	this->subWinSize = ts.subWinSize;
	this->subWinStep = ts.subWinStep;
		debugtrigger=ts.debugtrigger;
	return *this;
}

////////////2.2 public utilities	/////////////////////////////////////////////////////

template<class T, size_t cn> void Tensor<T,cn>::Display(int flag) const
{
	string tempName;
	CV_DbgAssert(cn == this->channels());
	Tensor<uchar,cn> tempTs;

	if (this->depth() != DataType<uchar>::depth)
		tempTs = Tensor<uchar,cn>(*this); //directly type conversion
	else
		tempTs = *this;
	
	if (cFileName.empty())
		tempName="Unknown";

	cv::namedWindow(tempName,flag|CV_WINDOW_FREERATIO|CV_GUI_EXPANDED);
	cvMoveWindow(tempName.c_str(),100,100);
	
	if (tsSize.depth==1)
	{
		cv::imshow(tempName,tempTs.GetFrame(0));
		cv::waitKey();
		cvDestroyWindow(tempName.c_str());
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
		cvDestroyWindow(tempName.c_str());
	}

}

template<class T, size_t cn> void Tensor<T,cn>::Display(int sec, int flag) const
{
	string tempName;
	CV_DbgAssert(cn == this->channels());

	Tensor<uchar,cn> tempTs;
	if (this->depth() != DataType<uchar>::depth)
		tempTs = Tensor<uchar,cn>(*this); //directly type conversion
	else
		tempTs = *this;
	

	if (cFileName.empty())
		tempName="Unknown";
	cv::namedWindow(tempName,flag|CV_WINDOW_KEEPRATIO|CV_GUI_EXPANDED);
	cvMoveWindow(tempName.c_str(),100,100);

	if (tsSize.depth==1)
	{
		cv::imshow(tempName,tempTs.GetFrame(0));
		cv::waitKey(sec);
		cvDestroyWindow(tempName.c_str());
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
		cvDestroyWindow(tempName.c_str());
	}

}

template<typename T, size_t cn> void Tensor<T,cn>::DisplayAll(std::vector<Tensor<T,cn>>& vts, int row, int col, bool save, string savename) const
{
  Size3 padsz(0,0,0);
  CV_DbgAssert(vts.size()<=row*col);
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
	cap>>mxFrame;//read next frame
	cv::Mat grayone;
	if (!mxFrame.data) //if cannot be read as a stream
	{
		cout<<"cannot be read by stream, use imread()"<<endl;
		if (cn == 3)
			mxFrame = cv::imread(cFileName,1);
		else if (cn == 1)
			mxFrame = cv::imread(cFileName,0);
		else if (cn==2) //complex
		{
			mxFrame = cv::imread(cFileName,0);
			Mat azMat = Mat::zeros(mxFrame.size(),mxFrame.type());
			vector<Mat> vcm;
			vcm.push_back(mxFrame);
			vcm.push_back(azMat);
			cv::merge(vcm,mxFrame);
		}
		else
			mxFrame = cv::imread(cFileName,-1);

		CV_Assert(!this->data && mxFrame.channels() == cn ); 
		const int sz[] = {1,mxFrame.rows,mxFrame.cols};
		mxFrame.convertTo(mxFrame,DataType<T>::depth);//new added
		this->create(3,sz,mxFrame.type());
		this->tsSize = Size3(sz[1],sz[2],sz[0]);
		SetFrame(0,mxFrame.clone());
		cout<<"size: "; 
		tsSize.Print();
		cout<<"type: "<<this->type();
		cout<<", channel: "<<this->channels()<<endl;
	}
	else
	{
		std::cout<<"read as stream\n";
		//cout<<mxFrame.channels()<<","<<mxFrame.type()<<endl;
		if (cn==1&&mxFrame.channels()==3)
		{
			cv::cvtColor(mxFrame,grayone,CV_RGB2GRAY);
			mxFrame = grayone;
		//	cout<<mxFrame.type();
		}
		if (cn==2) //complex
		{
			cv::cvtColor(mxFrame,grayone,CV_RGB2GRAY);
			Mat azMat = Mat::zeros(grayone.size(),grayone.type());
			vector<Mat> vcm;
			vcm.push_back(grayone);
			vcm.push_back(azMat);
			cv::merge(vcm,mxFrame);
		}
		int sz[] = {1,mxFrame.rows,mxFrame.cols};
    mxFrame.convertTo(mxFrame,DataType<T>::depth);//new added
    //cout<<mxFrame.channels()<<","<<mxFrame.type()<<endl;
		this->create(3,sz,mxFrame.type());
		this->tsSize = Size3(sz[1],sz[2],sz[0]);
		SetFrame(0,mxFrame.clone());
    //this->Print();
		cap>>mxFrame;
		//while(mxFrame.data)
		//{	
		//	if (cn==1&&mxFrame.channels()>1)
		//	{
		//		cv::cvtColor(mxFrame,grayone,CV_RGB2GRAY);
		//		mxFrame = grayone;	
		//	}
  //    Mat tempFrame(3,sz,mxFrame.type(),mxFrame.data);
		//	this->push_back(tempFrame.clone()); //untested
		//	cap>>mxFrame;
		//	tsSize.depth++; //measure the length of video
		//}
		*this = Tensor<T,cn>(*this);//new added
	}
	cap.release();
		debugtrigger=false;
	//this->Display();
}

template<class T, size_t cn> Size3 Tensor<T,cn>::size(void) const 
{
	return this->tsSize;
}


template<class T, size_t cn> typename Tensor<T,cn>::ref_type Tensor<T,cn>::operator() (int x, int y, int z)
{
	CV_DbgAssert( CV_MAKETYPE(DataType<T>::depth,cn) == this->type());
	return this->at<value_type>(z,x,y);

}

template<class T, size_t cn> typename Tensor<T,cn>::c_ref_type Tensor<T,cn>::operator() (int x, int y, int z) const
{
	int t = this->type();
	int t2 = CV_MAKETYPE(DataType<T>::depth,cn) ;
	CV_DbgAssert( CV_MAKETYPE(DataType<T>::depth,cn) == this->type());

	return this->at<value_type>(z,x,y);

}

template<class T, size_t cn> typename Tensor<T,cn>::ref_type Tensor<T,cn>::operator[] (const Point3i& pos)
{
	return this->operator()(pos.x,pos.y,pos.z);
}
template<class T, size_t cn> typename Tensor<T,cn>::c_ref_type Tensor<T,cn>::operator[](const Point3i& pos) const
{
	return this->operator()(pos.x,pos.y,pos.z);
}

template< class T, size_t cn> Mat Tensor<T,cn>::GetFrame(int i)
{
	if (dims <=2 )
		return *this;
	else
	{
		CV_DbgAssert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
    
		//mxFrame = this->row(i);
    mxFrame =  Mat(*this, Range(i, i+1), Range::all());
    //mylib::DisplayMat(mxFrame);
		Mat rst = Mat(2,mxFrame.size.p+1,mxFrame.type(),mxFrame.data,mxFrame.step.p+1);
    //if (rst.size().height<10)
    //mylib::DisplayMat(rst);
		return rst;
	}
}

template< class T, size_t cn> Mat& Tensor<T,cn>::GetFrameRef(int i) 
{
	if (dims <=2 )
		return *this;
	else
	{
		CV_DbgAssert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
		//mxFrame = this->row(i);
     //Mat* arrays[] = {this};
    Mat* arrays[] = {(Mat*)this};
    Mat planes[1];
    NAryMatIterator it((const Mat**)arrays, planes,1);
    mxFrame = it.planes[i].reshape(cn,this->size().height);
    //if (mxFrame.size().height<0)
    //mylib::DisplayMat(mxFrame);
    return mxFrame;
	}
}

template< class T, size_t cn> void Tensor<T,cn>::GetFrameRef(int i, Mat& rst) const 
{
	if (dims <=2 )
		rst= *this;
	else
	{
		CV_DbgAssert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
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
		return *this;
	else
	{
		CV_DbgAssert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
    Mat tempMat = Mat(*this, Range(i, i+1), Range::all());
	  //Mat tempMat = this->row(i);
		Mat rst = Mat(2, tempMat.size.p+1,tempMat.type(),tempMat.data,tempMat.step.p+1);
		return rst;
	}
}

template< class T, size_t cn> Mat Tensor<T,cn>::operator[](int i)
{
	return GetFrame(i);
}

template< class T, size_t cn> const Mat Tensor<T,cn>::operator[](int i) const
{
	return GetFrame(i);
}

template< class T, size_t cn> Point3i Tensor<T,cn>::offset() const
{
	return this->tsOffset;
}

template< class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Clone(void) const
{
  std::lock_guard<std::mutex> lock(tsmutex);
	//tsmutex.lock();
  Tensor<T,cn> rst(this->Mat::clone());
	rst.SetFileName(this->cFileName);
	rst.SetOffset(this->tsOffset);
	rst.SetSubWinSize(subWinSize);
	rst.SetSubWinStep(subWinStep); 
	rst.lightTag = lightTag;
	rst.lightCan = lightCan; 

	return rst;
}

template< class T, size_t cn> void Tensor<T,cn>::SetFrame(int i, const Mat& frm)
{
	//use to modify(or assignment by copy) one frame
	if (dims <=2 )
		*this = frm;
	else
	{
		CV_DbgAssert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
		Mat mxFrame = GetFrame(i); //get a reference to the frame
		frm.copyTo(mxFrame);
		//make sure frm is nonempty
        //CV_DbgAssert(frm.total()!=0&&frm.isContinuous() && mxFrame.isContinuous());
        //size_t planeSize = frm.elemSize() * frm.rows * frm.cols;
		//memcpy(mxFrame.data,frm.data,planeSize);
	}
}


template<class T, size_t cn> void Tensor<T,cn>::Print(const string& fname, bool tofileonly) const
{

	fstream logfile;
	logfile.open(".\\"+fname+".txt",ios::out); 
	//bool tofileonly=false;
	if (!tofileonly)
	{
		cout<<"File: "<<fname<<" Display "<<size().height <<" "<<size().width<<" Data type: "<<type()<<endl;
		cout<<"cn="<<cn<<endl;
	}
	//logfile<<"File: "<<cFileName<<"Display "<<size().height <<" "<<size().width<<" Data type: "<<type()<<endl;
	
	for (int d=0; d< size().depth; d++) 
	{
		if (!tofileonly)	
			cout<<"frame "<<d<<" :"<<endl;
		//logfile<<"frame "<<d<<" :"<<endl;
		for (int c=0; c< cn;c++)
		{
			if (!tofileonly)
				cout<<"channel "<<c<<" :"<<endl;
			//logfile<<"channel "<<c<<" :"<<endl;
			for (int x = 0; x < size().height; x++)
			{
				for (int y = 0; y < size().width; y++)
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

//template <class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::GetBlockRef(const Cube& roi)
//{
//		CV_DbgAssert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
//		//mxFrame = this->row(i);
//     //Mat* arrays[] = {this};
//    Mat* arrays[] = {(Mat*)this};
//    Mat planes[1];
//    NAryMatIterator it((const Mat**)arrays, planes,1);
//    mxFrame = it.planes[i].reshape(cn,this->size().height);
//    //mylib::DisplayMat(mxFrame);
//    return mxFrame(Rect(roi.y,roi.x,roi.width,roi.height));
//}
template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::GetBlock(const Cube& roi)
{
	Range r[] = {Range(roi.z,roi.z+roi.depth),Range(roi.x,roi.x+roi.height),Range(roi.y,roi.y+roi.width)};	
	Tensor<T,cn> rst = this->Mat::operator()(r);
	rst.SetFileName(this->cFileName);
	rst.SetOffset(roi.offset());
	rst.SetSubWinSize(this->subWinSize);
	rst.SetSubWinStep(this->subWinStep);
	return rst;
}
template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::GetBlock(const Cube& roi) const
{
	Range r[] = {Range(roi.z,roi.z+roi.depth),Range(roi.x,roi.x+roi.height),Range(roi.y,roi.y+roi.width)};	
	Tensor<T,cn> rst = this->Mat::operator()(r);
	rst.SetFileName(this->cFileName);
	rst.SetOffset(roi.offset()+this->offset());
	rst.SetSubWinSize(this->subWinSize);
	rst.SetSubWinStep(this->subWinStep);
	return rst;
}

template<class T, size_t cn> void Tensor<T,cn>::Ref(const Cube& roi, Tensor<T,cn>& dst) const
{
	Range r[] = {Range(roi.z,roi.z+roi.depth),Range(roi.x,roi.x+roi.height),Range(roi.y,roi.y+roi.width)};	
	dst.Mat::operator=(this->Mat::operator()(r));
  dst.tsSize = roi.size();
  dst.tsOffset = roi.offset();
	dst.SetFileName(this->cFileName);
	dst.SetOffset(roi.offset()+this->offset());
	dst.SetSubWinSize(this->subWinSize);
	dst.SetSubWinStep(this->subWinStep);
}

template <class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::SetBlock(const Point3i& pos, const Tensor<T,cn>& ts)
{
	CV_DbgAssert(ts.type() == type());
	AssertRange(pos,ts.size());
	
	Cube roi(pos,ts.size());
	cv::Mat tempMat = GetBlock(roi);

	ts.copyTo(tempMat);
	return *this;
}
template <class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::SetBlock(const Tensor<T,cn>& ts)
{
	return SetBlock(Point3i(0,0,0),ts);
}
template<class T, size_t cn> void Tensor<T,cn>::AssertRange(const Point3i& pos, const Size3& sz) const
{
	CV_DbgAssert( pos.x + sz.height <= tsSize.height &&
				  pos.y + sz.width  <= tsSize.width &&
				  pos.z + sz.depth  <= tsSize.depth);
}

template<class T, size_t cn> void Tensor<T,cn>::SaveBlock(const string& cFileName, bool isGray)
{
	if (size().depth>1)
	{
		cv::VideoWriter vd(cFileName+".avi",CV_FOURCC('D','I','B',' '),30, cv::Size(size().width,size().height),true);
		for (int i=0;i<size().depth;i++)
		{

			if (isGray)
			{
				Mat grayone = GetFrame(i);
				cv::cvtColor(grayone,grayone,CV_RGB2GRAY);
				vd<<grayone;
			}
			else
				vd<<GetFrame(i);

		}
	}
	else
	{
		if (isGray)
		{
			Mat grayone = GetFrame(0);
			cv::cvtColor(grayone,grayone,CV_RGB2GRAY);
			cv::imwrite(cFileName,grayone);
		}
		else
			cv::imwrite(cFileName,GetFrame(0));
	}
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Crop(const Point3i& pos, const Size3& sz) const
{
  return GetBlock(Cube(pos,sz)).Clone();
}

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

template<class T, size_t cn> bool Tensor<T,cn>::IsInside(const Point3i& pos) const
{
	if ( pos.x >= tsOffset.x && pos.x < tsOffset.x + tsSize.height &&
		pos.y >= tsOffset.y && pos.y < tsOffset.y + tsSize.width  &&
		pos.z >= tsOffset.z && pos.z < tsOffset.z + tsSize.depth )
		return true;
	else
		return 0;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Row(int i)
{
	Tensor<T,cn> rst(Size3(1,tsSize.width,tsSize.depth));
	if (tsSize.depth == 1)
		rst = this->row(i);
	else
	{
		Range r[] = {Range::all(),Range(i,i+1),Range::all()};
		rst = this->Mat::operator()(r);
	}
	//this->Ref(Cube(i,0,0,1,tsSize.width,tsSize.depth), *rst);
	return rst;
}


template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator()(const Cube& roi)
{	
	//return GetBlock(roi);
	AssertRange(roi.offset(),roi.size());
	Tensor<T,cn> dst(roi.size());
	for (int i=0; i< roi.depth; i++)
	{
		//if ( i < dst.size().depth )
		dst.SetFrame(i,this->GetFrame(i+roi.z)(Rect(roi.y,roi.x,roi.width,roi.height)));
		//else
		//	dst.vdStream.push_back((*this)[i](Rect(roi.y,roi.x,roi.width,roi.height)));
	}
	dst.tsOffset = roi.offset()+this->offset();
	dst.tsSize = roi.size();
	dst.cFileName = this->cFileName;	
	dst.SetSubWinSize(this->subWinSize);
	dst.SetSubWinStep(this->subWinStep);
	return dst; 
}

////////2.4 arithmatic operations //////////////////////////////////////////

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator+(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) + Tensor<T,cn>(size().height,size().width,size().depth,s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator-(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) - Tensor<T,cn>(size().height,size().width,size().depth,s);
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

	for (int z =0; z < size().depth; z++)
	{
		if (cn==2 || cn == 6)
			cv::mulSpectrums(rst.GetFrame(z),ts.GetFrame(z),rst.GetFrame(z),DFT_ROWS,false); // complex mul
		else
    {
			rst.SetFrame(z,GetFrame(z).mul(ts.GetFrame(z)));
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
			cv::mulSpectrums(ts.GetFrame(i),ts.GetFrame(i),tempMat,DFT_ROWS,true);
			//mylib::DisplayMat(tempMat);
			vector<Mat> temp;
			cv::split(tempMat, temp);
			temp[1] = temp[0];
			merge(temp,tempMat);
			//mylib::DisplayMat(tempMat);
			cv::mulSpectrums(rst.GetFrame(i),ts.GetFrame(i),rst.GetFrame(i),DFT_ROWS,true);
			//mylib::DisplayMat(rst.GetFrame(i));
			//Vec<T,cn> a= rst(0,0,0);
			//Vec<T,cn> b = tempMat.at<Vec<T,cn>> (0,0);
			rst.GetFrame(i) /= tempMat;
			//mylib::DisplayMat(rst.GetFrame(i));
		}
	}
	else
	{
		for (int i = 0; i< tsSize.depth; i++)
			rst.GetFrame(i) /= ts.GetFrame(i);
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
	return (*this)+All(s);
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator-(T s) const
{
	return (*this)-All(s);
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator*(T s) const
{
	return (*this)*All(s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator/(T s) const
{
	return (*this)/All(s);
}

template<class T, size_t cn> typename Tensor<T,cn>::value_type Tensor<T,cn>::All(T val) const
{
	return value_type::all(val);
}

template<class T, size_t cn> typename Tensor<T,cn> Tensor<T,cn>::Abs(void) const
{
	
	Tensor<T,cn> rst(this->size());
	if (cn == 2) // if data type is not float/double this will throw type exception
	{
		for (int i=0; i< size().depth; i++)
		{
			cv::mulSpectrums(this->GetFrame(i),this->GetFrame(i),rst.GetFrame(i),DFT_ROWS,true);
			cv::pow(rst.GetFrame(i),0.5,rst.GetFrame(i));
		}
	}
	else
	{
		for (int i=0; i< size().depth; i++)
			rst.SetFrame(i,cv::abs(this->GetFrame(i)));
	}
	return rst;
}


template<class T, size_t cn> typename Tensor<T,cn> Tensor<T,cn>::AbsDiff(c_ref_type s) const
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
			cv::absdiff(this->GetFrame(i),v,rst.GetFrame(i));
	}
	return rst;
}


template<class T, size_t cn> typename Tensor<T,cn> Tensor<T,cn>::AbsDiff(const Tensor& ts) const
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
			cv::absdiff(this->GetFrame(i),ts.GetFrame(i),rst.GetFrame(i));
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
					divide( 1., this->GetFrame(i), rst.GetFrame(i));
				if( ipower == -1 )
					return rst;
				ipower = -ipower;
			}
			else
				rst = this->Clone();

			switch( ipower )
			{
			case 0:
				rst = Tensor<T,cn>(tsSize.height,tsSize.width,tsSize.depth,Vec<T,cn>(1,0));
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
  cv::log(*this,rst);
  return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Exp(void) const
{
  Tensor<T,cn> rst(this->size());
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
	return mylib::VecDiv(Sum(),T(tsSize.area()));
}

template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::Var(void) const
{
	Vec<T,cn> mu, var;
	//this->Print();
	cv::meanStdDev(*this,mu,var);
	return mylib::VecMul(var,var);
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
		rst.subWinSize=subWinSize;
		rst.subWinStep=subWinStep;
	}
	//	(Size3(this->size().width,this->size().height,this->size().depth));
	for (int i=0; i< size().depth; i++)
		cv::transpose(GetFrame(i),rst[i]);
	return rst;
}


template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::ExtendBoundary(Size3 extSz, typename Tensor<T,cn>::value_type val) const
{
	Tensor<T,cn> rst(size()+(extSz*2),val);
	rst.tsOffset = Point3i(extSz.height,extSz.width,extSz.depth);
	rst.SetBlock(rst.tsOffset,*this);
	rst.SetSubWinSize(subWinSize);
	rst.SetSubWinStep(subWinStep);
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
	rst.SetSubWinSize(subWinSize);
	rst.SetSubWinStep(subWinStep);
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalMean(const Mat& ker,const Size3& subWinStep) const
{
	return Filter2D(ker).SubSample(subWinStep);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalMean(const Size3& subWinSize,const Size3& subWinStep) const
{
	//tsSize.operator-(subWinSize);
	Size3 rstSize = (Size3_<double>(tsSize-subWinSize)/Size3_<double>(subWinStep)).Ceil()+Size3(1,1,1);
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
	auto temp =  ((*this)*(this->Conjugate())).Filter2D(ker).SubSample(subWinStep)- (mu*mu.Conjugate());
  return temp.Abs();
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalVariance( const Tensor<T,cn>& mu, const Size3& subWinSize, const Size3& subWinStep) const
{
	return ((*this)*(this->Conjugate())).LocalMean(subWinSize,subWinStep)- (mu*mu.Conjugate());
}



////////2.5 comparison and metric //////////////////////////////////////////////


template <class T, size_t cn> bool Tensor<T,cn>::ComputeAIM(const Tensor& ts, double thrd) const
{
	//important! this Adaptive Interpolation Metric only works valid for 2D image!
	//this->AbsDiff(ts).Print(true);
  //ts.Print();
  //this->Print();
  //this->AbsDiff(ts).Print();
	Tensor<uchar,1> tempMap = this->AbsDiff(ts).CompareElement(thrd);
  //tempMap.Print();
	Tensor<uchar,1> extMap = tempMap.ExtendBoundary(Size3(1,1,0),255);
  //extMap.Print();
	//extMap.Print(true);
	for (int z = 0 ; z < size().depth; z++)
	{
		for (int x = 1; x < size().height+1; x++)
			for (int y=1; y < size().width+1; y++)
			{
				if (extMap(x,y,z)[0]== 0) //test 1 failed
				{
					if (extMap(x+1,y,z)[0]/255 + extMap(x,y+1,z)[0]/255
						+ extMap(x-1,y,z)[0]/255 + extMap(x,y-1,z)[0]/255 < 3)
						return false;
				}
			}
	}
	return true;
}

template<class T, size_t cn> Mat Tensor<T,cn>::CompareElement(const Mat& inMat, double thrd, int flag) const
{
	cv::Mat tempMat = Mat::ones(inMat.size(),CV_8U);
	vector<Mat> A,B(inMat.channels());
	cv::split(inMat,A);
	for (int c = 0; c < inMat.channels(); c++)
	{
		compare(A[c],thrd,B[c],flag);
    //mylib::DisplayMat(B[c]);
		tempMat = tempMat.mul(B[c]/255);
   // mylib::DisplayMat(tempMat);
	}
	tempMat = tempMat*255;
  //mylib::DisplayMat(tempMat);
	return tempMat;
}


template<class T, size_t cn> Tensor<uchar,1> Tensor<T,cn>::CompareElement(double thrd, int flag) const
{
	Tensor<uchar,1> rst(size());
	for (int i=0; i< size().depth; i++)
	{
		rst.SetFrame(i, CompareElement(GetFrame(i),thrd,flag));
   // mylib::DisplayMat(rst[i]);
	}
  //rst.Print();
	return rst;
}


template <class T, size_t cn> double Tensor<T,cn>::ComputeMSE(const Tensor& ts) const
{
	double distance=0;
	CV_DbgAssert(ts.channels() == cn);
	CV_DbgAssert(ts.size()==tsSize);
	//if (this->offset().x == DEBUG_X && this->offset().y== DEBUG_Y && this->size().height==DEBUG_SIZE)
	if (ts.debugtrigger==true||this->debugtrigger==true)
	{
		this->Print("bound tar",true);
		ts.Print("bound cand",true);
	}
	Scalar temp; 
	Mat v;
  Mat v1,v2;
  Mat v1_dft, v2_dft, v_dft;
  //int pnum = 4;
  //vector<thread> threads;
  //vector<Vec<T,cn>> sums(pnum);
	for (int z = 0; z < tsSize.depth; z++)
	{
		//cv::absdiff(GetFrame(z),ts[z],v);
   /* for (int p=0; p<pnum; p++)
    {
      threads.push_back(thread([&](int p, int z, int H, int W){
        Vec<T,cn> tt=0;
        for (int i=p; i<H; i+=pnum)
          for (int j=p; j<W; j+=pnum)
          {
            tt =  this->operator()(i,j,z)-ts(i,j,z);
            sums[p] += mylib::VecMul<T,cn>(tt,tt);
          }
      },p,z,ts.size().height,ts.size().width));
    }*/
    v1 = this->GetFrame(z);
    ts.GetFrameRef(z,v2);
		v = v1-v2;
    cv::pow(v,2.0,v);
    //cv::multiply(v1,v2,v);
    //cv::dft(v1,v1_dft);
    //cv::dft(v2,v2_dft);
   /* for (int p=0; p<pnum;p++)
    {
      threads[p].join();
      temp+=sums[p];
    }*/
    //temp[0]+=cv::norm(v1,v2,NORM_L1);
		temp +=cv::sum(v);
	}
	for (int i =0; i < cn; i++)
		distance+= temp[i];
	return distance;///tsSize.volumn()/cn;
}

template <class T, size_t cn> double Tensor<T,cn>::ComputeSAD(const Tensor& ts) const
{
	double distance=0;
	CV_DbgAssert(ts.channels() == cn);
	CV_DbgAssert(ts.size()==tsSize);
	Scalar temp; 
	Mat v;
	for (int z = 0; z < tsSize.depth; z++)
	{
		cv::absdiff(GetFrame(z),ts[z],v);
		temp +=cv::sum(v);
	}
	for (int i =0; i < cn; i++)
		distance+= temp[i];
	return distance/tsSize.volumn()/cn;
}

template <class T, size_t cn> double Tensor<T,cn>::ComputePSNR(const Tensor& ts) const
{
	double distance = ComputeMSE(ts);
  distance/= double(ts.size().area());//do it since ComputeMSE does not normalize
	return 10*log10(std::pow(255.0,2.0)/distance);
}
//gj01072013_2 add LRI based metric to Tensor::computeLRI()
template <class T, size_t cn> double Tensor<T,cn>::ComputeLRI(const Tensor& ts) const
{
  LRI lri;
  double rst=0;
  for (int z=0; z<tsSize.depth; z++)
  {
    rst+=lri.computeNewMetric(this->GetFrame(z),ts.GetFrame(z));
  }
  return rst/tsSize.depth;
}

template <class T, size_t cn> double Tensor<T,cn>::ComputeSSIM(const Tensor<T,cn>& ts,const Size3& subWinSize, const Size3& subWinStep, int nLevel, int nDir, int boundary_cut, int stsim2_pool_type, int stsim2_modifer) const
{
	double rst=0;

  fstream debugfile;
  if (ts.debugtrigger)
  {
//    string PID = boost::lexical_cast<string>(_getpid());

    debugfile.open("ssim_terms.txt",ios::app);
    debugfile.precision(3);
    //debugfile.width(10);
    debugfile<<"=========inter-band========="<<endl;
    debugfile<<left<<"Band"<<"\t"<<"L"<<"\t"<<"C"<<"\t"<<"C01"<<"\t"<<"C10"<<"\t"<<"Pool"<<endl;
  }

  //20130718
  if (stsim2_modifer == STSIM2_LSE) //do LSE based classifer here, and ignore other part
  {
    //Metric mc;
    int coeffNum=0;
    vector<int> coeff_in_band;
    vector<vector<Tensor<double,2>>> stats = vector<vector<Tensor<double,2>>>(2);
    stats[0] = this->ComputeStatistics(subWinSize,subWinStep);
    stats[1] = ts.ComputeStatistics(subWinSize,subWinStep);
    for (auto t : stats[0])
    {
      coeffNum += t.size().area(); 
      coeff_in_band.push_back( t.size().area() );
    }
    Mat f = Mat(coeff_in_band[0],(nLevel*nDir+2)*4*2,CV_64F);
    Rect roi;
    for (int k=0; k<stats[0].size(); k++)
    {
      roi = Rect(k*2,0,2,coeff_in_band[k]);
      Mat temp=   (stats[0][k] - stats[1][k]).GetFrameRef(0).clone();
      if (k/(nLevel*nDir+2)==0)//means max 255
        temp = temp/255;
      else if (k/(nLevel+nDir+2)==1)//variance max is given by max((mu-min_pix)*(max-pix - mu)) about 128^2
        temp = temp/16384;
       vector<Mat> tt;
      cv::split(temp,tt);
      //abs
      tt[0]=cv::abs(tt[0]);
      tt[1]=cv::abs(tt[1]);
      //satruation
      //cv::threshold(tt[0],tt[0],1,1,cv::THRESH_TRUNC); 
      //cv::threshold(tt[1],tt[1],1,1,cv::THRESH_TRUNC); 
      cv::merge(tt,temp);
      //mylib::DisplayMat(temp);
      //inverse
      if (k/(nLevel*nDir+2)==0||k/(nLevel*nDir+2)==1)
        temp = Scalar(1,1) - temp;
      //mylib::DisplayMat(temp);
      //save to feature vector
      temp.reshape(1,stats[0][k].size().area())
          .copyTo(f(roi));
     /*  (stats[0][k] - stats[1][k]).Abs()
          .GetFrameRef(0).clone()
          .reshape(1,stats[0][k].size().area())
          .copyTo(f(roi));*/
      //mylib::DisplayMat(f(roi));
    }
    Mat A(1,coeffNum*2,CV_64F);
    Mat temp(coeffNum,2,CV_64F);
    int copypos =0;
    for (int j=0; j<coeff_in_band.size();j++)
    {
      Rect roi(j*2,0,2,coeff_in_band[j]);
      Rect toroi(0,copypos,2,coeff_in_band[j]);
      f(roi).copyTo(temp(toroi));
      copypos+=coeff_in_band[j];
    }
    //mylib::DisplayMat(temp);
    stsim2_lse_weight = Mat(A.cols,A.rows,CV_64F,(void*)&STSIM2_LSE_WEIGHT); 
    //mylib::DisplayMat(stsim2_lse_weight,"weight",true);
    //mylib::DisplayMat(A,"A",true);
    temp.reshape(1,1).copyTo(A.row(0));
    Mat scorerst = A*stsim2_lse_weight;
   // mylib::DisplayMat(scorerst);
    return scorerst.at<double>(0,0)/10;
  }
	const double C0 = 0.001;
	const double C1 = 0.001;
	const double C2 = 0.001;
	const double C3 = 0.001;
	CV_DbgAssert(cn==2);
	bool cutBoundary = true;
	Steerable spA, spB;
	//this->ToComplex().Print();
	//ts.ToComplex().Print();
  bool downsample = false;
	spA.buildSCFpyr(this->ToComplex(),nLevel,nDir,1,downsample);//A = this is orgExt
	spB.buildSCFpyr(ts.ToComplex(),nLevel,nDir,1,downsample);//B = ts is candExt
	vector<Tensor<double,2>>& pyrA = spA.getSpaceDomainPyr();
	vector<Tensor<double,2>>& pyrB = spB.getSpaceDomainPyr();
	if ( boundary_cut == FILTER_BOUND_HALF)//modify Dec 27 2011, not cut half but cut half + boundary
	{
		for (unsigned int i = 0; i< pyrA.size(); i++)
		{
			//pyrA[i].Print();
			//dec 27 2011, extend the block by size of block size + bounday (below and right)
			//int hh = subWinSize.height*(int)ceil(float(pyrA[i].size().height)/2.0/float(subWinSize.height));
			//int ww = subWinSize.width*(int)ceil(float(pyrA[i].size().width)/2.0/float(subWinSize.width));
			//int dd = subWinSize.depth*(int)ceil(float(pyrA[i].size().depth)/2.0/float(subWinSize.depth));
			//just half, dec 28 2012
			int hh = pyrA[i].size().height/2;
			int ww = pyrA[i].size().width/2;
			int dd = pyrA[i].size().depth;
			Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, hh,ww,dd);
			//Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, pyrA[i].size().height/2, pyrA[i].size().width/2, pyrA[i].size().depth);
			pyrA[i] = pyrA[i](roi);//no size changed
			//pyrA[i].Print();
			pyrB[i] = pyrB[i](roi);
			//pyrB[i].Print();
		}
	}
  else if ( boundary_cut == FILTER_BOUND_VALID)
  {
    for (unsigned int i = 0; i< pyrA.size(); i++)
		{
  		int hh = pyrA[i].size().height/2 + pyrA[i].size().height/4;
			int ww = pyrA[i].size().width/2  + pyrA[i].size().width/4;
			int dd = pyrA[i].size().depth/2  + pyrA[i].size().depth/4;
      dd<1?dd=1:dd=dd;
			//just half, dec 28 2012
			//int hh = pyrA[i].size().height/2;
			//int ww = pyrA[i].size().width/2;
			//int dd = pyrA[i].size().depth;
			Cube roi(pyrA[i].size().height/8, pyrA[i].size().width/8, 0, hh,ww,dd); // temporaly use this /8
			//Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, pyrA[i].size().height/2, pyrA[i].size().width/2, pyrA[i].size().depth);
			pyrA[i] = pyrA[i](roi);//no size changed
			//pyrA[i].Print();
			pyrB[i] = pyrB[i](roi);
			//pyrB[i].Print();
    }
  }
	vector<Tensor<double,2>> mu_A(pyrA.size());
	vector<Tensor<double,2>> mu_B(pyrB.size());
	vector<Tensor<double,2>> sigma2_A(pyrA.size());
	vector<Tensor<double,2>> sigma2_B(pyrA.size());
	vector<Tensor<double,1>> L(pyrA.size());
	vector<Tensor<double,1>> C(pyrB.size());
	vector<Tensor<double,1>> C01(pyrA.size());
	vector<Tensor<double,1>> C10(pyrB.size());
	int crossbandNum = nDir*(nLevel-1) + nLevel* (mylib::combination(nDir,2));
	vector<Tensor<double,1>> C00 (crossbandNum);
	int index = 0;
	Size3 initSize = pyrA[0].size();
	Size3 sz = pyrA[0].size();
	double weight = sz.height*sz.width;
	cv::Mat temp;
	Size3 subWinStepLv = subWinStep;//include PLC boundary 
	Size3 subWinSizeLv = subWinSize;//include PLC boundary
  //if (boundary_cut == FILTER_BOUND_VALID)
  //{
  //  subWinSizeLv = subWinSizeLv + subWinSizeLv/2; //2*overlap + size
  //  subWinStepLv = subWinStepLv + subWinStepLv/2;
  //}
  //Tensor<double,1> rstMat(pyrA[0].size()/subWinStep - subWinSize+Size3(1,1,1));
	Tensor<double,1> rstMat(Size3(Size3_<double>(pyrA[0].size()-subWinSizeLv)/Size3_<double>(subWinStepLv)) + Size3(1,1,1),0);//gj20130120 set inital value to 1 for productive pooling
                                                                                                                          //set 0 for addictive pooling
  //rstMat.Print();

	for (index = 0; index < (int)pyrA.size(); index++)
	{
		int lvl=0;
    if (downsample)
    {
		  if (index ==pyrA.size()-1) //HP part
		  {
		    lvl = 0;
		  }
		  else
		  {
		    lvl = index/nDir;
		  }
      subWinSizeLv.height = subWinSize.height>>lvl;
		  subWinSizeLv.width = subWinSize.width>>lvl;
		  subWinStepLv.height = subWinStep.height>>lvl;
		  subWinStepLv.width = subWinStep.width>>lvl;
      //if (boundary_cut == FILTER_BOUND_VALID)
      //{
      //  subWinSizeLv = subWinSizeLv + subWinSizeLv/2; //2*overlap + size
      //  subWinStepLv = subWinStepLv + subWinStepLv/2;
      //}
      sz.height = initSize.height>>lvl;
		  sz.width = initSize.width>>lvl;
      //mylib::DisplayMat(gaussKernel);
		  weight = subWinSizeLv.volumn();
    }
		cv::Mat flatKel=cv::Mat(subWinSizeLv,pyrA[index][0].type()-((pyrA[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1.0/weight));
		cv::Mat gaussKernel = mylib::GenGaussKer(subWinSizeLv.height,double(subWinSizeLv.height)/6.0,CV_64F);

		//pyrA[index].Print("pyrA");
   // Tensor<double,1>(gaussKernel).Print("gauss");
		mu_A[index] = pyrA[index].LocalMean(gaussKernel,subWinStepLv);
		//if (ts.debugtrigger)
		//	mu_A[index].Print("mua");
		mu_B[index] = pyrB[index].LocalMean(gaussKernel,subWinStepLv);
		//if (ts.debugtrigger)
		//	mu_B[index].Print("mub");
    //mu_A[index] = pyrA[index].LocalMean(flatKel,subWinStepLv);    //Dec 30 2012, use flat (usuall) average
    //mu_B[index] = pyrB[index].LocalMean(flatKel,subWinStepLv);    //Dec 30 2012, use flat (usuall) average
   /* if (sz.height==32)
    {
      Tensor<double,2>(mu_A[index]).Print();
      Tensor<double,2>(mu_B[index]).Print();
    }*/
    if (stsim2_modifer == STSIM2_BASELINE||stsim2_modifer==STSIM2_TUNE)
      L[index] = (mu_A[index].Abs() * mu_B[index].Abs() *2 + C0).Real() / (mu_A[index]*mu_A[index].Conjugate() + mu_B[index]*mu_B[index].Conjugate() + C0).Real();
    else if (stsim2_modifer == STSIM2_NEW_L1)
    {
     
      double clipThrd = (mu_A[index].Abs().Sum()/double(mu_A[index].size().volumn()))[0]*0.1;
			clipThrd>4?clipThrd=clipThrd:clipThrd=4;
     	clipThrd = clipThrd*clipThrd;
			L[index] = ((mu_A[index]-mu_B[index])*((mu_A[index]-mu_B[index]).Conjugate())).Real().MaxClip(clipThrd);
      L[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>(clipThrd)) - L[index];
	    L[index] = L[index]/clipThrd;	
      if(ts.debugtrigger&&index==4)
      {
        mu_A[index].Print("mu_A");
        mu_B[index].Print("mu_B");
        (mu_A[index]-mu_B[index]).Abs().Real().Print("muA-muB");
        L[index].Print("L");
      }
    }
    else if (stsim2_modifer == STSIM2_NEW_L2)
    {
		//L[index] = (One - (mu_A[index]-mu_B[index]).Abs()/510).Real();
		//L[index] = Tensor<double,1>(mu_A[index].size(),1)- (mu_A[index]-mu_B[index]).Abs().Real()/510;
		//L[index].Print();
		  if (index == 12)//DC band only
		  {
			  double clipThrd = (mu_A[index].Abs().Sum()/double(mu_A[index].size().volumn()))[0]*0.1;
			  clipThrd>4?clipThrd=clipThrd:clipThrd=4;
			  //mu_A[index].Print("mu A");
			  //mu_B[index].Print("mu B");
			  clipThrd = clipThrd*clipThrd;
			  L[index] = ((mu_A[index]-mu_B[index])*((mu_A[index]-mu_B[index]).Conjugate())).Real().MaxClip(clipThrd);
			  //L[index] = mu_A[index].AbsDiff(mu_B[index]).Real().MaxClip(clipThrd); 
			  //L[index].Print();
			  L[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>(clipThrd)) - L[index];
			  //L[index].Print();
			  L[index] = L[index]/clipThrd;
			  //if (ts.debugtrigger)
			  //	L[index].Print();
		  }
      else
			  L[index] = (mu_A[index].Abs() * mu_B[index].Abs() *2 + C0).Real() / (mu_A[index]*mu_A[index].Conjugate() + mu_B[index]*mu_B[index].Conjugate() + C0).Real();
    }
    else if (stsim2_modifer == STSIM2_NEW_L3)
    {
      if (index >=8 && index <= 12)//DC band only and lower band
		  {
			  double clipThrd =(mu_A[index].Abs().Sum()/double(mu_A[index].size().volumn()))[0]*0.1;
			  clipThrd>4?clipThrd=clipThrd:clipThrd=4;
			  //mu_A[index].Print("mu A");
			  //mu_B[index].Print("mu B");
			  clipThrd = clipThrd*clipThrd;
			  L[index] = ((mu_A[index]-mu_B[index])*((mu_A[index]-mu_B[index]).Conjugate())).Real().MaxClip(clipThrd);
			  //L[index] = mu_A[index].AbsDiff(mu_B[index]).Real().MaxClip(clipThrd); 
			  //L[index].Print();
			  L[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>(clipThrd)) - L[index];
			  //L[index].Print();
			  L[index] = L[index]/clipThrd;
			  //if (ts.debugtrigger)
			  //	L[index].Print();
		  }
      else
			  L[index] = (mu_A[index].Abs() * mu_B[index].Abs() *2 + C0).Real() / (mu_A[index]*mu_A[index].Conjugate() + mu_B[index]*mu_B[index].Conjugate() + C0).Real();

    }
	else if (stsim2_modifer == STSIM2_NEW_L4)
	{
		//because the brightness perception is based on the ration, not the differences, so try to use another metric
    double clipMin = mu_A[index].Abs()(0,0,0)[0]*0.3;
    double clipMax = 4;
			//clipThrd>4?clipThrd=clipThrd:clipThrd=4;
     	//clipThrd = clipThrd*clipThrd;
      L[index] = ((mu_A[index]-mu_B[index])/mu_B[index]).Abs().Real().MaxClip(clipMax).MinClip(clipMin);
 
      //L[index].Print();
			//L[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>(clipMax)) - L[index];
      //L[index].Print();
	    //L[index] = L[index]/(clipMax-clipMin);	
      double tempclip = log(clipMax-clipMin);
      for (int ii=0; ii<L[index].size().height;ii++)
        for (int jj=0; jj<L[index].size().width;jj++)
          L[index](ii,jj,0)[0]=log(L[index](ii,jj,0)[0])/tempclip;
     //L[index].Print();
	}
    //L[index].Print();
		//mu_A[index].Print();
		//mu_B[index].Print();
		//(mu_A[index]-mu_B[index]).Print();
		//(mu_A[index]-mu_B[index]).Conjugate().Print();
		//((mu_A[index]-mu_B[index])*((mu_A[index]-mu_B[index]).Conjugate())).Print();
		//pyrA[index].Abs().Real().Print("band A");
		//mu_A[index].Abs().Real().Print("mu_A abs");
		
		//L[index].Print();
		//.Max(Vec<double,2>(clipThrd,0));
		//L[index] = Tensor<double,clipThrd - L[index];
		//max((mu_A[index]-mu_b[inxedx])*(mu_A[index]-mu_b[inxedx])/

		//L[index].Print();
		//here computed the abs(pyrA.* pyrA) that is (actuall) not equal to the true variance
		//for the complex number, it will be a covariance mx instead of var
		//which one is better should be studied.
    //use gaussian kernel (traditional) 
    if (ts.debugtrigger&&index==12)
    {
      pyrA[index].Print("A");
      pyrB[index].Print("B");
    }
		sigma2_A[index] = pyrA[index].LocalVariance(mu_A[index],gaussKernel,subWinStepLv);
		sigma2_B[index] = pyrB[index].LocalVariance(mu_B[index],gaussKernel,subWinStepLv);
    if (ts.debugtrigger&&index==12)
    {
      auto temp2 = (sigma2_A[index].Real()*sigma2_B[index].Real()).Sqrt() + C1;
      for (int ii=0; ii< rstMat.size().height; ii++)
        for (int jj=0; jj<rstMat.size().width; jj++)
        {
          double m = std::numeric_limits<double>::min();
          double temp = (sigma2_A[index](ii,jj,0)[0]*sigma2_B[index](ii,jj,0)[0]+C1);
          cout<<left<<index<<"\t"<<temp2(ii,jj,0)[0]<<endl;
        }  
    }
		//use flat kernel Dec 30 2012
    //sigma2_A[index] = pyrA[index].LocalVariance(mu_A[index],flatKel,subWinStepLv);
		//sigma2_B[index] = pyrB[index].LocalVariance(mu_B[index],flatKel,subWinStepLv);
    //sigma2_B[index].Print();
		//Try ignore the C term in order to work on smooth/smooth or smooth/texture region, gjin aug 24, 2011
		//C[index] = Tensor<double,1>(mu_A[index].size(),Vec<double,1>::all(1));
		C[index] = ((sigma2_A[index].Real()*sigma2_B[index].Real()).Sqrt()*2 + C1)/(sigma2_A[index]+sigma2_B[index] +C1).Real();
    if(ts.debugtrigger&&index==12)
    {
       ((sigma2_A[index].Real()*sigma2_B[index].Real()).Sqrt()*2 + C1).Real().Print("num");   
       (sigma2_A[index]+sigma2_B[index] +C1).Real().Print("denum");
		  C[index].Print();
    }
    C01[index] = ComputeCrossTerm(pyrA[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth)),
									  pyrA[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth)),
									  pyrB[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth)),
									  pyrB[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth)),
									  subWinSizeLv-Size3(0,1,0), subWinStepLv);
		//C01[index].Print();
		C10[index] = ComputeCrossTerm(pyrA[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth)),
									  pyrA[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth)),
									  pyrB[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth)),
									  pyrB[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth)),
									  subWinSizeLv-Size3(1,0,0), subWinStepLv);
		//C10[index].Print();
    
    Tensor<double,1> tempRstMat;
    if (stsim2_modifer == STSIM2_TUNE)
    {
      tempRstMat = L[index]*TUNE_WEIGHT[index*4]+C[index]*TUNE_WEIGHT[index*4+1]+C01[index]*TUNE_WEIGHT[index*4+2]+C10[index]*TUNE_WEIGHT[index*4+3];
    }
    else
      tempRstMat=(L[index]*C[index]*C01[index]*C10[index]).Pow(0.25);//change weight for C gj20130120
	//tempRstMat.Print();
    //if (stsim2_modifer == STSIM2_BASELINE)
  	  //  tempRstMat = tempRstMat.Pow(0.25);// ((L[index]*C[index]*C01[index]*C10[index]).Pow(0.25));
    //else
    //{//    rstMat += ((L[index]*C[index]*C01[index]*C10[index]).Pow(0.25))*()
		//   if (index == 12)
		//	   tempRstMat = tempRstMat.Pow(0.25)*(14.0/3.0);//+= ((L[index]*C[index]*C01[index]*C10[index]).Pow(0.25))*(14.0/3);//(3.5)
		//   else
	  //     tempRstMat = tempRstMat.Pow(0.25)*(28.0/39.0);// += ((L[index]*C[index]*C01[index]*C10[index]).Pow(0.25))*(28.0/39.0);//;*(42.0/52.0);
    //}

    rstMat=rstMat+tempRstMat; //gj20130120 to emphasize every term, use product instead of sum
      /*if (index==12) //lp
        rstMat += tempRstMat*(7);
      else if (index ==13 ) //hp
        rstMat += tempRstMat*(0.875);
      else if (index >=8) //level 1
        rstMat += tempRstMat*(0.875);
      else if (index >=4) //level 2
        rstMat += tempRstMat*(0.4375);
      else //level 3
        rstMat += tempRstMat*(0.2188);
*/
    //}
    //else if (stsim2_modifer == STSIM2_NEW_L2)
   // {
      //
   // }
    //else if (stsim2_modifer == STSIM2_NEW_L3)
    //{
      //
   // }
		//if (ts.debugtrigger)
		//	rstMat.Print();
    //print paramters
    if (ts.debugtrigger)
    {
      for (int ii=0; ii< rstMat.size().height; ii++)
        for (int jj=0; jj<rstMat.size().width; jj++)
          debugfile<<left<<index<<"\t"<<L[index](ii,jj,0)[0]<<"\t"<<C[index](ii,jj,0)[0]<<"\t"<<C01[index](ii,jj,0)[0]<<"\t"<<C10[index](ii,jj,0)[0]<<"\t"<<tempRstMat(ii,jj,0)[0]<<endl;
    }
	}
 
  Tensor<double,1> rstMatBackup;
  if (stsim2_modifer==STSIM2_TUNE)
    rstMatBackup = rstMat;
  else
  {
    //rstMatBackup = (rstMat/double(pyrA.size())).Clone();///;.Pow(1.0/double(index)) gj20130120 use product, so do geometric normalization
    rstMatBackup = rstMat.Clone();
    rstMat = Tensor<double,1>(Size3(Size3_<double>(pyrA[0].size()-subWinSizeLv)/Size3_<double>(subWinStepLv)) + Size3(1,1,1),0);
  }
  if (ts.debugtrigger)
    rstMatBackup.Print();
  if (ts.debugtrigger)
  {
    debugfile<<"STSIM-1 score"<<"\t\t\t\t\t"<<rstMatBackup.Min()[0]/*diabled when using product pooling/index*/<<endl;
    debugfile<<"====== cross band ======="<<endl;
    debugfile<<"Band A"<<"\t"<<"Band B"<<"\t"<<"C00"<<endl;
  }

  if (stsim2_modifer!=STSIM2_TUNE)
  {
	index = 0;
	for (int dr = 0; dr < nDir; dr++)
		for (int lv = 0; lv < nLevel-1; lv++)
		{
			C00[index] = ComputeCrossTerm( pyrA[ lv*nDir + dr].Abs(), pyrA[ (lv+1)*nDir + dr].Abs(),
										   pyrB[ lv*nDir + dr].Abs(), pyrB[ (lv+1)*nDir + dr].Abs(),
										   subWinSizeLv,subWinStepLv);
			//rst+= C00[index].Mean();
			rstMat+= C00[index];
    //  C00[index].Print();
      if(ts.debugtrigger)
      {
        for (int ii=0; ii< rstMat.size().height; ii++)
          for (int jj=0; jj<rstMat.size().width; jj++)
            debugfile<<lv*nDir+dr<<"\t"<<(lv+1)*nDir+dr<<"\t"<<C00[index](ii,jj,00)[0]<<endl;
      }
			index++;
		}
		
	for (int lv=0; lv< nLevel; lv++)
		for (int dr = 0; dr < nDir-1; dr++)
			for ( int p = dr+1; p < nDir; p++)
			{
				C00[index] = ComputeCrossTerm(pyrA[lv*nDir + dr].Abs(),pyrA[lv*nDir + p].Abs(),
											  pyrB[lv*nDir + dr].Abs(),pyrB[lv*nDir + p].Abs(),
											  subWinSizeLv, subWinStepLv);
				//rst+=C00[index].Mean();
				rstMat+= C00[index];
        //C00[index].Print();
        if (ts.debugtrigger)
        {
          for (int ii=0; ii< rstMat.size().height; ii++)
            for (int jj=0; jj<rstMat.size().width; jj++)
              debugfile<<lv*nDir+dr<<"\t"<<(lv)*nDir+p<<"\t"<<C00[index](ii,jj,00)[0]<<endl;
        }
				index++;
			}
  //rstMat = rstMat/double(crossbandNum);
  if (ts.debugtrigger)
    rstMat.Print();
  double totalbands = double(pyrA.size()+crossbandNum);
  
  //gj20120120 use product in STSIM-1 but still use summation in C00, so the pooling is different//
  //rstMat = rstMatBackup*double(pyrA.size())/totalbands + rstMat*double(crossbandNum)/totalbands;

  rstMat = (rstMatBackup + rstMat)/totalbands;
  
  // turn off cross band, so set rstMat
	//rstMat = rstMatBackup; // Dec28 2012, turn off cross band
    
//  rstMat = rstMat/double(pyrA.size());
  //pooling STSIM terms
	//if (rstMat.size().area() == 4)
   // rstMat.Print();
  }
	if (stsim2_pool_type == STSIM2_POOL_AVE)
	{
		rst = rstMat.Mean()[0];
	}
	else
	{
		rst = rstMat.Min()[0];
	}  
  if(ts.debugtrigger)
  {
    debugfile<<"STSIM-2 score"<<"\t\t\t\t\t"<<rst<<endl;
	  //return rst/(crossbandNum+pyrA.size())/sz.depth;
    debugfile.close();
    rstMat.Print();
  }
	return rst;
}

template <class T, size_t cn> double Tensor<T,cn>::ComputeSVMMetric(const Tensor<T,cn>& ts, const Size3& subWinSize, const Size3& subWinStep) const
{

  //SVM s;
  //s.load("svm");
  //Rect roi;
  //build feature
  Metric mc;
  mc.subwinSize = subWinSize;
  mc.subwinStep = subWinStep;
  mc.loadParams();
  mc.loadClassifier("svm");
  double rst= mc.computeMetric(*this,ts);
 //Mat f(1,504,CV_32F);
  return rst;
}


template <class T, size_t cn> vector<Tensor<double,2> > Tensor<T,cn>::ComputeStatistics(const Size3& subWinSize, const Size3& subWinStep, bool subsample, int nLevel, int nDir, bool changeWin) const
{
  typedef vector<Tensor<double,2> > vT;
  Steerable sp;
  sp.buildSCFpyr(this->ToComplex(),3,4,1,subsample);
  vT& pyr = sp.getSpaceDomainPyr();
  int bands = pyr.size();
  Size3 initSize = pyr[0].size();
  Size3 sz = pyr[0].size();
  double weight = subWinSize.height*subWinSize.width;
  double scalar = weight/(weight-1);
	cv::Mat temp;
  int index=0;
	Size3 subWinStepLv = subWinStep;
	Size3 subWinSizeLv = subWinSize;
  int lvl=0;
  vT mu(bands);
  vT sigma2(bands);
  vT rho01(bands);
  vT rho10(bands);
  vT statistics(4*bands); 
  for (index = 0; index < (int)pyr.size(); index++)
	{
    //pyr[index].Print("pry.txt",true); //checked, the pyr are correct
    sz = pyr[index].size();
    if (subWinSizeLv.height>pyr[index].size().height||subWinSizeLv.width>pyr[index].size().width)
    {
      lvl = lvl+1;
      subWinSizeLv = pyr[index].size();
      subWinStepLv = pyr[index].size();
    }
    else if (index==pyr.size()-1) //highpass
    {
      lvl = 0;
      subWinSizeLv = subWinSize;
      subWinStepLv = subWinStep;
    }
    else
      lvl = index/nDir;
    
    if (changeWin) //change window size and not lowpass band
    {
      subWinSizeLv = Size3(subWinSize.height/(1<<lvl),subWinSize.width/(1<<lvl),subWinSize.depth);
      subWinStepLv = Size3(subWinStep.height/(1<<lvl),subWinStep.width/(1<<lvl),subWinStep.depth);
    }
    weight = subWinSizeLv.height*subWinSizeLv.width;
    scalar = weight/(weight-1);

    cv::Mat flatKel1=cv::Mat(subWinSizeLv,pyr[index][0].type()-((pyr[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1.0/weight));
    cv::Mat flatKel2 = cv::Mat(subWinSizeLv,pyr[index][0].type()-((pyr[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1.0/(weight-1)));
    //mylib::DisplayMat(flatKel1);
   // pyr[index].Print("pyr",true);
    mu[index] = pyr[index].LocalMean(flatKel1,subWinStepLv);
    //check mu
   // mu[index].Print();
    sigma2[index] = pyr[index].LocalVariance(mu[index]*scalar,flatKel2,subWinStepLv).Sqrt();//take the sqrt to get std
    //is rho01 rho10 a good choice?????? YES 20130722
    rho01[index] = ComputeRho(pyr[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth)),
									  pyr[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth)),
									  subWinSizeLv-Size3(0,1,0), subWinStepLv);
    rho10[index] = ComputeRho(pyr[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth)),
									  pyr[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth)),
									  subWinSizeLv-Size3(1,0,0), subWinStepLv);
    statistics[index]=mu[index].Clone();
    statistics[index+bands]=sigma2[index].Clone();
    statistics[index+2*bands] = rho01[index].Clone();
    statistics[index+3*bands] = rho10[index].Clone();
  }
  return statistics;
}
template <class T, size_t cn> double Tensor<T,cn>::ComputeMahalanobis(const Tensor<T,cn>& ts, Size3 subWinSize, Size3 subWinStep, const Mat& iMcovar) const
{
	double rst=0;
  const double C = 0.001;
	const double C0 = 0.001;
	const double C1 = 0.001;
	const double C2 = 0.001;
	const double C3 = 0.001;
	//CV_DbgAssert(cn==2);
	bool cutBoundary = true;
	Steerable spA, spB;
  int nDir = 4; int nLevel = 3; 
	int crossbandNum = nDir*(nLevel-1) + nLevel* (mylib::combination(nDir,2));
	spA.buildSCFpyr(this->ToComplex(),nLevel,nDir,1,false);
	spB.buildSCFpyr(ts.ToComplex(),nLevel,nDir,1,false);
	vector<Tensor<double,2>>& pyrA = spA.getSpaceDomainPyr();
	vector<Tensor<double,2>>& pyrB = spB.getSpaceDomainPyr();
  int boundary_cut = FILTER_BOUND_HALF;
	if ( boundary_cut == FILTER_BOUND_HALF)//modify Dec 27 2011, not cut half but cut half + boundary
	{
		for (unsigned int i = 0; i< pyrA.size(); i++)
		{
			//pyrA[i].Print();
			//dec 27 2011, extend the block by size of block size + bounday (below and right)
			int hh = subWinSize.height*(int)ceil(float(pyrA[i].size().height)/2.0/float(subWinSize.height));
			int ww = subWinSize.width*(int)ceil(float(pyrA[i].size().width)/2.0/float(subWinSize.width));
			int dd = subWinSize.depth*(int)ceil(float(pyrA[i].size().depth)/2.0/float(subWinSize.depth));
			Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, hh,ww,dd);
			//Cube roi(pyrA[i].size().height/4, pyrA[i].size().width/4, 0, pyrA[i].size().height/2, pyrA[i].size().width/2, pyrA[i].size().depth);
			pyrA[i] = pyrA[i](roi);//no size changed
			//pyrA[i].Print();
			pyrB[i] = pyrB[i](roi);
			//pyrB[i].Print();
		}
	}
	Size3 sz = pyrA[0].size();
	//Tensor<double,1> rstMat(pyrA[0].size()/subWinStep - subWinSize+Size3(1,1,1));
	Tensor<double,1> rstMat(Size3(Size3_<double>(pyrA[0].size()-subWinSize)/Size3_<double>(subWinStep)) + Size3(1,1,1));
	//rstMat.Print();
	vector<Tensor<double,2>> mu_A(pyrA.size());
	vector<Tensor<double,2>> mu_B(pyrB.size());
  Tensor<double,2> mu11,mu12,mu21,mu22;
	vector<Tensor<double,2>> sigma2_A(pyrA.size());
	vector<Tensor<double,2>> sigma2_B(pyrA.size());
  vector<Tensor<double,2>> rho01_A(pyrA.size());
  vector<Tensor<double,2>> rho01_B(pyrA.size());
  vector<Tensor<double,2>> rho10_A(pyrA.size());
  vector<Tensor<double,2>> rho10_B(pyrA.size());
  Mat featureA(pyrA.size()*4,1,CV_64F);
  Mat featureB(pyrA.size()*4,1,CV_64F);
  Tensor<double,2> sigma1_cross,sigma2_cross;
	Tensor<double,2> sigma11,sigma12,sigma21,sigma22;

	vector<Tensor<double,1>> C00 (crossbandNum);
	int index = 0;

	double weight = sz.height*sz.width;
	cv::Mat flatKel=cv::Mat(sz,pyrA[index][0].type()-((pyrA[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1/weight));
	cv::Mat gaussKernel = mylib::GenGaussKer(subWinSize.height,double(subWinSize.height)/6.0,CV_64F);
	cv::Mat temp;
	for (index = 0; index < (int)pyrA.size(); index++)
	{
		mu_A[index] = pyrA[index].LocalMean(gaussKernel,subWinStep);
		mu_B[index] = pyrB[index].LocalMean(gaussKernel,subWinStep);
		sigma2_A[index] = pyrA[index].LocalVariance(mu_A[index],gaussKernel,subWinStep);
		sigma2_B[index] = pyrB[index].LocalVariance(mu_B[index],gaussKernel,subWinStep);
    Tensor<double,2>& im11 = pyrA[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth));
    Tensor<double,2>& im12 = pyrA[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth));
    Tensor<double,2>& im21 = pyrB[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth));
		Tensor<double,2>& im22 = pyrB[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth)),
	  flatKel = Mat(Size(subWinSize.width-1,subWinSize.height),CV_64F,Scalar(1/double((subWinSize.width-1)*subWinSize.height)));
 		mu11 = im11.LocalMean(subWinSize,subWinStep);         //mu11.Print();
	  mu12 = im12.LocalMean(subWinSize,subWinStep);
	  mu21 = im21.LocalMean(subWinSize,subWinStep);
	  mu22 = im22.LocalMean(subWinSize,subWinStep);
	  sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
	  sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep); 
	  sigma21 = im21.LocalVariance(mu21,subWinSize,subWinStep); 
	  sigma22 = im22.LocalVariance(mu22,subWinSize,subWinStep); 
	  sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
	  sigma2_cross = ((im21 * im22.Conjugate()).LocalMean(subWinSize,subWinStep) - mu21 * mu22.Conjugate());// sigma2_cross.Print();
	  rho01_A[index] = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
	  rho01_B[index] = (sigma2_cross + Vec2d(C,0))/((sigma21*sigma22).Sqrt()+Vec2d(C,0));// rho2.Print();
    im11 = pyrA[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth));    
    im12 = pyrA[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth));
    im21 = pyrB[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth));
    im22 = pyrB[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth));
    flatKel = Mat(Size(subWinSize.width,subWinSize.height-1),CV_64F,Scalar(1/double(subWinSize.width*(subWinSize.height-1))));
    mu11 = im11.LocalMean(subWinSize,subWinStep);         //mu11.Print();
	  mu12 = im12.LocalMean(subWinSize,subWinStep);
	  mu21 = im21.LocalMean(subWinSize,subWinStep);
	  mu22 = im22.LocalMean(subWinSize,subWinStep);
	  sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
    sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep); 
    sigma21 = im21.LocalVariance(mu21,subWinSize,subWinStep); 
	  sigma22 = im22.LocalVariance(mu22,subWinSize,subWinStep); 
 //    sigma11.Print("sigma11");
 //   sigma12.Print("sigma12");
 //   sigma1_cross.Print("sigma1_cross");
    sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
    sigma2_cross = ((im21 * im22.Conjugate()).LocalMean(subWinSize,subWinStep) - mu21 * mu22.Conjugate());// sigma2_cross.Print();
	  rho10_A[index] = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
    rho10_B[index] = (sigma2_cross + Vec2d(C,0))/((sigma21*sigma22).Sqrt()+Vec2d(C,0));// rho2.Print();

   // rho10_A[index].Print();
   // rho10_A[index].Abs().Real().Print();
   // cout<<rho10_A[index].Abs().Real()(0,0,0)[0]<<endl;;
    featureA.at<double>(index*4,0) = mu_A[index].Abs().Real()(0,0,0)[0];
    featureA.at<double>(index*4+1,0) = sqrt(sigma2_A[index].Abs().Real()(0,0,0)[0]);
    featureA.at<double>(index*4+2,0) = rho01_A[index].Abs().Real()(0,0,0)[0];
    featureA.at<double>(index*4+3,0) = rho10_A[index].Abs().Real()(0,0,0)[0];
		featureB.at<double>(index*4,0) = mu_B[index].Abs().Real()(0,0,0)[0];
    featureB.at<double>(index*4+1,0) = sqrt(sigma2_B[index].Abs().Real()(0,0,0)[0]);
    featureB.at<double>(index*4+2,0) = rho01_B[index].Abs().Real()(0,0,0)[0];
    featureB.at<double>(index*4+3,0) = rho10_B[index].Abs().Real()(0,0,0)[0];
	}
	//return rst/(crossbandNum+pyrA.size())/sz.depth;
 // mylib::DisplayMat(featureA,"fa");
 // mylib::DisplayMat(featureB,"fb");
 // mylib::DisplayMat((featureA-featureB).t());
 // mylib::DisplayMat(iMcovar,"imcovar");
/*  Mat temp2(1,iMcovar.size().width,CV_64F,Scalar(0));
 // mylib::DisplayMat(temp2,"temp2");
  cv::gemm(featureA-featureB,iMcovar,1,featureA.t(),0,temp2,GEMM_1_T);
  Mat F = (featureA - featureB);
  mylib::DisplayMat(F,"F");
  rst = 0;
  for (int jj=0; jj<iMcovar.size().width; jj++)
  {
    temp2.at<double>(0,jj)=0;
    for (int ii=0; ii<iMcovar.size().height; ii++)
    {
      temp2.at<double>(0,jj) += F.at<double>(ii,0)*iMcovar.at<double>(ii,jj);
      cout<<"+"<<F.at<double>(ii,0)<<" x "<<iMcovar.at<double>(ii,jj)<<" = "<< temp2.at<double>(0,jj)<<endl;
     
    }
  }
  mylib::DisplayMat(temp2,"temp2");
  mylib::DisplayMat((featureA-featureB).t()*iMcovar);
  temp2 = (featureA - featureB).t()*iMcovar*(featureA-featureB);
  rst = temp2.at<double>(0,0);*/
  return cv::Mahalanobis(featureA,featureB,iMcovar);
}

template<class T, size_t cn> cv::Mat Tensor<T,cn>::EstimateVarForMahalanobis(Size3 wsize, Size3 stepsize)
{
  Mat covar;
  Mat mean;
  //compute 2D only
  CV_DbgAssert(wsize.depth==1);
  int numberofparam = 4*14;
  int numberofvector = (this->tsSize.height)/stepsize.height*(this->tsSize.width)/stepsize.width;
  //CV_MAKE_TYPE(CV_64F,cn);
  //CV_MAKETYPE(
  cv::Mat samples(numberofparam,numberofvector,CV_64F);
  Cube roi(Point3i(0,0,0),wsize);
  Tensor<T,cn> block;
  int nLevel = 3;
  int nDir = 4;
  int nband = nLevel*nDir +2;
  vector<Tensor<double,2>> mu_A(nband);
	vector<Tensor<double,2>> sigma2_A(nband);
  vector<Tensor<double,2>> rho01(nband), rho10(nband);

  int crossbandNum = nDir*(nLevel-1) + nLevel* (mylib::combination(nDir,2));
	vector<Tensor<double,1>> C00 (crossbandNum);
  int index = 0;
  Size3 subWinSize = wsize;
  Size3 SetSubWinStep = stepsize;
	cv::Mat flatKel;//=cv::Mat(sz,CV_64F2[0].type()-((pyrA[index][0].channels()-1)<<CV_CN_SHIFT),Scalar(1/weight));
	cv::Mat gaussKernel = mylib::GenGaussKer(wsize.height,double(wsize.height)/6.0,CV_64F);
	cv::Mat temp;
  double C = 0.001;
  int col =0;
  for (int i=0; i<this->tsSize.height; i+=wsize.height)
  {
    for (int j=0; j<this->tsSize.width; j+=wsize.width)
    {
      block = this->GetBlock(roi);
      //Tensor<double,1>(block).Display();
      Steerable spA;
    	spA.buildSCFpyr(Tensor<double,1>(block).ToComplex(),3,4,1,false);
      vector<Tensor<double,2>>& pyrA = spA.getSpaceDomainPyr();
      for (int index=0; index<pyrA.size(); index++)
      {
       // pyrA[index].Print();
        mu_A[index] = pyrA[index].LocalMean(gaussKernel,subWinStep);
        sigma2_A[index] = pyrA[index].LocalVariance(mu_A[index],gaussKernel,subWinStep);
        Size3 sz = pyrA[index].size();
        Tensor<double,2>& im11 = pyrA[index](Cube(0,0,0,sz.height,sz.width-1,sz.depth));
        Tensor<double,2>& im12 = pyrA[index](Cube(0,1,0,sz.height,sz.width-1,sz.depth));
	      flatKel = Mat(Size(wsize.width-1,wsize.height),CV_64F,Scalar(1/double((wsize.width-1)*wsize.height)));
      	Tensor<double,2> mu11 = im11.LocalMean(subWinSize,subWinStep);         //mu11.Print();
      	Tensor<double,2> mu12 = im12.LocalMean(subWinSize,subWinStep);
      	Tensor<double,2> sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
      	Tensor<double,2> sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep); 
      	Tensor<double,2> sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
	      rho01[index] = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
		    im11 = pyrA[index](Cube(0,0,0,sz.height-1,sz.width,sz.depth));    
        im12 = pyrA[index](Cube(1,0,0,sz.height-1,sz.width,sz.depth));
        flatKel = Mat(Size(subWinSize.width,subWinSize.height-1),CV_64F,Scalar(1/double(subWinSize.width*(subWinSize.height-1))));
	      sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
      	sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep); 
      	sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
	      rho10[index] = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
        //cout<< mu_A[index].Abs().Real()(0,0,0)[0]<<endl<<sigma2_A[index].Abs().Real()(0,0,0)[0]<<endl<< rho01[index].Abs().Real()(0,0,0)[0] <<endl<< rho10[index].Abs().Real()(0,0,0)[0]<<endl;
        samples.at<double>(index*4,col) = mu_A[index].Abs().Real()(0,0,0)[0];
        samples.at<double>(index*4+1, col) = sqrt(sigma2_A[index].Abs().Real()(0,0,0)[0]);
        samples.at<double>(index*4+2,col) = rho01[index].Abs().Real()(0,0,0)[0];
        samples.at<double>(index*4+3,col) = rho10[index].Abs().Real()(0,0,0)[0];
      }
      cout<<roi.x<<"," <<roi.y<<endl;
      roi.y +=wsize.width;
      col++;
    }
    roi.x += wsize.height;
    roi.y =0;
    //col++;
  }
  //mylib::DisplayMat(samples, "samples",true);
  cv::calcCovarMatrix(samples,covar,mean,CV_COVAR_NORMAL|CV_COVAR_COLS);
 // mylib::DisplayMat(covar,"covar",true);
  //mylib::DisplayMat(mean,"mean",true);
  cv::Mat icovar;
  cv::invert(covar,icovar,DECOMP_SVD);
  //Tensor<double,1>(icovar).Print("icovar",true);
  return icovar;
}

template<class T, size_t cn> double Tensor<T,cn>::Compare(const Tensor<T,cn>& ts, int criteria,  Printable_t param1, Printable_t param2, Printable_t param3, Printable_t param4, Printable_t param5, Printable_t param6) const // double param1, double param2, int param3,int param4, int param5, cv::Mat& param6) const
{
  /*va_list ap;
  va_start(ap,criteria);
  vector<string> argname;
  vector<Printable_t> argdata;
  Printable_t printable;
  int count=0;
  for (count=0;count<paramNum;count++)
  {
    char* s = va_arg(ap,char*);
    printable.s = s;// va_arg(ap,char*);
    argname.push_back(printable.s);
      //cout<<Printable.s<<endl;;
    if (criteria==COMPARE_CRITERIA_SSIM) //(!strcmp(printable.s,"nDir")||!strcmp(printable.s,"nLevel")||!strcmp(printable.s,"pool_type")||!strcmp(printable.s,"boundary_cut"))
        printable.i = va_arg(ap,int);
    else if (criteria==COMPARE_CRITERIA_MAHALANOBIS)
        printable.m = va_arg(ap,cv::Mat*);
    else
        printable.d = va_arg(ap,double);
    argdata.push_back(printable);
      //cout<<Printable.i<<endl;
  }
  va_end(ap);
  */
	if (criteria == COMPARE_CRITERIA_MSE)
  { 
    int modifier = param1.i;//get the modifer
    int bsize = param2.i;
    int osize = param3.i;
    if (modifier == SE_MSE)
    {
      Tensor<T,cn> A_low,A_right,B_low,B_right;
      A_low = this->GetBlock(Cube(osize+bsize,osize,0,osize,bsize,1));
      A_right = this->GetBlock(Cube(osize,osize+bsize,0,bsize,osize,1));
      B_low = ts.GetBlock(Cube(osize+bsize,osize,0,osize,bsize,1));
      B_right=ts.GetBlock(Cube(osize,osize+bsize,0,bsize,osize,1));
      return A_low.ComputeMSE(B_low)+A_right.ComputeMSE(B_right); 
    }
		return ComputeMSE(ts);
  }
	else if (criteria == COMPARE_CRITERIA_INTERP)
	{
    //CV_DbgAssert(paramNum==1);
    
		bool rst = ComputeAIM(ts,param1.d);

		if (rst)
			return 1.0;
		else
			return 0;
	}
	else if (criteria == COMPARE_CRITERIA_SSIM)//|STSIM2_BASELINE|SSIM2_DCT|STSIM2_ADT_NATIVE|STSIM2_NEW_L1|STSIM2_NEW_L2|STSIM2_NEW_L3|STSIM2_NEW_L4|STSIM2_SE_MSE)
	{
    //CV_DbgAssert(paramNum==5);
		if (subWinSize == Size3(0,0,0))
			CV_Error(CV_StsBadSize,"subWinSize is un-set");
    int nDir, nLevel, bd_cut, pool_type, modifier;
    nDir = param2.i;
    nLevel = param1.i;
    bd_cut = param3.i;
    pool_type = param4.i;
    modifier = param5.i;
    /*
    count=0;
    for (string name : argname)
    {
      if (!name.compare("nDir"))
      {
        nDir = argdata[count].i;
      }
      else if (!name.compare("nLevel"))
      {
        nLevel = argdata[count].i;
      }
      else if (!name.compare("boundary_cut"))
      {
        bd_cut=argdata[count].i;
      }
      else if (!name.compare("pool_type"))
      {
        pool_type = argdata[count].i;
      }
      else if (!name.compare("metric_modifier"))
      {
        modifier = argdata[count].i;
      }
      else
      {
        CV_DbgAssert(0);
      }
      count++;
    }
    */
		Tensor<T,2> T1 = this->ToComplex();
    Size3 tempSubSize,tempSubStep;
    if (T1.size().height<=subWinSize.height&&bd_cut==FILTER_BOUND_HALF)
    {
      tempSubSize = T1.size()/2+Size3(0,0,1);
      tempSubStep = T1.size()/2+Size3(0,0,1);
    }
    else
    {
      tempSubSize = subWinSize;
      if (modifier==STSIM2_LSE)
        tempSubStep = tempSubSize;
      else
        tempSubStep = Size3(subWinSize.height/4,subWinSize.width/4,1);/*subWinStep;*///20130521 use 1/4 of subwinSize as the step size in order to do blk + LU boundary
    }

	  T1.SetSubWinSize(tempSubSize);
	  T1.SetSubWinStep(tempSubStep);

		Tensor<T,2> T2 = ts.ToComplex(); 
		T2.debugtrigger = ts.debugtrigger;
    
		//if (param3) //gj01142013
			return T1.ComputeSSIM(T2,tempSubSize,tempSubStep,nLevel,nDir,/*FILTER_BOUND_HALF*/bd_cut, pool_type, modifier);
		//else
		//	return T1.ComputeSSIM(T2,subWinSize,subWinStep,(int)param1,(int)param2, FILTER_BOUND_VALID, param4, param5);
	}
	else if (criteria == COMPARE_CRITERIA_SUBJECTIVE)
	{
		//TBD!
		return 1.0;
	}
	else if (criteria== COMPARE_CRITERIA_SAD)
	{
		return ComputeSAD(ts);
	}
  else if (criteria== COMPARE_CRITERIA_MAHALANOBIS)
  {
    //CV_DbgAssert(count==2);
    return ComputeMahalanobis(ts,subWinSize,subWinStep,*param1.m);
  }
  else if (criteria== COMPARE_CRITERIA_LRI)
  {
    return ComputeLRI(ts);
  }
  else if (criteria == COMPARE_CRITERIA_SVM)
  {
    return ComputeSVMMetric(ts,subWinSize,subWinStep);
  }
	else
		return INT_MAX;

}


template<class T, size_t cn>
Tensor<double,1> Tensor<T,cn>::ComputeRho(const Tensor<double,2>& im11, const Tensor<double,2>& im12, const Size3& subWinSize, const Size3& subWinStep) const
{

	if (im11.size()!= im12.size()) 
	{
		CV_Error(CV_StsUnmatchedSizes,"im11,im12 must in same size");
	}
  const double C = 0.0001;
	//im11.Print();
	//im12.Print();
	//im21.Print();
	//im22.Print();
	Tensor<double,2> mu11,mu12;
	Tensor<double,2> sigma_cross;
	Tensor<double,2> rho;
	Tensor<double,2> sigma11,sigma12;
	Mat flatKel = Mat(Size(subWinSize.width,subWinSize.height),CV_64F,Scalar(1/double(subWinSize.width*subWinSize.height)));
	mu11 = im11.LocalMean(subWinSize,subWinStep);      //   mu11.Print();
	mu12 = im12.LocalMean(subWinSize,subWinStep);//  mu12.Print();
	sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
	sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep);// sigma12.Print();
	sigma_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
	rho = (sigma_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
	return rho;
}

template<class T, size_t cn>
Tensor<double,1> Tensor<T,cn>::ComputeCrossTerm(const Tensor<double,2>& im11, const Tensor<double,2>& im12,Tensor<double,2>& im21, const Tensor<double,2>& im22, const Size3& subWinSize, const Size3& subWinStep) const
{

	if (im11.size()!= im12.size() ||
		im11.size()!= im21.size() ||
		im11.size()!= im22.size())
	{
		CV_Error(CV_StsUnmatchedSizes,"im11,im12,im21,im22 must in same size");
	}
	//im11.Print();
	//im12.Print();
	//im21.Print();
	//im22.Print();
	const double C = 0.001;
	Tensor<double,2> mu11,mu12,mu21,mu22;
	Tensor<double,2> sigma1_cross,sigma2_cross;
	Tensor<double,2> rho1,rho2;
	Tensor<double,2> sigma11,sigma12,sigma21,sigma22;
	Tensor<double,1> rst; 
	Mat flatKel;
	//if (im11.size().height> im11.size().width)
	//	flatKel = Mat(Size(subWinSize.width-1,subWinSize.height),CV_64F,Scalar(1/double((subWinSize.width-1)*subWinSize.height)));
	//else if (im11.size().height < im11.size().width)
	//	flatKel = Mat(Size(subWinSize.width,subWinSize.height-1),CV_64F,Scalar(1/double(subWinSize.width*(subWinSize.height-1))));
	//else
		flatKel = Mat(Size(subWinSize.width,subWinSize.height),CV_64F,Scalar(1/double(subWinSize.width*subWinSize.height)));
	//mu11 = im11.LocalMean(flatKel,subWinStep); mu11.Print();
	mu11 = im11.LocalMean(subWinSize,subWinStep);      //   mu11.Print();
	//mu12 = im12.LocalMean(flatKel,subWinStep); //mu12.Print();
	//mu21 = im21.LocalMean(flatKel,subWinStep); //mu21.Print();
	//mu22 = im22.LocalMean(flatKel,subWinStep); //mu22.Print();
	mu12 = im12.LocalMean(subWinSize,subWinStep);//  mu12.Print();
	mu21 = im21.LocalMean(subWinSize,subWinStep);//  mu21.Print();
	mu22 = im22.LocalMean(subWinSize,subWinStep);//  mu22.Print();
	sigma11 = im11.LocalVariance(mu11,subWinSize,subWinStep); //sigma11.Print();
	//sigma11 = im11.LocalVariance(mu11,flatKel,subWinStep); sigma11.Print();
	//sigma12 = im12.LocalVariance(mu12,flatKel,subWinStep); //sigma12.Print();
	//sigma21 = im21.LocalVariance(mu21,flatKel,subWinStep); //sigma21.Print();
	//sigma22 = im22.LocalVariance(mu22,flatKel,subWinStep); //sigma22.Print();
	sigma12 = im12.LocalVariance(mu12,subWinSize,subWinStep);// sigma12.Print();
	sigma21 = im21.LocalVariance(mu21,subWinSize,subWinStep);// sigma21.Print();
	sigma22 = im22.LocalVariance(mu22,subWinSize,subWinStep);/// sigma22.Print();

	//sigma1_cross = ((im11 * im12.Conjugate()).Filter2D(flatKel).SubSample(subWinStep) - mu11 * mu12.Conjugate()); sigma1_cross.Print();
	sigma1_cross = ((im11 * im12.Conjugate()).LocalMean(subWinSize,subWinStep) - mu11 * mu12.Conjugate()); //sigma1_cross.Print();
	sigma2_cross = ((im21 * im22.Conjugate()).LocalMean(subWinSize,subWinStep) - mu21 * mu22.Conjugate()); //sigma2_cross.Print();
	rho1 = (sigma1_cross + Vec2d(C,0))/((sigma11*sigma12).Sqrt()+Vec2d(C,0));// rho1.Print();
	//(sigma11*sigma12).Sqrt().Print();
	rho2 = (sigma2_cross + Vec2d(C,0))/((sigma21*sigma22).Sqrt()+Vec2d(C,0)); //rho2.Print();
	rst = Tensor<double,1>(rho1.size(),Vec<double,1>::all(1)) - (rho1-rho2).Abs().Real()/2;
  //rst.Print();
	return rst;
}

template<class T, size_t cn> 
double Tensor<T,cn>::ComputeCrossTerm(const cv::Mat &im11,const cv::Mat &im12,const cv::Mat &im21,const cv::Mat &im22) const
{
	if (im11.size()!= im12.size() ||
		im11.size()!= im21.size() ||
		im11.size()!= im22.size())
	{
		CV_Error(CV_StsUnmatchedSizes,"im11,im12,im21,im22 must in same size");
	}

	if (im11.channels() != im12.channels())
	{
		CV_Error(CV_StsUnmatchedSizes,"channels are not same");
	}

	const double C = 0.001;
	const double size = double(im11.size().height * im11.size().width);

	complex<double> mu11,mu12,mu21,mu22,sigma1_cross,sigma2_cross,rho1,rho2;
	double sigma11,sigma12,sigma21,sigma22,rst; 
	cv::Mat flatKel(im11.size(),im11.type()-((im11.channels()-1)<<CV_CN_SHIFT),Scalar(1/size));
	int a = flatKel.channels();
	int b = im11.channels();
	//	if (im11.channels() == 2 || im11.channels() == 6)
	//		{
	mu11= mylib::WeightedSumComplex(im11,flatKel);
	mu12= mylib::WeightedSumComplex(im12,flatKel);
	mu21= mylib::WeightedSumComplex(im21,flatKel);
	mu22= mylib::WeightedSumComplex(im22,flatKel);
	sigma11 = abs(mylib::WeightedSumComplex(mylib::complexMul(im11,im11),flatKel) - mu11*conj(mu11));
	sigma12 = abs(mylib::WeightedSumComplex(mylib::complexMul(im12,im12),flatKel) - mu12*conj(mu12));
	sigma21 = abs(mylib::WeightedSumComplex(mylib::complexMul(im21,im21),flatKel) - mu21*conj(mu21));
	sigma22 = abs(mylib::WeightedSumComplex(mylib::complexMul(im22,im22),flatKel) - mu22*conj(mu22));
	sigma1_cross = mylib::WeightedSumComplex(mylib::complexMul(im11,im12),flatKel) - mu11*conj(mu12);
	sigma2_cross = mylib::WeightedSumComplex(mylib::complexMul(im21,im22),flatKel) - mu21*conj(mu22);
	rho1 = (sigma1_cross + C)/(sqrt(sigma11*sigma12)+C);
	rho2 = (sigma2_cross + C)/(sqrt(sigma21*sigma22)+C);
	rst = 1 - 0.5*abs(rho1-rho2);
	return rst;
}

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
		cv::Laplacian(this->GetFrame(z),rst.GetFrame(z),this->depth(),1,1,0,BORDER_REFLECT_101);
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
			cv::filter2D(this->GetFrame(z),rst.GetFrame(z),-1,ker);
		}
		if (boundary == FILTER_BOUND_EXTEND)
		{
			CV_Error(CV_StsNotImplemented,"unsupport extend boundary version");
			return rst;
		}
		else if (boundary == FILTER_BOUND_VALID)
		{
			Size3 validSize(tsSize.height - ker.size().height +1, tsSize.width - ker.size().width+1, 1);
			Point3i validOffset(ker.size().height/2, ker.size().width/2,0);
			return rst.Crop(validOffset,validSize);
		}
		else
			return rst;
	
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::SubSample(const Size3& subrate)
{
	Size3 rstSize = (Size3_<double>(tsSize)/Size3_<double>(subrate)).Ceil();
	
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
		cv::minMaxLoc(this->GetFrame(z),&temp);
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
		rst.SetFrame(i,mylib::FreqComplexFilter(GetFrame(i),kernel,conj));
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::DFT(void) const
{

	Tensor<T,cn> rst(size());
	CV_Assert(cn ==2 || cn== 6); 
	for (int i=0; i< size().depth; i++)
	{
		cv::dft(GetFrame(i),rst[i]);
	}
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::IDFT(void) const
{
	Tensor<T,cn> rst(size());
	CV_Assert(cn ==2 || cn== 6); 
	for (int i=0; i< size().depth; i++)
	{
		cv::dft(GetFrame(i),rst[i],CV_DXT_INVERSE|DFT_SCALE);
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
		rst.SetFrame(i,DFTShift(GetFrame(i)));
	}
	return rst;
}


template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Normalize(void) const
{
	CV_Assert(cn == 1);
	Tensor<T,cn> rst(size());
	for (int i=0; i< size().depth; i++)
	{
		cv::normalize(GetFrame(i),rst[i],0,255,NORM_MINMAX);
	}
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LSFitting(int order) const
{
	if ( order > 2 || order < 1)
	{
		CV_Error(CV_StsNotImplemented,"unsupport");
	}
	if (order == 1)
	{
		int count = 0;
		Tensor<T,cn> H(Size3(tsSize.volumn(),3,1));
		Tensor<T,cn> t(Size3(tsSize.volumn(),1,1));
		Tensor<T,cn> rst(Size3(3,1,1));
		for (int i=0; i< tsSize.height; i++)
			for (int j=0; j < tsSize.width; j++)
			{
				H(count,0,0) = saturate_cast<T>(i);
				H(count,1,0) = saturate_cast<T>(j);
				H(count,2,0) = 1;
				t(count,0,0) = this->operator()(i,j,0);
				count++;
			}
			for (int i=0; i< cn; i++)
			{
				rst[i] = (H[i].t() * H[i]).inv()*(H[i].t()*t[i]);
			}
			return rst;
	}
	else
	{
		Tensor<T,cn> rst(Size3(6,1,1));
		Tensor<T,cn> H(Size3(tsSize.volumn(),6,1));
		Tensor<T,cn> t(Size3(tsSize.volumn(),1,1));
		int count =0;
		for (int i=0; i< tsSize.height; i++)
			for (int j=0; j< tsSize.width; j++)
			{
				H(count,0,0)=saturate_cast<T>(i*i);
				H(count,1,0)=saturate_cast<T>(i*j);
				H(count,2,0)=saturate_cast<T>(j*j);
				H(count,3,0)=saturate_cast<T>(i);
				H(count,4,0)=saturate_cast<T>(j);
				H(count,5,0)=1;
				t(count,0,0)=this->operator()(i,j,0);
				count++;
			}
			for (int i=0; i< cn; i++)
			{
				rst[i] = (H[i].t() * H[i]).inv()*(H[i].t()*t[i]);
			}
			return rst; 
	}
}

template<class T, size_t cn>
Tensor<T,cn> Tensor<T,cn>::BuildLightingPlane(const Tensor<T,cn>& param, int order) const
{
	if ( order > 2 || order < 1)
	{
		CV_Error(CV_StsNotImplemented,"unsupport");
	}
	Tensor<T,cn> rst(tsSize);
	//param.Print();
	for (int k=0; k<tsSize.depth; k++)
		for (int i=0; i<tsSize.height; i++)
			for (int j=0; j< tsSize.width; j++)
			{
				if(order==1)
					rst(i,j,k)=(param(0,0,k)*i+param(1,0,k)*j+param(2,0,k))(0);
				else
					rst(i,j,k)=(param(0,0,k)*i*i+param(1,0,k)*i*j+param(2,0,k)*j*j+param(3,0,k)*i+param(4,0,k)*j+param(5,0,k))(0);
			}
			return rst;
}

template <class T, size_t cn>
Tensor<T,cn> Tensor<T,cn>::LightingCorrection(const Tensor<T,cn>& changeTo,bool saveCodeLength)
{
	//do this after quilting
	//change to is the original block (the lighting should be)
	//Tensor<T,cn> tempTagLighting = changeTo.LSFitting();
	//tempTagLighting.Print();
	//Tensor<T,cn> tempCanLighting = LSFitting();
	//tempCanLighting.Print();
	//for (int i=0; i< tempTagLighting.size().height; i++)
	//	lightTag.push_back(tempTagLighting(i,0,0));
	//for (int j=0; j< tempCanLighting.size().height; j++)
	//	lightCan.push_back(tempCanLighting(j,0,0));
	//// add DCT coding light plane
	//Tensor<T,cn> lighting = BuildLightingPlane(tempTagLighting-tempCanLighting);
	Tensor<T,cn> lighting  = changeTo;
	//lighting.Print();
	///quantize the ligting surface by dct 
	//refresh lighting array
	//lightingDCTCoeffStat.clear();
	int lf = 0;//number of LF 
	const float* tbl;
	if (size().height<=16)
	{
		tbl = base8;
		lf = 64;
	}
	else if (size().height<=32)
	{
		tbl = base16;
		lf = 6;
	}
	else
	{
		tbl = base32;
		lf = 4;
	}
	int index = 0;
	int tempLength =0;
	codeLength = 0;
	float qfactor = 5;
	Tensor<double,cn> flit = Tensor<double,cn>(changeTo); //target X
	Tensor<double,cn> fcand = Tensor<double,cn>(*this); //candidate Y
	//flit.Print();
	//flit = flit - 128;
	flit = flit - fcand;

	for (int i=0; i<flit.size().depth; i++)
	{
		//mylib::DisplayMat(flit[i]);
		//mylib::DisplayMat(fcand[i]);

		cv::dct(flit[i],flit[i]); //dct of target
	//	cv::dct(fcand[i],fcand[i]); //dct of candidate
	//	mylib::DisplayMat(flit[i]);
	//	mylib::DisplayMat(fcand[i]);
		for (int x=0; x<flit[i].size().height; x++)
			for (int y=0; y<flit[i].size().width; y++)
				for (int cc=0; cc< cn; cc++)
				{
					if ((x+y>2) || (lf == 4 && x>=2 && y>=2))
					{
						flit(x,y,i)[cc]=0;
					//	fcand(x,y,i)[cc]=0;
					}
					else
					{	
						if (x==0&& y ==0)
							index = 0;
						else if (x==0&&y ==1)
							index = 1;
						else if (x==1&&y ==0)
							index = 2;
						else if (x==2 &&y==0)
							index = 3;
						else if (x==1 && y == 1)
							index = 4;
						else if (x==0 && y ==2 )
							index = 5;
						else
							index = 6;
						
						flit(x,y,i)[cc] = floor(flit(x,y,i)[cc]/(*(tbl+index) * qfactor)+0.5);//quantize
						//fcand(x,y,i)[cc] = 0; //clear the LF component of candidate
						if(saveCodeLength)
						{

							if (flit(x,y,i)[cc] != 0)
								tempLength = int(floor(log(abs(flit(x,y,i)[cc]))/log(2.0)))+1;
							if (size().height<32)
								lightingDCTCoeffStat.push_back(tempLength);
							if (index==0)
								codeLength = codeLength + light_dc_length[tempLength]+tempLength;
							else
								codeLength = codeLength + light_ac_length[tempLength]+tempLength;
						}
						//flit(x,y,i)[cc] = flit(x,y,i)[cc]+0.5;
						flit(x,y,i)[cc] *= (*(tbl+index) *qfactor);
						//cout<<flit(x,y,i)[cc]<<endl;;
					
						
					}
				}
		//mylib::DisplayMat(flit[i]);
	//	mylib::DisplayMat(fcand[i]);
		cv::idct(flit[i],flit[i],DFT_SCALE);
		//cv::idct(fcand[i],fcand[i],DFT_SCALE);
	}
	//flit.Print();
	lighting = Tensor<T,cn>(flit);
	//this->SetBlock(Tensor<T,cn>(Tensor<double,cn>(*this) - fcand + flit + 128));
	//this->SetBlock(Tensor<T,cn>(fcand + flit + 128)); //target lighting + (target - target lighting)
	this->SetBlock(Tensor<T,cn>(Tensor<double,cn>(*this) + flit));
	return lighting;
	//this->SetBlock((*this) + BuildLightingPlane(tempTagLighting - tempCanLighting));
}

template <class T, size_t cn>
Tensor<T,cn> Tensor<T,cn>::LightingCorrection(const Tensor<T,cn>& changeTo, const Tensor<T,cn>& VQCodebook)
{
	//do this after quilting
	//change to is the original block (the lighting should be)
	Tensor<T,cn> tempTagLighting = changeTo.LSFitting();
	Tensor<T,cn> tempCanLighting = LSFitting();
	for (int i=0; i< tempTagLighting.size().height; i++)
		lightTag.push_back(tempTagLighting(i,0,0));
	for (int j=0; j< tempCanLighting.size().height; j++)
		lightCan.push_back(tempCanLighting(j,0,0));
	//tempTagLighting.Print();
	//SearchCodeword(tempTagLighting,VQCodebook).Print();
	//tempCanLighting.Print();
	//SearchCodeword(tempCanLighting,VQCodebook).Print();
	//search for vq codebook and use the quantized version to recover lighting.
	Tensor<T,cn> lighting = BuildLightingPlane(SearchCodeword(tempTagLighting,VQCodebook) - SearchCodeword(tempCanLighting,VQCodebook));
	//BuildLightingPlane(tempTagLighting-tempCanLighting).Print();
	//lighting.Print();
	this->SetBlock((*this) + lighting);
	return lighting;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::SearchCodeword(const Tensor<T,cn>& val, const Tensor<T,cn>& VQCodeBook)
{
	int label=-1;
	T dist = std::numeric_limits<T>::max();
	//T dist=(T)DBL_MAX;
	T temp_dist = 0;
	Vec<T,cn> dist_vec;
	if(VQCodeBook.size().width!= val.size().height)
	{
		CV_Error(CV_StsUnmatchedSizes,"the input polynomial has different size as codebook dim");
	}
	for (int i=0; i< VQCodeBook.size().height; i++)
	{

		dist_vec =Vec<T,cn>::all(0);
		temp_dist = 0;
		for(int j=0; j< VQCodeBook.size().width; j++)
			dist_vec+= mylib::VecPow(VQCodeBook(i,j,0) - val(j,0,0),2.0);
		for (int c=0; c<cn; c++)
			temp_dist+= dist_vec[c];
		if (temp_dist <= dist)
		{
			dist = temp_dist;
			label = i;
		}
	}
	return Tensor<T,cn>(VQCodeBook[0].row(label)).Transpose();
}
template< class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::MicroShift(Size3 win, double maxdev)
{
   int     i, j, ii, jj, k, l;
    int     hh, ww, hs, ws, hss, wss;
    int      nchan;
    long    seed;
    double   ampl;
    Tensor<T,cn> rst = this->Clone();
    seed = 125982;
    ampl = maxdev;			/* max pixel deviation */
    wss = win.width;			/* window size */
    hss = win.height;			/* window size */

    ws = 2*wss+1;
    hs = 2*hss+1;

    hh = size().height;
    ww = size().width;
    nchan = this->channels();
    srand(seed);		/* initializes noise generator */

    for(j = 0; j < hh; j++)
    {
	    for(i = 0; i < ww; i++)
	    {

	      //for(n = 0; n < nchan; n++)
	      //{
//		      pt1 = SignalData(ptsig1,n).c;
	//	      pt2 = SignalData(ptsig2,n).c;
          rst((j*hh+i)/ww,(j*hh+i)%ww,0) = this->operator()(j,i,0);
		      //*(pt2+j*hh+i) = *(pt1+j*ww+i);
	      //}
	    }
    }
    //rst.Display();
    for(j = hss; j < hh-hss; j += hs)
    {
	    for(i = wss; i < ww-wss; i += ws)
	    {
	      int  xx, yy;
	      xx = int(ampl*(2.*double(rand())/RAND_MAX-1.)+0.5);
	      yy = int(ampl*(2.*double(rand())/RAND_MAX-1.)+0.5);

	      jj = j+yy;
	      if(jj < hss) jj = hss;
	      else if(jj > hh-hss-1) jj = hh-hss-1;

	      ii = i+xx;
	      if(ii < wss) ii = wss;
	      else if(ii > ww-wss-1) ii = ww-wss-1;

	      for(l = -hss; l <= hss; l++)
	      {
		      for(k = -wss; k <= wss; k++)
		      {
		        //for(n = 0; n < nchan; n++)
		        //{
			      //pt1 = SignalData(ptsig1,n).c;
			      //pt2 = SignalData(ptsig2,n).c;
            rst(((j+l)*hh+i+k)/ww,((j+l)*hh+i+k)%ww,0) = this->operator()(jj+l,ii+k,0);
			      //*(pt2+(j+l)*hh+(i+k)) = *(pt1+(jj+l)*ww+(ii+k));
		      }
		    }
	    }
	  }
    return rst;
}

template< class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::MicroRotate(Size3 win, double maxdev)
{
    int     i, j, ii, jj;
    int     hh, ww, hs, ws, hss, wss;
    long    seed;
    double   ampl;
    Tensor<T,cn> rst = this->Clone();
    seed = 125982;
    ampl = maxdev;			/* max pixel deviation */
    wss = win.width;			/* window size */
    hss = win.height;			/* window size */

    ws = 2*wss+1;
    hs = 2*hss+1;

    hh = this->size().height;//SignalHeight(ptsig1,0);
    ww = this->size().width;//SignalWidth(ptsig1,0);

    srand(seed);		/* initializes noise generator */

    for(j = 0; j < hh; j++)
    {
	    for(i = 0; i < ww; i++)
	    {
        rst((j*hh+i)/ww,(j*hh+i)%ww,0) = this->operator()(j,i,0);
	    }
    }
    for(j = hss; j < hh-hss; j += hs)
    {
	    for(i = wss; i < ww-wss; i += ws)
	    {
	      double  phi;
	      phi = ampl*(2.*double(rand())/RAND_MAX-1.);
	      //printf("j=%d,i=%d:phi=%f\n",j,i,phi);

	      for(jj = 0; jj < hs; jj++)
	      {
		      for(ii = 0; ii < ws; ii++)
		      {
            rst(j+jj,i+ii,0) = this->Bilinear(ii,jj,i,j,ww,hs,ws,hs,ws,phi);//i is ---> j is vvvvv here
      			//*(pt2+(j+jj)*ww+i+ii) = bilinear(ii,jj,pt1+j*ww+i,ww,hs,ws,hs,ws,phi);
		      }
	      }
	    }
    }
    return rst;
}

template < class T, size_t cn>  typename Tensor<T,cn>::value_type Tensor<T,cn>::Bilinear(int i,int j, int offset_i, int offset_j, int ww, int hei1, int wid1, int hei2, int wid2, double phi)
{
  
 /* 
 * img1: original floating image ( hei1 X wid1 )
 *
 * imgg: floating value of rotated image ( hei2 X wid2 )
 *       at location [j][i]
 *
 */

    double  alpha, beta, gamma, delta;
    double  y, x;
    int    xp, xn, yp, yn;
    int    outp;
    Tensor<T,cn>::value_type ucoutpix,outpix;

/*
 *  (y,x) is the location of the current pixel (j,i) in image1 coordinates.
 *  xp: is the previous horizontal location
 *  xn: is the next horizontal location
 *  yp: is the previous vertical location
 *  yn: is the next vertical location
 *
 */
    y =  (j-hei2/2.+.5)*cos(phi)+(i-wid2/2.+.5)*sin(phi)+hei1/2.-.5;
    x = -(j-hei2/2.+.5)*sin(phi)+(i-wid2/2.+.5)*cos(phi)+wid1/2.-.5;

    //    if(y < 0. || y > hei1-1 || x < 0. || x > wid1-1) return((float) color);

    xp = (int) x;
    xn = xp+1;
    yp = (int) y;
    yn = yp+1;

    alpha = x-xp;
    beta  = yn-y;
    gamma = xn-x;
    delta = y-yp;

    if(xn >= wid1) xn = xp;
    if(yn >= hei1) yn = yp;
    outpix = this->operator()(yn+offset_j,xp+offset_i,0)*delta*gamma;
    //outpix  = *(img1+yn*ww+xp)*delta*gamma;
    outpix += this->operator()(yp+offset_j,xp+offset_i,0)*beta*gamma;
    outpix += this->operator()(yp+offset_j,xn+offset_i,0)*alpha*beta;
    outpix += this->operator()(yn+offset_j,xn+offset_i,0)*alpha*delta;
    for (int z=0; z<cn; z++)
    {
      outp = int(outpix[z]+.5);
      if(outp < 0) ucoutpix[z] = 0;
      else if(outp > 255) ucoutpix[z] = 255;
      else ucoutpix[z] = (T)outp;
    }
    return(ucoutpix);
}

template< class T, size_t cn> void Tensor<T,cn>::RecodeLighting(void)
{
	fstream lightfile;
	lightfile.open(".\\lighting.txt",ios::app);

	for (unsigned int i=0; i< lightingDCTCoeffStat.size(); i++)
		lightfile<<lightingDCTCoeffStat[i]<<endl;
	lightfile.close();
}

template< class T, size_t cn> int Tensor<T,cn>::GetLightingCodeLength(void) const
{
	return this->codeLength;
}
template< class T, size_t cn> void Tensor<T,cn>::SetLightingCodeLength(int l)
{
	this->codeLength=l;
}

template< class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::ComputeSAT(void) const
{
	//2D only
	Tensor<T,cn> S = this->ExtendHalfBoundary();
	S = this->ComputeSAT(S,Point3i(0,0,0),tsSize.Point3()-Point3i(1,1,1));
//	for ( int i=1; i< this->tsSize.height+1; i++)
//		for (int j=1; j<this->tsSize.width+1; j++)
//		{
//			S(i,j,0)=S(i-1,j,0)+S(i,j-1,0)-S(i-1,j-1,0)+S(i,j,0);
//		}	
	return S;
}
template<class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::ComputeSAT(Tensor<T,cn>& SAT, Point3i& sPos, Point3i& ePos) const
{
	//remember S is extended boundary by 1 to the left and up (zeros)
	for (int i = sPos.x+1; i <= ePos.x+1; i++)
		for (int j=sPos.y+1; j <= ePos.y +1; j++)
		{
			SAT(i,j,0)=SAT(i-1,j,0)+SAT(i,j-1,0)-SAT(i-1,j-1,0)+SAT(i,j,0);
		}
	return SAT;

}


template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::ComputeTPSS(double p) const
{
  CV_Assert(cn==1);//only implement gray scale
  Tensor<T,cn> G=this->Clone();
  ThinPlateSpline tps;		
  Mat temp;
  for (int t=0; t<this->size().depth; t++)
  {
    
    tps.load(G.GetFrame(t));
    G.SetFrame(t,tps.solve(p));
  }
  return G;
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

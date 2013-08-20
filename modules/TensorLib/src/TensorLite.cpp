#include "TensorLite.h"
//20130815 decoule Tensor to high level operations
//#include "Steerable.h"
//#include "LRI.h"
//#include "Metric.h"
namespace tensor{
//20130815 template<class T, size_t cn> Mat Tensor<T,cn>::stsim2_lse_weight = Mat(); 
/////////////2.1 constructor, copy  assignment and type conversion ////////////////////////
BufferGPU gbuf;
template<class T, size_t cn> Tensor<T,cn>::Tensor(void):Mat()
{
	this->tsOffset = Point();
	this->cFileName = "unknown";
	this->mxFrame = Mat();
	//subWinSize = Size3();
	//subWinStep = Size3();
	this->debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(int height, int width, typename Tensor<T,cn>::c_ref_type val): Mat()
{
   
	*this = Mat(height,width, CV_MAKETYPE(DataType<T>::depth,cn), Scalar(val));
	this->cFileName = "unknown";
	this->tsOffset = Point();
	this->mxFrame = Mat(height,width,CV_MAKETYPE(DataType<T>::depth,cn));
	this->debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(const Size& size, const Vec<T,cn>& val): Mat()
{
	*this = Mat(size,  CV_MAKETYPE(DataType<T>::depth,cn),Scalar(val));
	this->cFileName = "unknown";
	this->tsOffset = Point();
	this->mxFrame = Mat(size.height,size.width,CV_MAKETYPE(DataType<T>::depth,cn));
	//subWinSize = Size3();
	//subWinStep = Size3();
	this->debugtrigger=false;
}
template<class T, size_t cn> Tensor<T,cn>::Tensor(const Size& size):Mat()
{
	*this = Mat(size, CV_MAKETYPE(DataType<T>::depth,cn));
	this->cFileName = "unknown";
	this->tsOffset = Point();
	this->mxFrame = Mat(size.height,size.width,CV_MAKETYPE(DataType<T>::depth,cn));
	//subWinSize = Size3();
	//subWinStep = Size3();
	this->debugtrigger=false;
}
template<class T, size_t cn> Tensor<T,cn>::Tensor(const Tensor<T,cn>& ts): Mat(ts)
{
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
	this->cFileName = "unknown";
	this->tsOffset = Point();
	if (mt.type() != CV_MAKETYPE(DataType<T>::depth,cn))
	{
		Tensor<T,cn> rst(this->size());
		this->convertTo(rst,DataType<T>::depth);
		*this = rst;
	}
	this->debugtrigger=false;
}

template<class T, size_t cn> Tensor<T,cn>::Tensor(const string cFileName)//:subTensor(Tensor_<T,cn>())
{

	this->Load(cFileName);
	this->tsOffset = Point();

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
	rst.SetFrame(mylib::toComplex(*this));
	rst.debugtrigger = this->debugtrigger;
	return rst;
}

template<class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::operator= (const Tensor<T,cn>& ts)
{
//	CV_Assert(ts.type() == CV_MAKETYPE(this->depth(), this->channels()));
	this->Mat::operator=(ts);
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

	cv::namedWindow(tempName,flag|WINDOW_FREERATIO/*|CV_GUI_EXPANDED*/);//if Qt enabled, you can uncomment this
	cv::moveWindow(tempName.c_str(),100,100);
	
		cv::imshow(tempName,tempTs);
		cv::waitKey();
		cv::destroyWindow(tempName.c_str());

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
	cv::namedWindow(tempName,flag|WINDOW_KEEPRATIO/*|CV_GUI_EXPANDED*/);//if Qt enabled, you can uncomment this
	cv::moveWindow(tempName.c_str(),100,100);

		cv::imshow(tempName,tempTs);
		cv::waitKey(sec);
		cv::destroyWindow(tempName.c_str());

}
//20130815 move this to utility module inside Tensor (not mylib)
/*
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
*/

template<class T, size_t cn> void Tensor<T,cn>::Load(string cFileName) 
{
	this->cFileName = cFileName;
	//cv::VideoCapture cap;
	//cap.open(cFileName);
	//cap>>this->mxFrame;//read next frame
	//cv::Mat grayone;
	//*if (!this->mxFrame.data) //if cannot be read as a stream
//	{
	//	cout<<"cannot be read by stream, use imread()"<<endl;
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
		this->mxFrame.convertTo(this->mxFrame,DataType<T>::depth);//new added
                *this = mxFrame.clone();
		cout<<"size: "<<this->size(); 
		cout<<", type: "<<this->type();
		cout<<", channel: "<<this->channels()<<endl;
	/*}
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
		int sz[] = {1,this->mxFrame.rows,this->mxFrame.cols};
		cout<<this->mxFrame.size()<<endl;
    		this->mxFrame.convertTo(this->mxFrame,DataType<T>::depth);//new added
    		//cout<<mxFrame.channels()<<","<<mxFrame.type()<<endl;
		this->create(3,sz,this->mxFrame.type());
		this->size() = Size3(sz[1],sz[2],sz[0]);
		this->SetFrame(0,this->mxFrame.clone());
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
	*/
	this->debugtrigger=false;
	//this->Display();
}

//////////////////////// 2.3 accessors

template<class T, size_t cn> typename Tensor<T,cn>::ref_type Tensor<T,cn>::operator() (int x, int y)
{
	CV_Assert( CV_MAKETYPE(DataType<T>::depth,cn) == this->type());
	return this->at<value_type>(x,y);

}

template<class T, size_t cn> typename Tensor<T,cn>::c_ref_type Tensor<T,cn>::operator() (int x, int y) const
{
	//int t = this->type();
	//int t2 = CV_MAKETYPE(DataType<T>::depth,cn) ;
	CV_Assert( CV_MAKETYPE(DataType<T>::depth,cn) == this->type());

	return this->at<value_type>(x,y);

}

template<class T, size_t cn> typename Tensor<T,cn>::ref_type Tensor<T,cn>::operator[] (const Point& pos)
{
	return this->operator()(pos.x,pos.y);
}
template<class T, size_t cn> typename Tensor<T,cn>::c_ref_type Tensor<T,cn>::operator[](const Point& pos) const
{
	return this->operator()(pos.x,pos.y);
}





template< class T, size_t cn> Point Tensor<T,cn>::offset() const
{
	return this->tsOffset;
}

template< class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Clone(void) const
{
  	std::lock_guard<std::mutex> lock(this->tsmutex);
	//tsmutex.lock();
  	Tensor<T,cn> rst(this->Mat::clone());
	rst.SetFileName(this->cFileName);
	rst.SetOffset(this->tsOffset);
	//rst.SetSubWinSize(subWinSize);
	//rst.SetSubWinStep(subWinStep); 
	//rst.lightTag = lightTag;
	//rst.lightCan = lightCan; 

	return rst;
}

template< class T, size_t cn> void Tensor<T,cn>::SetFrame( const Mat& frm)
{
	//use to modify(or assignment by copy) one frame
		*this = frm;
}


template<class T, size_t cn> void Tensor<T,cn>::Print(const string& fname, bool tofileonly) const
{

	fstream logfile;
	logfile.open("./"+fname+".txt",ios::out); 
	//bool tofileonly=false;
	if (!tofileonly)
	{
		cout<<"File: "<<fname<<" Display "<<this-> size().height <<" "<<this-> size().width<<" Data type: "<<this-> type()<<endl;
		cout<<"cn="<<cn<<endl;
	}
	//logfile<<"File: "<<cFileName<<"Display "<<size().height <<" "<<size().width<<" Data type: "<<type()<<endl;
	
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
						cout<< double(this->operator()(x,y)[c]) <<"\t";
					logfile<< double(this->operator()(x,y)[c])<<"\t";
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
	logfile.close();
}



template <class T, size_t cn> void Tensor<T,cn>::SetFileName(const string& name)
{
	this->cFileName = name;
}

template <class T, size_t cn> void Tensor<T,cn>::SetOffset(const Point& pos)
{
	this->tsOffset = pos;
}

//template <class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::GetBlockRef(const Cube& roi)
//{
//		CV_Assert(  this->dims == 3 && (unsigned)i < (unsigned)tsSize.depth);
//		//mxFrame = this->row(i);
//     //Mat* arrays[] = {this};
//    Mat* arrays[] = {(Mat*)this};
//    Mat planes[1];
//    NAryMatIterator it((const Mat**)arrays, planes,1);
//    mxFrame = it.planes[i].reshape(cn,this->size().height);
//    //mylib::DisplayMat(mxFrame);
//    return mxFrame(Rect(roi.y,roi.x,roi.width,roi.height));
//}
template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::GetBlock(const Rect& roi)
{
	Tensor<T,cn> rst = this->Mat::operator()(roi);
	rst.SetFileName(this->cFileName);
	rst.SetOffset(roi.offset());
	return rst;
}
template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::GetBlock(const Rect& roi) const
{
	Tensor<T,cn> rst = this->Mat::operator()(roi);
	rst.SetFileName(this->cFileName);
	rst.SetOffset(roi.offset()+this->offset());
	return rst;
}

template<class T, size_t cn> void Tensor<T,cn>::Ref(const Rect& roi, Tensor<T,cn>& dst) const
{
	dst.Mat::operator=(this->Mat::operator()(roi));
  	dst.tsOffset = roi.offset();
	dst.SetFileName(this->cFileName);
	dst.SetOffset(roi.offset()+this->offset());
	//dst.SetSubWinSize(this->subWinSize);
	//dst.SetSubWinStep(this->subWinStep);
}

template <class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::SetBlock(const Point& pos, const Tensor<T,cn>& ts)
{
	CV_Assert(ts.type() == type());
	AssertRange(pos,ts.size());
	
        Rect roi(pos, ts.size());
	cv::Mat tempMat = GetBlock(roi);

	ts.copyTo(tempMat);
	return *this;
}
template <class T, size_t cn> Tensor<T,cn>& Tensor<T,cn>::SetBlock(const Tensor<T,cn>& ts)
{
	return SetBlock(Point(0,0),ts);

}

template<class T, size_t cn> void Tensor<T,cn>::AssertRange(const Point& pos, const Size& sz) const
{
	CV_Assert( pos.x + sz.height <= this-> size().height &&
				  pos.y + sz.width  <= this-> size().width );
}

template<class T, size_t cn> void Tensor<T,cn>::SaveBlock(const string& cFileName, bool isGray)
{
		if (isGray)
		{
			Mat grayone = *this;
			cv::cvtColor(grayone,grayone,COLOR_RGB2GRAY);
			cv::imwrite(cFileName,grayone);
		}
		else
			cv::imwrite(cFileName,*this);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Crop(const Point& pos, const Size& sz) const
{
  return this-> GetBlock(Rect(pos,sz)).Clone();
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
template<class T, size_t cn> bool Tensor<T,cn>::IsInside(const Point& pos) const
{
	if ( pos.x >= tsOffset.x && pos.x < tsOffset.x + this-> size().height &&
		pos.y >= tsOffset.y && pos.y < tsOffset.y + this-> size().width  )
		return true;
	else
		return 0;
}


template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator()(const Rect& roi)
{	
	//return GetBlock(roi);
	AssertRange(roi.tl(),roi.size());
	Tensor<T,cn> dst(roi.size());
		//if ( i < dst.size().depth )
		dst.SetFrame(this->operator()(roi));
		//else
		//	dst.vdStream.push_back((*this)[i](Rect(roi.y,roi.x,roi.width,roi.height)));
	dst.tsOffset = roi.offset()+this->offset();
	dst.cFileName = this->cFileName;	
	//dst.SetSubWinSize(this->subWinSize);
	//dst.SetSubWinStep(this->subWinStep);
	return dst; 
}

////////2.4 arithmatic operations //////////////////////////////////////////

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator+(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) + Tensor<T,cn>(size().height,size().width,s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator-(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) - Tensor<T,cn>(size().height,size().width,s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator*(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) * Tensor<T,cn>(size().height,size().width, s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator/(typename Tensor<T,cn>::c_ref_type s) const
{
	return (*this) / Tensor<T,cn>(size().height,size().width, s);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator*(const Tensor<T,cn>& ts) const
{

	Tensor<T,cn> rst = Clone();

		if (cn==2 || cn == 6)
			cv::mulSpectrums(rst,ts,rst,DFT_ROWS,false); // use CCS strucrue for complex mul
		else
    		{
                        rst = rst.mul(ts);
    		}
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::operator/(const Tensor<T,cn>& ts) const
{
	Tensor<T,cn> rst = Clone();
	if (cn==2|| cn==6) //complex div
	{
			Mat tempMat;
			//mylib::DisplayMat(rst.GetFrame(i));
			cv::mulSpectrums(ts,ts,tempMat,DFT_ROWS,true);
			//mylib::DisplayMat(tempMat);
			vector<Mat> temp;
			cv::split(tempMat, temp);
			temp[1] = temp[0];
			merge(temp,tempMat);
			//mylib::DisplayMat(tempMat);
			cv::mulSpectrums(rst,ts,rst,DFT_ROWS,true);
			//mylib::DisplayMat(rst.GetFrame(i));
			//Vec<T,cn> a= rst(0,0,0);
			//Vec<T,cn> b = tempMat.at<Vec<T,cn>> (0,0);
			rst /= tempMat;
			//mylib::DisplayMat(rst.GetFrame(i));
	}
	else
	{
			rst /= ts;
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

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Abs(void) const
{
	
	Tensor<T,cn> rst(this->size());
	if (cn == 2) // if data type is not float/double this will throw type exception
	{
			cv::mulSpectrums(*this,*this,rst,DFT_ROWS,true);
			cv::pow(rst,0.5,rst);
	}
	else
	{
                        rst.SetFrame(cv::abs(*this));
	}
	return rst;
}


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
		cv::absdiff(*this,v,rst);
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
		cv::absdiff(*this,ts,rst);
	}
	return rst;
}

template<class T, size_t cn> Tensor<T,1> Tensor<T,cn>::Real(void) const
{
	if (cn==2)
	{
		vector<Mat> temp;
		cv::split(*this,temp);
                Tesnor<T,1> rst(temp[0]);
		return rst;
	}
	else
		return this->Clone();
}
template<class T, size_t cn> Tensor<T,1> Tensor<T,cn>::Imag(void) const
{
	if (cn==2)
	{
		vector<Mat> temp;
		cv::split(*this,temp);
		Tensor<T,1> rst(temp[1]);
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
				divide( 1., *this, rst);
				if( ipower == -1 )
					return rst;
				ipower = -ipower;
			}
			else
				rst = this->Clone();

			switch( ipower )
			{
			case 0:
				rst = Tensor<T,2>(this-> tsSize.height,this-> tsSize.width,Vec<T,2>(1,0));
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
	T a = T(this->size().area());
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
				tempmean[kk]+=(this->operator()(ii,jj)[kk]);
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
				tempmean[kk]+=(this->operator()(ii,jj)[kk]);
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
			{
				vecMat.at<Vec<T,cn>>(0,count)= this->operator()(i,j);
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
		rst = Tensor<T,cn>(tsSize.width,tsSize.height);
		rst.tsOffset = tsOffset;
		//rst.subWinSize=subWinSize;
		//rst.subWinStep=subWinStep;
	}
	//	(Size3(this->size().width,this->size().height,this->size().depth));
		cv::transpose(*this,rst);
	return rst;
}


template <class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::ExtendBoundary(Size extSz, typename Tensor<T,cn>::value_type val) const
{
	Tensor<T,cn> rst(size()+(extSz*2),val);
	rst.tsOffset = Point(extSz.height,extSz.width);
	rst.SetBlock(rst.tsOffset,*this);
	//rst.SetSubWinSize(subWinSize);
	//rst.SetSubWinStep(subWinStep);
	return rst;
}

template <class T, size_t cn>
Tensor<T,cn> Tensor<T,cn>::ExtendHalfBoundary(Size extSz, typename Tensor<T,cn>::value_type val, bool which_side) const
{
	Tensor<T,cn> rst(size()+(extSz),val);
	if (which_side)
		rst.tsOffset = Point();
	else
		rst.tsOffset = Point(extSz.height,extSz.width);
	rst.SetBlock(rst.tsOffset,*this);
	//rst.SetSubWinSize(subWinSize);
	//rst.SetSubWinStep(subWinStep);
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalMean(const Mat& ker,const Size& subWinStep) const
{
	return Filter2D(ker).SubSample(subWinStep);
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalMean(const Size& subWinSize,const Size& subWinStep) const
{
	//tsSize.operator-(subWinSize);
	Size_<double> temp1 = tsSize-subWinSize;
        Size_<double> temp2 = subWinStep;
	//(Size3_<double>(Size3(tsSize-subWinSize))/Size3_<double>(Size3(subWinStep))).Ceil()+Size3(1,1,1);
	Size temp3 = (temp1/temp2)//.Ceil();
	Size rstSize = temp3 + Size(1,1);
	Tensor<T,cn> rst(rstSize);
	for (int x=0; x< rstSize.height; x++)
		for (int y=0; y< rstSize.width; y++)
			{
				rst(x,y)= GetBlock(Rect(x*subWinStep.height,y*subWinStep.width,
					x*subWinStep.height+subWinSize.height<tsSize.height?subWinSize.height:tsSize.height-x*subWinStep.height,
					y*subWinStep.width +subWinSize.width <tsSize.width ?subWinSize.width :tsSize.width -y*subWinStep.width
					)).Mean();
			}
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalVariance( const Tensor<T,cn>& mu, const Mat& ker, const Size& subWinStep) const
{
	auto temp =  ((*this)*(this->Conjugate())).Filter2D(ker).SubSample(subWinStep)- (mu*mu.Conjugate());
  return temp.Abs();
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::LocalVariance( const Tensor<T,cn>& mu, const Size& subWinSize, const Size& subWinStep) const
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
	cv::Laplacian(*this,rst,this->depth(),1,1,0,BORDER_REFLECT_101);
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

			cv::filter2D(*this,rst,-1,ker);
		if (boundary == (int)FilterBoundary::FILTER_BOUND_EXTEND)
		{
			CV_Error(CV_StsNotImplemented,"unsupport extend boundary version");
			return rst;
		}
		else if (boundary == (int)FilterBoundary::FILTER_BOUND_VALID)
		{
			Size validSize(tsSize.height - ker.size().height +1, tsSize.width - ker.size().width+1);
			Point validOffset(ker.size().height/2, ker.size().width/2);
			return rst.Crop(validOffset,validSize);
		}
		else
			return rst;
	
}
//this is slow
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::SubSample(const Size& subrate)
{
	Size_<double> temp1 = tsSize;
	Size_<double> temp2 = subrate;
	Size rstSize = (temp1/temp2);//.Ceil();
	
	if (rstSize.height==0 || rstSize.width==0 )
		CV_Error(CV_StsBadSize, "sample rate too high");
	Tensor<T,cn> rst(rstSize);
		for (int i =0; i*subrate.height< tsSize.height; i++)
			for (int j=0; j*subrate.width< tsSize.width; j++)
			{
				rst(i,j) = this->operator()(i*subrate.height,j*subrate.width);
			}
	return rst;
}
template<class T, size_t cn> Vec<T,cn> Tensor<T,cn>::Min() const
{
	if (cn != 1)
		CV_Error(CV_StsNotImplemented,"do not work for multi-channel yet");
	double p=DBL_MAX;
	double temp=0;
		cv::minMaxLoc(*this,&temp);
		if (temp < p)
			p= temp;
	Vec<T,cn> rst;
	rst = cv::saturate_cast<T>(p);
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::MaxClip(const Vec<T,cn>& value) const
{
	Tensor<T,cn> rst = this->Clone();
	for (int i=0; i< this->size().height; i++)
		for (int j=0; j<this->size().width; j++)
			{
				//Vec<T,cn> temp = rst(i,j,k);
				//temp = (T)cv::norm(rst(i,j,k));
				//temp = (T)cv::norm(value);
				if (cv::norm(rst(i,j)) > cv::norm(value))
					rst(i,j)=value;
			}
	//rst.Print();
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::MinClip(const Vec<T,cn>& value) const
{
	Tensor<T,cn> rst = this->Clone();
	for (int i=0; i< this->size().height; i++)
		for (int j=0; j<this->size().width; j++)
			{
				//Vec<T,cn> temp = rst(i,j,k);
				//temp = (T)cv::norm(rst(i,j,k));
				//temp = (T)cv::norm(value);
				if (cv::norm(rst(i,j)) < cv::norm(value))
					rst(i,j)=value;
			}
	//rst.Print();
	return rst;
}
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Flip(bool byX, bool byY) const
{
	Tensor<T,cn> rst = this->Clone();
	if (byX)
	{
	for (int y=0; y< this->size().width; y++)
		for (int x=0; x< this->size().height/2; x++)
		{	
			rst(x,y) = this->operator()(tsSize.height-x-1,y);
			rst(tsSize.height-x-1,y)=this->operator()(x,y);
		}
	}
	if (byY)
	{
	for (int x=0; x< this->size().height; x++)
		for (int y=0; y< this->size().width/2; y++)
		{	
			rst(x,y) = this->operator()(x,tsSize.width-y-1);
			rst(x,tsSize.width-y-1)=this->operator()(x,y);
		}
	}
	
	return rst;
}	
template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::FreqComplexFilter(cv::Mat& kernel, bool conj) const
{
	Tensor<T,cn> rst(size());
	CV_Assert(this->channels() == kernel.channels() && this->channels()== 2); 
		rst = mylib::FreqComplexFilter(*this,kernel,conj);
	return rst;
}

template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::DFT(void) const
{

	Tensor<T,cn> rst(size());
	CV_Assert(cn ==2 || cn== 6); 
	#ifdef USE_GPU
		  gbuf.gI1.upload(*this);
		  gbuf.gI1.convertTo(gbuf.t1, CV_32FC2, gbuf.stream);
             	  gpu::dft(gbuf.t1, gbuf.t2, this->size(),0,gbuf.stream);
		  gbuf.t2.convertTo(gbuf.gI2, CV_64FC2,gbuf.stream);//make this adaptive
		  gbuf.gI2.download(rst); 
	#else
	cv::dft(*this,rst);
	#endif

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
		gbuf.gI1.upload(*this);
                gbuf.gI1.convertTo(gbuf.t1, CV_32FC2, gbuf.stream);
		gpu::dft(gbuf.t1, gbuf.t2, this->size(), DFT_INVERSE, gbuf.stream);
		gbuf.t2.convertTo(gbuf.gI2, CV_64FC2,gbuf.stream);
                gpu::split(gbuf.gI2, gbuf.v, gbuf.stream); 
                gpu::divide(gbuf.v[0], double(rst[i].size().height*rst[i].size().width), gbuf.v[0],1, -1, gbuf.stream);
                gpu::divide(gbuf.v[1], double(rst[i].size().height*rst[i].size().width), gbuf.v[1],1, -1, gbuf.stream); 
                gpu::merge(gbuf.v, gbuf.gI2, gbuf.stream);
		gbuf.gI2.download(rst);
		//vector<Mat> temp;
                //cv::split(rst[i],temp);
                //Tensor<double,1>(temp[0]/temp[0].size().height/temp[0].size().width).Display(); 
                //mylib::DisplayMat(temp[0]);
                //cv::divide(this->size().area(), rst[i], rst[i], -1);
		#else	
		cv::dft(*this,rst,DFT_INVERSE|DFT_SCALE);
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
		rst=DFTShift(*this);
	return rst;
}


template<class T, size_t cn> Tensor<T,cn> Tensor<T,cn>::Normalize(void) const
{
	CV_Assert(cn == 1);
	Tensor<T,cn> rst(size());
		cv::normalize(*this,rst,0,255,NORM_MINMAX);
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

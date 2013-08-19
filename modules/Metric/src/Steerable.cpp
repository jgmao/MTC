
#include "Steerable.h"

namespace metric{
Steerable::Steerable(void)
{
	max_ht=0;
	order=0;
	twidth=0;
	nbands=0;
	ctr=Point3i();
}


Steerable::~Steerable(void)
{
}

Steerable::Steerable(Steerable::c_real_ref im)
{
	
	max_ht = int(floor(log(double(im.size().height < im.size().width ? im.size().height : im.size().width))/log(2.0))-2);
	order = 0; //by default
	nbands = 0;
	twidth = 1; // by default
	dims = im.size();
	ctr.x = (int)ceil((double(dims.height)+0.5)/2.0); //= ceil((dims+0.5)/2);
	ctr.y = (int)ceil((double(dims.width)+0.5)/2.0);
	ctr.z = 0;

}


vector<Steerable::data_type>& Steerable::buildSCFpyr(Steerable::c_data_ref im, int nLevel, int nDir, int twidth, bool subsample)
{
	//step 0, make sure im is the size of 2's powers
	//step 1, get the dft
	//vector<STBlock> rst;//rst should be a series of complex image
	if (im.channels() == 1 || im.channels() == 3)
	{
		//unsupport
		exit(-1);
	}
	else if (im.size().height != im.size().width)
	{
		exit(-2);
	}
	else if (im.size().depth != 1)
	{
		exit(-3);
	}

	this->order = nDir -1;
	this->twidth = twidth;
	this->nbands = nDir;
	this->max_ht = nLevel;
	cv::Mat hi0mask = raisedCosine(im.size().height, 1, 0.5, twidth); // first level half band high
	cv::Mat lo0mask = raisedCosine(im.size().height, 0.5, 1, twidth); // first level half band low;
	cv::Mat c_hi0mask = toComplex(hi0mask);
        cv::Mat c_lo0mask = toComplex(lo0mask);

	data_type imdft = im.DFT().DFTShift();
	data_type lo_imdft = imdft.FreqComplexFilter(c_lo0mask);
        data_type hi_imdft = imdft.FreqComplexFilter(c_hi0mask);
	//imdft.Print();
	//STBlock hipart = mylib::FreqComplexFilter(imdft,toComplex(hi0mask));
	buildSCFpyrLevs(lo_imdft,pyr_freq,nLevel,nDir,twidth, subsample);
	pyr_freq.push_back(hi_imdft);
	pyr_space.clear();
	
	for (unsigned int i=0; i< pyr_freq.size(); i++)
	{
   // pyr_freq[i].Abs().Real().Display();
    //pyr_freq[i].Print("freq");
		pyr_space.push_back(pyr_freq[i].DFTShift().IDFT());
  //  pyr_space[i].Print("space");
		//pyr_space.push_back(IDFT(DFTShift(pyr_freq[i])));
	}
	return pyr_freq;
}
void Steerable::buildSCFpyrLevs(Steerable::data_ref loDft, vector<Steerable::data_type>& rst, int nLevel, int nDir, int twidth, bool subsample)
{
	//the input loDft must be a complex image!!
	int next_twidth;
	if (nLevel == 0) // stop and save the low-pass band
	{
		rst.push_back(loDft);// do not use Low pass band //gji Feb 17 2012, change it back on Aug 20, 2012
		return;
	}
	//int order = nDir -1; 
	//double aConst = (pow(2.0,2*order))*(pow(mylib::factorial(order),2.0))/(nDir* mylib::factorial(2*order));
	//loDft.Print();
	cv::Mat hi1mask = raisedCosine(loDft.size().height,0.5,0.25,twidth+1);
//	mylib::DisplayMat(hi1mask);
	hi1mask = toComplex(hi1mask);
	cv::Mat ag0mask;
	cv::Mat shiftfactor(loDft.size().height,loDft.size().width,loDft.type(),mylib::complexToScalar(pow(std::complex<double>(0,-1),nDir-1))); 
	data_type imbp = loDft.FreqComplexFilter(hi1mask);

	for (int i=0; i< nDir; i++)
	{
		ag0mask = angleFilter(loDft.size().height,nDir,i);
		//mylib::DisplayMat(ag0mask);
		ag0mask = toComplex(ag0mask);
		ag0mask = mylib::FreqComplexFilter(ag0mask,shiftfactor);
	
		rst.push_back(imbp.FreqComplexFilter(ag0mask));
	}

	if (subsample)
		next_twidth = twidth;
	else
		next_twidth = twidth+1;

	if (subsample)//normalize /4 on Dec 30 2012
		loDft = loDft.Crop(cv::Point3i(loDft.size().height/4, loDft.size().width/4,0),loDft.size()/Size3(2,2,1))/4;
	
	cv::Mat lo1mask = raisedCosine(loDft.size().height,0.5,1,next_twidth);
  //mylib::DisplayMat(lo1mask,"lowmask");
	lo1mask = toComplex(lo1mask);
	data_type lo_loDft = loDft.FreqComplexFilter(lo1mask);
	buildSCFpyrLevs(lo_loDft,rst,nLevel-1,nDir,next_twidth,subsample);
	//STBlock alpha = mylib::mod(xrCos + PI, 2*PI) - PI;
	//STBlock yCosn = 2*sqrt(aConst) * (mylib::cos(xCos
	//return vector<STBlock>();
}

cv::Mat Steerable::raisedCosine(int maskSize, double passPosition, double stopPosition, int octave, int channels)
{
	//firstly first, only consider squared blocks!
	//2D raised cosine filter
	//double step = 2*CV_PI/maskSize;
	double posX, posY, radius;
	double passPosDB = log(passPosition)/log(2.0);
	double stopPosDB = log(stopPosition)/log(2.0);
	double width = abs(stopPosition - passPosition);
	double scale = 1/(width*4)/pow(2.0,octave-1);
	int center = (int)ceil((double(maskSize)-0.5)/2.0);
	cv::Mat mask(maskSize,maskSize,6);
	for (int x=0; x < maskSize; x++)
		for (int y=0; y < maskSize; y++)
		{

			posX = double(x - center)/double(maskSize)/scale;
			posY = double(y - center)/double(maskSize)/scale;
			radius =log(posX*posX + posY*posY)/log(2.0)/2.0;	
			
			if (passPosition < stopPosition)// low pass filter
			{
				if (x == center && y == center)
				{
					mask.at<double>(y,x) = 1.0;
					continue;
				}
				if (radius <= passPosDB) //eg. -1 dB (PI/2)
					mask.at<double>(y,x) = 1.0;
				else if (radius >= stopPosDB) //eg. 0dB (PI)
					mask.at<double>(y,x) = 0;
				else
  				{
					double temp = sqrt(1.0 - pow(cos(CV_PI*(radius-stopPosDB)/2),2.0));
					mask.at<double>(y,x) = temp;//sqrt(1.0 - pow(cos(PI*(radius-stopPosDB)/4/width),2.0));
				}
			}
			else // high pass filter
			{
				if (x == center && y == center)
				{
					mask.at<double>(y,x) = 0.0;
					continue;
				}
				if (radius >= passPosDB)
					mask.at<double>(y,x) = 1.0;
				else if (radius <= stopPosDB)
					mask.at<double>(y,x) = 0;
				else
				{
					double temp = cos(CV_PI*(radius-passPosDB)/2);
					mask.at<double>(y,x) = temp;
				}
			}
		}

		return mask;

}
//
cv::Mat Steerable::angleFilter(int maskSize, int nBands, int direction)
{
	int order  = nBands -1 ;
	double temp;
	double angle;
	int center = (int)ceil((double(maskSize)-0.5)/2.0);
	//c = (pow(2.0,2*order))*(pow(mylib::factorial(order),2.0))/(nDir* mylib::factorial(2*order));
	double a = mylib::factorial(order);
	double b = mylib::factorial(2*order);
	double c = pow(2.0,order)*a / sqrt((double)nBands*b);
	//double c = pow(2.0, (double)order) * (double)mylib::factorial((long)order)/sqrt(double(nBands*mylib::factorial(long(order*2))));
	cv::Mat mask(maskSize,maskSize,CV_64F);
 	for (int x=0; x < maskSize; x++)
		for (int y=0; y < maskSize; y++)
		{
			//posX = double(x - center)/double(maskSize);
			//posY = double(y - center)/double(maskSize);
			angle = atan2(double(y- center),double(x - center));
			angle = angle -CV_PI*(direction)/nBands;
		
			if ( abs(angle) < CV_PI/2 || angle < -3*CV_PI /2) 
			{
				temp = 2*c*pow(cos(angle),order); 
				mask.at<double>(y,x) = temp;
			}
			else
				mask.at<double>(y,x) = 0.0;
		}
	return mask;

}
//
cv::Mat Steerable::getChannel(cv::Mat& im, int n)
{
	cv::Mat rst;
	if (im.channels() < n)
	{	
		cout<<"the input image do not have "<<n<< " channels, exit."<<endl;
		return rst;
	}
	else
	{
		vector<cv::Mat> temp;
		cv::split(im,temp);
		rst = temp[n];
		return rst;
	}

}
//
cv::Mat Steerable::toComplex(cv::Mat& A, cv::Mat& B)
{
	cv::Mat rst;
	vector<cv::Mat> temp;
	if (A.channels()==1 && B.channels()==1)
	{
		temp.push_back(A);
		temp.push_back(B);
	}
	else if (A.channels()==3 && B.channels()==3)
	{
		temp.push_back(getChannel(A,0));
		temp.push_back(getChannel(B,0));
		temp.push_back(getChannel(A,1));
		temp.push_back(getChannel(B,1));
		temp.push_back(getChannel(A,2));
		temp.push_back(getChannel(B,2));
	}
	else
	{
		cout<<"A and B channel are different, or they are neither 1 channel or 3 channels image\n";
		return rst;
	}
	cv::merge(temp,rst);
	return rst;

}
//
cv::Mat Steerable::toComplex(cv::Mat& A)
{
	cv::Mat AZim = cv::Mat::zeros(A.rows,A.cols,A.type());
	return toComplex(A,AZim);
}
//
//STBlock Steerable::toComplex(STBlock& A)
//{
//	STBlock rst(A.size,A.GetDataType()+((A.GetChannel()*2-1)<<CV_CN_SHIFT),A.GetChannel()*2);
//	int i=0;
//	for (A.iter = A.GetBeginning();A.iter!= A.GetEnding(); A.iter++)
//	{
//		rst.SetFrame(i,toComplex(*A.iter));
//		i++;
//	}
//	return rst;
//}
//
cv::Mat Steerable::getMagnitude(cv::Mat& A)
{
	cv::Mat rst;
	if (A.channels()==1 || A.channels() == 3)
		rst = cv::abs(A);
	else if (A.channels() == 2)
	{
		cv::magnitude(getChannel(A,0),getChannel(A,1),rst);
	}
	else if (A.channels() == 6)
	{
		vector<cv::Mat> temp;
		cv::Mat a,b,c;
		cv::magnitude(getChannel(A,0),getChannel(A,1),a);
		cv::magnitude(getChannel(A,2),getChannel(A,3),b);
		cv::magnitude(getChannel(A,4),getChannel(A,5),c);
		temp.push_back(a);
		temp.push_back(b);
		temp.push_back(c);
		cv::merge(temp,rst);
	}
	else
	{
		cout<<"undefined behavior, the channel of image is neither, 1,2,3,or 6\n";
	}
	return rst;
}
//
//STBlock Steerable::getMagnitude(STBlock & A)
//{
//	STBlock rst(A.size,A.GetDataType());
//	int i=0;
//	for (A.iter = A.GetBeginning(); A.iter!= A.GetEnding(); A.iter++)
//	{
//		rst.SetFrame(i,getMagnitude(*A.iter));
//		i++;
//	}
//	return rst;
//}
//
cv::Mat Steerable::DFT(cv::Mat& A)
{
	cv::Mat rst;
	cv::Mat temp;
	temp = convertTo(A, CV_64F + ((A.channels()-1)<<CV_CN_SHIFT));
	if (A.channels() == 1 || A.channels() == 3)
		temp = toComplex(temp);
	cv::dft(temp,rst);
	return rst;
}

cv::Mat Steerable::IDFT(cv::Mat & A)
{
	cv::Mat rst;
	if (A.channels()== 1 || A.channels()==3)
		cv::dft(toComplex(A),rst,CV_DXT_INVERSE|DFT_SCALE);
	else
		cv::dft(A,rst,CV_DXT_INVERSE|DFT_SCALE);
	return rst;
}
//
//
//STBlock Steerable::DFT(STBlock & A)
//{
//	//only implemented 1 frame (2D)
//	// if depth >1, will use 3D DFT instead
//	//STBlock rst;
//	cv::Mat rst;
//	if (A.size.depth == 1)
//	{
//		A.iter = A.GetBeginning();
//		rst = DFT(*A.iter);
//		//rst = temp;
//	}
//	else
//	{
//		cout<<"behavior dose not define, the depth > 1\n";
//	}
//	return rst;
//}
//
//STBlock Steerable::IDFT(STBlock &A)
//{
//	//STBlock rst;
//	cv::Mat rst;
//	if (A.size.depth==1)
//	{
//		A.iter = A.GetBeginning();
//		rst =IDFT(*A.iter);
//		//rst =temp;
//	}
//	else
//	{
//		cout<<"behavior dose not define , the depth >1\n";
//	}
//	return rst;
//}
//
//STBlock Steerable::DisplayLevel(int n, vector<STBlock>& pyr)
//{
//
//	if (n>=(int)pyr.size() )
//		exit(1);
//	else
//	{
//		int step =5;
//		BlockSize<int> size = pyr[n*nbands].size;
//		
//		//int dataType = pyr[n*nbands].GetDataType()- ((pyr[n*nbands].GetChannel()-1)<<3);
//		if ( nbands == 4)
//		{
//			STBlock temp((size*BlockSize<int>(2,2,1)),CV_8U);
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands]))),cv::Point3i(0,0,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+1]))), cv::Point3i(0,size.width,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+2]))), cv::Point3i(size.height,0,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+3]))), cv::Point3i(size.height,size.width,0));
//			temp.Display(0);
//			return temp;
//		}
//		else if (nbands == 8)
//		{
//			STBlock temp(size*BlockSize<int>(2,4,1),CV_8U);
//   			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands]))),cv::Point3i(0,0,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+1]))), cv::Point3i(0,size.width,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+2]))), cv::Point3i(size.height,0,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+3]))), cv::Point3i(size.height,size.width,0));
//      		temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+4]))), cv::Point3i(size.height*2,0,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+5]))), cv::Point3i(size.height*2,size.width,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+6]))), cv::Point3i(size.height*3,0,0));
//			temp.SetBlock(convertTo(normalize(getMagnitude(pyr[n*nbands+7]))), cv::Point3i(size.height*3,size.width,0));
//			temp.Display(0);
//			return temp;
//		}
//		else
//		{
//			//undefined
//			exit(1);
//		}
//	}
//
//}
//
//STBlock Steerable::DisplayLowpass(vector<STBlock>& pyr)
//{
//	
//	STBlock temp( convertTo(normalize(getMagnitude(*(pyr.end()-2)))));
//	temp.Display(0);
//	return temp;
//
//}
//
//STBlock Steerable::DisplayHighpass(vector<STBlock>& pyr)
//{
//	STBlock temp(convertTo(normalize(getMagnitude(*(pyr.end()-1)))));
//	temp.Display(0);
//	return temp;
//}
//
//
//
cv::Mat Steerable::normalize(cv::Mat & A)
{
	cv::Mat temp;
	cv::normalize(A,temp,0,255,NORM_MINMAX);
	return temp;
}
//
//STBlock Steerable::normalize(STBlock & A)
//{
//	STBlock rst(A.size,A.GetDataType());
//	int i=0;
//	for (A.iter=A.GetBeginning(); A.iter!= A.GetEnding(); A.iter++)
//	{
//		rst.SetFrame(i,normalize(*A.iter));
//		i++;
//	}
//	return rst;
//}
//
cv::Mat Steerable::convertTo(cv::Mat &A, int type)
{
	cv::Mat rst;
	A.convertTo(rst,type);
	return rst;
}

//STBlock Steerable::convertTo(STBlock &A, int type)
//{
//	STBlock rst(A.size,type);
//	int i=0;
//	cv::Mat temp;
//	for (A.iter=A.GetBeginning(); A.iter!= A.GetEnding(); A.iter++)
//	{
//		A.iter->convertTo(temp,type);
//		rst.SetFrame(i,temp);
//		i++;
//	}
//	return rst;
//}
//
vector<Steerable::data_type>& Steerable::getSpaceDomainPyr(void)
{
	return this->pyr_space;
}
//
//
cv::Mat Steerable::DFTShift(cv::Mat & A)
{
	cv::Mat rst = cv::Mat::zeros(A.size(),A.type());
	//only works when size is even
	cv::Rect roi1(0,0,A.size().width/2, A.size().height/2);
	cv::Rect roi2(A.size().height/2,0,A.size().width/2, A.size().height/2);
	cv::Rect roi3(0,A.size().width/2,A.size().width/2, A.size().height/2);
	cv::Rect roi4(A.size().height/2, A.size().width/2,A.size().width/2, A.size().height/2);
	A(roi4).copyTo(rst(roi1));
	//mylib::DisplayMat(rst,"rst");
	A(roi3).copyTo(rst(roi2));
	A(roi2).copyTo(rst(roi3));
	A(roi1).copyTo(rst(roi4));
	return rst;

}
//
//STBlock Steerable::DFTShift(STBlock & A)
//{
//	STBlock rst(A);
//	int i=0;
//	for (A.iter = A.GetBeginning(); A.iter != A.GetEnding(); A.iter++)
//	{
//		rst.SetFrame(i, DFTShift(*A.iter));
//		i++;
//
//	}
//	return rst;
//}
}

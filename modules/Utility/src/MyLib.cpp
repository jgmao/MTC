#define  _CRT_SECURE_NO_WARNINGS
#include "MyLib.h"
void mylib::DisplayMat(const cv::Mat &m, const string & name, bool tofileonly)
{
	fstream logfile;
	logfile.open("./temp/"+name+"_Ts.txt",ios::out);
	if (!tofileonly)
		std::cout<<"Display "<<m.size.p[0] <<" "<<m.size.p[1]<<" Matrix:"<<name<<" Data type: "<<m.type()<<endl;
	//logfile<<"Display "<<m.size.p[0] <<" "<<m.size.p[1] <<" Matrix:"<<name<<" Data type: "<<m.type()<<endl;

	if (m.channels()==1)
	{
		for (int i=0; i<m.size.p[0];i++)
		{
			for (int j=0;j<m.size.p[1];j++)
			{
				if (m.type()==5)
				{

					if (!tofileonly)
						std::cout<<m.at<float>(i,j)<<" ";
					logfile<<m.at<float>(i,j)<<" ";
				}
				if (m.type()==CV_8U)
				{
					if (!tofileonly)
					  std::cout<<(int)m.at<uchar>(i,j)<<" ";
					logfile<<(int)m.at<uchar>(i,j)<<" ";
				}
				if (m.type()==CV_64F)
				{
					if (!tofileonly)
						std::cout<<m.at<double>(i,j)<<" ";
					logfile<<m.at<double>(i,j)<<", ";
				}
				if (m.type()==CV_32S)
				{
					if (!tofileonly)
						std::cout<<m.at<int>(i,j)<<" ";
					logfile<<m.at<int>(i,j)<<" ";
				}
				if (m.type()==CV_16S)
				{

					if (!tofileonly)
						std::cout<<m.at<short>(i,j)<<" ";
					logfile<<m.at<short>(i,j)<<" ";
				}
			}
			logfile<<";"<<endl;
			if (!tofileonly)
				std::cout<<";"<<endl;
		}
	}
	else
	{
		vector<cv::Mat> plane;
		cv::split(m,plane);
		for (int c=0; c< m.channels();c++)
		{


			if (!tofileonly)
				std::cout<<"channel "<<c<<" :"<<endl;
			logfile<<"channel "<<c<<" :"<<endl;
			for (int i=0;i<plane[c].size().height;i++)
			{
				for (int j=0;j<plane[c].size().width;j++)
				{
					//int temp= plane[c].type();
					if (plane[c].type()==5)
					{
						if (!tofileonly)
							std::cout<<plane[c].at<float>(i,j)<<"\t";
						logfile<<plane[c].at<float>(i,j)<<"\t";
					}
					else if (plane[c].type()==0)
					{
						if (!tofileonly)
							std::cout<<(int)plane[c].at<uchar>(i,j)<<"\t";
						logfile<<(int)plane[c].at<uchar>(i,j)<<"\t";
					}
					else if (plane[c].type()==CV_64F)
					{

						if (!tofileonly)
							std::cout<<plane[c].at<double>(i,j)<<"\t";
						logfile<<plane[c].at<double>(i,j)<<"\t";
					}
				}
				if (!tofileonly)
					std::cout<<"-----------------------------"<<endl;
				logfile<<";"<<endl;
			}
		}

	}
	logfile.close();
}

cv::Mat mylib::Filter(cv::Mat src, cv::Mat kernel, int flag)
{
	cv::Mat rst;
	cv::filter2D(src,rst,-1,kernel,cv::Point(-1,-1),0.0,flag);
	return rst;
}

//STBlock mylib::FreqComplexFilter(STBlock& src, cv::Mat& kernel)
//{
//	STBlock rst(src.size, src.GetDataType(),src.GetChannel());
//	if ( src.GetChannel() != kernel.channels())
//		cout<<"error, the src and kerenl are not with same channels\n";
//	else if (src.GetChannel()== 2 || src.GetChannel()==6)
//	{
//		cv::Mat temp;
//		int i=0;
//		for (src.iter = src.GetBeginning(); src.iter!= src.GetEnding(); src.iter++)
//		{
//			cv::mulSpectrums(kernel,*src.iter,temp,CV_DXT_ROWS);
//			rst.SetFrame(i, temp);
//			i++;
//		}
//	}
//	else
//		cout<<"error, the input are not complex image\n";
//	return rst;
//}

cv::Mat mylib::FreqComplexFilter(const cv::Mat& src, cv::Mat& kernel, bool conj)
{
	//must make sure the both of src and kernel are complex
	cv::Mat rst(src.size(),src.type());
	if (src.channels()!= kernel.channels())
		cout<<"error, the src and kerenl are not with same channels\n";
	else if (src.channels() == 2 || src.channels() ==6)
		cv::mulSpectrums(kernel,src,rst,CV_DXT_ROWS,conj);
	else
		cout<<"error, the input are not complex image\n";
	return rst;
}

cv::Mat mylib::BinomialKernel(int size, int nDataType)
{
	cv::Mat rst(size,size,nDataType);
	if (nDataType == CV_32F)
	{	
		vector<float> kernel(2,0.5);
		vector<float> temp(2,0.5);
		for (int n=0;n<size-2;n++)
		{
			temp = conv(temp,kernel,0);
		}
		rst = mylib::VecToMat(crossproduct(temp,temp),nDataType);
		return rst;
	}
	else
	{
		vector<double> kernel(2,0.5);
		vector<double> temp(2,0.5);

		for (int n=0;n<size-2;n++)
		{
			temp = conv(temp,kernel,0);
		}
		rst = mylib::VecToMat(crossproduct(temp,temp),nDataType);
	//mylib::DisplayMat(rst,"binomial kernal");	
		return rst;
	}
}
//
//double mylib::ComparePSNR(STBlock & blockA, STBlock &blockB)
//{
//	double temp = blockA.Compare(blockB,0);
//	return 10*log10(std::pow(255.0,2.0)/temp);
//}

cv::Mat mylib::CrossProduct(const cv::Mat &colVec, const cv::Mat &rowVec)
{
	//only workable for depth 1 vectors with type CV_32F now.

	if (colVec.type()!=rowVec.type())
		exit(0);

	cv::Mat result(colVec.size().height,rowVec.size().width,colVec.type());
	for (int i=0; i< colVec.size().height; i++)
		for (int j=0; j< rowVec.size().width; j++)
		{
			if (colVec.type() == CV_32F)
				result.at<float>(i,j) = colVec.at<float>(i) * rowVec.at<float>(j);
			else if (colVec.type() == CV_64F)
				result.at<double>(i,j) = colVec.at<double>(i)* rowVec.at<double>(j);
			else 
				result.at<uchar>(i,j) = colVec.at<uchar>(i) * rowVec.at<uchar>(j);
		}
	return result;
}

double mylib::WeightedSum(const cv::Mat& A, const cv::Mat& W)
{
	if (A.size() != W.size())
		exit(0);
	//int a = A.type();
	//int b = W.type();
	double result=0.0;
	cv::Scalar s;
	vector<Mat> Wv;
	Mat W3;
	for (int i=0; i< A.channels(); i++)
	{
		Wv.push_back(W);
	}
	cv::merge(Wv,W3);
	W3 = A.mul(W3);
	s = cv::sum(W3);
	for (int i=0; i<A.channels(); i++)
		result += s[i];
	return result / A.channels();


}

std::complex<double> mylib::WeightedSumComplex(const cv::Mat & A, const cv::Mat& W)
{
	if (A.size()!=W.size())
		exit(1);
	//int a = A.type();
	//int b = W.type();
	cv::Mat AA;
	if (A.channels() !=2 && A.channels()!=6)
		AA = mylib::toComplex(A);
	else
		AA = A;

	complex<double> result;

	cv::Scalar s;
	vector<Mat> Wv;
	Mat W3;
	for (int i=0; i< AA.channels(); i++)
	{
		Wv.push_back(W);
	}
	cv::merge(Wv,W3);
	W3 = AA.mul(W3);
	s = cv::sum(W3);
	//may be not good for color complex image
	result = complex<double>(s[0],s[1]);
	return result;

}

cv::Mat mylib::GenGaussKer(int size, double sigma, int type)
{
	Mat gaussKernel = cv::getGaussianKernel(size,sigma,type);
	Mat gaussKernel_T;
	cv::transpose(gaussKernel,gaussKernel_T); //stupid
	gaussKernel = mylib::CrossProduct(gaussKernel,gaussKernel_T);
	return gaussKernel;
}


double mylib::AverageScalar3(cv::Scalar s)
{
	return (s[0]+s[1]+s[2])/3;
}

unsigned long mylib::factorial(unsigned long n)
{
	if (n == 1)
		return 1;
	else 
		return n*factorial(n-1);
}


double mylib::factorial(int n)
{
	if (n == 1)
		return 1;
	else 
		return n*factorial(n-1);
}

cv::Mat mylib::complexMul(const cv::Mat& A, const cv::Mat& B)
{
	if (A.channels() != B.channels())
		exit(1);
	else
	{
		cv::Mat AA,BB;
		//int ddd = A.channels();
		if (A.channels()!=2 && A.channels()!=6)
		{
			AA = toComplex(A);
			BB = toComplex(B);
		}
		else
		{
			AA = A;
			BB = B;
		}
		//int a = AA.type();
		//int b = A.type();
		//int c = BB.type();
		//int d = B.type();
		cv::Mat rst;
		cv::mulSpectrums(AA,BB,rst,DFT_ROWS,true);
		return rst;
	}
}


cv::Mat mylib::toComplex(const cv::Mat& A, const cv::Mat& B)
{
	cv::Mat rst;
	vector<cv::Mat> temp;
	//int a = A.channels();
	//int b = B.channels();
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

cv::Mat mylib::toComplex(const cv::Mat& A)
{
	if (A.channels() == 2) //already complex
		return A;
	else
	{
		cv::Mat AZim = cv::Mat::zeros(A.rows,A.cols,A.type());
		return toComplex(A,AZim);
	}
}

cv::Mat mylib::getChannel(const cv::Mat& im, int n)
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

int mylib::combination(int fromNum, int selectNum)
{
	double r; 
	r = factorial(fromNum)/factorial(selectNum)/factorial(fromNum - selectNum);
	return (int)r;
}
//temp disable 20130813 to avoid use Tensor, implement it using Mat later
//#include "Tensor.h"
/*void mylib::CombineImage(const vector<string>& infilenames, string& outfilename, int step)
{
	int count = infilenames.size();
	vector<Tensor<uchar,3>> t(count);
	int maxheight = 0;
	int totalwidth = 0;
	for (int i=0; i< count; i++)
	{
		t[i].Load(infilenames[i]);
		if (t[i].size().height > maxheight)
			maxheight = t[i].size().height;
		totalwidth+= t[i].size().width;
	}
	Tensor<uchar,3> r(Size3(maxheight,totalwidth + (count-1)*step,1));
	for (int i=0; i< count; i++)
		r.SetBlock(Point3i(0,i*step+i*t[i].size().width,0),t[i]);
	r.SaveBlock(outfilename);
	r.Display(3000,1);

}
*/
cv::Mat mylib::readMatFromTxt(string filename, int H, int W)
{
  CV_DbgAssert(H!=0&&W!=0);
  cv::Mat rst(H,W,CV_64F);
  ifstream fp(filename,ios::in);
  double temp;
  int count=0;
  while(count<H*W)
  {
    fp>>temp;
    rst.at<double>(count/W,count%W)=temp;
    count++;
  }
  return rst;
}

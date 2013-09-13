#include "Lighting.h"
namespace lighting{
Lighting::Lighting()
{


}



Mat Lighting::LSFitting(const Mat& im, int order) const
{
	if ( order > 2 || order < 1)
	{
		CV_Error(CV_StsNotImplemented,"unsupport");
	}
	const unsigned int cn = im.channels();
	CV_Assert(cn==1);
	Tensor<double,1> ts(im);
	if (order == 1)
	{
		int count = 0;
		Tensor<double,1> H(Size3(im.size().area(),3,1));
		Tensor<double,1> t(Size3(im.size().area(),1,1));
		Tensor<double,1> rst(Size3(3,1,1));
		for (int i=0; i< im.size().height; i++)
			for (int j=0; j < im.size().width; j++)
			{
				H(count,0,0) = saturate_cast<double>(i);
				H(count,1,0) = saturate_cast<double>(j);
				H(count,2,0) = 1;
				t(count,0,0) = ts(i,j,0);
				count++;
			}
			for (unsigned int i=0; i< cn; i++)
			{
				rst[i] = (H[i].t() * H[i]).inv()*(H[i].t()*t[i]);
			}
			return rst;
	}
	else
	{
		Tensor<double,1> rst(Size3(6,1,1));
		Tensor<double,1> H(Size3(im.size().area(),6,1));
		Tensor<double,1> t(Size3(im.size().area(),1,1));
		int count =0;
		for (int i=0; i< im.size().height; i++)
			for (int j=0; j< im.size().width; j++)
			{
				H(count,0,0)=saturate_cast<double>(i*i);
				H(count,1,0)=saturate_cast<double>(i*j);
				H(count,2,0)=saturate_cast<double>(j*j);
				H(count,3,0)=saturate_cast<double>(i);
				H(count,4,0)=saturate_cast<double>(j);
				H(count,5,0)=1;
				t(count,0,0)=ts(i,j,0);
				count++;
			}
			for (unsigned int i=0; i< cn; i++)
			{
				rst[i] = (H[i].t() * H[i]).inv()*(H[i].t()*t[i]);
			}
			return rst;
	}
}

Mat Lighting::BuildLightingPlane(const Mat& im, const Mat& param, int order) const
{
	if ( order > 2 || order < 1)
	{
		CV_Error(CV_StsNotImplemented,"unsupport");
	}
	CV_Assert(im.channels()==1);
	Tensor<double,1> rst;
	Tensor<double,1> p(param);
	//param.Print();
		for (int i=0; i<im.size().height; i++)
			for (int j=0; j< im.size().width; j++)
			{
				if(order==1)
					rst(i,j,0)=(p(0,0,0)*i+p(1,0,0)*j+p(2,0,0))(0);
				else
					rst(i,j,0)=(p(0,0,0)*i*i+p(1,0,0)*i*j+p(2,0,0)*j*j+p(3,0,0)*i+p(4,0,0)*j+p(5,0,0))(0);
			}
			return rst;
}

Mat Lighting::LightingCorrection(Mat& changeFrom, const Mat& changeTo,bool saveCodeLength)
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

  CV_Assert(changeFrom.channels()==changeTo.channels());
  CV_Assert(changeFrom.channels()==1);
  unsigned int cn = changeFrom.channels();
  Tensor<double,1> lighting(changeTo);
  //Tensor<double,cn> org(changeFrom);
	//lighting.Print();
	///quantize the ligting surface by dct
	//refresh lighting array
	//lightingDCTCoeffStat.clear();
	int lf = 0;//number of LF
	const float* tbl;
	if (lighting.size().height<=16)
	{
		tbl = base8;
		lf = 64;
	}
	else if (lighting.size().height<=32)
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

	Tensor<double,1> flit = Tensor<double,1>(changeTo); //target X
	Tensor<double,1> fcand = Tensor<double,1>(changeFrom); //candidate Y
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
				for (unsigned int cc=0; cc< cn; cc++)
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
							if (lighting.size().height<32)
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
  lighting = Tensor<double,1>(flit);

	//this->SetBlock(Tensor<T,cn>(Tensor<double,cn>(*this) - fcand + flit + 128));
	//this->SetBlock(Tensor<T,cn>(fcand + flit + 128)); //target lighting + (target - target lighting)
  changeFrom = changeFrom + flit;
        //SetBlock(Tensor<double,cn>(changeFrom) + flit);
	return lighting;
	//this->SetBlock((*this) + BuildLightingPlane(tempTagLighting - tempCanLighting));
}

Mat Lighting::LightingCorrection(Mat& changeFrom, const Mat& changeTo, const Tensor<double,1>& VQCodebook)
{
	//do this after quilting
	//change to is the original block (the lighting should be)
	Tensor<double,1> tempTagLighting = LSFitting(changeFrom);
	Tensor<double,1> tempCanLighting = LSFitting(changeTo);
	for (int i=0; i< tempTagLighting.size().height; i++)
		this->lightTag.push_back(tempTagLighting(i,0,0));
	for (int j=0; j< tempCanLighting.size().height; j++)
		this->lightCan.push_back(tempCanLighting(j,0,0));
	//tempTagLighting.Print();
	//SearchCodeword(tempTagLighting,VQCodebook).Print();
	//tempCanLighting.Print();
	//SearchCodeword(tempCanLighting,VQCodebook).Print();
	//search for vq codebook and use the quantized version to recover lighting.
	Tensor<double,1> lighting = BuildLightingPlane(changeFrom, SearchCodeword(tempTagLighting,VQCodebook) - SearchCodeword(tempCanLighting,VQCodebook));
	//BuildLightingPlane(tempTagLighting-tempCanLighting).Print();
	//lighting.Print();
	changeFrom = changeFrom + lighting;
	//this->SetBlock((*this) + lighting);
	return lighting;
}

Mat Lighting::SearchCodeword(const Mat& val, const Mat& VQCodeBook)
{
	int label=-1;
	double dist = std::numeric_limits<double>::max();
	const int cn = val.channels();
	CV_Assert(cn==1);
	double temp_dist = 0;
	Vec<double,1> dist_vec;
	if(VQCodeBook.size().width!= val.size().height)
	{
		CV_Error(CV_StsUnmatchedSizes,"the input polynomial has different size as codebook dim");
	}

	Tensor<double,1> tVQCodeBook(VQCodeBook);
	Tensor<double,1> tVal(val);

	for (int i=0; i< tVQCodeBook.size().height; i++)
	{

		dist_vec =Vec<double, 1>::all(0);
		temp_dist = 0;
		for(int j=0; j< tVQCodeBook.size().width; j++)
		{
		    Vec<double,1> temp  = tVQCodeBook(i,j,0) - tVal(j,0,0);
		    dist_vec+= mylib::VecPow<double,1>(temp,2.0);

		}
		for (int c=0; c<cn; c++)
			temp_dist+= dist_vec[c];
		if (temp_dist <= dist)
		{
			dist = temp_dist;
			label = i;
		}
	}
	return Tensor<double,1>(tVQCodeBook[0].row(label)).Transpose();
}




int Lighting::GetLightingCodeLength(void) const
{
	return this->codeLength;
}
void Lighting::SetLightingCodeLength(int l)
{
	this->codeLength=l;
}




void Lighting::RecordLighting(void)
{
  fstream lightfile;
  lightfile.open("./lighting.txt",ios::app);

  for (unsigned int i=0; i< lightingDCTCoeffStat.size(); i++)
    lightfile<<lightingDCTCoeffStat[i]<<endl;
  lightfile.close();
}

vector<Vec<double,1> > Lighting::GetTagLighting() const
{
  return this->lightTag;
}
vector<Vec<double,1> > Lighting::GetCanLighting() const
{
  return this->lightCan;
}

Mat ComputeTPSS(const Mat& im, double p)
{
  CV_Assert(im.channels()==1);//only implement gray scale

  Tensor<double,1> G(im.clone());
  ThinPlateSpline tps;
  for (int t=0; t<G.size().depth; t++)
  {
    Mat temp = G.GetFrame(t);
    tps.load(temp);
    G.SetFrame(t,tps.solve(p));
  }
  return G;
}

Mat ComputeSAT(const Mat& im)
{
  //2D only
  CV_Assert(im.channels()==1);
  Tensor<double,1> temp(im);
  Tensor<double,1> S = temp.ExtendHalfBoundary();

  Point3i pos(im.size().height-1, im.size().width-1,0);

  S = ComputeSAT(S,Point3i(0,0,0),pos);
//  for ( int i=1; i< this->tsSize.height+1; i++)
//    for (int j=1; j<this->tsSize.width+1; j++)
//    {
//      S(i,j,0)=S(i-1,j,0)+S(i,j-1,0)-S(i-1,j-1,0)+S(i,j,0);
//    }
  return S;
}

Mat ComputeSAT(const Mat& S, const Point3i& sPos, const Point3i& ePos)
{
  //remember S is extended boundary by 1 to the left and up (zeros)
  Mat SAT=S;
  for (int i = sPos.x+1; i <= ePos.x+1; i++)
    for (int j=sPos.y+1; j <= ePos.y +1; j++)
    {
      SAT.at<double>(i,j)=SAT.at<double>(i-1,j)+SAT.at<double>(i,j-1)-SAT.at<double>(i-1,j-1)+SAT.at<double>(i,j);
    }
  return SAT;

}
}



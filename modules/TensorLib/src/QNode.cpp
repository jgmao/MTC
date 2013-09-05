#include "QNode.h"

////////////////////////////// Implementation ///////////////////////////
namespace tensor{
template< class T, size_t cn>
QNode<T,cn>::~QNode(void)
{

}
template< class T, size_t cn>
QNode<T,cn>::QNode(void):Tensor<T,cn>()
{
	overlapSize= Size3();
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING;
	qsize = 0;
}
template< class T, size_t cn>
QNode<T,cn>::QNode(const cv::Mat &mt):Tensor<T,cn>(mt)
{
	//overlapSize= Size3();
	bits=0;
	predicted_method=CodingMethodNames::CODING_JPEG;
	qsize = 0;
	footPos = Point3i(this->size().height-1,this->size().width-1,this->size().depth-1);
	nodeFoot = this->operator[](footPos);
	overlapSize = this->size()/4;//by default
}
template< class T, size_t cn>
QNode<T,cn>::QNode(const QNode& nd):Tensor<T,cn>(nd)
{
	overlapSize = nd.overlapSize;
	bits = nd.bits;
	predicted_method = nd.predicted_method;
	footError = nd.footError;
	nodeFoot = nd.nodeFoot;
	approxFoot = nd.approxFoot;
	footPos = nd.footPos;
	quanBit = nd.quanBit;
	boundary = nd.boundary;
	qsize = nd.qsize;
	leftBound = nd.leftBound;
	upBound = nd.upBound;
	lowBound = nd.lowBound;
	rightBound = nd.rightBound;
	rightBound = nd.rightBound;
	bounds.clear();
	bounds.push_back(&upBound);
	bounds.push_back(&leftBound);
  leftFoot = nd.leftFoot;
  upFoot = nd.upFoot;
  eNode = nd.eNode;
  this->feet = nd.feet;
}
/*
template< class T, size_t cn>
QNode<T,cn>::QNode(const QTree& qt):QTree<T,cn>(qt)
{
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING
	qsize = 0;
	footPos = Point3i(size().height-1,size().width-1,size().depth-1);
	nodeFoot = this->operator[](footPos);
	overlapSize = size()/4;//by default
}
*/
template< class T, size_t cn>
QNode<T,cn>::QNode(const Size3& sz, const Vec<T,cn>& val):Tensor<T,cn>(sz,val)
{
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING;
	qsize = 0;
	footPos = Point3i(this->size().height-1,this->size().width-1,this->size().depth-1);
	nodeFoot = this->operator[](footPos);
	overlapSize = this->size()/4;//by default
}
/*
template< class T, size_t cn>
QNode<T,cn>::QNode(const QTree& qt, const Size3& overlapSize):QTree<T,cn>(qt)
{
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING
	qsize = 0;
	footPos = Point3i(size().height-1,size().width-1,size().depth-1);
	nodeFoot = this->operator[](footPos);
	//overlapSize = size()/4;
}
*/
template< class T, size_t cn>
QNode<T,cn>::QNode(const Size3& sz, const Size3& overlapSize):Tensor<T,cn>(sz)
{
	
	this->overlapSize = overlapSize;
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING;
	qsize = 0;
	footPos = Point3i(this->size().height-1,this->size().width-1,this->size().depth-1);
	nodeFoot = this->operator[](footPos);
}

template< class T, size_t cn>
QNode<T,cn>::QNode(const Size3& sz, double overlapRatio):Tensor<T,cn>(sz)
{
		
	this->overlapSize.height = int((double)sz.height*overlapRatio);
	this->overlapSize.width = int((double)sz.width*overlapRatio);
	this->overlapSize.depth = int((double)sz.depth*overlapRatio);
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING;
	qsize = 0;
	footPos = Point3i(this->size().height-1,this->size().width-1,this->size().depth-1);
	nodeFoot = this->operator[](footPos);
}

template< class T, size_t cn>
QNode<T,cn>& QNode<T,cn>::operator= (const QNode& nd)
{
	//deep copy
	Tensor<T,cn>::operator=(nd);
	overlapSize = nd.overlapSize;
	bits = nd.bits;
	predicted_method = nd.predicted_method;
	footError = nd.footError;
	nodeFoot = nd.nodeFoot;
	approxFoot = nd.approxFoot;
	footPos = nd.footPos;
	quanBit = nd.quanBit;
	boundary = nd.boundary;
	leftBound = nd.leftBound;
	upBound = nd.upBound;
	lowBound = nd.lowBound;
	rightBound = nd.rightBound;
	bounds.clear();
	bounds.push_back(&upBound);
	bounds.push_back(&leftBound);
	qsize = nd.qsize;
  leftFoot = nd.leftFoot;
  upFoot = nd.upFoot;
  eNode = nd.eNode;
  this->feet = nd.feet;
	return *this;
}
template< class T, size_t cn>
QNode<T,cn>::QNode(const Tensor<T,cn>& ts, const Size3& sz, const Point3i& pos, const Size3& boundarySize)
{

	Size3 bsize = boundarySize;
	if (boundarySize.height == 0)
		bsize.height = 1;
	if (boundarySize.width == 0)
		bsize.width = 1;
	Cube roi(pos,sz);
	ts.Ref(roi,*this);
	//*this = QTree(ts.Crop(pos,sz));
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING;
	qsize = 0;
	footPos = this->tsSize.Point3() - Point3i(1,1,1);//Point3i(size().height-1,size().width-1,size().depth-1);
	nodeFoot = ts[footPos+roi.offset()];
	overlapSize = bsize;
	leftBound = Tensor<T,cn>(Size3(bsize.height+sz.height,bsize.width,sz.depth));
	upBound = Tensor<T,cn>(Size3(bsize.height,bsize.width+sz.width,sz.depth));
	lowBound = Tensor<T,cn>(Size3(bsize.height,bsize.width+sz.width,sz.depth));
	rightBound = Tensor<T,cn>(Size3(bsize.height+sz.height,bsize.width,sz.depth));
  upFoot = QNode::value_type::all(0);
  leftFoot = QNode::value_type::all(0);
  Size3 extSize = sz;
  Point3i extPos=pos;
	if (pos.x !=0 && pos.y==0)
	{
		ts.Ref(Cube(Point3i(pos.x-bsize.height, pos.y, pos.z - bsize.depth),Size3(bsize.height,sz.width,sz.depth)),upBound);
		//leftBound = Tensor<T,cn>(
		//boundary.push_back(ts.Crop( Point3i(pos.x-bsize.height, pos.y, pos.z - bsize.depth),
		//						Size3(bsize.height,sz.width,sz.depth)));
    extPos.x = pos.x-bsize.height;
    extPos.y = 0;
    extSize.height+=bsize.height;
	}
	else if (pos.y !=0 && pos.x==0)
	{
		ts.Ref(Cube(Point3i(pos.x, pos.y-bsize.width,pos.z - bsize.depth), Size3(sz.height,bsize.width,sz.depth)),leftBound);
		//boundary.push_back(ts.Crop( Point3i(pos.x, pos.y-bsize.width,pos.z - bsize.depth),
		//						Size3(sz.height,bsize.width,sz.depth)));
    extPos.x=0;
    extPos.y = pos.y-bsize.width;
    extSize.width+=bsize.width;
	}
	else if (pos.x!=0 && pos.y!=0)
	{
		ts.Ref(Cube(Point3i(pos.x-bsize.height, pos.y-bsize.width,pos.z - bsize.depth), Size3(bsize.height,sz.width+bsize.width,sz.depth)),upBound);
		ts.Ref(Cube(Point3i(pos.x-bsize.height, pos.y-bsize.width,pos.z - bsize.depth),	Size3(sz.height+bsize.height,bsize.width,sz.depth)),leftBound);	
	  extPos.x=pos.x - bsize.height;
    extPos.y=pos.y - bsize.width;
    extSize.height+=bsize.height;
    extSize.width+=bsize.width;
  }

	if (pos.x + sz.height+ bsize.height -1 < ts.size().height) // have lower bound
	{
		if (pos.y-bsize.width < 0)
			ts.Ref(Cube(pos.x+sz.height,0, pos.z - bsize.depth,bsize.height, sz.width+bsize.width-pos.y,sz.depth),lowBound);
		else if (pos.y+sz.width+bsize.width > ts.size().width)
			ts.Ref(Cube(pos.x+sz.height,pos.y-bsize.width,0,bsize.height,sz.width+bsize.width,sz.depth),lowBound);
		else
			ts.Ref(Cube(pos.x+sz.height,pos.y-bsize.width,0,bsize.height,sz.width+2*bsize.width,sz.depth),lowBound);
    extSize.height+=bsize.height;
	}
	if (pos.y + sz.width + bsize.width -1 < ts.size().width)
	{
		if (pos.x - bsize.height < 0)
			ts.Ref(Cube(0, pos.y + sz.width, pos.z - bsize.depth, sz.height+bsize.height, bsize.width, sz.depth),rightBound);
		else if (pos.x+sz.height+bsize.height > ts.size().height)
			ts.Ref(Cube(pos.x - bsize.height, pos.y + sz.width, pos.z - bsize.depth, sz.height+bsize.height, bsize.width, sz.depth),rightBound);
		else
			ts.Ref(Cube(pos.x - bsize.height, pos.y + sz.width, pos.z - bsize.depth, sz.height+2*bsize.height, bsize.width, sz.depth),rightBound);
    extSize.width+=bsize.width;
	}
	bounds.push_back(&upBound);
	bounds.push_back(&leftBound);
  ts.Ref(Cube(extPos,extSize),eNode);
}
template< class T, size_t cn>
Size3 QNode<T,cn>::overlap(void) const
{
	return overlapSize;
}

template< class T, size_t cn>
QNode<T,cn>::QNode(string cFileName):Tensor<T,cn>(cFileName)
{
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING;
	qsize = 0;
	footPos = Point3i(this->size().height-1,this->size().width-1,this->size().depth-1);
	nodeFoot = this->operator[](footPos);
	overlapSize = this->size()/4;//by default
}
template< class T, size_t cn>
void QNode<T,cn>::SetOverlapSize(const Size3& overlapSize)
{

	this->overlapSize = overlapSize;

}
////////////////// Interpolation ///////////////////////
template< class T, size_t cn>
int QNode<T,cn>::BitLength(void) const
{
	return (int)quanBit.length();
}

template<class T,size_t cn>
string QNode<T,cn>::GetBits(void) const
{
	return quanBit;
}

template< class T, size_t cn>
QNode<T,cn>& QNode<T,cn>::LinearInterp(double qsize, bool computefoot, Directions direction)
{
	this->qsize = qsize;
	bits=0;
	if (this->size().height!= this->size().width)
	{
		//unsupport
		exit(-1);
	}
	if (this->size().depth>1)
	{
		cout<<"error, Node need to be 2D image at this time\n";
		exit(1);
	}
	//this->qsize = qsize;
	InitBound(qsize);//initail boundary when reach image bound
	if (computefoot)
		ComputeFoot(qsize, direction);
	InterpLastCol();
  //eNode.Display();
  //Print();
//	upBound.Print();
	InterpLastRow();
  //eNode.Display();
//	Print();
//	leftBound.Print();
	InterpCenter();
 // eNode.Display();
//  Print();
	this->predicted_method = CodingMethodNames::CODING_PQI;
	return *this;
}
template< class T, size_t cn>
QNode<T,cn>& QNode<T,cn>::LinearInterpGeneric(double qsize)
{
	//this->qsize = qsize;
	//InitBound(qsize);
	//Vec<T,cn> foot = this->operator()(this->GetFootPos());

	//Vec<T,cn> upPixel = this->GetBoundaryUp(this->GetBoundaryUp().size().Point3()-Point3i(1,1,1));
	//Vec<T,cn> leftPixel = this->GetBoundaryLeft(this->GetBoundaryLeft().size().Point3()-Point3i(1,1,1));


	return *this;

}

template<class T, size_t cn>
void QNode<T,cn>::InitBound(double qsize)
{
	if (this->offset().x==0 && this->offset().y==0)
	{
		//make up the very first node
		MakeVeryFirstNode(qsize);
	}
	else if (this->offset().x==0 && this->offset().y!=0)
	{
		//make the 1st row
		//boundary[0].Print();
		InterpZeroRow();
	//	extNode.GetBoundaryUp().Print();
	//	extNode.Print();
	}
	else if (this->offset().y==0 && this->offset().x!=0)
	{
		//boundary[0].Print();
		InterpZeroCol();
		//make the 1st column
		//extNode.GetBoundaryLeft().Print();
		//extNode.Print();
	}
}
template<class T, size_t cn>
Vec<T,cn> QNode<T,cn>::ComputePLCFoot(const Tensor<T,cn>& ref, const Point3i& pos, const Size3& sz, int caseNo, int method)
{
	Tensor<T,cn> temp;
	Point3i posnew;
	Size3 sznew;
	switch (caseNo) 
	{
	case 1://same pos, sz as block
		posnew = pos;
		sznew = sz;
		break;
	case 2: //pos at center of block, same sz as block
		posnew = pos+ Point3i(sz.height/2,sz.width/2,sz.depth/2);
		sznew = sz;
		break;
	case 3: //pos at center of block, half size as block
		posnew = pos + Point3i(sz.height/2,sz.width/2,sz.depth/2);
		sznew = Size3(sz.height/2,sz.width/2,sz.depth);
		break;
	default:
		cout<<"wrong case No : must from 1 to 3"<<endl; 
		CV_Assert(false);
	}

	//deal with exceed border problem 
	if (posnew.x<0)
	{
		sznew.height-=posnew.x;
		posnew.x=0;
	}
	if (posnew.y<0)
	{
		sznew.width-=posnew.y;
		posnew.y=0;
	}
	if (posnew.x+sznew.height>ref.size().height)
	{
		sznew.height= ref.size().height-posnew.x;
	}
	if (posnew.y+sznew.width>ref.size().width)
	{
		sznew.width = ref.size().width - posnew.y;
	}

	ref.Ref(Cube(posnew,sznew),temp);
	Vec<T,cn> rst;
  //i//
	//temp.Print();
	switch (method)
	{
		case 1://mean
			rst = temp.Mean();
			break;
		case 2://lower right triagle mean
			rst = temp.LRMean();
			break;
		case 3: //percentage mean
			rst = temp.PercentageMean(0.4,0.6);
			break;
		case 4: //median
			rst = temp.PercentageMean(0.5,0.5);
			break;
		case 5: //upper left triagnle mean
			rst = temp.URMean();
			break;
		default:
			cout<<"wrong method No, must from 1 to 4"<<endl;
	}
	//Vec<T,cn> rst = temp.PercentageMean(0.4,0.6);
	//Vec<T,cn> rst = temp.Mean();
	//Vec<T,cn> rst = temp.LRMean();
	return rst;

}
template<class T, size_t cn>
void QNode<T,cn>::ComputeFoot(double qsize, Directions direction, SrcCodingMethod coding_method)
{
	Point3i bPos = overlapSize.Point3() - Point3i(1,1,0);
	Point3i tempFootPos = footPos;
	Vec<double,cn> PredictedFoot;
	// i think this is not necessary anymore
	//if (offset() == Point3i(0,0,0))
	//	PredictedFoot = Vec<T,cn>::all(127);
	//else
	if (direction == Directions::DIRECTION_VERTICAL)
	{
		PredictedFoot = (Vec<double,cn>::all(0.5).mul(upBound[bPos+Point3i(0,this->tsSize.width,0)] + this->operator[](footPos)));//now foot is at center of last col
		tempFootPos = footPos - Point3i(this->size().height/2,0,0);
	}
	else if (direction == Directions::DIRECTION_HORIZONTAL)
	{
		PredictedFoot = (Vec<double,cn>::all(0.5).mul(this->operator[](footPos) + leftBound[bPos+Point3i(this->tsSize.height,0,0)])); // foot is at center of last col
		tempFootPos = footPos - Point3i(0,this->size().width/2,0);
	}
	else if (direction == Directions::DIRECTION_CENTER)
	{
		PredictedFoot = (Vec<double,cn>::all(0.25).mul(upBound[bPos+Point3i(0,this->tsSize.width/2,0)] +
															this->operator[](footPos-Point3i(0,this->tsSize.width/2,0))+
															this->operator[](footPos-Point3i(this->tsSize.height/2,0,0))+
															leftBound[bPos+Point3i(this->tsSize.height/2,0,0)]));

		tempFootPos = footPos - Point3i(this->size().height/2,this->size().width/2,0);
	}
	else
	{
		PredictedFoot = (Vec<double,cn>::all(0.5).mul(upBound[bPos+Point3i(0,this->tsSize.width,0)] + leftBound[bPos+Point3i(this->tsSize.height,0,0)]));
	}
	//Print();
	footError = Vec<double,cn>(this->operator[](tempFootPos)) - PredictedFoot;
	approxFoot = Quantize(footError,qsize,coding_method) + PredictedFoot;
	(*this)[tempFootPos] = approxFoot; // copy foot val 
}



template<class T, size_t cn>
void QNode<T,cn>::MakeVeryFirstNode(double qsize) 
{
	//Vec<T,cn> PredictedFoot = Vec<T,cn>::all(127);
	//footError = nodeFoot - PredictedFoot;
	//approxFoot = Quantize(footError,qsize) + PredictedFoot;
  Size3 upSize,leftSize;
  if (upBound.size()==Size3(0,0,0))
  {
    upSize = this->overlap()+Size3(0,this->tsSize.width,1);
  }
  else
    upSize = upBound.size();
  if (leftBound.size()==Size3(0,0,0))
    leftSize = this->overlap()+Size3(this->tsSize.height,0,1);
  else
    leftSize = leftBound.size();
	upBound = Tensor<T,cn>(upSize,Vec<T,cn>::all(127));
	leftBound = Tensor<T,cn>(leftSize,Vec<T,cn>::all(127));
	//(*this)[footPos] = approxFoot; // copy foot val to the SE of tempNode
}


template<class T, size_t cn>
void QNode<T,cn>::InterpZeroRow(void)
{

	Tensor<T,cn> tempBound(leftBound.size()+Size3(overlapSize.height,0,0));
	tempBound.SetBlock(Point3i(overlapSize.height,0,0),leftBound);
	leftBound = tempBound;
	Point3i bPos = overlapSize.Point3() - Point3i(1,1,0);
	leftBound[bPos] = leftBound[bPos+Point3i(1,0,0)];
  //if (upBound.size()==Size3(0,0,0))
    upBound = Tensor<T,cn>(overlap()+Size3(0,this->tsSize.width,1));
	upBound[bPos] = leftBound[bPos+Point3i(1,0,0)];
	upBound[bPos+Point3i(0,this->tsSize.width,0)] = leftBound[bPos+Point3i(this->tsSize.height,0,0)];
	for(int z=0; z < this->size().depth; z++)
	{
		int width = this->tsSize.width+1;
		int stepY = width>>1;
		for (; stepY>0; stepY=stepY>>1)
		{
			for (int i=1; i*stepY<width-1; i++)
			{
				if (i%2==0)
					continue;
				else
				{
					upBound[bPos+Point3i(0,i*stepY,z)] = Vec<double,cn>::all(0.5).mul(upBound[bPos+Point3i(0,(i-1)*stepY,z)]+ upBound[bPos+Point3i(0,(i+1)*stepY,z)]);
				}
			}
		}
	}
	return;
}

template<class T, size_t cn>
void QNode<T,cn>::InterpZeroCol(void)
{
	
	Tensor<T,cn> tempBound(upBound.size()+Size3(0,overlapSize.width,1));
	tempBound.SetBlock(Point3i(0,overlapSize.width,0),upBound);
	upBound = tempBound;
	Point3i bPos = overlapSize.Point3()-Point3i(1,1,0);
	upBound[bPos] = upBound[bPos+Point3i(0,1,0)];
  //if (leftBound.size()==Size3(0,0,0))
    leftBound = Tensor<T,cn>(overlap()+Size3(this->tsSize.height,0,1));
	leftBound[bPos]=upBound[bPos+Point3i(0,1,0)];
	leftBound[bPos+Point3i(this->tsSize.height,0,0)] = upBound[bPos+Point3i(0,this->tsSize.width,0)];
	for(int z=0; z< this->size().depth; z++)
	{
		int height = this->tsSize.height+1;
		int stepX = height>>1;
		for (; stepX>0; stepX=stepX>>1)
		{
			for (int i=1; i*stepX<height-1; i++)
			{
				if (i%2==0)
					continue;
				else
				{
					leftBound[bPos+ Point3i(i*stepX,0,z)] = Vec<double,cn>::all(0.5).mul
						(leftBound[bPos+Point3i((i-1)*stepX,0,z)]+
						leftBound[bPos+Point3i((i+1)*stepX,0,z)]);
				}
			}		
		}
	}
	return;
}

template<class T, size_t cn>
bool QNode<T,cn>::AssertExtend(const QNode& extNode,int n) const
{
	if (extNode.size() - Size3(n,n,0) == this->size())
		return true;
	else
		return false;
}
template<class T, size_t cn>
void QNode<T,cn>::InterpLastRow(void)
{
	Tensor<T,cn> extRow(Size3(1,this->tsSize.width+1,this->tsSize.depth));
	Point3i bPos = overlapSize.Point3() - Point3i(1,1,0);
	extRow(0,0,0)= leftBound[bPos + Point3i(this->tsSize.height,0,0)];
	extRow(0,this->tsSize.width,0) = this->operator[](footPos);//approxFoot; modify gjin Sep 15, 2011
	for(int z=0; z < this->size().depth; z++)
	{
		int width = extRow.size().width;
		int stepY = width>>1;
		for (; stepY>0; stepY=stepY>>1)
		{
			for (int i=1; i*stepY<width-1; i++)
			{
				if (i%2==0)
					continue;
				else
				{
					extRow(0,i*stepY,z) = Vec<double,cn>::all(0.5).mul(extRow(0,(i-1)*stepY,z)+extRow(0,(i+1)*stepY,z));
				}
			}
		}
	}
	this->SetBlock(Point3i(this->tsSize.height-1,0,0),extRow.Crop(Point3i(0,1,0),Size3(1,this->tsSize.width,1)));
	return;
}

template<class T, size_t cn>
void QNode<T,cn>::InterpLastCol(void)
{
	Tensor<T,cn> extCol(Size3(this->tsSize.height+1,1,this->tsSize.depth));
	Point3i bPos = overlapSize.Point3() - Point3i(1,1,0);
	extCol(0,0,0)= upBound[bPos + Point3i(0,this->tsSize.width,0)];
	extCol(extCol.size().height-1,0,0) = this->operator[](footPos);//approxFoot; modify gjin Sep 15, 2011
	for(int z=0; z< this->size().depth; z++)
	{
		int height = extCol.size().height;
		int stepX = height>>1;
		for (; stepX>0; stepX=stepX>>1)
		{
			for (int i=1; i*stepX<height-1; i++)
			{
				if (i%2==0)
					continue;
				else
				{
					extCol(i*stepX,0,z) = Vec<double,cn>::all(0.5).mul(extCol((i-1)*stepX,0,z)+extCol((i+1)*stepX,0,z));
				}
			}
		}
	}
  //extCol.Print();
        this->SetBlock(Point3i(0,this->tsSize.width-1,0),extCol.Crop(Point3i(1,0,0),Size3(this->tsSize.height,1,1)));
  //this->Print();
  //eNode.Print();
}
template<class T, size_t cn>
void QNode<T,cn>::InterpCenter(void)
{
	Tensor<T,cn> extNode = this->ExtendHalfBoundary();
	Point3i bPos = overlapSize.Point3() - Point3i(1,1,0);
	
	extNode.SetBlock(upBound.Crop(bPos,Size3(1,extNode.size().width,1)));
	extNode.SetBlock(leftBound.Crop(bPos,Size3(extNode.size().height,1,1)));
	//extNode.Display(0);

	for (int z=0; z< this->size().depth; z++)
	{
		int index=0;
		int height = extNode.size().height;//-qNode.overlapSize.height;
		int width = extNode.size().width;//-qNode.overlapSize.width;
		int stepX = height>>1;
		int stepY = width>>1;
		while (stepX>0 && stepY>0)
		{
			if (index%2==0)
			{
				for (int i=1; i*stepX<height-1; i++)
					for (int j=1; j*stepY<width-1; j++)
					{
						if ((i+j)%2==1 || index==0)	
						{
							extNode[Point3i(i*stepX, j*stepY,z)] = Vec<double,cn>::all(0.25).mul 
								(extNode[Point3i(i*stepX, (j-1)*stepY,z)]+
								 extNode[Point3i(i*stepX, (j+1)*stepY,z)]+
								 extNode[Point3i((i+1)*stepX, j*stepY,z)]+
								 extNode[Point3i((i-1)*stepX, j*stepY,z)]);
#ifdef INT_RST
				extNode[Point3i(i*stepX, j*stepY,z)]  = Vec<T,cn>(Vec<int,cn>(extNode[Point3i(i*stepX, j*stepY,z)] ));
#endif
						}
						else
							continue;
					}
				stepX=stepX>>1;
				stepY=stepY>>1;
			}
			else
			{
				for (int i=1; i*stepX<height-1; i++)
					for (int j=1; j*stepY<width-1; j++)
					{
						if (i%2==0||j%2==0)
							continue;
						else
						{
							extNode[Point3i(i*stepX, j*stepY,z)] = Vec<double,cn>::all(0.25).mul
								(extNode[Point3i((i+1)*stepX, (j-1)*stepY,z)]+ 
								 extNode[Point3i((i-1)*stepX, (j-1)*stepY,z)]+
								 extNode[Point3i((i+1)*stepX, (j+1)*stepY,z)]+
								 extNode[Point3i((i-1)*stepX, (j+1)*stepY,z)]);
#ifdef INT_RST
				extNode[Point3i(i*stepX, j*stepY,z)]  = Vec<T,cn>(Vec<int,cn>(extNode[Point3i(i*stepX, j*stepY,z)] ));
#endif
		
						}
					}
			}
			index++;
		}
	}
	this->SetBlock(extNode.Crop(Point3i(1,1,0),this->tsSize));
}
template<class T, size_t cn>
Point3i QNode<T,cn>::GetFootPos(void) const
{
	return this->offset()  + this->footPos;//this->footPos;
}
template<class T, size_t cn>
Point3i QNode<T,cn>::GetFootPos(void) 
{
	return this->offset()  + this->footPos;//this->footPos;
}
template<class T, size_t cn>
Vec<double,cn> QNode<T,cn>::Quantize(const Vec<double,cn>& val, double qsize, SrcCodingMethod coding_method)
{

	QNode<T,cn>::value_type rst;
	int index;
	
	for (int i=0; i< cn; i++)
	{
		index = (int)floor( (val[i]/qsize) + 0.5);
		if (coding_method == SrcCodingMethod::UNARY_CODE)
			bits+= UnaryCoding(index);
		else
			bits += VariableLengthCoding(index); 
		
#ifdef INT_QUAN
		rst[i] =T(int(T(index)*qsize)); //the quantizer output only int values
#else
		rst[i] = T((double)index*qsize);
#endif
	}
	return rst;
}


template <class T, size_t cn>
Vec<int,cn> QNode<T,cn>::UniformQuantize(const typename QNode<T,cn>::value_type& val, double qsize)
{
	Vec<int,cn> rst;
	for (int i=0; i<cn;i++)
	{
		rst[i]=(int)floor((val[i]/qsize)+0.5);
	}
	return rst;
}
template<class T, size_t cn>
int QNode<T,cn>::VariableLengthCoding(int index)
{
	const int MaxLevel = 50;
	quanBit = "";
	if (index > MaxLevel)
		index = MaxLevel;
	else if (index < -MaxLevel)
		index = -MaxLevel;
	if (index ==0)
	{
		quanBit += "00";
		return 2;
	}
	else if (index == 1)
	{
		quanBit += "01";
		return 2;
	}
	else if (index == -1)
	{
		quanBit += "10";
		return 2;
	}
	else
	{
		int footindex = int(abs(2*index)+((index>0)?(-0.5):0));
		for (int i=1; i< footindex; i++)
			quanBit+="1";
		quanBit+="0";
		return footindex;
	}
}
template<class T, size_t cn>
int QNode<T,cn>::UnaryCoding(int index)
{
	const int MaxLevel = 50;
	if (index > MaxLevel)
		index = MaxLevel;
	else if (index < -MaxLevel)
		index = -MaxLevel;
	//sign bit
	if (index >= 0)
		quanBit = "0";
	else 
	{
		quanBit = "1";
		index = -index;
	}
	while (index >1 )
	{
		quanBit+="1";
		index --;
	}
	quanBit+="0"; // end bit

	return int(quanBit.length());
	

}
//////////////////////////////////// Quilting ////////////////////


template<class T, size_t cn>
vector<vector<int>>& QNode<T,cn>::Quilting(QNode<T,cn>& ts) //quilting must be rechecked / rewrite to debug!! 03032013
{
	QNode<T,cn>::Matrix3Di path;
	path = FindPath(ts);
	seam = Blending(path,ts);
	//leftBound.Print();
	//upBound.Print();
	this->SetBlock(ts);//use predicted to substitude qNode
	return seam;
}
template<class T, size_t cn>
vector<vector<int>>& QNode<T,cn>::Quilting(vector<Tensor<T,cn>>& boundA, vector<Tensor<T,cn>>& boundB,BoundDir dir)
{
	QNode<T,cn>::Matrix3Di path;
	//Tensor<T,cn>::DisplayAll(boundA,1,1);
  //Tensor<T,cn>::DisplayAll(boundB,1,1);
	path = FindPath(boundA,boundB,dir);
	seam = Blending(path,boundA,boundB,dir);
  //Tensor<T,cn>::DisplayAll(boundA,1,1);
  //Tensor<T,cn>::DisplayAll(boundB,1,1);
	//SetBlock(ts);
	return seam;
}
template<class T, size_t cn>
vector<vector<int>> QNode<T,cn>::Blending(Matrix3Di path, vector<Tensor<T,cn>>& boundA, vector<Tensor<T,cn>>& boundB,BoundDir dir)
{
	cv::Mat tempMat, seamRst,crossRst,tempRst, tempOrg, temp;
	cv::Mat kernel;
	vector<cv::Mat> maskRGB;
	vector<vector<int>> seam(path.size());
	cv::Mat mask;
	//cv::Point3i pos;
	Vec3b tempElm(255,0,0);
	int position;
	Matrix3Di::iterator iter = path.begin();
	//pos = offset();
	//int index =0;
	int filterLength = min(boundA[0].size().height,boundA[0].size().width);
	filterLength = filterLength/2+1;
	for (unsigned int s=0; s<path.size(); s++)
	{
		tempMat = mylib::VecToMat(path[s],DataType<T>::depth);
		kernel = mylib::BinomialKernel(filterLength,DataType<T>::depth);
		//tempMat = mylib::Filter(tempMat,kernel,BORDER_REFLECT);//remove this 080413
		
    if (cn==3)
		{
			maskRGB.clear();
			maskRGB.push_back(tempMat);
			maskRGB.push_back(tempMat);
			maskRGB.push_back(tempMat);
			cv::merge(maskRGB,mask);
		}
		else
			mask = tempMat;
		//mylib::DisplayMat(mask,"1st");
		cv::multiply(boundB[s][0],mask,temp); 
		cv::subtract(1,tempMat,tempMat);
		if (cn==3)
		{
			maskRGB.clear();
			maskRGB.push_back(tempMat);
			maskRGB.push_back(tempMat);
			maskRGB.push_back(tempMat);
			cv::merge(maskRGB,mask);
		}
		else
			mask = tempMat;
		//mylib::DisplayMat(mask,"2nd");
		cv::multiply(boundA[s][0],mask,boundA[s][0]); 
		cv::add(boundA[s][0],temp,boundA[s][0]);
		vector<vector<int>>& tempBoundary=*iter;
		if (s==1 || boundA[s].offset().x == 0 || dir == BoundDir::LEFT)
		{
			for (unsigned int i=0; i< iter->size(); i++)
			{
				position = int(iter->begin()->size()-1);
				for (unsigned int j=0; j< iter->begin()->size()-1; j++)
				{
					if (tempBoundary[i][j]==1 && tempBoundary[i][j+1]==0)
						{
							position = j;
							break;
						}
				}
				seam[s].push_back(position);
		//			seamRst.at<Vec3b>(i,position)=tempElm;
		//			//Vec3b & mmm=tempRst.at<cv::Vec3b>(i,position);
		//			//STBlock(seamRst).Display(0);
			}
			
				
		}
		else
		{
			for (unsigned int j=0; j< iter->begin()->size(); j++)
			{
				position = int(iter->size()-1);
				for (unsigned int i=0; i< iter->size()-1; i++)
				{
					if (tempBoundary[i][j]==1 && tempBoundary[i+1][j]==0)
					{
						position = i;
						break;
					}
				}
				seam[s].push_back(position);
				
		//	seamRst.at<Vec3b>(position,j)=tempElm;
		//			//STBlock(seamRst).Display(0);

		//			//Vec3b & mmm=tempRst.at<cv::Vec3b>(position,j);
			}			
		}
		iter++;
		//	//STBlock(tempRst).Display(0);
		//	//STBlock(seamRst).Display(0);
		//	//copy to the seam version
		//	this->rst.SetBlock(STBlock(tempRst),pos - startPos[index]);	
		//	//this->rst.Display(0);
		//	this->rst_with_seam.SetBlock(STBlock(seamRst),pos-startPos[index]);			
		//	//rst_with_seam.Display(0);
		//	//this->rst.Display(0);
		//	//end lable the seam/////////////////////////////
		//	orgIter++;
		//	rstIter++;
		//	index++;
		}
		return seam;
		////copy the rest part 

}
template<class T, size_t cn>
vector<vector<int>> QNode<T,cn>::Blending(Matrix3Di path,QNode<T,cn>& ts, int criteria)
{
	cv::Mat tempMat, seamRst,crossRst,tempRst, tempOrg, temp;
	cv::Mat kernel;
	vector<cv::Mat> maskRGB;
	vector<vector<int>> seam(path.size());
	cv::Mat mask;
	cv::Point3i pos;
	Vec3b tempElm(255,0,0);
	int position;
	Matrix3Di::iterator iter = path.begin();
	//vector<cv::Point3i> startPos = GetBoundOffset(refNode,1);
	if (criteria ==0)
	{
		pos = this->offset();
		//int index =0;
		int filterLength = min(overlapSize.height,overlapSize.width);//big block will blur too much! 
		filterLength = filterLength/2+1;
		for (unsigned int s=0; s<path.size(); s++)
		{
			tempMat = mylib::VecToMat(path[s],DataType<T>::depth);
			kernel = mylib::BinomialKernel(filterLength,DataType<T>::depth);
			
			//	mylib::DisplayMat(kernel);
			tempMat = mylib::Filter(tempMat,kernel,BORDER_REFLECT);
			if (cn==3)
			{
				maskRGB.clear();
				maskRGB.push_back(tempMat);
				maskRGB.push_back(tempMat);
				maskRGB.push_back(tempMat);
				cv::merge(maskRGB,mask);
			}
			else
				mask = tempMat;
			
			cv::multiply(bounds[s]->GetFrame(0),mask,temp); // keep the left and up as rst
			cv::subtract(1,tempMat,tempMat);
			if (cn==3)
			{
				maskRGB.clear();
				maskRGB.push_back(tempMat);
				maskRGB.push_back(tempMat);
				maskRGB.push_back(tempMat);
				cv::merge(maskRGB,mask);
			}
			else
				mask = tempMat;
			cv::multiply(ts.bounds[s]->GetFrame(0),mask,bounds[s]->GetFrame(0)); // change right and below as the predicted
			cv::add(bounds[s]->GetFrame(0),temp,bounds[s]->GetFrame(0));
			vector<vector<int>>& tempBoundary=*iter;
			if (s==1 || this->offset().x == 0)
			{
				for (unsigned int i=0; i< iter->size(); i++)
				{
					position = int(iter->begin()->size()-1);
					for (unsigned int j=0; j< iter->begin()->size()-1; j++)
					{
						if (tempBoundary[i][j]==1 && tempBoundary[i][j+1]==0)
						{
							position = j;
							break;
						}
					}
				seam[s].push_back(position);
		//			seamRst.at<Vec3b>(i,position)=tempElm;
		//			//Vec3b & mmm=tempRst.at<cv::Vec3b>(i,position);
		//			//STBlock(seamRst).Display(0);
				}
			
				
			}
			else
			{
				for (unsigned int j=0; j< iter->begin()->size(); j++)
				{
					position = (int)iter->size()-1;
					for (unsigned int i=0; i< iter->size()-1; i++)
					{
						if (tempBoundary[i][j]==1 && tempBoundary[i+1][j]==0)
						{
							position = i;
							break;
						}
					}
				seam[s].push_back(position);
				
		//	seamRst.at<Vec3b>(position,j)=tempElm;
		//			//STBlock(seamRst).Display(0);

		//			//Vec3b & mmm=tempRst.at<cv::Vec3b>(position,j);
				}			
			}
			iter++;
		//	//STBlock(tempRst).Display(0);
		//	//STBlock(seamRst).Display(0);
		//	//copy to the seam version
		//	this->rst.SetBlock(STBlock(tempRst),pos - startPos[index]);	
		//	//this->rst.Display(0);
		//	this->rst_with_seam.SetBlock(STBlock(seamRst),pos-startPos[index]);			
		//	//rst_with_seam.Display(0);
		//	//this->rst.Display(0);
		//	//end lable the seam/////////////////////////////
		//	orgIter++;
		//	rstIter++;
		//	index++;
		}
		return seam;
		////copy the rest part 


	}
	else
	{
		//undef
		exit(-1);
		return seam;
	}
}



template<class T, size_t cn>
typename QNode<T,cn>::Matrix3Di QNode<T,cn>::FindPath(const vector<Tensor<T,cn>>& boundA,vector<Tensor<T,cn>>& boundB,BoundDir dir)
{
	if (boundA.size()!= boundB.size())
		CV_Error(CV_StsBadSize,"the two input boundary vectors do not have same number of elements");
	unsigned int bSize =(unsigned int) boundA.size();

	for (unsigned int i=0; i < boundA.size(); i++)
	{
		if (boundA[i].size()!= boundB[i].size())
			CV_Error(CV_StsBadSize,"difference boundary block size");
	}
	QNode<T,cn>::Matrix3Di path;
	cv::Mat dpmx,emx;
	cv::Mat temp; 
	T tempSum=0;
	Vec<T,cn> elem;
	//int count=0;
	for (unsigned int s=0; s< bSize; s++)
	{
		Size3 sz = boundA[s].size();
		emx = cv::Mat::zeros(sz.height, sz.width,DataType<T>::depth);
		dpmx = -1*cv::Mat::ones(sz.height, sz.width,DataType<T>::depth);
		path.push_back(vector<vector<int>>(sz.height,vector<int>(sz.width)));
		for ( int i = 0; i< sz.depth; i++)
		{
			cv::pow(cv::abs(boundA[s][i]-boundB[s][i]),/*2*/3,temp);//try different pow, gj01072013
			for (int i = 0; i<emx.size().height; i++)
				for (int j =0 ;j<emx.size().width;j++)
				{
					// I verified in this case, -----------------> x
					//							|
					//							|
					//							V
					//							y
					elem = temp.at<Vec<T,cn>>(i,j); 
					tempSum=0;
					for (int k=0; k<cn; k++)
					{
						tempSum+= elem[k];
					}
					emx.at<T>(i,j) = (T)std::sqrt(double(tempSum/cn));
				}
			if ( s==1 || boundA[s].offset().x==0 || dir==BoundDir::LEFT) //vertical
			{	
				emx.row(0).copyTo(dpmx.row(0));
				for (int l = 0; l<boundA[s].size().width; l++)
					ShortestPathDP(emx.size().height-1,l,emx,dpmx,*(path.end()-1),Directions::DIRECTION_VERTICAL);
				MakePath(*(path.end()-1),dpmx,Directions::DIRECTION_VERTICAL);
			}
			else // horizontal
			{	emx.col(0).copyTo(dpmx.col(0));
				for (int l=0; l<boundA[s].size().height; l++)
					ShortestPathDP(l,emx.size().width-1,emx,dpmx,*(path.end()-1),Directions::DIRECTION_HORIZONTAL);
				MakePath(*(path.end()-1),dpmx,Directions::DIRECTION_HORIZONTAL);
			}

		}	
			
	}
	if ( bSize  == 2)
		TreatCrossRegion(path,CornerPos::CORNER_SE);
	return path;

}







template<class T, size_t cn>
typename QNode<T,cn>::Matrix3Di QNode<T,cn>::FindPath(const QNode<T,cn>& ts) 
{
	//only work for 2D image
	QNode<T,cn>::Matrix3Di path;
	cv::Mat dpmx,emx;
	cv::Mat temp; 
	T tempSum=0;
	Vec<T,cn> elem;
	if (leftBound.size()!= ts.leftBound.size() || upBound.size()!=ts.upBound.size())
	{
		CV_Error(CV_StsBadSize,"The boundaries size are different!");
	}
	//int count=0;
	//for safty reconstruct bounds 
	bounds.clear();
	bounds.push_back(&upBound);
	bounds.push_back(&leftBound);
	for (unsigned int s=0; s< bounds.size(); s++)
	{
		Size3 sz = bounds[s]->size();
		emx = cv::Mat::zeros(sz.height, sz.width,DataType<T>::depth);
		dpmx = -1*cv::Mat::ones(sz.height, sz.width,DataType<T>::depth);
		path.push_back(vector<vector<int>>(sz.height,vector<int>(sz.width)));
		for ( int i = 0; i< sz.depth; i++)
		{
			cv::pow(cv::abs(bounds[s]->GetFrame(i)-ts.bounds[s]->GetFrame(i)),2/*3*/,temp);//try different pow, gj01072013

			for (int i = 0; i<emx.size().height; i++)
				for (int j =0 ;j<emx.size().width;j++)
				{
					// I verified in this case, -----------------> x
					//							|
					//							|
					//							V
					//							y
					elem = temp.at<Vec<T,cn>>(i,j); 
					tempSum=0;
					for (int k=0; k<cn; k++)
					{
						tempSum+= elem[k];
					}
					emx.at<T>(i,j) = (T)std::sqrt(double(tempSum/cn));
				}
			if ( s==1 )//changed gj03032013 , || bounds[s]->offset().x==0 )
			{	
				emx.row(0).copyTo(dpmx.row(0));
				for (int l = 0; l<overlapSize.width; l++)
					ShortestPathDP(emx.size().height-1,l,emx,dpmx,*(path.end()-1),Directions::DIRECTION_VERTICAL);
				MakePath(*(path.end()-1),dpmx,Directions::DIRECTION_VERTICAL);
			}
			else
			{	emx.col(0).copyTo(dpmx.col(0));
				for (int l=0; l<overlapSize.height; l++)
					ShortestPathDP(l,emx.size().width-1,emx,dpmx,*(path.end()-1),Directions::DIRECTION_HORIZONTAL);
				MakePath(*(path.end()-1),dpmx,Directions::DIRECTION_HORIZONTAL);
			}

		}	
			
	}
	if ( bounds.size()  == 2)
		TreatCrossRegion(path);
	return path;
}

template<class T, size_t cn>
double QNode<T,cn>::ShortestPathDP(int i, int j, const cv::Mat &emx, cv::Mat &dpmx, vector<vector<int>> &pathPlane, Directions direction)
{
	if (i<0 || j<0)
	{
		//std::cout<<"error i or j \n";
		return FLT_MAX;
	}
	else if (i>= emx.size().height || j >= emx.size().width)
	{
		//std::cout<<"out of boundary \n";
		return FLT_MAX;
	}
	else if (dpmx.at<T>(i,j)>=0)
	{
		return dpmx.at<T>(i,j);
	}
	else
	{
		if (direction == Directions::DIRECTION_VERTICAL)
		{
			T a=(T)ShortestPathDP(i-1,j-1,emx,dpmx,pathPlane,Directions::DIRECTION_VERTICAL);
			T b=(T)ShortestPathDP(i-1,j,emx,dpmx,pathPlane, Directions::DIRECTION_VERTICAL);
			T c=(T)ShortestPathDP(i-1,j+1,emx,dpmx,pathPlane, Directions::DIRECTION_VERTICAL);
			int tempPath=-1;
			dpmx.at<T>(i,j)= emx.at<T>(i,j);
			if (a<b)
			{
				if (a<c)
				{
					dpmx.at<T>(i,j)+=a;
					tempPath=j-1;
				}
				else
				{
				
					dpmx.at<T>(i,j)+=c;
					tempPath=j+1;
				}
			}
			else
			{
				if (b<c)
				{
					dpmx.at<T>(i,j)+=b;
					tempPath=j;
				}
				else
				{
					dpmx.at<T>(i,j)+=c;
					tempPath=j+1;
				}
			}
			pathPlane[i][j]=tempPath;
		}
		else
		{
			T a=(T)ShortestPathDP(i-1,j-1,emx,dpmx,pathPlane,Directions::DIRECTION_HORIZONTAL);
			T b=(T)ShortestPathDP(i,j-1,emx,dpmx,pathPlane,Directions::DIRECTION_HORIZONTAL);
			T c=(T)ShortestPathDP(i+1,j-1,emx,dpmx,pathPlane,Directions::DIRECTION_HORIZONTAL);
			int tempPath=-1;
			dpmx.at<T>(i,j)= emx.at<T>(i,j);
			if (a<b)
			{
				if (a<c)
				{
					dpmx.at<T>(i,j)+=a;
					tempPath=i-1;
				}
				else
				{
				
					dpmx.at<T>(i,j)+=c;
					tempPath=i+1;
				}
			}
			else
			{
				if (b<c)
				{
					dpmx.at<T>(i,j)+=b;
					tempPath=i;
				}
				else
				{
					dpmx.at<T>(i,j)+=c;
					tempPath=i+1;
				}
			}
			pathPlane[i][j]=tempPath;
		}

		//consider memerize the path at this point!!!!!
		/////////////////////////////////////////

		//std::cout<<"dpmx("<<i<<","<<j<<")="<<dpmx.at<float>(i,j)<<std::endl;
		//std::cout<<"emx("<<i<<","<<j<<")="<<emx.at<float>(i,j)<<std::endl;	
		return dpmx.at<T>(i,j);
	}
}

template<class T, size_t cn>
void QNode<T,cn>::TreatCrossRegion(Matrix3Di & path, CornerPos corner)
{
	QNode<T,cn>::Matrix3Di backup = path;
	unsigned int m,n,M,N;
	m = (unsigned int)path[0].size();
	n = (unsigned int)path[1][0].size();
	M = (unsigned int)path[0][0].size();
	N = (unsigned int)path[1].size(); //actually it is the overlap size
		for (unsigned int i=0; i< path[0].size();i++)
			for (unsigned int j=0; j<path[1][0].size();j++)
			{
				if (corner == CornerPos::CORNER_NW)
				{
					if (backup[0][i][j]==1 && backup[1][i][j] ==0)
						path[1][i][j]=1;
					if (backup[0][i][j]==0 && backup[1][i][j] ==1)
						path[0][i][j]=1;
				}
				else if (corner ==CornerPos::CORNER_SE)
				{

					if (backup[0][m-i-1][N-j-1]==1 && backup[1][M-i-1][n-j-1] ==0)
						path[1][M-i-1][n-j-1]=1;
					if (backup[0][m-i-1][N-j-1]==0 && backup[1][M-i-1][n-j-1] ==1)
						path[0][m-i-1][N-j-1]=1;
				}
				else
				{
					CV_Error(CV_StsBadFunc,"unsupport");
				}
			}
	
}

template<class T, size_t cn>
void QNode<T,cn>::MakePath(vector<vector<int>>&pathPlane, cv::Mat &dpmx, Directions direction)
{

	//mylib::DisplayMat(dpmx);
	T temp=saturate_cast<T>(FLT_MAX);
	int startPoint=0;
	int tempPoint=0;
	int height = (int) pathPlane.size();
	int width  = (int) pathPlane[0].size();
	
	if (direction == Directions::DIRECTION_VERTICAL)
	{
		for (int i=0;i<width;i++)
		{
			if (temp > dpmx.at<T>(height-1,i))
			{
				temp = dpmx.at<T>(height-1,i);
				tempPoint = i;//pathPlane[height-1][i];
			}
		}
		for (int i=height; i>0;i--)
		{
			startPoint = tempPoint;
			tempPoint= pathPlane[i-1][startPoint];
			for (int j=0;j<width;j++)
			{
				if (j<=startPoint)
					pathPlane[i-1][j]=1;
				else
					pathPlane[i-1][j]=0;
			}
		//	pathPlane[i-1][startPoint]=1;
		}
		
	}
	else
	{
		for (int i=0;i<height;i++)
		{
			if (temp > dpmx.at<T>(i,width-1))
			{
				temp = dpmx.at<T>(i,width-1);
				tempPoint =  i;//pathPlane[i][width-1];
			}
		}
		for (int i=width; i>0;i--)
		{
			startPoint = tempPoint;
			tempPoint= pathPlane[startPoint][i-1];
			for (int j=0;j<height;j++)
			{
				if (j<=startPoint)
					pathPlane[j][i-1]=1;
				else
					pathPlane[j][i-1]=0;
			}
			//	pathPlane[startPoint][i-1]=1;
		}
	}

}

//////////////////////////////////// Boundary operations /////////////////
template<class T, size_t cn> 
Tensor<T,cn> QNode<T,cn>::GetExtendTensor(bool up, bool left , bool down, bool right) const
{
	int sz_v=0,sz_h=0;
	if (up)
		sz_v++;
	if (left)
		sz_h++;
	if (down)
		sz_v++;
	if (right)
		sz_h++;

	Size3 s1 = this->size();
	Size3 s2(this->overlap().height*sz_v,this->overlap().width*sz_h,0);
	Tensor<T,cn> extNode(s1+s2);

	cv::Point3i center(this->overlap().height,this->overlap().width,0);
	cv::Point3i l(0,0,0),d(this->size().height+this->overlap().height,0,0),u(0,0,0),r(0,this->size().width+this->overlap().width,0);
	extNode.SetBlock(center,*this);
	if (up)
	{
		extNode.SetBlock(u,this->upBound);
	}
	if (left)
	{
		extNode.SetBlock(l, this->leftBound);
	}
	if (down)
	{
		extNode.SetBlock(d, this->lowBound);
	}
	if (right)
	{
		extNode.SetBlock(r, this->rightBound);
	}
	//extNode.SetSubWinSize(this->GetSubWinSize());
	//extNode.SetSubWinStep(this->GetSubWinStep());
	return extNode;
}
template<class T, size_t cn> 
Tensor<T,cn> QNode<T,cn>::GetExtendTensor( Vec<T,cn> padding, bool up, bool left , bool down, bool right) const
{
	int sz_v=0,sz_h=0;
	if (up)
		sz_v++;
	if (left)
		sz_h++;
	if (down)
		sz_v++;
	if (right)
		sz_h++;

	Size3 s1 = this->size();
	Size3 s2(this->overlap().height*sz_v,this->overlap().width*sz_h,0);
	Tensor<T,cn> extNode(s1+s2);

	cv::Point3i center(this->overlap().height,this->overlap().width,0);
	cv::Point3i l(0,0,0),d(this->size().height+this->overlap().height,0,0),u(0,0,0),r(0,this->size().width+this->overlap().width,0);
	extNode.SetBlock(center,*this);
	if (up)
	{
		extNode.SetBlock(u,Tensor<T,cn>(this->upBound.size(),padding));
	}
	if (left)
	{
		extNode.SetBlock(l, Tensor<T,cn>(this->leftBound.size(),padding));
	}
	if (down)
	{
		extNode.SetBlock(d, Tensor<T,cn>(this->lowBound.size(),padding));
	}
	if (right)
	{
		extNode.SetBlock(r, Tensor<T,cn>(this->rightBound.size(),padding));
	}
	//extNode.SetSubWinSize(this->GetSubWinSize());
	//extNode.SetSubWinStep(this->GetSubWinStep());
	return extNode;
}
template<class T, size_t cn>
Tensor<T,cn>& QNode<T,cn>::GetBoundaryLeft(void)
{
	return this->leftBound;

}


template<class T, size_t cn>
const Tensor<T,cn>& QNode<T,cn>::GetBoundaryLeft(void) const
{
	return this->leftBound;

}

template<class T, size_t cn>
Vec<T,cn>& QNode<T,cn>::GetBoundaryLeft(const Point3i& pos)
{

	return this->leftBound[pos];

}

template<class T, size_t cn>
Tensor<T,cn>& QNode<T,cn>::GetBoundaryUp(void)
{
	return this->upBound;

}

template<class T, size_t cn>
const Tensor<T,cn>& QNode<T,cn>::GetBoundaryUp(void) const
{
	return this->upBound;

}

template<class T, size_t cn>
Vec<T,cn>& QNode<T,cn>::GetBoundaryUp(const Point3i&pos)
{
	return this->upBound[pos];

}

template<class T, size_t cn>
Tensor<T,cn>& QNode<T,cn>::GetBoundary(int i)
{
	if (this->bounds.size() <= (unsigned int)i)
	{
		//oor
		exit(-1);
	}
	return *(this->bounds[i]);
}


template<class T, size_t cn>
const Tensor<T,cn>& QNode<T,cn>::GetBoundary(int i) const
{
	if (this->bounds.size() <= (unsigned int)i)
	{
		//oor
		exit(-1);
	}
	return *(this->bounds[i]);
}

template<class T, size_t cn>
int QNode<T,cn>::GetBoundarySize(void) const
{
	return int(this->bounds.size());
}
template<class T, size_t cn>
void QNode<T,cn>::Ref(const Tensor<T,cn>& ts, const Cube& roi, const Size3& boundarySize)
{
	Size3 bsize = boundarySize;
	if (boundarySize.height == 0)
		bsize.height = 1;
	if (boundarySize.width == 0)
		bsize.width = 1;
	Size3 sz = roi.size();
	Point3i pos = roi.offset();
	ts.Ref(roi,*this);
	//*this = QTree(ts.Crop(pos,sz));
	bits=0;
	predicted_method=CodingMethodNames::NO_CODING;
	qsize = 0;
	footPos = this->tsSize.Point3() - Point3i(1,1,1);//Point3i(size().height-1,size().width-1,size().depth-1);
	nodeFoot = ts[footPos+roi.offset()];
	overlapSize = boundarySize;
    Size3 extSize = sz;
  Point3i extPos=pos;
	//leftBound = Tensor<T,cn>(Size3(bsize.height+sz.height,bsize.width,sz.depth));
	//upBound = Tensor<T,cn>(Size3(bsize.height,bsize.width+sz.width,sz.depth));
	//lowBound = Tensor<T,cn>(Size3(bsize.height,bsize.width+sz.width,sz.depth));
	//rightBound = Tensor<T,cn>(Size3(bsize.height+sz.height,bsize.width,sz.depth));
	if (pos.x !=0 && pos.y==0)
	{
		ts.Ref(Cube(Point3i(pos.x-bsize.height, pos.y, pos.z - bsize.depth),Size3(bsize.height,sz.width,sz.depth)),upBound);
		//leftBound = Tensor<T,cn>(
		//boundary.push_back(ts.Crop( Point3i(pos.x-bsize.height, pos.y, pos.z - bsize.depth),
		//						Size3(bsize.height,sz.width,sz.depth)));
		//bounds.push_back(&upBound);
        extPos.x = pos.x-bsize.height;
    extPos.y = 0;
    extSize.height+=bsize.height;
	}
	else if (pos.y !=0 && pos.x==0)
	{
		ts.Ref(Cube(Point3i(pos.x, pos.y-bsize.width,pos.z - bsize.depth), Size3(sz.height,bsize.width,sz.depth)),leftBound);
		//boundary.push_back(ts.Crop( Point3i(pos.x, pos.y-bsize.width,pos.z - bsize.depth),
		//						Size3(sz.height,bsize.width,sz.depth)));
		//bounds.push_back(&leftBound);
        extPos.x=0;
    extPos.y = pos.y-bsize.width;
    extSize.width+=bsize.width;
	}
	else if (pos.x!=0 && pos.y!=0)
	{
		ts.Ref(Cube(Point3i(pos.x-bsize.height, pos.y-bsize.width,pos.z - bsize.depth), Size3(bsize.height,sz.width+bsize.width,sz.depth)),upBound);
		ts.Ref(Cube(Point3i(pos.x-bsize.height, pos.y-bsize.width,pos.z - bsize.depth),	Size3(sz.height+bsize.height,bsize.width,sz.depth)),leftBound);	
		//bounds.push_back(&upBound);
		//bounds.push_back(&leftBound);
    	  extPos.x=pos.x - bsize.height;
    extPos.y=pos.y - bsize.width;
    extSize.height+=bsize.height;
    extSize.width+=bsize.width;

	}
	if (pos.x + sz.height+ bsize.height -1 < ts.size().height) // have lower bound
	{
		if (pos.y-bsize.width < 0)
			ts.Ref(Cube(pos.x+sz.height,0, pos.z - bsize.depth,bsize.height, sz.width+bsize.width-pos.y,sz.depth),lowBound);
		else if (pos.y+sz.width+bsize.width > ts.size().width)
			ts.Ref(Cube(pos.x+sz.height,pos.y-bsize.width,0,bsize.height,sz.width+bsize.width,sz.depth),lowBound);
		else
			ts.Ref(Cube(pos.x+sz.height,pos.y-bsize.width,0,bsize.height,sz.width+2*bsize.width,sz.depth),lowBound);
        extSize.height+=bsize.height;
	}
	if (pos.y + sz.width + bsize.width -1 < ts.size().width)
	{
		if (pos.x - bsize.height < 0)
			ts.Ref(Cube(0, pos.y + sz.width, pos.z - bsize.depth, sz.height+bsize.height, bsize.width, sz.depth),rightBound);
		else if (pos.x+sz.height+bsize.height > ts.size().height)
			ts.Ref(Cube(pos.x - bsize.height, pos.y + sz.width, pos.z - bsize.depth, sz.height+bsize.height, bsize.width, sz.depth),rightBound);
		else
			ts.Ref(Cube(pos.x - bsize.height, pos.y + sz.width, pos.z - bsize.depth, sz.height+2*bsize.height, bsize.width, sz.depth),rightBound);
      extSize.width+=bsize.width;
	}
	bounds.push_back(&upBound);
	bounds.push_back(&leftBound);
  ts.Ref(Cube(extPos,extSize),eNode);

}

template <class T, size_t cn>
void QNode<T,cn>::SetLeftFoot(const typename QNode<T,cn>::value_type & left)
{
  this->leftFoot = left;
}

template <class T, size_t cn>
void QNode<T,cn>::SetUpFoot(const typename QNode<T,cn>::value_type & up)
{
  this->upFoot = up;
}


template <class T, size_t cn>
void QNode<T,cn>::SetFoot(const Vec<T,cn>& ft)
{
	this->nodeFoot = ft;
}

template <class T, size_t cn>
void QNode<T,cn>::AddFoot(const pair<Point3i,FootItem>& f)
{
	this->feet.push_back(f);
}


template <class T, size_t cn>
void QNode<T,cn>::clearFoot(void)
{
	this->feet.clear();
}
template <class T, size_t cn>
vector<pair<Point3i,FootItem>> QNode<T,cn>::getFeetCopy(void) const
{
  return this->feet;
}
template <class T, size_t cn>
int QNode<T,cn>::getFeetNumber(void) const
{
  return this->feet.size();
}

template <class T, size_t cn>
void QNode<T,cn>::SetFoot(const Tensor<T,cn>& src)
{
	this->nodeFoot = src[this->tsOffset+this->footPos];
}
template <class T, size_t cn>
Vec<T,cn> QNode<T,cn>::GetFoot() const
{
	return this->nodeFoot;//return a copy
}
template <class T, size_t cn>
pair<Point3i,FootItem> QNode<T,cn>::GetFoot(int i) const
{
	return this->feet[i];//return a copy
}
template <class T, size_t cn>
void QNode<T,cn>::SetFoot(const pair<Point3i,FootItem>& f,int i) 
{
  this->feet[i]=f;
}
template <class T, size_t cn>
Vec<T,cn> QNode<T,cn>::GetApproxFoot() const
{
  return this->approxFoot;
}
template <class T, size_t cn>
QNode<T,cn>& QNode<T,cn>::PoissonLightingCorrection(const QNode<T,cn>& changeTo,const Tensor<T,cn>& ref, const Tensor<T,cn>& rec, Size3 boundSize, int footRegion, int footMethod, double qsize)
{
  //ref is original image 
  //rec is reconstructed image
	//boundSize is the size of bound which is also lighting corrected
  this->qsize = qsize;
  vector<FootType> feetToPb; //the feet in PLC is different, it just contain x,y,value, no coding information

	//this->Print();
	//changeTo.Print();
  Tensor<T,cn> candExt = rec.GetBlock(Cube(this->offset()-this->overlap().Point3()-Point3i(1,1,0),this->size()+Size3(1,1,0)*2 + this->overlap().Point3()*2)).ExtendBoundary(boundSize-Size3(1,1,0),0);
 /* candExt.Print();*/
	//Tensor<T,cn> candExt = this->GetExtendTensor(1,1,1,1).ExtendBoundary(boundSize,0);
	Tensor<T,cn> orgExt = ref.GetBlock(Cube(changeTo.offset()-changeTo.overlap().Point3()*2,candExt.size()));//there was a bug here, it use tar instead of org to compute lighting feet
                                                                                                           //changeTo.GetExtendTensor(1,1,1,1);
  //candExt.Display();
  Tensor<T,cn> tarExt = rec.GetBlock(Cube(rec.offset()+changeTo.offset()-changeTo.overlap().Point3()*2, orgExt.size()));
  //tarExt.Display();
  //tarExt.Print();
	CV_DbgAssert(boundSize.height<=this->overlap().height && boundSize.width<=this->overlap().width);
	//T ave = changeTo.GetFoot()[0];
  //mask it self is the size of block + 2*boundSize (twice the blk size)
  //but the non-zero value is limited to block+1*boundSize
	Tensor<T,1> mask(this->size()+boundSize /* *2 */,255);//gj 2013-05-15 since only do forward (L,U) border, change boundSize to 1x
  Tensor<T,1> mask2 = mask.ExtendHalfBoundary(/* 20130605 this->overlap()-*/boundSize,0,1);
  //mask2.Display();
	//mask = mask.ExtendBoundary(boundSize,255);
	Tensor<T,1> maskExt = mask2.ExtendBoundary(boundSize,0);
	//maskExt.Print();
	//compute foot
	//Tensor<T,cn> temp ;
	//ref.Ref(Cube(changeTo.offset()+(changeTo.size()/2).Point3(),changeTo.size()/Size3(2,2,1)),temp);
	//ref.Ref(Cube(changeTo.offset()-(changeTo.size()/2).Point3(),changeTo.size()*Size3(2,2,1)),temp);
		//ref.Crop(changeTo.offset(),changeTo.size()); // I don't do TP on the boundary of image, so it is safe to do crop this way
					/*double tempmean=0;
					for (int ii=0;ii<temp.size().height; ii++)
						for (int jj=temp.size().height-ii;jj<temp.size().width;jj++)
							tempmean+=temp(ii,jj,0)[0];
					tempmean/=(temp.size().area()/2);*/
	//temp.Print();
	//remember both qNode and changeTo are bordersize doubled = to size/2
	int caseNo = footRegion;
	int method = footMethod;
	//case 1://mean
	//		rst = temp.Mean();
	//		break;
	//	case 2://lower right triagle mean
	//		rst = temp.LRMean();
	//		break;
	//	case 3: //percentage mean
	//		rst = temp.PercentageMean(0.4,0.6);
	//		break;
	//	case 4: //median
	//		rst = temp.PercentageMean(0.5,0.5);
	//		break;
	//	case 5: //upper left triagnle mean
  this->feet = changeTo.getFeetCopy();
  vector<pair<Point3i,FootItem>> extraFoot;
  //the feet position get from RetriveFeet does not include any boundary foot
  //we need up and left foot on the blending boundary + (if any) foot not on the boundary
  //I must consider more than 3 feet case!
 // cout<<"size of feet:"<<feet.size();
  //cout<<"org:"<<changeTo.offset().x<<","<<changeTo.offset().y<<endl;

    
  for (auto& f : this->feet)
  {
    //  cout<<"foot"<<f.first.x<<","<<f.first.y<<endl;
    if (f.second.value<-255)//no previous foot , use the foot on the boundary
    {
      if ((f.first.x == changeTo.GetFootPos().x && f.first.y == changeTo.GetFootPos().y) //except the one on the SE corner (true foot)
          ||(f.first.x != changeTo.GetFootPos().x -  this->size().height && f.first.y == changeTo.GetFootPos().y) //on the right edge
          ||(f.first.x == changeTo.GetFootPos().x && f.first.y != changeTo.GetFootPos().y - this->size().width)//on the bottom edge
        )
      {
        f.second.value = this->ComputePLCFoot(ref,rec.offset()+Point3i(f.first.x+1,f.first.y+1,0)-changeTo.overlap().Point3(),changeTo.overlap()*2+Size3(0,0,1),caseNo,method)[0]; 
        f.first.x++;
        f.first.y++;
      }
      else if (f.first.y == changeTo.GetFootPos().y) //up foot or others on the right bd
      {
        f.second.value = this->ComputePLCFoot(rec,rec.offset()+Point3i(f.first.x+1,f.first.y+1,0)-Point3i(changeTo.overlap().height*3,changeTo.overlap().width*2,0),changeTo.overlap()*2+Size3(0,0,1),caseNo,method)[0]; 
        f.first.x -= changeTo.overlap().height;
        //cout<<"osize"<<changeTo.overlap().height;
        f.first.y++;
      }
      else if (f.first.x == changeTo.GetFootPos().x) //left foot or others on the bottom bd
      {
        f.second.value = this->ComputePLCFoot(rec,rec.offset()+Point3i(f.first.x+1,f.first.y+1,0)-Point3i(changeTo.overlap().height*2,changeTo.overlap().width*3,0),changeTo.overlap()*2+Size3(0,0,1),caseNo,method)[0]; 
        f.first.y -= changeTo.overlap().width;
        f.first.x++;
      }
      else //error there shouldn't be anything else
      {
        CV_Error(CV_StsNotImplemented,"not implment error");
      }

    }
    else //exist privors coded foot
    {
      //keep it but add more
      //cout<<"warning, found extra foot\n";
      if (f.first.y == changeTo.GetFootPos().y && f.first.x != changeTo.GetFootPos().x) //up foot or others on the right bd
      {
        cout<<"warning extra foot available\n";
        Point3i pos = rec.offset()+Point3i(f.first.x-changeTo.overlap().height,f.first.y+1,0);
        FootItem ff;
        ff.value = this->ComputePLCFoot(rec,rec.offset()+Point3i(f.first.x+1,f.first.y+1,0)-Point3i(changeTo.overlap().height*3,changeTo.overlap().width*2,0),changeTo.overlap()*2+Size3(0,0,1),caseNo,method)[0]; 
        f.first.y++;
        f.first.x++;
        extraFoot.push_back(f);
        f = pair<Point3i,FootItem>(pos,ff);

      }
      else if (f.first.y != changeTo.GetFootPos().y && f.first.x == changeTo.GetFootPos().x) //left foot or others on the bottom bd
      {
        cout<<"warning extra foot available\n";
        Point3i pos = rec.offset()+Point3i(f.first.x+1,f.first.y-changeTo.overlap().width,0);
        FootItem ff;
        ff.value = this->ComputePLCFoot(rec,rec.offset()+Point3i(f.first.x+1,f.first.y+1,0)-Point3i(changeTo.overlap().height*2,changeTo.overlap().width*3,0),changeTo.overlap()*2+Size3(0,0,1),caseNo,method)[0]; 
        f.first.x++;
        f.first.y++;
        extraFoot.push_back(f);
        f = pair<Point3i,FootItem>(pos,ff);
      }
      else
      {
        f.first.x++;
        f.first.y++;
      }
    }
    //  cout<<"foot"<<f.first.x<<","<<f.first.y<<endl;
      //if (f.x==this->GetFootPos().x && f.y==this->GetFootPos.y)
    //f.first.x ++;//= changeTo.overlap().height/2+1;
    //f.first.y ++;//= changeTo.overlap().width/2+1;

  }
 
  string _quanBitString="";
  if (this->feet.size()>0)
  {
  int level = int(log(double(this->feet.size()-1)/2)/log(2.0));
  int stride = (this->feet.size()-1)/2;
  //qantize foot
 

  this->nodeFoot = ((this->feet.end()-1)->second.value);
  if (changeTo.GetFoot(this->feet.size()-1).second.value<-255)
  {
    Vec<double,cn> pred_foot = (feet[0].second.value  + feet[(2<<level)].second.value)/2;
 	  this->footError = this->nodeFoot - pred_foot;
          this->approxFoot = Quantize(footError,qsize/2,SrcCodingMethod::UNARY_CODE) + pred_foot;
    (this->feet.end()-1)->second.value = approxFoot[0];
    (this->feet.end()-1)->second.bits = quanBit.size();
    (this->feet.end()-1)->second.bitstring = quanBit;
    _quanBitString+=quanBit;
 }
  //quan other foot
  ///horizontal and vertical
  for (int step = (1<<(level-1)); step >0; step/=2)
  {
    int count=0;
    for (int i=step; i < (1<<level); i+=step)
    {
      if (count%2==0)
      {
        int next_idx = i+step;
        next_idx==stride?next_idx=2*stride:next_idx;
        Vec<double,cn> pred_foot = (feet[i-step].second.value + feet[next_idx].second.value)/2;
        feet[i].second.value = (Quantize((feet[i].second.value-pred_foot[0]),qsize/2,SrcCodingMethod::UNARY_CODE) + pred_foot)[0];
        feet[i].second.bits = quanBit.size();
        feet[i].second.bitstring = quanBit;
        _quanBitString+=quanBit;
        pred_foot = (feet[stride+i-step].second.value + feet[stride+i+step].second.value)/2;
        feet[stride+i].second.value = (Quantize(feet[stride+i].second.value-pred_foot[0],qsize/2,SrcCodingMethod::UNARY_CODE)+pred_foot)[0];
        feet[stride+i].second.bits = quanBit.size();
        feet[stride+i].second.bitstring = quanBit;
        _quanBitString+=quanBit;
      }
      count++;
    }

  }
  
  /*
   //compute foot if not available
  if (changeTo.leftFoot[0]<-255) 
  {

	  leftFoot = this->ComputePLCFoot(rec,rec.offset()+changeTo.offset()-Point3i(0,changeTo.size().width,0),changeTo.size(),caseNo,method);//+changeTo.overlap());

  }
  else
    leftFoot = changeTo.leftFoot;
  if (changeTo.upFoot[0]<-255)
  {

	  upFoot = this->ComputePLCFoot(rec,rec.offset()+changeTo.offset()-Point3i(changeTo.size().height,0,0),changeTo.size(),caseNo,method);//+changeTo.overlap());
  }
  else
    upFoot = changeTo.upFoot;

  if (changeTo.GetFoot()[0]<-255)
  {
	  nodeFoot = this->ComputePLCFoot(ref,ref.offset()+changeTo.offset()-(changeTo.overlap()/2).Point3(),changeTo.size()+changeTo.overlap(),caseNo,method);//+changeTo.overlap());
    Vec<T,cn> pred_foot = (upFoot  + leftFoot)/2;
 	  footError = nodeFoot - pred_foot;
	  approxFoot = Quantize(footError,qsize/2,UNARY_CODE) + pred_foot;
    _quanBitString+=quanBit;
  }
  else
  {
    approxFoot = changeTo.GetFoot();
  }
  */
 /* upave = Quantize(upave,qsize,UNARY_CODE);
  _quanBitString+=quanBit;
  leftave=Quantize(leftave,qsize,UNARY_CODE);
  _quanBitString+=quanBit;*/
  //temp.PercentageMean(0.4,0.6);

	//ref.Ref(Cube(changeTo.offset()-Point3i(changeTo.size().height,0,0)+(changeTo.size()/2).Point3(),changeTo.size()/Size3(2,2,1)),temp);
	//Vec<T,cn> upave = temp.PercentageMean(0.4,0.6);//.LRMean();
	//ref.Ref(Cube(changeTo.offset()-Point3i(0,changeTo.size().width,0)+(changeTo.size()/2).Point3(),changeTo.size()/Size3(2,2,1)),temp);
	//Vec<T,cn> leftave = temp.PercentageMean(0.4,0.6);//.LRMean();
/*
					tar2.SetFoot(temp.Mean());	
					temp = ensemble.Crop(Point3i(qNode.offset().x-qNode.size().height/2, qNode.offset().y+qNode.size().width/2,0),qNode.size());
					tar2.upNodeFoot = temp.Mean();
					temp = ensemble.Crop(Point3i(qNode.offset().x+qNode.size().height/2, qNode.offset().y-qNode.size().width/2,0),qNode.size());
					tar2.leftNodeFoot =temp.Mean();*/
  //tarExt.Print();
  //candExt.Print();

  //adjust feet to help poissonsolver, since poissonsolver bases its coordinate system on blocksize + 2* overlapsize (quarter) + 2
  //qNode feet based on the referece (or constructed rec) image coordinate, so the offset of referece must be compensated first
  int count=0;
  for (auto& f : this->feet)
  {
    feetToPb.push_back(FootType(f));
  }
  
  //feetToPb = this->feet;
  for (auto ef : extraFoot)
  {
    feetToPb.push_back(FootType(ef));
  }
  for (auto& f : feetToPb)
  {
    //cout<<"foot"<<f.x<<","<<f.y<<endl;
    f.x -= changeTo.offset().x - changeTo.overlap().height-1;
    f.y -= changeTo.offset().y - changeTo.overlap().width-1;
    //if (count<stride)
    //  f.y++;//move up foot right 1 pixel to help pb
    //else 
    //  f.x++;//move left foot down 1 pixel to help pb
    //if (count==stride*2)
    //  f.y++;
    count++;
  }
  }//end if feet.size()>0 
  else//no feet 
  {
    for (auto& f : this->feet)
    {
      feetToPb.push_back(FootType(f));
    }
  }

	//PoissonSolver pb(candExt.GetFrame(0),tarExt.GetFrame(0),maskExt.GetFrame(0),orgExt.GetFrame(0), approxFoot[0],upFoot[0],leftFoot[0]);//changeTo.upNodeFoot[0],changeTo.leftNodeFoot[0]);
  PoissonSolver pb(candExt.GetFrame(0),tarExt.GetFrame(0),maskExt.GetFrame(0),orgExt.GetFrame(0), feetToPb);
  //PoissonSolver pb(candExt.GetFrame(0),tarExt.GetFrame(0),maskExt.GetFrame(0),ave[0],upave[0],leftave[0]);//changeTo.upNodeFoot[0],changeTo.leftNodeFoot[0]);
	cv::Mat dst;
	//mylib::DisplayMat(dst,"dst_before_LC");
	//candExt.Print();
	pb.poissonLightCorrection(dst,changeTo.offset().x,changeTo.offset().y);
  //Tensor<T,cn>(dst).Display();
	//pb.gradientStiching(128);
	//mylib::DisplayMat(dst,"dst_after_LC");
  //candExt.Print();
  auto feetBak = this->feet;
  if (this->feet.size()>0)
  {
  int _bits = this->bits;
  Vec<T,cn> _foot = this->approxFoot;
	*this = QNode<T,cn>(Tensor<T,cn>(dst),this->size(),this->overlap().Point3()*2,this->overlap());
  this->feet = feetBak;
	//this->SetBlock(dst(cv::Rect(this->overlap().width,this->overlap().height,this->size().width,this->size().height)));
	//this->GetBoundaryLeft().SetBlock(dst(cv::Rect(0,0,this->overlap().width,this->overlap().height+this->size().height)));
	//this->GetBoundaryUp().SetBlock(dst(cv::Rect(0,0,this->overlap().width+this->size().width,this->overlap().height)));
	//this->rightBound.SetBlock(dst(cv::Rect(this->overlap().width+this->size().width,0,this->overlap().width,this->overlap().height*2+this->size().height)));
	//this->lowBound.SetBlock(dst(cv::Rect(0,this->overlap().height+this->size().height,this->overlap().width*2+this->size().width,this->overlap().height)));
  this->quanBit = _quanBitString;
  this->bits = _bits;
  this->approxFoot = _foot;
  }
	return *this;
}

template <class T, size_t cn>
void QNode<T,cn>::GradientStitching(const QNode<T,cn>& changeTo, BlendingLocation blendPos, Size3 boundSize, const Tensor<T,1>& tempmask)
{
	Tensor<T,1>  maskExt;
	CV_DbgAssert(overlap().height>=boundSize.height && overlap().width>=boundSize.width && overlap().depth>=boundSize.depth);
	Tensor<T,1> mask(this->size()+boundSize*2);
	if (blendPos == BlendingLocation::CUSTOM_BLENDING)
	{
		mask = tempmask;
	}
	else if (blendPos == BlendingLocation::FORWARD_BLENDING)
	{
		Tensor<T,1> mask2(this->size(),0);
		mask2 = mask2.ExtendHalfBoundary(boundSize,255);
		mask.SetBlock(mask2);
		//mask.Print();
	}
	else if (blendPos == BlendingLocation::POST_BLENDING_RIGHT)
	{
		Tensor<T,1> mask2(this->size().height, boundSize.width,1, 255);
		mask.SetBlock(Point3i(0,boundSize.width+this->size().width,0),mask2);
	}
	else if (blendPos == BlendingLocation::POST_BLENDING_LOW)
	{
		Tensor<T,1> mask2(boundSize.width,this->size().width,1,255);
		mask.SetBlock(Point3i(boundSize.height+this->size().height,0,0),mask2);
	}
	maskExt = mask.ExtendBoundary(this->overlap()-boundSize,0);
	//maskExt.Print();
	//the order of tarExt and candExt are changed, since the blending is the target gradueatly change to (light corrected) candidate. 
	Tensor<T,cn> tarExt = this->GetExtendTensor(1,1,1,1);
	Tensor<T,cn> candExt = changeTo.GetExtendTensor(1,1,1,1);
	//tarExt.Print();
	//candExt.Print();

	PoissonSolver pb(candExt.GetFrame(0),tarExt.GetFrame(0),maskExt.GetFrame(0));
	cv::Mat dst;
	pb.gradientStiching(dst,this->offset().x, this-> offset().y);
	//mylib::DisplayMat(dst,"dst");
	QNode<T,cn> temp = QNode<T,cn>(Tensor<T,cn>(dst),this->size(),this->overlap().Point3(),this->overlap());
	
	if (blendPos == BlendingLocation::FORWARD_BLENDING)
	{
		this->leftBound.SetBlock(temp.GetBoundaryLeft());
		this->upBound.SetBlock(temp.GetBoundaryUp());
		this->SetBlock(temp);
	}
	else if (blendPos == BlendingLocation::POST_BLENDING_RIGHT)
	{
		this->rightBound.SetBlock(temp.rightBound);
		//this->lowBound.SetBlock(temp.lowBound);
	}
	else if (blendPos == BlendingLocation::POST_BLENDING_LOW)
	{
		this->lowBound.SetBlock(temp.lowBound);
	}
	else
	{
		CV_Error(CV_StsNotImplemented,"not supoort yet\n");
	}
}


template <class T, size_t cn>
Tensor<T,cn> QNode<T,cn>::LightingCorrection(const Tensor<T,cn>& changeTo, bool saveCodeLength)
{
	//this is different to the lighting correction in tensor because qNode also need to include 4 boundaries
	Tensor<T,cn>  candid= this->ExtendBoundary(Size3(overlapSize.height,overlapSize.width,overlapSize.depth));
	candid.SetBlock(upBound);
	candid.SetBlock(leftBound);
	candid.SetBlock(Point3i(0,this->tsSize.width+overlapSize.width,0),rightBound);
	candid.SetBlock(Point3i(this->tsSize.height+overlapSize.height,0,0),lowBound);
	//candid.Display(0);
  // all lighting information will be saved into this->light
  Tensor<T,cn> lighting = this->lt.LightingCorrection(candid,changeTo,saveCodeLength);
	//this->lt.SetLightingCodeLength(candid.lt.GetLightingCodeLength());
	lt.RecordLighting();
	//candid.Display(0);
	this->SetBlock(candid.Crop(overlapSize.Point3(),this->tsSize));
	this->upBound.SetBlock(candid.Crop(Point3i(0,0,0),upBound.size()));
	this->leftBound.SetBlock(candid.Crop(Point3i(0,0,0), leftBound.size()));
	this->lowBound.SetBlock(candid.Crop(Point3i(this->tsSize.height+overlapSize.height,0,0),lowBound.size()));
	this->rightBound.SetBlock(candid.Crop(Point3i(0,this->tsSize.width+overlapSize.width,0),rightBound.size()));
	return lighting;
}

template <class T, size_t cn>
Tensor<T,cn> QNode<T,cn>::LightingCorrection(const Tensor<T,cn>& changeTo, const Tensor<T,cn>& VQCodebook)
{
	//this is different to the lighting correction in tensor because qNode also need to include 4 boundaries
	Tensor<T,cn>  candid= this->ExtendBoundary(Size3(overlapSize.height,overlapSize.width,overlapSize.depth));
	candid.SetBlock(upBound);
	candid.SetBlock(leftBound);
	candid.SetBlock(Point3i(0,this->tsSize.width+overlapSize.width,0),rightBound);
	candid.SetBlock(Point3i(this->tsSize.height+overlapSize.height,0,0),lowBound);
	//candid.Display(0);
	Tensor<T,cn> lighting = lt.LightingCorrection(candid, changeTo,VQCodebook);
	//candid.Display(0);
	this->SetBlock(candid.Crop(overlapSize.Point3(),this->tsSize));
	this->upBound.SetBlock(candid.Crop(Point3i(0,0,0),upBound.size()));
	this->leftBound.SetBlock(candid.Crop(Point3i(0,0,0), leftBound.size()));
	this->lowBound.SetBlock(candid.Crop(Point3i(this->tsSize.height+overlapSize.height,0,0),lowBound.size()));
	this->rightBound.SetBlock(candid.Crop(Point3i(0,this->tsSize.width+overlapSize.width,0),rightBound.size()));
	return lighting;

}

template <class T, size_t cn>
void QNode<T,cn>::LightingCorrection2(const Tensor<T,cn>& fromLight, const Tensor<T,cn>& toLight)
{
	/// only work for 2D case
	Tensor<T,cn>  candid= this->ExtendBoundary(Size3(overlapSize.height,overlapSize.width,overlapSize.depth));
	candid.SetBlock(upBound);
	candid.SetBlock(leftBound);
	candid.SetBlock(Point3i(0,this->tsSize.width+overlapSize.width,0),rightBound);
	candid.SetBlock(Point3i(this->tsSize.height+overlapSize.height,0,0),lowBound);
	CV_Assert(fromLight.size().width == toLight.size().width && fromLight.size().height == toLight.size().height);
	CV_Assert(fromLight.size().width == candid.size().width && fromLight.size().height == candid.size().height);
	CV_Assert(candid.type()!=0); // if use unsigned type, wrong answers!
	//candid = candid.ExtendHalfBoundary();
	//candid.Display(0);
	//fromLight.Display(0);
	//toLight.Display(0);
	//QNode<T,cn> extNode(
	candid = (candid - fromLight) + toLight;
	//candid.Display(0);
	this->SetBlock(candid.Crop(overlapSize.Point3(),this->tsSize));
	this->upBound.SetBlock(candid.Crop(Point3i(0,0,0),upBound.size()));
	this->leftBound.SetBlock(candid.Crop(Point3i(0,0,0), leftBound.size()));
	this->lowBound.SetBlock(candid.Crop(Point3i(this->tsSize.height+overlapSize.height,0,0),lowBound.size()));
	this->rightBound.SetBlock(candid.Crop(Point3i(0,this->tsSize.width+overlapSize.width,0),rightBound.size()));
	return;
}
template <class T, size_t cn>
QNode<T,cn>& QNode<T,cn>::LightingCorrection2(const QNode<T,cn>& changeTo)
{
	//this is different to the lighting correction in tensor because qNode also need to include 4 boundaries
	Tensor<T,cn>  candid= this->ExtendBoundary(Size3(overlapSize.height,overlapSize.width,overlapSize.depth));
	candid.SetBlock(upBound);
	candid.SetBlock(leftBound);
	candid.SetBlock(Point3i(0,this->tsSize.width+overlapSize.width,0),rightBound);
	candid.SetBlock(Point3i(this->tsSize.height+overlapSize.height,0,0),lowBound);
	//Tensor<T,cn> lighting = candid.LightingCorrection(changeTo);
	candid = candid.ExtendHalfBoundary();
	candid.Display(0);
	Tensor<T,cn> target = changeTo.ExtendBoundary(Size3(overlapSize.height,overlapSize.width,overlapSize.depth));
	target.SetBlock(changeTo.GetBoundaryUp());
	target.SetBlock(changeTo.GetBoundaryLeft());
	target = target.ExtendHalfBoundary();
	target.Display(0);

	this->SetBlock(candid.Crop(overlapSize.Point3(),this->tsSize));
	this->upBound.SetBlock(candid.Crop(Point3i(0,0,0),upBound.size()));
	this->leftBound.SetBlock(candid.Crop(Point3i(0,0,0), leftBound.size()));
	this->lowBound.SetBlock(candid.Crop(Point3i(this->tsSize.height+overlapSize.height,0,0),lowBound.size()));
	this->rightBound.SetBlock(candid.Crop(Point3i(0,this->tsSize.width+overlapSize.width,0),rightBound.size()));
	return *this;
}
template <class T, size_t cn> QNode<T,cn> QNode<T,cn>::Clone() const
{

	//QNode<T,cn> rst(this->tsSize);
	QNode<T,cn> rst(this->Tensor<T,cn>::Clone());
	//rst = Tensor<T,cn>::Clone();
	rst.overlapSize = overlapSize;
	rst.tsOffset = this->tsOffset;
	rst.bits = bits;
	rst.predicted_method = predicted_method;
	rst.footError = footError;
	rst.nodeFoot = nodeFoot;
	rst.approxFoot = approxFoot;
	rst.footPos = footPos;
	rst.quanBit = quanBit;
	rst.boundary = boundary;
	rst.qsize = qsize;
	rst.bounds.clear();	
	rst.leftBound = this->leftBound.Clone();
	rst.rightBound = this->rightBound.Clone();
	rst.upBound = this->upBound.Clone();
	rst.lowBound = this->lowBound.Clone();
	rst.bounds.push_back(&rst.upBound);
	rst.bounds.push_back(&rst.leftBound);
	//rst.SetSubWinSize(this->subWinSize);
	//rst.SetSubWinStep(this->subWinStep);
  rst.feet = this->feet;
	//rst.SetBlock(*this);
	return rst;
}
//template <class T, size_t cn> void QNode<T,cn>::SetFootRegion(int region)
//{
//	this->footComputeRegion= region;
//}
//template <class T, size_t cn> void QNode<T,cn>::SetFootComputeMethod(int method)
//{
//	this->footComputeMethod = method;
//}
///explicit instantiation
template class QNode<double,1>;
template class QNode<double,2>;
template class QNode<double,3>;
template class QNode<uchar,1>;
template class QNode<uchar,2>;
template class QNode<uchar,3>;
template class QNode<float,1>;
template class QNode<float,2>;
template class QNode<float,3>;
}

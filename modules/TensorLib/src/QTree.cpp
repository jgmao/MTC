#include "QTree.h"
namespace tensor{
//////////////////////////////////////////////////////////////////////////////
/////////////////////Implementation /////////////////////////////////////////
template<class T, size_t cn>
void QTree<T,cn>::InitEntries(void)
{
	treeLevel=0;
	treePos=0;
	subTrees = vector<QTree<T,cn>>(0);
	subTreeIter = subTrees.begin();
	splitMark=false;
	parentTree = NULL;
	nextTree = NULL;
	branchNum=0;
}

template<class T, size_t cn>
void QTree<T,cn>::InitEntries(const QTree & qt)
{
	treeLevel = qt.treeLevel;
	treePos = qt.treePos;
	splitMark = qt.splitMark;
	parentTree = qt.parentTree;
	nextTree = qt.nextTree;
	branchNum = qt.branchNum;
	subTrees = qt.subTrees;
}


template<class T, size_t cn>
QTree<T,cn>::QTree(void):QNode<T,cn>()
{
	InitEntries();
}

template<class T, size_t cn>
QTree<T,cn>::~QTree(void)
{
}

template<class T, size_t cn>
QTree<T,cn>::QTree(const cv::Mat& mt):QNode<T,cn>(mt)
{
	InitEntries();
}

template<class T, size_t cn>
QTree<T,cn>::QTree(const Tensor<T,cn>& ts):QNode<T,cn>(ts)
{
	InitEntries();
}

template<class T, size_t cn>
QTree<T,cn>::QTree(const QTree& qt):QNode<T,cn>(qt)
{
	InitEntries(qt);
}

template<class T, size_t cn>
QTree<T,cn>::QTree(string cFileName):QNode<T,cn>(cFileName)
{
	InitEntries();
}

template<class T, size_t cn>
QTree<T,cn>::QTree(const Size3& sz, const Vec<T,cn>& val):QNode<T,cn>(sz,val)
{
	InitEntries();
}

template< class T, size_t cn>
QTree<T,cn>::QTree(const QNode<T,cn>& qt):QNode<T,cn>(qt)
{
	InitEntries();
}

template<class T, size_t cn>
QTree<T,cn>& QTree<T,cn>::operator= (const QTree& qt)
{
	QNode<T,cn>::operator=(qt);
	InitEntries(qt);
	return *this;
}


template<class T, size_t cn> 
void QTree<T,cn>::Split(void)
{
	int size = 1<<(int)ceil(log((double)max(this->size().width,this->size().height))/log(2.0));
	int subsize = size/2;
  Size3 suboverlap = this->overlap()/2;
	int count = 0;
	Size3 tempSize;
	Cube roi;
  Cube extroi;
	subTrees = vector<QTree<T,cn>>(4);
  //Tensor<T,cn> test = this->operator()(Cube(-8,-8,0,16,16,1));
  //test.Display();
 	for (int t=0; t<this->size().depth; t=(int)ceil((float)t+float(this->size().depth)/2))
		for (int x=0; x<this->size().height; x+= subsize)
			for (int y=0; y<this->size().width; y+=subsize)
			{	
				tempSize.width = subsize;
				tempSize.height = subsize;
				tempSize.depth = int(ceil(float(this->size().depth)/2));

				if (x+tempSize.height>this->size().height)
					tempSize.height = this->size().height - x;
				if (y+tempSize.width> this->size().width)
					tempSize.width = this->size().width - y;
				if (t+tempSize.depth > this->size().depth)
					tempSize.depth = this->size().depth - t;

				roi = Cube(x,y,t,tempSize.height,tempSize.width,tempSize.depth);
	      extroi = Cube(Point3i(x,y,t)+this->overlap().Point3(),tempSize);
	if (this->eNode.offset().x==0)//deal with the border
	  extroi.x-=this->overlap().height;
	if (this->eNode.offset().y==0)
	  extroi.y-=this->overlap().width;
       
        subTrees[count].Ref(this->eNode,extroi,suboverlap);
        subTrees[count].tsOffset = this->tsOffset+roi.offset();
				//subTrees[count] = this->operator()(roi);
        //subTrees[count].upBound = this->upBound(Cube(suboverlap.Point3(),suboverlap+Size3(0,subsize,1)));
				//subTrees.push_back(Crop(cv::Point3i(x,y,t),tempSize));
				subTrees[count].SetTreeLevel(GetTreeLevel()+1);
				subTrees[count].SetTreePos(count);
				subTrees[count].parentTree = this;
				//	subTrees[count].cFileName = cFileName;
				count++;

			}
			//assign pointers 
			for(subTreeIter= GetSubTreeBeginning(); subTreeIter!=GetSubTreeEnding()-1; subTreeIter++)
			{
				subTreeIter->nextTree = &*(subTreeIter+1);
			}

			(GetSubTreeEnding()-1)->nextTree = this->nextTree;

			branchNum=count;
			subTreeIter = GetSubTreeBeginning();
			splitMark = true;
}

template< class T, size_t cn>
void QTree<T,cn>::Split(bool direction)
{
	//use specially to split 2x2 to two rectangular blocks in PQI
	subTrees = vector<QTree<T,cn>>(2);
	Cube roi;
	if (direction == 0) //vertical
	{
		roi = Cube(0,0,0,2,1,1);
		//subTrees.push_back(Crop(cv::Point3i(0,0,0),Size3(2,1,1)));
		subTrees[0] = this->operator()(roi);
		subTrees[0].SetTreeLevel(GetTreeLevel()+1);
		subTrees[0].SetTreePos(0);
		subTrees[0].parentTree = this;
		roi = Cube(0,1,0,2,1,1);
		//subTrees.push_back(Crop(cv::Point3i(0,1,0),Size3(2,1,1)));
		subTrees[1] = this->operator()(roi);
		subTrees[1].SetTreeLevel(GetTreeLevel()+1);
		subTrees[1].SetTreePos(1);
		subTrees[1].parentTree = this;

	}
	else // horizantal
	{
		roi = Cube(0,0,0,1,2,1);
		subTrees[0] = this->operator()(roi);
		//subTrees.push_back(Crop(cv::Point3i(0,0,0),Size3(1,2,1)));
		subTrees[0].SetTreeLevel(GetTreeLevel()+1);
		subTrees[0].SetTreePos(0);
		subTrees[0].parentTree = this;
		roi = Cube(1,0,0,1,2,1);
		subTrees[1] = this->operator()(roi);
		//subTrees.push_back(Crop(cv::Point3i(1,0,0),Size3(1,2,1)));
		subTrees[1].SetTreeLevel(GetTreeLevel()+1);
		subTrees[1].SetTreePos(1);
		subTrees[1].parentTree = this;

	}

	for(subTreeIter= GetSubTreeBeginning(); subTreeIter!=GetSubTreeEnding()-1; subTreeIter++)
	{
		subTreeIter->nextTree = &*(subTreeIter+1);
	}

	(GetSubTreeEnding()-1)->nextTree = this->nextTree;

	branchNum=2;
	subTreeIter = GetSubTreeBeginning();
	splitMark = true;

}

//template< class T, size_t cn> QTree<T,cn> QTree<T,cn>::Clone() const
//{
//	QTree<T,cn> rst(tsSize);
//	rst.treeLevel = treeLevel;
//	rst.treePos = treePos;
//	rst.splitMark = splitMark;
//	rst.parentTree = parentTree;
//	rst.rst.nextTree = nextTree;
//	rst.branchNum = branchNum;
//	rst.subTrees = subTrees;
//	rst.SetBlock(this->Tensor<T,cn>::Clone());
//}

template< class T, size_t cn>
QTree<T,cn> * QTree<T,cn>::NextLeaf(void)
{

//	QTree<T,cn> * tempTree;
  /*
	if (this->splitMark)
	{	tempTree= &*GetSubTree(0);
	if (tempTree->splitMark)
		return tempTree->NextLeaf();
	else
	{

		return tempTree;
	}
	}
	else 
  */
  if (nextTree == NULL)
		return NULL;
	else if ( nextTree->splitMark)
	{
		return nextTree->NextLeaf();
	}
	else
		return nextTree;
}

template< class T, size_t cn>
typename QTree<T,cn>::qiter QTree<T,cn>::GetSubTree(int pos) 
{
	if (pos >= (int)subTrees.size())
	{
		cout<<"error, pos exceed sub tree size"<<endl;
		exit(0);
	}
	else 
	{
		return subTrees.begin()+pos;
	}
}


template< class T, size_t cn>
int QTree<T,cn>::GetTreeLevel(void) const
{
	return this->treeLevel;
}

template< class T, size_t cn>
int QTree<T,cn>::GetTreePos(void) const
{
	return this->treePos;
}

template< class T, size_t cn>
typename QTree<T,cn>::qiter QTree<T,cn>::GetSubTreeBeginning(void) 
{
	return subTrees.begin();
}

template< class T, size_t cn>
typename QTree<T,cn>::qiter QTree<T,cn>::GetSubTreeEnding(void) 
{
	return subTrees.end();
}

template< class T, size_t cn>
void QTree<T,cn>::SetTreeLevel(int lv)
{
	this->treeLevel = lv;
}

template< class T, size_t cn>
void QTree<T,cn>::SetTreePos(int pos)
{
	this->treePos = pos;
}
template< class T, size_t cn>
int QTree<T,cn>::GetSubTreeSize(void) const
{
	return int(subTrees.size());
}

template< class T, size_t cn>
int QTree<T,cn>::GetPeerTreeSize(void) const
{
	if (parentTree!=NULL)
		return parentTree->GetSubTreeSize();
	else
		exit(0);
}
template <class T, size_t cn>
int QTree<T,cn>::CollectBits(void) const
{
  int rst=0;
  if (this->splitMark)
    rst+=bitcount.GetDecBitN();
  else
  {
    rst+=bitcount.GetDecBitN();
    rst+=bitcount.GetFootBitN();
    rst+=bitcount.GetCodingBitN();
  }
  if (this->splitMark)
  {
    for (int i=0; i<this->GetSubTreeSize();i++)
      rst+= this->subTrees[i].CollectBits();
  }
  return rst;
}

template <class T, size_t cn> 
int QTree<T,cn>::CollectTreeBits(void) const
{
  int rst=0;
  if (this->splitMark)
  {
    rst+=bitcount.GetDecBitN();
    for (int i=0; i<this->GetSubTreeSize();i++)
      rst+= this->subTrees[i].CollectTreeBits();
  }
  return rst;
}
template <class T, size_t cn> 
int QTree<T,cn>::CollectModeBits(void) const
{
  int rst=0;
  if (!this->splitMark)
    rst+=bitcount.GetDecBitN();
  else
  {
    for (int i=0; i<this->GetSubTreeSize();i++)
      rst+= this->subTrees[i].CollectModeBits();
  }
  return rst;
}
template <class T, size_t cn> 
vector<int> QTree<T,cn>::CollectMTCBits(int level) const
{
  vector<int> rst(int(log(this->tsSize.height)/log(2.0))-3);
  if (this->splitMark)
  {

    for (int i=0; i<this->GetSubTreeSize();i++)
    {
      vector<int> temp= this->subTrees[i].CollectMTCBits(level+1);
      for (int j=0; j< temp.size(); j++)
      {
        rst[j+1]+= temp[j];
      }
    }
  }
  else
  {
    if (this->bitcount.GetCodingMethod()==CodingMethodNames::CODING_MTC)
      rst[level]+=bitcount.GetCodingBitN();
  }
  return rst;
}

template <class T, size_t cn> 
vector<int> QTree<T,cn>::CollectJPEGBits(void) const
{
  vector<int> rst(2);
  if (this->splitMark)
  {
    for (int i=0; i<this->GetSubTreeSize();i++)
    {
      vector<int> temprst = this->subTrees[i].CollectJPEGBits();
      rst[0]+=temprst[0];
      rst[1]+=temprst[1];
    }
  }
  else
  {
    if (this->bitcount.GetCodingMethod()==CodingMethodNames::CODING_JPEG)
      rst[0]+=bitcount.GetCodingBitN();
    else if(this->bitcount.GetCodingMethod()==CodingMethodNames::CODING_JPEG_DEGRADE)
      rst[1]+=bitcount.GetCodingBitN();
  }
  return rst;
}
template <class T, size_t cn> 
int QTree<T,cn>::CollectPQIBits(void) const
{
  int rst=0;
  if (this->splitMark)
  {
    for (int i=0; i<this->GetSubTreeSize();i++)
    {
      rst+=this->subTrees[i].CollectPQIBits();
    }
  }
  else
  {
    if (this->bitcount.GetCodingMethod()==CodingMethodNames::CODING_PQI)
      rst+=bitcount.GetCodingBitN();
  }
  return rst;
}
template <class T, size_t cn> 
vector<int> QTree<T,cn>::GetMTCN(int level) const
{
  vector<int> rst(1);
  if (this->splitMark)
  {
    rst.push_back(0);
    for (int i=0; i<this->GetSubTreeSize();i++)
    {
      vector<int> temp= this->subTrees[i].GetMTCN(level+1);
      for (int j=0; j< temp.size(); j++)
      {
        rst[j+1]+= temp[j];
      }
    }
  }
  else
  {
    if (this->bitcount.GetCodingMethod()==CodingMethodNames::CODING_MTC)
      rst[level]++;
  }
  return rst;
}

template <class T, size_t cn> 
vector<int> QTree<T,cn>::GetJPEGN(void) const
{
  vector<int> rst(2);
  if (this->splitMark)
  {
    for (int i=0; i<this->GetSubTreeSize();i++)
    {
      vector<int> temprst = this->subTrees[i].GetJPEGN();
      rst[0]+=temprst[0];
      rst[1]+=temprst[1];
    }
  }
  else
  {
    if (this->bitcount.GetCodingMethod()==CodingMethodNames::CODING_JPEG)
      rst[0]++;
    else if(this->bitcount.GetCodingMethod()==CodingMethodNames::CODING_JPEG_DEGRADE)
      rst[1]++;
  }
  return rst;
}
template <class T, size_t cn> 
vector<int> QTree<T,cn>::GetPQIN(int level) const
{
  vector<int> rst(1);
  if (this->splitMark)
  {
    rst.push_back(0);
    for (int i=0; i<this->GetSubTreeSize();i++)
    {
      vector<int> temp= this->subTrees[i].GetMTCN(level+1);
      for (int j=0; j< temp.size(); j++)
      {
        rst[j+1]+= temp[j];
      }
    }
  }
  else
  {
    if (this->bitcount.GetCodingMethod()==CodingMethodNames::CODING_PQI)
      rst[level]++;
  }
  return rst;
}
//explicit instantiation
template class QTree<double,1>;
template class QTree<double,2>;
template class QTree<double,3>;
template class QTree< uchar,1>;
template class QTree< uchar,2>;
template class QTree< uchar,3>;
template class QTree<float,1>;
template class QTree<float,2>;
template class QTree<float,3>;
}

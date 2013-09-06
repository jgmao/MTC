#include "QGrid.h"
#include <future>
/////////////////////////////// implementation ////////////////

namespace tensor{
template< class T, size_t cn>
QGrid<T,cn>::~QGrid(void)
{
}


template <class T, size_t cn> 
QGrid<T,cn>::QGrid(void )
{
	cFileName = "Unknown";
}
template <class T, size_t cn>
QGrid<T,cn>::QGrid(const string& cFileName)
{
	this->cFileName = cFileName;
	ensemble.Load(cFileName);
}
template<class T, size_t cn>
QGrid<T,cn>& QGrid<T,cn>::operator=(const QGrid& g)
{
	cFileName = g.cFileName;
	blockSize = g.blockSize;
	overlapSize = g.overlapSize;
	gridSize = g.gridSize;
	ensemble = g.ensemble;
	causalMap = g.causalMap;
	grid = g.grid;
	searchStep = g.searchStep;
	candidNum = g.candidNum;
	return *this;
}
template<class T, size_t cn>
QGrid<T,cn>::QGrid(const QGrid& g)
{
	*this = g;
}


template<class T, size_t cn>
QGrid<T,cn>::QGrid(const string& cFileName, const Size3& blockSize, const Size3& overlapSize)
{
	ensemble.Load(cFileName);
	this->cFileName = cFileName;
	this->blockSize = blockSize;
	this->overlapSize = overlapSize;
	searchStep = Size3(4,4,1);
	candidNum = 4;
	InitGrid();
}


template<class T, size_t cn>
QGrid<T,cn>::QGrid(const string& cFileName, const Size3& blockSize, const Size3_<double>& overlapRatio)
{
	this->overlapSize = Size3(Size3_<double>(blockSize)*overlapRatio);
	ensemble.Load(cFileName);
	this->cFileName = cFileName;
	this->blockSize = blockSize;
	searchStep = Size3(4,4,1);
	candidNum = 4;
	InitGrid();
}

template<class T, size_t cn>
void QGrid<T,cn>::InitGrid(void)
{
  this->gridSize = (Size3_<double>(ensemble.size())/Size3_<double>(blockSize)).Ceil();
	grid = Grid(gridSize.depth,vector<vector<QTree<T,cn>>>(gridSize.height));
	rst = Tensor<T,cn>(ensemble.size());
  /*
  rstExt = Tensor<T,cn>(overlapSize*4 + ensemble.size());
  for ( int z = 0; z < rst.size().depth; z++)
  {
    Mat temp;
    cv::copyMakeBorder(rst.GetFrame(z),temp,overlapSize.height*2,overlapSize.height*2,overlapSize.width*2,overlapSize.width*2,cv::BORDER_REFLECT);
    rstExt.SetFrame(z,temp);
  }
  rstExt.Display(); 
  rstExt.SetSubWinSize(ensemble.GetSubWinSize());
	rstExt.SetSubWinStep(ensemble.GetSubWinStep());
  */
        //rst.SetSubWinSize(ensemble.GetSubWinSize());
        //rst.SetSubWinStep(ensemble.GetSubWinStep());

	rstSqr = Tensor<T,cn>(ensemble.size());
	SetVarThrd();
	grid = vector<vector<vector<QTree<T,cn>>>>(gridSize.depth);
	for (int t=0; t<gridSize.depth; t++)
		grid[t]= vector<vector<QTree<T,cn>>>(gridSize.height,vector<QTree<T,cn>>(gridSize.width));
	localVarMap = Tensor<double,cn>(ensemble.size()/searchStep);//*searchStep.height,gridSize.width*searchStep.width,gridSize.depth*searchStep.depth);
	//cout<<"done 4\n";
	for (int t=0;t<gridSize.depth;t++)
		for (int x=0;x<gridSize.height;x++)
			for (int y=0;y<gridSize.width;y++)
			{
				localVarMap(x,y,t)=0;
				//QNode<T,cn>(rst,
				grid[t][x][y].Ref(rst,Cube(Point3i(x*blockSize.height,y*blockSize.width,t*blockSize.depth),blockSize),overlapSize);
			}
			//grid[0][8][8].SetBlock(Tensor_<T,1>(blockSize,Vec<T,cn>(255)));
			//rst.Display();
			causalMap = Tensor<T,1>(ensemble.size());
			return;
}

template<class T, size_t cn>
void QGrid<T,cn>::ReInitGrid(void)
{
	InitGrid();
}

template<class T, size_t cn>
void QGrid<T,cn>::SetCausalMap(const QNode<T,cn>& nd)
{
	Tensor<T,1> temp(nd.size(),Vec<T,1>::all(255));
	causalMap.SetBlock(nd.offset(),temp);
}

template<class T, size_t cn>
bool QGrid<T,cn>::IsInsideCausalRegion(const Point3i& curPos, const QNode<T,cn>& nd, const int multiBound) const
{
	
	//must gurantee do not exceed CausalMap before use it.
	if (curPos.x<0||curPos.y<0||curPos.z<0)
		return false;
	Size3 absBoundary = ensemble.size();
	//for post blending, overlap*2
	//to avoid boundary overlap add additional 1 overlap size

	Size3 temp = nd.overlap()*multiBound+ nd.size() + curPos - Size3(1,1,1);
	Point3i tempPos(temp.height,temp.width,temp.depth);
	if (tempPos.x>=absBoundary.height ||
		tempPos.y>=absBoundary.width||
		tempPos.z>=absBoundary.depth)
		return false;
	//this method of check SE corner only work in Raster scan order
	//cout<<tempPos;
	if (causalMap[tempPos][0] > 0)
		return true;
	else
		return false;

}
template<class T, size_t cn>
QTree<T,cn>& QGrid<T,cn>::GetNode(const Point3i& pos) {
	return grid[pos.z][pos.x][pos.y];
}

template<class T, size_t cn>
void QGrid<T,cn>::DisplayCasualMap(int flag)
{
	causalMap.Display(flag);
}

template<class T, size_t cn>
void QGrid<T,cn>::Display(int flag)
{
	rst.Display(flag);
}

template<class T, size_t cn>
vector<Point3i> QGrid<T,cn>::BoundaryMatching(QNode<T,cn>& qNode, MatchingMethod matching_method, double matching_thrd)
{
	//qNode.Display();
	//qNode.leftBound.Display();
	//qNode.upBound.Display();
	double diff=0.0;
	double N = 0.0;
	//this->causalMap.Display();
	//bool findIndicator = false;
	//vector<STBlock> fromRef = GetRefBoundary(refNode,0);
	//vector<STBlock> fromTag;
	//vector<STBlock>::iterator refIter, tagIter;

  ///todo: 20130604//change here, make sure have sufficent/necessary gap /////////////////
  std::mutex queue_lock;

	//for post blending, use extra boundary: overlap *2
	//I should change the gap, since I am not going to use 2x such big block size in PLC
  Size3 searchRegion = rst.size() - qNode.size() - qNode.overlap()*3 + Size3(1,1,1);
  //Size3 searchRegion = rst.size() - qNode.size() - qNode.overlap()*2+ Size3(1,1,1);//20130604 changed , not necessary such big gap?????
  // for limited search region, change here
  double ratio = 5;
  int offsetUp = qNode.offset().x - int(ratio*double(qNode.size().height));
  offsetUp > 1 ? offsetUp = offsetUp: offsetUp=1;
  int offsetLeft = qNode.offset().y - int(ratio*double(qNode.size().width));
  offsetLeft>1? offsetLeft = offsetLeft: offsetLeft = 1;
  int offsetRight = qNode.offset().y + int(ratio*double(qNode.size().width));
  offsetRight>searchRegion.width-1?offsetRight = searchRegion.width-1:offsetRight = offsetRight;
  int offsetDown = qNode.offset().x + int(2*double(qNode.size().height));
  offsetDown>searchRegion.height-1?offsetDown = searchRegion.height-1:offsetDown  = offsetDown;
  LinkStruct* queue;
  //if (candidNum>0)
  //  queue = new LinkQueue(candidNum);
	  //LinkQueue queue(candidNum);//gj15012013  set possibility of dynamic number of candidate
  //else
    queue = new LinkArray();

  LinkQueue varQueue(100);
  vector<Point3i> varList;
	Tensor<T,cn> target;
  Tensor<T,cn> candUp, candLeft;
	Tensor<double,cn> tempMu,tempVar;
	Tensor<double,cn> tempMap;
   vector<Point3i> sideMatchAddr;
  fstream logfile;
  logfile.open("./temp/matching.txt",ios::app);

	if (qNode.offset().x ==416 && qNode.offset().y==352)
		tempMap = Tensor<double,cn>(rst.size());
	if (matching_method == MatchingMethod::MATCHING_MSE || matching_method ==  MatchingMethod::MATCHING_SAD|| matching_method == MatchingMethod::MATCHING_MSE_CONSTRAINT)
	{
	//old method, template matching include boundarys x 3

	if (qNode.offset().x == DEBUG_X && qNode.offset().y == DEBUG_Y && qNode.size().height==DEBUG_SIZE)
	{
		rst.SaveBlock("rst.tif");
		//logfile<<"\n ======== output debug info, template matching ====== \n";
		qNode.debugtrigger = true;
	}
  Tensor<T,cn> tarUp , tarLeft;
  qNode.upBound.Ref(Cube(0,qNode.overlap().width,0,qNode.overlap().height,qNode.size().width,qNode.size().depth),tarUp);
  tarLeft = qNode.leftBound;
  if (qNode.offset().x == DEBUG_X && qNode.offset().y == DEBUG_Y && qNode.size().height==DEBUG_SIZE)
	{
    tarUp.Print("up");
    tarLeft.Print("left");
  }
  Tensor<T,cn> TT(tarUp.size().width+tarLeft.size().height,tarUp.size().height,1);
  
  Vec<T,cn> Mt=0;
  double Vt0=0,Vt1=0,Vt=0;
	double lambda = 10;
  
  //substract mean before matching, 11212012
  //estimate mean by sliding windows 05232013
  Tensor<T,cn> upMean, leftMean;
#ifdef SUBSTRACT_BOUND_MEAN
  ComputeBoundMean(tarUp,tarLeft,qNode.offset(),qNode.size(),qNode.overlap(),upMean,leftMean);
#else
  upMean = Tensor<T,cn>(tarUp.size(),0);
  leftMean = Tensor<T,cn>(tarLeft.size(),0);
#endif
  
#ifndef PARALLEL_ENABLED
  Tensor<T,cn> tarUp_norm = tarUp - upMean;
  Tensor<T,cn> tarLeft_norm = tarLeft- leftMean;
  if (qNode.offset().x == DEBUG_X && qNode.offset().y == DEBUG_Y && qNode.size().height==DEBUG_SIZE)
	{
    tarUp_norm.Print("up_norm");
    tarLeft_norm.Print("left_norm");
  }
#endif

  if (matching_method ==  MatchingMethod::MATCHING_MSE_CONSTRAINT)
  {
    //before doing that, up, left and anti feets must be prepared. 
    TT.SetBlock(tarUp.Transpose());
    TT.SetBlock(Point3i(tarUp.size().width,0,0),tarLeft);
    //Mt = TT.Mean();
    Vt = TT.Var()[0];
  }

	//Tensor<T,cn> T0 =  qNode.GetBoundary(0).Crop(Point3i(0,qNode.overlap().width,0),qNode.GetBoundary(0).size()-Size3(0,qNode.overlap().width,0));
	//Tensor<T,cn> T1 = 	 qNode.GetBoundary(1).Clone();
	//Tensor<T,cn> TT(Size3(T0.size().height,T0.size().width+T1.size().width,1));
	//Tensor<T,cn> VV(Size3(T0.size().height,T0.size().width+T1.size().width,1));
	//TT.SetBlock(T1.Transpose());
	//TT.SetBlock(Point3i(0,T1.size().width,0),T0);
	//cv::Mat gaussKernel = mylib::GenGaussKer(3,double(3)/6.0,CV_64F);
	//Tensor<T,cn> T0_low=T0.Filter2D(gaussKernel,FILTER_BOUND_VALID);
	//Tensor<T,cn> T1_low=T1.Filter2D(gaussKernel,FILTER_BOUND_VALID);
	//Tensor<T,cn> T0_Lap=T0.Laplacian();
	//Tensor<T,cn> T1_Lap=T1.Laplacian();
	Tensor<T,cn> filterBD,lowBD;
	//TT.Print();
	/*if (qNode.debugtrigger)
	{
		T0.Print();
		T1.Print();
	}*/

  
  
#ifndef PARALLEL_ENABLED
  /// comment them when multi-thread /////////////// 
  Tensor<T,cn> candUp_norm, candLeft_norm;
  Tensor<T,cn> VV(tarUp.size().width+tarLeft.size().height,tarUp.size().height,1);
      double Vb0=0,Vb1=0,Vb=0;
      Vec<T,cn> Mb=0; 
  ///// comment above when multi-thread //////
#endif



  //do hihgpass before matching
  Mat ker = (Mat_<double>(3,3) <<-0.0113437365584951,	-0.0838195058022106,	-0.0113437365584951,
                                 -0.0838195058022106,	0.380652969442823,	-0.0838195058022106,
                                 -0.0113437365584951,	-0.0838195058022106,	-0.0113437365584951);
  //T0 = T0.Filter2D(ker,cv::BORDER_REFLECT);
  //T1 = T1.Filter2D(ker,cv::BORDER_REFLECT);
  //
  //std::thread threads[1];//[(offsetDown-offsetUp)/searchStep.height];
#ifdef PARALLEL_ENABLED
  int pnum = 4;
  vector<thread> threads;
#else
  int pnum =1;// thread::hardware_concurrency();
#endif
  //int blines = (offsetDown-offsetUp)/searchStep.height/pnum;
  int brows = (offsetRight - offsetLeft)/searchStep.width/pnum;
  vector<LinkStruct*> queues(pnum);
  //for (auto& q : queues)
  //  q = new LinkArray();
	for (int t=0; t< searchRegion.depth;t+=searchStep.depth)
  {
    for (int p=0; p<pnum;p++)
    {
      double localdiff=0.0;//make it local
      double localN=0.0;
      queues[p]=new LinkArray();
#ifdef PARALLEL_ENABLED
      threads.push_back(thread([&](int p, int t, Tensor<T,cn> localTarUp, Tensor<T,cn> localTarLeft, LinkStruct* myqueue ){
     // cout<<"work in p: "<< p<<","<<  std::this_thread::get_id()<<endl;
      std::mutex rst_mutex;
      //auto localrst = rst.Clone();//
      Tensor<T,cn> localCandUp, localCandLeft; //make them local
      Tensor<T,cn> candUp_norm, candLeft_norm; //make them local
      Tensor<T,cn> VV(tarUp.size().width+tarLeft.size().height,tarUp.size().height,1);
      QNode<T,cn> localCandid;
      double Vb0=0,Vb1=0,Vb=0;
      Vec<T,cn> Mb=0;
      Tensor<T,cn> tarUp_norm = localTarUp-upMean;
      Tensor<T,cn> tarLeft_norm = localTarLeft- leftMean;
#endif
		  //for (int x=offsetUp+p*blines; x< min(offsetDown,offsetUp+(p+1)*blines); x+=searchStep.height)
      for (int y = offsetLeft +p*brows; y<min(offsetRight, offsetLeft+(p+1)*brows);y+=searchStep.width)
      {
          for (int x=offsetUp; x<= offsetDown; x+=searchStep.height)
			    //for (int y=offsetLeft; y< offsetRight; y+=searchStep.width)
			    {
				    if (IsInsideCausalRegion(cv::Point3i(x,y,t),qNode,3/*3*/))//20130605  change from 3 to 2
				    {
              //rearrage this part, it looks so arkward. 
              if (qNode.GetBoundarySize()!=2)
                CV_Error(CV_StsBadSize,"boundary missing , must hav both up and left\n");
              //localCandLeft = localrst.GetBlock(Cube(Point3i(x,y,t),qNode.leftBound.size()));
              //localCandUp = localrst.GetBlock(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)));
              Tensor<T,cn> candUpMean, candLeftMean;
#ifdef PARALLEL_ENABLED
              {
              mymux.lock();
              localCandLeft = rst.Crop(Point3i(x,y,t),qNode.leftBound.size());
              localCandUp = rst.Crop(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
              mymux.unlock();
              //lock_guard<std::mutex> lock(rst_mutex);
              //rst.Ref(Cube(Point3i(x,y,t),qNode.leftBound.size()),localCandLeft);
              //rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),localCandUp);
#ifdef SUBSTRACT_BOUND_MEAN
              ComputeBoundMean(localCandUp,localCandLeft,qNode.offset(),qNode.size(),qNode.overlap(),candUpMean,candLeftMean);
#else
              candUpMean = Tensor<T,cn>(localCandUp.size(),0);
              candLeftMean = Tensor<T,cn>(localCandLeft.size(),0);
#endif
  
              }
#else
              rst.Ref(Cube(Point3i(x,y,t),qNode.leftBound.size()),candLeft);
              rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),candUp);
              if (qNode.offset().x == DEBUG_X && qNode.offset().y == DEBUG_Y && qNode.size().height==DEBUG_SIZE)
	            {
                if (x==29&&y==52)
                {
                candUp.Print("cand_up");
                candLeft.Print("cand_left");
                }
              }
#ifdef SUBSTRACT_BOUND_MEAN
              ComputeBoundMean(candUp,candLeft,qNode.offset(),qNode.size(),qNode.overlap(),candUpMean,candLeftMean);
#else
              candUpMean = Tensor<T,cn>(candUp.size(),0);
              candLeftMean = Tensor<T,cn>(candLeft.size(),0);
#endif
  
#endif
                                            if (matching_method ==  MatchingMethod::MATCHING_MSE_CONSTRAINT)
              {
#ifdef PARALLEL_ENABLED
                VV.SetBlock(localCandUp.Transpose());
                VV.SetBlock(Point3i(tarUp.size().width,0,0),localCandLeft);
#else
                VV.SetBlock(candUp.Transpose());
                VV.SetBlock(Point3i(tarUp.size().width,0,0),candLeft);
#endif
                //Mb = VV.Mean();
                Vb = VV.Var()[0];
              }              
#ifdef PARALLEL_ENABLED
              candUp_norm = localCandUp - candUpMean;
              candLeft_norm = localCandLeft - candLeftMean;
#else
              candUp_norm = candUp - candUpMean;
              candLeft_norm = candLeft - candLeftMean;
#endif
	      if (matching_method ==  MatchingMethod::MATCHING_MSE || matching_method== MatchingMethod::MATCHING_MSE_CONSTRAINT)
					    {  //diff+= (T1_low.Compare(lowBD,COMPARE_CRITERIA_MSE)+T1_Lap.Compare(filterBD,COMPARE_CRITERIA_MSE));
						    localdiff = metric::Compare(tarUp_norm,candUp_norm,CompareCriteria::COMPARE_CRITERIA_MSE) +
			   metric::Compare(tarLeft_norm, candLeft_norm, CompareCriteria::COMPARE_CRITERIA_MSE);

                //else if (matching_method == MATCHING_MSE_CONSTRAINT) //	diff+= qNdoe.GetBoundary(i).Compare(target,COMPARE_CRITERIA_MSE_CONSTRAINT);
              }
              else if (matching_method == MatchingMethod::MATCHING_SAD)
					    {
		    localdiff+= metric::Compare(tarUp_norm,candUp_norm,CompareCriteria::COMPARE_CRITERIA_SAD)+
		   metric::Compare(tarLeft_norm,candLeft_norm,CompareCriteria::COMPARE_CRITERIA_SAD);
              }
              else
						    CV_Error(CV_StsUnsupportedFormat,"unsupport criteria\n");
					    localN= candUp_norm.size().volumn()+candLeft_norm.size().volumn();

					    //target.debugtrigger=false;
					    //localdiff= localdiff/localN;
	      localdiff/=65025;//normalize by 255^2
					    if (matching_method == MatchingMethod::MATCHING_MSE_CONSTRAINT)
					    {
						    //double delta = ((Vt0 + Vt1)/2 - (Vb0+Vb1)/2);
						    //VV.Print();
						    //double delta = ((Vt0-Vb0)*(Vt0-Vb0)+(Vt1-Vb1)*(Vt1-Vb1))/2;
						    //double delta = abs(Vb-Vt);
						    //lambda = 1.0/log(delta);
						    //lambda = 10;
                //double d0 = abs(Vb0-Vt0);
                //double d1 = abs(Vb1-Vt1);
						    double delta = abs(log(Vb/Vt));
                //normalize mse, since mse is bounded to 255^2, done above gj20130116
                //diff = diff/255/255;
                //clip delta and normalize it
                double delta_thrd = 100;
                delta > delta_thrd? delta = 1 : delta=delta/delta_thrd;
                //d1 > delta_thrd? d1 = 1 : d1 = d1/delta_thrd;
                //d0 > delta_thrd? d0 = 1 : d0 = d0/delta_thrd;
                lambda = 0.5;
                double lambda2 = 0.5;
                double masking=0;
                if (Mt[0]<50)
                  masking = 0.5;
                else if (Mt[0]>200)
                  masking =0.5;
                else
                  masking = 1;
                  //double constraint = lambda*(d0*d0 + d1*d1);//may be constraint by distribution similarity
						    double constraint = lambda*delta;// + masking*lambda2 ;
            
                localdiff = localdiff + constraint ;
					    }
              /// debug here !!!!!!!!!!!!!
              //protect shared data here
              #ifdef PARALLEL_ENABLED
              mymux.lock();
#endif
              if (localdiff<=matching_thrd)
					      /*myqueue*/queue->compareInsert(-localdiff,cv::Point3i(x,y,t));
              if (/*myqueue*/queue->getLength()>candidNum&&candidNum>0)
                /*myqueue*/queue->pop();
              #ifdef PARALLEL_ENABLED
              mymux.unlock();
#endif
				    }
				    else
				    {
					    //cout<<x<<","<<y<<","<<t<<endl;
				    }
			    }
        }
        //queues[n]=localqueue;
        //p.set_value(localqueue);
#ifdef PARALLEL_ENABLED
      },p,t,tarUp, tarLeft, queues[p]));
#endif
      //cout<<"push at : "<<std::this_thread::get_id()<<endl;
      //threads.push_back(searchaline);
    }
  }
#ifdef PARALLEL_ENABLED
  for (auto& thread : threads){
    thread.join();
  }
#endif
	if (qNode.offset().x == DEBUG_X && qNode.offset().y == DEBUG_Y && qNode.size().height==DEBUG_SIZE)
	{
				//logfile.close();
	}
	qNode.debugtrigger=false;

 /* logfile<<"PRINT ALL IN THE PARALLEL QUEUES: "<<qNode.offset().x<<","<<qNode.offset().y<<qNode.size().height<<endl;
  for (auto& q : queues)
  {
    if (q->getLength()>0)
    {
      logfile<<"    ";
      for (int qi=0; qi<q->getLength(); qi++)
        logfile<<q->getData(qi)<<",";
      logfile<<endl;
    } 
  }*/
  ///combine all the queues
  //for (int iii=0; iii<candidNum; iii++)
  //{
  //int ct=0;
  //int cidx=-1;
  //double tempd = DBL_MAX;
  //CandidateRecord temp;
  //CandidateRecord itis;
  ////print all in the queues
  //
  //for (auto& q : queues)
  //{
  //  if (q->getLength()>0)
  //  {
  //   temp = q->getMin();
  //   if(temp.data<tempd)
  //   {
  //     itis = temp;
  //     cidx =ct;
  //     tempd = itis.data;
  //   }
  //  }
  //  ct++;
  //}
  //if (cidx>=0)
  //{
  //  queues[cidx]->pop();
  //  queue->compareInsert(itis.data,itis.addr);
  //}
  //}
  
  for (auto&q : queues)
  {
    delete q;
  }
    //  for (int iq=0; iq<q->getLength(); iq++)
  //  {
  //   queue->compareInsert(q->getData(iq),q->GetAddress(iq));
  //   if (queue->getLength()>candidNum&&candidNum>0)
  //              queue->pop();
  //  }
  


	}
	//old method finish here 
	else if (matching_method == MatchingMethod::MATCHING_SAT)
	{

	//qNode.Print();
	}
	else if (matching_method == MatchingMethod::MATCHING_VAR)
	{
		Tensor<T,cn> tempBlk = ensemble.Crop(qNode.offset(),qNode.size());
		tempMu = tempBlk.LocalMean(searchStep,searchStep);
		tempVar = tempBlk.LocalVariance(tempMu,searchStep,searchStep);
		localVarMap.SetBlock(Point3i(qNode.offset().x/searchStep.height,qNode.offset().y/searchStep.width,qNode.offset().z/searchStep.depth),tempVar);
	
	//localVarMap.Print();
	//coarse seaching 32x32 by variance with Z8
	for (int t=0; t< searchRegion.depth;t+=searchStep.depth)
		for (int x=offsetUp; x< offsetDown; x+=searchStep.height)
			for (int y=offsetLeft; y< offsetRight; y+=searchStep.width)
			{
				//because next step will be seaching in 8x8 neighbor , so half of searchStep back is used to make sure all possible position will be searched
				if (IsInsideCausalRegion(cv::Point3i(x-searchStep.height/2,y-searchStep.width/2,t-searchStep.depth/2),qNode,0))
				{
					diff = (localVarMap(Cube(Point3i(x/searchStep.height,y/searchStep.width,
						t/searchStep.depth),qNode.size()/searchStep)) - tempVar).Abs().Sum()[0];
					
					queue->compareInsert(diff,Point3i(x,y,t));
				}
			}
	
			//only keep the valid position;
			while(true)
			{
				if (queue->getLength()>0)
				{
					if (queue->GetAddress()[0]==Point3i(-1,-1,-1))
						queue->pop();
					else
						break;
				}
				else
					break;
			}

			//finer search around the valid coarse positions
			Point3i tempPos;
			double tempDiff;
			for (int k=0; k<queue->getLength();k++)
			{
				tempPos = queue->GetAddress()[k];
				tempDiff = DBL_MAX;
				queue->set(k,tempDiff,Point3i(-1,-1,-1));//eliminate first to guantee no non-causal candidate exist
				for (int t=0; t<1; t++)
					for (int x=(int)-searchStep.height/2; x<(int)searchStep.height/2; x++)
						for (int y=(int)-searchStep.width/2; y<(int)searchStep.width/2;y++)
						{
							diff=0;
							N=0;
							if (IsInsideCausalRegion(tempPos+Point3i(x,y,t)-qNode.overlap().Point3(),qNode,3))
							{
									
								for (int i = 0; i< qNode.GetBoundarySize(); i++)
								{
									rst.Ref(Cube(Point3i(x,y,t)+tempPos-qNode.overlap().Point3(),qNode.GetBoundary(i).size()),target);
									diff+= metric::Compare(qNode.GetBoundary(i),target,CompareCriteria::COMPARE_CRITERIA_MSE);
									N+= qNode.GetBoundary(i).size().volumn();
								}
								if (diff<tempDiff)
									queue->set(k,diff,tempPos+Point3i(x,y,t)-qNode.overlap().Point3());
								//diff= diff/N;
								
							}
						}
				if (diff > varThrd1)//use varThrd1 instead of matching_thrd
					queue->set(k,diff,Point3i(-1,-1,-1));
				//else
					//queue->setData(k,diff);
			}
			//filter again the invalid candidate
	}
	else if (matching_method == MatchingMethod::MATCHING_HIERARCHY)
	{
	  Tensor<T,cn> T_left,B_left; 
	  Tensor<T,cn> T_up, B_up;
    Tensor<T,cn> TT(Size3(qNode.overlap().height,2*qNode.size().width+qNode.overlap().width,1));
	  Tensor<T,cn> VV(Size3(qNode.overlap().height,2*qNode.size().width+qNode.overlap().width,1));
    double var_blk,var_tag;
    bool has_left=false, has_up=false;
    double var_tag_left=-1,var_tag_up=-1, var_blk_left=-1, var_blk_up=-1;


    Vec<T,cn> muT,muB;
    Tensor<T,cn> T_left_norm, T_up_norm, B_up_norm,B_left_norm;
    if (qNode.offset().x > 0  && qNode.offset().y>0)
    {
      T_up =  qNode.upBound.Crop(Point3i(0,qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
      T_left = 	 qNode.leftBound.Clone();
      TT.SetBlock(T_left.Transpose());
	    TT.SetBlock(Point3i(0,T_left.size().width,0),T_up);
      var_tag = TT.Var()[0];
      var_tag_left = T_left.Var()[0];
      var_tag_up = T_up.Var()[0];
      has_left = true;
      has_up = true;
      muT = TT.Mean();
      T_up_norm = T_up-muT;
      T_left_norm = T_left-muT;
    }
    else if (qNode.offset().x ==0&&qNode.offset().y>0) //only have left
    {
      T_left = qNode.leftBound.Crop(Point3i(qNode.overlap().height,0,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0));
      var_tag = T_left.Var()[0];
      muT = T_left.Mean();
      var_tag_left = var_tag;
      has_left = true;
      T_left_norm = T_left-muT;
    }
    else if (qNode.offset().x>0&&qNode.offset().y==0) // only have up 
    {
      T_up = qNode.upBound.Crop(Point3i(0,qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
      var_tag = T_up.Var()[0];
      var_tag_up = var_tag;
      has_up = true;
      muT = T_up.Mean();
      T_up_norm=T_up-muT;
    }
    else
    {
      CV_Error(CV_StsNotImplemented,"no boundary for matching"); 
      return sideMatchAddr;
    }
	  //cv::Mat gaussKernel = mylib::GenGaussKer(3,double(3)/6.0,CV_64F);
	  //Tensor<T,cn> T0_low=T0.Filter2D(gaussKernel,FILTER_BOUND_VALID);
	  //Tensor<T,cn> T1_low=T1.Filter2D(gaussKernel,FILTER_BOUND_VALID);
	  //Tensor<T,cn> T0_Lap=T0.Laplacian();
	  //Tensor<T,cn> T1_Lap=T1.Laplacian();
	  //Tensor<T,cn> filterBD,lowBD;
	for (int t=0; t< searchRegion.depth;t+=searchStep.depth)
		for (int x=offsetUp; x< offsetDown; x+=searchStep.height)
			for (int y=offsetLeft; y< offsetRight; y+=searchStep.width)
			  {
				  if (IsInsideCausalRegion(cv::Point3i(x,y,t),qNode,3))
				  {
					  diff=0;
					  N=0;	
		
            if (has_left && has_up)
            {
            	rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
						  rst.Ref(Cube(Point3i(x,y,t),qNode.leftBound.size()),B_left);
              VV.SetBlock(B_left.Transpose());
	            VV.SetBlock(Point3i(0,T_left.size().width,0),B_up);
              var_blk = VV.Var()[0];
              var_blk_up= B_up.Var()[0];
              var_blk_left = B_left.Var()[0];
              muB = VV.Mean();
              B_up_norm=B_up-muB;
              B_left_norm = B_left - muB;
            }
            else if (has_left)
            {
              rst.Ref(Cube(Point3i(x+qNode.overlap().height,y,t),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0)),B_left);
              var_blk = B_left.Var()[0];
              var_blk_left = var_blk;
              muB =B_left.Mean();
              B_left_norm = B_left - muB;
            }
            else if (has_up)
            {
             	rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
						  var_blk = B_up.Var()[0];
              var_blk_up = var_blk;
              muB = B_up.Mean();
              B_up_norm = B_up - muB;
            }
            //diff = abs(var_blk-var_tag);
            //diff = abs(log10(var_blk/var_tag));//change to this, Nov 10, 2012
            double diff_left=0,diff_up=0;
            if (has_left)
              diff_left = (log10((var_blk_left-var_tag_left)/var_tag_left));
            if (has_up)
              diff_up = (log10((var_blk_up-var_tag_up)/var_tag_up));
            if (diff_left < varThrd1  && diff_up < varThrd1 ) //varThrd1) //change thred, Nov 10, 2012
            {//gj01132013 change threahold from 0.5 to 0.2, gj01142013, change thred co
              //varqueue->compareInsert(diff,cv::Point3i(x,y,t));
              //varList.push_back(cv::Point3i(x,y,t));
              //round 2, search for nbs
                
		          for (int nbx=-searchStep.height/2; nbx<= searchStep.height/2; nbx+=1) //change from < to <=, Nov 10,2012
			          for (int nby=-searchStep.width/2; nby<= searchStep.width/2; nby+=1)
                {
                 
                  if (IsInsideCausalRegion(cv::Point3i(nbx+x,nby+y,t),qNode,3))
                  {
                    double mse_diff=DBL_MAX;
                    N=0;
                    if (has_left && has_up)
                    {
                      mse_diff = 0;
            	        rst.Ref(Cube(Point3i(x+nbx,y+nby+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
						          rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.leftBound.size()),B_left);

                      mse_diff+= metric::Compare(T_up_norm,B_up_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
                      mse_diff+= metric::Compare(T_left_norm,B_left_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
                      N+= T_up.size().volumn();
                      N+= T_left.size().volumn();
                    }
                    else if (has_left)
                    {
                      mse_diff = 0;
                      rst.Ref(Cube(Point3i(x+nbx+qNode.overlap().height,y+nby,t),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0)),B_left);
                      mse_diff+= metric::Compare(T_left_norm, B_left_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
                      N+= T_left.size().volumn();
                    }
                    else if (has_up)
                    {
                      mse_diff=0;
             	        rst.Ref(Cube(Point3i(x+nbx,y+nby+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
                      mse_diff+= metric::Compare(T_up_norm,B_up_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
                      N+= T_up.size().volumn();
                    }
                    mse_diff/=double(N);
                    mse_diff/=65025;//normalize 255^2
                    if (mse_diff<=matching_thrd)
                      queue->compareInsert(-mse_diff,Point3i(x+nbx,y+nby,t));
                    if (queue->getLength()>candidNum&&candidNum>0)
                      queue->pop();
                  }
                }
            }
          }
        }
	}
			if (queue->getLength()>0)
			{
        logfile<<"Tar: ("<<qNode.offset().x<<","<<qNode.offset().y<<")"<<endl;
        cout<<"Tar: ("<<qNode.offset().x<<","<<qNode.offset().y<<")"<<endl;
				for (int k=0; k<queue->getLength();k++)
				{
					if (queue->GetAddress()[k]!=Point3i(-1,-1,-1)) 
					{
            sideMatchAddr.push_back(queue->GetAddress()[k]);
            logfile<<"       <<<<<<<<  Cand: ("<<sideMatchAddr[k].x<<","<<sideMatchAddr[k].y<<") Dist: "<<queue->getData(k)<<endl;
            cout<<"       <<<<<<<<  Cand: ("<<sideMatchAddr[k].x<<","<<sideMatchAddr[k].y<<") Dist: "<<queue->getData(k)<<endl;
          }
				}
			}
      delete queue;
      logfile.close();
			return sideMatchAddr;
}
template<class T, size_t cn>
void QGrid<T,cn>::SetNode(const Point3i& gridPos, const Tensor<T,cn>& ts, const Cube& roi, const Size3& overlapSize)
{

	grid[gridPos.z][gridPos.x][gridPos.y].SetBlock(ts.Crop(roi.offset(),roi.size()));
}
template<class T, size_t cn>
Size3 QGrid<T,cn>::GetGridSize(void) const
{
	return gridSize;
}
template<class T, size_t cn>
void QGrid<T,cn>::SetVarThrd(double t1)
{
	this->varThrd1 = t1;
}

template<class T, size_t cn>
void QGrid<T,cn>::UpdateRstSqr(const Point3i& sPos, const Point3i& ePos)
{
  Cube temp = Cube(sPos,Size3(ePos.x-sPos.x+1,ePos.y-sPos.y+1,ePos.z-sPos.z+1));
        this->UpdateRstSqr(temp);
}


template<class T, size_t cn>
void QGrid<T,cn>::UpdateRstSqr(const Cube& cb)
{
	this->rstSqr.SetBlock(cb.offset(),(this->rst(cb)*this->rst(cb)));
}

/*
template<class T, size_t cn>
void QGrid<T,cn>::SetSubWinSize(const Size3& sz)
{
	this->ensemble.SetSubWinSize(sz);
	this->rst.SetSubWinSize(sz);
}

template<class T, size_t cn>
void QGrid<T,cn>::SetSubWinStep(const Size3& sz)
{
	this->ensemble.SetSubWinStep(sz);
	this->rst.SetSubWinStep(sz);
}


template<class T, size_t cn>
Size3 QGrid<T,cn>::GetSubWinSize(void) const
{
	return ensemble.GetSubWinSize();
}


template<class T, size_t cn>
Size3 QGrid<T,cn>::GetSubWinStep(void) const
{
	return ensemble.GetSubWinStep();
}
*/

template<class T, size_t cn>
void QGrid<T,cn>::ComputeBoundMean(const Tensor<T,cn>& tarUp,  const Tensor<T,cn>& tarLeft, const Point3i& offset, const Size3& size, const Size3& overlap, Tensor<T,cn>& upMean, Tensor<T,cn>& leftMean)
{
    leftMean = Tensor<T,cn>(tarLeft.size());
    upMean = Tensor<T,cn>(tarUp.size());
    Point3i up = offset + Point3i(-size.height/2,size.width/2,0);
    Point3i left = offset + Point3i(size.height/2,-size.width/2,0);
    Point3i anti = offset - overlap.Point3()*2;
    Size3 cropSize = Size3(size.height/2, size.width/2,1);
    auto u = rst.Crop(up,cropSize).Mean();

    auto l = rst.Crop(left,cropSize).Mean();
    auto a = rst.Crop(anti,cropSize).Mean();
    //interpolate upMean
    for (int ii=0; ii <upMean.size().width; ii++)
    {
      auto interpVal = ((upMean.size().width - ii-1)*u + (ii+1)*a)/upMean.size().width;
      for (int jj=0; jj<upMean.size().height; jj++)
        upMean(jj,ii,0) = interpVal;
    }
    for (int jj=overlap.height; jj<leftMean.size().height; jj++)
    {
      auto interpVal = ((leftMean.size().height-overlap.height-jj-1)*l + (jj-overlap.height+1)*a)/(leftMean.size().height-overlap.height);
      for (int ii=0; ii<leftMean.size().width; ii++)
        leftMean(jj,ii,0) = interpVal;
    }

    for (int jj=0; jj<overlap.height; jj++)
      for (int ii=0; ii<overlap.width;ii++)
      {
        leftMean(jj,ii,0) = (a + upMean(jj,0,0) + leftMean(overlap.height,ii,0))/3;
      }
}
template class QGrid<double,1>;
template class QGrid<double,2>;
template class QGrid<double,3>;
template class QGrid<uchar,1>;
template class QGrid<uchar,2>;
template class QGrid<uchar,3>;
template class QGrid<float,1>;
template class QGrid<float,2>;
template class QGrid<float,3>;


}

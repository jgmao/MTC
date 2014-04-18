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
    L1_train_len = 0;
    lenH0=0;
    lenH1=0;
    train1 = true;
    train2 = true;
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
  vector<pair<Point3i,double> > QGrid<T,cn>::BoundaryMatching(QNode<T,cn>& qNode, MatchingMethod matching_method, double matching_thrd, Size3 subWinSize)
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
    double ratio = 3;//20140128 from 5-->3
    int offsetUp = qNode.offset().x - int(ratio*double(qNode.size().height));
    offsetUp > 1 ? offsetUp = offsetUp: offsetUp=1;
    int offsetLeft = qNode.offset().y - int(ratio*double(qNode.size().width));
    offsetLeft>1? offsetLeft = offsetLeft: offsetLeft = 1;
    int offsetRight = qNode.offset().y + int(ratio*double(qNode.size().width));
    offsetRight>searchRegion.width-1?offsetRight = searchRegion.width-1:offsetRight = offsetRight;
    int offsetDown = qNode.offset().x + int(2*double(qNode.size().height));
    offsetDown>searchRegion.height-1?offsetDown = searchRegion.height-1:offsetDown  = offsetDown;
    LinkArray* queue;
    //if (candidNum>0)
    //  queue = new LinkQueue(candidNum);
    //LinkQueue queue(candidNum);//gj15012013  set possibility of dynamic number of candidate
    //else
    queue = new LinkArray();

    LinkQueue varQueue(100);
    vector<Point3i> varList;
    Tensor<T,cn> target;
    Tensor<T,cn> candUp, candLeft;
    QNode<T,cn> cand;
    Tensor<double,cn> tempMu,tempVar;
    Tensor<double,cn> tempMap;
    //vector<Point3i> sideMatchAddr;
    vector<pair<Point3i,double> > sideMatchAddr;
    Tensor<T,cn> org = ensemble.Crop(qNode.offset(),qNode.size());
    fstream logfile;
    logfile.open("./temp/matching.txt",ios::app);

    if (qNode.offset().x ==32 && qNode.offset().y==64)
      tempMap = Tensor<double,cn>(rst.size());
    //! 20130916 use opencv template matching
    Tensor<T,cn> tarUp , tarLeft;
    qNode.upBound.Ref(Cube(0,qNode.overlap().width,0,qNode.overlap().height,qNode.size().width,qNode.size().depth),tarUp);
    tarLeft = qNode.leftBound;
    // tarUp.Print();
    // tarLeft.Print();
    Tensor<T,cn> tarDevH, tarDevV;
    Vec<T,cn> tv;//total variance
    if (matching_method == MatchingMethod::MATCHING_OPENCV)
      {


        Mat matchUp, matchLeft,match;
        int cv_matchmethod =  TM_SQDIFF;
        Mat flUp,flLeft;
        Mat flRst;
        rst.convertTo(flRst,CV_32F);
        tarUp.convertTo(flUp,CV_32F);
        tarLeft.convertTo(flLeft,CV_32F);
        cv::matchTemplate(flRst,flUp,matchUp,cv_matchmethod);
        cv::matchTemplate(flRst,flLeft,matchLeft,cv_matchmethod);
        cout<<matchUp.size()<<endl;
        cout<<matchLeft.size()<<endl;
        //match = matchUp + matchLeft;
        normalize(matchUp, matchUp, 0, 255, NORM_MINMAX, -1, Mat() );
        normalize(matchLeft, matchLeft, 0, 255, NORM_MINMAX, -1, Mat() );
        Tensor<double,1>(matchUp).Display();
        Tensor<double,1>(matchLeft).Display();

      }
    else if (matching_method == MatchingMethod::MATCHING_MSE || matching_method ==  MatchingMethod::MATCHING_SAD|| matching_method == MatchingMethod::MATCHING_MSE_CONSTRAINT|| matching_method == MatchingMethod::MATCHING_DIRECT)
      {
        //old method, template matching include boundarys x 3

        if (qNode.offset().x == DEBUG_X && qNode.offset().y == DEBUG_Y && qNode.size().height==DEBUG_SIZE)
          {
            rst.SaveBlock("rst.tif");
            //logfile<<"\n ======== output debug info, template matching ====== \n";
            qNode.debugtrigger = true;
          }

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
#if SUBSTRACT_BOUND_MEAN
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
            //compute total variance here
            tarDevH = qNode.ComputeBoundHorDev();
            tarDevV = qNode.ComputeBoundVerDev();
            tv = tarDevH.Pow(2).Sum()+tarDevV.Pow(2).Sum();
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
                Tensor<T,cn> candDevH, candDevV;
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
#if SUBSTRACT_BOUND_MEAN
                                  ComputeBoundMean(localCandUp,localCandLeft,qNode.offset(),qNode.size(),qNode.overlap(),candUpMean,candLeftMean);
#else
                                  candUpMean = Tensor<T,cn>(localCandUp.size(),0);
                                  candLeftMean = Tensor<T,cn>(localCandLeft.size(),0);
#endif

                                }
#else
                                rst.Ref(Cube(Point3i(x,y,t),qNode.leftBound.size()),candLeft);
                                rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),candUp);

                                //! 20131008 use qNode instead of tensor to compute boundary total variance
                                //rst.Ref(Cube(Point3i(x+qNode.overlap().height,y+qNode.overlap().width,t),qNode.size()),cand);
                                cand = QNode<T,cn>(rst,qNode.size(),Point3i(x+qNode.overlap().height,y+qNode.overlap().width,t),qNode.overlap());
                                if (qNode.offset().x == DEBUG_X && qNode.offset().y == DEBUG_Y && qNode.size().height==DEBUG_SIZE)
                                  {
                                    if (x==29&&y==52)
                                      {
                                        candUp.Print("cand_up");
                                        candLeft.Print("cand_left");
                                      }
                                  }
#if SUBSTRACT_BOUND_MEAN
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
                                    candDevH = cand.ComputeBoundHorDev();
                                    candDevV = cand.ComputeBoundVerDev();
                                    tv += candDevH.Pow(2).Sum()+candDevV.Pow(2).Sum();
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
				    //! 20130913
				    //tarUp_norm.Print();
				    //candUp_norm.Print();
				    //tarLeft_norm.Print();
				    //candLeft_norm.Print();
				    localdiff = metric::Compare(tarUp_norm,candUp_norm,CompareCriteria::COMPARE_CRITERIA_MSE) +
					metric::Compare(tarLeft_norm, candLeft_norm, CompareCriteria::COMPARE_CRITERIA_MSE);
				    // cout<<"mse="<<localdiff<<endl;
				    //else if (matching_method == MATCHING_MSE_CONSTRAINT) //	diff+= qNdoe.GetBoundary(i).Compare(target,COMPARE_CRITERIA_MSE_CONSTRAINT);
				  }
				else if (matching_method == MatchingMethod::MATCHING_SAD)
				  {
				    localdiff+= metric::Compare(tarUp_norm,candUp_norm,CompareCriteria::COMPARE_CRITERIA_SAD)+
					metric::Compare(tarLeft_norm,candLeft_norm,CompareCriteria::COMPARE_CRITERIA_SAD);
				  }
				else if (matching_method == MatchingMethod::MATCHING_DIRECT)
				  {
				    localdiff= metric::Compare(org,cand,CompareCriteria::COMPARE_CRITERIA_SSIM,Size3(16,16,1),Size3(16,16,1),3,4,(int)FilterBoundary::FILTER_BOUND_FULL,(int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_NEW_L1,0,false);

				  }
				else
				  CV_Error(CV_StsUnsupportedFormat,"unsupport criteria\n");
				localN= candUp_norm.size().volumn()+candLeft_norm.size().volumn();

				//target.debugtrigger=false;
				//localdiff= localdiff/localN;
                //! 20131008 remove normilzation localdiff/=65025;//normalize by 255^2
                //! 20131008 use sqrt
                localdiff = sqrt(localdiff);
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
				    double delta_thrd = log(100);
				    delta > delta_thrd? delta = 1 : delta=delta/delta_thrd;
				    //d1 > delta_thrd? d1 = 1 : d1 = d1/delta_thrd;
				    //d0 > delta_thrd? d0 = 1 : d0 = d0/delta_thrd;
				    double lambda2 = 0.5;
				    double masking=0;
				    if (Mt[0]<50)
				      masking = 0.5;
				    else if (Mt[0]>200)
				      masking =0.5;
				    else
				      masking = 1;
				    //double constraint = lambda*(d0*d0 + d1*d1);//may be constraint by distribution similarity
                   //lambda2 = 1/sqrt( ((TT-VV)*(TT-VV)).Sum()[0]);
                    //!20131008 use other by guoxin double constraint = lambda*delta;// + masking*lambda2 ;
                   // cout<<"tv="<<tv[0]<<endl;
                   // cout<<"xyH="<<(2*(tarDevH*candDevH).Sum())[0]<<endl;
                   // cout<<"xyV="<<(2*(tarDevV*candDevV).Sum())[0]<<endl;
                   double constraint = sqrt((tv-2*(tarDevH*candDevH).Sum()-2*(tarDevV*candDevV).Sum())[0]);
                   constraint /= (2*sqrt(2)*256);
                    //!20131011 retray use variance
                    localdiff/=256;
                   // cout<<"constr="<<constraint<<endl;
           // cout<<"diff="<<localdiff<<endl;
                    localdiff = 1*localdiff +0.5*delta +0*constraint ;
          }
				/// debug here !!!!!!!!!!!!!
				//protect shared data here
#ifdef PARALLEL_ENABLED
				mymux.lock();
#endif
//cout<<matching_thrd<<endl;

				if (matching_method == MatchingMethod::MATCHING_DIRECT)
		  queue->compareInsert(localdiff,cv::Point3i(x,y,t));
				else{
                    if (localdiff<=matching_thrd||matching_thrd==0)//0 means accept all
                                      queue->compareInsert(-localdiff,cv::Point3i(x,y,t));
				  }
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
      else if (matching_method == MatchingMethod::MATCHING_STAT)
      {
        Tensor<T,cn> T_left,B_left;
        Tensor<T,cn> T_up, B_up;
        Tensor<T,cn> TT(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
#if !PARALLEL_MATCHING
        Tensor<T,cn> VV(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
        Tensor<T,cn> candLeftLocal,candUpLocal;
#endif
        double var_blk,var_tag;
        bool has_left=false, has_up=false;
        double diff_up=1000, diff_left=1000;
        double var_tar_left=-1,var_tar_up=-1, var_blk_left=-1, var_blk_up=-1;
        Size3 oSize(max(qNode.overlap().height,subWinSize.height),max(qNode.overlap().width,subWinSize.width),1);
        if (offsetUp + qNode.overlap().height-oSize.height<0)
          offsetUp += oSize.height-qNode.overlap().height;
        if (offsetLeft+qNode.overlap().width -oSize.width<0)
          offsetLeft += oSize.width-qNode.overlap().width;
        int shift_x = 0;
        int shift_y = 0;
        if (qNode.overlap().height<subWinSize.height)
          shift_x = qNode.overlap().height-oSize.height;
        if (qNode.overlap().width<subWinSize.width)
          shift_y = qNode.overlap().width-oSize.width;

        TT.SetBlock(tarUp.Transpose());
        TT.SetBlock(Point3i(tarUp.size().width,0,0),tarLeft);
        //Mt = TT.Mean();
        double Vt = TT.Var()[0];
        int t=0;
        Vec<T,cn> muT,muB;
        ofstream thrdfile("./thrdfile.txt",std::ios::out|std::ios::app);
        Tensor<T,cn> T_left_norm, T_up_norm, B_up_norm,B_left_norm;
//        if (qNode.offset().x > 0  && qNode.offset().y>0)
        {
            T_up =  qNode.upBound.Clone();// qNode.upBound.Crop(Point3i(0,qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
            T_left = 	 qNode.leftBound.Clone();
//            TT.SetBlock(T_left.Transpose());
//            TT.SetBlock(Point3i(0,T_left.size().width,0),T_up);
//            var_tag = TT.Var()[0];
            var_tar_left = T_left.Var()[0];
            var_tar_up = T_up.Var()[0];
            //has_left = true;
           // has_up = true;
           // muT = TT.Mean();
           // T_up_norm = T_up-muT;
           // T_left_norm = T_left-muT;
         }
//        else if (qNode.offset().x ==0&&qNode.offset().y>0) //only have left
//          {
//            T_left = qNode.leftBound.Crop(Point3i(qNode.overlap().height,0,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0));
//            var_tag = T_left.Var()[0];
//            muT = T_left.Mean();
//            var_tag_left = var_tag;
//            has_left = true;
//            T_left_norm = T_left-muT;
//          }
//        else if (qNode.offset().x>0&&qNode.offset().y==0) // only have up
//          {
//            T_up = qNode.upBound.Crop(Point3i(0,qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
//            var_tag = T_up.Var()[0];
//            var_tag_up = var_tag;
//            has_up = true;
//            muT = T_up.Mean();
//            T_up_norm=T_up-muT;
//          }
//        else
//          {
//            CV_Error(CV_StsNotImplemented,"no boundary for matching");
//            return sideMatchAddr;
//          }
        //cv::Mat gaussKernel = mylib::GenGaussKer(3,double(3)/6.0,CV_64F);
        //Tensor<T,cn> T0_low=T0.Filter2D(gaussKernel,FILTER_BOUND_VALID);
        //Tensor<T,cn> T1_low=T1.Filter2D(gaussKernel,FILTER_BOUND_VALID);
        //Tensor<T,cn> T0_Lap=T0.Laplacian();
        //Tensor<T,cn> T1_Lap=T1.Laplacian();
        //Tensor<T,cn> filterBD,lowBD;

#if PARALLEL_MATCHING
        int pnum = 4;
        vector<thread> threads;
#else
        Tensor<T,cn> tarSide(oSize);
        Tensor<T,cn> canSide(oSize);
        int pnum =1;// thread::hardware_concurrency();
#endif
        vector<LinkArray*> localqueue(pnum);
        int brows = (offsetDown - offsetUp)/searchStep.height/pnum;
        for (int p=0; p<pnum; p++)
        {
#if PARALLEL_MATCHING

          threads.push_back(thread([&](int p, int t){
          Tensor<T,cn> tarSide(oSize);
          Tensor<T,cn> canSide(oSize);
          Tensor<T,cn> VV(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
          Tensor<T,cn> candLeftLocal,candUpLocal;
          localqueue[p] = new LinkArray();
#else
          localqueue[p] = queue;
#endif
          for (int x=offsetUp+p*brows; x< min(offsetDown,offsetUp+(p+1)*brows); x+=searchStep.height)
            for (int y=offsetLeft; y< offsetRight; y+=searchStep.width)
              {
                if (IsInsideCausalRegion(cv::Point3i(x,y,0),qNode,3))
                  {
                    diff=0;
                    N=0;
                    //20140127 use varratio
                    rst.Ref(Cube(Point3i(x,y,t),qNode.leftBound.size()),B_left);//candLeftLocal);
                    //rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),candUpLocal);
                    rst.Ref(Cube(Point3i(x,y,t),qNode.upBound.size()),B_up);//candUpLocal);

                    //VV.SetBlock(candUpLocal.Transpose());
                    //VV.SetBlock(Point3i(tarUp.size().width,0,0),candLeftLocal);
                    //get local variance
                    //Tensor<T,cn> tempblk;
                    //rst.Ref(Cube(Point3i(x,y,t),qNode.size()+qNode.overlap()),tempblk);
                    //double Vb = VV.Var()[0];
                    //double Vb = tempblk.Var()[0];//local variance is close to target side variance
                    //double delta = abs(1-Vb/Vt);
                    //clip
                   // double delta_thrd = 20;//10*log10(100);
                    //delta > delta_thrd? delta = 1 : delta=delta/delta_thrd;
                    //double delta = abs(Vb - Vt);
                    //cout<<"delta is "<<delta<<endl;
//                    if (has_left && has_up)
                      {
//                        rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                        rst.Ref(Cube(Point3i(x,y,0),qNode.leftBound.size()),B_left);
//                        VV.SetBlock(B_left.Transpose());
//                        VV.SetBlock(Point3i(0,T_left.size().width,0),B_up);
//                        var_blk = VV.Var()[0];
                          var_blk_up= B_up.Var()[0];
                          var_blk_left = B_left.Var()[0];
                          //muB_up = B_up.Mean();
                         // muB_left = B_left.Mean();
//                        muB = VV.Mean();
                         // B_up_norm=B_up-muB_up;
                         // B_left_norm = B_left - muB_left;
                          //diff_up = abs(var_blk_up - var_tar_up);
                          //diff_left = abs(var_blk_left-var_tar_left);
                          diff_up = (var_blk_up < var_tar_up)?abs(1-var_blk_up/var_tar_up):abs(1-var_tar_up/var_blk_up);
                          diff_left = (var_blk_left<var_tar_left)? abs(1-var_blk_left/var_tar_left):abs(1-var_tar_left/var_blk_left);
                          //for the small variance, I want to make the diff less sensitive 20140130
                          double adj_left = log10(var_tar_left);
                          double adj_up = log10(var_tar_up);
                          adj_left *= adj_left;
                          adj_left = (adj_left>1)?1:adj_left;
                          adj_up *= adj_up;
                          adj_up = (adj_up>1)?1:adj_up;
                          diff_up*=adj_up;
                          diff_left*=adj_left;
                      }
//                    else if (has_left)
//                      {
//                        rst.Ref(Cube(Point3i(x+qNode.overlap().height,y,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0)),B_left);
//                        var_blk = B_left.Var()[0];
//                        var_blk_left = var_blk;
//                        muB =B_left.Mean();
//                        B_left_norm = B_left - muB;
//                      }
//                    else if (has_up)
//                      {
//                        rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                        var_blk = B_up.Var()[0];
//                        var_blk_up = var_blk;
//                        muB = B_up.Mean();
//                        B_up_norm = B_up - muB;
//                      }
//                    //diff = abs(var_blk-var_tag);
//                    //diff = abs(log10(var_blk/var_tag));//change to this, Nov 10, 2012
//                    double diff_left=0,diff_up=0;
//                    if (has_left)
//                    {
//                        if (abs(var_blk_left-var_tag_left)<0.1)
//                            diff_left =0;
//                        else
//                            diff_left = abs(log10((var_blk_left-var_tag_left)/var_tag_left));
//                    }
//                    if (has_up)
//                    {
//                        if (abs(var_blk_up-var_tag_up)<0.1)
//                            diff_left = 0;
//                        else
//                            diff_up = abs(log10((var_blk_up-var_tag_up)/var_tag_up));
//                    }
                    //cout<<diff_left<<","<<diff_up<<endl;
                    //if (delta<=varThrd1) //20140127 first trim (coarest) discard large var ratio
                    if (true) //varThrd1) //change thred, Nov 10, 2012
                    {//gj01132013 change threahold from 0.5 to 0.2, gj01142013, change thred co
                        //varqueue->compareInsert(diff,cv::Point3i(x,y,t));
                        //varList.push_back(cv::Point3i(x,y,t));
                        //round 2, search for nbs

                        int nbx = 0;
                        int nby = 0;
                        {
                              if (x+nbx<offsetUp||x+nbx>offsetDown||y+nby<offsetLeft||y+nby>offsetRight) //20131227, make sure searching in neighborhood do not voilate search range constrain
                                continue;
                              if (x+nbx<offsetUp+p*brows||x+nbx>=offsetUp+(p+1)*brows) //no duplicate search in different threads
                                continue;
                              if (IsInsideCausalRegion(cv::Point3i(nbx+x,nby+y,0),qNode,3))
                                {
                                  thrdfile<<qNode.offset().x<<","<<qNode.offset().y<<","<<x<<","<<y<<","<<var_tar_left<<","<<var_tar_up<<","<<var_blk_left<<","<<var_blk_up<<","<< diff_left <<","<<diff_up<<",";

                                  rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.leftBound.size()),B_left);//candLeftLocal);
                                  rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.upBound.size()),B_up);//candUpLocal);
                                  double mse_up = metric::ComputeMSE(T_up,B_up);
                                  double mse_left = metric::ComputeMSE(T_left,B_left);//metric::ComputeMSE(VV,TT);
                                  //mse /= 65025; //normalize
                                  thrdfile<<mse_left<<","<<mse_up<<",";
                                  N=0;
                                  double score=0;
                                  double final_score=1;
                                  int count=0;
                                  for (int m = 0; m<qNode.leftBound.size().height-shift_x; m+=oSize.height)
                                  {
                                    ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(m,0,1),oSize),tarSide);
                                    rst.Ref(Cube(Point3i(x+nbx+m+shift_x,y+nby+shift_y,t),oSize),canSide);
                                    score=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize,oSize, 3, 1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                                    if (score<final_score)
                                      final_score= score;
                                    count++;
                                  }
                                  for (int n = oSize.width; n<qNode.upBound.size().width-shift_y; n+=oSize.width)
                                  {
                                      ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(0,n,1),oSize),tarSide);
                                      rst.Ref(Cube(Point3i(x+nbx+shift_x,y+nby+n+shift_y,t),oSize),canSide);
                                      score=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize, oSize,3, 1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                                      if (score<final_score)
                                        final_score= score;
                                      count++;
                                  }
                                  //score = score/count;//average;
                                  //got the min of all scores
                                  thrdfile<< score<<",";
                                  Tensor<T,cn> cand;
                                  Tensor<T,cn> org;
                                  ensemble.Ref(Cube(qNode.offset(),qNode.size()),org);
                                  rst.Ref(Cube(Point3i(x+nbx+qNode.overlap().height, y+nby+qNode.overlap().width,0),qNode.size()),cand);  
                                  double temp = metric::Compare(org,cand,CompareCriteria::COMPARE_CRITERIA_SSIM,Size3(16,16,1), Size3(16,16,1),3,4,(int)FilterBoundary::FILTER_BOUND_EXTEND/*true*/,(int)FeaturePoolType::FEATURE_POOL_MIN,(int)MetricModifier::STSIM2_BASELINE,0,false);
                                  thrdfile<<temp<<endl;                  
                                  localqueue[p]->compareInsert(final_score,cv::Point3i(x+nbx,y+nby,t));
                                  if (localqueue[p]->getLength()>candidNum/pnum&&candidNum>0)
                                    localqueue[p]->pop();
//                                  if (has_left && has_up)
//                                    {
//                                      mse_diff = 0;
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby,0),qNode.leftBound.size()),B_left);
//#if DEBUG_20131226
//                                      B_up.Print();
//                                      B_left.Print();
//#endif
//                                      double diff_up =   metric::Compare(T_up,B_up,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      double diff_left = metric::Compare(T_left,B_left,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      //Compare normalize the value, it's not correct
//                                      //so fix it
//                                      diff_up*=(double)T_up.size().area();
//                                      diff_left*=(double)T_left.size().area();
//                                      mse_diff = diff_up+diff_left;
//                                      N+= T_up.size().volumn();
//                                      N+= T_left.size().volumn();
//                                    }
//                                  else if (has_left)
//                                    {
//                                      mse_diff = 0;
//                                      rst.Ref(Cube(Point3i(x+nbx+qNode.overlap().height,y+nby,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0)),B_left);
//                                      mse_diff+= metric::Compare(T_left_norm, B_left_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      N+= T_left.size().volumn();
//                                    }
//                                  else if (has_up)
//                                    {
//                                      mse_diff=0;
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                                      mse_diff+= metric::Compare(T_up_norm,B_up_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      N+= T_up.size().volumn();
//                                    }
//                                  mse_diff/=double(N);
//#if DEBUG_20131226
//                                  cout<<x+nbx<<","<<y+nby<<endl;
//                                  cout<<mse_diff<<endl;
//                                  T_up.Print();
//                                  T_left.Print();
//#endif
//                                  mse_diff/=65025;//normalize 255^2
//#if DEBUG_20131226
//                                  cout<<mse_diff<<endl;
//#endif
//                                  if (mse_diff<=matching_thrd)
//                                    queue->compareInsert(-mse_diff,Point3i(x+nbx,y+nby,0));
//                                  if (queue->getLength()>candidNum&&candidNum>0)
//                                    queue->pop();
                                }
                            }
                      }
                  }
              }
#if PARALLEL_MATCHING
       },p,t));
#endif
     }
#if PARALLEL_MATCHING
     for (int p=0; p<pnum; p++)
     {
       auto& thread = threads[p];
       thread.join();
       for (int ii=0; ii<localqueue[p]->getLength(); ii++)
       {
           queue->compareInsert(localqueue[p]->getData(ii),localqueue[p]->GetAddress(ii));
       }
       delete localqueue[p];
     }
#endif
     #if     OUTPUT_THRDFILE
        thrdfile.close();
#endif
   }

    else if (matching_method == MatchingMethod::MATCHING_HIERARCHY)
      {
        Tensor<T,cn> T_left,B_left;
        Tensor<T,cn> T_up, B_up;
        Tensor<T,cn> TT(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
#if !PARALLEL_MATCHING
        Tensor<T,cn> VV(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
        Tensor<T,cn> candLeftLocal,candUpLocal;
#endif
        double var_blk,var_tag;
        bool has_left=false, has_up=false;
        double diff_up=1000, diff_left=1000;
        double var_tar_left=-1,var_tar_up=-1, var_blk_left=-1, var_blk_up=-1;
        Size3 oSize(max(qNode.overlap().height,subWinSize.height),max(qNode.overlap().width,subWinSize.width),1);
        if (offsetUp + qNode.overlap().height-oSize.height<0)
          offsetUp += oSize.height-qNode.overlap().height;
        if (offsetLeft+qNode.overlap().width -oSize.width<0)
          offsetLeft += oSize.width-qNode.overlap().width;
        int shift_x = 0;
        int shift_y = 0;
        if (qNode.overlap().height<subWinSize.height)
          shift_x = qNode.overlap().height-oSize.height;
        if (qNode.overlap().width<subWinSize.width)
          shift_y = qNode.overlap().width-oSize.width;

        TT.SetBlock(tarUp.Transpose());
        TT.SetBlock(Point3i(tarUp.size().width,0,0),tarLeft);
        //Mt = TT.Mean();
        double Vt = TT.Var()[0];
        int t=0;
        Vec<T,cn> muT,muB;
#if     OUTPUT_THRDFILE
        ofstream thrdfile("./thrdfile.txt",std::ios::out|std::ios::app);
#endif
        Tensor<T,cn> T_left_norm, T_up_norm, B_up_norm,B_left_norm;
//        if (qNode.offset().x > 0  && qNode.offset().y>0)
        {
            T_up =  qNode.upBound.Clone();// qNode.upBound.Crop(Point3i(0,qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
            T_left = 	 qNode.leftBound.Clone();
//            TT.SetBlock(T_left.Transpose());
//            TT.SetBlock(Point3i(0,T_left.size().width,0),T_up);
//            var_tag = TT.Var()[0];
            var_tar_left = T_left.Var()[0];
            var_tar_up = T_up.Var()[0];
            //has_left = true;
           // has_up = true;
           // muT = TT.Mean();
           // T_up_norm = T_up-muT;
           // T_left_norm = T_left-muT;
#if     OUTPUT_THRDFILE
           thrdfile << "org: "<<qNode.offset().x<<","<<qNode.offset().y<<"var: left "<<var_tar_left<<", up "<<var_tar_up<<"----------------------------"<<endl;
#endif
         }
//        else if (qNode.offset().x ==0&&qNode.offset().y>0) //only have left
//          {
//            T_left = qNode.leftBound.Crop(Point3i(qNode.overlap().height,0,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0));
//            var_tag = T_left.Var()[0];
//            muT = T_left.Mean();
//            var_tag_left = var_tag;
//            has_left = true;
//            T_left_norm = T_left-muT;
//          }
//        else if (qNode.offset().x>0&&qNode.offset().y==0) // only have up
//          {
//            T_up = qNode.upBound.Crop(Point3i(0,qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
//            var_tag = T_up.Var()[0];
//            var_tag_up = var_tag;
//            has_up = true;
//            muT = T_up.Mean();
//            T_up_norm=T_up-muT;
//          }
//        else
//          {
//            CV_Error(CV_StsNotImplemented,"no boundary for matching");
//            return sideMatchAddr;
//          }
        //cv::Mat gaussKernel = mylib::GenGaussKer(3,double(3)/6.0,CV_64F);
        //Tensor<T,cn> T0_low=T0.Filter2D(gaussKernel,FILTER_BOUND_VALID);
        //Tensor<T,cn> T1_low=T1.Filter2D(gaussKernel,FILTER_BOUND_VALID);
        //Tensor<T,cn> T0_Lap=T0.Laplacian();
        //Tensor<T,cn> T1_Lap=T1.Laplacian();
        //Tensor<T,cn> filterBD,lowBD;

#if PARALLEL_MATCHING
        int pnum = 4;
        vector<thread> threads;
#else
        Tensor<T,cn> tarSide(oSize);
        Tensor<T,cn> canSide(oSize);
        int pnum =1;// thread::hardware_concurrency();
#endif
        vector<LinkArray*> localqueue(pnum);
        int brows = (offsetDown - offsetUp)/searchStep.height/pnum;
        for (int p=0; p<pnum; p++)
        {
#if PARALLEL_MATCHING

          threads.push_back(thread([&](int p, int t){
          Tensor<T,cn> tarSide(oSize);
          Tensor<T,cn> canSide(oSize);
          Tensor<T,cn> VV(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
          Tensor<T,cn> candLeftLocal,candUpLocal;
          localqueue[p] = new LinkArray();
#else
          localqueue[p] = queue;
#endif
          for (int x=offsetUp+p*brows; x< min(offsetDown,offsetUp+(p+1)*brows); x+=searchStep.height)
            for (int y=offsetLeft; y< offsetRight; y+=searchStep.width)
              {
                if (IsInsideCausalRegion(cv::Point3i(x,y,0),qNode,3))
                  {
                    diff=0;
                    N=0;
                    //20140127 use varratio
                    rst.Ref(Cube(Point3i(x,y,t),qNode.leftBound.size()),B_left);//candLeftLocal);
                    //rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),candUpLocal);
                    rst.Ref(Cube(Point3i(x,y,t),qNode.upBound.size()),B_up);//candUpLocal);

                    //VV.SetBlock(candUpLocal.Transpose());
                    //VV.SetBlock(Point3i(tarUp.size().width,0,0),candLeftLocal);
                    //get local variance
                    //Tensor<T,cn> tempblk;
                    //rst.Ref(Cube(Point3i(x,y,t),qNode.size()+qNode.overlap()),tempblk);
                    //double Vb = VV.Var()[0];
                    //double Vb = tempblk.Var()[0];//local variance is close to target side variance
                    //double delta = abs(1-Vb/Vt);
                    //clip
                   // double delta_thrd = 20;//10*log10(100);
                    //delta > delta_thrd? delta = 1 : delta=delta/delta_thrd;
                    //double delta = abs(Vb - Vt);
                    //cout<<"delta is "<<delta<<endl;
//                    if (has_left && has_up)
                      {
//                        rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                        rst.Ref(Cube(Point3i(x,y,0),qNode.leftBound.size()),B_left);
//                        VV.SetBlock(B_left.Transpose());
//                        VV.SetBlock(Point3i(0,T_left.size().width,0),B_up);
//                        var_blk = VV.Var()[0];
                          var_blk_up= B_up.Var()[0];
                          var_blk_left = B_left.Var()[0];
                          //muB_up = B_up.Mean();
                         // muB_left = B_left.Mean();
//                        muB = VV.Mean();
                         // B_up_norm=B_up-muB_up;
                         // B_left_norm = B_left - muB_left;
                          //diff_up = abs(var_blk_up - var_tar_up);
                          //diff_left = abs(var_blk_left-var_tar_left);
                          diff_up = (var_blk_up < var_tar_up)?abs(1-var_blk_up/var_tar_up):abs(1-var_tar_up/var_blk_up);
                          diff_left = (var_blk_left<var_tar_left)? abs(1-var_blk_left/var_tar_left):abs(1-var_tar_left/var_blk_left);
                          //for the small variance, I want to make the diff less sensitive 20140130
                          double adj_left = log10(var_tar_left);
                          double adj_up = log10(var_tar_up);
                          adj_left *= adj_left;
                          adj_left = (adj_left>1)?1:adj_left;
                          adj_up *= adj_up;
                          adj_up = (adj_up>1)?1:adj_up;
                          //20140205 tempdisable
                          //diff_up*=adj_up;
                          //diff_left*=adj_left;
#if     OUTPUT_THRDFILE
                          thrdfile<<"-- can: "<<x<<", "<<y<<"var: "<<var_blk_left<<", "<<var_blk_up<<" -- layer1: "<< diff_left <<", "<<diff_up<<endl;
#endif
                      }
//                    else if (has_left)
//                      {
//                        rst.Ref(Cube(Point3i(x+qNode.overlap().height,y,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0)),B_left);
//                        var_blk = B_left.Var()[0];
//                        var_blk_left = var_blk;
//                        muB =B_left.Mean();
//                        B_left_norm = B_left - muB;
//                      }
//                    else if (has_up)
//                      {
//                        rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                        var_blk = B_up.Var()[0];
//                        var_blk_up = var_blk;
//                        muB = B_up.Mean();
//                        B_up_norm = B_up - muB;
//                      }
//                    //diff = abs(var_blk-var_tag);
//                    //diff = abs(log10(var_blk/var_tag));//change to this, Nov 10, 2012
//                    double diff_left=0,diff_up=0;
//                    if (has_left)
//                    {
//                        if (abs(var_blk_left-var_tag_left)<0.1)
//                            diff_left =0;
//                        else
//                            diff_left = abs(log10((var_blk_left-var_tag_left)/var_tag_left));
//                    }
//                    if (has_up)
//                    {
//                        if (abs(var_blk_up-var_tag_up)<0.1)
//                            diff_left = 0;
//                        else
//                            diff_up = abs(log10((var_blk_up-var_tag_up)/var_tag_up));
//                    }
                    //cout<<diff_left<<","<<diff_up<<endl;
                    //if (delta<=varThrd1) //20140127 first trim (coarest) discard large var ratio
                    if (diff_left < varThrd1  && diff_up < varThrd1 ) //varThrd1) //change thred, Nov 10, 2012
                    {//gj01132013 change threahold from 0.5 to 0.2, gj01142013, change thred co
                        //varqueue->compareInsert(diff,cv::Point3i(x,y,t));
                        //varList.push_back(cv::Point3i(x,y,t));
                        //round 2, search for nbs


                        for (int nbx=-searchStep.height/2; nbx< searchStep.height/2; nbx+=1) //change from < to <=, Nov 10,2012, wrong ! 20140130
                          for (int nby=-searchStep.width/2; nby< searchStep.width/2; nby+=1)
                            {
                              if (x+nbx<offsetUp||x+nbx>offsetDown||y+nby<offsetLeft||y+nby>offsetRight) //20131227, make sure searching in neighborhood do not voilate search range constrain
                                continue;
                              if (x+nbx<offsetUp+p*brows||x+nbx>=offsetUp+(p+1)*brows) //no duplicate search in different threads
                                continue;
                              if (IsInsideCausalRegion(cv::Point3i(nbx+x,nby+y,0),qNode,3))
                                {
                                  rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.leftBound.size()),B_left);//candLeftLocal);
                                  rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.upBound.size()),B_up);//candUpLocal);
                                  double mse_up = metric::ComputeMSE(T_up,B_up);
                                  double mse_left = metric::ComputeMSE(T_left,B_left);//metric::ComputeMSE(VV,TT);
                                  //mse /= 65025; //normalize
                                  #if     OUTPUT_THRDFILE
                                  thrdfile<<"---layer 2:"<<mse_left<<", "<<mse_up<<endl;
                                   #endif
                                  if (mse_left>matching_thrd&&mse_up>matching_thrd) //2nd trim: discard large side mse candidate
                                    continue;
                                  N=0;
                                  double score=0;
                                  double final_score=1;
                                  int count=0;
                                  for (int m = 0; m<qNode.leftBound.size().height-shift_x; m+=oSize.height)
                                  {
                                    ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(m,0,1),oSize),tarSide);
                                    rst.Ref(Cube(Point3i(x+nbx+m+shift_x,y+nby+shift_y,t),oSize),canSide);
                                    score=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize,oSize, 3, 1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                                    if (score<final_score)
                                      final_score= score;
                                    count++;
                                  }
                                  for (int n = oSize.width; n<qNode.upBound.size().width-shift_y; n+=oSize.width)
                                  {
                                      ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(0,n,1),oSize),tarSide);
                                      rst.Ref(Cube(Point3i(x+nbx+shift_x,y+nby+n+shift_y,t),oSize),canSide);
                                      score=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize, oSize,3, 1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                                      if (score<final_score)
                                        final_score= score;
                                      count++;
                                  }
                                  //score = score/count;//average;
                                  //got the min of all scores
                                  #if     OUTPUT_THRDFILE
                                  thrdfile<<"---layer 3:"<<x+nbx<<", "<<y+nby<<" score: "<< score<<endl;
                                  if (qNode.offset().x == DEBUG_X && qNode.offset().y==DEBUG_Y&&x+nbx==18 && y+nby==292)
                                    cout<<"x,y "<<x<<", "<<y<<endl;
                                  #endif
                                  localqueue[p]->compareInsert(final_score,cv::Point3i(x+nbx,y+nby,t));
                                  if (localqueue[p]->getLength()>candidNum/pnum&&candidNum>0)
                                    localqueue[p]->pop();
//                                  if (has_left && has_up)
//                                    {
//                                      mse_diff = 0;
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby,0),qNode.leftBound.size()),B_left);
//#if DEBUG_20131226
//                                      B_up.Print();
//                                      B_left.Print();
//#endif
//                                      double diff_up =   metric::Compare(T_up,B_up,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      double diff_left = metric::Compare(T_left,B_left,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      //Compare normalize the value, it's not correct
//                                      //so fix it
//                                      diff_up*=(double)T_up.size().area();
//                                      diff_left*=(double)T_left.size().area();
//                                      mse_diff = diff_up+diff_left;
//                                      N+= T_up.size().volumn();
//                                      N+= T_left.size().volumn();
//                                    }
//                                  else if (has_left)
//                                    {
//                                      mse_diff = 0;
//                                      rst.Ref(Cube(Point3i(x+nbx+qNode.overlap().height,y+nby,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0)),B_left);
//                                      mse_diff+= metric::Compare(T_left_norm, B_left_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      N+= T_left.size().volumn();
//                                    }
//                                  else if (has_up)
//                                    {
//                                      mse_diff=0;
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                                      mse_diff+= metric::Compare(T_up_norm,B_up_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      N+= T_up.size().volumn();
//                                    }
//                                  mse_diff/=double(N);
//#if DEBUG_20131226
//                                  cout<<x+nbx<<","<<y+nby<<endl;
//                                  cout<<mse_diff<<endl;
//                                  T_up.Print();
//                                  T_left.Print();
//#endif
//                                  mse_diff/=65025;//normalize 255^2
//#if DEBUG_20131226
//                                  cout<<mse_diff<<endl;
//#endif
//                                  if (mse_diff<=matching_thrd)
//                                    queue->compareInsert(-mse_diff,Point3i(x+nbx,y+nby,0));
//                                  if (queue->getLength()>candidNum&&candidNum>0)
//                                    queue->pop();
                                }
                            }
                      }
                  }
              }
#if PARALLEL_MATCHING
       },p,t));
#endif
     }
#if PARALLEL_MATCHING
     for (int p=0; p<pnum; p++)
     {
       auto& thread = threads[p];
       thread.join();
       for (int ii=0; ii<localqueue[p]->getLength(); ii++)
       {
           queue->compareInsert(localqueue[p]->getData(ii),localqueue[p]->GetAddress(ii));
       }
       delete localqueue[p];
     }
#endif
     #if     OUTPUT_THRDFILE
        thrdfile.close();
#endif
   }
  else if (matching_method == MatchingMethod::MATCHING_HIERARCHY2)
      {
  //matching stsim first then do constraint matching
        Tensor<T,cn> T_left,B_left;
        Tensor<T,cn> T_up, B_up;
        Tensor<T,cn> TT(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
#if !PARALLEL_MATCHING
        Tensor<T,cn> VV(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
        Tensor<T,cn> candLeftLocal,candUpLocal;
#endif
        double var_blk,var_tag;
        bool has_left=false, has_up=false;
        double diff_up=1000, diff_left=1000;
        double var_tar_left=-1,var_tar_up=-1, var_blk_left=-1, var_blk_up=-1;
        Size3 oSize(max(qNode.overlap().height,subWinSize.height),max(qNode.overlap().width,subWinSize.width),1);
        if (offsetUp + qNode.overlap().height-oSize.height<0)
          offsetUp += oSize.height-qNode.overlap().height;
        if (offsetLeft+qNode.overlap().width -oSize.width<0)
          offsetLeft += oSize.width-qNode.overlap().width;
        int shift_x = 0;
        int shift_y = 0;
        if (qNode.overlap().height<subWinSize.height)
          shift_x = qNode.overlap().height-oSize.height;
        if (qNode.overlap().width<subWinSize.width)
          shift_y = qNode.overlap().width-oSize.width;

        TT.SetBlock(tarUp.Transpose());
        TT.SetBlock(Point3i(tarUp.size().width,0,0),tarLeft);
        //Mt = TT.Mean();
        double Vt = TT.Var()[0];
        int t=0;
        Vec<T,cn> muT,muB;
#if     OUTPUT_THRDFILE
        ofstream thrdfile("./thrdfile.txt",std::ios::out|std::ios::app);
#endif
        Tensor<T,cn> T_left_norm, T_up_norm, B_up_norm,B_left_norm;
        {
            T_up =  qNode.upBound.Clone();
            T_left = 	 qNode.leftBound.Clone();
            var_tar_left = T_left.Var()[0];
            var_tar_up = T_up.Var()[0];
#if     OUTPUT_THRDFILE
           thrdfile << "org: "<<qNode.offset().x<<","<<qNode.offset().y<<"var: left "<<var_tar_left<<", up "<<var_tar_up<<"----------------------------"<<endl;
#endif
         }
#if PARALLEL_MATCHING
        int pnum = 4;
        vector<thread> threads;
#else
        Tensor<T,cn> tarSide(oSize);
        Tensor<T,cn> canSide(oSize);
        int pnum =1;// thread::hardware_concurrency();
#endif
        vector<LinkArray*> localqueue(pnum);
        int brows = (offsetDown - offsetUp)/searchStep.height/pnum;
        for (int p=0; p<pnum; p++)
        {
#if PARALLEL_MATCHING

          threads.push_back(thread([&](int p, int t){
          Tensor<T,cn> tarSide(oSize);
          Tensor<T,cn> canSide(oSize);
          Tensor<T,cn> VV(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
          Tensor<T,cn> candLeftLocal,candUpLocal;
          localqueue[p] = new LinkArray();
#else
          localqueue[p] = queue;
#endif
          for (int x=offsetUp+p*brows; x< min(offsetDown,offsetUp+(p+1)*brows); x+=searchStep.height)
            for (int y=offsetLeft; y< offsetRight; y+=searchStep.width)
              {
                if (IsInsideCausalRegion(cv::Point3i(x,y,0),qNode,3))
                  {
                    diff=0;
                    N=0;
                    double score;
                    double final_score=1;
                    int count = 0;
                    for (int m = 0; m<qNode.leftBound.size().height-shift_x; m+=oSize.height)
                    {
                          ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(m,0,1),oSize),tarSide);
                          rst.Ref(Cube(Point3i(x+m+shift_x,y+shift_y,t),oSize),canSide);
                          score=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize,oSize, 3, 1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                          if (score<final_score)
                            final_score= score;
                          count++;
                    }
                    for (int n = oSize.width; n<qNode.upBound.size().width-shift_y; n+=oSize.width)
                    {
                            ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(0,n,1),oSize),tarSide);
                            rst.Ref(Cube(Point3i(x+shift_x,y+n+shift_y,t),oSize),canSide);
                            score=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize, oSize,3, 1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                            if (score<final_score)
                              final_score= score;
                            count++;
                    }
                    if (final_score>matching_thrd) 
                    {
                        for (int nbx=-searchStep.height/2; nbx< searchStep.height/2; nbx+=1) //change from < to <=, Nov 10,2012, wrong ! 20140130
                          for (int nby=-searchStep.width/2; nby< searchStep.width/2; nby+=1)
                            {
                              if (x+nbx<offsetUp||x+nbx>offsetDown||y+nby<offsetLeft||y+nby>offsetRight) //20131227, make sure searching in neighborhood do not voilate search range constrain
                                continue;
                              if (x+nbx<offsetUp+p*brows||x+nbx>=offsetUp+(p+1)*brows) //no duplicate search in different threads
                                continue;
                              if (IsInsideCausalRegion(cv::Point3i(nbx+x,nby+y,0),qNode,3))
                                {
                                  rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.leftBound.size()),B_left);//candLeftLocal);
                                  rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.upBound.size()),B_up);//candUpLocal);
                                  double mse_up = metric::ComputeMSE(T_up,B_up);
                                  double mse_left = metric::ComputeMSE(T_left,B_left);//metric::ComputeMSE(VV,TT);
                                  var_blk_left = B_left.Var()[0];
                                  var_blk_up = B_up.Var()[0]; 
                                  diff_up = abs(var_blk_up - var_tar_up);
                                  diff_left = abs(var_blk_left<var_tar_left);
                                   
                                  diff = (mse_up + 0.5*diff_up) + (mse_left + 0.5*diff_left);
                                  localqueue[p]->compareInsert(diff,cv::Point3i(x+nbx,y+nby,t));
                                  if (localqueue[p]->getLength()>candidNum/pnum&&candidNum>0)
                                    localqueue[p]->pop();
                                }
                            }
                      }
                  }
              }
#if PARALLEL_MATCHING
       },p,t));
#endif
     }
#if PARALLEL_MATCHING
     for (int p=0; p<pnum; p++)
     {
       auto& thread = threads[p];
       thread.join();
       for (int ii=0; ii<localqueue[p]->getLength(); ii++)
       {
           queue->compareInsert(localqueue[p]->getData(ii),localqueue[p]->GetAddress(ii));
       }
       delete localqueue[p];
     }
#endif
     #if     OUTPUT_THRDFILE
        thrdfile.close();
#endif
   }
else if (matching_method == MatchingMethod::MATCHING_HIERARCHY3)
      {
        Tensor<T,cn> T_left,B_left;
        Tensor<T,cn> T_up, B_up;
        Tensor<T,cn> TT(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
#if !PARALLEL_MATCHING
        Tensor<T,cn> VV(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
        Tensor<T,cn> candLeftLocal,candUpLocal;
#endif
        bool dec1,dec2;
        double lambda_l1=-1;
        double lambda_l2=-1;
        double lambda_u1 = -1;
        double lambda_u2=-1;
        double eta1=100,eta2=100;
        double likelihood1,likelihood2;
        double likelihood1_l,likelihood1_u;
        double likelihood2_l,likelihood2_u;
        double var_blk,var_tag;
        bool has_left=false, has_up=false;
        double diff_up=1000, diff_left=1000;
        double Ph0, Ph1;
        double ratio_up=1, ratio_left=1;
        double var_tar_left=-1,var_tar_up=-1, var_blk_left=-1, var_blk_up=-1;
        Size3 oSize(max(qNode.overlap().height,subWinSize.height),max(qNode.overlap().width,subWinSize.width),1);
        if (offsetUp + qNode.overlap().height-oSize.height<0)
          offsetUp += oSize.height-qNode.overlap().height;
        if (offsetLeft+qNode.overlap().width -oSize.width<0)
          offsetLeft += oSize.width-qNode.overlap().width;
        int shift_x = 0;
        int shift_y = 0;
        if (qNode.overlap().height<subWinSize.height)
          shift_x = qNode.overlap().height-oSize.height;
        if (qNode.overlap().width<subWinSize.width)
          shift_y = qNode.overlap().width-oSize.width;
        Mat var_ratio(1,2,CV_64F);
        Mat logmse(1,2,CV_64F);
        TT.SetBlock(tarUp.Transpose());
        TT.SetBlock(Point3i(tarUp.size().width,0,0),tarLeft);
        //Mt = TT.Mean();
        double Vt = TT.Var()[0];
        int t=0;
        Vec<T,cn> muT,muB;
        Tensor<T,cn> cand,org;
#if     OUTPUT_THRDFILE
        ofstream thrdfile("./thrdfile.txt",std::ios::out|std::ios::app);
#endif
        Tensor<T,cn> T_left_norm, T_up_norm, B_up_norm,B_left_norm;
        //int tempPh1=0;
        //for (auto iter=this->queueLen.begin(); iter!=this->queueLen.end();iter++)
       // {
        //  if (*iter)
        //    tempPh1++;
        //}
        //Ph1 = double(tempPh1)/double(this->queueLen.size());
        //Ph0 = 1 - Ph1;
        //Ph1 = this->L1Model[1].size();
        //Ph0 = this->L1Model[0].size();
        //Ph1 = Ph1/(Ph1+Ph0);
        //Ph0 = 1 - Ph1;
        int recordcount = 0;
        if (qNode.offset().x==80&& qNode.offset().y==128&&qNode.size().height==16)
          cout<<"debug here"<<endl;
        BayesianRecord tempRecord1,tempRecord2;
        for (int x=offsetUp; x< offsetDown; x++)
          for (int y=offsetLeft; y< offsetRight; y++)
          {
            unsigned long key =x*ensemble.size().width+y;
            //key = key*ensemble.size().area()+x*ensemble.size().width+y;
            auto iter1 = L1Record.find(key);
            if (iter1!=L1Record.end())
            {
              tempRecord1.addData(iter1->second.data0,false);
              tempRecord1.addData(iter1->second.data1,true);
            }
            auto iter2 = L2Record.find(key);
            if (iter2!=L2Record.end())
            {
              tempRecord2.addData(iter2->second.data0,false);
              tempRecord2.addData(iter2->second.data1,true);
            }
          }

        //tempRecord.print();
        L1Model[0].clear();
        L1Model[1].clear();
        L1Model[0].addData(tempRecord1.data0);
        L1Model[1].addData(tempRecord1.data1);

        Ph0 = L1Model[0].size();
        Ph1 = L1Model[1].size();
        //if (qNode.size().height>16) //only consider 16x16 for Bayesian
        //  train1 = true;
        //else
        if
        (L1Model[0].size()<2||L1Model[1].size()<2||Ph1<10)
        {
            train1 = true;
        }
        else
        {
          train1 = false;
          Ph1 = Ph1/(Ph0+Ph1);
          Ph0 = 1-Ph1;
          L1Model[0].mle1D();
          L1Model[1].mle1D();
          //L1Model[0].mle();
          //cout<<L1Model[0].mu<<endl;
          //cout<<L1Model[0].igamma<<endl;
          //L1Model[1].mle();
          if (L1Model[0].lambda==0||L1Model[1].lambda==0)//low rank
            train1 = true;
          //cout<<L1Model[1].mu<<endl;
          //cout<<L1Model[1].igamma<<endl;
          eta1 = double(Ph0)/double(Ph1)/2;
        }

        L2Model[0].clear();
        L2Model[1].clear();
        L2Model[0].addData(tempRecord2.data0);
        L2Model[1].addData(tempRecord2.data1);
        //tempRecord2.print();
        Ph0 = L2Model[0].size();
        Ph1 = L2Model[1].size();
        //20140213 LR 16 only
        //if (qNode.size().height>16)
        //  train2=true;
        //else if
        if(L2Model[0].size()<2||L2Model[1].size()<2||Ph1<10)
        {
            train2 = true;
        }
        else
        {
                  train2 = false;
                  Ph1 = Ph1/(Ph0+Ph1);
                  Ph0 = 1-Ph1;
                  L2Model[0].mle1D();
                  //cout<<L1Model[0].mu<<endl;
                  //cout<<L1Model[0].igamma<<endl;
                  L2Model[1].mle1D();
                  //cout<<L1Model[1].mu<<endl;
                  //cout<<L1Model[1].igamma<<endl;
                  eta2 = (double)Ph0/double(Ph1);
                  if (L2Model[0].detSigma==0||L2Model[1].detSigma==0)//low rank
                    train2 = true;
        }
//        Ph0 = 0.5462;
//        Ph1 = 0.4538;
//        L1Model[0].lambda =3.39851394298828;//3.6371;
//        L1Model[1].lambda = 5.29223331976362;//1.4592;
//        L1Model[0].mu = (Mat_<double>(1,2)<<-0.0175409670043977,0.111436732926657);//0.3133,0.4540);
//        L1Model[1].mu = (Mat_<double>(1,2)<<0.440592517010730,0.391815187735022);//-0.0374,-0.0625);
//        L1Model[0].gamma=(Mat_<double>(2,2)<<1.17058747370365,0.768749840589956,0.768749840589956,1.35912638153674);//1.31827907937583,0.787828165735224,0.787828165735224,1.22938552547847);
//        L1Model[1].gamma=(Mat_<double>(2,2)<<1.16387250740165,0.773056143087246,0.773056143087246,1.37267251370308);//1.21482731227115,0.546032825848818,0.546032825848818,1.06858961252486);
//        L1Model[0].igamma = L1Model[0].gamma.inv();
//        L1Model[1].igamma = L1Model[1].gamma.inv();
//        eta1 = Ph0/Ph1/2;
//        if (qNode.offset().x > 0  && qNode.offset().y>0)
        {
            T_up =  qNode.upBound.Clone();// qNode.upBound.Crop(Point3i(0,qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
            T_left = 	 qNode.leftBound.Clone();
//            TT.SetBlock(T_left.Transpose());
//            TT.SetBlock(Point3i(0,T_left.size().width,0),T_up);
//            var_tag = TT.Var()[0];
            var_tar_left = T_left.Var()[0];
            var_tar_up = T_up.Var()[0];
            if (var_tar_left<0.01)
              var_tar_left = 0.01;
            if (var_tar_up<0.01)
              var_tar_up = 0.01;
            //has_left = true;
           // has_up = true;
           // muT = TT.Mean();
           // T_up_norm = T_up-muT;
           // T_left_norm = T_left-muT;
#if     OUTPUT_THRDFILE
           thrdfile << "org: "<<qNode.offset().x<<","<<qNode.offset().y<<"var: left "<<var_tar_left<<", up "<<var_tar_up<<"----------------------------"<<endl;
#endif
         }
//        else if (qNode.offset().x ==0&&qNode.offset().y>0) //only have left
//          {
//            T_left = qNode.leftBound.Crop(Point3i(qNode.overlap().height,0,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0));
//            var_tag = T_left.Var()[0];
//            muT = T_left.Mean();
//            var_tag_left = var_tag;
//            has_left = true;
//            T_left_norm = T_left-muT;
//          }
//        else if (qNode.offset().x>0&&qNode.offset().y==0) // only have up
//          {
//            T_up = qNode.upBound.Crop(Point3i(0,qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0));
//            var_tag = T_up.Var()[0];
//            var_tag_up = var_tag;
//            has_up = true;
//            muT = T_up.Mean();
//            T_up_norm=T_up-muT;
//          }
//        else
//          {
//            CV_Error(CV_StsNotImplemented,"no boundary for matching");
//            return sideMatchAddr;
//          }
        //cv::Mat gaussKernel = mylib::GenGaussKer(3,double(3)/6.0,CV_64F);
        //Tensor<T,cn> T0_low=T0.Filter2D(gaussKernel,FILTER_BOUND_VALID);
        //Tensor<T,cn> T1_low=T1.Filter2D(gaussKernel,FILTER_BOUND_VALID);
        //Tensor<T,cn> T0_Lap=T0.Laplacian();
        //Tensor<T,cn> T1_Lap=T1.Laplacian();
        //Tensor<T,cn> filterBD,lowBD;

#if PARALLEL_MATCHING
        int pnum = 4;
        vector<thread> threads;
#else
        Tensor<T,cn> tarSide(oSize);
        Tensor<T,cn> canSide(oSize);
        int pnum =1;// thread::hardware_concurrency();
#endif
        vector<LinkArray*> localqueue(pnum);
        int brows = (offsetDown - offsetUp)/searchStep.height/pnum;
        ensemble.Ref(Cube(qNode.offset(),qNode.size()),org);
        for (int p=0; p<pnum; p++)
        {
#if PARALLEL_MATCHING

          threads.push_back(thread([&](int p, int t){
          Tensor<T,cn> tarSide(oSize);
          Tensor<T,cn> canSide(oSize);
          Tensor<T,cn> VV(Size3(qNode.size().height*2+qNode.overlap().height,qNode.overlap().width,1));
          Tensor<T,cn> candLeftLocal,candUpLocal;
          localqueue[p] = new LinkArray();
#else
          localqueue[p] = queue;
#endif
          for (int x=offsetUp+p*brows; x< min(offsetDown,offsetUp+(p+1)*brows); x+=searchStep.height)
            for (int y=offsetLeft; y< offsetRight; y+=searchStep.width)
              {
                if (IsInsideCausalRegion(cv::Point3i(x,y,0),qNode,3))
                  {
                    diff=0;
                    N=0;
                    //20140127 use varratio
                    rst.Ref(Cube(Point3i(x,y,t),qNode.leftBound.size()),B_left);//candLeftLocal);
                    //rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,t),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),candUpLocal);
                    rst.Ref(Cube(Point3i(x,y,t),qNode.upBound.size()),B_up);//candUpLocal);

                    //VV.SetBlock(candUpLocal.Transpose());
                    //VV.SetBlock(Point3i(tarUp.size().width,0,0),candLeftLocal);
                    //get local variance
                    //Tensor<T,cn> tempblk;
                    //rst.Ref(Cube(Point3i(x,y,t),qNode.size()+qNode.overlap()),tempblk);
                    //double Vb = VV.Var()[0];
                    //double Vb = tempblk.Var()[0];//local variance is close to target side variance
                    //double delta = abs(1-Vb/Vt);
                    //clip
                   // double delta_thrd = 20;//10*log10(100);
                    //delta > delta_thrd? delta = 1 : delta=delta/delta_thrd;
                    //double delta = abs(Vb - Vt);
                    //cout<<"delta is "<<delta<<endl;
//                    if (has_left && has_up)
                      {
//                        rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                        rst.Ref(Cube(Point3i(x,y,0),qNode.leftBound.size()),B_left);
//                        VV.SetBlock(B_left.Transpose());
//                        VV.SetBlock(Point3i(0,T_left.size().width,0),B_up);
//                        var_blk = VV.Var()[0];
                          var_blk_up= B_up.Var()[0];
                          var_blk_left = B_left.Var()[0];
                          if (var_blk_left<0.01)
                            var_blk_left = 0.01;
                          if (var_blk_up<0.01)
                            var_blk_up = 0.01;
                          diff_up = var_blk_up - var_tar_up;
                          diff_left = var_blk_left - var_tar_left;
                          ratio_left = std::log(var_blk_left/var_tar_left);
                          ratio_up = std::log(var_blk_up/var_tar_up);
                          if (ratio_left<-10)
                            ratio_left = -10;
                          if (ratio_left>10)
                            ratio_left = 10;
                          if (ratio_up<-10)
                            ratio_up = -10;
                          if (ratio_up>10)
                            ratio_up = 10;
                          //ratio_up = (var_blk_up < var_tar_up)?abs(1-var_blk_up/var_tar_up):abs(1-var_tar_up/var_blk_up);
                          //ratio_left = (var_blk_left<var_tar_left)? abs(1-var_blk_left/var_tar_left):abs(1-var_tar_left/var_blk_left);
                          dec1 = false;
                          //if (lenH0>=L1_train_max&& lenH1>=L1_train_max && train==true)
                          //if (lenH0>50&&lenH1>50&&lenH0+lenH1>this->L1_train_max2)
                          //if (train&&L1Model[0].size()>=10&&L1Model[1].size()>=10&&L1Model[1].size()+L1Model[0].size()>=this->L1_train_max)
                          //{
                            //compute mle
                          //    train = false;
                          //    L1Model[0].mle();//H0 left
                          //    L1Model[1].mle();//H0 up
                              //L1Model[2].mle();
                             // L1Model[3].mle();
                          //}
                          var_ratio.at<double>(0,0) = ratio_left;
                          var_ratio.at<double>(0,1) = ratio_up;
                          likelihood1 = -1;//mark as not used
                          if (!train1)
                          {
                            //compute llr
                            //likelihood1 =  (L1Model[2].pdf(diff_left))/(L1Model[0].pdf(diff_left));
                            //likelihood2 =  (L1Model[3].pdf(diff_left))/(L1Model[1].pdf(diff_left));
                            //eta1 = double(L1Model[0].N)/double(L1Model[0].N+L1Model[2].N);
                            //eta1 = double(this->lenH0)/double(this->lenH1)/2;


                            //cout<<var_ratio<<endl;
                            dec1 = false;
                            likelihood1_l = L1Model[1].pdf1D(ratio_left)/L1Model[0].pdf1D(ratio_left);
                            likelihood1_u = L1Model[1].pdf1D(ratio_up)/L1Model[0].pdf1D(ratio_up);
                            //likelihood1 = L1Model[1].pdf(var_ratio)/L1Model[0].pdf(var_ratio);
                            //if (likelihood1>eta1)
                            if (likelihood1_l>eta1&&likelihood1_u>eta1)
                                dec1 = true;
                            /*
                            double pmiss = 0.4;
                            dec1 = false;
                            if (L1Model[0].mu<L1Model[2].mu) // H0 on the left of H1
                            {
                              lambda_l1 = L1Model[2].cdfinv(pmiss);
                              if (diff_left >=lambda_l1)
                                dec1 = true;
                            }
                            else
                            {
                              lambda_l2 = L1Model[2].cdfinv(1-pmiss);
                              if (diff_left < lambda_l2)
                                dec1 = true;
                            }

                            //lambda_l1 = L1Model[2].cdfinv(pmiss/2);//left
                            //lambda_l2 = L1Model[2].cdfinv(1-pmiss/2);//right
                            //if (diff_left<lambda_l1||diff_left>=lambda_l2)
                            //  dec1 = true;

                            //eta1 = double(lenH0)/double(lenH1)/2;
                            //eta1 = 0.5;
                            //if (likelihood1>=lambda_l)
                            //  dec1 = true;
                            //lambda_u1 = L1Model[3].cdfinv(pmiss/2);//left
                            //lambda_u2 = L1Model[3].cdfinv(1-pmiss/2);//right

                            if (L1Model[1].mu<L1Model[3].mu)//H0 on left of H1
                            {
                                lambda_u1 = L1Model[3].cdfinv(pmiss);
                                if(diff_up >= lambda_u1)
                                  dec1 = true&dec1;
                                else
                                  dec1 = false&dec1;
                            }
                            else
                            {
                                lambda_u2 = L1Model[3].cdfinv(1-pmiss);
                                if (diff_up<lambda_u2)
                                  dec1 = true&dec1;
                                else
                                  dec1 = false&dec1;
                            }
                            */
                            //if (diff_up<lambda_u1||diff_up>=lambda_u2)
                            //  dec1 = dec1&&true;
                            //else
                            //  dec1 = dec1&&false;
                            //if (likelihood2>=lambda_l)
                            //  dec1 = true&&dec1;
                          }
                          //muB_up = B_up.Mean();
                         // muB_left = B_left.Mean();
//                        muB = VV.Mean();
                         // B_up_norm=B_up-muB_up;
                         // B_left_norm = B_left - muB_left;
                          //diff_up = abs(var_blk_up - var_tar_up);
                          //diff_left = abs(var_blk_left-var_tar_left);
                          //diff_up = (var_blk_up < var_tar_up)?abs(1-var_blk_up/var_tar_up):abs(1-var_tar_up/var_blk_up);
                          //diff_left = (var_blk_left<var_tar_left)? abs(1-var_blk_left/var_tar_left):abs(1-var_tar_left/var_blk_left);
                          //for the small variance, I want to make the diff less sensitive 20140130
                          //20140205 tempdisable
                          //diff_up*=adj_up;
                          //diff_left*=adj_left;
#if     OUTPUT_THRDFILE
                          thrdfile<<"-- can: "<<x<<", "<<y<<"var: "<<var_blk_left<<", "<<var_blk_up<<" -- layer1: "<< var_ratio.at<double>(0,0) <<", "<<var_ratio.at<double>(0,1)<<"train1? "<<train1<<" lr1: "<<likelihood1<<" eta1: "<<eta1<<" dec1: "<<dec1<<endl;
#endif
                      }
//                    else if (has_left)
//                      {
//                        rst.Ref(Cube(Point3i(x+qNode.overlap().height,y,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0)),B_left);
//                        var_blk = B_left.Var()[0];
//                        var_blk_left = var_blk;
//                        muB =B_left.Mean();
//                        B_left_norm = B_left - muB;
//                      }
//                    else if (has_up)
//                      {
//                        rst.Ref(Cube(Point3i(x,y+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                        var_blk = B_up.Var()[0];
//                        var_blk_up = var_blk;
//                        muB = B_up.Mean();

//                        B_up_norm = B_up - muB;
//                      }
//                    //diff = abs(var_blk-var_tag);
//                    //diff = abs(log10(var_blk/var_tag));//change to this, Nov 10, 2012
//                    double diff_left=0,diff_up=0;
//                    if (has_left)
//                    {
//                        if (abs(var_blk_left-var_tag_left)<0.1)
//                            diff_left =0;
//                        else
//                            diff_left = abs(log10((var_blk_left-var_tag_left)/var_tag_left));
//                    }
//                    if (has_up)
//                    {
//                        if (abs(var_blk_up-var_tag_up)<0.1)
//                            diff_left = 0;
//                        else
//                            diff_up = abs(log10((var_blk_up-var_tag_up)/var_tag_up));
//                    }
                    //cout<<diff_left<<","<<diff_up<<endl;
                    //if (delta<=varThrd1) //20140127 first trim (coarest) discard large var ratio
                    //test, make case 1 always fail and only test case 2
                    //dec1 = false;train1 = true;
                    if (dec1||(train1&&abs(ratio_left)<log(3.0)&&abs(ratio_up)<log(3.0))) //varThrd1) //change thred, Nov 10, 2012
                    {//gj01132013 change threahold from 0.5 to 0.2, gj01142013, change thred co
                        //varqueue->compareInsert(diff,cv::Point3i(x,y,t));
                        //varList.push_back(cv::Point3i(x,y,t));
                        //round 2, search for nbs


                        for (int nbx=-searchStep.height/2; nbx< searchStep.height/2; nbx+=1) //change from < to <=, Nov 10,2012, wrong ! 20140130
                          for (int nby=-searchStep.width/2; nby< searchStep.width/2; nby+=1)
                            {
                              if (x+nbx<offsetUp||x+nbx>offsetDown||y+nby<offsetLeft||y+nby>offsetRight) //20131227, make sure searching in neighborhood do not voilate search range constrain
                                continue;
                              if (x+nbx<offsetUp+p*brows||x+nbx>=offsetUp+(p+1)*brows) //no duplicate search in different threads
                                continue;
                              if (IsInsideCausalRegion(cv::Point3i(nbx+x,nby+y,0),qNode,3))
                                {
                                  rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.leftBound.size()),B_left);//candLeftLocal);
                                  rst.Ref(Cube(Point3i(x+nbx,y+nby,t),qNode.upBound.size()),B_up);//candUpLocal);
                                  double mse_up = metric::ComputeMSE(T_up,B_up);
                                  double mse_left = metric::ComputeMSE(T_left,B_left);//metric::ComputeMSE(VV,TT);
                                  likelihood2 = -1;
                                  if (!train2)
                                  {
                                      mse_up = log(mse_up);
                                      mse_left = log(mse_left);
                                      logmse.at<double>(0,0) = mse_left;
                                      logmse.at<double>(0,1) = mse_up;
                                      likelihood2_l = L2Model[1].pdf1D(mse_left)/L2Model[0].pdf1D(mse_left);
                                      likelihood2_u = L2Model[1].pdf1D(mse_up)/L2Model[0].pdf1D(mse_up);
                                      //likelihood2 = L2Model[1].pdf(logmse)/L2Model[0].pdf(logmse);
                                      dec2 = false;
                                      //if (likelihood2>eta2)
                                      if (likelihood2_l>eta2&&likelihood2_u>eta2)
                                        dec2 = true;
                                  }



                                  //mse /= 65025; //normalize

                                  /*
                                  if (train)
                                  {
                                    rst.Ref(Cube(Point3i(x+nbx+qNode.overlap().height,y+nby+qNode.overlap().width,0),qNode.size()),cand);
                                    double temp = metric::Compare(org,cand,CompareCriteria::COMPARE_CRITERIA_SSIM, Size3(16,16,1),Size3(16,16,1), 3, 4, (int)FilterBoundary::FILTER_BOUND_EXTEND, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_NEW_L1,0,false);
                                    if (temp<=varThrd1)//here varThrd1 = stsim2 score
                                    {
                                      L1Model[0].addData(diff_left);
                                      L1Model[1].addData(diff_up);
                                    }
                                    else
                                    {
                                      L1Model[2].addData(diff_left);
                                      L1Model[3].addData(diff_up);
                                    }
                                    this->L1_train_len++;
                                  }
                                  */
                                  if (!train2){
                                    //if (mse_left>matching_thrd||mse_up>matching_thrd) //2nd trim: discard large side mse candidate
                                    if (!dec2)
                                      continue;
                                  }
                                  else
                                  {
                                    if (mse_left>(matching_thrd*2)||mse_up>(matching_thrd*2)) //if no successfull MTC in searching area, relax the layer2
                                      continue;

                                  }
#if     OUTPUT_THRDFILE
                                  thrdfile<<"---layer 2:"<<mse_left<<", "<<mse_up<<"train2? "<<train2<<" lr2: "<<likelihood2<<" eta2: "<<eta2<<" dec2: "<<dec2<<endl;
#endif
                                  N=0;
                                  double score=0;
                                  double final_score=1;
                                  int count=0;
                                  for (int m = 0; m<qNode.leftBound.size().height-shift_x; m+=oSize.height)
                                  {
                                    ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(m,0,1),oSize),tarSide);
                                    rst.Ref(Cube(Point3i(x+nbx+m+shift_x,y+nby+shift_y,t),oSize),canSide);
                                    score=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize,oSize, 3, 1, (int)FilterBoundary::FILTER_BOUND_EXTEND, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                                    if (score<final_score)
                                      final_score= score;
                                    count++;
                                  }
                                  for (int n = oSize.width; n<qNode.upBound.size().width-shift_y; n+=oSize.width)
                                  {
                                      ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(0,n,1),oSize),tarSide);
                                      rst.Ref(Cube(Point3i(x+nbx+shift_x,y+nby+n+shift_y,t),oSize),canSide);
                                      score=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM,
                                                            oSize, oSize,3, 1, (int)FilterBoundary::FILTER_BOUND_EXTEND,
                                                            (int)FeaturePoolType::FEATURE_POOL_MIN,
                                                            (int)MetricModifier::STSIM2_PART,0,false);
                                      if (score<final_score)
                                        final_score= score;
                                      count++;
                                  }
                                  //score = score/count;//average;
                                  //got the min of all scores
                                  #if     OUTPUT_THRDFILE
                                  thrdfile<<"---layer 3:"<<x+nbx<<", "<<y+nby<<" score: "<< score<<endl;
                                  if (qNode.offset().x == DEBUG_X && qNode.offset().y==DEBUG_Y&&x+nbx==18 && y+nby==292)
                                    cout<<"x,y "<<x<<", "<<y<<endl;
                                  #endif
                                  localqueue[p]->compareInsert(final_score,cv::Point3i(x+nbx,y+nby,t));
                                  //if (!train){//no limit of candidate when training
                                  if (localqueue[p]->getLength()>candidNum/pnum&&candidNum>0)
                                    localqueue[p]->pop();
                                  //}
//                                  if (has_left && has_up)
//                                    {
//                                      mse_diff = 0;
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby,0),qNode.leftBound.size()),B_left);
//#if DEBUG_20131226
//                                      B_up.Print();
//                                      B_left.Print();
//#endif
//                                      double diff_up =   metric::Compare(T_up,B_up,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      double diff_left = metric::Compare(T_left,B_left,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      //Compare normalize the value, it's not correct
//                                      //so fix it
//                                      diff_up*=(double)T_up.size().area();
//                                      diff_left*=(double)T_left.size().area();
//                                      mse_diff = diff_up+diff_left;
//                                      N+= T_up.size().volumn();
//                                      N+= T_left.size().volumn();
//                                    }
//                                  else if (has_left)
//                                    {
//                                      mse_diff = 0;
//                                      rst.Ref(Cube(Point3i(x+nbx+qNode.overlap().height,y+nby,0),qNode.leftBound.size()-Size3(qNode.overlap().height,0,0)),B_left);
//                                      mse_diff+= metric::Compare(T_left_norm, B_left_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      N+= T_left.size().volumn();
//                                    }
//                                  else if (has_up)
//                                    {
//                                      mse_diff=0;
//                                      rst.Ref(Cube(Point3i(x+nbx,y+nby+qNode.overlap().width,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)),B_up);
//                                      mse_diff+= metric::Compare(T_up_norm,B_up_norm,CompareCriteria::COMPARE_CRITERIA_MSE);
//                                      N+= T_up.size().volumn();
//                                    }
//                                  mse_diff/=double(N);
//#if DEBUG_20131226
//                                  cout<<x+nbx<<","<<y+nby<<endl;
//                                  cout<<mse_diff<<endl;
//                                  T_up.Print();
//                                  T_left.Print();
//#endif
//                                  mse_diff/=65025;//normalize 255^2
//#if DEBUG_20131226
//                                  cout<<mse_diff<<endl;
//#endif
//                                  if (mse_diff<=matching_thrd)
//                                    queue->compareInsert(-mse_diff,Point3i(x+nbx,y+nby,0));
//                                  if (queue->getLength()>candidNum&&candidNum>0)
//                                    queue->pop();
                                }
                            }
                      }
                  }
              }
#if PARALLEL_MATCHING
       },p,t));
#endif
     }
#if PARALLEL_MATCHING
     for (int p=0; p<pnum; p++)
     {
       auto& thread = threads[p];
       thread.join();
       for (int ii=0; ii<localqueue[p]->getLength(); ii++)
       {
           queue->compareInsert(localqueue[p]->getData(ii),localqueue[p]->GetAddress(ii));
       }
       delete localqueue[p];
     }
#endif
     #if     OUTPUT_THRDFILE
        thrdfile.close();
#endif
   }

   else if (matching_method==MatchingMethod::MATCHING_STSIM)
   {
      Size3 oSize(max(qNode.overlap().height,subWinSize.height),max(qNode.overlap().width,subWinSize.width),1);
      if (offsetUp + qNode.overlap().height-oSize.height<0)
        offsetUp += oSize.height-qNode.overlap().height;
      if (offsetLeft+qNode.overlap().width -oSize.width<0)
        offsetLeft += oSize.width-qNode.overlap().width;
      int shift_x = 0;
      int shift_y = 0;
      if (qNode.overlap().height<subWinSize.height)
        shift_x = qNode.overlap().height-oSize.height;
      if (qNode.overlap().width<subWinSize.width)
        shift_y = qNode.overlap().width-oSize.width;
#if PARALLEL_MATCHING
        int pnum = 4;
        vector<thread> threads;
#else
        Tensor<T,cn> tarSide(oSize);
        Tensor<T,cn> canSide(oSize);
        int pnum =1;// thread::hardware_concurrency();
#endif
      int brows = (offsetRight - offsetLeft)/searchStep.width/pnum;
      for (int t=0; t< searchRegion.depth;t+=searchStep.depth)
      {
        for (int p=0; p<pnum; p++)
        {
#if PARALLEL_MATCHING
       threads.push_back(thread([&](int p, int t){
        Tensor<T,cn> tarSide(oSize);
        Tensor<T,cn> canSide(oSize);
#endif
        for (int y = offsetLeft +p*brows; y<min(offsetRight, offsetLeft+(p+1)*brows);y+=searchStep.width)
        {
          for (int x=offsetUp; x<= offsetDown; x+=searchStep.height)
          {
            if (IsInsideCausalRegion(cv::Point3i(x,y,t),qNode,3/*3*/))//20130605  change from 3 to 2
            {
              double score=0;
              int count=0;
              for (int m = 0; m<qNode.leftBound.size().height-shift_x; m+=oSize.height)
              {
               ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(m,0,1),oSize),tarSide);
                rst.Ref(Cube(Point3i(x+m+shift_x,y+shift_y,t),oSize),canSide);
                score+=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize,oSize, 3,1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_NEW_L1,0,false);
                count++;
              }
              for (int n = oSize.width; n<qNode.upBound.size().width-shift_y; n+=oSize.width)
              {
                  ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(0,n,1),oSize),tarSide);
                  rst.Ref(Cube(Point3i(x+shift_x,y+n+shift_y,t),oSize),canSide);
                  score+=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize, oSize,3,1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_NEW_L1,0,false);
                  count++;
              }
              score = score/count;//average;
              queue->compareInsert(score,cv::Point3i(x,y,t));
              if (/*myqueue*/queue->getLength()>candidNum&&candidNum>0)
                /*myqueue*/queue->pop();
            }
         }
       }
#if PARALLEL_MATCHING
       },p,t));
#endif
        }
      }
#if PARALLEL_MATCHING
      for (auto& thread : threads){
        thread.join();
      }
#endif
    }
    else if (matching_method==MatchingMethod::MATCHING_STSIM_PART)
    {
        Size3 oSize(max(qNode.overlap().height,subWinSize.height),max(qNode.overlap().width,subWinSize.width),1);
        if (offsetUp + qNode.overlap().height-oSize.height<0)
          offsetUp += oSize.height-qNode.overlap().height;
        if (offsetLeft+qNode.overlap().width -oSize.width<0)
          offsetLeft += oSize.width-qNode.overlap().width;
        int shift_x = 0;
        int shift_y = 0;
        if (qNode.overlap().height<subWinSize.height)
          shift_x = qNode.overlap().height-oSize.height;
        if (qNode.overlap().width<subWinSize.width)
          shift_y = qNode.overlap().width-oSize.width;
  #if PARALLEL_MATCHING
          int pnum = 4;
          vector<thread> threads;
  #else
          Tensor<T,cn> tarSide(oSize);
          Tensor<T,cn> canSide(oSize);
          int pnum =1;// thread::hardware_concurrency();
  #endif
        int brows = (offsetRight - offsetLeft)/searchStep.width/pnum;
        for (int t=0; t< searchRegion.depth;t+=searchStep.depth)
        {
          for (int p=0; p<pnum; p++)
          {
  #if PARALLEL_MATCHING
          threads.push_back(thread([&](int p, int t){
          Tensor<T,cn> tarSide(oSize);
          Tensor<T,cn> canSide(oSize);
  #endif
          for (int y = offsetLeft +p*brows; y<min(offsetRight, offsetLeft+(p+1)*brows);y+=searchStep.width)
          {
            for (int x=offsetUp; x<= offsetDown; x+=searchStep.height)
            {
              if (IsInsideCausalRegion(cv::Point3i(x,y,t),qNode,3/*3*/))//20130605  change from 3 to 2
              {
                double score=0;
                int count=0;
                for (int m = 0; m<qNode.leftBound.size().height-shift_x; m+=oSize.height)
                {
                  ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(m,0,1),oSize),tarSide);
                  rst.Ref(Cube(Point3i(x+m+shift_x,y+shift_y,t),oSize),canSide);
                  score+=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize,oSize, 3, 1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                  count++;
                }
                for (int n = oSize.width; n<qNode.upBound.size().width-shift_y; n+=oSize.width)
                {
                    ensemble.Ref(Cube(qNode.offset()-oSize.Point3()+Point3i(0,n,1),oSize),tarSide);
                    rst.Ref(Cube(Point3i(x+shift_x,y+n+shift_y,t),oSize),canSide);
                    score+=metric::Compare(tarSide,canSide,CompareCriteria::COMPARE_CRITERIA_SSIM, oSize, oSize,3, 1, (int)FilterBoundary::FILTER_BOUND_FULL, (int)FeaturePoolType::FEATURE_POOL_MIN, (int)MetricModifier::STSIM2_PART,0,false);
                    count++;
                }
                score = score/count;//average;
                queue->compareInsert(score,cv::Point3i(x,y,t));
                if (/*myqueue*/queue->getLength()>candidNum&&candidNum>0)
                  /*myqueue*/queue->pop();
              }
            }
          }
  #if PARALLEL_MATCHING
         },p,t));
  #endif
          }
        }
  #if PARALLEL_MATCHING
        for (auto& thread : threads){
          thread.join();
        }
  #endif
    }

    if (queue->getLength()>0)
      {
        logfile<<"Tar: ("<<qNode.offset().x<<","<<qNode.offset().y<<")"<<endl;
        cout<<"Tar: ("<<qNode.offset().x<<","<<qNode.offset().y<<")"<<endl;
        for (int k=0; k<queue->getLength();k++)
          {
            if (queue->GetAddress()[k]!=Point3i(-1,-1,-1))
              {
                sideMatchAddr.push_back(pair<Point3i,double>(queue->GetAddress()[k],queue->getData(k)));
                logfile<<"       <<<<<<<<  Cand: ("<<sideMatchAddr[k].first.x<<","<<sideMatchAddr[k].first.y<<") Dist: "<<queue->getData(k)<<endl;
                cout<<"       <<<<<<<<  Cand: ("<<sideMatchAddr[k].first.x<<","<<sideMatchAddr[k].first.y<<") Dist: "<<queue->getData(k)<<endl;
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

#include "Tester.h"
#include <stdio.h>
//#include "Steerable.h"
#include <iostream>
#include <sstream>
#include <bitset>
//#include "PoissonSolver.h"
#include <sys/types.h>
//#include <process.h>
//#include <Metric.h>

using namespace std;
using namespace tensor;
using namespace metric;
Tester::Tester(void)
{
}


Tester::~Tester(void)
{
}
/*
void Tester::TestInterp(string filename)
{
	typedef double T;
	const int cn=1;
	Tensor<T,cn> ensemble(filename);
	QNode<T,cn> node(ensemble,Size3(16,16,1),Point3i(16,16,0),Size3(4,4,0));
	node.LinearInterp(8);
	node.Display(0);
}

void Tester::TestRoi(string filename)
{
	typedef double T;
	const int cn=1;
	Tensor<T,cn> ensemble(filename);
	Tensor<T,cn> blk(Size3(20,4,1));
	Tensor<T,cn> part;

	part = blk(Cube(4,0,0,16,4,1));
	//part.SetBlock(Tensor_<T,cn>(part.size(),Vec<T,cn>(255)));
	blk.Display();
	part = ensemble(Cube(30,30,0,16,4,1));
	blk.Display(0);

	QNode<T,cn> node;
	//node.Ref(ensemble,Cube(16,16,0,16,16,1),Size3(4,4,0));
	//node.SetBlock(Tensor_<T,cn>(node.size()));
	//ensemble.Display();
}

void Tester::TestQuilting(string filename)
{
	typedef double T;
	const int cn=1;
	Tensor<T,cn> ensemble(filename);
	QNode<T,cn> A(ensemble,Size3(16,16,1),Point3i(30,30,0),Size3(4,4,0));
	QNode<T,cn> B(ensemble,Size3(16,16,1),Point3i(80,80,0),Size3(4,4,0));
	A.Compare(B,COMPARE_CRITERIA_SSIM);
	vector<vector<int>> seam = A.Quilting(B);
	//ensemble.Print();


}
void Tester::TestQuilting(string f1, string f2)
{
  Tensor<double,1> A(f1);
  Tensor<double,1> B(f2);
  Tensor<double,1> two(A.size().height,A.size().width*2,A.size().depth);//Size3(0,0,1));
  two.SetBlock(A);
  two.SetBlock(Point3i(0,A.size().width,0),B);
  CV_Assert(A.size()==B.size());
  Size3 osize = A.size()/8;
  Size3 bsize = A.size()-osize*2;
  QNode<double,1> q1(A,bsize,osize.Point3(),osize);
  QNode<double,1> q2(B,bsize,osize.Point3(),osize);
  vector<Tensor<double,1>> b1,b2;
  Tensor<double,1> leftBound;
  b1.push_back(q1.rightBound);
  B.Ref(Cube(q2.leftBound.offset(),Size3(B.size().height,osize.width,1)),leftBound);
  b2.push_back(leftBound);
  vector<Tensor<double,1>> all;
  vector<vector<int>> seam = q2.Quilting(b2,b1,LEFT);
  Tensor<double,1> rst(A.size().height,A.size().width*2-osize.width, A.size().depth);
  rst.SetBlock(A(Cube(Point3i(0,0,0),Size3(A.size().height,A.size().width-osize.width,A.size().depth))));
  rst.SetBlock(Point3i(0,A.size().width-osize.width,0),B);
  all.push_back(two);
  all.push_back(rst);
  A.DisplayAll(all,2,1);
  two.SaveBlock("../tempresults/two.tiff");
  rst.SaveBlock("../tempresults/combin_two.tiff");
}


void Tester::TestEmbedding(string f1, string f2)
{
  Tensor<double,1> A(f1);
  Tensor<double,1> B(f2);
  Tensor<double,1> C= A.Clone();
  CV_Assert(A.size()==B.size());
  Size3 osize = A.size()/4;
  Size3 bsize = A.size()-osize*2;
  QNode<double,1> q1(A,bsize,osize.Point3(),osize);
  QNode<double,1> q2(B,bsize,osize.Point3(),osize);
  QNode<double,1> q3(C,bsize,osize.Point3(),osize);
  vector<Tensor<double,1>> all;
  all.push_back(A.Clone());
  all.push_back(B.Clone());
  q1.Quilting(q2);
  q2.Quilting(q3);
 // all.push_back(q1.GetExtendTensor().Clone());
 // all.push_back(q2.GetExtendTensor().Clone());
  //  A.DisplayAll(all,2,2);
  vector<Tensor<double,1>> b1,b2,b3;
  b1.push_back(q1.rightBound(Cube(osize.height,0,0,bsize.height,osize.width,1)));
  b2.push_back(q2.rightBound(Cube(osize.height,0,0,bsize.height,osize.width,1)));
  b3.push_back(q3.rightBound(Cube(osize.height,0,0,bsize.height,osize.width,1)));
  q1.Quilting(b1,b2,LEFT);
  q1.rightBound.SetBlock(Point3i(osize.height,0,0),b1[0]);
  q2.Quilting(b2,b3,LEFT);
  q2.rightBound.SetBlock(Point3i(osize.height,0,0),b2[0]);
  // all.push_back(q1.GetExtendTensor().Clone());
  //all.push_back(q2.GetExtendTensor().Clone());
   //A.DisplayAll(all,3,2);
  b1.pop_back();
  b2.pop_back();
  b3.pop_back();
  b1.push_back(q1.lowBound(Cube(0,osize.width,0,osize.height,bsize.height,1)));
  b2.push_back(q2.lowBound(Cube(0,osize.width,0,osize.height,bsize.height,1)));
  b3.push_back(q3.lowBound(Cube(0,osize.width,0,osize.height,bsize.height,1)));
  q1.Quilting(b1,b2,DOWN);
  q1.lowBound.SetBlock(Point3i(0,osize.width,0),b1[0]);
  q2.Quilting(b2,b3,DOWN);
  q2.lowBound.SetBlock(Point3i(0,osize.width,0),b2[0]);
  all.push_back(q1.GetExtendTensor().Clone());
  all.push_back(q2.GetExtendTensor().Clone());
  A.DisplayAll(all,2,2);
  q1.GetExtendTensor().SaveBlock("../tempresults/q1.tiff");
  q2.GetExtendTensor().SaveBlock("../tempresults/q2.tiff");

}



void Tester::TestInitGrid(void)
{
	typedef double T;
	const int cn=1;
	//Tensor_<T,cn> ensemble("baboon.pgm");
	//vector<vector<vector<QNode<T,cn>>>> vg(1);
	//vg[0]=vector<vector<QNode<T,cn>>>(32,32);

	//for (int i=0; i< 32; i++)
	//	for (int j=0; j<32; j++)
	//		vg[0][i][j].Ref(ensemble,Cube(i*16,j*16,0,16,16,1),Size3(4,4,0));
	//
	//
	//vg[0][8][8].SetBlock(Tensor_<T,cn>(Size3(16,16,1)));
	//vg[0][12][16].SetBlock(Tensor_<T,cn>(Size3(16,16,1),Vec<T,cn>(255)));
	//ensemble.Display();

	QGrid<T,cn> g("baboon.pgm",Size3(16,16,1),Size3(4,4,0));

}

void Tester::TestBoundMatch(void)
{
	typedef double T;
	const int cn=1;
	QGrid<T,cn> g("baboon.pgm",Size3(16,16,1),Size3(4,4,0));
	Tensor<T,cn> ts("baboon.pgm");
	for (int i=0; i< 6; i++)
		for (int j=0; j< 32; j++)
		{
			g.SetNode(Point3i(i,j,0),ts,Cube(i*16,j*16,0,16,16,1),Size3(4,4,0));
			g.SetCausalMap(g.GetNode(Point3i(i,j,0)));
		}
	g.SetNode(Point3i(6,0,0),ts,Cube(6*16,0,0,16,16,1),Size3(4,4,0));
	g.SetCausalMap(g.GetNode(Point3i(6,0,0)));
	g.Display();
	vector<Point3i> candid = g.BoundaryMatching(g.GetNode(Point3i(6,1,0)));
}



void Tester::TestTPC(void)
{
	string filename = "baboon_part2.tif";
	double q[] = {13,10,8,6.5, 5};
	double t[] = {0.92,0.93,0.94,0.95};
	int m[] = {CODING_MODE_PQI,CODING_MODE_TPQI};

	float r = cv::fastAtan2(-4,-4);

	//another more general setup model;
	TPC tpc(filename,CODING_MODE_JPEG);
	//tpc.SetInitQStep(8);
	//tpc.seti
	//TPC tpc(filename,0.93,Size3(32,32,1),Size3(1,1,1),4,8);
	tpc.Coding();
	//
}


void Tester::TestVarLenFunction(void)
{
  Tensor<double,1> A(32,32,1,255);
  Tensor<double,1> B(32,32,1,128);
  double d = A.Compare(B,COMPARE_CRITERIA_MSE);
}
void Tester::TestSTSIMBorder(void)
{
  int tar_x =64; 
  int tar_y = 480;
  int can_x = 5;
  int can_y = 719;
  //ifstream fp("best_candExt_PLC.txt",ios::in);
  ifstream fp("cand_0.txt",ios::in);
  double temp;
  int count=0;
  int dbsize=64;
  int bsize=dbsize/2;
  int osize=bsize/4;
  Tensor<double,1> cand(dbsize,dbsize,1,0);
  Tensor<double,1> cand_ctr;// = cand.Clone();
  Tensor<double,1> cand_br;//.. =cand.Clone();
  Tensor<double,1> im("../images/woman_g.tiff");
  //Tensor<double,1> rst("../rst.tif");
  //Tensor<double,1> cand_from_rst =rst.GetBlock(Cube(can_x-bsize/2,can_y-bsize/2,0,dbsize,dbsize,1));
  while(count<dbsize*dbsize)
  {
    if(count%dbsize==0)
      cout<<endl;
    fp>>temp;
    cout<<temp<<", ";
    //cout<<endl
    cand(count/dbsize,count%dbsize,0)=temp;
    count++;
  }
  
  //cand.Display();
  //cand.Print();
  cand_ctr=cand(Cube(bsize/2,bsize/2,0,bsize,bsize,1));
  cand_br=cand(Cube(bsize/4,bsize/4,0,bsize+2*osize,bsize+2*osize,1));
  //cand.Display();
  //cand_br.Display();
  //cand_ctr.Display();
 // Tensor<double,1> im("../images/woman_g.tiff");
  Tensor<double,1> tar = im.GetBlock(Cube(tar_x-bsize/2,tar_y-bsize/2,0,dbsize,dbsize,1));
  Tensor<double,1> tar_ctr;
  Tensor<double,1> tar_br;
  tar_ctr = tar(Cube(bsize/2,bsize/2,0,bsize,bsize,1));
  tar_br = tar(Cube(bsize/4,bsize/4,0,bsize+2*osize,bsize+2*osize,1));
  //tar.Display();
  //tar_br.Display();
  //tar_ctr.Display();
  Size3 subWinSz(16,16,1);
  Size3 subWinSp(8,8,1);
  Tensor<double,1> cand_ctr_pad(dbsize,dbsize,1,0);
  Tensor<double,1> tar_ctr_pad(dbsize,dbsize,1,0);
  Tensor<double,1> tar_br_pad(dbsize,dbsize,1,0);
  Tensor<double,1> tar_br_pad2((osize*2+bsize)*2,(osize*2+bsize)*2,1,0);
  Tensor<double,1> cand_br_pad2 = tar_br_pad2.Clone();
  Tensor<double,1> cand_br_pad(dbsize,dbsize,1,0);
  cand_ctr_pad.SetBlock(Point3i(bsize,bsize,0),cand_ctr);
  cand_br_pad.SetBlock(Point3i(osize,osize,0),cand_br);
  tar_ctr_pad.SetBlock(Point3i(bsize,bsize,0),tar_ctr);
  tar_br_pad.SetBlock(Point3i(osize,osize,0),tar_br);
  cand_br_pad2.SetBlock(Point3i(osize+bsize/2,osize+bsize/2,0),cand_br);
  tar_br_pad2.SetBlock(Point3i(osize+bsize/2,osize+bsize/2,0),tar_br);
  //cand_from_rst.Display();
  //cand_br_pad.Display();
  //tar_br_pad.Display();
  vector<Tensor<double,1>> vts;
  //vts.push_back(cand_from_rst);
  vts.push_back(cand_br_pad);
  vts.push_back(tar_br_pad);
  cand.DisplayAll(vts,1,2);
  Tensor<double,1> cand_pad(dbsize*2,dbsize*2,1,0);
  cand_pad.SetBlock(Point3i(dbsize,dbsize,0),cand);
  Tensor<double,1> tar_pad(dbsize*2,dbsize*2,1,0);
  tar_pad.SetBlock(Point3i(dbsize,dbsize,0),tar);
  //tar.debugtrigger=true;
  double full = cand.ToComplex().ComputeSSIM(tar.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_FULL,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double half = cand.ToComplex().ComputeSSIM(tar.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_HALF,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double valid = cand.ToComplex().ComputeSSIM(tar.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_VALID,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double ctr_full = tar_ctr.ToComplex().ComputeSSIM(cand_ctr.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_FULL,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  //double ctr_half = tar_ctr.ToComplex().ComputeSSIM(cand_ctr.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_HALF,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double br_full = tar_br.ToComplex().ComputeSSIM(cand_br.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_FULL,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  //double br_valid = tar_br.ToComplex().ComputeSSIM(cand_br.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_VALID,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double pad_full = cand_pad.ToComplex().ComputeSSIM(tar_pad.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_FULL,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double pad_half = cand_pad.ToComplex().ComputeSSIM(tar_pad.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_HALF,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double ctr_pad_full = cand_ctr_pad.ToComplex().ComputeSSIM(tar_ctr_pad.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_FULL,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double ctr_pad_half = cand_ctr_pad.ToComplex().ComputeSSIM(tar_ctr_pad.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_HALF,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double br_pad_full = cand_br_pad.ToComplex().ComputeSSIM(tar_br_pad.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_FULL,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double br_pad_half = cand_br_pad.ToComplex().ComputeSSIM(tar_br_pad.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_HALF,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  tar_br_pad.debugtrigger=true;
  auto tester = tar_br_pad.ToComplex();
  double br_pad_valid = cand_br_pad.ToComplex().ComputeSSIM(tar_br_pad.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_VALID,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double br_pad2_full = cand_br_pad2.ToComplex().ComputeSSIM(tar_br_pad2.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_FULL,STSIM2_POOL_MIN,STSIM2_NEW_L1);
  double br_pad2_half = cand_br_pad2.ToComplex().ComputeSSIM(tar_br_pad2.ToComplex(),subWinSz,subWinSp,3,4,FILTER_BOUND_HALF,STSIM2_POOL_MIN,STSIM2_NEW_L1);
 
}
void Tester::TestMatching()
{
  int tar_x =64;
  int tar_y = 480;
  int candNum = 8;
  int bsize = 32;
  int osize = bsize/4;
  int dbsize = bsize*2;
  Tensor<double,1> dummy(dbsize,dbsize,1,0);
  vector<Tensor<double,1>> cands,bds;
  Tensor<double,1> im("../images/woman_g.tiff");
  Tensor<double,1> tar = im.GetBlock(Cube(tar_x-bsize/2,tar_y-bsize/2,0,dbsize,dbsize,1));
  dummy.SetBlock(Point3i(osize,osize,0),tar.Crop(Point3i(osize,osize,0),Size3(bsize+2*osize,bsize+2*osize,1)));
  cands.push_back(dummy.Clone());
  cands.push_back(dummy.Clone());
  Tensor<double,1> tarLeft = tar.Crop(Point3i(osize,osize,0),Size3(osize+bsize,osize,1));//
  Tensor<double,1> tarUp = tar.Crop(Point3i(osize,osize+osize,0),Size3(osize,bsize,1));
  bds.push_back(tarLeft);
  bds.push_back(tarUp);
  vector<double> varsUp,varsLeft;
  varsLeft.push_back(tarLeft.Var()[0]);
  varsUp.push_back(tarUp.Var()[0]);
    cout<<" Tar Left var: "<<varsLeft[0]<<", Up var: "<<varsUp[0]<<endl;
  for (int i=0;i<candNum;i++)
  {
    string no = boost::lexical_cast<string>(i);
    Tensor<double,1> temp = mylib::readMatFromTxt("cand_PLC_"+no+".txt",dbsize,dbsize);
    dummy.SetBlock(Point3i(osize,osize,0),temp.Crop(Point3i(osize,osize,0),Size3(bsize+2*osize,bsize+2*osize,1)));
    cands.push_back(dummy.Clone());
    Tensor<double,1> candLeft = temp.Crop(Point3i(osize,osize,0),Size3(osize+bsize,osize,1));//
    Tensor<double,1> candUp = temp.Crop(Point3i(osize,osize+osize,0),Size3(osize,bsize,1));
    bds.push_back(candLeft.Clone());
    bds.push_back(candUp.Clone());
    varsLeft.push_back(candLeft.Var()[0]);
    varsUp.push_back(candUp.Var()[0]);
    cout<<"    >>>>>> Cand Left var: "<<varsLeft[i+1]<<", Up var: "<<varsUp[i+1]<<", delta= " << 
      abs(varsLeft[0]-varsLeft[i+1])/varsLeft[0]<< " and " <<abs(varsUp[0]-varsUp[i+1])/varsUp[0]<<
      " pooling: "<< sqrt((abs(varsLeft[0]-varsLeft[i+1])/varsLeft[0])*(abs(varsUp[0]-varsUp[i+1])/varsUp[0]))<<
      ", NMSE = "<<((candLeft-tarLeft).Pow(2.0).Sum()+(candUp-tarUp).Pow(2.0).Sum())[0]/double(candUp.size().area()+candLeft.size().area())/255/255<<endl;

  }

  dummy.DisplayAll(cands,5,2);
  dummy.DisplayAll(bds,9,2);
}




void Tester::TestLighting(void)
{
	Tensor<double,1> A("BBB_gray_p1.tif");
	Tensor<double,1> B("BBB_gray_p2.tif");
	A.LSFitting().Print();
	A.BuildLightingPlane(A.LSFitting()).SaveBlock(".\\Mar11\\T_L.tif");
	B.LSFitting().Print();
	B.BuildLightingPlane(B.LSFitting()).SaveBlock(".\\Mar11\\C_L.tif");

	(A.LSFitting()-B.LSFitting()).Print();
	Tensor<double,1> Diff(A.size());
	Diff = Diff.BuildLightingPlane(A.LSFitting() - B.LSFitting());
	Diff.SaveBlock(".\\Mar11\\D_L.tif");
	B = B+Diff;
	B.SaveBlock(".\\Mar11\\corrected_candid.tif");
 	B.LightingCorrection(A);
	B.Display(0);
//	Tensor_<double,1> rA  = A.LSFitting(1);
}

void Tester::TestBoundary(void)
{
	Tensor<double,1> im("baboon.pgm");
	QNode<double,1> q(im,Size3(32,32,1),Point3i(64,64,0),Size3(8,8,0));
	q.GetBoundaryLeft().Print();
	q.GetBoundaryLeft().SetBlock(Tensor<double,1>(q.GetBoundaryLeft().size()));
	q.GetBoundary(1).Print();
	q.GetBoundaryLeft().Print();
	im.Display();
}
*/
void Tester::TestSSIM(string file1, string file2)
{
	typedef double T;
	const int cn=1;
	Tensor<T,cn> A(file1);
	A.Print();
	Tensor<T,cn> B(file2);
	B.Print();
	//A.SetSubWinSize(Size3(16,16,1));
	//A.SetSubWinStep(Size3(1,1,1));
	ComputeStatistics(A,Size3(16,16,1),Size3(16,16,1));
	double ss = Compare(A,B,CompareCriteria::COMPARE_CRITERIA_SSIM,Size3(16,16,1),Size3(16,16,1),3,4,(int)FilterBoundary::FILTER_BOUND_HALF,(int)FeaturePoolType::FEATURE_POOL_MIN,(int)MetricModifier::STSIM2_BASELINE,false,true);
	cout<<ss<<endl;
}

void Tester::TestFilter2D(void)
{
	typedef double T;
	Tensor<T,1> im("/home/guoxin/Projects/MTC/data/baboon.pgm");
	im.Display();
	Tensor<T,1> part = im.Crop(Point3i(32,32,0),Size3(32,32,1));
	Steerable sp;
	sp.buildSCFpyr(part.ToComplex(),3,4,1,false);
	vector<Tensor<T,2>>& pyr = sp.getSpaceDomainPyr();
	cv::Mat gauss = mylib::GenGaussKer(8,8.0/6.0,CV_64F);
	Tensor<T,2> dst(part.size());
	cout<<pyr[0].size()<<", "<<pyr[0].channels()<<endl;
	auto temp = pyr[0];
	temp.Print();
	mylib::DisplayMat(gauss);
	//cv::flip(gauss,gauss,-1);
	cv::filter2D(pyr[0][0],dst[0],-1,gauss);
	dst= dst.Crop(Point3i(4,4,0),Size3(25,25,1));
	dst.Print();

	im.Filter2D(gauss).Display();

}
void Tester::TestMetric()
{
   Metric mc;
   mc.subwinSize = Size3(16,16,1);
   mc.subwinStep = Size3(16,16,1);
   //mc.subsample = true;
   //mc.changeWin=true;
   //mc.searchPath =  "H:/Code/Matlab/metric trainning/data/";
   string scorefilepath = "../../Size/subtestoutput_corbis64_inter.txt";
   mc.searchPath = "../../../data/totest/dist_corbis_64_inter/";
   mc.trainMetirc(mc.searchPath,scorefilepath);
}
void Tester::StudyMetricFeatures()
{
    Metric mc;
    mc.subwinSize = Size3(16,16,1);
    mc.subwinStep = Size3(16,16,1);
    //mc.subsample = true;
    //mc.changeWin=true;
    //mc.searchPath =  "H:/Code/Matlab/metric trainning/data/";
    string logfilepath = "./everything/everything.txt";
    mc.searchPath = "./everything/";
    //mc.trainSTSIM2Weights(mc.searchPath,scorefilepath);
    mc.studyMetricFeature(logfilepath);
}

void Tester::TestGranulateTrain()
{
  Metric mc;
  string scorefilepath = "/modules/Site/grantest_output.txt";
  mc.searchPath = "/data/totest/gran/";
  mc.trainGranularity(mc.searchPath,scorefilepath);
}


void Tester::TestEncoding(int argc, char* argv[])
{
	//CV_Assert(argc==2);


	cout<<getpid()<<endl;

	MTC mtc;
	mtc.ParseCfg(*(argv+1));
	//tpc.SetFootComputeRegion(boost::lexical_cast<int>(*(argv+2)));
	//tpc.SetFootComputeMethod(boost::lexical_cast<int>(*(argv+3)));
	//TPC tpc(*(argv+1));
//	string mode = *(argv+2);
//	cout<<"done step1\n";
//	//testing::InitGoogleTest(&argc, argv);
//	//RUN_ALL_TESTS();
//	if (string(*(argv+2)) == "TJPG")
//		tpc.SetCodingMode(CODING_MODE_TJPG);
//	else if (string(*(argv+2))=="JPEG")
//		tpc.SetCodingMode(CODING_MODE_JPEG)	;
//	else if (string(*(argv+2))=="PQI")
//		tpc.SetCodingMode(CODING_MODE_PQI);
//	else if (string(*(argv+2))=="TPQI")
//		tpc.SetCodingMode(CODING_MODE_TPQI);
//	else if (string(*(argv+2))=="PTP")
//		tpc.SetCodingMode(CODING_MODE_POST_TP);
//	else
//		CV_Error(CV_StsBadFlag,"Wrong coding mode");
#if NDEBUG
	mtc.SetTest(false);//test in debug
#else
	mtc.SetTest(false);
#endif
//	tpc.SetMatchingMethod(MATCHING_MSE);//=
//	tpc.SetJpegQFactor(50);//=
//	cout<<"done step 2.5\n";
//	tpc.SetInitBlockSize(Size3(32,32,1));//=
//	tpc.SetSTSIMSubWinSize(Size3(16,16,1));//=(20,20,1));//20x20 because the extended block in comparsion STSIM2 is 16+4 (boundary) or 32 +8 (boundary)
//	tpc.SetSTSIMSubWinStep(Size3(4,4,1));//=if equal to subwin size, means non_overlap, it must be a factor of block size!
//	tpc.SetQFactor(1.03);//=the threshold for 16x16 block is threshold*qfactor
//	tpc.SetSTSIMQualigyThrd(0.86);//=
//	tpc.SetInitQSize(8);//=
//	tpc.SetOverlapSizeByRatio(Size3_<double>(0.25,0.25,0)); //=
//	tpc.SetSearchStep(Size3(4,4,1));//=
//	tpc.SetCandidNum(4); //=
//	tpc.SetVarThrd(1000);//=
//	cout<<"done step 2.75\n";
//	tpc.SetPQIRectType(false);//=
//	tpc.SetLCType(POISSON_LC);//=
//	tpc.SetPBType(POST_BLENDING_ONLINE);//=
//	tpc.SetJpegQTblType(JPEG_ADAPTIVE);//=
//	tpc.SetMetricModifier(STSIM2_ADT_NAIVE);//=
//	tpc.SetSTSIM2PoolType(STSIM2_POOL_MIN);//=
//	tpc.SetPQILFBlockSize(16);//=
	cout<<"done step 3\n";
	mtc.UpdateParameters();
	cout<<"done step 4\n";
  //tpc.DebugTest();
  mtc.DummyCoding();
}

/*
void Tester::TestKmeans(void)
{
	ifstream tagfile, candfile;
	tagfile.open(".\\Apr14\\TJ_light_PB_0.900_30_32_4_1_woman_p1_TagLighting.txt");
	candfile.open(".\\Apr14\\TJ_light_PB_0.900_30_32_4_1_woman_p1_CanLighting.txt");
	char str[100];
	cv::Mat tag(282,6,CV_32F);
	cv::Mat can(282,6,CV_32F);
	for (int i=0; i < 282; i++)
		for (int j=0; j< 6; j++)
		{
			if(!tagfile.eof())
			{
				tagfile.getline(str,100,',');

				if (string(str).empty())
				{
					tagfile.getline(str,100,',');
				}
				std::istringstream s(str);
				
				s>>tag.at<float>(i,j);
				float x = tag.at<float>(i,j);
			}
			if(!candfile.eof())
			{
				candfile.getline(str,100,',');
				if(string(str).empty())
				{
					candfile.getline(str,100,',');
				}
				std::istringstream s(str);
				s>>can.at<float>(i,j);
				
			}
		}
	Mat Light = tag-can;
	Mat labels, codebook;
	cv::kmeans(Light,64,labels,TermCriteria(CV_TERMCRIT_EPS+CV_TERMCRIT_ITER,100,0.001),10,KMEANS_RANDOM_CENTERS,codebook);
	mylib::DisplayMat(codebook);
	mylib::DisplayMat(labels);
	union
	{
		float input;   // assumes sizeof(float) == sizeof(int)
		int   output;
	}data;

	for (int i=0; i< 64; i++)
		for (int j=0; j<6; j++)
		{
			data.input = codebook.at<float>(i,j);
			std::bitset<sizeof(float) * CHAR_BIT>   bits(data.output);
			string str = bits.to_string();
		}
}



void Tester::TestCovar(void)
{
  //test complex operation
  const double c = 0.001;
  Tensor<double,2> A(1,1,1);
  Tensor<double,2> B(1,1,1);
  Tensor<double,2> C(1,1,1);
  A(0,0,0) = Vec2d(3.05252,0);
  B(0,0,0) = Vec2d(-2.21077,0);
  C(0,0,0) = Vec2d(2.64,-5e-16);
  Tensor<double,2> D = C+Vec2d(c,0);
  D.Print("D");
  Tensor<double,2> E = (A*B);
  E.Print("E");
  E.Sqrt().Print("sqrt(E)");
  (D/E).Print("D/E");
}
void Tester::TestPoisson(void)
{
	cv::Mat cand = cv::imread("cbd.tif",0);
	cv::Mat mask = cv::imread("mask.tif",0);
	cv::Mat tag = cv::imread("tbd.tif",0);
	cv::Mat org = cv::imread("obd.tif",0);
	double ave = cv::mean(org)[0];
	PoissonSolver p(cand,tag,mask,org,ave);
	cv::Mat dst;
	p.poissonLightCorrection(dst,0,0);
	mylib::DisplayMat(dst,"z");
	cv::imwrite("merge.tif",dst);
}
 void Tester::TestMedian(void)
 {
	 Tensor<double,1> ref("woman_g.tiff");
	 int x = 336, y= 592;
	 Point3i offset(x,y,0);
	 Size3 overlap(8,8,0);
	 Size3 sz(16,16,1);
	 int caseNo = 3;
	 int method = 3;
	QNode<double,1> changeTo(ref,sz,offset,overlap);
	Vec<double,1> ave = changeTo.ComputePLCFoot(ref,changeTo.offset()-(changeTo.overlap()/2).Point3(),changeTo.size()+changeTo.overlap(),caseNo,method);//+changeTo.overlap());
	Vec<double,1> upave = changeTo.ComputePLCFoot(ref,changeTo.offset()-(changeTo.overlap()/2).Point3()-Point3i(changeTo.size().height,0,0),changeTo.size()+changeTo.overlap(),caseNo,method);//+changeTo.overlap());
	Vec<double,1> leftave = changeTo.ComputePLCFoot(ref,changeTo.offset()-(changeTo.overlap()/2).Point3()-Point3i(0,changeTo.size().width,0),changeTo.size()+changeTo.overlap(),caseNo,method);//+changeTo.overlap());
	

 }
 
 void Tester::TestDistortion(string fname)
 {
   vector<Tensor<double,1>> all;
   Tensor<double,1> A(fname);
   Mat a = A.GetFrame(0);
   Mat b = a.clone();
   Mat c = a.clone();
   Size win(5,5);
   long seed = 125982;
   srand(seed);
   Point center(win.height/2,win.width/2); 
   Mat rotate;
   double  phi;
   double ampl=3.;
	 for (int i=0; i<a.rows-win.height; i=i+win.height)
     for (int j=0; j<a.cols-win.width; j=j+win.width)
     {
        phi = 0.3*180*(2.*double(rand())/RAND_MAX-1.)/3.1415926;
        rotate = getRotationMatrix2D(center,phi,1);
        Mat suba = a(cv::Rect(i,j,win.height,win.width));
        Mat subb = b(cv::Rect(i,j,win.height,win.width));
        Mat subc = c(cv::Rect(i,j,win.height,win.width));
        int xx = ampl*(2.*double(rand())/RAND_MAX-1.)+0.5;
	      int yy = ampl*(2.*double(rand())/RAND_MAX-1.)+0.5;
        Mat shift = (Mat_<double>(2,3)<< 1, 0,xx, 0, 1, yy); 
        mylib::DisplayMat(shift);
        warpAffine(suba,subb,shift,suba.size(),INTER_LINEAR,IPL_BORDER_WRAP);
        warpAffine(suba,subc,rotate,suba.size(),INTER_LINEAR,IPL_BORDER_WRAP);
     }
   all.push_back(A);
   all.push_back(A);
   all.push_back(Tensor<double,1>(b));
   Tensor<double,1> B = A.MicroShift(Size3(2,2,0),3.0);
   all.push_back(B);
   all.push_back(Tensor<double,1>(c));
   Tensor<double,1> C = A.MicroRotate(Size3(2,2,0),0.3);
   all.push_back(C);
   A.DisplayAll(all,3,2,1);

 }
//
//

//TEST(TestTensor, TestPlusMinusMultiDiv)
//{
//	typedef float T;
//	Tensor<float,3> ts1(4,4,2,Vec<float,3>::all(9));
//	Tensor<float,3> ts2(4,4,2,Vec<float,3>::all(6));
//	Tensor<float,2> tc1(3,3,2,Vec<float,2>(3,2));
//	Tensor<float,2> tc2(3,3,2,Vec<float,2>(1,4));
//
//
//	//Exception Assert
//	EXPECT_NO_THROW(ts1+ts2);
//	EXPECT_NO_THROW(ts1-ts2);
//	EXPECT_NO_THROW(ts1*ts2);
//	ASSERT_NO_THROW(ts1/ts2);
//	EXPECT_NO_THROW(tc1+tc2);
//	EXPECT_NO_THROW(tc1-tc2);
//	EXPECT_NO_THROW(tc1*tc2);
//	EXPECT_NO_THROW(tc1/tc2);
//
//	//type check
//	EXPECT_EQ(ts1.type(),(ts1+ts2).type())<<"+"; 
//	EXPECT_EQ(ts1.type(),(ts1-ts2).type())<<"-"; 
//	EXPECT_EQ(ts1.type(),(ts1*ts2).type())<<"*"; 
//	EXPECT_EQ(ts1.type(),(ts1/ts2).type())<<"/"; 
//	
//	EXPECT_EQ(tc1.type(),(tc1+tc2).type());
//	EXPECT_EQ(tc1.type(),(tc1-tc2).type());
//	EXPECT_EQ(tc1.type(),(tc1*tc2).type());
//	EXPECT_EQ(tc1.type(),(tc1/tc2).type());
//	
//	///test result
//
//	EXPECT_EQ(15,(ts1+ts2)(2,2,1)[1]);
//	EXPECT_EQ(3,(ts1-ts2)(2,3,1)[2]);
//	EXPECT_EQ(54,(ts1*ts2)(3,2,1)[0]);
//	EXPECT_EQ(1.5,(ts1/ts2)(1,1,0)[1]);
//
//	
//	EXPECT_EQ(4,(tc1+tc2)(2,2,1)[0]);
//	EXPECT_EQ(6,(tc1+tc2)(2,2,1)[1]);
//	EXPECT_EQ(2,(tc1-tc2)(2,1,1)[0]);
//	EXPECT_EQ(-2,(tc1-tc2)(2,1,1)[1]);
//
//	EXPECT_EQ(-5,(tc1*tc2)(1,1,1)[0]);
//	EXPECT_EQ(14,(tc1*tc2)(1,1,1)[1]);
//	
//	EXPECT_NEAR(float(0.8462),(tc1/tc2)(1,1,0)[0],0.0001);
//	EXPECT_NEAR(float(-0.7692),(tc1/tc2)(1,1,0)[1],0.0001);
//
//	//test scalr arithmatic
//	Vec<float,3> a(2,3,6);
//	Vec<float,2> b(1,4); 
//	EXPECT_NO_THROW(ts1+a);
//	EXPECT_NO_THROW(ts1-a);
//	EXPECT_NO_THROW(ts1*a);
//	EXPECT_NO_THROW(ts1/a);
//	EXPECT_NO_THROW(tc1+b);
//	EXPECT_NO_THROW(tc1-b);
//	EXPECT_NO_THROW(tc1*b);
//	EXPECT_NO_THROW(tc1/b);
//	EXPECT_EQ(12,(ts1+a)(1,2,1)[1]);
//	EXPECT_EQ(3,(ts1-a)(0,0,0)[2]);
//	EXPECT_EQ(18,(ts1*a)(0,1,0)[0]);
//	EXPECT_EQ(4.5,(ts1/a)(1,1,1)[0]);
//	
//	EXPECT_EQ(4,(tc1+b)(2,2,1)[0]);
//	EXPECT_EQ(6,(tc1+b)(2,2,1)[1]);
//	EXPECT_EQ(2,(tc1-b)(2,1,1)[0]);
//	EXPECT_EQ(-2,(tc1-b)(2,1,1)[1]);
//
//	EXPECT_EQ(-5,(tc1*b)(1,1,1)[0]);
//	EXPECT_EQ(14,(tc1*b)(1,1,1)[1]);
//	
//	EXPECT_NEAR(float(0.8462),(tc1/b)(1,1,0)[0],0.0001);
//	EXPECT_NEAR(float(-0.7692),(tc1/b)(1,1,0)[1],0.0001);
//
//
//}
//
//
//TYPED_TEST(TestTensor, TestABS)
//{
//	if(DataType<TypeParam>::depth == 0) // uchar
//	{
//		EXPECT_EQ(0,(this->t23 - this->t13).Abs()(0,0,0)[1]);
//		EXPECT_THROW((this->t12 - this->t22).Abs(),cv::Exception);
//		EXPECT_EQ(3,this->t23.AbsDiff(v13)(0,0,0)[1]);
//		EXPECT_THROW(this->t12.AbsDiff(v22),cv::Exception);
//		EXPECT_EQ(3,this->t23.AbsDiff(t13)(0,0,0)[1]);
//		EXPECT_THROW(this->t12.AbsDiff(t22),cv::Exception);
//	}
//	else
//	{
//		EXPECT_EQ(3,(this->t23 - this->t13).Abs()(0,0,0)[1]);
//		EXPECT_EQ(TypeParam(sqrt(8.0)),(this->t12 - this->t22).Abs()(0,0,0)[0]);
//		EXPECT_EQ(3,this->t23.AbsDiff(v13)(0,0,0)[1]);
//		EXPECT_EQ(TypeParam(sqrt(8.0)), this->t12.AbsDiff(v22)(0,0,0)[0]);
//		EXPECT_EQ(3, this->t23.AbsDiff(t13)(0,0,0)[1]);
//		EXPECT_EQ(TypeParam(sqrt(8.0)),this->t12.AbsDiff(t22)(0,0,0)[0]);
//	}
//}
//
//TYPED_TEST(TestTensor, TestSumSqrtMeanRI)
//{
//	EXPECT_EQ(3,this->t12.Real()(0,0,0)[0]);
//	EXPECT_EQ(2,this->t12.Imag()(0,1,0)[0]);
//	EXPECT_EQ(64,this->t11.Sum()[0]);
//	EXPECT_EQ(4,this->t21.Mean()[0]);
//	if (DataType<TypeParam>::depth == 0) //uchar
//	{
//		EXPECT_THROW(this->t21.Sqrt()(0,0,0)[0],cv::Exception);
//	}
//	else
//	{
//		EXPECT_EQ(2,this->t21.Sqrt()(0,0,0)[0]);
//	}
//}
//
//TYPED_TEST(TestTensor, TestPowConjTrans)
//{
//
//	EXPECT_EQ(8,this->t11.Pow(3)(0,0,0)[0]);
//	if (DataType<TypeParam>::depth == 0) //uchar
//	{
//		EXPECT_THROW(this->t21.Pow(1.5)(0,0,0)[0],cv::Exception);
//		EXPECT_THROW(this->t12.Pow(2)(0,0,0)[0],cv::Exception);
//		EXPECT_THROW(this->t12.Pow(2)(0,0,0)[1],cv::Exception);
//		EXPECT_EQ(3,this->t12.Conjugate().Real()(0,0,0)[0]);
//		EXPECT_EQ(0,this->t12.Conjugate().Imag()(0,0,0)[0]);
//	
//	}
//	else
//	{
//		EXPECT_EQ(TypeParam(16.0),this->t21.Pow(2)(0,0,0)[0]);
//		EXPECT_EQ(5,this->t12.Pow(2)(0,0,0)[0]);
//		EXPECT_EQ(12, this->t12.Pow(2)(0,0,0)[1]);
//		EXPECT_EQ(3,this->t12.Conjugate().Real()(0,0,0)[0]);
//		EXPECT_EQ(-2,this->t12.Conjugate().Imag()(0,0,0)[0]);
//	}
//	
//
//
//	EXPECT_NO_THROW(t23.Transpose());
//
//}
//
//TYPED_TEST(TestTensor, TestExtendBoundary)
//{
//	Tensor<TypeParam,1> temp1;
//	Tensor<TypeParam,2> temp2;
//	temp1 = this->t11.ExtendBoundary(Size3(1,1,1), this->v21);
//	temp2 = this->t22.ExtendBoundary(Size3(2,2,0));
//	EXPECT_EQ(4,temp1(0,0,0)[0]);
//	EXPECT_EQ(2,temp1(1,1,1)[0]);
//	EXPECT_EQ(Size3(6,6,4),temp1.size());
//	EXPECT_EQ(Size3(7,7,2),temp2.size());
//	EXPECT_EQ(0,temp2(1,1,1)[0]);
//	EXPECT_EQ(0,temp2(0,0,0)[1]);
//	EXPECT_EQ(1,temp2(2,2,0)[0]);
//	EXPECT_EQ(4,temp2(2,2,1)[1]);
//	//test half extend
//	temp1 = this->t21.ExtendHalfBoundary(Size3(1,1,0),this->v11);
//	temp2 = this->t12.ExtendHalfBoundary(Size3(2,2,0));
//	EXPECT_EQ(Size3(5,5,2),temp1.size());
//	EXPECT_EQ(Size3(5,5,2),temp2.size());
//	EXPECT_EQ(this->v11,temp1(0,0,0));
//	EXPECT_EQ(this->v21,temp1(1,1,0));
//	Vec<TypeParam,2> temp(0,0);
//	EXPECT_EQ(temp,temp2(1,1,0));
//	EXPECT_EQ(this->v12,temp2(2,2,0));
//}
*/
void Tester::TestTensor(void)
{

  Tensor<double,1> ts("/home/guoxin/Projects/MTC/data/texture1.png");
        cout<<ts.size()<<endl;
  ts.Display();
        ts.Print();
  Tensor<double,2> cts = ts.ToComplex();
  double time = (double)getTickCount();
   Tensor<double,2> f;
  for (int i=0; i<1000; i++)
  {
    f = cts.DFT();
     cout<<i<<endl;
  }
  f = f.IDFT();

  time = 1000*((double)getTickCount() - time)/getTickFrequency();
  cout<<"time = " <<time<<"ms"<<endl;
  Mat temp = f.GetFrame(0);
  vector<Mat> tv;
  cv::split(temp,tv);
  Tensor<double,1>(tv[0]).Display();
  Cube c(0,0,0,10,10,1);
  cout<<c.area()<<endl;
  Steerable sp;
  sp.buildSCFpyr(f,3,4,1,false);
  vector<Tensor<double,2> >& pyr = sp.getSpaceDomainPyr();
  cv::split(pyr[12],tv);
  Tensor<double,1>(tv[0]).Display();
   Tensor<double,1> ts2("/home/guoxin/Projects/MTC/data/texture2.png");
        ts2.Print();
   cout<<ComputeMSE(ts,ts2)<<endl;

}

void Tester::TestAlgorithms(void)
{
  Tensor<double,1> ts1("/home/guoxin/Projects/MTC/data/texture1.png");
  Tensor<double,1> ts2("/home/guoxin/Projects/MTC/data/texture2.png");
        ts1.Print();
  cout<<"MSE:"<<ComputeMSE(ts1,ts2)<<endl;
  cout<<"Compare two:\n";
         auto temp = CompareElement(ts1,ts2,CMP_LE);
        temp.Print();
        temp.Display();
        cout<<"thresholding:\n";
        temp = CompareElement(ts1,Scalar(128),CMP_LE);
        temp.Print();
        temp.Display();
        cout<<"AIM:"<<ComputeAIM(ts1,ts2,50)<<endl;

}

void Tester::TestGranulateGen()
{
  Granulate g;
  g.readFiles("/home/guoxin/Projects/MTC/data/totest/","tiff");
  g.generateGrid();
}

void Tester::debugBlock(int x, int y, int sz)
{
    if (x<0||y<0||sz<0)
        return;
    cout<<"debuging ("<<x<<","<<y<<") size of "<<sz<<endl;
    string path = "./everything/";
    cout<<"loading path "<<path<<endl;
    char posstr[20];
    sprintf(posstr,"_%d_(%d_%d)",sz,x,y);
    char patt[20];
    sprintf(patt,"_%d_\\(%d_%d\\)",sz,x,y);
    Tensor<double,1> org(path+"org"+posstr+".png");
    Point3i pos(sz/4,sz/4,0);
              Size3 bsz(sz,sz,1);
              org = org.Crop(pos,bsz);
    vector<Tensor<double,1>> cands,cands_plc;
    DIR *dir;
    struct dirent *ent;
    dir = opendir(path.c_str());
    boost::regex rx("(cand"+string(patt)+")(.*)");
    boost::regex rx2("(cand_plc"+string(patt)+")(.*)");
    cout<<rx<<endl;
    cout<<rx2<<endl;
    while ((ent = readdir(dir)) != NULL)
    {
      if (!strcmp(ent->d_name,"." )) continue;
      if (!strcmp(ent->d_name,"..")) continue;
      string file_name = ent->d_name;
      cout<<file_name<<endl;
      boost::smatch res;
      boost::regex_search(file_name,res,rx);
      cout<<res[0]<<","<<res[1]<<endl;
      if (boost::regex_match(file_name,rx))
      {
          Tensor<double,1> temp = Tensor<double,1>(path+file_name);
      cout<<file_name<<endl;
          cands.push_back(temp.Crop(pos,bsz));

      }
      boost::regex_search(file_name,res,rx2);
            cout<<res[0]<<","<<res[1]<<endl;
            if (boost::regex_match(file_name,rx2))
            {
                 cout<<file_name<<endl;
                Tensor<double,1> temp = Tensor<double,1>(path+file_name);
                         cands_plc.push_back(temp.Crop(pos,bsz));

            }

    }

    closedir(dir);
    for (unsigned int i=0; i<cands.size();i++)
    {
        cout<<"Computing STSIM-2 for candidate "<<i<<endl;
        metric::ComputeSTSIM2(org,cands[i],Size3(16,16,1),Size3(16,16,1),3,4,false,FilterBoundary::FILTER_BOUND_FULL,FeaturePoolType::FEATURE_POOL_MIN,MetricModifier::STSIM2_NEW_L1,true);
    }


}

#define  _CRT_SECURE_NO_WARNINGS
//#include "QGrid.h"
#include "Tester.h"
#include "Steerable.h"
//#include "TPC.h"

using namespace cv;
using namespace std;
int main(int argc, char* argv[])
{
	//string s="ATS\n012\1\\";
	//int i =strlen("ATS\n012\1\\");
	Tester test;
  cout<<CV_MAJOR_VERSION<<endl;
  cout<<CV_MINOR_VERSION<<endl;
	//Tensor<double,1> ts;
	//ts.Load("woman_g.tiff");
	//ts(Cube(0,0,0,16,16,1)).Print();
	//ts.Display();
	//Mat_<cv::Complex<double>> a(Size(1,1));
	//a.at<cv::Complex<double>>(0,0)=cv::Complex<double>(0.2,0.8);
	////mylib::DisplayMat(a);
	//Mat_<cv::Complex<double>> b(Size(1,1));
	//b.at<cv::Complex<double>>(0,0)=cv::Complex<double>(0.4,0.5);
	////Mat c= a/b;
	////mylib::DisplayMat(c);
	//Tensor_<double,2> A(a);
	//Tensor_<double,2> B(b);
	//Tensor_<double,2> C = A/B;

	//C.Abs().Print();
 	//test.TestSSIM("BBB_gray_p1.tif","BBB_gray_p2.tif");
	//vector<string> filenames;
	//filenames.push_back("H:\\MyStudy\\Research Reports\\Mar11\\TJ_0.93_4_seam.tif");
	//filenames.push_back("H:\\MyStudy\\Research Reports\\Mar11\\TJ_0.93_8_seam.tif");
	//filenames.push_back("H:\\MyStudy\\Research Reports\\Mar11\\TJ_0.93_16_seam.tif");
	////filenames.push_back("H:\\MyStudy\\Research Reports\\Mar11\\TJ_0.93_32_nolighting.tif");
	//////filenames.push_back(".//Mar6//TJ_0.9300_10.4_32_4_1_baboon.tif");
	//////filenames.push_back(".//Mar6//J_0.9300_10.4_32_4_1_baboon.tif");
	//mylib::CombineImage(filenames,string("H:\\MyStudy\\Research Reports\\Mar11\\combine_seam.tif"));
	//test.TestBoundary();
	//test.TestTensor();
	//test.TestFilter2D();
	//test.TestAlgorithms();
	//test.TestSSIM("/home/guoxin/Projects/MTC/data/texture1.png","/home/guoxin/Projects/MTC/data/texture2.png");
	test.TestGranulate();
	//test.TestKmeans();
	//test.TestHuffTree();
	//test.TestJPEG();
	//test.Test3DMat();
	//test.TestPoisson();
	//test.TestMedian();
  //test.TestCovar();
  //test.TestSTSIMBorder();
  //test.TestVarLenFunction();
 //	test.TestMatching();
  //test.TestDistortion("..\\images\\test.tiff");
  //test.TestQuilting("..\\images\\training\\8064_1_g.tiff","..\\images\\training\\42-17730995_0_d_g.tiff");
  //test.TestEmbedding("..\\images\\training\\8064_1_g.tiff","..\\images\\training\\42-17730995_0_d_g.tiff");
  
  //test.TestEncoding(argc,argv); 
	//testing::InitGoogleTest(&argc, argv);
	//RUN_ALL_TESTS();
	//test.TestTPC();
	//test.TestLighting();
	//test.TestRoi("baboon.pgm"); 
	//test.TestInterp("baboon.pgm");
	//test.TestQuilting("baboon.pgm");
	//test.TestInitGrid();
	//test.TestBoundMatch();
	return 0;

}

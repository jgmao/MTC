#pragma once
#define _CRT_SECURE_NO_WARNINGS
//#include "QGrid.h"

#include "MyLib.h"
#include "MTC.h"
//#include "CubePlus.h"
//#include "HuffTree.h"
//#include <gtest\gtest.h>
#include "TensorLite.h"
#include "algorithms.h"
#include <Granulate.h>
#include <opencv2/opencv.hpp>
#include <Metric.h>
using namespace std;
using namespace tensor;
using namespace metric;
using namespace mtc;
class Tester
{
public:
	Tester(void);
	~Tester(void);
	//void TestInterp(string filename);
	//void TestRoi(string filename);
	//void TestQuilting(string filename);
  	//void TestQuilting(string f1, string f2);
  	//void TestEmbedding(string f1,string f2);
	//void TestInitGrid(void);
	//void TestBoundMatch(void);
	//void TestTPC(void);
	void TestSSIM(string file1, string file2);
	//void TestJPEG(void);
	//void TestLighting(void);
	//void TestBoundary(void);
	//void TestMedian(void);
	void TestFilter2D(void);
	//void TestKmeans(void);
	//void TestHuffTree(void);
	void TestEncoding(int argc, char* argv[]);
	//void Test3DMat(void);
	//void TestPoisson(void);
  	//void TestCovar(void);
  //void TestSTSIMBorder(void);
  //void TestVarLenFunction(void);
  //void TestMatching(void);
  //void TestDistortion(string fname);
  void TestTensor();
  void TestAlgorithms();
  void TestGranulateGen();
  void TestMetric();
  void StudyMetricFeatures();
  void TestGranulateTrain();
  void debugBlock(int x, int y, int sz);
  void TestSpeed(void);

};/*

template <class T> class TestTensor : public testing::Test
{
public:
	Tensor<T,1> t11;
	Tensor<T,2> t12;
	Tensor<T,3> t13;
	Tensor<T,1> t21;
	Tensor<T,2> t22;
	Tensor<T,3> t23;
	Vec<T,1> v11, v21;
	Vec<T,2> v12, v22;
	Vec<T,3> v13, v23;
	TestTensor()
	{
		v11 = Vec<T,1>::all(2);
		v21 = Vec<T,1>::all(4);
		v12 = Vec<T,2>(3,2);
		v22 = Vec<T,2>(1,4);
		v13 = Vec<T,3>::all(9);
		v23 = Vec<T,3>::all(6);
		t11 = Tensor<T,1>(4,4,2,v11);
		t21 = Tensor<T,1>(4,4,2,v21);
		t12 = Tensor<T,2>(3,3,2,v12);
		t22 = Tensor<T,2>(3,3,2,v22);
		t13 = Tensor<T,3>(4,4,2,v13);
		t23 = Tensor<T,3>(4,4,2,v23);
	}
	~TestTensor()
	{
	}
};*/

//typedef testing::Types<float, double,UINT8> TestTypes;
//TYPED_TEST_CASE(TestTensor,TestTypes);

#pragma once
#define  _CRT_SECURE_NO_WARNINGS
#ifndef MTC_H
#define MTC_H
//#define _USE_DATABASE
#include<thread>
#include<cstdint>
#ifdef _USE_DATABASE
#include <SQLAPI.h>
#define DB_HOST "localhost\\OMEGA@MTC"
#define DB_USR "jgmao"
#define DB_PWD "040806"
#define DB_NAME " Test_MTC_Tbl "
#endif 

#include "QGrid.h"
#include "Statistics.h"
//#include <process.h>
#include <boost/lexical_cast.hpp>
//#include <boost/filesystem.hpp>
#include "PostBlendSet.h"
//#include <gtest/gtest.h>
#include <boost/dynamic_bitset.hpp>
#include "CoderDef.h"
#include "NTest.h"
#include <sys/stat.h>
#include <TensorHelper.h>
#include "HuffTree.h"
using namespace std;
namespace mtc{

class MTC :public QGrid<double,1>
{

public:
	typedef double T;
  //Metric mc;
	static const int cn = 1;
	EXPORTLIB MTC(void);
	EXPORTLIB MTC(const string& imageName, CodingMode codemode=CodingMode::CODING_MODE_MTC);
	EXPORTLIB MTC(const string& imageName, double qualityThrd , const Size3& initSize = Size3(16,16,1), Size3 searchStep = Size3(4,4,1), int candidNum = 4, double initQSize = 8, const Size3_<double>& overlapRatio = Size3_<double>(0.25,0.25,0));
	EXPORTLIB void Coding(CodingMode codemode); //uniform coding interface, manage file names and coding modes;
	EXPORTLIB void Coding(void);//coding use preassigned code mode
	EXPORTLIB void Encode(void); //generate bit streams
	EXPORTLIB void Decode(ofstream& fh);//decode image
	EXPORTLIB ~MTC(void);
	EXPORTLIB void DummyCoding(void);
	EXPORTLIB void Init(const string& imageName, CodingMode codemode=CodingMode::CODING_MODE_MTC);
	Statistics stat;
	//////////// MTC parameters set
	EXPORTLIB void SetInitBlockSize(const Size3& initBlockSize);
	//EXPORTLIB void SetPredictQualityThrd(double qualityThrd);
	EXPORTLIB void SetSearchStep(const Size3& step);//search setp of matching
	EXPORTLIB void SetSTSIMSubWinSize(const Size3& subWinSize);//=
	EXPORTLIB void SetSTSIMSubWinStep(const Size3& subWinStep);//=
	EXPORTLIB void SetCandidNum(unsigned int candidNum);
	EXPORTLIB void SetQFactor(double qFactor);//quality test/quantize step incresing ratio
	EXPORTLIB void SetSTSIMQualigyThrd(double qThrd);
  EXPORTLIB void SetMSEThrd(double mseThrd);
        EXPORTLIB void SetInitQSize(double qSize);
        EXPORTLIB void SetCodingMode(CodingMode qmode);
        EXPORTLIB void SetMetricModifier(MetricModifier modifier);
        EXPORTLIB void SetOverlapSizeByRatio(const Size3_<double>& ratio);
        EXPORTLIB void SetOverlapSize(const Size3& overlap);
        EXPORTLIB void SetPQIRectType(bool rectType);
        EXPORTLIB void SetPBType(PostBlendingType pbType); // set postblending type
        EXPORTLIB void SetLCType(LightingCorrectionType lcType); // set Lighting correction type
        EXPORTLIB void SetJpegQFactor(int q); //between 0 - 100
        EXPORTLIB void SetJpegQTblType(JPEGType quanType); // set the quantizer table type
        EXPORTLIB void UpdateParameters(void); //must call this after change any parameters.
        EXPORTLIB void SetTest(bool test);
        EXPORTLIB void SetMatchingMethod(MatchingMethod m_method);
        EXPORTLIB void SetSTSIM2PoolType(FeaturePoolType pool_type);
        EXPORTLIB void SetPQILFBlockSize(int pqiSize);
        EXPORTLIB void ParseCfg(const string& cfgname);
        EXPORTLIB void SetBlendingMethod(BlendingMethod blend_method=BlendingMethod::SHORTEST_PATH_BLENDING);
        EXPORTLIB void SetFootComputeRegion(int region=0);
        EXPORTLIB void SetFootComputeMethod(int method=0);
        EXPORTLIB void DebugTest(void);
        EXPORTLIB Size3 GetSubWinSize() const;
        EXPORTLIB Size3 GetSubWinStep() const;
  protected:
	//i/o parameters: include path, filename and storage modifier
	string prefix; //prefix of output
	string path; //output path
	string baseName;//basename (without .ext)
	fstream outputfile,cbfile;
	fstream logfile;
  fstream candPosLog;
        string count_exp;
	int footComputeRegion;
	int footComputeMethod;
#ifdef _USE_DATABASE
	SAConnection conn;//Database connector 
	SACommand cmd;//Database command executor
#endif
	//setup parameters: need to be initialized before coding
	double qualityThrd;
  double mseThrd;
	double orgQualityThrd; ///backup
	double qfactor;
	int max_bits;
	int levels;
	bool rectType; //false is rect split type 2, true is tpype 1 (type 1 has 1x1 blocks)
	CodingMode mode;//coding mode
	JPEGType jpegType;
	MetricModifier metricModifier;
	PostBlendingType postBlendType;
	LightingCorrectionType lightCorrectionType;
	bool unitTest;
	FeaturePoolType stsim2PoolType;
	BlendingMethod blend_method;
  int acount;
	//other finer defined parameters: 
	//Size3 stsimSubWinSize;
	//Size3 stsimSubWinStep;
	//Size3 stsimSubWinSizeRatio;
	//Size3 stsimSubWinStepRatio;
  Mat lse_weight;
  Size3 subSize;
  Size3 subStep;
	int jpegQNum;
	Size3 initSize;
	double initQSize;
	vector<double> qstep;
	vector<double> thrd;
	MatchingMethod matching_method;
	//intermediat or final data structures
	Tensor<T,cn> dctEnsemble;
	Tensor<T,cn> edgeEnsemble;
	Tensor<T,cn> LFImage;
	vector<vector<int>> edgeMap;//0: T/T, 1: S/S, 2: T, 3: T/S, 4: S, 5: T/T/S
	Tensor<T,3> rst_with_seam;
	Tensor<T,1> pred_rst;
  Tensor<T,1> pred_afterLC;
  Tensor<T,1> pred_postLC;
  Tensor<T,1> plc_rst;
	Tensor<T,1> footmap;
	vector<uint16_t> JpegQuanTbl,JpegQuanTbl_alt;
	vector<vector<uint8_t>> JpegHuffTbl_DC,JpegHuffTbl_AC;
	int totalBits;
	PBSet pbSet;
	//use for record the block offset and size to be used for post blending
	vector<vector<Vec<T,cn>>> lightTagVec;
	vector<vector<Vec<T,cn>>> lightCanVec;
	Tensor<T,cn> VQCodeBook;
	vector<int> VQLabels;
	vector<Tensor<T,cn>> lightPlanes;
	//unit test
	NTest<string> nts;
	NTest<int> nti;
	NTest<T> nt;
	int nPQILFSize;
	Tensor<T,cn> SAT;
	bool acceptDirect;
  cv::Mat iMahaCovar;
  //list<FootItem> FootTable;
  std::mutex foot_lock;
	std::map<Point3i,FootItem,ComparePoint3i> FootTable;
  int fLevel;
  struct Option
	{
		string name;
		void (MTC::*pfuncInt)(int);
		void (MTC::*pfuncSize)(const Size3&);
		void (MTC::*pfuncPt)(const Point3i&);
	}option;

public:
	std::map<int,string> Huffman_DC_Map, Huffman_AC_Map;
	HuffTree Huffman_AC_Tree, Huffman_DC_Tree;
   public:
    Tensor<T,cn>& PQICodingLFComponent(void);

private:
	void CodingMTC(QTree<T,cn> &qNode, int qLevel=0); //this is online coding
	bool IsAcceptPredict(QNode<T,cn>& qNode,  CompareCriteria criteria, double param1, double param2 = 4);
	int IsAcceptPredict(const vector<Point3i>& matchCandid,  QTree<T,cn>& qNode,MetricModifier metricModifier, int level=-1);//! 20130913 level -1 means no start with 0 foot , no LC
	//=COMPARE_CRITERIA_SSIM, double param1=3, double param2 = 4);
	void TreatAcceptedInterp(QTree<T,cn> &qNode, int qLevel);
	int TreatRectSplit(QTree<T,cn>& qNode, int qLevel);
	void PostBlending(const vector<TPrecord>& tprecs);
	void PostBlending(const PBRecord& rec, CubePlus& roi);
  Tensor<T,cn> PostPLC(Tensor<T,cn>& changeFrom,Tensor<T,cn>& changeTo);
	bool TexturePrediction(QTree<T,cn>& qNode, int qLevel=0); //routine do TP and modify results
  bool AdaptiveTPSS(const QTree<T,cn>& qNode, int qLevel=0); 
	void PostCoding(QTree<T,cn> &qNode, int qLevel=0); //this is offline coding
	void CodingJPEG(QTree<T,cn>& qNode,int qLevel=0); //the size of qNode must be 8x8
	//bool QuadTreeIterpCoding(QNode<T,cn>& qNode, int qLevel=0);
	void InitJpegParam(int quality=50);
	string ComputeJpegHuff(const Vec<T,cn>& val, std::map<int,string>& huffMap);
	string ComputeJpegHuff(const Vec<T,cn>& val, uint8_t runlength, std::map<int,string>& huffMap);
	void UpdatePBSet(const QNode<T,cn> &qNode, Point3i matchPos);
	void UpdatePBSet(const QNode<T,cn> &qNode, const QNode<T,cn>& candid);
	//void SaveLighting(void);
	void CollectLighing(void);
  void CollectBits(void);
	void UpdateLog(void);
	void BuildVQCodebook(vector<Tensor<T,cn>>& lightPlanes,Tensor<T,cn>& coeffs, int bits=6);
	//int SearchCodeword(vector<Vec<T,cn>>& poly);

	Tensor<T,cn>& LoadPreDefLighting(void);
	QTree<T,cn>& ComputeJPEG(QTree<T,cn>& qNode);
	QTree<T,cn>& CodingCore(QTree<T,cn>& qTree);
	Tensor<T,cn>& ComputeEnsembleDCT(int bsize);
	void ModifyMetric(const Point3i& blockPos);
	void UpdateSAT(const Point3i& sPos,const Size3& sz);
	QNode<T,cn> GetValidCandid(Point3i matchPos, const QNode<T,cn>& qNode);
  QNode<T,cn> GetValidCandidSimple(Point3i matchPos, const QNode<T,cn>& qNode);
  QNode<T,cn> GetValidNode(Point3i matchPos, const QNode<T,cn>& qNode, const Tensor<T,cn>& ref);
	void ScanFile(istream& in);
	void ScanLine(string& line);
  vector<pair<cv::Point3i, FootItem>>  RetrieveFeet(QNode<double,1>& qNode,int level=0);
  auto RetrieveFoot(QNode<double,1>& qNode, int level = 0) -> decltype( MTC::FootTable.begin()); //the return value indiate which foot is found
  auto RetrieveFootsNew(QNode<double,1>& qNode, int level = 0) -> decltype( MTC::FootTable.begin()); //the return value indiate which foot is found
  void UpdateFoots(const QNode<double,1>& org, const QNode<double,1>& cand);
  void UpdateFeet(const vector<pair<Point3i, FootItem>>& feet);
  void RemoveFeet(const vector<pair<Point3i, FootItem>>& feet);
  CompareCriteria ParseCriteria(MetricModifier metricModifier);
  void InitExt(void);
private:
	void InitHuffCoeff(void);
};
}
#endif

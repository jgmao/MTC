#pragma once
#define  _CRT_SECURE_NO_WARNINGS
#include<vector>
#include<string>
//#include<cv.h>
#include <opencv2/opencv.hpp>
#include "Size3.h"
#include <future>
#include "Cube.h"
namespace mtc{
using namespace std;
using namespace tensor;
class TPrecord
{
public:
	cv::Point3i matchPos; //raw position, do not include boundary 
	Cube targetROI;
	Size3 overlapSize;
	TPrecord(const cv::Point3i& match, const Cube& target, const Size3& ovSize)
	{
		this->matchPos = match;
		this->targetROI = target;
		this->overlapSize = ovSize;
	}
};

class Statistics
{
public:
	Statistics(void);
	Statistics(int levels, int init_size);
	~Statistics(void);
public:
	vector<string> tempstr;
	vector<string> footbits;
	vector<string> treebits;
	vector<cv::Point3i> footPos;
	vector<Size3> blockSize;
	vector<int> numTestInLevel;
	vector<int> numBlocksInLevel;
	vector<int> footBitInLevel;
	vector<int> treeBitInLevel;
	vector<double> psnrInLevel;
	vector<double> bppInLevel;
	vector<int> pixInLevel;
	string outputBit;
	///use for texture prediction
	vector<int> numTextureInLevel;
	vector<double> psnrTextureInLevel;
	vector<int> textureBitInLevel;
	vector<TPrecord> tprecs;
	//map<int, int> numTestInLevel;
	//map<int, int> numBlocksInLevel;
	int levels;
	int init_size;
	int footbitLength;
	int treebitLength;
	int baseCoderBitLength;
  int adaptiveJpegBitLength;
  int plcBitLength;
	double PSNR;
	double bpp;
	double SSIM;
	int lightingCodeLength32;
	int lightingCodeLength16;
	int tpBlockNum32;
	int tpBlockNum16;
	vector<int> tpBlockNum;
	vector<int> tpBits;
	vector<int> jBlockNum;
	vector<int> jBits;
	int beginH,endH,beginM,endM,beginS,endS;
	double duration;
	time_t beginT,endT;
	//for record of lighting correction
	vector<Cube> lightROI;
  double lightingDiff;
public:
	EXPORTLIB int getTotalLength(void);
	EXPORTLIB int getTotalBlocks(void);
	EXPORTLIB double getPercentage(int qlevel);
	EXPORTLIB TPrecord& searchTPRecord();
};

}

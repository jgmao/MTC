#define  _CRT_SECURE_NO_WARNINGS
#include "Statistics.h"
namespace mtc{
//#define MY_MAKELEVEL(max_size, cur_size) log(double(max_size/cur_size))/log(2.0)
//#define MY_TEXTLV(max_size) log(double(max_size / 8))/log(2.0)
Statistics::Statistics(void)
{
	init_size=16;
	levels = 5;
	footbitLength=0;
	treebitLength=0;
	baseCoderBitLength = 0;
	PSNR = 0.0;
	numTestInLevel = vector<int>(levels);
	numBlocksInLevel = vector<int>(levels);
	SSIM =0;
	//blockSize.clear();
	lightingCodeLength32= 0;
	tpBlockNum32 = 0;
	lightingCodeLength32 = 0;
	tpBlockNum32 = 0;
  adaptiveJpegBitLength = 0;
  plcBitLength =0;
}


Statistics::Statistics(int levels, int init_size)
{
	this->init_size= init_size;
	this->levels =  levels;
	footbitLength=0;
	treebitLength=0;
	baseCoderBitLength = 0;
	PSNR = 0.0;
	bpp = 0.0;
	SSIM = 0;
	numTestInLevel = vector<int>(levels);
	numBlocksInLevel = vector<int>(levels);
	footBitInLevel = vector<int>(levels);
	treeBitInLevel = vector<int>(levels);
	psnrInLevel = vector<double>(levels);
 	bppInLevel = vector<double>(levels);
	pixInLevel = vector<int>(levels);
	numTextureInLevel = vector<int>(levels);
	psnrTextureInLevel = vector<double>(levels);
	textureBitInLevel = vector<int>(levels);
	lightingCodeLength32= 0;
	tpBlockNum32 = 0;
	lightingCodeLength16 = 0;
	tpBlockNum16 = 0;
	for (int i=0; i<levels -2;i++)
		pixInLevel[i] = (init_size>>i)*(init_size>>i);
	pixInLevel[levels-2] = 2;
	pixInLevel[levels-1] = 1;
	//blockSize.clear();
	outputBit = "";
  adaptiveJpegBitLength = 0;
  plcBitLength = 0;
  lightingDiff=0;
}


Statistics::~Statistics(void)
{
}

int Statistics::getTotalLength(void)
{
	return footbitLength+treebitLength;
}

int Statistics::getTotalBlocks(void)
{
	return blockSize.size();
}

double Statistics::getPercentage(int qLevel)
{
	return  100.0*(double)numBlocksInLevel[qLevel]/(double)numTestInLevel[qLevel];
}
}

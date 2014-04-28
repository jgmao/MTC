
#include "Metric.h"
#include <string>
#include <iostream>
#include <Steerable.h>
#include <regex>
using namespace cv;
using namespace std;
using namespace metric;
int main(void)
{
  Metric mc;


  mc.subwinSize = Size3(16,16,1);
  mc.subwinStep = Size3(16,16,1);
  //mc.subsample = true;
  //mc.changeWin=true;
  //mc.searchPath =  "H:/Code/Matlab/metric trainning/data/";
  string scorefilepath = "./modules/Site/subtestoutput_coding64_inter.txt";
  mc.searchPath = "./data/totest/dist_coding_64_inter/";
  //mc.trainSTSIM2Weights(mc.searchPath,scorefilepath);
  mc.trainMetirc(mc.searchPath,scorefilepath);
  return 0;
}

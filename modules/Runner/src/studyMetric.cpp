
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
  string logfilepath = "./everything/everything.txt";
  mc.searchPath = "./everything/";
  //mc.trainSTSIM2Weights(mc.searchPath,scorefilepath);
  mc.studyMetricFeature(logfilepath);
  return 0;
}


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
  string scorefilepath = "/modules/Site/grantest_output.txt";
  mc.searchPath = "/data/totest/gran/";
  mc.subsample = true;
  mc.changeWin = false;
  mc.trainGranularity(mc.searchPath,scorefilepath);

  return 0;
}



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
  //64x64
  mc.subwinSize = Size3(64,64,1);
  mc.subwinStep = Size3(64,64,1);
  //mc.subsample = true;
  //mc.changeWin=true;
  mc.searchPath =  "/home/guoxin/Projects/MTC/data/coding/coding64/D1/";
  mc.searchExt = ".tif";
  mc.readFiles(mc.searchPath,"([a-zA-Z0-9_-]+)_(\\d+).tif","");
  mc.searchPath = "/home/guoxin/Projects/MTC/data/coding/coding64/D2/";
  mc.readFiles(mc.searchPath,"([a-zA-Z0-9_-]+)_(\\d+).tif","");
  mc.searchPath = "/home/guoxin/Projects/MTC/data/curet/curet64/";
  mc.readFiles(mc.searchPath,"([a-zA-Z0-9_-]+)__(\\d+).tif","");
  mc.searchPath = "/home/guoxin/Projects/MTC/data/corbis/p1_64/";
  mc.readFiles(mc.searchPath,"([a-zA-Z0-9_-]+)_(\\d+).tif","");
  mc.computeFeatures();


  return 0;
}

#include "LightingVerify.h"
//#include "PoissonSolver.h"
namespace lighting{
LightingVerify::LightingVerify(void)
{
}


LightingVerify::~LightingVerify(void)
{
}

void LightingVerify::readfile(string filepath, string ext)
{
  
  this->fpath = filepath;
  Metric mc;
  this->fnames = mc.readFiles(filepath,"","",ext);

}

vector<Tensor<double,1>>& LightingVerify::estimateLight(vector<Tensor<double,1>>& images, vector<Tensor<double,1>>& rst, int blksize, int bdsize, int fsize)
{
  for (Tensor<double,1>& A : images)
  {
    Tensor<double,1> temp= lt.ComputeTPSS(A,0.001);
    rst.push_back(temp.Clone());
  }
  return rst;
}

void LightingVerify::estimateLight(Tensor<double,1>& im, Tensor<double,1>& rst, int blksize, int bdsize, int fsize)
{
  rst = lt.ComputeTPSS(im,0.001);
}
vector<double>& LightingVerify::doPLC(int blksize, int bdsize, int fsize )
{
  this->blksize = blksize;
  this->bdsize = bdsize;
  this->fsize=fsize;
  Mat mask = Mat::ones(Size(blksize,blksize), CV_64F)*255;
  Mat extmask;
  //cv::copyMakeBorder(mask,extmask,bdsize,bdsize,bdsize,bdsize,BORDER_CONSTANT,0);
  cv::copyMakeBorder(mask,extmask,bdsize,bdsize,bdsize,bdsize,BORDER_CONSTANT,0);
  Vec<double,1> upfoot,leftfoot,foot;
  this->orignalWithBorder.clear();
  this->MSEs.clear();
  this->originalImages.clear();
  this->correctedImages.clear();
  this->lightingOrg.clear();
  this->lightingRst.clear();
  this->feet.clear();
  for (auto& fname : this->fnames)
  {
    Tensor<double,1> A(fpath+fname);
    this->orignalWithBorder.push_back(A.Clone());
    Tensor<double,1> blk = A.Crop(Point3i(bdsize,bdsize,0),Size3(blksize,blksize,1));
    this->originalImages.push_back(blk);
    foot = A(Cube(blksize+bdsize-fsize/2,blksize+bdsize-fsize/2,0,fsize+1,fsize+1,1)).Mean();
    upfoot = A(Cube(bdsize-fsize/2,blksize+bdsize-fsize/2,0,fsize+1,fsize+1,1)).Mean();
    leftfoot = A(Cube(blksize+bdsize-fsize/2,bdsize-fsize/2,0,fsize+1,fsize+1,1)).Mean();
    feet.push_back(FootType(upfoot[0],0,blksize+1));
    feet.push_back(FootType(leftfoot[0],blksize+1,0));
    feet.push_back(FootType(foot[0],blksize+1,blksize+1));
    PoissonSolver pb(A.GetFrame(0),A.GetFrame(0),extmask,A.GetFrame(0),feet);
    Mat rst;
    pb.poissonLightCorrection(rst);
    correctedImages.push_back(rst(Rect(bdsize,bdsize,blksize,blksize)));
  }

  estimateLight(originalImages,lightingOrg,blksize,bdsize,fsize);
  estimateLight(correctedImages,lightingRst,blksize,bdsize,fsize); 

  for (int i=0; i< fnames.size(); i++)
  {
    double temp = ComputeMSE(lightingOrg[i],lightingRst[i]);
    MSEs.push_back(temp);///double(lightingOrg[i].size().area())
  }
  return MSEs;
}

vector<double>& LightingVerify::doAdaptivePLC(int blksize, int bdsize, int fsize,double thrd )
{
  this->blksize = blksize;
  this->bdsize = bdsize;
  this->fsize=fsize;
  Mat mask = Mat::ones(Size(blksize,blksize), CV_64F)*255;
  Mat extmask;
  cv::copyMakeBorder(mask,extmask,bdsize,bdsize,bdsize,bdsize,BORDER_CONSTANT,0);
  Vec<double,1> upfoot,leftfoot,foot;
  int count=0;
  this->orignalWithBorder.clear();
  this->originalImages.clear();
  this->correctedImages.clear();
  this->lightingOrg.clear();
  this->lightingRst.clear();
  this->MSEs.clear();
  this->feet.clear();
  for (auto& fname : this->fnames)
  {
    Tensor<double,1> A(fpath+fname);
    this->orignalWithBorder.push_back(A.Clone());
    Tensor<double,1> blk = A.Crop(Point3i(bdsize,bdsize,0),Size3(blksize,blksize,1));
    this->originalImages.push_back(blk);
    foot = A(Cube(blksize+bdsize-fsize/2,blksize+bdsize-fsize/2,0,fsize+1,fsize+1,1)).Mean();
    upfoot = A(Cube(bdsize,blksize+bdsize-fsize/2,0,fsize+1,fsize+1,1)).Mean();
    leftfoot = A(Cube(blksize+bdsize-fsize/2,bdsize,0,fsize+1,fsize+1,1)).Mean();
    feet.push_back(FootType(upfoot[0],0,blksize+1));
    feet.push_back(FootType(leftfoot[0],blksize+1,0));
    feet.push_back(FootType(foot[0],blksize+1,blksize+1));
    PoissonSolver pb(A.GetFrame(0),A.GetFrame(0),extmask,A.GetFrame(0),feet);
    Mat rst;
    pb.poissonLightCorrection(rst);
    correctedImages.push_back(rst(Rect(bdsize,bdsize,blksize,blksize)));
    lightingOrg.push_back(originalImages[count]);
    lightingRst.push_back(correctedImages[count]);
    estimateLight(originalImages[count],lightingOrg[count],blksize,bdsize,fsize);
    estimateLight(correctedImages[count],lightingRst[count],blksize,bdsize,fsize);

    double mse = ComputeMSE(lightingOrg[count],lightingRst[count]);///double(lightingOrg[count].size().area());
    if (mse > thrd)
    {
      pb.addFoot(FootType(A(Cube(bdsize+blksize/2-fsize/2.0,blksize+bdsize-fsize/2,0,fsize+1,fsize+1,1)).Mean()[0],blksize/2+1,blksize+1));//right center 
      pb.addFoot(FootType(A(Cube(blksize+bdsize-fsize/2,bdsize+blksize/2-fsize/2.0,0,fsize+1,fsize+1,1)).Mean()[0],blksize+1,blksize/2+1));//bottom center
      pb.poissonLightCorrection(rst);
      correctedImages[count]=(rst(Rect(bdsize,bdsize,blksize,blksize)));
      //Tensor<double,1> templt = lightingRst[count].Clone();
      estimateLight(correctedImages[count],lightingRst[count],blksize,bdsize,fsize);//lighting rst seems doesn't update
      //double tempmse = ComputeMSE(templt, lightingRst[count]);
      mse = ComputeMSE(lightingOrg[count],lightingRst[count]);///double(lightingOrg[count].size().area());
    }
    MSEs.push_back(mse);
    count++;
  }
  return MSEs;
}
void LightingVerify::saveResults(string folder)
{
  string path = fpath+folder;
  ofstream fp;
# ifdef _WIN32
  _mkdir(path.c_str());
# else
  mkdir(path.c_str(),0755);
# endif
  fp.open(path+"mse.txt",ios::out);
  for (int i=0; i< fnames.size();i++)
  {
    vector<Tensor<double,1>> dummy;
    Tensor<double,1> temp = this->originalImages[i];    
    dummy.push_back(temp);
    dummy.push_back(lightingOrg[i]);//.ExtendBoundary(Size3(bdsize,bdsize,0),255));
    Tensor<double,1> x;
    x.DisplayAll(dummy,2,1,true,path+"org_"+fnames[i]);
    dummy.clear();
    //temp.SetBlock(Point3i(bdsize,bdsize,0),correctedImages[i]);
    temp=correctedImages[i];
    dummy.push_back(temp);
    dummy.push_back(lightingRst[i]);//.ExtendBoundary(Size3(bdsize,bdsize,0),255));
    x.DisplayAll(dummy,2,1,true,path+"test_"+fnames[i]);
    fp<<i<<":"<<MSEs[i]<<endl;
  }
  fp.close();
}

}

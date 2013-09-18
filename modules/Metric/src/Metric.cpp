#include "Metric.h"
#include <dirent.h>
#include <sys/stat.h>
namespace metric{


Metric::Metric()
{
  featureNum = 4;
  subwinSize  = Size3(16,16,1);
  subwinStep = Size3(16,16,1);
  Nor = 4;
  Nsc = 3;
  gamma0 = vector<vector<Mat_<double> > >(Nor*Nsc+2,vector<Mat_<double>>(featureNum));
  gamma1 = vector<vector<Mat_<double> > >(Nor*Nsc+2,vector<Mat_<double>>(featureNum));
  lambda1 = vector<vector<double> >(Nor*Nsc+2,vector<double>(featureNum));
  lambda0 = vector<vector<double> >(Nor*Nsc+2,vector<double>(featureNum));
  subsample = false;
  changeWin = false;
  sameCount = 0;
  diffCount = 0;
  coeffNum = 0;
  coeffcount0=0;
  coeffcount1=0;
}

bool Metric::readFiles(string path, string searchPattern, string searchExt)
{
  //check if path end with /
  char endwith = *path.rbegin();
  if (endwith!='/')
    path=path+"/";
  this->clustercount=0;
  string filename;
#ifdef WIN32
  string fullSearchPath = path + searchPattern+"*"+searchExt;
  WIN32_FIND_DATA FindData;
  HANDLE hFind;
  hFind = FindFirstFile( fullSearchPath.c_str(), &FindData );
  if( hFind == INVALID_HANDLE_VALUE )
  {
    cout << "Error searching directory\n";
    return false;
  }
#else
  DIR *dir;
  struct dirent *ent;
  dir = opendir(path.c_str());
#endif
  do
  {
#ifdef WIN32
      filename = FindData.cFileName;
      filenames.push_back(FindData.cFileName);
#else
      if (!strcmp(ent->d_name,"." )) continue;
      if (!strcmp(ent->d_name,"..")) continue;
      string file_name = ent->d_name;
      unsigned long found = file_name.rfind("."+searchExt);
      if (found!=std::string::npos)
      {
          filename = path + "/" + file_name;
          filenames.push_back(filename);
          cout<<"Got: "<<file_name<<endl;
      }
#endif
    std::smatch res;
    std::regex rx(searchPattern+"_([^]*)_(\\d+)"+searchExt);
    std::regex_search(filename,res,rx);
    std::cout << "class: "<<res[1] << " , num: " << res[2] << "\n";
    if (clusternames.find(res[1])==clusternames.end())
    {
      clusternames.insert(pair<string,int>(res[1],clustercount));
      clustercount++;
    }
  }
#ifdef WIN32
  while( FindNextFile(hFind, &FindData) > 0 );
  FindClose(path);
#else
  while ((ent = readdir(dir)) !=NULL);
  closedir(dir);
#endif

#ifdef WIN32
  if( GetLastError() != ERROR_NO_MORE_FILES )
  {
    cout << "Something went wrong during searching\n";
    return false;
  }
  else
  {
#endif
    pairNum = filenames.size()*(filenames.size()-1)/2;
    label = Mat(pairNum,1,CV_8S);
    return true;
#ifdef WIN32
  }
#endif
}

vector<string> Metric::readFiles(string path, string searchPrefix, string searchPattern, string searchExt)
{
  char last = *path.rbegin();

  if (last!='/')
    path = path+"/";
  last = *searchExt.begin();
  if (last!='.')
    searchExt="."+searchExt;

  this->clustercount=0;

#ifdef WIN32
  string fullSearchPath = path + searchPrefix+"*"+searchExt;
  WIN32_FIND_DATA FindData;
  HANDLE hFind;
  wstring wfullSearchPath;
  wfullSearchPath.assign(fullSearchPath.begin(),fullSearchPath.end());
  hFind = FindFirstFile( wfullSearchPath.c_str(), &FindData );
  if( hFind == INVALID_HANDLE_VALUE )
  {
    cout << "Error searching directory\n";
    return filenames;
  }
#else
  DIR *dir;
  struct dirent *ent;
  dir = opendir(path.c_str());
#endif


  do
  {
    string filename;
#ifdef WIN32
    wstring wfilename= FindData.cFileName;
    filename.assign(wfilename.begin(),wfilename.end());
    filenames.push_back(filename);
#else
    if (!strcmp(ent->d_name,"." )) continue;
    if (!strcmp(ent->d_name,"..")) continue;
    string file_name = ent->d_name;
    unsigned long found = file_name.rfind("."+searchExt);
    if (found!=std::string::npos)
    {
        filename = path + "/" + file_name;
        filenames.push_back(filename);
        cout<<"Got: "<<file_name<<endl;
    }
#endif

    std::smatch res;
    //std::regex rx(searchPattern+"_([^]*)_(\\d+)"+searchExt);
    std::regex rx(searchPrefix+searchPattern+searchExt);
    std::regex_search(filename,res,rx);
    std::cout << "class: "<<res[1] << " , num: " << res[2] << ", d? "<<res[3]<<"\n";
    string cluster = string(res[1])+string(res[3]);
    if (clusternames.find(cluster)==clusternames.end())
    {
      clusternames.insert(pair<string,int>(cluster,clustercount));
      clustercount++;
    }
  }
#ifdef WIN32
  while( FindNextFile(hFind, &FindData) > 0 );
  FindClose(path);
#else
  while ((ent=readdir(dir))!=NULL);
  closedir(dir);
#endif

#ifdef WIN32
  if( GetLastError() != ERROR_NO_MORE_FILES )
  {
    cout << "Something went wrong during searching\n";
  }
  else
  {
#endif
    pairNum = filenames.size()*(filenames.size()-1)/2;
    label = Mat(pairNum,1,CV_8S);
    return filenames;
#ifdef WIN32
  }
#endif
}

bool Metric::readFileData(string path, string searchPrefix, string searchPattern, string searchExt, int flag)
{
  char last = *path.rbegin();

  if (last!='/')
    path = path+"/";
  last = *searchExt.begin();
  if (last!='.')
    searchExt="."+searchExt;

  this->clustercount=0;
#ifdef WIN32
  WIN32_FIND_DATA FindData;
  HANDLE hFind;
  string fullSearchPath = path + searchPrefix+"*"+searchExt;
  wstring wfullSearchPath;
  wfullSearchPath.assign(fullSearchPath.begin(),fullSearchPath.end());
  hFind = FindFirstFile( wfullSearchPath.c_str(), &FindData );

  if( hFind == INVALID_HANDLE_VALUE )
  {
    cout << "Error searching directory\n";
    return false;
  }
#else
  DIR *dir;
  struct dirent *ent;
  dir = opendir(path.c_str());
#endif

  FNode* clusterHead, *clusterTail;
  clusterHead = ftable;
  clusterTail = ftable;
  int inClusterCount=0;
  int imagecount=0;
  do
  {
    string filename;
#ifdef WIN32
    wstring wfilename= FindData.cFileName;
    filename.assign(wfilename.begin(),wfilename.end());
    filenames.push_back(filename);
#else
    if (!strcmp(ent->d_name,"." )) continue;
    if (!strcmp(ent->d_name,"..")) continue;
    string file_name = ent->d_name;
    unsigned long found = file_name.rfind("."+searchExt);
    if (found!=std::string::npos)
    {
        filename = path + "/" + file_name;
        filenames.push_back(filename);
        cout<<"Got: "<<file_name<<endl;
    }
#endif
    std::smatch res;
    //std::regex rx(searchPattern+"_([^]*)_(\\d+)"+searchExt);
    std::regex rx(searchPrefix+searchPattern+searchExt);
    std::regex_search(filename,res,rx);
    std::cout << "class: "<<res[1] << " , num: " << res[2] << ", d? "<<res[3]<<"\n";
    string cluster = string(res[1])+string(res[3]);
    if (clusternames.find(cluster)==clusternames.end())
    {
      clusternames.insert(pair<string,int>(cluster,clustercount));
      if (clusterHead==NULL)
      {
        clusterHead = new FNode(cluster,0);
        ftable = clusterHead;
      }
      else
      {
        clusterTail->nextNode = NULL;
        clusterHead->nextClass = new FNode(cluster,0);
        clusterHead->clusterNumber = inClusterCount;
        clusterHead = clusterHead->nextClass;
      }
      clusterTail = clusterHead;
      clustercount++;
      inClusterCount=0;
    }
    clusterTail->nextNode = new FNode(path+filename,stoi(string(res[2])),flag);
    clusterTail->nextNode->clusterNumber = inClusterCount;
    clusterTail->nextClass = clusterHead->nextClass;
    clusterTail = clusterTail->nextNode;
    inClusterCount++;
    imagecount++;
  }
#ifdef WIN32
  while( FindNextFile(hFind, &FindData) > 0 );
  FindClose(path);
#else
  while ((ent = readdir(dir)) !=NULL);
  closedir(dir);
#endif
  clusterHead->nextClass=NULL;
  clusterTail->nextNode=NULL;
  cout<<"No. of clusters:"<<clustercount<<endl;
  cout<<"No. of images:"<<imagecount<<endl;
#ifdef WIN32
  if( GetLastError() != ERROR_NO_MORE_FILES )
  {
    cout << "Something went wrong during searching\n";
    return false;
  }
  else
  {
#endif
    pairNum = filenames.size()*(filenames.size()-1)/2;
    label = Mat(pairNum,1,CV_8S);
    return true;
#ifdef WIN32
  }
#endif
}


int Metric::computeStats(string path, string searchPattern, string searchExt)
{
  stats = vector<vector<Tensor<double,2>>>(filenames.size());
  int count = 0;
  coeffNum = 0;
  for (string filename : this->filenames)
  {
    Tensor<double,1> temp(path+filename);
    //temp.Print();
    vector<Tensor<double,2> > stat = ComputeStatistics(temp,subwinSize,subwinStep);
   // stat[0+1].Print();
    //stat[14+1].Print();
    //stat[28+1].Print();
    //stat[42+1].Print();
    stats[count] = stat;
    count++;
    //if (count > 2) break;
  }
  //by product compute number of coeffients (in case downsampling)
  for (auto t : stats[0])
  {
    coeffNum += t.size().area();
    coeff_in_band.push_back( t.size().area() );
  }
  llr = Mat(pairNum,coeffNum,CV_64F);
  dec = Mat(pairNum,coeffNum,CV_64F);
  return stats.size();
}
 int Metric::computeStats(const Tensor<double,1>& im1, const Tensor<double,1>& im2, FeaturePoolType pooltype)
 {
   pairNum = 1;
   coeffNum = 0;
   stats = vector<vector<Tensor<double,2> > >(2);
   stats[0] = ComputeStatistics(im1, subwinSize,subwinStep,this->subsample,changeWin,3,4,FilterBoundary::FILTER_BOUND_FULL,pooltype);
   stats[1] = ComputeStatistics(im2, subwinSize,subwinStep,this->subsample,changeWin,3,4,FilterBoundary::FILTER_BOUND_FULL,pooltype);
   coeffNum=0;
   coeff_in_band.clear();
  for (auto t : stats[0])
  {
    coeffNum += t.size().area();
    coeff_in_band.push_back( t.size().area() );
  }
  llr = Mat(pairNum,coeffNum,CV_64F);
  dec = Mat(pairNum,coeffNum,CV_64F);
  return stats.size();
 }

int Metric::computeStats(const Tensor<double,1>& im, FeaturePoolType pooltype)
{
  pairNum=1;
  coeffNum=0;
  stats = vector<vector<Tensor<double,2>>>(1);
  stats[0] = ComputeStatistics(im, this->subwinSize,this->subwinStep,this->subsample,changeWin,3,1,FilterBoundary::FILTER_BOUND_EXTEND,pooltype,false);
  coeffNum=0;
  coeff_in_band.clear();
  for (auto t : stats[0])
  {
     // t.Print();
    coeffNum += t.size().area();
    coeff_in_band.push_back( t.size().area() );
  }
  llr = Mat(pairNum,coeffNum,CV_64F);
  dec = Mat(pairNum,coeffNum,CV_64F);
  return stats.size();
}
int Metric::computeFeatures(string searchPattern, string searchExt)
{
  std::smatch res;
  string clusterA,clusterB;
  std::regex rx(searchPattern+"_([^]*)_(\\d+)"+searchExt);
  int paircount=0;
  sameCount = 0;
  diffCount=0;
  int posCount =0;
  Mat temp,tempf;
  Tensor<double,2> tempt;
  //vector<Tensor<double,2>> f((Nor*Nsc+2)*featureNum);
  f0 = Mat(pairNum*coeff_in_band[0],(Nsc*Nor+2)*featureNum*2,CV_64F);
  f1 = Mat(pairNum*coeff_in_band[0],(Nsc*Nor+2)*featureNum*2,CV_64F);
  Mat* pf;
  for (unsigned int i=0; i<filenames.size();i++)
  {
    std::cout<<"i= "<<i<<endl;
    for (unsigned int j=i+1; j< filenames.size(); j++)
    {
      string s = filenames[i];
      string p = filenames[j];
      std::regex_search(s,res,rx);
      clusterA = res[1];
      std::regex_search(p,res,rx);
      clusterB = res[1];
      //vector<Tensor<double,2>> f;
      CV_Assert(stats[i].size()==stats[j].size());
      int count=0;
      //bool inverse = 0;
      Rect roi;
      if (~clusterA.compare(clusterB)) //same cluster
      {
        posCount = sameCount;
        pf = &f1;
        sameCount ++;
        //sameClusterFeature.push_back(f);
        //coeffcount1 += ount;
        label.at<INT8>(paircount,0)=1;
      }
      else
      {
        pf = &f0;
        posCount = diffCount;
        diffCount ++;
        //diffClusterFeature.push_back(f);
        //coeffcount0 += count;
        label.at<INT8>(paircount,0)=0;
      }
      for (unsigned int k=0; k<stats[i].size(); k++)
      {

        roi = Rect(k*2,posCount*coeff_in_band[k],2,coeff_in_band[k]);

        count+=coeff_in_band[k];

        (stats[i][k] - stats[j][k])
          .GetFrameRef(0)
          .reshape(1,stats[i][k].size().area())
          .copyTo(pf->operator()(roi));
        //mylib::DisplayMat(pf->operator()(roi));
        count+=coeff_in_band[k];
      }
      //mylib::DisplayMat(pf->operator()(Rect(0,posCount*9,2,9)),"frow",true);

      paircount++;
    }
  }
  CV_Assert(paircount==pairNum);


  //f0 = Mat::zeros(coeffNum*diffClusterFeature.size(),2*featureNum*(Nsc*Nor+2),CV_64F);//duff
  //f1 = Mat::zeros(coeffNum*sameClusterFeature.size(),2*featureNum*(Nsc*Nor+2),CV_64F);
  //int count = 0;
  //int old_count = count;
  //
  //for (int m = 0; m< sameClusterFeature.size();m++)
  //{
  //  for (int k=0; k< sameClusterFeature[m].size(); k++)
  //  {
  //    count = old_count;
  //    for (int i=0; i< sameClusterFeature[m][k].size().height; i++)
  //      for (int j=0; j<sameClusterFeature[m][k].size().width;j++)
  //      {
  //        f1.at<double>(count,k*2) = sameClusterFeature[m][k](i,j,0)[0];
  //        f1.at<double>(count,k*2+1) = sameClusterFeature[m][k](i,j,0)[1];
  //        count++;
  //      }
  //  }
  //  old_count = count;
  //}
  //for (auto c :diffClusterFeature)
  //{
  //  int count=0;
  //  for (int k=0; k< c.size(); k++)
  //  {
  //    int fidx = k%(Nsc*Nor+2);
  //    count = k/(Nsc*Nor+2)*diffClusterFeature.size()*c[k].size().area();
  //    for (int i=0; i< c[k].size().height; i++)
  //      for (int j=0; j<c[k].size().width;j++)
  //      {
  //        f0.at<double>(fidx*diffClusterFeature.size()*c[k].size().area()+count,0) = c[k](i,j,0)[0];
  //        f0.at<double>(fidx*diffClusterFeature.size()*c[k].size().area()+count,1) = c[k](i,j,0)[0];
  //        count++;
  //      }
  //  }
  //}
  coeffcount0 = diffCount*coeffNum/featureNum/(Nsc*Nor+2);
  coeffcount1 = sameCount*coeffNum/featureNum/(Nsc*Nor+2);

  f0 = f0(Rect(0,0,f0.size().width,coeffcount0));
  f1 = f1(Rect(0,0,f1.size().width,coeffcount1));
  mylib::DisplayMat(getFeature(f1,0,1),"f1var",true);

 // mylib::DisplayMat(f1(Rect(0,0,2,coeffcount1)),"f1",true);
  return sameCount+diffCount;
}

int Metric::computeFeatures(const vector<vector<Tensor<double,2>>>& stats)
{
 // cout<<pairNum<<","<<featureNum<<","<<coeff_in_band[0]<<endl;
  //cout<<"number of features are"<<  stats[0].size()<<endl;
  f = Mat(pairNum*coeff_in_band[0],(Nsc*Nor+2)*featureNum*2+2*(mylib::combination(Nor,2)*Nsc+(Nsc-1)*Nor),CV_64F);
  //cout<<(Nsc*Nor+2)*featureNum*2+2*(mylib::combination(Nor,2)*Nsc+(Nsc-1)*Nor)<<endl;
  Rect roi;
   for (unsigned int k=0; k<stats[0].size(); k++)
   {
      roi = Rect(k*2,0,2,coeff_in_band[k]);
      //stats[0][k].Print();
      //stats[1][k].Print();
      //abs difference
      //(stats[0][k] - stats[1][k]).Print();
      //0 is org, 1 is cand
      //do Weber's law (cand - org)/org
      Mat logit;
      //Mat temp=   (stats[1][k] - stats[0][k]).GetFrameRef(0).clone();
      //Tensor<double,2> num(temp);
      if (stats.size()==2)
      {
        Tensor<double,2> num(stats[1][k]-stats[0][k]);
        Tensor<double,2> dem(stats[0][k]);
        Tensor<double,2> web = num.Abs()/dem.Abs();
        logit = web.GetFrameRef(0).clone();
      }
      else if (stats.size()==1)
      {
        logit = stats[0][k].GetFrameRef(0).clone();
      }
      else
      {
          CV_Assert(stats.size()==1||stats.size()==2);
      }
      /*try simple div only 20130805
      //num.Print();
      //dem.Print();
      //web.Print();
      Mat abstemp = web.Abs().GetFrameRef(0).clone();
      //normailze

     // mylib::DisplayMat(abstemp);
     // Mat exptemp;
      //mylib::DisplayMat(temp);
      //mylib::DisplayMat(abstemp);
      vector<Mat> tt;
      cv::split(abstemp,tt);
      cv::exp(-tt[0],exptemp);
      //mylib::DisplayMat(-tt[0]);
      //mylib::DisplayMat(exptemp);
      exptemp=exptemp+1;
      cv::divide(2.0,exptemp,tt[0]);
      tt[0] = 2 - tt[0];
      cv::merge(tt,logit);
      //mylib::DisplayMat(logit);
      //mylib::DisplayMat(temp);
      //mylib::DisplayMat(exptemp);
      //temp = exptemp/(exptemp+1);
      //Mat logis = 1 - abs(temp);
      //mylib::DisplayMat(temp);
     /*/ // try simple div only
      //mylib::DisplayMat(temp);
      //normalize each feature
      //if (k/(Nsc*Nor+2)==0)//means max 255
      //  temp = temp/255;
      //else if (k/(Nsc+Nor+2)==1)//variance max is given by max((mu-min_pix)*(max-pix - mu)) about 128^2
      //  temp = temp/16384;
      //for rho term, which are crosscoefficient , is range from 0 ~ 1 no need to normalize
     // mylib::DisplayMat(temp);

      //abs
      //tt[0]=logis;
      //tt[1]=cv::Mat::zeros(temp.size(),temp.type());//cv::abs(tt[1]);
      //satruation
      //cv::threshold(tt[0],tt[0],1,1,cv::THRESH_TRUNC);
      //cv::threshold(tt[1],tt[1],1,1,cv::THRESH_TRUNC);
      //cv::merge(tt,temp);
      //mylib::DisplayMat(temp);
      //inverse
      //if (k/(Nsc*Nor+2)==0||k/(Nsc*Nor+2)==1)
      //  temp = Scalar(1,1) - temp;
      //mylib::DisplayMat(temp);
      //save to feature vector
     // cout<<roi<<endl;
     // cout<<f.size()<<endl;
     // cout<<logit.size()<<endl;
      logit.reshape(1,stats[0][k].size().area())
          .copyTo(f(roi));
     //mylib::DisplayMat(f,"f",true);
   }

   return 1;
}

Mat Metric::getFeature(const Mat& f, int subband, int feature)
{
  int fnum = ((Nsc*Nor+2)*feature + subband)*2;
  //int fnum = stats[0][(Nor*Nsc+2)*feature+subband].size().area();//get the featurenum in band

  Rect roi(fnum,0,2,f.size().height);
  return f(roi);
}


Mat& Metric::estimateML(const Mat& f, int subband, int feature, vector<vector<Mat_<double>>>& gamma, vector<vector<double>>& lambda)
{
  //f is the input feature (all the features)
  //subband is the subband index (0-13)
  //feature is the feature index: 0 - mu, 1 - sigma , 2 rho01 , 3 rho10

  Mat tempmean,covar;
  Mat tempf = getFeature(f,subband,feature);
  //mylib::DisplayMat(tempf,"c.txt",true);
  cv::calcCovarMatrix(getFeature(f,subband,feature),covar,tempmean,CV_COVAR_NORMAL|CV_COVAR_ROWS);
  covar = covar/(tempf.size().height-1);
  //mylib::DisplayMat(covar);
  double dr = cv::determinant(covar);
  dr = sqrt(dr);
 // if (subband==13&&feature==3)
  //{
  //  cout<<dr<<endl;
 // }
  if (dr<1e-6)
    lambda[subband][feature] = sqrt(abs(covar.at<double>(0,0)));
  else
    lambda[subband][feature] = dr;
  gamma[subband][feature] =  (covar/lambda[subband][feature]).inv(DECOMP_SVD);
  //if (subband==13&&feature==3)
  //{
  //  mylib::DisplayMat(tempf,"c.txt",true);
   // mylib::DisplayMat(covar/lambda[subband][feature]);
  //}
  //out put gamma is pinv of gamma
  //std::cout<<lambda[subband][feature]<<endl;
  //mylib::DisplayMat(gamma[subband][feature]);
  return gamma[subband][feature];
}
void Metric::trainMetirc(string path,string scorefilepath)
{
  ifstream scorefile;
  if (scorefilepath.empty())
    scorefile.open(path+"/subtestoutput.txt");
  else
    scorefile.open(scorefilepath);

  //std::regex rx(searchPattern+"_([^]*)_(\\d+)"+searchExt);
  char temp[1024];
  //string temp;
  /*int count=0;
  while(scorefile.getline(temp,1024))
    count++;
  scorefile.close();*/
  //scorefile.open(path+"/subtestoutput.txt");
  vector<Mat> features;
  vector<double> scores;
  while(scorefile.getline(temp,1024))
  {
    string str(temp);
    std::smatch res;
    std::regex rx("\\b([^ ]*),");
    std::regex_search(str,res,rx);
    cout<<"org: "<<res[1]<<", ";
    Tensor<double,1> org(path+"/"+string(res[1]));
    org = org.Crop((org.size()/4).Point3(),org.size()/2+Point3i(0,0,1));
    str = res.suffix().str();
    std::regex_search(str,res,rx);
    cout<<"cand "<<res[1]<<endl;
    Tensor<double,1> cand(path+"/"+string(res[1]));
    cand = cand.Crop((cand.size()/4).Point3(),cand.size()/2+Point3i(0,0,1));
    str = res.suffix().str();
    std::regex_search(str,res,rx);
    scores.push_back((double)std::stoi(res[1]));
    computeStats(org,cand);
    computeFeatures(stats);
    //do LS
    features.push_back(f);
    //mylib::DisplayMat(f);
  }
  cout<<features.size()<<endl;
  Mat A(features.size(),coeffNum*2,CV_64F);
  Mat b(features.size(),1,CV_64F);
  Mat bl(features.size(),1,CV_64F);
  cout<<features[0].at<double>(0,0)<<","<<features[1].at<double>(0,0)<<endl;
  for (unsigned int i=0; i< features.size(); i++)
  {
    Mat temp(coeffNum,2,CV_64F);
    int copypos =0;
    //mylib::DisplayMat(features[i],"featurei",true);
    for (unsigned int j=0; j<coeff_in_band.size();j++)
    {
      Rect roi(j*2,0,2,coeff_in_band[j]);
      Rect toroi(0,copypos,2,coeff_in_band[j]);
      features[i](roi).copyTo(temp(toroi));
    //  mylib::DisplayMat(features[i](roi));
    //  mylib::DisplayMat(temp,"temp",true);
      copypos+=coeff_in_band[j];
    }
    //mylib::DisplayMat(temp);
    temp.reshape(1,1).copyTo(A.row(i));
    b.at<double>(i,0) = scores[i]/10;
    //if (b.at<double>(i,0)>0.99)
    //  b.at<double>(i,0)=0.99;
    //if (b.at<double>(i,0)<0.01)
    //  b.at<double>(i,0)=0.01;
    //bl.at<double>(i,0) = log(b.at<double>(i,0)/(1 - b.at<double>(i,0)));
  }
  mylib::DisplayMat(A,"A",true);
  mylib::DisplayMat(b,"b",true);
  Mat x;
  //cv::solve(A,b,x,DECOMP_SVD);
  //x = (A.t()*A).inv()*A*b;
  cv::invert(A,x,cv::DECOMP_SVD);
  x = x*b;
  weights = x;
  mylib::DisplayMat(x,"lse",true);
  scorefile.close();
  //compute err
  //Mat tempmat = A*x;
  //Mat expmat;
  //exp(tempmat,expmat);
  Mat rst  = A*x;//expmat/(expmat+1);

  mylib::DisplayMat(rst-b,"err",true);
  double err = cv::sum(cv::abs(rst-b))[0];
  cout<<"abs error is: "<<err<<endl;
  ofstream coefffile;
  coefffile.open("coeff_in_band.txt");
  for (unsigned int i=0; i< coeff_in_band.size();i++)
  {
    coefffile<<coeff_in_band[i]<<";"<<endl;
  }
  coefffile.close();
  ofstream weightfile;
  weightfile.open("weight.txt");
  for (int i=0; i<x.size().height; i++)
  {
    weightfile<<x.at<double>(i,0)<<",";
  }
  weightfile.close();
}
void Metric::trainMetric(string path, string searchPattern, string searchExt)
{
  if(readFiles(path,searchPattern, searchExt)!=true)
  {
    cout<<"read database error\n";
    exit(-1);
  }
  computeStats(path,searchPattern,searchExt);
  computeFeatures(searchPattern,searchExt);
  //int count = 0;
  //Rect roi;
  dec = Mat(coeffNum,pairNum,CV_64F);
  for (int subband = 0; subband<Nor*Nsc+2; subband++)
  {
    cout<<subband<<endl;
    for (int feature = 0; feature < featureNum; feature++)
    {
      estimateML(f0,subband,feature,gamma0,lambda0);
      estimateML(f1,subband,feature,gamma1,lambda1);
      Mat llr1 = llr(Rect(0,0,coeffNum,sameCount));
      Mat llr0 = llr(Rect(0,sameCount,coeffNum,diffCount));
      computeLLR(f1,llr1,subband,feature, gamma0,gamma1,lambda0,lambda1);
      computeLLR(f0,llr0,subband,feature,gamma0,gamma1,lambda0,lambda1);
    }
  }
  saveGamma("gamma0",gamma0);
  loadGamma("gamma0",gamma0);
  saveGamma("gamma1",gamma1);
  loadGamma("gamma1",gamma1);
  saveLambda("lambda0",lambda0);
  loadLambda("lambda0",lambda0);
  saveLambda("lambda1",lambda1);
  loadLambda("lambda1",lambda1);
  //mylib::DisplayMat(llr,"llr",true);
  makeDec(llr,dec,log(1.0),DEC_FLAG_SOFT);

  label = -1*cv::Mat::ones(pairNum,1,CV_32F);
  label(Rect(0,0,1,sameCount))=1;
  params.svm_type = SVM::NU_SVR;
  params.kernel_type = SVM::LINEAR;
  params.term_crit = TermCriteria(CV_TERMCRIT_ITER,100,1e-6);
  params.nu = 0.5;
  Mat floatdec;
  dec.convertTo(floatdec,CV_32F);
  //mylib::DisplayMat(floatdec,"dec.txt",true);
  s.train_auto(floatdec,label,Mat(),Mat(),params);
  Mat lb(label.size(),label.type());
  s.save("svm");
  SVM s2;
  s2.load("svm");
  for (int i=0; i < pairNum;i++)
  {

    lb.at<float>(i,0) = s2.predict(floatdec.row(i));
  }

  cv::compare(lb,0,lb,CV_CMP_GE);
  lb.convertTo(lb,CV_32F);
  lb = lb/255*2.0 -1.0;
  mylib::DisplayMat(lb,"lb",true);
  mylib::DisplayMat(label,"label",true);
  double rst = cv::sum((lb-label))[0];
  cout<<rst<<endl;

}

Mat& Metric::computeLLR(const Mat& f, Mat& llr, int subband, int feature, vector<vector<Mat_<double>>>& gamma0, vector<vector<Mat_<double>>>& gamma1, vector<vector<double>>& lambda0, vector<vector<double>>& lambda1)
{
  CV_Assert(f.size().height/coeff_in_band[0]==llr.size().height);
   Mat tempf = getFeature(f,subband,feature);

   Mat temp0 = tempf*gamma0[subband][feature];
   Mat temp1 = tempf*gamma1[subband][feature];
   Mat temp;
   mylib::DisplayMat(tempf,"f",true);
   mylib::DisplayMat(gamma0[subband][feature]);
   mylib::DisplayMat(gamma1[subband][feature]);
   mylib::DisplayMat(temp0,"temp0",true);
   mylib::DisplayMat(temp1,"temp1",true);
   double p0,p1,normz0,normz1;
  // llr = Mat(coeff_in_band[subband],1,CV_64F);
   for (int i=0; i<temp0.size().height; i++)
   {
     int findex = i%coeff_in_band[subband]+subband*coeff_in_band[subband]+feature*coeff_in_band[subband]*(Nsc*Nor+2);
     int sindex = i/coeff_in_band[subband];
      cv::multiply(temp0.row(i),tempf.row(i),temp);
      normz0 = sqrt(sum(temp)[0]);
      cv::multiply(temp1.row(i),tempf.row(i),temp);
      normz1 = sqrt(sum(temp)[0]);
      p0 = 2*cyl_bessel_k(0,sqrt(2*normz0/lambda0[subband][feature]))/(2*pi)/lambda0[subband][feature];
      p1 = 2*cyl_bessel_k(0,sqrt(2*normz1/lambda1[subband][feature]))/(2*pi)/lambda1[subband][feature];
      llr.at<double>(sindex,findex) = log(p1/p0);
     // cout<<llr.at<double>(sindex,findex)<<endl;
   }

   return llr;
}

Mat& Metric::makeDec(const Mat& LLR, Mat& dec, double thred, int flag )
{
  if (flag == DEC_FLAG_HARD)
  {
    cv::compare(LLR,thred,dec,CMP_GT);
  }
  else
  {
    dec = LLR - thred;
  }
  return dec;
}

bool Metric::loadClassifier(string filename)
{
  try
  {
  this->s.load(filename.c_str());
  }
  catch (Exception e)
  {
    cout<<"load classifier error\n";
    return false;
  }
  return true;
}

float Metric::predict(const Mat& f)
{
  return s.predict(f);
}

void Metric::saveGamma(string filename, const vector<vector<Mat_<double>>>& gamma)
{
  cv::FileStorage fs(filename,cv::FileStorage::WRITE);
  int count=0;
  for (auto gline : gamma)
  {
    for (auto g: gline)
    {
      fs << filename+std::to_string(count)<<g;
      count++;
    }
  }
  fs.release();
}

void Metric::loadGamma(string filename, vector<vector<Mat_<double>>>& gamma)
{
  cv::FileStorage fs(filename,cv::FileStorage::READ);
  int count = 0;
  for (int i=0; i< (Nor*Nsc+2); i++)
    for (int j=0; j< featureNum; j++)
    {
      fs[filename+std::to_string(count)]>>gamma[i][j];
      count++;
    }

  fs.release();
}

void Metric::saveLambda(string filename, const vector<vector<double>>& lambda)
{
  //convert lambda to mat
  cv::Mat L(lambda.size(),lambda[0].size(),CV_64F);
  for (unsigned int i=0; i< lambda.size();i++)
    for (unsigned int j=0; j< lambda[0].size(); j++)
      L.at<double>(i,j) = lambda[i][j];
  cv::FileStorage fs(filename,cv::FileStorage::WRITE);
  fs << filename<<L;
  fs.release();
}

void Metric::loadLambda(string filename, vector<vector<double>>& lambda)
{
  cv::FileStorage fs(filename,cv::FileStorage::READ);
  Mat_<double> L ;
  fs[filename]>>L;
  for (unsigned int i=0; i< lambda.size();i++)
    for (unsigned int j=0; j< lambda[0].size(); j++)
      lambda[i][j] = L.at<double>(i,j);
  fs.release();
}

double Metric::computeMetric(const Tensor<double,1>& im1, const Tensor<double,1>& im2)
{
  computeStats(im1,im2);
  computeFeatures(this->stats);
  for (int i=0; i< Nsc*Nor+2; i++)
    for (int j=0; j< featureNum; j++)
      computeLLR(f,this->llr,i,j,gamma0,gamma1,lambda0,lambda1);
  makeDec(llr,dec,log(1.0),DEC_FLAG_SOFT);
  Mat floatdec;
  dec.convertTo(floatdec,CV_32F);
  double rst = s.predict(floatdec);
    return rst;
}

double Metric::computeMetricLSE(const Tensor<double,1>& im1, const Tensor<double,1>& im2)
{
  computeStats(im1,im2);
  computeFeatures(this->stats);
    Mat A(1,coeffNum*2,CV_64F);
    Mat temp(coeffNum,2,CV_64F);
    int copypos =0;
   // mylib::DisplayMat(f,"f",true);
    for (unsigned int j=0; j<coeff_in_band.size();j++)
    {
      Rect roi(j*2,0,2,coeff_in_band[j]);
      Rect toroi(0,copypos,2,coeff_in_band[j]);
      f(roi).copyTo(temp(toroi));
      copypos+=coeff_in_band[j];
    }
    temp.reshape(1,1).copyTo(A.row(0));
    //mylib::DisplayMat(A,"A_test",true);
    Mat scorerst = A*weights;
    mylib::DisplayMat(scorerst);
    return scorerst.at<double>(0,0)/10;
}

void Metric::loadParams(void)
{
  loadGamma("gamma0",gamma0);
  loadGamma("gamma1",gamma1);
  loadLambda("lambda0",lambda0);
  loadLambda("lambda1",lambda1);
  loadClassifier("svm");
}

void Metric::trainGranularity(string path, string scorefilepath,FeaturePoolType pooltype)
{
  cout<<this->subwinSize<<endl;
  char curpath[FILENAME_MAX];
  GetCurDir(curpath,sizeof(curpath));
  cout<<"reading scorefile: "<<string(curpath)+scorefilepath<<endl;
  ifstream scorefile;
  scorefile.open(string(curpath)+scorefilepath,std::ios::in);
  char temp[1024];
  vector<Mat> features;
  vector<double> scores;
  vector<int> sizes;
  while(scorefile.getline(temp,1024))
  {
    string str(temp);
    cout<<str<<endl;
    boost::smatch res;
    boost::regex rx("\\b(.*)_org\\.png");//find everyword end with _org.png
    boost::regex_search(str,res,rx);
    cout<<"org: "<<res[0]<<", "<<res[1]<<", "<<endl;
    cout<<"path"<<path<<endl;
    Tensor<double,1> org(string(curpath)+path+"/"+string(res[0]));
    rx.set_expression("\\b(\\d),");
    boost::regex_search(str,res,rx);
    int yesSize = Granulate::sz_idx[std::stoi(res[1])];
    scores.push_back(1.0);
    this->subwinSize=Size3(16,16,1);
    this->subwinStep=Size3(16,16,1);
    computeStats(org,FeaturePoolType::FEATURE_POOL_ALL);
   // cout<<"stat size"<<stats[0].size()<<endl;
    computeFeatures(stats);
    features.push_back(f);
    sizes.push_back(yesSize);
  }
  //after reading done
  //1. ouptut statistcs for viewing
  Mat A(features.size(),coeffNum*2,CV_64F);
  Mat b(features.size(),1,CV_64F);
  Mat bl(features.size(),1,CV_64F);
  Mat ms(features.size(),1,CV_32F);
  //cout<<features[0].at<double>(0,0)<<","<<features[1].at<double>(0,0)<<endl;
  cout<<A.size()<<endl;
  for (unsigned int i=0; i< features.size(); i++)
  {
    Mat temp(coeffNum,2,CV_64F);
    int copypos =0;
    //mylib::DisplayMat(features[i],"featurei",true);
    for (unsigned int j=0; j<coeff_in_band.size();j++)
    {
      Rect roi(j*2,0,2,coeff_in_band[j]);
      Rect toroi(0,copypos,2,coeff_in_band[j]);
      features[i](roi).copyTo(temp(toroi));
    //  mylib::DisplayMat(features[i](roi));
    //  mylib::DisplayMat(temp,"temp",true);
      copypos+=coeff_in_band[j];
    }
    //mylib::DisplayMat(temp);
    temp.reshape(1,1).copyTo(A.row(i));
    b.at<double>(i,0) = scores[i];
    if (scores[i]==1)
      bl.at<double>(i,0)= 1;
    else
      bl.at<double>(i,0)=-1;
    //if (b.at<double>(i,0)>0.99)
    //  b.at<double>(i,0)=0.99;
    //if (b.at<double>(i,0)<0.01)
    //  b.at<double>(i,0)=0.01;
    //bl.at<double>(i,0) = log(b.at<double>(i,0)/(1 - b.at<double>(i,0)));

   ms.at<float>(i,0)=float(sizes[i]);
  }
  mylib::DisplayMat(A,"A",true);
  mylib::DisplayMat(b,"b",true);
  mylib::DisplayMat(bl,"bl",true);
  mylib::DisplayMat(ms,"ms",true);
  ofstream coefffile;
  coefffile.open("./temp/coeff_in_band.txt");
  for (unsigned int i=0; i< coeff_in_band.size();i++)
  {
    coefffile<<coeff_in_band[i]<<";"<<endl;
  }
  coefffile.close();
  cout<<"feature exaction done"<<endl;
  //2. make a SVM
  bl.convertTo(label,CV_32F);
  params.svm_type = SVM::NU_SVR;
  params.kernel_type = SVM::RBF;
  params.term_crit = TermCriteria(CV_TERMCRIT_ITER,100,1e-6);
  params.nu = 0.5;
  Mat floatA;
  A.convertTo(floatA,CV_32F);
  //mylib::DisplayMat(floatdec,"dec.txt",true);
  s.train_auto(floatA,label,Mat(),Mat(),params);
  s.save("./temp/svm");
  Mat lb(label.size(),label.type());
  SVM s2;
  s2.load("./temp/svm");
  for (int i=0; i < label.size().height;i++)
  {

    lb.at<float>(i,0) = s2.predict(floatA.row(i));
  }

  cv::compare(lb,0,lb,CV_CMP_GE);
  lb.convertTo(lb,CV_32F);
  lb = lb/255*2.0 -1.0;
  mylib::DisplayMat(lb,"lb",true);
  mylib::DisplayMat(label,"label",true);
  double rst = cv::sum((lb-label))[0];
  cout<<rst<<endl;

  cout<<"svm done"<<endl;

}
}

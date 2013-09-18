#define  _CRT_SECURE_NO_WARNINGS
#define SQL_NOUNICODEMAP
//#include <Windows.h>
#ifdef WIN32
#include <direct.h>
#endif

#include <MTC.h>
#include  "jpeglib.h"
#include "JpegAddon.h"
#include <bitset>
//#include <Windows.h>
//#include <sql.h>
//#include <sqltypes.h>
//#include <sqlext.h>
#include <string>
#include <cmath>
#include <boost/filesystem.hpp>
#include <sys/dir.h>
#include <Metric.h>
namespace mtc  {

  bool debug_switch = false;

  ////comment out for none debug run
  //#ifndef DEBUG
  //#define DEBUG
  //#endif
  ////
  extern const int zig_zag_x[];
  extern const int zig_zag_y[];

  MTC::MTC(void)
  {
  }


  MTC::~MTC(void)
  {
  }
  MTC::MTC(const string& imageName, CodingMode codemode):QGrid<MTC::T,MTC::cn>(imageName)
  {
    this->acount=0;
    if (codemode == CodingMode::CODING_MODE_JPEG || codemode ==CodingMode::CODING_MODE_MTC)
      {
        SetJpegQFactor(30);//default value, dont change
        InitJpegParam(jpegQNum);
      }
    SetInitBlockSize(Size3(32,32,1));
    SetSTSIMSubWinSize(Size3(16,16,1));
    SetSTSIMSubWinStep(this->subSize);//means non_overlap
    SetQFactor(1.3);
    SetSTSIMQualigyThrd(0.86);
    SetInitQSize(8);
    SetCodingMode(codemode);
    SetOverlapSizeByRatio(Size3_<double>(0.25,0.25,0));
    SetSearchStep(Size3(4,4,1));
    SetCandidNum(4);
    max_bits = 70;//internal parameters, does not change frequently
    //rst_with_seam = rst.CvtColor<3>(CV_GRAY2RGB);
    SetPQIRectType(false);
    SetLCType(LightingCorrectionType::HAS_LIGHTING_CORRECTION);
    SetPBType(PostBlendingType::POST_BLENDING_ONLINE);
    //UpdateParameters();
    acceptDirect=false;
    mseThrd = 0;//default accept all
  }
  void MTC::Init(const string& imageName, CodingMode codemode)
  {
    this->acount=0;
    this->cFileName = imageName;
    ensemble.Load(cFileName);
    //ensemble = ensemble.Log();//20130708
    if (codemode == CodingMode::CODING_MODE_JPEG || codemode == CodingMode::CODING_MODE_MTC)
      {
        SetJpegQFactor(30);//default value, dont change
        InitJpegParam(jpegQNum);
      }
    SetInitBlockSize(Size3(32,32,1));
    SetSTSIMSubWinSize(Size3(16,16,1));
    SetSTSIMSubWinStep(this->subSize);//means non_overlap
    SetQFactor(1.3);
    SetSTSIMQualigyThrd(0.86);
    SetInitQSize(8);
    SetCodingMode(codemode);
    SetOverlapSizeByRatio(Size3_<double>(0.25,0.25,0));
    SetSearchStep(Size3(4,4,1));
    SetCandidNum(4);
    max_bits = 70;//internal parameters, does not change frequently
    //rst_with_seam = rst.CvtColor<3>(CV_GRAY2RGB);
    SetPQIRectType(false);
    SetLCType(LightingCorrectionType::HAS_LIGHTING_CORRECTION);
    SetPBType(PostBlendingType::POST_BLENDING_ONLINE);
    //UpdateParameters();
    acceptDirect=false;
    //InitGrid();
  }


  MTC::MTC(const string& imageName, double qualityThrd , const Size3& initSize, Size3 searchStep, int candidNum, double initQSize, const Size3_<double>& overlapRatio):
    QGrid<MTC::T,MTC::cn>(imageName,initSize,overlapRatio)
  {
    this->acount=0;
    this->initSize = initSize;
    this->initQSize = initQSize/ pow(1.3,(double)initSize.height/16.0 - 1);
    this->qualityThrd = qualityThrd;
    this->qfactor = 1.3;
    this->max_bits = 70;
    this->searchStep = searchStep;
    this->candidNum = candidNum;
    this->cFileName = imageName;
    //as in tcy code
    /*double arr[] = {6,9,11,13,15,15,15};
        this->qstep = vector<double>(arr, arr+ sizeof(arr)/sizeof(double));
        this->thrd = qstep;
*/
    //UpdateParameters();
    //rst = Tensor_<T,cn>(ensemble.size());
    //rst_with_seam = rst.CvtColor<3>(CV_GRAY2RGB);
    rectType = false;
    stat = Statistics(levels,initSize.height);
    footmap = Tensor<T,1>(ensemble.size(),Vec<T,1>::all(255));
    dctEnsemble = Tensor<T,cn>(ensemble.size());
    acceptDirect=false;
  }
  void MTC::DebugTest(void)
  {

    //i want to test if change node in the grid changes the rst: yes
    QNode<double,1> qNode = GetNode(Point3i(2,2,0));
    qNode.SetBlock(Tensor<double,1>(qNode.size(),255));
    //now test if boundary also changes rst: yes
    qNode.GetBoundaryLeft().SetBlock(Tensor<double,1>(qNode.leftBound.size(),128));
    qNode.GetBoundaryUp().SetBlock(Tensor<double,1>(qNode.upBound.size(),128));
    qNode.rightBound.SetBlock(Tensor<double,1>(qNode.rightBound.size(),128));
    qNode.lowBound.SetBlock(Tensor<double,1>(qNode.lowBound.size(),128));
    //now change rst, to see whether qNode changes:
    rst.SetBlock(Point3i(58,58,0),Tensor<double,1>(qNode.upBound.size()/2+Size3(0,0,1),200));
    rst.SetBlock(Point3i(72,72,0),Tensor<double,1>(Size3(16,16,1),200));
    rst.Display();
    qNode.upBound.Display();
    qNode.Display();
    qNode.eNode.Display();
    qNode.eNode.SetBlock(Point3i(20,20,0),Tensor<double,1>(16,16,1,0));
    rst.Display();
    //test tree
    QTree<double,1> qTree= GetNode(Point3i(2,2,0));
    qTree.Split();
    qTree.GetSubTree(0)->Display();
    //qTree.GetSubTree(0)->leftBound.Display();
    //qTree.GetSubTree(0)->upBound.Display();
    qTree.GetSubTree(0)->GetExtendTensor().Display();
    qTree.GetSubTree(3)->eNode.Display();
    qTree.GetSubTree(3)->SetBlock(Point3i(8,8,0),Tensor<double,1>(Size3(8,8,1),128));
    rst.Display();
    // qNode.bitcount.SetCodeBit(10);
    //qNode.bitcount.SetCodeBit(10);
    qTree.bitcount.SetDecBit("2");

  }


  void MTC::DummyCoding(void)
  {

    if (unitTest)
      {
        nts.AssertEQ(this->cFileName,"baboon.pgm","Filename");
      }

    if (lightCorrectionType == LightingCorrectionType::PQI_LF_ENCODING)
      {
        this->PQICodingLFComponent();
      }
    else if (lightCorrectionType ==  LightingCorrectionType::PREDEF_LIGHTING)
      {
        this->LoadPreDefLighting();
      }
    if (this->mode != CodingMode::CODING_MODE_TPSS)
      ComputeEnsembleDCT(8);
    if (this->metricModifier == MetricModifier::MAHALANOBIS_DIST) //MAHALANOBIS
      {
        this->iMahaCovar = metric::EstimateVarForMahalanobis(ensemble,this->GetSubWinSize(),this->GetSubWinStep());
        mylib::DisplayMat(this->iMahaCovar);
      }
    this->SAT = rst.ExtendHalfBoundary();
    //SAT.Display();
    //SAT.Display();
    //rst = ensemble;
    //this->UpdateRstSqr(Cube(ensemble.offset(),ensemble.size()));
    //this->rstSqr(Cube(0,0,0,16,16,1)).Print();
    acount=0;
    //this->SAT = rstSqr.ComputeSAT();
    //SAT(Cube(0,0,0,16,16,1)).Print();
    int tempTotalBits=0, tempPQIBits=0, tempTreeBits=0, tempModeBits=0;
    vector<int> tempJPEGBits(2), tempMTCBits(int(log(blockSize.height/8)/log(2))+1);
    cout<<stat.plcBitLength<<endl;
    for (int t=0; t< gridSize.depth; t++)
      for (int i= 0; i < gridSize.height; i++)
        {
          for (int j=0; j< gridSize.width; j++)
            {
              if (mode == CodingMode::CODING_MODE_POST_TP)
                {//unmaintained
                  QTree<T,cn> &tempNode = GetNode(Point3i(i,j,0));
                  PostCoding(tempNode);//<T,cn>(rst,tempNode.size(),tempNode.offset(),tempNode.overlap()));
                }
              else
                {
                  QTree<T,cn> qTree = GetNode(Point3i(i,j,0));
                  //qNode.SetBlock(ensemble.Crop(qNode.offset(),qNode.size()));
                  //(rst).Display();
                  //ModifyMetric(Point3i(i,j,t));
                  //not using oracle
                  cout<<"("<<qTree.offset().x<<","<<qTree.offset().y<<")\n";
                  if (i==2 && j==25)
                    cout<<"debug now"<<endl;
                  CodingCore(qTree);
                  cout<<stat.plcBitLength<<endl;
                  //rst.Display();
                  if (i*qTree.size().height==DEBUG_DISP_X && j*qTree.size().width==DEBUG_DISP_Y)
                    {
                      rst_with_seam.SaveBlock("rst_seam_debug.tiff");
                      rst_with_seam.Display();
                      //rst_with_seam.Print("rst_seam_temp",true);
                      //rst.Display();
                    }
                  ////(rst+128).Display(500,1);
                  //(rst_with_seam+128).Display(500,1);
                  tempTotalBits+=qTree.CollectBits();
                  vector<int> temp = qTree.CollectJPEGBits();
                  tempJPEGBits[0]+=temp[0];
                  tempJPEGBits[1]+=temp[1];
                  //          temp = qTree.CollectMTCBits();
                  //        for (int j=0; j<temp.size();j++)
                  //        tempMTCBits[j]+=temp[j];
                  tempModeBits+=qTree.CollectModeBits();
                  tempTreeBits+=qTree.CollectTreeBits();
                  tempPQIBits+=qTree.CollectPQIBits();
                }


            }
        }
    //save rsts to a video
#ifdef RECORD_EVERYTHING

    //save rsts to a video
    cv::VideoWriter writer;
    cv::VideoCapture cap("./everything/rst_%04d.bmp");
#if CV_MINOR_VERSION < 6
    writer.open("rst.avi",CV_FOURCC('Y','U','V','9'),1,Size(rst.size().width,rst.size().height),true);
#else
    writer.open("rst.avi",cv::VideoWriter::fourcc('Y','U','V','9'),1,Size(rst.size().width,rst.size().height),true);
#endif
    cv::Mat frame;
    cv::Mat frame2;
    for(;;)
      {
        cap>>frame;
        writer.write(frame);
        if(frame.empty())
          break;
      }
    writer.release();
    cap.release();
#endif
    //rst.Display();
    if (this->mode == CodingMode::CODING_MODE_TPSS)
      std::cout<<endl<<"lighting error PSNR"<<10.0*log10f(255.0*255.0/(stat.lightingDiff/double(rst.size().area())))<<endl;
    std::cout<<endl<<"pre compute bits: "<<tempTotalBits;
    std::cout<<endl<<"pre compute bpp: "<<double(tempTotalBits)/double(ensemble.size().area());
    std::cout<<endl<<"PQI bits: "<<double(tempPQIBits);
    std::cout<<endl<<"done~~~~~~~~~~~!"<<endl;
    Tensor<T,cn> tempRS = (rst*rst);
    Vec<T,cn> vRs = (tempRS-rstSqr).Sum();

    //std::cout<<"square: "<<vRs[0]<<endl;
    //(rst*rst)(Cube(12,92,0,20,20,1)).Print();
    //rstSqr(Cube(12,92,0,20,20,1)).Print();
    //Vec<T,cn> vSS = (rst*rst).Sum()-SAT(SAT.size().height-1,SAT.size().width-1,SAT.size().depth-1);
    //Tensor<T,cn> tempSS = (rst*rst).ComputeSAT();
    //Vec<T,cn> vSS = (tempSS-SAT).Sum();
    //std::cout<<"SAT: "<<vSS[0]<<endl;
    //SAT(Cube(12,92,0,21,21,1)).Print();
    //rst = rst.Exp(); //20130708
    UpdateLog();

  }
  void MTC::CollectBits(void)
  {
    Statistics stat2;
    //  int temp;
    vector<int> tempvec;
    for (int t=0; t< gridSize.depth; t++)
      for (int i= 0; i < gridSize.height; i++)
        for (int j=0; j< gridSize.width; j++)
          {
            QTree<T,cn> qTree = GetNode(Point3i(i,j,0));

          }
  }
  void MTC::UpdateLog(void)
  {
    time_t t = time(&stat.endT);
    tm* tt = localtime(&t);
    stat.duration = double(stat.endT - stat.beginT);
    //stat.duration = difftime(stat.endT,stat.beginH);
    logfile<<endl;
    //stat.duration = (stat.endH - stat.beginH)*3600+(stat.endM-stat.beg
    std::cout<<rst.size().height<<","<<rst.size().width<<endl;
#ifdef WIN32
    string PID = boost::lexical_cast<string>(_getpid());
#else
    string PID = boost::lexical_cast<string>(getpid());
#endif
    //comment out on Jan 7
    if (lightCorrectionType==LightingCorrectionType::NO_LIGHTING_CORRECTION|| lightCorrectionType==LightingCorrectionType::POISSON_LC)
      (rst).SaveBlock(path+count_exp+"_"+PID+"_rst.tif");
    else
      (rst+LFImage).SaveBlock(path+count_exp+"_"+PID+"_rst.tif");
    rst_with_seam.SaveBlock(path+count_exp+"_"+PID+"_seam.tif");
    pred_rst.SaveBlock(path+count_exp+"_"+PID+"_pred.tif");
    plc_rst.SaveBlock(path+count_exp+"_"+PID+"_plc.tif");
    pred_afterLC.SaveBlock(path+count_exp+"_"+PID+"_predAfterLC.tif");
    pred_postLC.SaveBlock(path+count_exp+"_"+PID+"_predPostLC.tif");
    ensemble.SaveBlock(path+count_exp+"_"+PID+"_org.tif");
    logfile<<"PID :"<<PID<<endl;
    logfile<<"Start time: "<<stat.beginH<<":"<<stat.beginM<<":"<<stat.beginS<<endl;
    logfile<<"End time: "<< tt->tm_hour<<":"<<tt->tm_min<<":"<<tt->tm_sec<<endl;
    //std::cout<<white<<"Duration: "<<stat.duration<<" second."<<endl;
    std::cout<<white<<"Duration: "<<int(stat.duration)/60<<" min, "<<int(stat.duration)%60<<" second."<<endl;
    std::cout<<"PID :"<<PID<<endl;
    if (stat.duration>=60)
      logfile<<"Duration: "<<int(stat.duration)/60<<"min, "<<int(stat.duration)%60<<" second."<<endl;
    else
      logfile<<"Duration: "<<stat.duration<<" second."<<endl;
    logfile<<"Coding Setup:"<<endl;
    if (mode!=CodingMode::CODING_MODE_JPEG||mode!=CodingMode::CODING_MODE_PQI)
      logfile<<"STSIM2 threshold "<<qualityThrd<<endl;
    if (mode!=CodingMode::CODING_MODE_PQI||mode!=CodingMode::CODING_MODE_TPQI)
      logfile<<"JPEG original factor "<<jpegQNum<<endl;
    if (jpegType == JPEGType::JPEG_ADAPTIVE)
      logfile<<"JPEG degrade factor "<<(double)jpegQNum*0.4<<endl;
    logfile<<"JPEG type"<<(int)jpegType<<endl;
    logfile<<"Post blending type"<<(int)postBlendType<<endl;
    logfile<<"Lighting Correction Type:"<<(int)lightCorrectionType<<endl;
    logfile<<"Side Matching Method:"<<(int)matching_method<<endl;
    logfile<<"STSIM2 thrd multiplier "<<qfactor<<endl;
    logfile<<"Candidate Num "<<candidNum<<endl;
    logfile<<"Boundary Matching Threshold"<<mseThrd<<endl;
    logfile<<"Variance Matching thrd "<<varThrd1<<endl;
    logfile<<"Init. Block Size:"<<"("<<initSize.height<<"x"<<initSize.width<<"x"<<initSize.depth<<")"<<endl;
    logfile<<"Init. Q Size:"<<initQSize<<endl;
    logfile<<"Seaching Step:"<<"("<<searchStep.height<<"x"<<searchStep.width<<"x"<<searchStep.depth<<")"<<endl;
    logfile<<"Overlap Size for 32x32:"<<"("<<overlapSize.height<<"x"<<overlapSize.width<<"x"<<overlapSize.depth<<")"<<endl;
    logfile<<"Overlap Size for 16x16:"<<"("<<overlapSize.height/2<<"x"<<overlapSize.width/2<<"x"<<overlapSize.depth/2<<")"<<endl;
    logfile<<"STSIM2 Subwin Size:"<<"("<<this->GetSubWinSize().height<<"x"<<this->GetSubWinSize().width<<"x"<<this->GetSubWinSize().depth<<")"<<endl;
    logfile<<"STSIM2 Subwin Step:"<<"("<<this->GetSubWinStep().height<<"x"<<this->GetSubWinStep().width<<"x"<<this->GetSubWinStep().depth<<")"<<endl;
    logfile<<"STSIM2 pooling method"<<(int)stsim2PoolType<<endl;
    logfile<<"Metric Modifer type"<<(int)metricModifier<<endl;
    long totalMTCPel=0;//compute the total pixels coded by MTC
    for (unsigned int l=0; l<stat.tpBlockNum.size();l++)
      {
        totalMTCPel+=stat.tpBlockNum[l]*(16<<l)*(16<<l);
        std::cout<<"tp block "<<l<<": "<<stat.tpBlockNum[l]<<endl;
      }
    long totalBasePel = 0; //compute the total pixels coded by base coder
    for (unsigned int l=0; l<stat.jBlockNum.size();l++)
      {
        totalBasePel+=stat.jBlockNum[l]*64;
        std::cout<<"j block"<<l<<": "<<stat.jBlockNum[l]<<endl;
      }
    if (this->mode != CodingMode::CODING_MODE_TPSS)
      {
        CV_Assert(totalBasePel+totalMTCPel == ensemble.size().volumn());
      }
    if (this->mode != CodingMode::CODING_MODE_JPEG)
      {
        if(mode == CodingMode::CODING_MODE_TPSS)
          {
            logfile<<"PLC foot bits bpp"<<double(stat.plcBitLength)/double(rst.size().volumn())<<endl;
            std::cout<<"PLC foot bits bpp"<<double(stat.plcBitLength)/double(rst.size().volumn())<<endl;
          }
        else if (lightCorrectionType==LightingCorrectionType::POISSON_LC)
          {

            logfile<<"PLC foot bits bpp"<<double(stat.plcBitLength)/double(totalMTCPel)<<endl;
            std::cout<<"PLC foot bits bpp"<<double(stat.plcBitLength)/double(totalMTCPel)<<endl;
          }
        else if (lightCorrectionType==LightingCorrectionType::PQI_LF_ENCODING)
          {
            std::cout<<white<<"PQI LF Encoding bpp = "<<(double)stat.footbitLength/(double)totalMTCPel<<endl;
            logfile<<white<<"PQI LF Encoding bpp = "<<(double)stat.footbitLength/(double)totalMTCPel<<endl;

          }
        logfile<<"Foot Compute Region"<<this->footComputeRegion<<endl;
        logfile<<"Foot Compute Method"<<this->footComputeMethod<<endl;
        std::cout<<"Foot Compute Region"<<this->footComputeRegion<<endl;
        std::cout<<"Foot Compute Method"<<this->footComputeMethod<<endl;

        //std::cout<<white<<"PQI LF Encoding bpp = "<<(double)stat.footbitLength/(double)ensemble.size().volumn()<<endl;
        //logfile<<white<<"PQI LF Encoding bpp = "<<(double)stat.footbitLength/(double)ensemble.size().volumn()<<endl;
        if (metricModifier==MetricModifier::STSIM2_ADT_NAIVE)
          {
            //std::cout<<"Region Information bpp = "<<(double)3/(double)ensemble.size().volumn()<<endl;
            //logfile<<"Region Information bpp = "<<(double)3/(double)ensemble.size().volumn()<<endl;
            //totalBits+= 3*ensemble.size().volumn()/initSize.volumn();
          }
      }
    if (mode != CodingMode::CODING_MODE_TPSS)
      {
        //std::cout<<"Base Coder BPP: "<<(double)stat.baseCoderBitLength/(double)ensemble.size().volumn()<<endl;
        std::cout<<"Base Coder BPP: "<<(double)stat.baseCoderBitLength/(double)totalBasePel<<endl;
        std::cout<<"Original JPEG block ct: "<<stat.jBlockNum[0]<<" with BPP: "<< double(stat.jBits[0])/double(stat.jBlockNum[0]*8*8) << " => "<<double(stat.jBlockNum[0]*8*8*100)/double(ensemble.size().area())<<"% of image"<<endl;
        //logfile<<"Base Coder BPP: "<<(double)stat.baseCoderBitLength/(double)ensemble.size().volumn()<<endl;
        logfile<<"Base Coder BPP: "<<(double)stat.baseCoderBitLength/(double)totalBasePel<<endl;
        logfile<<"Original JPEG block ct: "<<stat.jBlockNum[0]<<" with BPP: "<< double(stat.jBits[0])/double(stat.jBlockNum[0]*8*8) << " => "<<double(stat.jBlockNum[0]*8*8*100)/double(ensemble.size().area())<<"% of image"<<endl;

        if (jpegType == JPEGType::JPEG_ADAPTIVE)
          {
            std::cout<<"Degrade JPEG block ct: "<<stat.jBlockNum[1]<<" with BPP: "<< double(stat.jBits[1])/double(stat.jBlockNum[1]*8*8) << "=> " <<double(stat.jBlockNum[1]*8*8*100)/double(ensemble.size().area())<<"% of image"<<endl;
            logfile<<"Degrade JPEG block ct: "<<stat.jBlockNum[1]<<" with BPP: "<< double(stat.jBits[1])/double(stat.jBlockNum[1]*8*8) << "=> " <<double(stat.jBlockNum[1]*8*8*100)/double(ensemble.size().area())<<"% of image"<<endl;
            std::cout<<"Bpp for indicate adaptive jpeg: "<<double(stat.adaptiveJpegBitLength)/double((stat.jBlockNum[0]+stat.jBlockNum[1])*64)<<endl;///bpp here
            logfile<<"Bits for indicate adaptive jpeg: "<<stat.adaptiveJpegBitLength/double((stat.jBlockNum[0]+stat.jBlockNum[1])*64)<<endl;
          }
      } //endif mode != TPSS

    std::cout<<"Tree Bit BPP: "<<(double)stat.treebitLength/(double)ensemble.size().volumn()<<endl;
    logfile<<"Tree Bit BPP: "<<(double)stat.treebitLength/(double)ensemble.size().volumn()<<endl;
    //std::cout<<"32x32 block ct: "<<stat.tpBlockNum[1]<<" with BPP: "<<double(stat.tpBits[1])/double(stat.tpBlockNum[1]*32*32) << "=> "<<double(stat.tpBlockNum[1]*32*32*100)/double(ensemble.size().area())<<"% of image"<<endl;
    //std::cout<<"16x16 block ct: "<<stat.tpBlockNum[0]<<" with BPP: "<<double(stat.tpBits[0])/double(stat.tpBlockNum[0]*16*16) << "=> "<<double(stat.tpBlockNum[0]*16*16*100)/double(ensemble.size().area())<<"% of image"<<endl;

    for (int l=blockSize.height; l >8; l/=2)
      {

        std::cout<<l<<"x"<<l<<" block ct: "<<stat.tpBlockNum[l>>5]<<" with BPP: "<<double(stat.tpBits[l>>5])/double(stat.tpBlockNum[l>>5]*l*l) << "=> "<<double(stat.tpBlockNum[l>>5]*l*l*100)/double(ensemble.size().area())<<"% of image"<<endl;
        logfile<<l<<"x"<<l<<" block ct: "<<stat.tpBlockNum[l>>5]<<" with BPP: "<<double(stat.tpBits[l>>5])/double(stat.tpBlockNum[l>>5]*l*l) << "=> "<<double(stat.tpBlockNum[l>>5]*l*l*100)/double(ensemble.size().area())<<"% of image"<<endl;
      }
    //logfile<<"32x32 block ct: "<<stat.tpBlockNum[1]<<" with BPP: "<<double(stat.tpBits[1])/double(stat.tpBlockNum[1]*32*32) << "=> "<<double(stat.tpBlockNum[1]*32*32*100)/double(ensemble.size().area())<<"% of image"<<endl;
    //logfile<<"16x16 block ct: "<<stat.tpBlockNum[0]<<" with BPP: "<<double(stat.tpBits[0])/double(stat.tpBlockNum[0]*16*16) << "=> "<<double(stat.tpBlockNum[0]*16*16*100)/double(ensemble.size().area())<<"% of image"<<endl;

    if (matching_method==MatchingMethod::MATCHING_VAR)
      {
        std::cout<<"Motion vector BPP"<< double(20*(stat.tpBlockNum[0]+stat.tpBlockNum[1]))/double(16*16*stat.tpBlockNum[0]+32*32*stat.tpBlockNum[1])<<endl;
        logfile<<"Motion vector BPP"<< double(20*(stat.tpBlockNum[0]+stat.tpBlockNum[1]))/double(16*16*stat.tpBlockNum[0]+32*32*stat.tpBlockNum[1])<<endl;
        totalBits+=20*(stat.tpBlockNum[0]+stat.tpBlockNum[1]);
      }
    CV_Assert((unsigned int)totalBits == stat.outputBit.length());
    std::cout<<"Total BPP: "<< double(totalBits)/double(ensemble.size().volumn())<<endl;
    logfile<<"Total BPP: "<< double(totalBits)/double(ensemble.size().volumn())<<endl;
    std::cout<<"Encoding finished!"<<endl;
    logfile<<"Encoding finished!"<<endl;
    logfile<<"================================================="<<endl;
#if defined(_USE_DATABASE)
    try
    {
      //remember add "space" after each entry
      SAString dbstr = "update ";
      dbstr+=DB_NAME;
      dbstr+=" set [End Time]=GETDATE()";
      dbstr+=" ,[Code Mode] = :3";
      dbstr+=" ,[STSIM2 Base Thrd] = :4";
      dbstr+=" ,[STSIM2 Thrd Multiplier] = :5";
      dbstr+=",[JPEG Type] = :6";
      dbstr+=" ,[JPEG Base Quality] =:7";
      dbstr+=" ,[JPEG Degrade Quality] = :8";
      dbstr+=" ,[PB Type] = :9";
      dbstr+=",[LC Type] = :10";
      dbstr+=" ,[Side Matching Type] = :11";
      dbstr+=" ,[Var Matching Thrd ]= :12";
      dbstr+=" ,[Candidate Number] = :13";
      dbstr+=" ,[Init Block Size Height] = :14h";
      dbstr+=" ,[Init Block Size Width] = :14w";
      dbstr+=" ,[Init Block Size Depth]= :14d";
      dbstr+=" ,[Init Q Size] = :15";
      dbstr+=" ,[Search Step X]= :16x";
      dbstr+=" ,[Search Step Y]= :16y";
      dbstr+=" ,[Search Step Z]= :16z";
      dbstr+=" ,[Overlap Ratio] = :17";
      dbstr+=" ,[STSIM2 SubWin Size]= :18";
      dbstr+=" ,[STSIM2 SubWin Step]= :19";
      dbstr+=" ,[STSIM2 Pooling Type]= :20";
      dbstr+=" ,[Metric Modifier Type] = :21";
      dbstr+=" ,[LF BPP] = :22";
      //dbstr+=" ,[Region Detector BPP]= :23";//not necessary
      dbstr+=" ,[Base Coder BPP] = :24";
      dbstr+=" ,[Base JPEG %]= :25";
      dbstr+=" ,[Base JPEG BPP] = :26";
      dbstr+=" ,[Degrade JPEG %]= :27";
      dbstr+=" ,[Degrade JPEG BPP]= :28";
      dbstr+=" ,[Tree BPP]= :29";
      dbstr+=" ,[32x32 block %]= :30";
      dbstr+=" ,[32x32 block BPP]= :31";
      dbstr+=" ,[16x16 block %]= :32";
      dbstr+=" ,[16x16 block BPP]= :33";
      dbstr+=" ,[Total BPP]= :34";
      dbstr+=" ,[Result Image Path] = :35";
      dbstr+=" ,[LF Image Path]= :36L";
      //dbstr+=" ,[HF Image Path]= :36H";
      dbstr+=" ,[Seam Path]= :37";
      dbstr+=" ,Complete = :38";
      dbstr+=" ,[Test Image] = :39";
      dbstr+= " where [Begin Time]=(select max([Begin Time]) from ";
      dbstr+= DB_NAME;
      dbstr+=")";

      cmd.setCommandText(dbstr);
      cmd.Execute();

      cmd.Param("3").setAsShort()=this->mode;
      cmd.Param("4").setAsDouble()=this->qualityThrd;
      cmd.Param("5").setAsDouble()=this->qfactor;
      cmd.Param("6").setAsShort() = this->jpegType;
      cmd.Param("7").setAsShort() = this->jpegQNum;
      cmd.Param("8").setAsShort() = short(double(this->jpegQNum)*0.4);
      cmd.Param("9").setAsShort() = this->postBlendType;
      cmd.Param("10").setAsShort() = this->lightCorrectionType;
      cmd.Param("11").setAsShort() = this->matching_method;
      cmd.Param("12").setAsDouble() = this->varThrd1;
      cmd.Param("13").setAsShort() = this->candidNum;
      cmd.Param("14h").setAsShort() = this->initSize.height;
      cmd.Param("14w").setAsShort() = this->initSize.width;
      cmd.Param("14d").setAsShort() = this->initSize.depth;
      cmd.Param("15").setAsDouble() = this->initQSize;
      cmd.Param("16x").setAsShort() = this->searchStep.height;
      cmd.Param("16y").setAsShort() = this->searchStep.width;
      cmd.Param("16z").setAsShort() = this->searchStep.depth;
      cmd.Param("17").setAsDouble() = double(this->initSize.height)/double(this->overlapSize.height);
      cmd.Param("18").setAsShort() = this->this->GetSubWinSize().height;
      cmd.Param("19").setAsShort() = this->this->GetSubWinStep().height;
      cmd.Param("20").setAsShort() = this->stsim2PoolType;
      cmd.Param("21").setAsShort() = this->metricModifier;
      cmd.Param("22").setAsDouble() = (double)stat.footbitLength/(double)ensemble.size().volumn();
      //cmd.Param("23").setAsDouble() = //region bpp
      cmd.Param("24").setAsDouble() = (double)stat.baseCoderBitLength/(double)ensemble.size().volumn();//base coder bpp
      cmd.Param("25").setAsDouble() = double(stat.jBlockNum[0]*8*8*100)/double(ensemble.size().area());//base j %
      cmd.Param("26").setAsDouble() = double(stat.jBits[0])/double(stat.jBlockNum[0]*8*8);//base j bpp
      cmd.Param("27").setAsDouble() = double(stat.jBlockNum[1]*8*8*100)/double(ensemble.size().area());//de j %
      cmd.Param("28").setAsDouble() = double(stat.jBits[1])/double(stat.jBlockNum[1]*8*8);//de j bpp
      cmd.Param("29").setAsDouble() = (double)stat.treebitLength/(double)ensemble.size().volumn();//tree bpp
      cmd.Param("30").setAsDouble() = double(stat.tpBlockNum[1]*32*32*100)/double(ensemble.size().area());//32 %
      cmd.Param("31").setAsDouble() = (stat.tpBlockNum[1]==0)?0:double(stat.tpBits[1])/double(stat.tpBlockNum[1]*32*32);
      cmd.Param("32").setAsDouble() = double(stat.tpBlockNum[0]*16*16*100)/double(ensemble.size().area());//16 %
      cmd.Param("33").setAsDouble() = (stat.tpBlockNum[0]==0)?0:double(stat.tpBits[0])/double(stat.tpBlockNum[0]*16*16);//16 bpp
      cmd.Param("34").setAsDouble() =  double(totalBits)/double(ensemble.size().volumn());//total bpp
      //SAString tempstr = SAString((path+count_exp+"_rst.tif").c_str());
      TCHAR s[100];
      char sc[100];
      DWORD b = GetCurrentDirectory(100, s);
      SAString sas(s);
      SAString sa_num((path.substr(1,path.length()-1)+count_exp).c_str());
      cmd.Param("35").setAsString() = sas + sa_num +"_rst.tif";
      cmd.Param("36L").setAsString() = sas+sa_num +"_LFImage.tif";
      //cmd.Param("36H").setAsString() = sas+sa_num +"_HFImage.tif";
      cmd.Param("37").setAsString() = sas+sa_num+"_seam.tif";
      cmd.Param("38").setAsBool() = true;
      cmd.Param("39").setAsString() = this->cFileName.c_str();
      cmd.Execute();
      conn.Commit();
      conn.Disconnect();
      std::cout<<green<<endl<<"DB done!"<<white<<endl;
    }
    catch (SAException &x)
    {
      try
      {
        conn.Rollback();
      }
      catch(SAException &)
      {

      }
      printf("%s\n",(const char*)x.ErrText());
    }

#endif
    candPosLog.close();
    logfile.close();
    //copy matching and ssim_terms to folder
    string topath = path+count_exp+"_"+PID+"_matching.txt";
    //wstring wstemp(topath.begin(),topath.end());
    boost::filesystem::copy_file("./temp/matching.txt", topath);
    //CopyFile(L"matching.txt",wstemp.c_str(),0);
    ifstream my_file("./temp/ssim_terms.txt");
    if (my_file.good())
      {
        topath = path+count_exp+"_"+PID+"_ssim_terms.txt";
        //wstemp = wstring(topath.begin(),topath.end());
        boost::filesystem::copy_file("./temp/ssim_terms.txt",topath);
        //CopyFile(L"ssim_terms.txt",wstemp.c_str(),0);
      }
  }
  Tensor<MTC::T,MTC::cn>& MTC::LoadPreDefLighting(void)
  {
    unsigned found = this->cFileName.find_last_of(".");
    string filepath = this->cFileName.substr(0,found);
    LFImage.Load(filepath+"_smooth.tif");
    LFImage.SaveBlock(path+count_exp+"_LFImage.tif");
    (ensemble - LFImage+128).SaveBlock(path+count_exp+"_HFImage.tif");
    //comment out on Jan 7, nolonger using only HF
    ensemble = (ensemble - LFImage); //get the HF
    std::cout<<white<<"Load Predef lighting done! bpp = N/A"<<endl;
    LFImage.Display(3000,1);
    //this->ReInitGrid();
    return LFImage;
  }
  Tensor<MTC::T,MTC::cn>& MTC::PQICodingLFComponent(void)
  {
    LFImage = ensemble.Clone();
    QNode<T,cn> qNode,tempAve;
    QTree<T,cn> qTree,temp;
    QTree<T,cn>* tempTree;
    int count = 0;
    double qsize = 10;
    Tensor<T,cn> tempMx(ensemble.size());
    //compute the average of each 16x16 center
    for (int x=0; x<ensemble.size().height-2*nPQILFSize+1; x+=nPQILFSize)
      for (int y=0; y<ensemble.size().width-2*nPQILFSize+1; y+=nPQILFSize)
        for (int t=0; t<ensemble.size().depth; t++)
          {
            temp = ensemble(Cube(x,y,t,2*nPQILFSize,2*nPQILFSize,1));
            LFImage(x+nPQILFSize-1,y+nPQILFSize-1,t) = temp.Mean();
            //tempMx(x+nPQILFSize,y+nPQILFSize,t) = temp.Mean();
            if (unitTest)
              {
                if (x==0&&y==0)
                  nt.ExpectNR(LFImage(x+7,y+7,t)[0],94.2539,"Mean vaule compute")<<"("<<x<<","<<y<<")";
                if (x==280 && y== 296)
                  nt.ExpectNR(LFImage(x+7,y+7,t)[0],180.7734,"Mean vaule compute")<<"("<<x<<","<<y<<")";
              }
            //ensemble(Cube(x,y,t,16,16,1)).Print();
            //LFImage(Cube(x,y,t,16,16,1)).Print();
          }
    //tempMx.Display();
    LFImage.SaveBlock(".\\LFImage_ave.tif");
    for (int x=0; x<GetGridSize().height;x++)
      for (int y=0; y<GetGridSize().width; y++)
        for (int t=0; t<GetGridSize().depth; t++)
          {
            qNode = GetNode(Point3i(x,y,t));
            temp.Ref(LFImage,Cube(Point3i(x*blockSize.height,y*blockSize.width,t*blockSize.depth),qNode.size()),qNode.overlap());
            tempTree = &temp;
            //std::cout<<x<<","<<y<<endl;

            while (tempTree != NULL)
              {
                //count++;
                if (tempTree->size().height> nPQILFSize && tempTree->size().width> nPQILFSize)
                  {
                    if (tempTree->size().height == 2*nPQILFSize) //16x16, compute 3 special feet
                      {
                        temp.Ref(LFImage,Cube(tempTree->offset(),tempTree->size()),Size3(1,1,0));
                        //the true values of the feet are the average around 16x16 block
                        temp.InitBound(qsize);//prepare boundary when reach the image bound
                        temp.ComputeFoot(qsize,Directions::DIRECTION_OTHER,SrcCodingMethod::UNARY_CODE);  // compute foot at (w,h)
                        stat.outputBit+=temp.GetBits();

                        stat.footBitInLevel[2]+=temp.BitLength();
                        stat.footbitLength+=temp.BitLength();

                        temp.ComputeFoot(qsize,Directions::DIRECTION_HORIZONTAL,SrcCodingMethod::UNARY_CODE); // compute foot at (w/2,h)
                        stat.outputBit+=temp.GetBits();
                        stat.footBitInLevel[2]+=temp.BitLength();
                        stat.footbitLength+=temp.BitLength();

                        temp.ComputeFoot(qsize,Directions::DIRECTION_VERTICAL,SrcCodingMethod::UNARY_CODE);//compute foot at (w,h/2)
                        stat.outputBit+=temp.GetBits();
                        stat.footbitLength+=temp.BitLength();
                        stat.footBitInLevel[2]+= temp.BitLength();
                        temp.ComputeFoot(qsize,Directions::DIRECTION_CENTER,SrcCodingMethod::UNARY_CODE);
                        stat.outputBit+=temp.GetBits();
                        stat.footBitInLevel[2]+=temp.BitLength();
                        stat.footbitLength+=temp.BitLength();

                      }

                    tempTree->Split();
                    //tempTree = tempTree->NextLeaf();
                    tempTree = &*tempTree->GetSubTree(0);
                  }
                else
                  {
                    //if (tempTree->size().height == nPQILFSize)
                    //{
                    temp.Ref(LFImage,Cube(tempTree->offset(),tempTree->size()),Size3(1,1,0));
                    temp.InitBound(qsize);
                    temp.ComputeFoot(qsize,Directions::DIRECTION_OTHER,SrcCodingMethod::UNARY_CODE);
                    stat.outputBit+= temp.GetBits();
                    stat.footBitInLevel[1]+=temp.BitLength();
                    stat.footbitLength+= temp.BitLength();
                    temp.LinearInterp(qsize,false);
                    tempTree = tempTree->NextLeaf();
                    //}

                    //temp.Ref(LFImage,Cube(tempTree->offset(),tempTree->size()),Size3(1,1,0));
                    //	//temp.upBound.Print();
                    //	//temp.leftBound.Print();
                    //	//temp.Print();
                    //temp.LinearInterp(qsize,false);//use quantization size 10
                    //	//temp.upBound.Print();
                    //	//temp.leftBound.Print();
                    //	//temp.Print();NextLeaf();
                  }
              }
          }
    //LFImage.Display();
    if (unitTest)
      {
        //nt.ExpectNR(LFImage(0,0,0)[0],126.375,"Test PQI LF.")<<"@(0,0)";
        nt.ExpectNR(LFImage(394,55,0)[0],132.3221,"Test PQI LF.")<<"@(394,55)";
        nt.ExpectNR(LFImage(363,354,0)[0],117.7321,"Test PQI LF.")<<"@(363,354)\n";
      }
    LFImage.SaveBlock(path+count_exp+"_LFImage.tif");
    (ensemble - LFImage+128).SaveBlock(path+count_exp+"_HFImage.tif");
    //comment out on Jan 7, nolonger using only HF
    ensemble = (ensemble - LFImage); //get the HF
    std::cout<<white<<"PQI LF Encoding done! bpp = "<<(double)stat.footbitLength/(double)ensemble.size().volumn()<<endl;
    if (unitTest)
      {
        //		nt.ExpectNR(ensemble(394,55,0)[0],3.6779,"Test HP Image.")<<"@(394,55)";
        //		nt.ExpectNR(ensemble(101,280,0)[0],44.7125,"Test HP Image.")<<"@(101,280)";
        nt.ExpectNR((double)stat.footbitLength/(double)ensemble.size().volumn(),0.034,"Test LF BPP.\n\n");
      }
    this->totalBits += stat.footbitLength;
    ////comput 8x8 dct store in dctEnsemble
    //for (int i= 0; i < ensemble.size().height; i+=8)
    //{
    //	for (int j=0; j< ensemble.size().width; j+=8)
    //	{
    //		QNode<T,cn> org = ensemble.Crop(Point3i(i,j,0),Size3(8,8,1));
    //		//org.Print();
    //		/////////// if using Watson quan table , convert to YCbCr first
    //		Tensor<float,cn> floatVerNode(org.size());
    //		//floatVerNode.Print();
    //		if (jpegType == JPEG_WASTON)
    //		{
    //			Tensor<float,3> colorVerNode = org.CvtColor<3>(CV_GRAY2RGB);
    //			colorVerNode = colorVerNode.CvtColor<3>(CV_RGB2YCrCb);
    //			//colorVerNode.Print();
    //			vector<Mat> planes;
    //			cv::split(colorVerNode[0],planes);
    //			floatVerNode[0] = planes[0]-128;
    //		}
    //		else
    //			floatVerNode = Tensor<float,cn>(org-128);
    //		//floatVerNode.Print();
    //		cv::dct(floatVerNode[0],floatVerNode[0]);
    //		//floatVerNode.Print();
    //		dctEnsemble.SetBlock(Point3i(i,j,0),Tensor<double,cn>(floatVerNode));
    //	}
    //}
    //dctEnsemble.SaveBlock(".\\dctensemble.tif");
    //ensemble.Display(3000,1);
    LFImage.Display(3000,1);
    return LFImage;
  }
  Tensor<MTC::T,MTC::cn>& MTC::ComputeEnsembleDCT(int bsize)
  {
    //ensemble.Display(3000,1);
    if (unitTest)
      {
        nt.ExpectNR(ensemble(65,266,0)[0],-14.3251,"Test HP Image.")<<"@(65,266)";
        nt.ExpectNR(ensemble(101,280,0)[0],44.7125,"Test HP Image.")<<"@(101,280)";
      }
    //comput 8x8 dct store in dctEnsemble
    for (int i= 0; i < ensemble.size().height; i+=bsize)
      {
        for (int j=0; j< ensemble.size().width; j+=bsize)
          {
            QNode<T,cn> org = ensemble.Crop(Point3i(i,j,0),Size3(bsize,bsize,1));
            //org.Print();
            /////////// if using Watson quan table , convert to YCbCr first
            Tensor<float,cn> floatVerNode(org.size());
            //floatVerNode.Print();
            if (jpegType == JPEGType::JPEG_WASTON)
              {
#if CV_MINOR_VERSION < 5
		Tensor<float,3> colorVerNode = org.CvtColor<3>(CV_GRAY2RGB);
		colorVerNode = colorVerNode.CvtColor<3>(CV_RGB2YCrCb);
#else
		Tensor<float,3> colorVerNode = org.CvtColor<3>(COLOR_GRAY2RGB);
		colorVerNode = colorVerNode.CvtColor<3>(COLOR_RGB2YCrCb);
#endif

		//colorVerNode.Print();
		vector<Mat> planes;
		cv::split(colorVerNode[0],planes);
		floatVerNode[0] = planes[0];
	      }
	    else
	      floatVerNode = Tensor<float,cn>(org);
	    //floatVerNode.Print();
	    cv::dct(floatVerNode[0],floatVerNode[0]);
	    //cv::idct(floatVerNode[0],floatVerNode[0],DFT_SCALE);

	    //floatVerNode.Print();
	    dctEnsemble.SetBlock(Point3i(i,j,0),Tensor<double,cn>(floatVerNode));
	  }
      }
    if (unitTest)
      {
        nt.ExpectNR(dctEnsemble(65,266,0)[0],10.03,"Test DCT.")<<"@(65,266)";
      }
    return dctEnsemble;
  }


  //tempTree = tempTree->
  QTree<MTC::T,MTC::cn>& MTC::CodingCore(QTree<T,cn>& qTree)
  {
    //qNode.Display();
    //pre coding
    /*
  if (qTree.upBound.Var()[0]<50&&qTree.leftBound.Var()[0]<50&&qTree.size().height>=8)//QNode<T,cn>(ensemble,qTree.size(),qTree.offset(), qTree.overlap()).Var()[0]<50 && qTree.size().height>=8) //do PQI
  {
    //qTree = QNode<T,cn>(rst,qTree.size(),qTree.offset(), qTree.overlap());
    Tensor<T,cn> org = ensemble(Cube(qTree.offset(),qTree.size()));
    //qNode.SetFoot(ensemble);
    qTree(qTree.size().height-1,qTree.size().width-1,0) = ensemble[qTree.GetFootPos()];
    //qNode.SetBlock(ensemble(Cube(qNode.offset(),qNode.size())));
    //qNode.GetExtendTensor().Display();
    //rst.SaveBlock("rst_before_PQI.tiff");
    qTree.LinearInterp(qTree.size().height);

    //ensemble.SaveBlock("org_after_PQI.tiff");
    //qNode.Print();
    //rst.Display();
    //rst.SaveBlock("rst.tiff");
    if(qTree.Compare(org,COMPARE_CRITERIA_INTERP,thrd[(int)log(double(blockSize.height)/double(qTree.size().height))/log(2.0)])>0)
    {
      qTree.bitcount.SetCodeMethod(CODING_PQI);
      qTree.bitcount.AddCodeBit(   qTree.GetBits());
      cout<<"BPP of PQI block:"<<qTree.offset()<<" "<<double(qTree.GetBits().length())/double(qTree.size().area())<<endl;
      return qTree;
    }
    else
    {
      qTree.SetBlock(Tensor<T,cn>(qTree.size(),0)); //reset to 0;
      //rst.Display();
    }
  }
  */
    if (!TexturePrediction(qTree)&&!AdaptiveTPSS(qTree,0))
      {
        //qNode.Display();
        if (qTree.size().height==8 && qTree.size().width==8)
          {
            ComputeJPEG(qTree);
            //rst.Display();
            //std::cout<<qNode.offset()<<endl;
            //rst(Cube(qNode.offset(),qNode.size()+Size3(2,2,0))).Print();
            //this->UpdateRstSqr(Cube(qNode.offset(),qNode.size()));
            //this->SAT.SetBlock(qNode.offset()+Point3i(1,1,0),rstSqr(Cube(qNode.offset(),qNode.size())));
            //this->SAT = rstSqr.ComputeSAT(this->SAT,qNode.offset(),qNode.offset()+qNode.size().Point3()-Point3i(1,1,1));
            Point3i pos = qTree.offset();
            Size3 sz = qTree.size();
            UpdateSAT(pos,sz);
            //rstSqr(Cube(qNode.offset(),qNode.size()+Size3(2,2,0))).Print();
            //
            //SAT(Cube(qNode.offset(),qNode.size()+Size3(1,1,0)+Size3(2,2,0))).Print();

          }
        else
          {
            qTree.Split();
            if(this->mode == CodingMode::CODING_MODE_MTC)//split output 1bit to indicate
              {
                stat.outputBit+="1";
                outputfile<<"1";
                totalBits++;
                stat.treebitLength++;
                qTree.bitcount.SetDecBit("1");
              }
            QTree<T,cn> *tempTree = &*qTree.GetSubTree(0);
            //tempTree = qTree.NextLeaf();////

            while(true)
              {
                if(tempTree == NULL)
                  break;
                else
                  {
                    CodingCore(*tempTree);//QNode<T,cn>(rst,tempTree->size(),tempTree->offset(),qTree.overlap()/2));
                  }
                if(tempTree->GetTreePos()==tempTree->GetPeerTreeSize()-1)
                  {
                    if (postBlendType == PostBlendingType::POST_BLENDING_ONLINE&&this->mode == CodingMode::CODING_MODE_MTC)
                      {

                        PBSet::iterator it = pbSet.find(PBRecord(qTree.offset(),BoundDir::LEFT));
                        if (it!=pbSet.end()&&it->second.height <= qTree.size().height)
                          {
                            cout<<"offset is: "<<qTree.offset()<<endl;
                            cout<<"pb offset is: "<<it->first.offset<<"which  points to: "<<it->second.offset()<<endl;
                            PostBlending(it->first,it->second);//pb for different size may goes wrong!
                            pbSet.erase(it->first	,it->second);
                            //rst_with_seam.Display(1000,1);
                            UpdateSAT(qTree.offset(),qTree.size());
                          }
                        it = pbSet.find(PBRecord(qTree.offset(),BoundDir::UP));
                        if (it!=pbSet.end()&&it->second.width <= qTree.size().width)
                          {
                            PostBlending(it->first,it->second);
                            pbSet.erase(it->first,it->second);
                            //rst_with_seam.Display(1000,1);
                            UpdateSAT(qTree.offset(),qTree.size());
                          }

                      }
                    break;
                  }
                else
                  tempTree=tempTree->NextLeaf();
              }
          }
      }

    return qTree;
  }
  void MTC::ModifyMetric(const Point3i& blockPos)
  {
    if (metricModifier==MetricModifier::STSIM2_ADT_NAIVE)
      {
        switch (edgeMap[blockPos.x][blockPos.y])
          {
          case 0:
            qualityThrd = qfactor * orgQualityThrd;
            break;
          case 1:
            qualityThrd = qfactor * orgQualityThrd;
            break;
          case 2:
            qualityThrd = orgQualityThrd;
            break;
          case 3:
            qualityThrd = qfactor*qfactor * orgQualityThrd;
            break;
          case 4:
            qualityThrd = orgQualityThrd;
            break;
          case 5:
            qualityThrd = qfactor * orgQualityThrd;
            break;
          default:
            qualityThrd = orgQualityThrd;
            break;
          }
      }
    else
      {
        qualityThrd = orgQualityThrd;
      }

  }

  QNode<MTC::T,MTC::cn> MTC::GetValidCandid(Point3i matchPos, const QNode<T,cn>& qNode)
  {
    QNode<T,cn> candid;
    if (matchPos.x + qNode.overlap().height - qNode.size().height/2 < 0 ||
        matchPos.y + qNode.overlap().width - qNode.size().width/2 <0 ||
        matchPos.x + qNode.overlap().height + 3 * qNode.size().height/2 > ensemble.size().height ||
        matchPos.y + qNode.overlap().width + 3*qNode.size().width/2 > ensemble.size().width)
      {
        Point3i pos1,pos2;
        pos1 = Point3i(matchPos.x + qNode.overlap().height - qNode.size().height/2, matchPos.y + qNode.overlap().width - qNode.size().width/2 , 0);
        Size3 sz(qNode.size()*Size3(2,2,1));
        Tensor<T,cn> temp1(sz,0);
        //temp1.SetSubWinSize(qNode.GetSubWinSize());
        //temp1.SetSubWinStep(qNode.GetSubWinStep());
        if (matchPos.x+qNode.overlap().height - qNode.size().height/2 < 0 )
          {
            pos1.x=0;
            pos2.x= qNode.size().height/2 -( matchPos.x+qNode.overlap().height );
            sz.height -= pos2.x;
          }
        if (	matchPos.y + qNode.overlap().width - qNode.size().width/2 <0 )
          {
            pos1.y=0;
            pos2.y= qNode.size().width/2 - (matchPos.y + qNode.overlap().width);
            sz.width -= pos2.y;
          }
        if (matchPos.x + qNode.overlap().height + 3 * qNode.size().height/2 > ensemble.size().height)
          {
            pos1.x=matchPos.x + qNode.overlap().height ;
            pos2.x=0;
            sz.height-= qNode.size().height/2 - (matchPos.x + qNode.overlap().height + 3 * qNode.size().height/2 - rst.size().height);
          }
        if (matchPos.y + qNode.overlap().width + 3*qNode.size().width/2 > ensemble.size().width)
          {
            pos1.y=matchPos.y + qNode.overlap().width;
            pos2.y=0;
            sz.width-= qNode.size().width/2 - (matchPos.y + qNode.overlap().width + 3*qNode.size().width/2 -  rst.size().width);
          }
        Tensor<T,cn> temp2 = rst.Crop( pos1, sz);
        //temp2.Print();
        temp1.SetBlock(pos2,temp2);
        Tensor<T,cn> reflect_temp,reflect;
        if (sz.height<temp1.size().height )//&& sz.width==temp1.size().width)
          {
            if (pos2.x>0) // empty in upper border
              {
                reflect_temp=temp2.Crop(Point3i(0,0,0),Size3(pos2.x,sz.width,sz.depth));
                //reflect_temp.Print();
                //cv::flip(reflect_temp,reflect,0);//vertical flip
                reflect = reflect_temp.Flip(1,0);
                //reflect.Print();
                temp1.SetBlock(Point3i(0,0,0),reflect);
                /*temp2.Print();
                                reflect_temp.Print();
                                reflect.Print();
                                temp1.Print();*/
              }
            else if (pos1.x>0) //empty in lower border
              {
                reflect_temp=temp2.Crop(Point3i(sz.height-pos1.x,0,0),Size3(pos1.x,sz.width,sz.depth));
                //cv::flip(reflect_temp,reflect,0);//vertical flip
                reflect= reflect_temp.Flip(1,0);
                temp1.SetBlock(Point3i(sz.height,0,0),reflect);
              }
            else
              {
                std::cout<<"wrong!!!!!!!!!!!!!!!!!!! in reflect border\n";
                CV_DbgAssert(false);
              }

          }

        if (sz.width < temp1.size().width )// && sz.height==temp1.size().height )
          {
            if (pos2.y>0) // empty in left border
              {
                reflect_temp=temp2.Crop(Point3i(0,0,0),Size3(sz.height,pos2.y,sz.depth));
                //cv::flip(reflect_temp,reflect,1);//horizontal flip
                reflect = reflect_temp.Flip(0,1);
                temp1.SetBlock(reflect);

              }
            else if (pos1.y>0) //empty in right border
              {
                reflect_temp=temp2.Crop(Point3i(0,sz.height-pos1.y,0),Size3(sz.width,pos1.y,sz.depth));
                //cv::flip(reflect_temp,reflect,1);//horizontal flip
                reflect = reflect_temp.Flip(0,1);
                temp1.SetBlock(Point3i(0,sz.width,0),reflect);
                //temp1.Print();
              }
            else
              {
                std::cout<<"wrong!!!!!!!!!!!!!!!!!!! in reflect border\n";
                CV_DbgAssert(false);
              }

          }

        candid = QNode<T,cn>(temp1,qNode.size(),(qNode.size()/2).Point3(),qNode.size()/2).Clone();
        //old version, use same boundary
        //
        //	candid = rst.Crop( qNode.offset() - (qNode.size()/2).Point3(),qNode.size()*Size3(2,2,1));
        //	Point3i pos = (qNode.size() * Size3(1,3,1) / 2).Point3();
        //	Size3 sz = qNode.size() * Size3(3,1,1) / Size3(2,2,1);
        //	candid.SetBlock(pos,ensemble.Crop( qNode.offset() - (qNode.size()/2).Point3() + pos, sz));
        ////candid.Print();
        //	sz = qNode.size()* Size3(1,3,1)/ Size3(2,2,1);
        //	pos = (qNode.size() * Size3(3,0,1) / 2).Point3();
        //	candid.SetBlock(pos,ensemble.Crop( qNode.offset() - (qNode.size()/2).Point3() + pos, sz));
        ////candid.Print();
        //	pos = (qNode.size()/2).Point3();
        //	candid.SetBlock(pos,rst.Crop( matchCandid[i] + qNode.overlap().Point3(), qNode.size()));
      }
    else
      {
        //candid.Print();
        //new version, use the true boundary of node
        //matchCandid get the position of boundary matching, which + overlap size is the candid position
        //candid = rst.Crop(matchCandid[i] + qNode.overlap().Point3() - (qNode.size()/2).Point3(),
        //				  qNode.size()*Size3(2,2,1));
        candid = QNode<T,cn>(rst,qNode.size(),matchPos+qNode.overlap().Point3(),
                             qNode.size()/2).Clone();
      }
    return candid;
  }

  QNode<MTC::T,MTC::cn> MTC::GetValidCandidSimple(Point3i matchPos, const QNode<T,cn>& qNode)
  {
    //20130605 : for get rid of extra boundary
    QNode<T,cn> candid;
    //20130605 if (matchPos.x + qNode.overlap().height - qNode.size().height/2 < 0 ||
    //	matchPos.y + qNode.overlap().width - qNode.size().width/2 <0 ||
    //	matchPos.x + qNode.overlap().height + 3 * qNode.size().height/2 > ensemble.size().height ||
    //	matchPos.y + qNode.overlap().width + 3*qNode.size().width/2 > ensemble.size().width)
    if (matchPos.x + qNode.overlap().height*2 + qNode.size().height < ensemble.size().height ||
        matchPos.y + qNode.overlap().width*2 + qNode.size().width < ensemble.size().width)
      {

        candid = QNode<T,cn>(rst,qNode.size(),matchPos+qNode.overlap().Point3(),
                             qNode.overlap()).Clone();
      }
    else
      {
        std::cout<<"wrong!!!!!!!!!!!!!!!!!!! in get valid candidate\n";
        CV_DbgAssert(false);
      }
    return candid;
  }


  QNode<MTC::T,MTC::cn> MTC::GetValidNode(Point3i matchPos, const QNode<T,cn>& qNode, const Tensor<T,cn>& ref)
  {
    //match position is the block offset - block overlap
    QNode<T,cn> candid;
    if (matchPos.x + qNode.overlap().height - qNode.size().height/2 < 0 ||
        matchPos.y + qNode.overlap().width - qNode.size().width/2 <0 ||
        matchPos.x + qNode.overlap().height + 3 * qNode.size().height/2 > ensemble.size().height ||
        matchPos.y + qNode.overlap().width + 3*qNode.size().width/2 > ensemble.size().width)
      {
        Point3i pos1,pos2,offsetAdjust;
        pos1 = Point3i(matchPos.x + qNode.overlap().height - qNode.size().height/2, matchPos.y + qNode.overlap().width - qNode.size().width/2 , 0);
        offsetAdjust = pos1;
        Size3 sz(qNode.size()*Size3(2,2,1));
        Tensor<T,cn> temp1(sz,0);
        //temp1.SetSubWinSize(qNode.GetSubWinSize());
        //temp1.SetSubWinStep(qNode.GetSubWinStep());
        if (matchPos.x+qNode.overlap().height - qNode.size().height/2 < 0 )
          {
            pos1.x=0;
            pos2.x= qNode.size().height/2 -( matchPos.x+qNode.overlap().height );
            sz.height -= pos2.x;
          }
        if (	matchPos.y + qNode.overlap().width - qNode.size().width/2 <0 )
          {
            pos1.y=0;
            pos2.y= qNode.size().width/2 - (matchPos.y + qNode.overlap().width);
            sz.width -= pos2.y;
          }
        if (matchPos.x + qNode.overlap().height + 3 * qNode.size().height/2 > ensemble.size().height)
          {
            pos1.x=matchPos.x - qNode.overlap().height ;
            pos2.x=0;
            sz.height-=  (matchPos.x + qNode.overlap().height + 3 * qNode.size().height/2 - ref.size().height);
          }
        if (matchPos.y + qNode.overlap().width + 3*qNode.size().width/2 > ensemble.size().width)
          {
            pos1.y=matchPos.y - qNode.overlap().width;
            pos2.y=0;
            sz.width-=  (matchPos.y + qNode.overlap().width + 3*qNode.size().width/2 - ref.size().width);
          }
        Tensor<T,cn> temp2 = ref.Crop( pos1, sz);
        //temp2.Print();
        temp1.SetBlock(pos2,temp2);
        Tensor<T,cn> reflect_temp,reflect;
        if (sz.height<temp1.size().height )//&& sz.width==temp1.size().width)
          {
            if (pos2.x>0) // empty in upper border
              {
                reflect_temp=temp2.Crop(Point3i(0,0,0),Size3(pos2.x,sz.width,sz.depth));
                //reflect_temp.Print();
                //cv::flip(reflect_temp,reflect,0);//vertical flip
                reflect = reflect_temp.Flip(1,0);
                //reflect.print();
                temp1.SetBlock(Point3i(0,pos2.y,pos2.z),reflect);
                /*temp2.Print();
                                reflect_temp.Print();
                                reflect.Print();
                                temp1.Print();*/
              }
            else if (pos1.x>0) //empty in lower border
              {
                reflect_temp=temp2.Crop(Point3i(qNode.size().height,0,0),Size3(temp1.size().height-sz.height,sz.width,sz.depth));
                //cv::flip(reflect_temp,reflect,0);//vertical flip
                reflect= reflect_temp.Flip(1,0);
                temp1.SetBlock(Point3i(sz.height,pos2.y,0),reflect);
              }
            else
              {
                std::cout<<"wrong!!!!!!!!!!!!!!!!!!! in reflect border\n";
                CV_DbgAssert(false);
              }

          }

        if (sz.width < temp1.size().width )// && sz.height==temp1.size().height )
          {
            if (pos2.y>0) // empty in left border
              {
                reflect_temp=temp2.Crop(Point3i(0,0,0),Size3(sz.height,pos2.y,sz.depth));
                //cv::flip(reflect_temp,reflect,1);//horizontal flip
                reflect = reflect_temp.Flip(0,1);
                temp1.SetBlock(Point3i(pos2.x,0,pos2.z),reflect);

              }
            else if (pos1.y>0) //empty in right border
              {
                reflect_temp=temp2.Crop(Point3i(0,qNode.size().width,0),Size3(sz.height,temp1.size().width-sz.width,sz.depth));
                //cv::flip(reflect_temp,reflect,1);//horizontal flip
                reflect = reflect_temp.Flip(0,1);
                temp1.SetBlock(Point3i(pos2.x,sz.width,0),reflect);
                //temp1.Print();
              }
            else
              {
                std::cout<<"wrong!!!!!!!!!!!!!!!!!!! in reflect border\n";
                CV_DbgAssert(false);
              }

          }
        if (qNode.offset().x == 0 && qNode.offset().y==992)
          temp1.Display();
        if (sz.height<temp1.size().height&& sz.width < temp1.size().width)
          {
            if (pos2.y>0 && pos2.x >0) //NW corner
              {
                reflect_temp=temp2.Crop(Point3i(0,0,0),Size3(pos2.x,pos2.y,sz.depth));
                reflect = reflect_temp.Flip(1,1);
                temp1.SetBlock(reflect);
              }
            if (pos1.y>0 && pos1.x>0) //SE corner
              {
                reflect_temp=temp2.Crop(Point3i(qNode.size().height,qNode.size().width,0),Size3(temp1.size().height-sz.height, temp1.size().width - sz.width,sz.depth));
                //cv::flip(reflect_temp,reflect,1);//horizontal flip
                reflect = reflect_temp.Flip(1,1);
                temp1.SetBlock(Point3i(sz.height,sz.width,0),reflect);
              }
            if (pos2.x>0 && pos1.y>0) //NE corner
              {
                reflect_temp=temp2.Crop(Point3i(0,qNode.size().width,0),Size3(temp1.size().height - sz.height,temp1.size().width-sz.width,sz.depth));
                //cv::flip(reflect_temp,reflect,1);//horizontal flip
                reflect = reflect_temp.Flip(1,1);
                temp1.SetBlock(Point3i(0,sz.width,0),reflect);
              }
            if (pos2.y>0 && pos1.x>0) //SW corner
              {
                reflect_temp=temp2.Crop(Point3i(qNode.size().height,0,0),Size3(temp1.size().height-sz.height, temp1.size().width-sz.width,sz.depth));
                //cv::flip(reflect_temp,reflect,1);//horizontal flip
                reflect = reflect_temp.Flip(1,1);
                temp1.SetBlock(Point3i(sz.height,0,0),reflect);
              }
          }
        candid = QNode<T,cn>(temp1,qNode.size(),(qNode.size()/2).Point3(),qNode.size()/2).Clone();
        candid.SetOffset(candid.offset()+offsetAdjust);
        //old version, use same boundary
        //
        //	candid = rst.Crop( qNode.offset() - (qNode.size()/2).Point3(),qNode.size()*Size3(2,2,1));
        //	Point3i pos = (qNode.size() * Size3(1,3,1) / 2).Point3();
        //	Size3 sz = qNode.size() * Size3(3,1,1) / Size3(2,2,1);
        //	candid.SetBlock(pos,ensemble.Crop( qNode.offset() - (qNode.size()/2).Point3() + pos, sz));
        ////candid.Print();
        //	sz = qNode.size()* Size3(1,3,1)/ Size3(2,2,1);
        //	pos = (qNode.size() * Size3(3,0,1) / 2).Point3();
        //	candid.SetBlock(pos,ensemble.Crop( qNode.offset() - (qNode.size()/2).Point3() + pos, sz));
        ////candid.Print();
        //	pos = (qNode.size()/2).Point3();
        //	candid.SetBlock(pos,rst.Crop( matchCandid[i] + qNode.overlap().Point3(), qNode.size()));
      }
    else
      {
        //candid.Print();
        //new version, use the true boundary of node
        //matchCandid get the position of boundary matching, which + overlap size is the candid position
        //candid = rst.Crop(matchCandid[i] + qNode.overlap().Point3() - (qNode.size()/2).Point3(),
        //				  qNode.size()*Size3(2,2,1));
        candid = QNode<T,cn>(ref,qNode.size(),matchPos+qNode.overlap().Point3(),
                             qNode.size()/2).Clone();
      }
    // candid.GetExtendTensor().Display();
    return candid;
  }
  void MTC::UpdateSAT(const Point3i& sPos, const Size3& sz)
  {
    //if (sPos == Point3i(12,92,0))
    //	rst(Cube(sPos,sz)).Print();
    Cube temp(sPos,sz);
    this->UpdateRstSqr(temp);
    //if (sPos == Point3i(12,92,0))
    //{
    //	rstSqr(Cube(sPos,sz)).Print();
    //		this->SAT(Cube(sPos,sz+Size3(1,1,0))).Print();
    //	}
    this->SAT.SetBlock(sPos+Point3i(1,1,0),rstSqr(Cube(sPos,sz)));
    //	if (sPos == Point3i(12,92,0))
    //		this->SAT(Cube(sPos,sz+Size3(1,1,0))).Print();
    this->SAT = lighting::ComputeSAT(this->SAT,sPos,sPos+sz.Point3()-Point3i(1,1,1));
    //	if (sPos == Point3i(12,92,0))
    //		this->SAT(Cube(sPos,sz+Size3(1,1,0))).Print();
    return;
  }

  std::map<cv::Point3i, FootItem, ComparePoint3i>::iterator MTC::RetrieveFootsNew(QNode<double,1>& qNode, int level)
  {
    //notice: the input qNode size has boundary size half , (different to the quarter size)
    //since the getvalidNode  (or getvalidblock) give that size boundary
    std::lock_guard<std::mutex> lock(foot_lock);
    vector<Point3i> feetPos;
    Point3i upFootPos = qNode.offset()+Point3i(-qNode.overlap().height/2,qNode.overlap().width/2+qNode.size().width,0)-Point3i(1,1,0);
    //    Point3i(qNode.size().height,0,0)+Point3i(-qNode.overlap().height+1,qNode.overlap().width,0);
    Point3i leftFootPos = qNode.offset()+Point3i(qNode.overlap().height/2+qNode.size().height,-qNode.overlap().width/2,0)-Point3i(1,1,0);
    for (int i=0; i< qNode.size().height+qNode.overlap().height; i+= (qNode.size().height+qNode.overlap().height)/(1<<level))
      {
        // 0 : bsize/1
        // 1: bsize/2
        // 2: bsize/4
        feetPos.push_back(upFootPos + Point3i(i,0,0));
      }
    for (int j=0; j< qNode.size().width+qNode.overlap().width; j+= (qNode.size().width+qNode.overlap().width)/(1<<level))
      feetPos.push_back(leftFootPos + Point3i(0,j,0));
    feetPos.push_back(qNode.GetFootPos()+(qNode.overlap()/2).Point3());
    //foot from former coded block can also be usefull
    Point3i footFromUpBlock = upFootPos+Point3i(qNode.overlap().height/2,0,0);
    Point3i footFromLeftBlock = leftFootPos + Point3i(0,qNode.overlap().width/2,0);
    feetPos.push_back(footFromUpBlock);
    feetPos.push_back(footFromLeftBlock);
    //UINT8 ret=0;
    auto upFoot = FootTable.find(upFootPos);
    auto leftFoot = FootTable.find(leftFootPos);
    auto foot = FootTable.find(qNode.GetFootPos());
    if (upFoot!= FootTable.end())
      {
        qNode.SetUpFoot(upFoot->second.value);
        //ret +=1;
      }
    else
      qNode.SetUpFoot(-500);
    if (leftFoot!= FootTable.end())
      {
        qNode.SetLeftFoot(leftFoot->second.value);
        //ret +=2;
      }
    else
      qNode.SetLeftFoot(-500);
    if (foot!= FootTable.end())
      {
        qNode.SetFoot(foot->second.value);
        //ret +=4;

      }
    else
      qNode.SetFoot(-500);
    qNode.clearFoot();
    for (auto& f : feetPos)
      {
        auto temp = FootTable.find(f);

        FootItem footitem =FootItem();
        if (temp!=FootTable.end())
          {
            footitem.value = temp->second.value;
          }
        else
          {
            footitem.value = -500;
          }
        qNode.AddFoot(pair<Point3i,FootItem>(Point3i(f.x,f.y,0),footitem));

      }

    return foot;
  }
  vector<pair<cv::Point3i, FootItem>>  MTC::RetrieveFeet(QNode<double,1>& qNode,int level)
  {
                                       std::lock_guard<std::mutex> lock(foot_lock);
                                       vector<Point3i> feetPos;
                                       Point3i upFootPos = qNode.GetFootPos()-Point3i(qNode.size().height,0,0);
                                       Point3i leftFootPos = qNode.GetFootPos()-Point3i(0,qNode.size().width,0);
                                       Point3i antiFootPos = qNode.offset()-Point3i(1,1,0);//the upper-left corner
                                       vector<Point3i> feetPosTobeCode;
                                       vector<pair<cv::Point3i, FootItem>>  feet;
                                       for (int i=0; i< qNode.size().height; i+= qNode.size().height/(1<<level))
  {
    // 0 : bsize/1
    // 1: bsize/2
    // 2: bsize/4
    if (i!=0)
      feetPosTobeCode.push_back(upFootPos + Point3i(i,0,0));
    feetPos.push_back(upFootPos + Point3i(i,0,0));
  }
  for (int j=0; j< qNode.size().width; j+= qNode.size().width/(1<<level))
  {
    if (j!=0)
      feetPosTobeCode.push_back(leftFootPos + Point3i(0,j,0));
    feetPos.push_back(leftFootPos + Point3i(0,j,0));
  }
  feetPos.push_back(qNode.GetFootPos());
  feetPosTobeCode.push_back(qNode.GetFootPos());
  //UINT8 ret=0;
  auto upFoot = FootTable.find(upFootPos);
  auto leftFoot = FootTable.find(leftFootPos);
  auto foot = FootTable.find(qNode.GetFootPos());

  if (upFoot!= FootTable.end())
  {
    qNode.SetUpFoot(upFoot->second.value);
    //ret +=1;
  }
  else
  qNode.SetUpFoot(-500);
  if (leftFoot!= FootTable.end())
  {
    qNode.SetLeftFoot(leftFoot->second.value);
    //ret +=2;
  }
  else
  qNode.SetLeftFoot(-500);
  if (foot!= FootTable.end())
  {
    qNode.SetFoot(foot->second.value);
    //ret +=4;

  }
  else
  qNode.SetFoot(-500);
  qNode.clearFoot();
  for (auto& f : feetPos)
  {
    auto temp = FootTable.find(f);
    FootItem footitem =FootItem();
    if (temp!=FootTable.end())
      {
        footitem.value = temp->second.value;
      }
    else
      {
        footitem.value = -500;
      }
    qNode.AddFoot(pair<Point3i,FootItem>(Point3i(f.x,f.y,0),footitem));

  }
  for (auto& f : feetPosTobeCode)
  {
    auto temp = FootTable.find(f);
    FootItem footitem =FootItem();
    if (temp!=FootTable.end())
      {
        footitem.value = temp->second.value;
      }
    else
      {
        footitem.value = -500;
      }
    //feet.push_back(FootTable.find(f));
    feet.push_back(pair<Point3i,FootItem>(f,footitem));
  }
  return feet;
}
std::map<cv::Point3i, FootItem, ComparePoint3i>::iterator MTC::RetrieveFoot(QNode<double,1>& qNode, int level)
{
  std::lock_guard<std::mutex> lock(foot_lock);
  if (level<0)
    return FootTable.end();
  vector<Point3i> feetPos;
  Point3i upFootPos = qNode.GetFootPos()-Point3i(qNode.size().height,0,0);
  Point3i leftFootPos = qNode.GetFootPos()-Point3i(0,qNode.size().width,0);
  Point3i antiFootPos = qNode.offset()-Point3i(1,1,0);//the upper-left corner
  for (int i=0; i< qNode.size().height; i+= qNode.size().height/(1<<level))
    {
      // 0 : bsize/1
      // 1: bsize/2
      // 2: bsize/4
      feetPos.push_back(upFootPos + Point3i(i,0,0));
    }
  for (int j=0; j< qNode.size().width; j+= qNode.size().width/(1<<level))
    {
      feetPos.push_back(leftFootPos + Point3i(0,j,0));
    }
  feetPos.push_back(qNode.GetFootPos());
  //UINT8 ret=0;
  auto upFoot = FootTable.find(upFootPos);
  auto leftFoot = FootTable.find(leftFootPos);
  auto foot = FootTable.find(qNode.GetFootPos());

  if (upFoot!= FootTable.end())
    {
      qNode.SetUpFoot(upFoot->second.value);
      //ret +=1;
    }
  else
    qNode.SetUpFoot(-500);
  if (leftFoot!= FootTable.end())
    {
      qNode.SetLeftFoot(leftFoot->second.value);
      //ret +=2;
    }
  else
    qNode.SetLeftFoot(-500);
  if (foot!= FootTable.end())
    {
      qNode.SetFoot(foot->second.value);
      //ret +=4;

    }
  else
    qNode.SetFoot(-500);
  qNode.clearFoot();
  for (auto& f : feetPos)
    {
      auto temp = FootTable.find(f);
      FootItem footitem =FootItem();
      if (temp!=FootTable.end())
        {
          footitem = temp->second;
          //footitem.value = temp->second.value;
        }
      else
        {
          footitem.value = -500;
        }
      qNode.AddFoot(pair<Point3i,FootItem>(Point3i(f.x,f.y,0),footitem));

    }

  return foot;
}
void MTC::UpdateFoots(const QNode<double,1>& org, const QNode<double,1>& cand)
{
  lock_guard<mutex> lock(foot_lock);
  FootItem foot;
  foot.value = cand.GetApproxFoot()[0];
  foot.bits = cand.GetBits().length();
  foot.bitstring = cand.GetBits();
  if (foot.value >-255) //makesure it is a valid foot
    FootTable.insert(pair<cv::Point3i, FootItem>(org.GetFootPos(),foot));
}
void MTC::UpdateFeet(const vector<pair<cv::Point3i, FootItem>>& feet)   //const QNode<double,1>& org, const QNode<double,1>& cand)
{
  lock_guard<mutex> lock(foot_lock);
  for (auto& f  : feet)
    {
      //foot.value = cand.GetApproxFoot()[0];
      //foot.bits = cand.GetBits().length();
      //foot.bitstring = cand.GetBits();
      if (f.second.value >-255 &&f.second.bits>0) //makesure it is a valid foot
        {
          auto it = FootTable.find(f.first);
          if (it == FootTable.end())
            FootTable.insert(f);
          //FootTable.insert(pair<cv::Point3i, FootItem>(org.GetFootPos(),foot));
        }
    }
}
void MTC::RemoveFeet(const vector<pair<cv::Point3i, FootItem>>& feet)   //const QNode<double,1>& org, const QNode<double,1>& cand)
{
  lock_guard<mutex> lock(foot_lock);
  for (auto& f  : feet)
    {
      //foot.value = cand.GetApproxFoot()[0];
      //foot.bits = cand.GetBits().length();
      //foot.bitstring = cand.GetBits();
      if (f.second.value >-255 &&f.second.bits>0) //makesure it is a valid foot
        {
          auto it = FootTable.find(f.first);
          FootTable.erase(it);
          //FootTable.insert(pair<cv::Point3i, FootItem>(org.GetFootPos(),foot));
        }
    }
}
QTree<MTC::T,MTC::cn>& MTC::ComputeJPEG(QTree<MTC::T,MTC::cn>& qNode)
{
  QNode<T,cn> org = ensemble.Crop(qNode.offset(),qNode.size());
  /////////// if using Watson quan table , convert to YCbCr first
  Tensor<float,cn> floatVerNode(org.size());

  ////floatVerNode.Print();
  //if (jpegQuanTblType == JPEG_WASTON)
  //{
  //	Tensor<float,3> colorVerNode = org.CvtColor<3>(CV_GRAY2RGB);
  //	colorVerNode = colorVerNode.CvtColor<3>(CV_RGB2YCrCb);
  //	//colorVerNode.Print();
  //	vector<Mat> planes;
  //	cv::split(colorVerNode[0],planes);
  //	floatVerNode[0] = planes[0]-128;
  //}
  //else
  //	floatVerNode = Tensor<float,cn>(org-128);
  ////floatVerNode.Print();
  //cv::dct(floatVerNode[0],floatVerNode[0]);
  ////cv::idct(floatVerNode[0],floatVerNode[0],DFT_SCALE);
  ////floatVerNode.Print();
  floatVerNode = dctEnsemble.Crop(qNode.offset(),floatVerNode.size()); //for the new coming block, can use pre-stored DCT table
  Tensor<T,cn> quanRst(qNode.size()); //index of quan
  //floatVerNode.Print();
  //dctEnsemble.SetBlock(qNode.offset(),Tensor<double,cn>(floatVerNode));
  //qNode.Print();
  //org.Print();
  vector<UINT16>* currentTbl;
  int high_low; //1 low , 0 high
  if (jpegType == JPEGType::JPEG_ADAPTIVE)
    {
      Tensor<T,cn> var = org.LocalVariance(org.LocalMean(qNode.size(),qNode.size()),qNode.size(),qNode.size());
      //var.Print();
      totalBits++;
      stat.adaptiveJpegBitLength++;
      if (var(0,0,0)[0]>50)
        {
          currentTbl = &JpegQuanTbl_alt;
          high_low = 1; //low
          stat.outputBit+="1";
          outputfile<<"1";
          qNode.bitcount.SetCodeMethod(CodingMethodNames::CODING_JPEG_DEGRADE);
          qNode.bitcount.SetDecBit("1");

        }
      else
        {
          currentTbl = &JpegQuanTbl;
          high_low = 0; //high
          stat.outputBit+="0";
          outputfile<<"0";
          qNode.bitcount.SetCodeMethod(CodingMethodNames::CODING_JPEG);
          qNode.bitcount.SetDecBit("0");
        }
    }
  else
    {
      currentTbl = &JpegQuanTbl;
      high_low = 0; //always high
      qNode.bitcount.SetCodeMethod(CodingMethodNames::CODING_JPEG);
    }
  for (int x=0; x < qNode.size().height;x++)
    for (int y=0; y < qNode.size().width;y++)
      {

        quanRst(x,y,0) = qNode.UniformQuantize(floatVerNode(x,y,0),(*currentTbl)[x*qNode.size().width+y]);
        floatVerNode(x,y,0) = quanRst(x,y,0).mul(Vec<T,cn>::all((*currentTbl)[x*qNode.size().width+y]));
        //quanRst(x,y,0)/JpegQuanTbl[x*qNode.size().width+y];
      }
  //quanRst.Print();
  //floatVerNode.Print();
  ////coding DC component
  ///treat different in different position, necessary recompute the dct
  Vec<T,cn> prevBlockDC;
  if (qNode.offset()==Point3i(0,0,0))
    {
      prevBlockDC = 0;
    }
  else if(qNode.offset().y==0)
    {
      //take the upper one block
      /*Tensor<T,cn> tempDCTBlk = rst.Crop(qNode.offset()-Point3i(qNode.size().height,0,0),qNode.size());
                                if (tempDCTBlk.LocalVariance(tempDCTBlk.LocalMean(qNode.size()),qNode.size())(0,0,0)[0] > 50)
                                {
                                        cv::dct(tempDCTBlk[0],tempDCTBlk[0]);
                                        prevBlockDC = qNode.UniformQuantize(tempDCTBlk(0,0,0),JpegQuanTbl_alt[0]);
                                }
                                else
                                {
                                        cv::dct(tempDCTBlk[0],tempDCTBlk[0]);
                                        prevBlockDC = qNode.UniformQuantize(tempDCTBlk(0,0,0),JpegQuanTbl[0]);
                                }*/
      prevBlockDC = qNode.UniformQuantize(dctEnsemble[qNode.offset()-Point3i(qNode.size().height,0,0)],(*currentTbl)[0]);
    }
  else if (qNode.offset().x == 0 || this->mode==CodingMode::CODING_MODE_JPEG)
    {
      prevBlockDC = qNode.UniformQuantize(dctEnsemble[qNode.offset()-Point3i(0,qNode.size().width,0)],(*currentTbl)[0]);
    }
  else
    {
      Tensor<T,cn> tempDCTBlk = rst.Crop(qNode.offset()-Point3i(0,qNode.size().width,0),qNode.size());
      if (tempDCTBlk.LocalVariance(tempDCTBlk.LocalMean(qNode.size(),qNode.size()),qNode.size(),qNode.size())(0,0,0)[0] > 50 &&  (jpegType==JPEGType::JPEG_ADAPTIVE))
        {
          cv::dct(tempDCTBlk[0],tempDCTBlk[0]);
          prevBlockDC = qNode.UniformQuantize(tempDCTBlk(0,0,0),JpegQuanTbl_alt[0]);
        }
      else
        {
          cv::dct(tempDCTBlk[0],tempDCTBlk[0]);
          prevBlockDC = qNode.UniformQuantize(tempDCTBlk(0,0,0),JpegQuanTbl[0]);
        }
      //prevBlockDC = qNode.UniformQuantize(dctEnsemble[qNode.offset()-Point3i(0,qNode.size().width,0)],(*currentTbl)[0]);
    }
  Vec<T,cn> diffDC = quanRst(0,0,0) - prevBlockDC;
  string tempstr = ComputeJpegHuff(diffDC,Huffman_DC_Map);
  qNode.bitcount.AddCodeBit(tempstr);
  int codelength = tempstr.length();
  //for (int cc=0; cc< cn; cc++)
  //	stat.outputBit+=Huffman_DC_Map[(int)diffDC[cc]];
  UINT8 runlength=0;
  //find EOB
  int eob=qNode.size().area()-2;
  for (int i=eob; i>=0; i--)
    {
      //find the last non zero coefficent position
      if (quanRst(zig_zag_x[i],zig_zag_y[i],0) != Vec<T,cn>::all(0))
        {
          eob = i+1;
          break;
        }
      else if (i==0)
        {
          eob = 1;
          break;
        }
    }
  for (int i=0; i<eob; i++)
    {
      if (quanRst(zig_zag_x[i],zig_zag_y[i],0)==Vec<T,cn>::all(0))
        {
          runlength++;
          if(runlength==16)
            {
              codelength+= 11; //ZRL 11111111001
              stat.outputBit+="11111111001";
              outputfile<<"11111111001";
              runlength =0;
              totalBits+=11;
              qNode.bitcount.AddCodeBit("11111111001");
            }
        }
      else
        {
          tempstr=ComputeJpegHuff(quanRst(zig_zag_x[i],zig_zag_y[i],0),runlength,Huffman_AC_Map);
          qNode.bitcount.AddCodeBit(tempstr);
          codelength+=tempstr.length();
          //	for (int cc=0; cc< cn; cc++)
          //		stat.outputBit+=Huffman_AC_Map[(int)quanRst(zig_zag_x[i],zig_zag_y[i],0)[cn]];
          runlength = 0;//lear runlength
        }
    }
  //add the eob bit in
  stat.outputBit+="1010";
  outputfile<<"1010";
  codelength+=4;
  totalBits+=4;
  qNode.bitcount.AddCodeBit("1010");
  stat.baseCoderBitLength += codelength;
  //reconstruct result
  //	floatVerNode.Print();
  stat.jBits[high_low] += codelength;
  stat.jBlockNum[high_low] ++;

  std::cout<<"BPP of JPEG "<<high_low<<" block: "<< qNode.offset()<< double(codelength)/qNode.size().area()<<endl;
  cv::idct(floatVerNode[0],floatVerNode[0],DFT_SCALE);
  rst.SetBlock(qNode.offset(),Tensor<T,cn>(floatVerNode));
  //floatVerNode.Print();
  SetCausalMap(qNode);
  //update localVarMap
  localVarMap.SetBlock(Point3i(qNode.offset().x/8,qNode.offset().y/8,qNode.offset().z),floatVerNode.LocalVariance(floatVerNode.LocalMean(Point3i(8,8,1),Point3i(8,8,1)),Point3i(8,8,1),Point3i(8,8,1)));
  //localVarMap(Cube(0,0,0,16,16,1)).Print();
#if CV_MINOR_VERSION <6
  rst_with_seam.SetBlock(qNode.offset(),Tensor<T,cn>(floatVerNode).CvtColor<3>(CV_GRAY2RGB));
#else
  rst_with_seam.SetBlock(qNode.offset(),Tensor<T,cn>(floatVerNode).CvtColor<3>(COLOR_GRAY2RGB));
#endif
  //stat.numBlocksInLevel[qLevel]++;
  //stat.footBitInLevel[qLevel]+=codelength;
  return qNode;
}

void MTC::Coding(void)
{
  Coding(this->mode);
}
void MTC::Coding(CodingMode codemode)
{
  this->mode = codemode;
  fstream fl;
  fl.open(".\\lighting.txt",ios::out);
  fl.close();
  ensemble.Display(3000,1);
  //2D only
  QNode<double,1> node;
  ///////init parameters;
  //UpdateParameters();//not necessary, gjin sep 11
  if (path.length()==0)
    {
      CV_Error(CV_StsError,"no path setted up");
    }
  else
    outputfile.open(path+prefix+"bits_"+baseName+".txt",ios::out);
  if (codemode == CodingMode::CODING_MODE_POST_TP)
    {
      rst = ensemble;
#if CV_MINOR_VERSION <6
      rst_with_seam = rst.CvtColor<3>(CV_GRAY2RGB);
#else
      rst_with_seam = rst.CvtColor<3>(COLOR_GRAY2RGB);
#endif

    }
  if (lightCorrectionType == LightingCorrectionType::PQI_LF_ENCODING)
    this->PQICodingLFComponent();
  //CollectLighing();
  //ensemble.Display();
  double orgQualityThrd = qualityThrd;
  for (int i= 0; i < gridSize.height; i++)
    {
      for (int j=0; j< gridSize.width; j++)
        {
          std::cout<<"("<<i<<","<<j<<")"<<endl;
          if (codemode == CodingMode::CODING_MODE_POST_TP)
            {
              QTree<T,cn> &tempNode = GetNode(Point3i(i,j,0));
              PostCoding(tempNode);//QNode<T,cn>(rst,tempNode.size(),tempNode.offset(),tempNode.overlap()));
            }
          else if (codemode == CodingMode::CODING_MODE_JPEG || codemode == CodingMode::CODING_MODE_MTC)
            {

              QTree<T,cn> qNode = GetNode(Point3i(i,j,0));
              qNode.SetBlock(ensemble.Crop(qNode.offset(),qNode.size()));
              //qNode.Print();
              switch (edgeMap[i][j])
                {
                case 0:
                  qualityThrd =orgQualityThrd;// qfactor * orgQualityThrd;
                  break;
                case 1:
                  qualityThrd = orgQualityThrd;//qfactor * orgQualityThrd;
                  break;
                case 2:
                  qualityThrd = orgQualityThrd;
                  break;
                case 3:
                  qualityThrd = orgQualityThrd;//qfactor*qfactor * orgQualityThrd;
                  break;
                case 4:
                  qualityThrd = orgQualityThrd;
                  break;
                case 5:
                  qualityThrd = orgQualityThrd;//qfactor * orgQualityThrd;
                  break;
                default:
                  qualityThrd = orgQualityThrd;
                  break;
                }
              /*if (i==6 && j==19)
                                {
                                        rst.SaveBlock("temp_rst.tif");
                                        causalMap.SaveBlock("temp_cuz.tif");
                                }*/
              CodingJPEG(qNode);
              std::cout<<float(stat.lightingCodeLength16)/float(stat.tpBlockNum16)/16/16<<endl;
              //rst.SaveBlock(".\\tempRst.tif");
              //	std::cout<<stat.lightingCodeLength16<<endl;
              //rst_with_seam.Display(1000,1);
            }
          else
            CodingMTC(GetNode(Point3i(i,j,0)));
        }
    }
  this->localVarMap.Print();
  stat.footbitLength = 0;
  stat.treebitLength = 0;
  stat.PSNR = metric::ComputePSNR(ensemble, rst);
  for (int k=0; k< stat.levels; k++)
    {
      stat.psnrInLevel[k] /= (double)stat.numBlocksInLevel[k];
      stat.psnrInLevel[k] = 10.0*log10(255.0*255.0/stat.psnrInLevel[k]);
      if (k<2)
        {
          stat.psnrTextureInLevel[k] /= (double)stat.numTextureInLevel[k];
          stat.psnrTextureInLevel[k] = 10.0*log10(255.0*255.0/stat.psnrTextureInLevel[k]);
        }
      stat.bppInLevel[k] = (double)stat.footBitInLevel[k]/(double)stat.numBlocksInLevel[k]/(double)stat.pixInLevel[k];
      stat.footbitLength += (stat.footBitInLevel[k] + stat.textureBitInLevel[k]);//*(int)(log((double)this->candidNum)/log(2.0)));
      if (k < 2)
        stat.treeBitInLevel[k] = stat.numTestInLevel[k] + stat.numTextureInLevel[k];
      else if (mode!= CodingMode::CODING_MODE_MTC)
        stat.treeBitInLevel[k] = stat.numTestInLevel[k];
      if (k== levels -3)
        stat.treeBitInLevel[k] = stat.numBlocksInLevel[k];
      if (k == levels - 2)
        stat.treeBitInLevel[k] = stat.numBlocksInLevel[k];//if no 1x1 block, each 2x2 block will cost 2 bit for rect split 10 and 11
      //if (k== levels-1)
      //	stat.treeBitInLevel[k] = stat.numBlocksInLevel[k]*2/4;

      stat.treebitLength+= stat.treeBitInLevel[k];
    }
  if (postBlendType==PostBlendingType::POST_BLENDING_OFFLINE) //offline post blending
    {
      rst.SaveBlock(path+prefix+baseName+"_no_post.tif");
      PostBlending(stat.tprecs);
    }
  //	if (lightCanVec.size()>0)
  //		SaveLighting();

  //add LF back
  rst = (rst)+LFImage;


  if (this->mode == CodingMode::CODING_MODE_JPEG)
    rst.SaveBlock(path+prefix+baseName+".jpg");
  else
    {
      rst.SaveBlock(path+prefix+baseName+".tif");
      rst_with_seam.SaveBlock(path+prefix+"seam_"+baseName+".tif");
      footmap.SaveBlock(path+prefix+"footmap_"+baseName+".tif");
      vector<string> filenames;
      filenames.push_back(cFileName);
      filenames.push_back(path+prefix+baseName+".tif");
      filenames.push_back(path+prefix+"seam_"+baseName+".tif");
      string tempfilename = path+"combin_"+prefix+baseName+".tif";
      tensor::CombineImage(filenames,tempfilename);
    }
  int lightingbits = 0;
  if (codemode==CodingMode::CODING_MODE_MTC&&lightCorrectionType==LightingCorrectionType::HAS_LIGHTING_CORRECTION)
    lightingbits= (lightCanVec.size()*lightCanVec[0].size()*6 + 64*6*32);
  //stat.bpp = double(lightingbits+stat.footbitLength+stat.treebitLength)/double(ensemble.size().height)/double(ensemble.size().width);
  stat.bpp = double(totalBits)/double(ensemble.size().area());
  /*if (rst.size().height != rst.size().width)
                std::cout<<"Do not compute SSIM because image is not square size\n";
        else
        {
                std::cout<<"Computing final SSIM..."<<endl;
                rst.SetSubWinSize(this->GetSubWinSize());
                stat.SSIM = rst.Compare(ensemble,COMPARE_CRITERIA_SSIM,3,4,false);
                std::cout<<"SSIM= "<<stat.SSIM<<endl;
        }*/
  fstream logfile;
  int wth = 10;
  logfile.open(path+prefix+"data_"+baseName+".txt",ios::out);
  logfile.precision(7);
  logfile.fill(' ');
  logfile<<"                                                  "<<"horizontal    vertical    pixels \n";
  logfile<<"                                    image size    "<<ensemble.size().width<<"           "<<ensemble.size().height<<"           "<<ensemble.size().area()<<endl<<endl;
  logfile<<"                                 block size N=    " <<"32            16             8             total"<<endl;
  logfile<<"                     pixels in block of size N    "<<32*32<<"           "<<16*16<<"           "<<8*8<<endl;
  logfile<<"                     blocks in image of size N    "<<(ensemble.size().width/32)*(ensemble.size().height/32)<<"          "
        <<(ensemble.size().width/16)*(ensemble.size().height/16)<<"          "
       <<(ensemble.size().width/8)*(ensemble.size().height/8)<<endl;
  logfile<<"                      blocks of size N encoded    "<<stat.tpBlockNum32<<"           "<<stat.tpBlockNum16<<"           "<<stat.numBlocksInLevel[2]<<"              "<<stat.tpBlockNum32+stat.tpBlockNum16+stat.numBlocksInLevel[2]<<endl;
  logfile<<"            pixels encoded in blocks of size N    "<<stat.tpBlockNum32*32*32<<"          "<<stat.tpBlockNum16*16*16<<"          "<<stat.numBlocksInLevel[2]*8*8<<"          "
        <<stat.tpBlockNum32*32*32+stat.tpBlockNum16*16*16+stat.numBlocksInLevel[2]*8*8<<endl;
  logfile<<"         % of image encoded in block of size N    "<<100*float(stat.tpBlockNum32*32*32)/float(ensemble.size().area())<<"    "
        <<100*float(stat.tpBlockNum16*16*16)/float(ensemble.size().area())<<"    "<<100*float(stat.numBlocksInLevel[2]*8*8)/float(ensemble.size().area())<<endl;
  logfile<<"          nontree coded bits for blocks size N    "<<stat.lightingCodeLength32<<"          "<<stat.lightingCodeLength16<<"          "<<stat.baseCoderBitLength<<"          "<<stat.lightingCodeLength16+stat.lightingCodeLength32+stat.baseCoderBitLength<<endl;
  //logfile<<"tree bits"<<"\t"<<stat.treebitLength<<"\t"<<"candid indicator bits"<<"\t"<<stat.tpBlockNum32*2 + stat.tpBlockNum16*2<<endl;

  logfile<<"                 tree bits for block of size N    "<<(ensemble.size().width/32)*(ensemble.size().height/32)<<"          "
        <<(ensemble.size().width/16)*(ensemble.size().height/16) - stat.tpBlockNum32*4<<"        0        "<<stat.treebitLength<<endl;
  logfile<<"            candidate indicator bits of size N    "<<stat.tpBlockNum32*(log(float(candidNum))/log(2.0))<<"            "<<stat.tpBlockNum16*(log(float(candidNum))/log(2.0))<<"       0       "<<(log(float(candidNum))/log(2.0))*(stat.tpBlockNum16+stat.tpBlockNum32)<<endl;
  logfile<<"              total # bits for block of size N    "<<stat.lightingCodeLength32+(ensemble.size().width/32)*(ensemble.size().height/32)+stat.tpBlockNum32*(log(float(candidNum))/log(2.0))<<"           "
        <<stat.lightingCodeLength16+(ensemble.size().width/16)*(ensemble.size().height/16) - stat.tpBlockNum32*4+stat.tpBlockNum16*(log(float(candidNum))/log(2.0))<<"         "
       <<stat.baseCoderBitLength<<"           "<<totalBits<<endl;
  logfile<<"                  non tree code rate of size N    "<<float(stat.lightingCodeLength32)/float(stat.tpBlockNum32*32*32)<<"    "<<float(stat.lightingCodeLength16)/float(stat.tpBlockNum16*16*16)<<"    "
        <<float(stat.baseCoderBitLength)/float(stat.numBlocksInLevel[2]*8*8)<<"    "<<float(stat.lightingCodeLength16+stat.lightingCodeLength32+stat.baseCoderBitLength)/float(ensemble.size().area())<<endl;
  logfile<<"                      tree code rate of size N    "<<float((ensemble.size().width/32)*(ensemble.size().height/32))/float(stat.tpBlockNum32*32*32)<<"    "
        <<float((ensemble.size().width/16)*(ensemble.size().height/16) - stat.tpBlockNum32*4)/float(stat.tpBlockNum16*16*16)<<"      0       "<<float(stat.treebitLength)/float(ensemble.size().area())<<endl;
  logfile<<"                candid ind. bit rate of size N    "<<float(stat.tpBlockNum32*(log(float(candidNum))/log(2.0)))/float(stat.tpBlockNum32*32*32)<<"    "
        <<float(stat.tpBlockNum16*(log(float(candidNum))/log(2.0)))/float(stat.tpBlockNum16*16*16)<<"      0       "<<float((log(float(candidNum))/log(2.0))*(stat.tpBlockNum16+stat.tpBlockNum32))/float(ensemble.size().area())<<endl;
  logfile<<"                     total code rate of size N    "<<float(stat.lightingCodeLength32+(ensemble.size().width/32)*(ensemble.size().height/32)+stat.tpBlockNum32*(log(float(candidNum))/log(2.0)))/float(stat.tpBlockNum32*32*32)<<"    "
        <<float(stat.lightingCodeLength16+(ensemble.size().width/16)*(ensemble.size().height/16) - stat.tpBlockNum32*4+stat.tpBlockNum16*(log(float(candidNum))/log(2.0)))/float(stat.tpBlockNum16*16*16)<<"    "
       <<float(stat.baseCoderBitLength)/float(stat.numBlocksInLevel[2]*8*8)<<"    "<<float(totalBits)/float(ensemble.size().area())<<endl;
  logfile<<"        bpp contrib. of nontree bits of size N    "<<float(stat.lightingCodeLength32)/float(ensemble.size().area())<<"    "<<float(stat.lightingCodeLength16)/float(ensemble.size().area())<<"    "
        <<float(stat.baseCoderBitLength)/float(ensemble.size().area())<<"    "<<float(stat.lightingCodeLength16+stat.lightingCodeLength32+stat.baseCoderBitLength)/float(ensemble.size().area())<<endl;
  logfile<<"      bpp contrib. of tree code rate of size N    "<<float((ensemble.size().width/32)*(ensemble.size().height/32))/float(ensemble.size().area())<<"    "
        <<float((ensemble.size().width/16)*(ensemble.size().height/16) - stat.tpBlockNum32*4)/float(ensemble.size().area())<<"      0       "<<float(stat.treebitLength)/float(ensemble.size().area())<<endl;
  logfile<<"bpp contrib. of candid ind. bit rate of size N    "<<float(stat.tpBlockNum32*(log(float(candidNum))/log(2.0)))/float(ensemble.size().area())<<"    "
        <<float(stat.tpBlockNum16*(log(float(candidNum))/log(2.0)))/float(ensemble.size().area())<<"      0      "<<float((log(float(candidNum))/log(2.0))*(stat.tpBlockNum16+stat.tpBlockNum32))/float(ensemble.size().area())<<endl;
  logfile<<"bpp contrib. of total bits for block of size N    "<<float(stat.lightingCodeLength32+(ensemble.size().width/32)*(ensemble.size().height/32)+stat.tpBlockNum32*(log(float(candidNum))/log(2.0)))/float(ensemble.size().area())<<"    "
        <<float(stat.lightingCodeLength16+(ensemble.size().width/16)*(ensemble.size().height/16) - stat.tpBlockNum32*4+stat.tpBlockNum16*(log(float(candidNum))/log(2.0)))/float(ensemble.size().area())<<"    "
       <<float(stat.baseCoderBitLength)/float(ensemble.size().area())<<"    "<<float(totalBits)/float(ensemble.size().area())<<endl;


  //
  logfile<<"PSNR: "<<stat.PSNR;
  logfile<<"\t BPP: " <<stat.bpp;
  ////logfile<<"\t SSIM: "<<stat.SSIM<<endl;
  //logfile<<"numBlocksInLevel:\n";
  //for (int i=0; i< this->levels; i++)
  //	logfile<<stat.numBlocksInLevel[i]<<"\t";
  //logfile<<endl;
  //logfile<<"footBitInLevel:\n";
  //for (int i=0; i< this->levels; i++)
  //	logfile<<stat.footBitInLevel[i]<<"\t";
  //logfile<<endl;
  //logfile<<"psnrInLevel:\n";
  //for (int i=0; i< this->levels; i++)
  //	logfile<<stat.psnrInLevel[i]<<"\t";
  //logfile<<endl;
  //logfile<<"textureBlockInLevel:\n";
  //for (int i=0; i< this->levels; i++)
  //	logfile<<stat.numTextureInLevel[i]<<"\t";
  //logfile<<endl;
  //logfile<<"texturePSNR:\n";
  //for (int i=0; i< this->levels; i++)
  //	logfile<<stat.psnrTextureInLevel[i]<<"\t";
  //logfile<<endl;
  //logfile<<"16x16 LF bpp"<<float(stat.lightingCodeLength16)/float(stat.tpBlockNum16)/16/16<<endl;
  //logfile<<"32x32 LF bpp"<<float(stat.lightingCodeLength32)/float(stat.tpBlockNum32)/32/32<<endl;
  ////output total number of blocks
  //logfile<<"32x32\t 16x16\t 8x8\t Total\n";
  //int tempTotalBlocks = stat.tpBlockNum16+ stat.tpBlockNum32+stat.numBlocksInLevel[2];
  //logfile<<stat.tpBlockNum32<<"\t"<<stat.tpBlockNum16<<"\t"<<stat.numBlocksInLevel[2]<<"\t"<<tempTotalBlocks<<endl;
  ////out put bit usage
  //logfile<<"32 lighting\t 16 lighting\t 8 JPEG \t TreeBits \t Total\n";
  //logfile<<stat.lightingCodeLength32<<"\t"<<stat.lightingCodeLength16<<"\t"<<stat.baseCoderBitLength<<"\t"<<stat.treebitLength<<"\t"<<this->totalBits<<endl;
  //logfile<<"bit percentage"<<endl;
  //logfile<<float(stat.lightingCodeLength32)/float(totalBits)<<"\t"<<float(stat.lightingCodeLength16)/float(totalBits)<<"\t"<<float(stat.baseCoderBitLength)/float(totalBits)<<"\t"<<float(stat.treebitLength)/float(totalBits)<<"\t"<<this->totalBits<<endl;
  ////out put bpp %
  //logfile<<"bpp contribution\n";
  //logfile<<float(stat.lightingCodeLength32)/float(stat.tpBlockNum32)/32/32<<"\t"<<float(stat.lightingCodeLength16)/float(stat.tpBlockNum16)/16/16<<"\t"<<float(stat.baseCoderBitLength)/float(stat.numBlocksInLevel[2])/8/8<<float(stat.treebitLength)/float(ensemble.size().area())<<endl;
  //logfile.close();
  ////output string version bits
  //outputfile.close();
  //std::cout<<"total Bits count "<<totalBits<<endl;
  //std::cout<<"bits in outputBits "<<stat.outputBit.length()<<endl;
  //logfile.open(path+prefix+"bits_cb_"+baseName+".txt",ios::out);
  //logfile<<stat.outputBit<<endl;
  logfile.close();
  rst.Display(3000,1);

}

bool MTC::AdaptiveTPSS(const QTree<T,cn>& qTree, int qLevel)
{
  bool tempAccept=true;
  vector<pair<Point3i,FootItem>> feet;
  QNode<T,cn> tar;
  if (this->mode != CodingMode::CODING_MODE_TPSS)
    {
      tempAccept = false;
    }
#ifdef ADAPTIVE_TPSS
  else if (qLevel>1) //at most level 1 -- 3 feet + up and left feet
    {
      tempAccept = false;
      //reject
      stat.outputBit+="1";
      outputfile<<"1";
      totalBits++;
      stat.treeBitInLevel[qLevel]++;
      stat.treebitLength++;
    }
#endif
  else
    {
      tar = this->GetValidNode(qTree.offset()-(qTree.overlap()).Point3(),qTree,ensemble);
      if (qLevel>=0)
        {
          RetrieveFeet(tar,qLevel);//retrieve foot from target possition
        }

      // TODO asfsafd
      //in retrieveFoot, make a list of feet required
      //and then in PossionLighintCorrection, compute each un-determined feet
      //interpolate by the list feet
      //update bits!!!!!!
      Tensor<double,1> orglt = lighting::ComputeTPSS(tar.GetExtendTensor(),0.001);
      tar.PoissonLightingCorrection(tar.Clone(),ensembleExt,ensembleExt,qTree.overlap(), footComputeRegion, footComputeMethod,this->initQSize);
      UpdateFeet(tar.getFeetCopy());
      Tensor<double,1> temp = lighting::ComputeTPSS(tar.GetExtendTensor(), 0.001);
      double diff = metric::ComputeMSE(orglt.Crop((tar.size()/2).Point3(),tar.size()), temp.Crop((tar.size()/2).Point3(),tar.size()));
#   ifdef ADAPTIVE_TPSS
      if (diff/double(tar.size().area())<qualityThrd)
        {
#   endif
          rst.SetBlock(qTree.offset(),temp(Cube((tar.size()/2).Point3(),tar.size())));
          stat.lightingDiff+=diff;
#   ifdef ADAPTIVE_TPSS
          if (qLevel == 0) //accept use 1 feet
            {
              stat.outputBit+="00";
              outputfile<<"00";
              totalBits+=2;
              stat.treeBitInLevel[qLevel]+=2;
              stat.treebitLength+=2;

            }
          else if (qLevel>0) //accept use more feet7
            {
              stat.outputBit+="01";
              outputfile<<"01";
              totalBits+=2;
              stat.treeBitInLevel[qLevel]+=2;
              stat.treebitLength+=2;
            }
#   endif
          for (auto& f : tar.getFeetCopy())
            {
              if (f.second.bits>0)
                {
                  stat.outputBit+=f.second.bitstring;
                  outputfile<<f.second.bitstring;
                  totalBits+=f.second.bits;
                  stat.plcBitLength+=f.second.bits;
                }
            }

#   ifdef ADAPTIVE_TPSS
        }
      else //increase number of feet
        {
          if (qTree.size().height==16 && qLevel==1) //get the minimum size, always accept but leave it blank
            {
              rst.SetBlock(qTree.offset(),Tensor<double,1>(qTree.size()));
              stat.lightingDiff+=diff;
              tempAccept = true;
              stat.outputBit+="1"; //reject but still use 3 feet
              outputfile<<"1";
              totalBits++;
              stat.treeBitInLevel[qLevel]++;
              stat.treebitLength++;
              for (auto& f : tar.getFeetCopy())
                {
                  if (f.second.bits>0)
                    {
                      stat.outputBit+=f.second.bitstring;
                      outputfile<<f.second.bitstring;
                      totalBits+=f.second.bits;
                      stat.plcBitLength+=f.second.bits;
                    }
                }
            }
          else
            {
              RemoveFeet(tar.getFeetCopy());
              tempAccept=AdaptiveTPSS(qTree,qLevel+1);
            }


        }
#   endif
    }
  return tempAccept;
}
bool MTC::TexturePrediction(QTree<T,cn>& qNode, int qLevel)
{
  bool tempAccept = true;
  //check boundary condition
  if (this->mode == CodingMode::CODING_MODE_TPSS)
    {
      tempAccept = false;
    }
  else if (qNode.size().height <= 8 && qNode.size().width <= 8)//force to set min blk size 32
    {
      tempAccept =false;
    }
  else if(qNode.offset().x == 0 && qNode.offset().y==0 && qNode.offset().z==0)
    {
      tempAccept = false;
    }
  else if(qNode.offset().x - qNode.overlap().height < 0 || qNode.offset().y - qNode.overlap().width < 0 || qNode.offset().z - qNode.overlap().depth < 0)
    {
      tempAccept = false;
    }
  else if(qNode.offset().x + qNode.size().height >= ensemble.size().height ||
          qNode.offset().y + qNode.size().width  >= ensemble.size().width )
    {
      tempAccept = false;
    }
  if (mode == CodingMode::CODING_MODE_PQI || mode== CodingMode::CODING_MODE_JPEG)//pure pqi, pure jpeg, turn off TP
    tempAccept = false;
  if (tempAccept == false)
    {
      stat.numTestInLevel[qLevel]++; //case 1: no match performed
    }
  else
    {
      //qNode.Print();
      //qNode.Print();
      //qNode.leftBound.Print();
      //qNode.upBound.Print();
      //		(rst+128).Display();
      if (qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y && qNode.size().height==DEBUG_SIZE)
        {
          qNode.GetExtendTensor().Print("tar_to_be_match",true);
          //temp pause here
          /*int temp;
                        cin>>temp;*/
          this->causalMap.SaveBlock("cmap.tif");
        }
      /* RetrieveFoot(qNode,fLevel);*/
      //cout<<"begin matching....\n";
      vector<cv::Point3i> matchCandid(candidNum,Point3i(-1,-1,-1));
      matchCandid = BoundaryMatching(qNode,matching_method,mseThrd);

      //qNode.bounds[0].Display();
      //qNode.bounds[1].Display();
      //if (matchCandid[candidNum-1] == Point3i(-1,-1,-1) ) // no match find
      if (qNode.offset().x==48&& qNode.offset().y==96)
        cout<<qNode.size().height<<endl;
      if (matchCandid.size()==0)
        {
          stat.numTestInLevel[qLevel]++; //case 2: matching, but no match found
          tempAccept = false;
        }
      else
        {
          //match light corrected version now!
          //(rst+128).Display();
          int index = -1;
          fLevel = -1;
          index = IsAcceptPredict(matchCandid,qNode,metricModifier,fLevel);
          //cout<<fLevel<<endl;
          /*
      if (metricModifier == M_DIST)
        index = IsAcceptPredict(matchCandid,qNode,COMPARE_CRITERIA_MAHALANOBIS,3,4);
      else if (metricModifier== LRI_METRIC)
        index = IsAcceptPredict(matchCandid,qNode,COMPARE_CRITERIA_LRI,0);
      else
        index = IsAcceptPredict(matchCandid,qNode,COMPARE_CRITERIA_SSIM,3,4);
        */

          if (index >= 0)
            {
              //record the tar and cand position for debugging

              //accept bit, indicate by "0" (no split) and store candid index
              stat.outputBit+="0";
              totalBits++;
              outputfile<<"0";
              string bset_to_string_str;
              boost::dynamic_bitset<> bset((int)ceil(log(double(/*candidNum*/matchCandid.size()))/log(2.0)),index);//gj 01162013
              boost::to_string(bset,bset_to_string_str);
              stat.outputBit+=bset_to_string_str;
              totalBits+= bset_to_string_str.length();
              outputfile<<bset_to_string_str;
              stat.numTextureInLevel[qLevel]++;
              stat.textureBitInLevel[qLevel]+= (int)ceil(log((double)this->candidNum)/log(2.0));
              qNode.bitcount.SetDecBit("0");
              qNode.bitcount.AddCodeBit(bset_to_string_str);
              qNode.bitcount.SetCodeMethod(CodingMethodNames::CODING_MTC);
              //update tpblock

              stat.tpBlockNum[qNode.size().height>>5]++;
              stat.tpBits[qNode.size().height>>5]+= (bset_to_string_str.length()+1);
              //if exceed the image boundary
              if ( matchCandid[index].x+2*qNode.overlap().height+qNode.size().height < ensemble.size().height &&
                   matchCandid[index].y+2*qNode.overlap().width+qNode.size().width < ensemble.size().width &&
                   qNode.GetFootPos().x+qNode.overlap().height < ensemble.size().height &&
                   qNode.GetFootPos().y+qNode.overlap().width  < ensemble.size().width )
                {
                  stat.tprecs.push_back(TPrecord(matchCandid[index]+qNode.overlap().Point3(),Cube(qNode.offset(),qNode.size()),qNode.overlap()));
                }
              //the candidate must be a copy from the best matching place
              QNode<T,cn> temp(rst,qNode.size(),matchCandid[index] + qNode.overlap().Point3(),qNode.overlap());
              QNode<T,cn> candid = temp.Clone();
#if CV_MINOR_VERSION <6
              Tensor<T,3> temp3= candid.CvtColor<3>(CV_GRAY2RGB);
#else
              Tensor<T,3> temp3= candid.CvtColor<3>(COLOR_GRAY2RGB);
#endif
	      Vec<T,3> red(0,255,0);
	      for (int ii=0; ii< temp3.size().height; ii++)
		for (int jj=0; jj< temp3.size().width; jj++)
		  {
		    if (ii==0||jj==0||ii==temp.size().height-1||jj==temp.size().width-1)
		      temp3(ii,jj,0)= red;
		  }
	      pred_rst.SetBlock(qNode.offset(),temp);

	      //ensemble(Cube(qNode.offset(),qNode.size())).Display();
	      //candid.Display();
	      //if (qNode.offset().x==32 && qNode.offset().y==336)
	      //	rst.SaveBlock("tempSave.tif");
	      //if (lightCorrectionType!=NO_LIGHTING_CORRECTION)
	      //{
	      //	Tensor<T,cn> fromLight = LFImage(Cube(matchCandid[index],candid.size()+(qNode.overlap()*2)));
	      //		Tensor<T,cn> toLight = LFImage(Cube(qNode.offset()-qNode.overlap().Point3(),candid.size()+(qNode.overlap()*2)));
	      //		//qNode.GetBoundaryUp(qNode.
	      //		candid.LightingCorrection2(fromLight,toLight);
	      //	}
	      //if (qNode.offset().x==32 && qNode.offset().y==336)
	      //{
	      //	candid.Print();
	      //}
	      //candid.SetBlock(Tensor<T,cn>(candid.size()));
	      //rst.Display();
	      /*	if (qNode.offset().x== 4*32 && qNode.offset().y==1*32)
				{
					rst.Display();
					rst_with_seam.Display();
					qNode.leftBound.Display(0);
					qNode.upBound.Display(0);
					std::cout<<matchCandid[index]<<endl;

				}*/

	      //if ( qNode.GetBoundaryUp().IsInside(target.GetFootPos()) ||
	      //	qNode.GetBoundaryLeft().IsInside(target.GetFootPos())||
	      //	qNode.GetBoundaryUp().IsInside(target.offset()+Point3i(target.size().height,0,0)) ||
	      //	qNode.GetBoundaryLeft().IsInside(target.offset()+Point3i(0,target.size().height,0)))
	      //{
	      //	std::cout<<"overlaped blocks offset"<<target.offset()<<endl;
	      //	//target = QNode<T,cn>(target); //if the target overlap with qNode's boundary, make a copy to avoid data loss during quilting
	      //}
	      //
	      //correct the lighting of the target, include the boundary(up, down, left, right)
	      QNode<T,cn> tar2 = QNode<T,cn>(rst,qNode.size(),qNode.offset(),qNode.overlap()).Clone();
	      QNode<T,cn> cand2 = GetValidCandidSimple(matchCandid[index],qNode);
	      cout<<"matchCandid"<<matchCandid[index]<<endl;
	      //cand2.Display();
	      QNode<T,cn> cand2Copy = cand2.Clone();
	      QNode<T,cn> org2 = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(),qNode.overlap());
	      if (qNode.size().height==DEBUG_SIZE&&qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
		{
		  tar2.GetExtendTensor(1,1,1,1).Print("tarExt",true);
		  cand2.GetExtendTensor(1,1,1,1).Print("best_candExt",true);
		  cout<<"=============="<<index<<endl;
		}
	      if(lightCorrectionType==LightingCorrectionType::HAS_LIGHTING_CORRECTION && acceptDirect == false)
		{
		  Tensor<T,cn> forLight = ensemble.Crop(qNode.offset()-qNode.overlap().Point3(),qNode.size()+(qNode.overlap()*2)); // I don't do TP on the boundary of image, so it is safe to do crop this way

		  //	rst.Display();
		  candid.LightingCorrection(forLight); //correct the lighing, include 4 boundaries to the candidate

		  //	rst.Display();
		  //approximate (guess) bits
		  totalBits+=candid.lt.GetLightingCodeLength(); //each dct coeff use 8 on average, and 10 significant coeffs.
		  if (candid.size().height>=32)
		    {
		      stat.lightingCodeLength32 += candid.lt.GetLightingCodeLength();
		      stat.tpBlockNum32++;
		    }
		  else
		    {
		      stat.lightingCodeLength16 += candid.lt.GetLightingCodeLength();
		      stat.tpBlockNum16++;
		    }
		  lightTagVec.push_back(candid.lt.GetTagLighting());
		  lightCanVec.push_back(candid.lt.GetCanLighting());
		}
	      else if (lightCorrectionType == LightingCorrectionType::POISSON_LC && acceptDirect==false)
		{
		  //Tensor<T,cn> temp = rst.Crop(candid.offset()-qNode.overlap().Point3(),candid.size()+(candid.overlap()*2)); // I don't do TP on the boundary of image, so it is safe to do crop this way
		  //temp.Print();
		  //= QNode<T,cn>(rst,qNode.size(),matchCandid[index]+qNode.overlap().Point3(),qNode.overlap()*2).Clone();

		  //Tensor<T,cn> temp = ensemble.Crop(qNode.offset()+(qNode.size()/2).Point3(),qNode.size()); // I don't do TP on the boundary of image, so it is safe to do crop this way
		  ////Tensor<T,cn> temp = ensemble.Crop(qNode.offset()-(qNode.overlap().Point3()),qNode.size()+qNode.overlap());
		  ////compute triangle mean
		  ///*double tempmean=0;
		  //for (int ii=0;ii<temp.size().height; ii++)
		  //	for (int jj=temp.size().height-ii;jj<temp.size().width;jj++)
		  //		tempmean+=temp(ii,jj,0)[0];
		  //tempmean/=(temp.size().area()/2);*/
		  //tar2.SetFoot(temp.Mean());
		  //temp = ensemble.Crop(Point3i(qNode.offset().x-qNode.size().height/2, qNode.offset().y+qNode.size().width/2,0),qNode.size());
		  //tar2.upNodeFoot = temp.Mean();
		  //temp = ensemble.Crop(Point3i(qNode.offset().x+qNode.size().height/2, qNode.offset().y-qNode.size().width/2,0),qNode.size());
		  //tar2.leftNodeFoot =temp.Mean();

		  //maybe store the uncorrected extended version of boundary here before do PLC
		  //this is the border + 1 pixel surrounded
		  auto footEntry = RetrieveFoot(org2,fLevel); //the block used as the key to RetrieveFoots, otherwise the initial value of foots will be wrong

		  // cand2.GetExtendTensor().Display();
		  //! PLC only apply to block+LU boundary
		  if (fLevel>=0)//if use PLC
		  {
		    cand2.PoissonLightingCorrection(/*org2 does not available at decoder*/org2,ensemble,rst,qNode.overlap(),footComputeRegion, footComputeMethod,this->initQSize);
		      // cand2.GetExtendTensor().Display();
		    if (footEntry==FootTable.end())
			  {
			  UpdateFoots(tar2,cand2);

			  }
		  }
		  //also need to update bits here.

		  if (qNode.size().height==DEBUG_SIZE&&qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
		    {
		      cand2.GetExtendTensor(1,1,1,1).Print("best_candExt_PLC",true);
		    }
		  candid = QNode<T,cn>(cand2.GetExtendTensor(1,1,1,1),qNode.size(),cand2.overlap().Point3(),qNode.overlap());
		  pred_afterLC.SetBlock(qNode.offset(),candid);
		  /*vector<Tensor<double,1>> dummy;
	  dummy.push_back(org2);
	  dummy.push_back(cand2);
	  cand2.DisplayAll(dummy,1,2);*/
		}
	      if (qNode.size().height==DEBUG_SIZE&&qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
		{
		  candid.GetExtendTensor(1,1,1,1).Print("best_cand_PLC",true);
		}
	      //UpdatePBSet(qNode,matchCandid[index]);//modify it to store the lighting correcjted candidate
	      if (lightCorrectionType == LightingCorrectionType::POISSON_LC)
		UpdatePBSet(qNode,cand2Copy);//cand2Copy contains uncorrected version of cand2
	      else
		UpdatePBSet(qNode,cand2);
	      //qNode.SaveBlock(path+prefix+"_tag("+boost::lexical_cast<string>(qNode.offset().x)+","+boost::lexical_cast<string>(qNode.offset().y)+".tif");
	      //candid.SaveBlock(path+prefix+"_can("+boost::lexical_cast<string>(qNode.offset().x)+","+boost::lexical_cast<string>(qNode.offset().y)+".tif");

	      //qNode.SetBlock(candid);
	      //qNode.GetBoundaryUp().SetBlock(candid.GetBoundaryUp());
	      //qNode.GetBoundaryLeft().SetBlock(candid.GetBoundaryLeft());
	      //==== tempoary disable Aug 14//
	      /* candid.upBound.Print();
	candid.leftBound.Print();*/
	      if (blend_method == BlendingMethod::SHORTEST_PATH_BLENDING)
		qNode.Quilting(candid); // do forward quilting
	      else if (blend_method == BlendingMethod::GRADIENT_BLENDING)
		{

		  //if (qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
		  //{
		  //	rst.Display();
		  //}
		  //cand2.GetExtendTensor().Print();
		  //tar2.GetExtendTensor().Print();
		  tar2.GradientStitching(cand2,BlendingLocation::FORWARD_BLENDING,qNode.overlap(),Tensor<double,1>());//remember the order is tar2 change to cand2

		  //if (qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
		  //{
		  //	rst.Display();
		  //}
		  //tar2.Print("tar2");
		  //copy to qNode next....
		  //qNode.Print("qNode before");
		  //cand2.GetExtendTensor().Print();
		  //tar2.GetExtendTensor().Print();

		  candid = QNode<T,cn>(cand2.GetExtendTensor(1,1,1,1),qNode.size(),cand2.overlap().Point3(),qNode.overlap());
		  QNode<T,cn> target = QNode<T,cn>(tar2.GetExtendTensor(),qNode.size(),tar2.overlap().Point3(),qNode.overlap());
		  //if (qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
		  //	{
		  //	qNode.SetBlock(Tensor<T,cn>(candid.size(),255));
		  //rst.Display();
		  //qNode.SetBlock(candid);
		  //rst.Display();
		  //	qNode.leftBound.SetBlock(Tensor<T,cn>(target.leftBound.size(),255));
		  //	qNode.upBound.SetBlock(Tensor<T,cn>(target.upBound.size(),255));
		  //	rst.Display();
		  //	}
		  qNode.leftBound.SetBlock(target.leftBound);
		  qNode.upBound.SetBlock(target.upBound);
		  //	if (qNode.offset().x==DEBUG_X&& qNode.offset().y==DEBUG_Y)
		  //	{
		  //	rst.Display();
		  //	}
		  qNode.SetBlock(cand2);
		  //	if (qNode.offset().x==DEBUG_X&& qNode.offset().y==DEBUG_Y)
		  //	{
		  //	rst.Display();
		  //	}
		  //rst.Display();
		  //qNode.Print("qNode after");
		  //////
		}
	      else
		{
		  qNode.SetBlock(candid);
		}
	      if (qNode.size().height==DEBUG_SIZE&&qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
		{
		  qNode.GetExtendTensor().SaveBlock("cand_blend.tif");
		}
	      //rst.SetBlock(qNode.offset(),candid);
	      //this->UpdateRstSqr(Cube(qNode.offset()-qNode.overlap().Point3(),qNode.overlap()+qNode.size()));
	      //this->SAT = rstSqr.ComputeSAT(this->SAT,qNode.offset()-qNode.overlap().Point3(),qNode.GetFootPos());
	      UpdateSAT(qNode.offset()-qNode.overlap().Point3(),qNode.size()+qNode.overlap());
	      /*	if (qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
				{
					rst(Cube(qNode.offset()-qNode.overlap().Point3(),qNode.overlap()+qNode.size())).Print();
					rstSqr(Cube(qNode.offset()-qNode.overlap().Point3(),qNode.overlap()+qNode.size())).Print();
					SAT(Cube(qNode.offset()-qNode.overlap().Point3(),qNode.overlap()+qNode.size()+Size3(1,1,0))).Print();
				}*/
	      SetCausalMap(qNode);
#if CV_MINOR_VERSION <6
	      rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(CV_GRAY2RGB));
#else
	      rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(COLOR_GRAY2RGB));
#endif
              // rst.SetBlock(qNode.offset(),qNode);
              //rst.Display();
              //rst_with_seam.Display();
              //deal with seam
              //qNode.GetBoundaryUp().Print();
              //qNode.GetBoundaryLeft().Print();
              ////// do post blending
              //=== begin temp disable ======== Aug 14
              if (blend_method == BlendingMethod::SHORTEST_PATH_BLENDING)
                {
#if CV_MINOR_VERSION <6
                  Tensor<T,3> tempBound = qNode.GetBoundaryUp().CvtColor<3>(CV_GRAY2RGB);
#else
                  Tensor<T,3> tempBound = qNode.GetBoundaryUp().CvtColor<3>(COLOR_GRAY2RGB);
#endif
		  Vec<T,3> tempElm(0,0,255);
		  for (unsigned int b=0; b < qNode.seam[0].size(); b++)
		    {
		      tempBound[Point3i(qNode.seam[0][b],b,0)] = tempElm;
		    }
		  rst_with_seam.SetBlock(qNode.offset() - qNode.overlap().Point3(),tempBound);
#if CV_MINOR_VERSION <6
		  tempBound  = qNode.GetBoundaryLeft().CvtColor<3>(CV_GRAY2RGB);
#else
		  tempBound  = qNode.GetBoundaryLeft().CvtColor<3>(COLOR_GRAY2RGB);
#endif

		  for (unsigned int b=0; b < qNode.seam[1].size(); b++)
		    {
		      tempBound[Point3i(b,qNode.seam[1][b],0)] = tempElm;
		    }
		  rst_with_seam.SetBlock(qNode.offset() - qNode.overlap().Point3(),tempBound);
		}
	      //	rst.Display();
	      // rst_with_seam.Display(0);
	      //qNode.Print();
	      stat.psnrTextureInLevel[qLevel] += metric::ComputeMSE(qNode, ensemble.Crop(qNode.offset(),qNode.size()));
	      //	rst_with_seam.Display(1000,1);
	      //rst.Display();
	      //===== end of temp disable

	      std::cout<<"MTC block: "<< qNode.offset()<< qNode.size().height<<"x"<<qNode.size().width<<endl;
	      if (qNode.size().height==DEBUG_SIZE&&qNode.offset().x==DEBUG_X && qNode.offset().y==DEBUG_Y)
		{
		  rst.SaveBlock("rst_after.tif");

		}
#ifdef RECORD_EVERYTHING
              string extra_path = "./everything/";
              char idxstr[20];
              sprintf(idxstr,"rst_%04d",acount);
              acount++;
              rst.SaveBlock(extra_path+idxstr+".bmp");
#endif
	    }
	  else
	    {
	      tempAccept = false;
	      stat.numTestInLevel[qLevel]++;
	    }
	}
    }
  return tempAccept;
}

int MTC::IsAcceptPredict(const vector<Point3i>& matchCandid, QTree<T,cn>& qNode,MetricModifier metricModifier, int level)//, double param1, double param2)
{
  double distance = INT_MIN, temp;
  CompareCriteria criteria = ParseCriteria(metricModifier);
  if (metricModifier == MetricModifier::MAHALANOBIS_DIST)
    distance = INT_MAX;
  int index = -1;
#ifndef METRIC_PARALLEL
  QNode<T,cn> candid;
  QNode<T,cn> org;
  Tensor<T,cn> tpss;
#endif
  Tensor<T,cn> org_mse;
  Tensor<T,cn> candid_mse;
  bool debugsignal=false;
  int bits_for_PLC_foot=0;
  string bits_stream_for_PLC_foot="";
  fstream scorefile;
  ofstream debugfile;

#ifdef RECORD_EVERYTHING
  ofstream everything;
  //double lightPSNR=0;
  double lightPSNRThrd = 10*log10(255*255/30); //db
  everything.open("./everything/everything.txt",ios::app);
#endif
  auto footEntry = FootTable.end();
  int N = (qNode.size().height*qNode.overlap().width+qNode.size().width*qNode.overlap().height);

  if (qNode.offset().x == DEBUG_X && qNode.offset().y==DEBUG_Y&&qNode.size().height==DEBUG_SIZE)
    {
      debugsignal = true;
      rst.SaveBlock("rst.tif");
      scorefile.open("scores.txt",ios::out);
    }


  if (criteria == CompareCriteria::COMPARE_CRITERIA_INTERP)
    {
      //unsupport
      exit(-2);
    }
  else
    {
      //get extra border for PLC
#ifndef METRIC_PARALLEL
      //if (lightCorrectionType == POISSON_LC)
      //{
      //    //if (criteria == COMPARE_CRITERIA_SSIM || criteria == COMPARE_CRITERIA_MAHALANOBIS||criteria==COMPARE_CRITERIA_LRI)
      //org = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(), qNode.size()/2);//no need for extra boundary
      //    org = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(), qNode.overlap());//no need for extra boundary
      //	//org = ensemble.Crop(qNode.offset() - (qNode.size()/2).Point3(),qNode.size()*Size3(2,2,1));
      //}
      //else
      //org = ensemble.Crop(qNode.offset(),qNode.size());
      org = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(),qNode.overlap());
      //org.Display();
      if (lightCorrectionType == LightingCorrectionType::POISSON_LC)
        {
          tpss = lighting::ComputeTPSS(org,0.01);
        }
      //org.Display();
      //org = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(),qNode.overlap());
#ifdef RECORD_EVERYTHING
      string extra_path = "./everything/";
      //_mkdir(extra_path.c_str());
      char posstr[20];
      sprintf(posstr,"_%d_(%d_%d)",qNode.size().height,qNode.offset().x,qNode.offset().y);
      //org.Print();
      org.GetExtendTensor(1,1,0,0)./*Crop(Point3i(org.size()/4),org.size()+org.size()/4).*/SaveBlock(extra_path+"org"+string(posstr)+".png");
      everything<<"org size "<<qNode.size().height<<" with LU:("<<org.offset().x<<","<<org.offset().y<<") matches:"<<endl;
      // org.GetExtendTensor(1,1,0,0).Crop(Point3i(org.size()/4),org.size()+org.size()/4).Display();
#endif
      if(debugsignal)
        {
#ifdef WIN32
          string PID = boost::lexical_cast<string>(_getpid());
#else
          string PID = boost::lexical_cast<string>(getpid());
#endif

          debugfile.open("./temp/ssim_terms.txt",ios::out);
          debugfile<<"org mean"<<org.Mean()[0]<<"\t"<<"org var"<<org.Var()[0]<<endl;
          debugfile.close();
        }

#else
#     ifdef RECORD_EVERYTHING
      string extra_path = "./everything/";
      mkdir(extra_path.c_str());
      char posstr[20];
      sprintf(posstr,"_%d_(%d_%d)",qNode.size().height,qNode.offset().x,qNode.offset().y);
      QNode<double,1> orgtemp;
      if (lightCorrectionType == POISSON_LC)
        {
          //if (criteria == COMPARE_CRITERIA_SSIM || criteria == COMPARE_CRITERIA_MAHALANOBIS||criteria==COMPARE_CRITERIA_LRI)
          orgtemp = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(), qNode.size()/2);
          //org = ensemble.Crop(qNode.offset() - (qNode.size()/2).Point3(),qNode.size()*Size3(2,2,1));
          orgtemp.GetExtendTensor(1,1,0,0).Crop(Point3i(qNode.size()/4),qNode.size()+qNode.size()/4).SaveBlock(extra_path+"org"+string(posstr)+".png");
        }
      else
        {
          //org = ensemble.Crop(qNode.offset(),qNode.size());
          orgtemp = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(),qNode.overlap());
          if (lightCorrectionType == PREDEF_LIGHTING)
            (orgtemp.GetExtendTensor(1,1,0,0).Crop(Point3i(0,0,0),qNode.size()+qNode.size()/4)+128).SaveBlock(extra_path+"org"+string(posstr)+".png");
          else
            orgtemp.GetExtendTensor(1,1,0,0).Crop(Point3i(0,0,0),qNode.size()+qNode.size()/4).SaveBlock(extra_path+"org"+string(posstr)+".png");
        }
      everything<<"org size "<<qNode.size().height<<" with LU:("<<orgtemp.offset().x<<","<<orgtemp.offset().y<<") matches:"<<endl;
#     endif
#endif
      vector<double> vdist(matchCandid.size());
      vector<double> ldist(matchCandid.size());
      vector<thread> threads;
      for (int i=0; i< (int)matchCandid.size(); i++)
        {

#ifdef METRIC_PARALLEL
	  threads.push_back(thread([&](const int i, vector<double>::iterator tempdist, vector<double>::iterator lightdist){
	      QNode<T,cn> candid;
	      QNode<T,cn> org;
	      Tensor<T,cn> tpss;
	      //get extra border for PLC
	      //if (lightCorrectionType == POISSON_LC)
	      //{
	      //  org = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(), qNode.size()/2);
	      //}
	      //else
	      org = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(),qNode.overlap());
	      tpss = org.ComputeTPSS(0.01);
	      //! 20130913 add tpss of un-corrected candid

#endif


	      //if (criteria == COMPARE_CRITERIA_SSIM ||criteria == COMPARE_CRITERIA_MAHALANOBIS||criteria==COMPARE_CRITERIA_LRI)
	      //!20130913                 if (lightCorrectionType == LightingCorrectionType::POISSON_LC)
	      //!20130913{
	      //candid= QNode<T,cn>(rst,qNode.size(),matchCandid[i]+qNode.overlap().Point3(),qNode.overlap());
	      if (metricModifier==MetricModifier::STSIM2_ADT_NAIVE)
		{
		  //verify??????
		  acceptDirect = false;
		  //candid.Print();
		  //add pre-compare process
		  org_mse= ensemble.Crop(qNode.offset(),qNode.size());
		  //org_mse.Print();
		  candid_mse = rst.Crop(matchCandid[i]+qNode.overlap().Point3(),qNode.size());
		  //candid_mse.Print();
		  /*	if (lightCorrectionType == HAS_LIGHTING_CORRECTION)
						candid_mse.LightingCorrection(org_mse,false);*/
		  //if (lightCorrectionType!=NO_LIGHTING_CORRECTION)
		  //{
		  //		Tensor<T,cn> fromLight = LFImage(Cube(candid_mse.offset(),candid_mse.size()));
		  //		Tensor<T,cn> toLight = LFImage(Cube(org_mse.offset(),org_mse.size()));
		  //		candid_mse= candid_mse - fromLight + toLight;//.LightingCorrection2(fromLight,toLight);
		  //	}
		  //qNode.GetBoundaryUp(qNode.

		  double mse = metric::Compare(org_mse, candid_mse,CompareCriteria::COMPARE_CRITERIA_MSE);
		  //if (mse < 15)
		  //{
		  //	printf("mse = %f < 15, accept directly\n",mse);
		  //	if (debugsignal==true)
		  //	{
		  //		char str[20];
		  //		sprintf(str,"cand_%d",i);
		  //		candid = this->GetValidCandid(matchCandid[i],qNode);
		  //		candid.GetExtendTensor().Print(string(str));
		  //	//org.GetExtendTensor().Print("org_before_PLC");
		  //	}
		  //	acceptDirect = true;
		  //	return i;
		  //}
		  //else
		  //{
		  double var = org.LocalVariance(org.LocalMean(org.size()),org.size())(0,0,0)[0];
		  if (var > varThrd1)
		    {
		      printf("var = %f > %f, increase stsim2 thrd\n",mse,varThrd1);
		      qualityThrd = orgQualityThrd*qfactor;
		    }
		  else
		    {
		      qualityThrd = orgQualityThrd;
		    }

		  //}
		}
	      //if (lightCorrectionType==POISSON_LC)
	      //{
	      //	;//do nothing
	      //}
	      //else
	      //{

	      //20130605 candid = this->GetValidCandid(matchCandid[i],qNode);
	      candid = this->GetValidCandidSimple(matchCandid[i],qNode);
	      //if (matchCandid[i].x + qNode.overlap().height - qNode.size().height/2 < 0 ||
	      //	matchCandid[i].y + qNode.overlap().width - qNode.size().width/2 <0 ||
	      //	matchCandid[i].x + qNode.overlap().height + 3 * qNode.size().height/2 > ensemble.size().height ||
	      //	matchCandid[i].y + qNode.overlap().width + 3*qNode.size().width/2 > ensemble.size().width)
	      //{
	      //	Point3i pos1,pos2;
	      //	pos1 = Point3i(matchCandid[i].x + qNode.overlap().height - qNode.size().height/2, matchCandid[i].y + qNode.overlap().width - qNode.size().width/2 , 0);
	      //	Size3 sz(qNode.size()*Size3(2,2,1));
	      //	Tensor<T,cn> temp1(sz,0);
	      //	if (matchCandid[i].x+qNode.overlap().height - qNode.size().height/2 < 0 )
	      //	{
	      //		pos1.x=0;
	      //		pos2.x= qNode.size().height/2 -( matchCandid[i].x+qNode.overlap().height );
	      //		sz.height -= pos2.x;
	      //	}
	      //	if (	matchCandid[i].y + qNode.overlap().width - qNode.size().width/2 <0 )
	      //	{
	      //		pos1.y=0;
	      //		pos2.y= qNode.size().width/2 - (matchCandid[i].y + qNode.overlap().width);
	      //		sz.width -= pos2.y;
	      //	}
	      //	if (matchCandid[i].x + qNode.overlap().height + 3 * qNode.size().height/2 > ensemble.size().height)
	      //	{
	      //		pos1.x=matchCandid[i].x + qNode.overlap().height ;
	      //		pos2.x=0;
	      //		sz.height-= qNode.size().height/2 - (matchCandid[i].x + qNode.overlap().height + 3 * qNode.size().height/2 - rst.size().height);
	      //	}
	      //	if (matchCandid[i].y + qNode.overlap().width + 3*qNode.size().width/2 > ensemble.size().width)
	      //	{
	      //		pos1.y=matchCandid[i].y + qNode.overlap().width;
	      //		pos2.y=0;
	      //		sz.width-= qNode.size().width/2 - (matchCandid[i].y + qNode.overlap().width + 3*qNode.size().width/2 -  rst.size().width);
	      //	}
	      //	Tensor<T,cn> temp2 = rst.Crop( pos1, sz);
	      //	//temp2.Print();
	      //	temp1.SetBlock(pos2,temp2);
	      //	candid = QNode<T,cn>(temp1,qNode.size(),(qNode.size()/2).Point3(),qNode.size()/2).Clone();
	      ////old version, use same boundary
	      ////
	      ////	candid = rst.Crop( qNode.offset() - (qNode.size()/2).Point3(),qNode.size()*Size3(2,2,1));
	      ////	Point3i pos = (qNode.size() * Size3(1,3,1) / 2).Point3();
	      ////	Size3 sz = qNode.size() * Size3(3,1,1) / Size3(2,2,1);
	      ////	candid.SetBlock(pos,ensemble.Crop( qNode.offset() - (qNode.size()/2).Point3() + pos, sz));
	      //////candid.Print();
	      ////	sz = qNode.size()* Size3(1,3,1)/ Size3(2,2,1);
	      ////	pos = (qNode.size() * Size3(3,0,1) / 2).Point3();
	      ////	candid.SetBlock(pos,ensemble.Crop( qNode.offset() - (qNode.size()/2).Point3() + pos, sz));
	      //////candid.Print();
	      ////	pos = (qNode.size()/2).Point3();
	      ////	candid.SetBlock(pos,rst.Crop( matchCandid[i] + qNode.overlap().Point3(), qNode.size()));
	      //}
	      //else
	      //{
	      //	//candid.Print();
	      ////new version, use the true boundary of node
	      ////matchCandid get the position of boundary matching, which + overlap size is the candid position
	      ////candid = rst.Crop(matchCandid[i] + qNode.overlap().Point3() - (qNode.size()/2).Point3(),
	      ////				  qNode.size()*Size3(2,2,1));
	      //candid = QNode<T,cn>(rst,qNode.size(),matchCandid[i]+qNode.overlap().Point3(),
	      //	qNode.size()/2).Clone();
	      //}
	      //}

	      //! 20130913 }
	      //! 20130913 else
	      //candid = rst.Crop(matchCandid[i] + qNode.overlap().Point3(), qNode.size());
	      //! 20130913 candid = QNode<T,cn>(rst,qNode.size(),matchCandid[i] + qNode.overlap().Point3(),qNode.overlap());
	      //if (qNode.size().height==16) //16x16
	      //{
	      //	org.SetSubWinSize(this->stsimSubWinSize/2);
	      //	org.SetSubWinStep(this->stsimSubWinStep/2);
	      //}
	      //else
	      //{
	      /*org.SetSubWinSize(this->stsimSubWinSize);
			org.SetSubWinStep(this->stsimSubWinStep);
			candid.SetSubWinSize(this->stsimSubWinSize);
			candid.SetSubWinStep(this->stsimSubWinStep);*/
	      //}
	      Tensor<T,cn> orgExt = org.GetExtendTensor(1,1,1,1);
	      QNode<T,cn> candidExt = candid.GetExtendTensor(1,1,1,1);
	      candidExt.ExtendBoundary(qNode.overlap(),0);//20130605
	      //orgExt.Print();
	      if (debugsignal)
		candidExt.Print();

	      //orgExt.SetSubWinSize(this->stsimSubWinSize);
	      //orgExt.SetSubWinStep(this->stsimSubWinStep);
	      //candidExt.SetSubWinSize(this->stsimSubWinSize);
	      //candidExt.SetSubWinStep(this->stsimSubWinStep);
	      if (debugsignal==true)
		{
		  char str[20];
		  sprintf(str,"cand_%d",i);
		  candid.GetExtendTensor().Print(string(str),true);
		  org.GetExtendTensor().Print("org_before_PLC",true);
		  //rst.Display();
		  //orgExt.debugtrigger=true;
		}
#ifdef RECORD_EVERYTHING
              //save candid before PLC
              char idxstr[20];
              sprintf(idxstr,"_(%d_%d)_(%d)",matchCandid[i].x+qNode.overlap().height,matchCandid[i].y+qNode.overlap().width,i);
              if (lightCorrectionType == LightingCorrectionType::PREDEF_LIGHTING)
                (candid.GetExtendTensor(1,1,0,0)+128)/*20130605.Crop(Point3i(org.size()/4),org.size()+org.size()/4)*/.SaveBlock(extra_path+"cand"+string(posstr)+string(idxstr)+".png");
              else
                candid.GetExtendTensor(1,1,0,0)/*20130605.Crop(Point3i(org.size()/4),org.size()+org.size()/4)*/.SaveBlock(extra_path+"cand"+string(posstr)+string(idxstr)+".png");
#endif

	  if (lightCorrectionType == LightingCorrectionType::HAS_LIGHTING_CORRECTION)
		  candidExt.LightingCorrection(orgExt,false);//do lighting correction before stsim2
	  else if (lightCorrectionType == LightingCorrectionType::POISSON_LC)
		{
		  //need to do foot recording here
		  //yes, update foot values, here, use a dynamic table to record each computed foot
		  // (x,y) value, bit
	    //for (int ll=-1; ll<2; ll++)
		  //{
		    //QNode<T,cn> tempCand = candid.Clone(); //!20130916 can clone foot information?
		    footEntry = RetrieveFoot(org,level);//this and compute foot, and footEntry  can be done outside loop, only once
		    if (level>=0) //no PLC if level <0
		    {
		      candid.PoissonLightingCorrection(org/*change to tar??,no*/,ensemble, rst,qNode.overlap(), footComputeRegion, footComputeMethod,this->initQSize);
		    }
		    // candid.GetExtendTensor(1,1,1,1).Display();
		    Tensor<T,cn> candTPSS = lighting::ComputeTPSS(candid, 0.01); //after correction
		    //vector<Tensor<T,cn>> templight;
		    /* templight.push_back(org);
	      templight.push_back(candid);
	      templight.push_back(tpss);
	      templight.push_back(candTPSS);
	      candTPSS.DisplayAll(templight,2,2);*/
	#ifdef METRIC_PARALLEL
	*lightdist = metric::ComputePSNR(tpss, candTPSS);
	#else
	ldist[i] = metric::ComputePSNR(tpss, candTPSS);
	#endif
	/*double lightPSNRBeforePLC = tpss.ComputePSNR(candTPSSBefore);
	if (lightPSNRBeforePLC>lightPSNR)
	lightPSNR = lightPSNRBeforePLC;*/
	#ifdef RECORD_EVERYTHING
	//save candidate after PLC
	candid.GetExtendTensor(1,1,0,0)/*20130605.Crop(Point3i(org.size()/4),org.size()+org.size()/4)*/.SaveBlock(extra_path+"cand_plc"+string(posstr)+string(idxstr)+"_("+to_string(level+1)+").png");
	//candid.GetExtendTensor(1,1,1,1).SaveBlock(extra_path+"cand_full_plc"+string(posstr)+string(idxstr)+".png");
	#endif
	if (footEntry==FootTable.end()&&level>=0)
	  UpdateFoots(org, candid);

        if (debugsignal==true)
        {
          char str[20];
          sprintf(str,"cand_PLC_%d",i);
          candid.GetExtendTensor().Print(string(str),true);//should change
          //57org.GetExtendTensor().Print("org_after_PLC",true);//shoube be the same as before_PLC
        }
        //candid.Display();
        //org.Display();
        if (criteria == CompareCriteria::COMPARE_CRITERIA_LRI)
        {
          candidExt = candid.GetExtendTensor(1,1,0,0).GetBlock(Cube(qNode.overlap().Point3(),qNode.size()+qNode.overlap()));
          orgExt = org.GetExtendTensor(1,1,0,0).GetBlock(Cube(qNode.overlap().Point3(),qNode.size()+qNode.overlap()));
        }
        else
        {
          Size3 bsize = candid.size();
          // candidExt = candid.GetExtendTensor(0,1,1,1,1);
          //candidExt.SetBlock((bsize/4).Point3(), candid.GetExtendTensor(1,1,1,1).GetBlock(Cube((bsize/4).Point3(),bsize+bsize/2)));// zero padding
          //20130606 only take U,L border+blk
          candidExt = candid.GetExtendTensor(1,1,0,0);
          orgExt = org.GetExtendTensor(1,1,0,0);//.Crop(qNode.overlap().Point3(),qNode.size()+qNode.overlap());
          //orgExt = org.GetExtendTensor(0,1,1,1,1).Crop((bsize/4).Point3(),bsize+qNode.overlap());
          // orgExt.SetBlock((bsize/4).Point3(),org.GetExtendTensor(1,1,1,1).GetBlock(Cube((bsize/4).Point3(),bsize+bsize/2)));
        }
        //	candidExt.Print();
        //temp = orgExt.Compare(candidExt,criteria,param1,param2,true,stsim2PoolType);              }
        //}//for levels
    }//if PLC

              //else
              /*
      Tensor<T,cn> orgPadding(orgExt.size());
      orgPadding.SetBlock((orgExt.size()/4).Point3(),orgExt(Cube((orgExt.size()/4).Point3(),qNode.size())));
      orgPadding.SetSubWinSize(orgExt.GetSubWinSize());
      orgPadding.SetSubWinStep(orgExt.GetSubWinStep());
      Tensor<T,cn> candPadding(candidExt.size());
      candPadding.SetBlock((candidExt.size()/4).Point3(),candidExt(Cube((candidExt.size()/4).Point3(),qNode.size())));
      if(debugsignal==true)
      {
        orgPadding.Display();
        candPadding.Display();
      }
        */
          if (debugsignal)
          {
            #ifdef WIN32
            string PID = boost::lexical_cast<string>(_getpid());
            #else
            string PID = boost::lexical_cast<string>(getpid());
            #endif
            debugfile.open("./temp/ssim_terms.txt",ios::app);
            debugfile<<" -------------"<<i<<"th candidate-----------------"<<endl;
            debugfile<<"candid mean\t"<<candid.Mean()[0]<<"\t"<<"candid var\t"<<candid.Var()[0]<<endl;
            debugfile.close();
            candidExt.debugtrigger=true;
          }
          //  std::thread tt(&orgExt.Compare,candidExt,criteria,param1,param2,true,stsim2PoolType,metricModifier, this->iMahaCovar);
          //zero padding
          if (metricModifier == MetricModifier::SE_MSE||metricModifier ==MetricModifier::STSIM2_SE_MSE)
            {
              Tensor<T,cn> candidBd = candidExt(Cube(qNode.overlap().Point3(),qNode.size()+qNode.overlap()*2));
              Tensor<T,cn> orgBd = orgExt(Cube(qNode.overlap().Point3(),qNode.size()+qNode.overlap()*2));
              temp = metric::Compare(orgBd, candidBd,criteria,this->subSize, this->subStep, (int)metricModifier, qNode.size().height,qNode.overlap().height);
              //compute PSNR
              temp = 10*log10(double(N)*65025.0/temp);
              //temp = double(qNode.size().height*qNode.overlap().width+qNode.size().width*qNode.overlap().height)/temp;
            }
          else if (criteria == CompareCriteria::COMPARE_CRITERIA_SSIM)
            {
              //20130516 since PLC only take care of L,U boundary + blk, SSIM cannot compare block with boundary
              //temp = orgExt.Compare(candidExt,criteria,3,4,FILTER_BOUND_VALID/*true*/,stsim2PoolType,metricModifier);
              //temp = orgExt.Compare(candidExt,criteria,3,4,FILTER_BOUND_HALF/*true*/,stsim2PoolType,metricModifier);
              //20130516 try take original
              //Tensor<T,cn> orgPd=org.ExtendBoundary(Size3(org.size().height/2, org.size().width/2,0));
              //Tensor<T,cn> canPd = candid.ExtendBoundary(Size3(org.size().height/2,org.size().width/2,0));
              //Tensor<T,cn> orgPd=org.GetExtendTensor(1,1,0,0).Crop(Point3i(org.size()/4),org.size()+org.size()/4).ExtendBoundary((qNode.size()+qNode.overlap())/2);
              //Tensor<T,cn> canPd=candid.GetExtendTensor(1,1,0,0).Crop(Point3i(org.size()/4),org.size()+org.size()/4).ExtendBoundary((qNode.size()+qNode.overlap())/2);
              //candidExt.Display();
              if (metricModifier==MetricModifier::STSIM3_LSE)
                {
                  //Metric mc;
                  Tensor<T,cn> orgPLC = orgExt.Crop((orgExt.size()/5).Point3(),orgExt.size()/5*4+Size3(0,0,1));
                  Tensor<T,cn> candPLC = candidExt.Crop((candidExt.size()/5).Point3(),candidExt.size()/5*4+Size3(0,0,1));
                  temp = metric::Compare(orgPLC, candPLC,criteria,this->subSize, this->subStep, 3,4,(int)FilterBoundary::FILTER_BOUND_FULL,(int)stsim2PoolType,(int)metricModifier,0,debugsignal);
                }
              else
                {
                  //! from 20130909 I decide to use core block (does not include boundary) to compute STSIM
                  //! no padding is used too.
                  // cout<<"orgExt size="<<orgExt.size()<<endl;
                  Tensor<T,cn> orgPd = orgExt.Crop((orgExt.size()/5).Point3(),orgExt.size()/5*4+Size3(0,0,1));
                  Tensor<T,cn> canPd = candidExt.Crop((orgExt.size()/5).Point3(),orgExt.size()/5*4+Size3(0,0,1));
                  //orgPd.Print("orgPd");
                  //canPd.Print("canPd");
                  temp = metric::Compare(orgPd,canPd,criteria,this->subSize, this->subStep,3,4,(int)FilterBoundary::FILTER_BOUND_FULL/*true*/,(int)stsim2PoolType,(int)metricModifier,0,debugsignal);
                  //cout<<"20130912 i="<<i<<", temp"<<temp<<endl;
                }
              //temp = org.Compare(candid,criteria,3,4,FILTER_BOUND_FULL/*true*/,stsim2PoolType,metricModifier);
              if (temp>1)//do it again when something goes wrong
                {
                  //candid.Print("cand",true);
                  //org.Print("org",true);
                  //candid.debugtrigger=true;
                  temp = metric::Compare(org, candid,criteria,this->subSize, this->subStep,3,4,(int)FilterBoundary::FILTER_BOUND_FULL/*true*/,(int)stsim2PoolType,(int)metricModifier,0,debugsignal);
                }
#ifdef METRIC_PARALLEL
              *tempdist = temp;
#else
              vdist[i]=temp;
#endif
            }
          else if (criteria == CompareCriteria::COMPARE_CRITERIA_SVM)
            {
              if (orgExt.size().height<=32)//only do 32x32
                temp = 0;
              else
                {
                  temp= metric::Compare(orgExt.GetBlock(Cube(8,8,0,48,48,1)), candidExt.GetBlock(Cube(8,8,0,48,48,1)),criteria);
                }
            }
          else
            temp = metric::Compare(orgExt, candidExt,criteria);
          //temp = orgPadding.Compare(candPadding,criteria,param1,param2,true,stsim2PoolType,metricModifier, this->iMahaCovar);
          if (debugsignal)
            {
              std::cout<<"cand"<<i<<"_PLC score:"<<temp<<endl;
              scorefile<<"cand"<<i<<"_PLC score:"<<temp<<endl;

            }
          if (((temp > distance)&&(criteria != CompareCriteria::COMPARE_CRITERIA_MAHALANOBIS)) ||((temp<distance)&&(criteria == CompareCriteria::COMPARE_CRITERIA_MAHALANOBIS)))
            {
              distance = temp;
              index = i;
              bits_stream_for_PLC_foot = candid.GetBits();
              bits_for_PLC_foot = bits_stream_for_PLC_foot.length();
            }
          //print all candidates
          vdist[i]=temp;
#ifdef METRIC_PARALLEL
            },i,vdist.begin()+i,ldist.begin()+i));
#endif
        }
#ifdef METRIC_PARALLEL
      for (auto& t:threads)
        t.join();
#endif



      //if ( distance >= this->qualityThrd)
      //adaptive quality
      //if (criteria == COMPARE_CRITERIA_SVM && qNode.size().height==32)
      //{
      //  vector<Tensor<double,1>> dummy;
      //    rst.Display();
      //  org_mse= ensemble.Crop(qNode.offset()-qNode.overlap().Point3(),qNode.size()+qNode.overlap()*2);
      //  dummy.push_back(org_mse);
      //  dummy.push_back(org_mse);
      //  for (int ii=0; ii< matchCandid.size(); ii++)
      //  {
      //    candid_mse = rst.Crop(matchCandid[ii],qNode.size()+qNode.overlap()*2);
      //    dummy.push_back(candid_mse.Clone());
      //    cout<<vdist[ii]<<endl;
      //    //candid_mse.Display();
      //  }
      //  qNode.DisplayAll(dummy,5,2);
      //}
#   ifdef RECORD_EVERYTHING
      for (unsigned int j=0; j<matchCandid.size();j++)
        {
          everything<<"candid "<<j<<" with LU: ("<<matchCandid[j].x+qNode.overlap().height<<","<<matchCandid[j].y+qNode.overlap().width<<"), "<<"flevel: "<<level<<" light: "<<ldist[j]<<" , score: "<<vdist[j]<<endl;
        }
      everything.close();
#   endif
      int accepted = -1;


#   ifdef METRIC_PARALLEL
      distance = vdist[0];
      index = 0;
      for (int idx = 1; idx< candidNum; idx++)
        {
          if (distance < vdist[idx])
            {
              distance = vdist[idx];
              index = idx;
            }
        }
#   endif

      if (criteria == CompareCriteria::COMPARE_CRITERIA_MAHALANOBIS)
        {

          if(qNode.size().height >16 && distance <= this->qualityThrd) //good
            accepted = index; //return index;
          else if (distance <= this->qualityThrd * this->qfactor)//increase the quality in small block
            accepted = index;//return index;
          else
            {
              if (lightCorrectionType == LightingCorrectionType::POISSON_LC) //reset the FootTable
                {
                  if (footEntry!=FootTable.end())
                    FootTable.erase(footEntry);
                }
              accepted = -1;//return -1;
            }
        }
      else if (criteria == CompareCriteria::COMPARE_CRITERIA_SSIM)
        {
          double thresholdAdaptor = 1;
          double lightAdaptor = 1;
          //if (qNode.upBound.Var()[0]<50 && qNode.leftBound.Var()[0]<50)
          //  thresholdAdaptor = 1.13;//gj20130121 smooth region increase thrd
          //#ifdef METRIC_PARALLEL
          //      QNode<T,cn> org = QNode<T,cn>(ensemble,qNode.size(),qNode.offset(),qNode.overlap());
          //#endif
          //thresholdAdaptor = (log10((1.0+org.Var()[0])/16834)+6)/4;
          // double up = 1;
          //double low = 0.5;
          ////double k = 0.8;
          ////double a = (1-k)/(up-low);
          ////double b = 1.0-a*up;
          //thresholdAdaptor>up?thresholdAdaptor=up:thresholdAdaptor=thresholdAdaptor;
          //thresholdAdaptor<low?thresholdAdaptor=low:thresholdAdaptor=thresholdAdaptor;
          //thresholdAdaptor = a*thresholdAdaptor+b;
          //  thresholdAdaptor = 1 - thresholdAdaptor*1;
          ////06272013 tpss light result affect threshold, twice the threshold if the light is not good
          //if (lightCorrectionType == POISSON_LC)
          //{
          //
          //  if (ldist[index] < lightPSNRThrd)
          //    lightAdaptor = 2;
          //}
          //if (metricModifier==STSIM3_LSE)
          //{
          //  if (abs(distance) < this->qualityThrd )
          //    accepted = index;
          //}

          //! 20130916 compute var
          //cout<<this->qfactor<<endl;
          double orgvar = org.Var()[0];
          if (orgvar > 10 && orgvar <200)//incease threaold
            thresholdAdaptor*=1.02;
          if (level<0 && orgvar < 100) //herustic set use PLC for smooth region
            accepted = -1;
          else if(qNode.size().height ==blockSize.height && distance >= this->qualityThrd*thresholdAdaptor*lightAdaptor) //good
            accepted = index; //return index;
          else if (distance >= this->qualityThrd * this->qfactor*thresholdAdaptor*lightAdaptor)//increase the quality in small block
            accepted = index;//return index;
          else //if not accepted, reset foottable
            {
              if (lightCorrectionType == LightingCorrectionType::POISSON_LC) //reset the FootTable
                {
                  if (footEntry!=FootTable.end())
                    FootTable.erase(footEntry);
                }
              accepted = -1;//return -1;
            }
        }
      else if (criteria==CompareCriteria::COMPARE_CRITERIA_MSE)
        {
          if (metricModifier == MetricModifier::SE_MSE)
            {
              //use PSNR
              if (distance>qualityThrd)
                accepted=index;
            }
          else if(metricModifier==MetricModifier::STSIM2_SE_MSE)
            {
              vector<Point3i> tempvp;
              tempvp.push_back(matchCandid[index]);
              if(IsAcceptPredict(tempvp,qNode,MetricModifier::STSIM2_NEW_L2)>=0)
                accepted=index;
            }
          else
            {
              accepted= index;
            }

        }
      else if (criteria ==CompareCriteria::COMPARE_CRITERIA_LRI)
        {
          if (distance>qualityThrd)
            accepted=index;
        }
      else if (criteria ==CompareCriteria::COMPARE_CRITERIA_SVM)
        {
          //if (distance>qualityThrd)
          if (distance>qualityThrd)
            accepted = index;
        }

      if (debugsignal)
        scorefile.close();
      if (accepted>=0)
        {
          if (qNode.offset().x == 664 && qNode.offset().y==24)
            {
              cout<<stat.outputBit<<endl;
              cout<<footEntry->first<<endl;
              cout<<footEntry->second.value<<endl;
              cout<<footEntry->second.bitstring;
            }

          if (lightCorrectionType==LightingCorrectionType::POISSON_LC)
            {
              if (footEntry==FootTable.end()) // make sure read the updated entry when only 1 candidate available
                //footEntry = RetrieveFoots(org);
                footEntry = RetrieveFoot(qNode,level);
              if (level>=0){
              stat.outputBit+=footEntry->second.bitstring;//bits_stream_for_PLC_foot;
              stat.plcBitLength+=footEntry->second.bitstring.length();
              this->totalBits += footEntry->second.bitstring.length();
              this->outputfile<<footEntry->second.bitstring;
              qNode.bitcount.SetFootBit(footEntry->second.bitstring);
              }
            }

          candPosLog<<"Tar: ("<<qNode.offset().x<<","<<qNode.offset().y<<") ==== Best Cand_"<<accepted<<": ("<<matchCandid[accepted].x+qNode.overlap().height<<","<<matchCandid[accepted].y+qNode.overlap().width<<") ====== Score: "<<distance<<endl;
          //if (criteria == COMPARE_CRITERIA_SSIM)
          //{
          for (int i=0; i<(int)matchCandid.size();i++)
            {
              if (i!=accepted)
                candPosLog<<"      >>>>>>>> Cand_"<<i<<": ("<<matchCandid[i].x+qNode.overlap().height<<","<<matchCandid[i].y+qNode.overlap().width<<") ====== Score: "<<vdist[i]<<endl;
            }
          candPosLog<<"---------------------------------"<<endl;
          //}

        }
      else //increase foot number and try again
        {
          if (level<1)//((qNode.size().height>>level)>8)
            {
              fLevel= level+1;
              accepted = IsAcceptPredict(matchCandid, qNode,metricModifier,fLevel);
            }
        }
      return accepted;
    }
}
void MTC::CodingJPEG(QTree<T,cn>& qTree, int qLevel)
{
  //qNode.Print();
  if ( !TexturePrediction(qTree,qLevel)) // if texture prediction fail
    {
      if (qTree.size().height==8 && qTree.size().width==8)
        {
          //ensemble.Display();
          QNode<T,cn> org = ensemble.Crop(qTree.offset(),qTree.size());
          /////////// if using Watson quan table , convert to YCbCr first
          Tensor<float,cn> floatVerNode(org.size());
          ////floatVerNode.Print();
          //if (jpegQuanTblType == JPEG_WASTON)
          //{
          //	Tensor<float,3> colorVerNode = org.CvtColor<3>(CV_GRAY2RGB);
          //	colorVerNode = colorVerNode.CvtColor<3>(CV_RGB2YCrCb);
          //	//colorVerNode.Print();
          //	vector<Mat> planes;
          //	cv::split(colorVerNode[0],planes);
          //	floatVerNode[0] = planes[0]-128;
          //}
          //else
          //	floatVerNode = Tensor<float,cn>(org-128);
          ////floatVerNode.Print();
          //cv::dct(floatVerNode[0],floatVerNode[0]);
          ////cv::idct(floatVerNode[0],floatVerNode[0],DFT_SCALE);
          ////floatVerNode.Print();
          qTree.bitcount.SetCodeMethod(CodingMethodNames::CODING_JPEG);
          floatVerNode = dctEnsemble.Crop(qTree.offset(),floatVerNode.size());
          Tensor<T,cn> quanRst(qTree.size()); //index of quan
          //floatVerNode.Print();
          //dctEnsemble.SetBlock(qNode.offset(),Tensor<double,cn>(floatVerNode));
          for (int x=0; x < qTree.size().height;x++)
            for (int y=0; y < qTree.size().width;y++)
              {
                quanRst(x,y,0) = qTree.UniformQuantize(floatVerNode(x,y,0),JpegQuanTbl[x*qTree.size().width+y]);
                floatVerNode(x,y,0) = quanRst(x,y,0).mul(Vec<T,cn>::all(JpegQuanTbl[x*qTree.size().width+y]));
                //quanRst(x,y,0)/JpegQuanTbl[x*qNode.size().width+y];
              }
          //quanRst.Print();
          //floatVerNode.Print();
          ////coding DC component
          ///treat different in different position
          Vec<T,cn> prevBlockDC;
          if (qTree.offset()==Point3i(0,0,0))
            {
              prevBlockDC = 0;
            }
          else if(qTree.offset().y==0)
            {
              prevBlockDC = qTree.UniformQuantize(dctEnsemble[qTree.offset()-Point3i(qTree.size().height,0,0)],JpegQuanTbl[0]);
            }
          else
            {
              prevBlockDC = qTree.UniformQuantize(dctEnsemble[qTree.offset()-Point3i(0,qTree.size().width,0)],JpegQuanTbl[0]);
            }
          Vec<T,cn> diffDC = quanRst(0,0,0) - prevBlockDC;
          string tempstr = ComputeJpegHuff(diffDC,Huffman_DC_Map);
          int codelength = tempstr.length();
          qTree.bitcount.AddCodeBit(tempstr);
          //for (int cc=0; cc< cn; cc++)
          //	stat.outputBit+=Huffman_DC_Map[(int)diffDC[cc]];
          UINT8 runlength=0;
          //find EOB
          int eob=qTree.size().area()-2;
          for (int i=eob; i>=0; i--)
            {
              //find the last non zero coefficent position
              if (quanRst(zig_zag_x[i],zig_zag_y[i],0) != Vec<T,cn>::all(0))
                {
                  eob = i+1;
                  break;
                }
              else if (i==0)
                {
                  eob = 1;
                  break;
                }
            }
          for (int i=0; i<eob; i++)
            {
              if (quanRst(zig_zag_x[i],zig_zag_y[i],0)==Vec<T,cn>::all(0))
                {
                  runlength++;
                  if(runlength==16)
                    {
                      codelength+= 11; //ZRL 11111111001
                      stat.outputBit+="11111111001";
                      outputfile<<"11111111001";
                      runlength =0;
                      totalBits+=11;
                      qTree.bitcount.AddCodeBit("11111111001");
                    }
                }
              else
                {
                  tempstr= ComputeJpegHuff(quanRst(zig_zag_x[i],zig_zag_y[i],0),runlength,Huffman_AC_Map);
                  codelength +=tempstr.length();
                  qTree.bitcount.AddCodeBit(tempstr);
                  //	for (int cc=0; cc< cn; cc++)
                  //		stat.outputBit+=Huffman_AC_Map[(int)quanRst(zig_zag_x[i],zig_zag_y[i],0)[cn]];
                  runlength = 0;//lear runlength
                }
            }
          //add the eob bit in
          stat.outputBit+="1010";
          outputfile<<"1010";
          codelength+=4;
          totalBits+=4;
          stat.baseCoderBitLength += codelength;
          qTree.bitcount.AddCodeBit("1010");
          //reconstruct result
          //	floatVerNode.Print();
          cv::idct(floatVerNode[0],floatVerNode[0],DFT_SCALE);
          floatVerNode = floatVerNode+128;
          //floatVerNode.Print();
          //floatVerNode.Print();
          rst.SetBlock(qTree.offset(),Tensor<T,cn>(floatVerNode));
          SetCausalMap(qTree);
          //update localVarMap

          localVarMap.SetBlock(Point3i(qTree.offset().x/searchStep.height,qTree.offset().y/searchStep.width,qTree.offset().z/searchStep.depth),floatVerNode.LocalVariance(floatVerNode.LocalMean(searchStep,searchStep),searchStep,searchStep));
          //localVarMap(Cube(0,0,0,16,16,1)).Print();
#if CV_MINOR_VERSION<6
          rst_with_seam.SetBlock(qTree.offset(),Tensor<T,cn>(floatVerNode).CvtColor<3>(CV_GRAY2RGB));
#else
          rst_with_seam.SetBlock(qTree.offset(),Tensor<T,cn>(floatVerNode).CvtColor<3>(COLOR_GRAY2RGB));
#endif
	  stat.numBlocksInLevel[qLevel]++;
	  stat.footBitInLevel[qLevel]+=codelength;
	}
      else
        {
          qTree.Split();
          if(this->mode == CodingMode::CODING_MODE_MTC)//split output 1bit to indicate
            {
              stat.outputBit+="1";
              outputfile<<"1";
              totalBits++;
              stat.treebitLength++;
              qTree.bitcount.SetDecBit("1");
            }
          QTree<T,cn> *tempTree;
          tempTree = qTree.NextLeaf();
          while(true)
            {
              if(tempTree == NULL)
                break;
              else
                {
                  CodingJPEG(*tempTree,qLevel+1);//QNode<T,cn>(rst,tempTree->size(),tempTree->offset(),qTree.overlap()/2),qLevel+1);
                }
              if(tempTree->GetTreePos()==tempTree->GetPeerTreeSize()-1)
                {
                  if (postBlendType == PostBlendingType::POST_BLENDING_ONLINE)
                    {
                      PBSet::iterator it = pbSet.find(PBRecord(qTree.offset(),BoundDir::LEFT));
                      if (it!=pbSet.end()&&it->second.height <= qTree.size().height)
                        {
                          PostBlending(it->first,it->second);
                          pbSet.erase(it->first	,it->second);
                          //rst_with_seam.Display(1000,1);
                        }
                      it = pbSet.find(PBRecord(qTree.offset(),BoundDir::UP));
                      if (it!=pbSet.end()&&it->second.width <= qTree.size().width)
                        {
                          PostBlending(it->first,it->second);
                          pbSet.erase(it->first,it->second);
                          //rst_with_seam.Display(1000,1);
                        }
                    }
                  break;
                }
              else
                tempTree=tempTree->NextLeaf();
            }
        }
    }
}

void MTC::Encode()
{
  Coding();
}

void MTC::Decode(ofstream& fh)
{

}
void MTC::CodingMTC(QTree<MTC::T,MTC::cn> &qNode, int qLevel)
{



  //std::cout<<"working on offset:"<<qNode.offset().x<<","<<qNode.offset().y<<endl;
  qNode.SetFoot(ensemble);
  QTree<double,1> *tempTree;
  bool tempAccept = true;
  QNode<double,1> interpNode;


  //replaced by method: TexturePrediction()
  //vector<cv::Point3i> matchCandid(candidNum,Point3i(-1,-1,-1));

  //qNode should be initilized in QGrid, only update boundary here!
  //UpdateBoundary(qNode);

  /////////////////////////////////////////
  //stat.numTestInLevel[qLevel]++;


  ////qNode.SetBlock(CropFromOrg(qNode.offset,qNode.GetSize()));
  //if (qNode.size().height < 16 && qNode.size().width < 16)
  //{
  //	tempAccept =false;
  //}
  //else if(qNode.offset().x == 0 && qNode.offset().y==0 && qNode.offset().z==0)
  //{
  //	tempAccept = false;
  //}
  //else if(qNode.offset().x - qNode.overlap().height < 0 || qNode.offset().y - qNode.overlap().width < 0 || qNode.offset().z - qNode.overlap().depth < 0)
  //{
  //	tempAccept = false;
  //}
  //else if(qNode.offset().x + qNode.size().height >= ensemble.size().height ||
  //	    qNode.offset().y + qNode.size().width  >= ensemble.size().width )
  //{
  //	tempAccept = false;
  //}
  //if (mode==CODING_MODE_PQI)
  //	tempAccept = false;//switch off TP
  //if ( tempAccept) // true means size > 16, not the very first block
  //{
  ////	//qNode.Display(1);

  //	matchCandid = BoundaryMatching(qNode);
  //	if (matchCandid[candidNum-1] == Point3i(-1,-1,-1) ) // no match find
  //		tempAccept = false;
  //	else
  //	{
  //		//rst_with_seam.WriteBlock("tempRst_seam.tif");
  //		int index = IsAcceptPredict(matchCandid,qNode,COMPARE_CRITERIA_SSIM,3,4);
  //		if (index > 0)
  //		{
  //
  //			stat.numTextureInLevel[qLevel]++;
  //			stat.textureBitInLevel[qLevel]+= (int)ceil(log((double)this->candidNum)/log(2.0));
  //		//	rst.SaveBlock("intermediate.tif");
  //			if ( matchCandid[index].x+2*qNode.overlap().height+qNode.size().height < ensemble.size().height &&
  //				 matchCandid[index].y+2*qNode.overlap().width+qNode.size().width < ensemble.size().width &&
  //				 qNode.GetFootPos().x+qNode.overlap().height < ensemble.size().height &&
  //				 qNode.GetFootPos().y+qNode.overlap().width  < ensemble.size().width )
  //				stat.tprecs.push_back(TPrecord(matchCandid[index]+qNode.overlap().Point3(),Cube(qNode.offset(),qNode.size()),qNode.overlap()));
  //
  //			QNode<T,cn> target(rst,qNode.size(),matchCandid[index] + qNode.overlap().Point3(),qNode.overlap());
  //

  //		/*	if (qNode.offset() == Point3i(96,32,0))
  //			{
  //				rst.Display();
  //				target.Display();
  //				target.GetBoundaryUp().Display();
  //				target.GetBoundaryLeft().Display();
  //				qNode.Display();
  //				qNode.GetBoundaryUp().Display();
  //				qNode.GetBoundaryLeft().Display();
  //			}*/
  //
  //			if ( qNode.GetBoundaryUp().IsInside(target.GetFootPos()) || qNode.GetBoundaryLeft().IsInside(target.GetFootPos()))
  //				target = QNode<T,cn>(target); //if the target overlap with qNode's boundary, make a copy to avoid data loss during quilting
  //
  //			qNode.Quilting(target);

  //		/*	if (qNode.offset() == Point3i(96,32,0))
  //			{
  //				rst.Display();
  //				target.Display();
  //				target.GetBoundaryUp().Display();
  //				target.GetBoundaryLeft().Display();
  //				qNode.Display();
  //				qNode.GetBoundaryUp().Display();
  //				qNode.GetBoundaryLeft().Display();
  //			}*/
  //			//qNode.Print();
  //			//rst.Display();
  //			//qNode.GetBoundaryUp().Print();
  //			//qNode.GetBoundaryLeft().Print();
  //			//rst.SetBlock(qNode.offset(),target);
  //			//for (int k = 0; k< qNode.GetBoundarySize(); k++)
  //			//	rst.SetBlock(qNode.offset() - qNode.overlap().Point3(),qNode.GetBoundary(k));
  //			SetCausalMap(qNode);
  //			rst_with_seam.SetBlock(qNode.offset(),target.CvtColor<3>(CV_GRAY2RGB));
  //			//deal with seam
  //			Tensor_<T,3> tempBound = qNode.GetBoundaryUp().CvtColor<3>(CV_GRAY2RGB);
  //			Vec<T,3> tempElm(0,0,255);
  //			for (unsigned int b=0; b < qNode.seam[0].size(); b++)
  //			{
  //				tempBound[Point3i(qNode.seam[0][b],b,0)] = tempElm;
  //			}
  //			rst_with_seam.SetBlock(qNode.offset() - qNode.overlap().Point3(),tempBound);
  //			tempBound  = qNode.GetBoundaryLeft().CvtColor<3>(CV_GRAY2RGB);
  //			for (unsigned int b=0; b < qNode.seam[1].size(); b++)
  //			{
  //				tempBound[Point3i(b,qNode.seam[1][b],0)] = tempElm;
  //			}
  //			rst_with_seam.SetBlock(qNode.offset() - qNode.overlap().Point3(),tempBound);
  //			stat.psnrTextureInLevel[qLevel] += qNode.ComputeMSE(ensemble.Crop(qNode.offset(),qNode.size()));
  //			//rst_with_seam.Display();
  //			return;
  //		}
  //		else
  //			tempAccept = false;
  //	}
  //}
  ///////////////////////////////////////////////////////
  //if tempAccept == false, that means either size < 16, qNode.offet = (0,0,0) or the texture preduction failed
  if ( !TexturePrediction(qNode,qLevel)) // if texture prediction fail, try smooth predict
    {
      //qNode.SetBlock(Point3i(0,0,0),ensemble.Crop(qNode.offset(),qNode.size()));
#ifdef DEBUG
      if(qLevel > 3 && debug_switch)
        {
          qNode.Print();
          qNode.GetBoundaryUp().Print();
          qNode.GetBoundaryLeft().Print();
        }
#endif

      qNode.LinearInterp(qstep[qLevel]);
#ifdef DEBUG
      if(qLevel > 3 && debug_switch)
        {
          qNode.Print();
          qNode.GetBoundaryUp().Print();
          qNode.GetBoundaryLeft().Print();
          ensemble.Crop(qNode.offset(),qNode.size()).Print();
        }
#endif
      //mylib::DisplayMat(interpNode[0],"interp");
      //LinearInterp Changed qNode!
      if (IsAcceptPredict(qNode,CompareCriteria::COMPARE_CRITERIA_INTERP,thrd[qLevel]))
        //if (qNode.Compare(CropFromOrg(qNode.offset,qNode.size),COMPARE_CRITERIA_INTERP_1,param.thrd[qLevel]))
        {
#ifdef DEBUG
	  if (qLevel == 0)
	    {
	      debug_switch = true;
	      IsAcceptPredict(qNode,CompareCriteria::COMPARE_CRITERIA_INTERP,thrd[qLevel]);
	    }
#endif
	  //Display();
	  TreatAcceptedInterp(qNode,qLevel);

	  //	rst.SaveBlock("intermediate_rst.pgm");
	  //this->DisplayCasualMap();
	  //	Display(0);
#ifdef DEBUG
	  if (debug_switch = true)
	    debug_switch = false;
#endif
	  return;
	}
      else //then both texture and smooth prediction failed,
        {
          int tempMark = 0;
          if (qNode.size().height==2 && qNode.size().width==2)// do rectangular spliting
            {
              tempMark =	TreatRectSplit(qNode, qLevel+1);
              if (tempMark < 2)
                {
                  ///	rst.SaveBlock("intermediate_rst.pgm");
                  //Display();
                  return;
                }
              //rst.SaveBlock("intermediate_rst.pgm");
              //DisplayCausalMap(0);
              //Display(0);

            }
#ifdef DEBUG
	  if (tempMark == 2 && debug_switch)
	    {
	      qNode.Print();
	      rst.Crop(qNode.offset() - Point3i(1,1,0),qNode.size()+Size3(1,1,0)).Print();
	      qNode.GetBoundaryUp().Print();
	      qNode.GetBoundaryLeft().Print();
	      ensemble.Crop(qNode.offset() - Point3i(1,1,0),qNode.size()+Size3(1,1,0)).Print();
	    }
#endif
	  qNode.Split();

	  if (tempMark == 2) // do 1x1
	    {
	      qLevel++;
	    }


	  //stat.treebits.push_back("1");
	  //stat.treebitLength++;
	  //stat.treeBitInLevel[qLevel]++;

	  tempTree = qNode.NextLeaf();
	  while (true)
	    {
	      if (tempTree == NULL)
		break;
	      else
		{
		  ///update overlap size information, by param.overlapRatio
		  CodingMTC(*tempTree,qLevel+1);//QNode<T,cn>(rst,tempTree->size(),tempTree->offset(),qNode.overlap()/2),qLevel+1);
		}
	      if (tempTree->GetTreePos()==tempTree->GetPeerTreeSize()-1)
		break;
	      else
		tempTree = tempTree->NextLeaf();
	    }
	}
    }
}
void MTC::PostCoding(QTree<T,cn> &qNode, int qLevel)
{
  QTree<T,cn> * tempTree;
  if (!TexturePrediction(qNode,qLevel))
    {
      if (qNode.size().height <= 16 && qNode.size().width <= 16)
        {
          SetCausalMap(qNode);
          //causalMap.Display(1000,1);
          return;
        }
      else
        {
          qNode.Split();
          tempTree = qNode.NextLeaf();
          while (true)
            {
              if (tempTree == NULL)
                break;
              else
                {
                  ///update overlap size information, by param.overlapRatio
                  PostCoding(*tempTree,qLevel+1);//QNode<T,cn>(rst,tempTree->size(),tempTree->offset(),qNode.overlap()/2),qLevel+1);
                }
              if (tempTree->GetTreePos()==tempTree->GetPeerTreeSize()-1)
                break;
              else
                tempTree = tempTree->NextLeaf();
            }
        }
    }
  else
    {
      //rst_with_seam.Display();
      //rst.Display();
    }

}
void MTC::UpdatePBSet(const QNode<T,cn>& qNode, Point3i matchPos )
{
  matchPos = matchPos+qNode.overlap().Point3();
  //pbq.push(Cube());
  // the new TB block coincide with existing record, erase it
  pbSet.erase(PBRecord(qNode.offset(),BoundDir::LEFT),Cube(matchPos+Point3i(0,qNode.size().width,0), qNode.leftBound.size() - Size3(qNode.overlap().height,0,0)));
  pbSet.erase(PBRecord(qNode.offset(),BoundDir::UP),Cube(matchPos+Point3i(qNode.size().height,0,0),qNode.upBound.size()-Size3(0,qNode.overlap().width,0)));
  //if the new TB block coincide with the offset and and direction, but differnet size
  pbSet.insert(pair<PBRecord,Cube>(PBRecord(qNode.offset()+Point3i(0,qNode.size().width,0),BoundDir::LEFT),Cube(matchPos+Point3i(0,qNode.size().width,0),qNode.leftBound.size() - Size3(qNode.overlap().height,0,0))));
  pbSet.insert(pair<PBRecord,Cube>(PBRecord(qNode.offset()+Point3i(qNode.size().height,0,0),BoundDir::UP),Cube(matchPos+Point3i(qNode.size().height,0,0),qNode.upBound.size() - Size3(0,qNode.overlap().width,0))));
}



void MTC::UpdatePBSet(const QNode<T,cn>& qNode, const QNode<T,cn>& candid )//the second part of pbSet store the actual canndidate post bourdary location and size
{
  Point3i matchPos = candid.offset();
  //cout<<"matchpos is: "<<matchPos<<endl;
  //cout<<"candid is: "<<candid.offset()<<endl;
  //cout<<"candid bound are: " <<candid.rightBound.offset()<<", "<<candid.lowBound.offset()<<endl;
  //pbq.push(Cube());
  // the new TB block coincide with existing record, erase it
  pbSet.erase(PBRecord(qNode.offset(),BoundDir::LEFT),Cube(matchPos+Point3i(0,qNode.size().width,0), qNode.rightBound.size() - Size3(qNode.overlap().height*2,0,0)));
  pbSet.erase(PBRecord(qNode.offset(),BoundDir::UP),Cube(matchPos+Point3i(candid.size().height,0,0),qNode.lowBound.size()-Size3(0,qNode.overlap().width*2,0)));
  //if the new TB block coincide with the offset and and direction, but differnet size
  //01052013 take the entire boundary will be usefull for post PLC.
  //05152013 not only extend the long boudnary, but also take extra short bournday
  //05152013 for example a 4x16 boundary will store the extra 6x26 boundary?????
  //Tensor<double,1> candExt = candid.GetExtendTensor();
  Tensor<double,1> aRecord = candid.rightBound.Crop(Point3i(candid.overlap().height,0,0),qNode.rightBound.size()-Size3(qNode.overlap().height*2,0,0));
  //cout<<"crop "<<candid.rightBound.Crop(Point3i(candid.overlap().height,0,0),qNode.rightBound.size()-Size3(qNode.overlap().height*2,0,0)).offset()<<endl;
  //cout<<"a record "<<aRecord.offset()<<aRecord.size()<<endl;
  CubePlus temp(aRecord);
  //cout<<"record building 1: "<<temp.offset()<<", "<<temp.size()<<endl;
  //  temp.SetExtraContent(candExt.Crop(Point3i(candid.overlap().height-1,candExt.size().width-candid.overlap().width-1,0),  qNode.rightBound.size()+Size3(2-qNode.overlap().height*2,2,0)));
  pbSet.insert(pair<PBRecord,CubePlus>(PBRecord(qNode.offset()+Point3i(0,qNode.size().width,0),BoundDir::LEFT),temp/*CubePlus(candid.rightBound.Crop(Point3i(qNode.overlap().height,0,0),candid.rightBound.size()-Size3(qNode.overlap().height*2,0,0)))*/));
  Tensor<double,1> bRecord = candid.lowBound.Crop(Point3i(0,candid.overlap().width,0),qNode.lowBound.size()-Size3(0,qNode.overlap().width*2,0));
  temp = CubePlus(bRecord);
  //temp.SetExtraContent(candExt.Crop(Point3i(candExt.size().height-candid.overlap().height-1, candid.overlap().width-1,0),qNode.lowBound.size()+Size3(2,2-qNode.overlap().width*2,0)));
  //cout<<"record building 2: "<<temp.offset()<<", "<<temp.size()<<endl;
  pbSet.insert(pair<PBRecord,CubePlus>(PBRecord(qNode.offset()+Point3i(qNode.size().height,0,0),BoundDir::UP),temp/*CubePlus(candid.lowBound.Crop(Point3i(0,qNode.overlap().width,0),candid.lowBound.size()-Size3(0,qNode.overlap().width*2,0)))*/));
}




bool MTC::IsAcceptPredict(QNode<T,cn>& qNode, CompareCriteria criteria, double param1, double param2)
{
  if (qNode.size().height == 1 && qNode.size().width ==1)
    return true;
  Tensor<T,cn> temp;
  ensemble.Ref(Cube(qNode.offset(),qNode.size()),temp);
  double distance = metric::Compare(qNode,temp,criteria,this->subSize, this->subStep, param1,param2);
  bool accept;
  if ( criteria != CompareCriteria::COMPARE_CRITERIA_INTERP)
    {
      if (distance <= this->qualityThrd)
        accept = true;
      else
        accept = false;
    }
  else
    {
      if (distance < 0.5)
        accept = false;
      else
        accept = true;
    }

#ifdef DEBUG
  if (debug_switch)
    {
      rst.SaveBlock("intermediate_rst.pgm");
      ensemble.Crop(qNode.offset(),qNode.size()).Print();
      qNode.GetBoundaryUp().Print();
      qNode.GetBoundaryLeft().Print();
      qNode.AbsDiff(ensemble.Crop(qNode.offset(),qNode.size())).CompareElement(param1).Print();
      qNode.ComputeAIM(ensemble.Crop(qNode.offset(),qNode.size()),param1);
    }
#endif
  return accept;
}

void MTC::TreatAcceptedInterp(QTree<T,cn>& qNode, int qLevel)
{
  //stat.treebits.push_back("0");
  //stat.treebitLength++;
  //stat.treeBitInLevel[qLevel]++;
  //stat.footbits.push_back(qNode.GetBits());//for blackwhite only!
  //stat.footbitLength += qNode.BitLength();
  stat.footBitInLevel[qLevel] += qNode.BitLength();
  //stat.blockSize.push_back(qNode.size);
  //stat.footPos.push_back(qNode.footPos);
  stat.numBlocksInLevel[qLevel] ++;
  stat.psnrInLevel[qLevel] += metric::ComputeMSE(qNode, ensemble.Crop(qNode.offset(),qNode.size()));
  //rst.SetBlock(qNode.offset(),qNode);
  //SetRst(qNode);
#if CV_MINOR_VERSION <6
  rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(CV_GRAY2RGB));
#else
  rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(COLOR_GRAY2RGB));

#endif
  footmap[qNode.GetFootPos()] =0;
  //rst_with_seam.Display();
  SetCausalMap(qNode);
  return;
}

int MTC::TreatRectSplit(QTree<T,cn>& qNode, int qLevel)
{

  if (rectType) //rect CASE 1 //method as in Teng's thesis
    {
      //qNode.Print();
      Tensor<T,cn> org = ensemble.Crop(qNode.offset(),qNode.size());
      Tensor<uchar,1> map = metric::CompareElement(qNode.AbsDiff(org), thrd[qLevel]);
      //	map.Print();//
      QNode<T,cn> newBlock(qNode);
      //		qNode.GetBoundaryUp().Print();
      //		qNode.GetBoundaryLeft().Print();
      stat.numTestInLevel[qLevel]+=2;
      if (map[Point3i(0,1,0)][0]>0 && map[Point3i(1,1,0)][0]>0)
        {
          //output 2x2 foot
          newBlock[Point3i(1,0,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryLeft(Point3i(2,0,0))+qNode[Point3i(1,1,0)]);
          newBlock[Point3i(1,0,0)] = newBlock[Point3i(1,0,0)] + newBlock.Quantize(org[Point3i(1,0,0)] - newBlock[Point3i(1,0,0)],qstep[qLevel]);
#ifdef INT_RST
	  //for all int value result
	  for (int i=0; i<cn; i++)
	    newBlock[Point3i(1,0,0)][i] = T(int(newBlock[Point3i(1,0,0)][i]));
#endif
	  newBlock[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryUp(Point3i(0,1,0))+newBlock[Point3i(1,0,0)]);
	  //newBlock.Print();
	  if (IsAcceptPredict(newBlock,CompareCriteria::COMPARE_CRITERIA_INTERP,thrd[qLevel]))
	    {
	      stat.footBitInLevel[qLevel] += qNode.BitLength();
	      //stat.treeBitInLevel[qLevel] ++;
	      stat.footBitInLevel[qLevel] += newBlock.BitLength();
	      //stat.treeBitInLevel[qLevel]++;
	      qNode = newBlock;
	      qNode.Split(0); // vertical
	      stat.numBlocksInLevel[qLevel]+=2;
	      stat.psnrInLevel[qLevel] += metric::ComputeMSE(qNode, ensemble.Crop(qNode.offset(),qNode.size()));
#if CV_MINOR_VERSION <6
	      rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(CV_GRAY2RGB));
#else
	      rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(COLOR_GRAY2RGB));
#endif
	      rst.SetBlock(qNode.offset(),qNode);
	      SetCausalMap(qNode);
	      footmap[qNode.GetFootPos()] = 0;
	      footmap[qNode.GetFootPos() - Point3i(0,1,0)] = 0;
	      return 0;
	    }
	}
      newBlock = qNode;
      if (map[Point3i(1,0,0)][0]>0 && map[Point3i(1,1,0)][0]>0)
        {
          newBlock[Point3i(0,1,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryUp(Point3i(0,2,0))+qNode[Point3i(1,1,0)]);
          newBlock[Point3i(0,1,0)] = newBlock[Point3i(0,1,0)] + newBlock.Quantize(org[Point3i(0,1,0)] - newBlock[Point3i(0,1,0)],qstep[qLevel]);
#ifdef INT_RST
	  //for all int value result
	  for (int i=0; i<cn; i++)
	    newBlock[Point3i(0,1,0)][i] = T(int(newBlock[Point3i(0,1,0)][i]));
#endif
	  newBlock[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryLeft(Point3i(1,0,0))+newBlock[Point3i(0,1,0)]);
	  //newBlock.Print();
	  if (IsAcceptPredict(newBlock, CompareCriteria::COMPARE_CRITERIA_INTERP,thrd[qLevel]))
	    {
	      stat.footBitInLevel[qLevel] += qNode.BitLength();
	      stat.footBitInLevel[qLevel] += newBlock.BitLength();
	      //stat.treeBitInLevel[qLevel]++;
	      //stat.treeBitInLevel[qLevel] ++;

	      qNode = newBlock;
	      qNode.Split(1); // horizontal
	      stat.numBlocksInLevel[qLevel]+=2;
	      stat.psnrInLevel[qLevel] += metric::ComputeMSE(qNode, ensemble.Crop(qNode.offset(),qNode.size()));
#if CV_MINOR_VERSION <6
	      rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(CV_GRAY2RGB));
#else
	      rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(COLOR_GRAY2RGB));
#endif
	      rst.SetBlock(qNode.offset(),qNode);
	      SetCausalMap(qNode);
	      footmap[qNode.GetFootPos()] = 0;
	      footmap[qNode.GetFootPos() - Point3i(1,0,0)] = 0;
	      return 1;
	    }
	}

      newBlock = qNode;
      {
        //manually do 4 point
        //newBlock[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryUp(Point3i(0,1,0)) + qNode.GetBoundaryLeft(Point3i(1,0,0)));
        //newBlock[Point3i(0,1,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryUp(Point3i(0,2,0)) + newBlock[Point3i(0,0,0)]);
        //newBlock[Point3i(1,0,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryUp(Point3i(2,0,0)) + newBlock[Point3i(0,0,0)]);
        //newBlock[Point3i(1,1,0)] = Vec<T,cn>::all(0.5).mul(newBlock[Point3i(0,1,0)] + newBlock[Point3i(1,0,0)]);
        //qNode = newBlock;
        //qNode.Split();
        return 2;
      }
    }
  else	 //rect CASE 2
    {


#ifdef BOUNDARY_ADT

      Tensor<T,cn> org;
      ensemble.Ref(Cube(qNode.offset(),qNode.size()),org);
      //org.Print();
      //qNode.Print();
      QNode<T,cn> newBlock1(qNode);
      stat.numTestInLevel[qLevel]+=2;
      //newBlock1.GetBoundaryLeft().Print();
      //newBlock1.GetBoundaryUp().Print();
      newBlock1[Point3i(1,0,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryLeft(Point3i(2,0,0))+qNode[Point3i(1,1,0)]);
      newBlock1[Point3i(1,0,0)] = newBlock1[Point3i(1,0,0)] + newBlock1.Quantize(org[Point3i(1,0,0)] - newBlock1[Point3i(1,0,0)],qstep[qLevel]);
      newBlock1[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryUp(Point3i(0,1,0))+newBlock1[Point3i(1,0,0)]);
      //newBlock1.Print();
      QNode<T,cn> newBlock2(qNode);
      //newBlock2.GetBoundaryLeft().Print();
      //newBlock2.GetBoundaryUp().Print();
      newBlock2[Point3i(0,1,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryUp(Point3i(0,2,0))+qNode[Point3i(1,1,0)]);
      newBlock2[Point3i(0,1,0)] = newBlock2[Point3i(0,1,0)] + newBlock2.Quantize(org[Point3i(0,1,0)] - newBlock2[Point3i(0,1,0)],qstep[qLevel]);
      newBlock2[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(qNode.GetBoundaryLeft(Point3i(1,0,0))+newBlock2[Point3i(0,1,0)]);
      //newBlock2.Print();
      double MSE1 = org.ComputeMSE(newBlock1);
      double MSE2 = org.ComputeMSE(newBlock2);
      if (MSE1 <= MSE2)
        {
          stat.footBitInLevel[qLevel] += qNode.BitLength();
          stat.footBitInLevel[qLevel] += newBlock1.BitLength();
          qNode.SetBlock(newBlock1);
          qNode.Split(0); // vertical
          stat.numBlocksInLevel[qLevel]+=2;
          stat.psnrInLevel[qLevel] += MSE1;
          rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(CV_GRAY2RGB));
          //rst.SetBlock(qNode.offset(),qNode);
          SetCausalMap(qNode);
          footmap[qNode.GetFootPos()] = 0;
          footmap[qNode.GetFootPos() - Point3i(0,1,0)] = 0;

          return 0;
        }
      else
        {
          stat.footBitInLevel[qLevel] += qNode.BitLength();
          stat.footBitInLevel[qLevel] += newBlock2.BitLength();
          qNode.SetBlock(newBlock2);
          qNode.Split(1); // horizontal
          stat.numBlocksInLevel[qLevel]+=2;
          stat.psnrInLevel[qLevel] += MSE2;
          rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(CV_GRAY2RGB));
          //rst.SetBlock(qNode.offset(),qNode);
          SetCausalMap(qNode);
          footmap[qNode.GetFootPos()] = 0;
          footmap[qNode.GetFootPos() - Point3i(1,0,0)] = 0;

          return 1;
        }
#else
      Tensor<T,cn> newBlock1(Size3(2,2,1));
      if ( qNode.offset().x+1 == ensemble.size().height || qNode.offset().y==0)
        newBlock1[Point3i(1,0,0)] = Vec<T,cn>::all(0.5).mul(Vec<T,cn>::all(128) + qNode[Point3i(1,1,0)]);
      else
        newBlock1[Point3i(1,0,0)] = Vec<T,cn>::all(0.5).mul(rst[Point3i(qNode.offset().x+1,qNode.offset().y-1,0)]+qNode[Point3i(1,1,0)]);
      if (qNode.offset().x != 0)
        newBlock1[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(rst[Point3i(qNode.offset().x-1,qNode.offset().y,0)]+newBlock1[Point3i(1,0,0)]);
      else
        newBlock1[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(Vec<T,cn>::all(128)+newBlock1[Point3i(1,0,0)]);


      newBlock1[Point3i(0,1,0)] = qNode[Point3i(0,1,0)];
      newBlock1[Point3i(1,1,0)] = qNode[Point3i(1,1,0)];
      //mylib::DisplayMat(newBlock1[0],"block 1");
      Tensor<T,cn> newBlock2(Size3(2,2,1));
      if (qNode.offset().x == 0 || qNode.offset().y+1 == rst.size().width)
        newBlock2[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(Vec<T,cn>::all(128)+newBlock2[Point3i(1,0,0)]);
      else
        newBlock2[Point3i(0,1,0)] = Vec<T,cn>::all(0.5).mul(rst[Point3i(qNode.offset().x-1,qNode.offset().y+1,0)]+qNode[Point3i(1,1,0)]);
      if (qNode.offset().y == 0)
        newBlock2[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(Vec<T,cn>::all(128)+newBlock2[Point3i(0,1,0)]);
      else
        newBlock2[Point3i(0,0,0)] = Vec<T,cn>::all(0.5).mul(rst[Point3i(qNode.offset().x,qNode.offset().y-1,0)]+newBlock2[Point3i(0,1,0)]);
      newBlock2[Point3i(1,1,0)] = qNode[Point3i(1,1,0)];
      newBlock2[Point3i(1,0,0)] = qNode[Point3i(1,0,0)];
      //mylib::DisplayMat(newBlock2[0], "block 2");
      double MSE1 = metric::ComputeMSE(newBlock1, ensemble.Crop(qNode.offset(),qNode.size()));
      double MSE2 = metric::ComputeMSE(newBlock2, ensemble.Crop(qNode.offset(),qNode.size()));
      if (MSE1 < MSE2) // vertical split
        {
          qNode.SetBlock(Point3i(0,0,0),newBlock1);
          qNode.Split(0);// 0 vertical  1 horizontal
          stat.treebits.push_back("110");
          stat.treebitLength+=2;
          stat.treeBitInLevel[qLevel+1]+=2;
          //treeStream<<"110";
          //BitsTree+=3;
          //stat.blockSize.push_back(BlockSize<int>(1,2,1));
          //stat.blockSize.push_back(BlockSize<int>(1,2,1));
          //stat.footbits.push_back(
        }
      else
        {
          qNode.SetBlock(Point3i(0,0,0),newBlock2);
          qNode.Split(1);
          stat.treebits.push_back("111");
          stat.treeBitInLevel[qLevel+1]+=2;
          //treeStream<<"111";
          stat.treebitLength+=2;
          //BitsTree+=3;
          //stat.blockSize.push_back(BlockSize<int>(2,1,1));
          //stat.blockSize.push_back(BlockSize<int>(2,1,1));
        }
#endif
    }
  if(qNode.GetSubTreeSize() <4)//if rect spliting
    {
      stat.numTestInLevel[qLevel+1]++;
      stat.footbits.push_back(qNode.GetBits());
      stat.footbitLength+= qNode.BitLength();
      stat.footBitInLevel[qLevel+1]+=qNode.BitLength();
      //stat.footPos.push_back(qNode.offset);
      stat.numBlocksInLevel[qLevel+1] +=2;
      //SetRst(qNode);
#if CV_MINOR_VERSION <6
      rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(CV_GRAY2RGB));
#else
      rst_with_seam.SetBlock(qNode.offset(),qNode.CvtColor<3>(COLOR_GRAY2RGB));

#endif
      rst.SetBlock(qNode.offset(),qNode);
      SetCausalMap(qNode);
    }
  return 0;

}


void MTC::PostBlending(const vector<TPrecord>& tprecs)
{
  for (unsigned int b=0; b<tprecs.size(); b++)
    {
      vector<Tensor<T,cn>> postBoundsA(2); //after target position, format, below, then right
      vector<Tensor<T,cn>> postBoundsB(2); //after predict position
      TPrecord rec = tprecs[b];
      rst.Ref(Cube(rec.targetROI.offset()+Point3i(rec.targetROI.height,0,0),Size3(rec.overlapSize.height,rec.targetROI.width+rec.overlapSize.width,1)),postBoundsA[0]);
      rst.Ref(Cube(rec.targetROI.offset()+Point3i(0,rec.targetROI.width,0),Size3(rec.targetROI.height+rec.overlapSize.height,rec.overlapSize.width,1)),postBoundsA[1]);
      rst.Ref(Cube(rec.matchPos+Point3i(rec.targetROI.height,0,0),Size3(rec.overlapSize.height,rec.targetROI.width+rec.overlapSize.width,1)),postBoundsB[0]);
      rst.Ref(Cube(rec.matchPos+Point3i(0,rec.targetROI.width,0),Size3(rec.targetROI.height+rec.overlapSize.height,rec.overlapSize.width,1)),postBoundsB[1]);
      QNode<T,cn> dummy;
      dummy.Quilting(postBoundsA,postBoundsB);
#if CV_MINOR_VERSION <6

      Tensor<T,3> tempBound = postBoundsA[0].CvtColor<3>(CV_GRAY2RGB);
#else
      Tensor<T,3> tempBound = postBoundsA[0].CvtColor<3>(COLOR_GRAY2RGB);
#endif

      Vec<T,3> tempElm(0,255,0);
      for (unsigned int k=0; k < dummy.seam[0].size(); k++)
        {
          tempBound[Point3i(dummy.seam[0][k],k,0)] = tempElm;
        }
      rst_with_seam.SetBlock(postBoundsA[0].offset(),tempBound);
#if CV_MINOR_VERSION <6

      tempBound  = postBoundsA[1].CvtColor<3>(CV_GRAY2RGB);
#else
      tempBound  = postBoundsA[1].CvtColor<3>(COLOR_GRAY2RGB);
#endif
      for (unsigned int k=0; k < dummy.seam[1].size(); k++)
        {
          tempBound[Point3i(k,dummy.seam[1][k],0)] = tempElm;
        }
      rst_with_seam.SetBlock(postBoundsA[1].offset(),tempBound);
    }
}
void MTC::PostBlending(const PBRecord& rec, CubePlus& roi)
{
  //i knew all the block size is a power of 2 such that I can do this, but this is not generic
  //Size3 sz(roi.height+roi.width,roi.height+roi.width,roi.depth);
  //Tensor<T,cn> candid = rst.Crop(rec.offset,sz);
  //Tensor<T,cn> target;
  //rst.Ref()
  vector<Tensor<T,cn>> tagBound(1);
  vector<Tensor<T,cn>> refBound(1);
  //tagBound[0] = rst.Crop(roi.offset(),roi.size());
  //rst.Ref(roi,tagBound[0]);//tag is the predicted block(candid)
  rst.Ref(Cube(rec.offset,roi.size()),refBound[0]);//ref is the block to be predict(target)
  tagBound[0] = roi.GetContent().Clone();

  //lighting correction
  //tagBound[0].LightingCorrection(refBound[0]);
  QNode<T,cn> dummy;
  //Tensor<T,cn> tempTarget,tempRef; //tempTarget is the MTC coded block + PLC corrected boundary , which will be correced
  //                                 //tempRef is the MTC coded block + JPEG coded boundary, which will be used to correct lighting
  Tensor<T,cn> changeFrom, changeTo;
  //20130605 testing ignore the lighting discontiuity , just force to do PLC on blk + bd
  cout<<rec.offset<<"size"<<roi.GetContent().size().height<<","<<roi.GetContent().size().width<<", dirc";
  cout<<(int)rec.direction<<endl;
  int bsize = max(roi.height,roi.width);
  Point3i pos;
  Size3 sz;
  //rst.Ref(Cube(rec.offset-Point3i(1,1,0),roi.GetContent().size()+Point3i(2,2,0)),changeTo);//this is PLC on post boundary only
  //do blk plc first
  if (rec.direction ==BoundDir::LEFT)
    {
      //pos = rec.offset - Point3i(1,1,0) - Point3i(0,bsize,0);
      sz  =  Size3(0,bsize,0);
    }
  else
    {
      //pos = rec.offset - Point3i(1,1,0) - Point3i(bsize,0,0);
      sz = Size3(bsize,0,0);
    }
  //rst.Ref(Cube(rec.offset - Point3i(1,1,0) - sz.Point3(), roi.size() + Point3i(2,2,0) + sz),changeTo);//this is PLC with blk + post bd


  //// //then do pb plc

  rst.Ref(Cube(rec.offset-Point3i(1,1,0),roi.GetContent().size()+Point3i(2,2,0)),changeTo);//this is PLC on post boundary only
  if (rec.direction == BoundDir::LEFT||rec.direction==BoundDir::RIGHT)
    {
      //fill the last row of changeTo as the last row of JPEG boundary
      for (int i=0; i<changeTo.size().width; i++)
        {
          changeTo(changeTo.size().height-1,i,0) = changeTo(changeTo.size().height-2,i,0);
        }
    }
  else
    {
      //fill the last row of changeTo as the last row of JPEG boundary
      for (int i=0; i<changeTo.size().height; i++)
        {
          changeTo(i,changeTo.size().width-1,0) = changeTo(i,changeTo.size().width-2,0);
        }
    }
  //changeFrom = changeTo.Clone();
  cout<<roi.offset()<<roi.GetContent().size()<<endl;
  rst.Ref(Cube(roi.offset() - Point3i(1,1,0),roi.GetContent().size()+Point3i(2,2,0)),changeFrom);
  if (lightCorrectionType == LightingCorrectionType::POISSON_LC)
    {
      Tensor<T,cn> postPLCRst = PostPLC(changeFrom,changeTo);//from xxx change to xxx
      tagBound[0] = postPLCRst.Crop(Point3i(1,1,0),tagBound[0].size());
      //   //tagBound[0] = postPLCRst.Crop(Point3i(1,1,0)+sz.Point3(),roi.size());
    }
  dummy.Quilting(refBound,tagBound,rec.direction);
  ////refBound[0].Print();

  Tensor<T,3> tempBound = refBound[0].CvtColor<3>(COLOR_GRAY2RGB);
  Vec<T,3> tempElm(0,255,0);
  if(rec.direction == BoundDir::UP)
    {
      for (unsigned int k=0; k < dummy.seam[0].size(); k++)
        {
          tempBound[Point3i(dummy.seam[0][k],k,0)] = tempElm;
        }
      rst_with_seam.SetBlock(refBound[0].offset(),tempBound);
    }
  if (rec.direction == BoundDir::LEFT)
    {
      for (unsigned int k=0; k < dummy.seam[0].size(); k++)
        {
          tempBound[Point3i(k,dummy.seam[0][k],0)] = tempElm;
        }
      rst_with_seam.SetBlock(refBound[0].offset(),tempBound);
    }
  //do relighting on the blk
  //!enable relighting
   rst.Ref(Cube(rec.offset - Point3i(1,1,0) - sz.Point3(), Size3(2,2,0) + Size3(bsize,bsize,1)),changeTo);//this is PLC with blk only
   //actual candidate location

   rst.Ref(Cube(roi.offset() - Point3i(1,1,0) - sz.Point3(), Size3(2,2,0) + Size3(bsize,bsize,1)),changeFrom);
   if (rec.direction == BoundDir::LEFT||rec.direction==BoundDir::RIGHT)
   {
     //fill the last row of changeTo as the last row of JPEG boundary
     for (int i=0; i<changeTo.size().width; i++)
     {
       changeTo(changeTo.size().height-1,i,0) = changeTo(changeTo.size().height-2,i,0);
       //changeFrom(changeTo.size().height-1,i,0) = changeTo(changeTo.size().height-2,i,0);
     }
     //fill the last col as the post boundary
     //for (int i=1; i<changeFrom.size().height-1; i++)
     //{
     //  changeFrom(i,changeFrom.size().width-1,0) = roi.GetContent()(i-1,0,0);
     //}
   }
   else
   {
      //fill the last col of changeTo as the last row of JPEG boundary
     for (int i=0; i<changeTo.size().height; i++)
     {
       changeTo(i,changeTo.size().width-1,0) = changeTo(i,changeTo.size().width-2,0);
       //changeFrom(i,changeTo.size().width-1,0) = changeTo(i,changeTo.size().width-2,0);
     }
     //fill the last row as the post boundary
     //for (int i=1; i<changeFrom.size().width-1; i++)
     //{
     //  changeFrom(changeFrom.size().height-1,i,0) = roi.GetContent()(0,i-1,0);
     //}
   }

   //change 20130616 the boundary will keep the same as target for pb
   //changeFrom = Tensor<T,cn>(changeTo.size(),0);
   //changeFrom.SetBlock(Point3i(1,1,0),changeTo.Crop(Point3i(1,1,0),changeTo.size()-Size3(2,2,0)));
   //changeFrom = changeTo.Clone();

   //changeFrom.SetBlock(Point3i(1,1,0),roi.GetContent());
   //changeFrom.SetBlock(Point3i(1,1,0)+ sz.Point3(), roi.GetContent()); //plc of blk+pb together
   //changeTo.Print("changeTo",true);
   //changeFrom.Print("changetFrom",true);

   /*
  // rst.Display();
   int bsize = max(roi.size().height,roi.size().width);
   int osize = min(roi.size().height,roi.size().width);
   int extSize = bsize+2*osize;
   if (rec.direction == LEFT)
   {
     rst.Ref(Cube(rec.offset-Point3i(0,bsize,0)-Point3i(1,1,0),Size3(roi.size().height+1,roi.size().height+roi.size().width+1,roi.size().depth)),tempRef);
     //tempRef.Print();
     tempTarget = Tensor<T,cn>(Size3(roi.size().height+1,roi.size().height+roi.size().width+1,roi.size().depth));
     //tempTarget.Print();
     tempTarget.SetBlock(rst.Crop(rec.offset-Point3i(0,bsize,0)-Point3i(1,1,0),Size3(bsize+1,bsize+1,1)));
     //tempTarget.Print();
     tempTarget.SetBlock(Point3i(0,bsize+1,0),roi.GetExtraContent().GetBlock(Cube(Point3i(osize-1,0,0),Size3(bsize+1,osize,1))).Clone());
    // tempTarget.Print();
    // tempRef.Print();
   }
   else if (rec.direction == UP)
   {
     rst.Ref(Cube(rec.offset-Point3i(bsize,0,0)-Point3i(1,1,0),Size3(roi.size().height+roi.size().width+1,roi.width+1,roi.size().depth)),tempRef);
     tempTarget = Tensor<T,cn>(Size3(roi.size().height+roi.size().width+1,roi.width+1,roi.size().depth));
     tempTarget.SetBlock(rst.Crop(rec.offset-Point3i(bsize,0,0)-Point3i(1,1,0),Size3(bsize+1,bsize+1,1)));
     tempTarget.SetBlock(Point3i(bsize+1,0,0),roi.GetExtraContent().GetBlock(Cube(Point3i(0,osize-1,0),Size3(osize,bsize+1,1))).Clone());
    // if (roi.GetExtraContent().size().height==0)
   //  {
    //   cout<<"\n==== PB record error ======"<<rec.offset.x<<", "<<rec.offset.y<<endl;
    //   roi.GetExtraContent().Print();
     //}
   }
   */
   //do PLC before blernding
   if (lightCorrectionType == LightingCorrectionType::POISSON_LC)
   {
     Tensor<T,cn> postPLCRst = PostPLC(changeFrom,changeTo);//from xxx change to xxx
     //rst.SetBlock(rec.offset - sz.Point3(), postPLCRst.Crop(Point3i(1,1,0),Size3(bsize,bsize,1))); //use for plc of blk + pb together
     rst.SetBlock(rec.offset-sz.Point3(), postPLCRst.Crop(Point3i(1,1,0),Size3(bsize,bsize,1)));
     /*
     if (rec.direction==UP)
     {
       rst.SetBlock(rec.offset-Point3i(bsize,0,0),postPLCRst.Crop(Point3i(1,1,0),Size3(bsize,bsize,1)));
       tagBound[0].SetBlock(postPLCRst.Crop(Point3i(1+bsize,1,0),refBound[0].size()));
     }
     else if (rec.direction == LEFT)
     {
       rst.SetBlock(rec.offset-Point3i(0,bsize,0),postPLCRst.Crop(Point3i(1,1,0),Size3(bsize,bsize,1)));
       tagBound[0].SetBlock(postPLCRst.Crop(Point3i(1,1+bsize,0),refBound[0].size()));
     }
     */
     //pred_postLC.SetBlock(tempRef.offset(),postPLCRst.Crop(Point3i(1,1,0),Size3(bsize,bsize,1)));
     //tagBound[0] = postPLCRst.Crop(Point3i(1,1,0),tagBound[0].size());
     //tagBound[0] = postPLCRst.Crop(Point3i(1,1,0)+sz.Point3(),roi.size());
     rst_with_seam.SetBlock(rec.offset-sz.Point3(),postPLCRst.Crop(Point3i(1,1,0),Size3(bsize,bsize,1)).CvtColor<3>(COLOR_GRAY2RGB));
     pred_postLC.SetBlock(rec.offset-sz.Point3(),postPLCRst.Crop(Point3i(1,1,0),Size3(bsize,bsize,1)));
   }

}

Tensor<MTC::T,MTC::cn> MTC::PostPLC(Tensor<T,cn>& changeFrom, Tensor<T,cn>& changeTo)
{
  Tensor<T,1> mask(changeTo.size(),0);
  mask.SetBlock(Point3i(1,1,0),Tensor<T,1>(changeTo.size()-Size3(2,2,0),255));
  //Tensor<T,cn> tarExt = target.ExtendBoundary(Size3(1,1,0));
  //Tensor<T,1> maskExt = mask.ExtendBoundary(Size3(1,1,0));
  //cout<<"("<<changeFrom.size().height<<","<<changeFrom.size().width<<")"
  //<<"=== ("<<mask.size().height<<","<<mask.size().width<<")"<<endl;
  PoissonSolver pb(changeFrom.GetFrame(0),changeTo.GetFrame(0),mask.GetFrame(0));
  Mat dst;
  pb.PostPLC(dst);
  return dst;
}
//void MTC::SetPredictQualityThrd(double qualityThrd)
//{
//	this->qualityThrd = qualityThrd;
//}
//
//void MTC::SetQFactor(double qFactor)
//{
//	this->qfactor = qFactor;
//}
//
//void MTC::SetInitQStep(double qStep)
//{
//	this->qstep = qStep;
//}
void MTC::SetSearchStep(const Size3& step)
{
  searchStep = step;
}



void MTC::InitJpegParam(int quality)
{
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  cinfo.image_width = ensemble.size().width;
  cinfo.image_height = ensemble.size().height;
  cinfo.input_components = ensemble.channels();
  if (ensemble.channels()==1)
    cinfo.in_color_space = JCS_GRAYSCALE;
  else
    cinfo.in_color_space = JCS_RGB;
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo,quality,TRUE);
  stat.jBits.push_back(0);
  stat.jBlockNum.push_back(0);
  //this is the baseline table
  if (jpegType == JPEGType::JPEG_BASELINE)
    JpegQuanTbl = vector<UINT16>(cinfo.quant_tbl_ptrs[0]->quantval,cinfo.quant_tbl_ptrs[0]->quantval+sizeof(cinfo.quant_tbl_ptrs[0]->quantval)/sizeof(UINT16));
  //this is the modified LF component table 1
  else if (jpegType == JPEGType::JPEG_ADAPTIVE)
    {
      JpegQuanTbl = vector<UINT16>(cinfo.quant_tbl_ptrs[0]->quantval,cinfo.quant_tbl_ptrs[0]->quantval+sizeof(cinfo.quant_tbl_ptrs[0]->quantval)/sizeof(UINT16));
      jpeg_set_quality(&cinfo,int((float)quality*0.4),TRUE);//change from 0.4 to 0.8 GJ, Nov 10 2012
      JpegQuanTbl_alt = vector<UINT16>(cinfo.quant_tbl_ptrs[0]->quantval,cinfo.quant_tbl_ptrs[0]->quantval+sizeof(cinfo.quant_tbl_ptrs[0]->quantval)/sizeof(UINT16));
      jpeg_set_quality(&cinfo,quality,TRUE);
      //in adaptive mode, 2 level jpegs
      stat.jBits.push_back(0);
      stat.jBlockNum.push_back(0);

    }
  else if (jpegType == JPEGType::JPEG_LF_3)
    {
      JpegQuanTbl = vector<UINT16>(cinfo.quant_tbl_ptrs[0]->quantval,cinfo.quant_tbl_ptrs[0]->quantval+sizeof(cinfo.quant_tbl_ptrs[0]->quantval)/sizeof(UINT16));
      for (int i=0; i< 8; i++)
        for (int j=0; j<8; j++)
          if (i+j>1)
            JpegQuanTbl[i*8 + j] = 10000;
    }
  else if (jpegType == JPEGType::JPEG_LF_6)
    {
      JpegQuanTbl = vector<UINT16>(cinfo.quant_tbl_ptrs[0]->quantval,cinfo.quant_tbl_ptrs[0]->quantval+sizeof(cinfo.quant_tbl_ptrs[0]->quantval)/sizeof(UINT16));
      for (int i=0; i< 8; i++)
        for (int j=0; j<8; j++)
          if (i+j>2)
            JpegQuanTbl[i*8 + j] = 10000;
    }
  else if (jpegType == JPEGType::JPEG_CUSTOM_TBL)
    {
      JpegQuanTbl = vector<UINT16>(64,10000);
      JpegQuanTbl[0]=80;
      JpegQuanTbl[1] = 160;
      JpegQuanTbl[8] = 160;
    }
  else
    {
      //this is watson's
      double defectNum = -0.2 * double(quality) + 20;
      UINT16 temptbl_y[] ={ 10 , 7,   7,   9,  13,  17,  22,  29,
                            7,   8,   6,   7,  10,  13,  16,  21,
                            7,   6,  11,  10,  12,  15,  19,  24,
                            9,   7,  10,  17,  17,  19,  23,  28,
                            13,  10, 12,  17,  24,  26,  29,  34,
                            17,  13, 15,  19,  26,  35,  38,  43,
                            22,  16, 19,  23,  29,  38,  53,  55,
                            29,  21, 24,  28,  34,  43,  55,  70};
      UINT16 temptbl_cr[] ={ 14,  14,  29,  36,  49,  65,  86, 112 ,
                             14,   9,  20,  29,  37,  49,  63,  82 ,
                             29,  20,  34,  41,  47,  57,  72,  91 ,
                             36,  29,  41,  63,  64,  73,  87, 106 ,
                             49,  37,  47,  64,  92,  98, 111, 129 ,
                             65,  49,  57,  73,  98, 135, 145, 162 ,
                             86,  63,  72,  87, 111, 145, 201, 207 ,
                             112,  82,  91, 106, 129, 162, 207, 255};
      UINT16 temptb_cb[] = { 28,  28,  57,  93, 124, 166, 219, 255,
                             28,  17,  38,  75,  95, 124, 162, 209,
                             57,  38,  67, 104, 120, 147, 184, 231,
                             93,  75, 104, 162, 165, 187, 223, 255,
                             124,  95, 120, 165, 235, 251, 255, 255,
                             166, 124, 147, 187, 251, 255, 255, 255,
                             219, 162, 184, 223, 255, 255, 255, 255,
                             255, 209, 231, 255, 255, 255, 255, 255};
      JpegQuanTbl =  vector<UINT16>(temptbl_y, temptbl_y + sizeof(temptbl_y) / sizeof(UINT8));
      if (defectNum >0)
        for (unsigned int i=0; i< JpegQuanTbl.size(); i++)
          JpegQuanTbl[i]  *= (unsigned short)defectNum;
    }
  //vector<UINT16> QTable2(cinfo.quant_tbl_ptrs[1]->quantval,cinfo.quant_tbl_ptrs[1]->quantval+sizeof(cinfo.quant_tbl_ptrs[1]->quantval)/sizeof(UINT16));
  JpegHuffTbl_DC = vector<vector<UINT8>>(sizeof(cinfo.dc_huff_tbl_ptrs[0]->bits)-1);
  JpegHuffTbl_AC = vector<vector<UINT8>>(sizeof(cinfo.ac_huff_tbl_ptrs[0]->bits)-1);

  //	,HTableAC2,HTableDC1,HTableDC2;
  int count = 0;
  for (unsigned int i=0; i < JpegHuffTbl_DC.size(); i++)
    {
      for (int j=0; j< cinfo.dc_huff_tbl_ptrs[0]->bits[i+1]; j++)
        {
          JpegHuffTbl_DC[i].push_back(cinfo.dc_huff_tbl_ptrs[0]->huffval[count]);
          count++;
        }

    }
  count =0;
  for (unsigned int i=0; i < JpegHuffTbl_AC.size(); i++)
    {
      for (int j=0; j< cinfo.ac_huff_tbl_ptrs[0]->bits[i+1]; j++)
        {
          JpegHuffTbl_AC[i].push_back(cinfo.ac_huff_tbl_ptrs[0]->huffval[count]);
          count++;
        }

    }
}

string MTC::ComputeJpegHuff(const Vec<T,cn>& val, std::map<int,string>& huffMap)
{
  return ComputeJpegHuff(val,0,huffMap);
}

string MTC::ComputeJpegHuff(const Vec<T,cn>& val, UINT8 runlength, std::map<int,string>& huffMap)
{
  UINT8 cata =0;
  UINT8 symb  =0;
  int rst=0;
  string rststr="";
  boost::dynamic_bitset<> bset;
  string codes,bcode;
  for (int i=0; i< cn; i++)
    {
      cata =  (UINT8)floor(log(abs(val[i]))/log(2.0));//this is the catagory
      symb = (runlength<<4)+cata+1;
      codes = huffMap[(int)symb]; //get the VLC part
      //get the VLI part
      if(val[i]>0)
        {
          bset = boost::dynamic_bitset<>(cata,(int)val[i] - (2>>cata));
          codes+="1";
        }
      else
        {
          bset = boost::dynamic_bitset<>(cata,(int)val[i] + (2>>(cata+1))-1);
          codes+="0";
        }
      boost::to_string(bset,bcode);
      stat.outputBit+=codes;
      stat.outputBit+=bcode;
      outputfile<<codes<<bcode;
      totalBits+= (cata+codes.length());
      rststr=(codes+bcode);
      //stat.outputBit+=//2's complement
      rst+= (cata+codes.length());
      //search for the catagory
      //for (unsigned int j=0; j< huffTbl.size(); j++)
      //	for (unsigned int k=0; k<huffTbl[j].size();k++)
      //	{
      //		if (huffTbl[j][k] == symb)
      //		{
      //			rst+= (j+1)+(cata+1);//two part, symbol length (j+1) and catalog(=number length)
      //			Huff
      //			j=huffTbl.size();//for jump out loop
      //			break;
      //		}
      //	}
    }
  return rststr;
}
//
//void MTC::SaveLighting(void)
//{
//	fstream lightfile;
//	lightfile.open(path+prefix+baseName+"_TagLighting.txt",ios::out);
//	for (unsigned int i=0; i< this->lightTagVec.size(); i++)
//	{
//		for (unsigned int j=0; j < this->lightTagVec[i].size(); j++)
//		{
//			for (int k=0; k < cn; k++)
//				lightfile<<","<<lightTagVec[i][j][k]<<"\t";
//		}
//		lightfile<<";\n";
//	}
//	lightfile.close();
//	lightfile.open(path+prefix+baseName+"_CanLighting.txt",ios::out);
//	for (unsigned int i=0; i< this->lightCanVec.size(); i++)
//	{
//		for (unsigned int j=0; j < this->lightCanVec[i].size(); j++)
//		{
//			for (int k=0; k < cn; k++)
//				lightfile<<","<<lightCanVec[i][j][k]<<"\t";
//		}
//		lightfile<<";\n";
//	}
//	lightfile.close();
//	Mat samples(lightTagVec.size(),lightTagVec[0].size(), CV_MAKE_TYPE(CV_32F,cn));
//	Mat labels;
//	Mat codebook;
//	for (unsigned int i=0; i< this->lightTagVec.size(); i++)
//		for (unsigned int j=0; j< this->lightTagVec[i].size(); j++)
//		{
//			samples.at<Vec<float,cn>>(i,j) = Vec<float,cn>(lightTagVec[i][j] - lightCanVec[i][j]);
//		}
//	cv::kmeans(samples,64,labels,TermCriteria(CV_TERMCRIT_EPS,100,0.001),10,KMEANS_RANDOM_CENTERS,&codebook);
//
//	//output codebook first
//	string cb=""; //code book string
//	union
//	{
//		float input;   // assumes sizeof(float) == sizeof(int)
//		int   output;
//	}data;
//
//	for (int i=0; i< codebook.rows; i++)
//		for (int j=0; j< codebook.cols; j++)
//			for (int k=0; k< cn; k++)
//			{
//				data.input = codebook.at<Vec<float,cn>>(i,j)[k];
//				std::bitset<sizeof(float) * CHAR_BIT>   bits(data.output);
//				cb+= bits.to_string();
//			}
//	std::cout<<"codebook length "<< cb.length()<<endl;
//
//	for (int i=0; i< labels.rows; i++)
//	{
//		cb+=std::bitset<6>(labels.at<int>(i,0)).to_string();
//	}
//	std::cout<<"label length "<< labels.rows<<endl;
//	cbfile.open(path+"cb_"+baseName+".txt",ios::out);
//	//out put cb
//	cbfile<<cb;
//
//	stat.outputBit = cb + stat.outputBit;
//}

void MTC::InitHuffCoeff(void)
{
  ifstream fh;
  fh.open("DC_hufftalbe.txt");
  std::cout<<"done 4\n";
  int val;
  char value[100];
  string vs;
  while(!fh.eof())
    {
      fh.getline(value,100,' '); //read value
      std::stringstream ss;
      ss << std::hex << value;
      ss >> val;
      fh.getline(value,100,'\n'); //read binary string
      vs= value;
      vs.erase(remove_if(vs.begin(), vs.end(), ::isspace), vs.end());
      MTC::Huffman_DC_Tree.Insert(vs,val);
      Huffman_DC_Map[val]=vs;
      std::cout<<"done 4.25\n";
    }
  fh.close();
  std::cout<<"done 4.5\n";
  fh.open("AC_hufftable.txt");
  while(!fh.eof())
    {
      fh.getline(value,100,' '); //read value
      std::stringstream ss;
      ss << std::hex << value;
      ss >> val;

      fh.getline(value,100,'\n'); //read binary string
      vs = value;
      vs.erase(remove_if(vs.begin(), vs.end(), ::isspace), vs.end());
      MTC::Huffman_AC_Tree.Insert(vs,val);
      Huffman_AC_Map[val]=vs;
    }
  fh.close();
}


///////////// parameter setter //////////////////////////
void MTC::SetInitBlockSize(const Size3& initBlockSize)
{
  this->initSize = initBlockSize;
  this->blockSize = initBlockSize;
}

void MTC::SetSTSIMSubWinSize(const Size3& subWinSize)
{
  //this->stsimSubWinSize = subWinSize;
  //this->SetSubWinSize(subWinSize);
  this->subSize = subWinSize;
}

void MTC::SetSTSIMSubWinStep(const Size3& subWinStep)
{
  //this->stsimSubWinStep = subWinStep;
  //this->ensemble.SetSubWinSize(subWinStep);
  //this->rst.SetSubWinStep(subWinStep);
  //this->SetSubWinStep(subWinStep);
  this->subStep = subWinStep;
}

Size3 MTC::GetSubWinSize() const
{
  return this->subSize;
}

Size3 MTC::GetSubWinStep() const
{
  return this->subStep;
}

void MTC::SetQFactor(double qFactor)
{
  this->qfactor = qFactor;
}

void MTC::SetSTSIMQualigyThrd(double qThrd)
{
  this->qualityThrd = qThrd;
  this->orgQualityThrd = qThrd;
}

void MTC::SetMSEThrd(double mseThrd)
{
  this->mseThrd = mseThrd;
}

void MTC::SetInitQSize(double qSize)
{
  this->initQSize = qSize;
}

void MTC::SetCodingMode(CodingMode qmode)
{
  this->mode = qmode;
}
void MTC::SetMetricModifier(MetricModifier modifier)
{
  this->metricModifier = modifier;
}
void MTC::SetTest(bool test)
{
  this->unitTest = test;
}
void MTC::SetMatchingMethod(MatchingMethod m_method)
{
  this->matching_method = m_method;
}

void MTC::SetSTSIM2PoolType(FeaturePoolType pool_type)
{
  this->stsim2PoolType = pool_type;
}
void MTC::SetBlendingMethod(BlendingMethod blend_method)
{
  this->blend_method = blend_method;
}
void MTC::SetPQILFBlockSize(int pqiSize)
{
  this->nPQILFSize = pqiSize;
}

void MTC::SetOverlapSizeByRatio(const Size3_<double>& ratio)
{
  this->overlapSize = Size3(Size3_<double>(blockSize)*ratio);
}

void MTC::SetOverlapSize(const Size3& overlap)
{
  this->overlapSize = overlap;
}

void MTC::SetCandidNum(unsigned int candidNum)
{
  this->candidNum = candidNum;
}

void MTC::SetPQIRectType(bool rectType)
{
  this->rectType = rectType;
}

void MTC::SetPBType(PostBlendingType pbType)
{
  this->postBlendType = pbType;
}

void MTC::SetLCType(LightingCorrectionType lcType)
{
  this->lightCorrectionType = lcType;
}

void MTC::SetJpegQFactor(int q)
{
  this->jpegQNum = q;
}

void MTC::SetJpegQTblType(JPEGType quanType)
{
  this->jpegType = quanType;
}

void MTC::UpdateParameters(void)
{
  //SQLHANDLE hdlEnv,hdlConn,hdlStmt;
  ////SQLCHAR* stmt;// = (SQLCHAR*)"SELECT * FROM Table5";
  //string stmt;
  //SQLCHAR retVal[256];
  //SQLINTEGER cbData;
  //RETCODE retcode = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &hdlEnv);
  //retcode = SQLSetEnvAttr(hdlEnv,SQL_ATTR_ODBC_VERSION,(void*)SQL_OV_ODBC3,0);
  //retcode = SQLAllocHandle(SQL_HANDLE_DBC, hdlEnv, &hdlConn);
  //retcode = SQLConnect(hdlConn, (SQLCHAR*)"MTC",SQL_NTS,(SQLCHAR*)"Guoxin Jin",SQL_NTS, (SQLCHAR*)"ti+PIsR5", SQL_NTS);
  //retcode = SQLAllocHandle(SQL_HANDLE_STMT, hdlConn, &hdlStmt);
  //retcode = SQLExecDirect(hdlStmt, (SQLCHAR*)"INSERT INTO Table5(A,C) VALUES('TJ','mm')", SQL_NTS);
  //if ( (retcode != SQL_SUCCESS) && (retcode != SQL_NEED_DATA) && (retcode != SQL_SUCCESS_WITH_INFO) ) {
  //  printf("SQLExecDirect Failed\n\n");
  //}
  //SQLExecDirect(hdlStmt, stmt, SQL_NTS);


  //while(SQL_SUCCEEDED(SQLFetch(hdlStmt)))
  //{
  //	retcode = SQLGetData(hdlStmt,1,SQL_C_CHAR,retVal,255,&cbData);
  //	std::cout << retVal << "\t";
  //	retcode = SQLGetData(hdlStmt,2,SQL_C_CHAR,retVal,255,&cbData);
  //	std::cout << retVal << std::endl;
  //	//SQLGetData(hdlStmt,3,SQL_C_CHAR,retVal,255,&cbData);
  //	//std::cout << retVal << std::endl;


  //}
#if defined(_USE_DATABASE)
  try
  {
    conn.Connect(DB_HOST,DB_USR,DB_PWD,SA_SQLServer_Client);
    cmd.setConnection(&conn);
    SAString dbstr("insert into ");
    dbstr += DB_NAME;
    dbstr += "([Begin Time]) values (GETDATE())";
    cmd.setCommandText(dbstr);
    cmd.Execute();
    std::cout<<green<<endl<<"DB Connected"<<white<<endl;
  }
  catch (SAException &x)
  {
    try
    {
      conn.Rollback();
    }
    catch(SAException &)
    {

    }
    printf("%s\n",(const char*)x.ErrText());
  }
#endif
  levels = (int)ceil(log10(double(initSize.height))/log10(double(2)))+2;
  qstep.push_back(this->initQSize);
  thrd.push_back(this->initQSize);
  for (int i=1; i<levels-2; i++)
    {
      thrd.push_back(qfactor*thrd[i-1]);
      qstep.push_back(qfactor*qstep[i-1]);
    }
  thrd.push_back(*(thrd.end()-1));
  qstep.push_back(*(thrd.end()-1));
  thrd.push_back(*(thrd.end()-1));
  qstep.push_back(*(thrd.end()-1));
  stat = Statistics(levels,initSize.height);

  if (mode == CodingMode::CODING_MODE_JPEG || mode == CodingMode::CODING_MODE_MTC)
    {
      if (ensemble.channels()!=1)
        CV_Error(CV_StsUnsupportedFormat,"dct can  only work on 1 channel image");

      InitJpegParam(jpegQNum);
    }
  double tempVar = this->varThrd1 ;
  std::cout<<"done 3.25\n";
  ReInitGrid();
  std::cout<<"done 3.375\n";
  this->SetVarThrd(tempVar);
  std::cout<<"done 3.45\n";
  InitHuffCoeff();
  std::cout<<"done 3.5\n";
  footmap = Tensor<T,1>(ensemble.size(),Vec<T,1>::all(255));
  dctEnsemble = Tensor<T,cn>(ensemble.size());
  rst_with_seam = rst.CvtColor<3>(COLOR_GRAY2RGB);
  pred_rst = rst.Clone();//.CvtColor<3>(CV_GRAY2RGB);
  plc_rst = rst.Clone();
  pred_afterLC = rst.Clone();
  pred_postLC = rst.Clone();
  //update output path
  std::cout<<"done 3.6\n";
  time_t t = time(&stat.beginT);
  tm* tt = localtime(&t);
  ///data base update from here
  //retcode = SQLExecDirect(hdlStmt, (SQLCHAR*)"INSERT INTO MTC_Table(Date) VALUES()", SQL_NTS);
  //stmt = "INSERT INTO MTC_Table(Date) VALUES("+ boost::lexical_cast<string>(tt->tm_mon+1) +"-"+ boost::lexical_cast<string>(tt->tm_mday) + "-" + boost::lexical_cast<string>(tt->tm_year+1900) + ")";
  //retcode = SQLExecDirect(hdlStmt,(SQLCHAR*)stmt.c_str(),SQL_NTS);
  //SQLExecDirect(hdlStmt,(SQLCHAR*)"SELECT Date FROM MTC_Table",SQL_NTS);
  //while(SQL_SUCCEEDED(SQLFetch(hdlStmt)))
  //	{
  //	retcode = SQLGetData(hdlStmt,1,SQL_C_CHAR,retVal,255,&cbData);
  //	std::cout << retVal << endl;
  //}

  stat.beginH = tt->tm_hour;
  stat.beginM = tt->tm_min;
  stat.beginS = tt->tm_sec;
  stat.tpBlockNum.clear();
  stat.tpBits.clear();
  std::cout<<"done 3.6\n";
  for (int l=initSize.height; l > 8; l/=2)
    {
      stat.tpBlockNum.push_back(0);
      stat.tpBits.push_back(0);
    }
  path = "./temp/";
  path+=boost::lexical_cast<string>(tt->tm_year+1900);//year start with 1900
  path+=boost::lexical_cast<string>(month[tt->tm_mon]);
  path+=boost::lexical_cast<string>(tt->tm_mday);
  path += "/";
  //if (postBlendType != POST_BLENDING_ONLINE && lightCorrectionType != NO_LIGHTING_CORRECTION||lightCorrectionType!=POISSON_LC)
  //	CV_Error(CV_StsError,"cannot do lighting correction without online post blending.");
  std::cout<<"done 3.7\n";
  mkdir(path.c_str(),0755);
  std::cout<<"done 3.75\n";
  string PID = boost::lexical_cast<string>(getpid());
  candPosLog.open(path+prefix+PID+"_TarCandPos.txt",ios::out);
  ifstream ifile(path+"log.txt");
  if (ifile) {
      ifile>>count_exp;
      cout<<"experiment count: "<<count_exp<<endl;
      ifile.close();
      count_exp=boost::lexical_cast<string>(boost::lexical_cast<int>(count_exp)+1);
      logfile.open(path+"log.txt",ios::out|ios::in);
      logfile.seekg(0,ios::beg);
      logfile.write(count_exp.c_str(),count_exp.length());
      logfile.close();
      logfile.open(path+"log.txt",ios::out|ios::app);
      logfile.seekg(0,ios::end);
      logfile<<endl<<"No."<<count_exp<<endl;
      logfile<<"PID: "<<getpid()<<endl;
      //logfile.close();
      // The file exists, and is open for input
      std::cout<<"done 3.76\n";
    }
  else
    {
      logfile.open(path+"log.txt",ios::out);
      logfile<<"1"<<endl;
      this->count_exp ="1";
      logfile<<endl<<"No."<<count_exp<<endl;
      logfile<<"PID: "<<getpid()<<endl;
      //logfile.close();
      std::cout<<"done 3.77\n";
    }
  std::cout<<"done 3.8\n";
  if (mode == CodingMode::CODING_MODE_PQI)
    prefix = "PQI_";
  else if (mode == CodingMode::CODING_MODE_POST_TP)
    prefix = "Post_";
  else if (mode == CodingMode::CODING_MODE_JPEG)
    prefix = "J_";
  else if (mode == CodingMode::CODING_MODE_MTC)
    prefix = "TJ_";
  else
    prefix = "MTC_";
  if (lightCorrectionType==LightingCorrectionType::HAS_LIGHTING_CORRECTION)
    prefix = prefix + "light_";
  else
    prefix = prefix + "nolgt_";
  if (postBlendType == PostBlendingType::NO_POST_BLENDING)
    prefix+= "noPB_";
  else if (postBlendType == PostBlendingType::POST_BLENDING_ONLINE)
    prefix+= "PB_";
  else
    prefix+= "offlinePB_";
  //add prefix indicate for min pooling stsim and new stsim boundary
  prefix +="min_";
  prefix +="nb_";
  //add prefix indicate the qfactor
  prefix += boost::lexical_cast<string>(this->qfactor).substr(0,4);//ratio of stsim threshold 16/32
  prefix += "_";

  string::size_type idx = cFileName.find_last_of(".");
  string::size_type idx2 = cFileName.find_last_of("/");
  baseName = cFileName.substr(0,idx).substr(idx2+1,idx-idx2);
  //edgeEnsemble.Load(baseName+"_aca.tiff");
  //edgeEnsemble.Display(1);
  prefix +=boost::lexical_cast<string>(this->qualityThrd).substr(0,5);
  //prefix += "pure_";
  prefix+="_";
  if (mode==CodingMode::CODING_MODE_MTC ||mode== CodingMode::CODING_MODE_JPEG)
    prefix += boost::lexical_cast<string>(this->jpegQNum);
  else
    prefix += boost::lexical_cast<string>(this->initQSize*this->qfactor).substr(0,4);
  prefix += "_";
  prefix += boost::lexical_cast<string>(this->initSize.height);
  prefix += "_";
  prefix += boost::lexical_cast<string>(this->candidNum);
  prefix += "_";
  prefix += boost::lexical_cast<string>(this->searchStep.height);
  prefix += "_";
  totalBits=0;

  //dctEnsemble.Print();
  //reading edge map
  //if use oracle, disable temp == aug 17 2012
  edgeMap = vector<vector<int>>(this->gridSize.height);
  //ifstream fh;
  //fh.open(baseName+"_oracle.csv");
  //string vs,vss;
  //int lineIdx = 0;
  //while(!fh.eof())
  //{
  //	fh>>vs; //read value
  //	while (vs.length()>0)
  //	{
  //		int comaPos = vs.find_first_of(",");
  //		if (comaPos<0)
  //		{
  //			vss = vs;
  //			vs="";
  //		}
  //		else
  //		{
  //			vss = vs.substr(0,comaPos);
  //			vs = vs.substr(comaPos+1,vs.length()-comaPos);
  //		}
  //		if (vss == "T")
  //			edgeMap[lineIdx].push_back(2);
  //		else if (vss == "S")
  //			edgeMap[lineIdx].push_back(4);
  //		else if (vss == "T/T")
  //			edgeMap[lineIdx].push_back(0);
  //		else if (vss == "T/S")
  //			edgeMap[lineIdx].push_back(3);
  //		else if (vss == "S/S")
  //			edgeMap[lineIdx].push_back(1);
  //		else
  //			edgeMap[lineIdx].push_back(5);
  //	}
  //	lineIdx++;
  //}
  //fh.close();
  fstream matchingfile;
  matchingfile.open("./temp/matching.txt",ios::out);
  matchingfile<<"starting recording matching .....\n";
  matchingfile.close();
#ifdef RECORD_EVERYTHING
  system("rm -rf ./everything");
  //std::remove("./everything");
  mkdir("./everything",0755);
  ofstream everything;
  everything.open("./everything/everything.txt",ios::out);
  everything.close();
#endif
  if (metricModifier == MetricModifier::STSIM3_LSE)
    {
      ifstream file ( "../Metric/weight.txt" );
      vector<string> value;

      while ( file.good() )
        {
          string temp;
          getline ( file, temp, ',' );
          if (!temp.empty())
            value.push_back(temp);
        }
      file.close();
      Tensor<double,1> dummy;
      lse_weight = Mat(value.size(),1,CV_64F);
      for (unsigned int i=0; i< value.size(); i++)
        {
          lse_weight.at<double>(i,0) = std::stod(value[i]);
        }

    }
}

void MTC::CollectLighing(void)
{
  for (int i=0; i < gridSize.height; i++)
    for (int j=0; j < gridSize.width; j++)
      {
        QTree<T,cn>& temp = GetNode(Point3i(i,j,0));
        QTree<T,cn> qNode(QNode<T,cn>(ensemble.Clone(),temp.size(),temp.offset(),temp.overlap()));
        //	qNode.Print();
        if(qNode.offset().x>0 && qNode.offset().y>0)
          {
            lightPlanes.push_back(qNode.LightingCorrection(Tensor<T,cn>(qNode.size())));
            lightCanVec.push_back(qNode.lt.GetCanLighting());

          }
        //then split
        qNode.Split();
        QTree<T,cn> *tempTree = qNode.NextLeaf();
        while(true)
          {
            if(tempTree == NULL)
              break;
            else
              {
                if(tempTree->offset().x > 0 && tempTree->offset().y>0)
                  {
                    QNode<T,cn> tempNode(this->ensemble.Clone(),tempTree->size(),tempTree->offset(),qNode.overlap()/2);
                    lightPlanes.push_back(tempNode.LightingCorrection(Tensor<T,cn>(tempTree->size())));
                    lightCanVec.push_back(tempNode.lt.GetCanLighting());
                  }
              }
            if(tempTree->GetTreePos()==tempTree->GetPeerTreeSize()-1)
              {
                break;
              }
            else
              tempTree=tempTree->NextLeaf();
          }
      }

  //build VQ

  Mat samples(lightCanVec.size(),lightCanVec[0].size(), CV_MAKE_TYPE(DataType<T>::depth,cn));
  Mat labels;
  Mat codebook;
  for (unsigned int i=0; i< this->lightCanVec.size(); i++)
    for (unsigned int j=0; j< this->lightCanVec[i].size(); j++)
      {
        samples.at<Vec<T,cn>>(i,j) = Vec<float,cn>(lightCanVec[i][j]);
      }
  Tensor<T,cn> ss(samples);
  //	mylib::DisplayMat(samples);
  BuildVQCodebook(lightPlanes,ss,6);
  //cv::kmeans(samples,64,labels,TermCriteria(CV_TERMCRIT_EPS,100,0.001),10,KMEANS_RANDOM_CENTERS,&codebook);
  ////mylib::DisplayMat(codebook);
  //VQCodeBook = Tensor<T,cn>(Tensor<float,cn>(codebook));
  ////output codebook bits
  //string cb=""; //code book string
  //union
  //{
  //	float input;   // assumes sizeof(float) == sizeof(int)
  //	int   output;
  //}data;
  //for (int i=0; i< codebook.rows; i++)
  //	for (int j=0; j< codebook.cols; j++)
  //		for (int k=0; k< cn; k++)
  //		{
  //			data.input = codebook.at<Vec<float,cn>>(i,j)[k];
  //			std::bitset<sizeof(float) * CHAR_BIT>   bits(data.output);
  //			cb+= bits.to_string();
  //		}
  //std::cout<<"codebook length "<< cb.length()<<endl;
  //totalBits+=cb.length();
  //outputfile<<cb;
  lightCanVec.clear();//flush the data
  //VQCodeBook.Print();
}
//
//int MTC::SearchCodeword(vector<Vec<T,cn>>& poly)
//{
//	int label=-1;
//	double dist=INFINITE;
//	double temp_dist = 0;
//	Vec<T,cn> dist_vec;
//	if(VQCodeBook.size().width!= poly.size())
//	{
//		CV_Error(CV_StsUnmatchedSizes,"the input polynomial has different size as codebook dim");
//	}
//	for (int i=0; i< VQCodeBook.size().height; i++)
//	{
//		dist_vec =Vec<T,cn>::all(0);
//		temp_dist = 0;
//		for(int j=0; j< VQCodeBook.size().width; j++)
//			dist_vec+= mylib::VecPow(VQCodeBook(i,j,0) - poly[j],2.0);
//		for (int c=0; c<cn; c++)
//			temp_dist+= dist_vec[c];
//		if (temp_dist <= dist)
//		{
//			dist = temp_dist;
//			label = i;
//		}
//	}
//	return label;
//}

void MTC::BuildVQCodebook(vector<Tensor<T,cn>>& lightPlanes, Tensor<T,cn>& coeffs, int bits)
{
  cv::RNG rng(time(0));
  //init guess
  VQCodeBook = Tensor<T,cn>(Size3((1<<bits),coeffs.size().width,coeffs.size().depth));
  //vector<Tensor<T,cn>*> revCenterLight_large(1<<bits,Tensor<T,cn>(initSize));
  //vector<Tensor<T,cn>*> revCenterLight_small(1<<bits,Tensor<T,cn>(initSize/Size3(2,2,1)));
  VQLabels = vector<int>(lightPlanes.size(),-1);
  double totalerror=0;
  vector<int> bins(1<<bits);
  Tensor<T,cn> tempCoeff(Size3(1,VQCodeBook.size().width,1));
  //init
  for (int j=0; j< (1<<bits); j++)
    {

      Tensor<T,cn> temp = coeffs.Row(rng.uniform(0,coeffs.size().height));
      VQCodeBook.SetBlock(Point3i(j,0,0),temp);
      //revCenterLight_large[j] = revCenterLight_large[j].BuildLightingPlane(*temp);
      //revCenterLight_small[j] = revCenterLight_small[j].BuildLightingPlane(*temp);
      //delete temp;
    }
  //VQCodeBook.Print();
  for (int loop = 0; loop < 90; loop++)
    {
      totalerror = 0;
      std::cout<<"loop:"<<loop<<endl;
      for (int i=0; i<coeffs.size().height; i++)
        {
          //lightPlanes[i].Print();
          int label = -1;
          double dist = std::numeric_limits<double>::max();
          double temp_dist =0;
          for (int k=0; k< VQCodeBook.size().height; k++)
            {

              Tensor<T,cn> temp = VQCodeBook.Row(k);
              QNode<T,cn> light(lightPlanes[i].size());
              Tensor<T,cn> plane = light.lt.BuildLightingPlane(light, temp.Transpose());
              temp_dist  = plane.AbsDiff(lightPlanes[i]).Pow(2.0).Mean()[0];
              if (temp_dist < dist)
                {
                  dist = temp_dist;
                  VQLabels[i]=k;
                }

            }
          totalerror +=dist;
          //	sampleLT.push_back(lightPlane);
        }
      //flush codebook
      //VQCodeBook = Tensor<T,cn>(Size3((1<<bits),coeffs.size().width,coeffs.size().depth));
      bins = vector<int>(1<<bits);
      for (unsigned int i=0; i< VQLabels.size(); i++)
        {
          Tensor<T,cn> tempCBRow  = VQCodeBook.Row(VQLabels[i]);
          Tensor<T,cn> tempCoeffRow = coeffs.Row(VQLabels[i]);
          if(bins[VQLabels[i]] == 0)
            {
              (tempCBRow) = (tempCoeffRow);
            }
          else
            {
              tempCBRow = (tempCBRow) + (tempCoeffRow);
            }
          VQCodeBook.SetBlock(Point3i(VQLabels[i],0,0),tempCBRow);
          bins[VQLabels[i]]++;
        }
      for (unsigned int j=0; j<bins.size(); j++)
        {
          if (bins[j]>0)
            {
              Tensor<T,cn> tempCBRow  = VQCodeBook.Row(j);
              (tempCBRow) = (tempCBRow)/bins[j];
              VQCodeBook.SetBlock(Point3i(j,0,0),tempCBRow);
            }

        }
      //VQCodeBook.Print();
      //std::cout<<totalerror<<endl;
    }
}

//some backup codes
//if (jpegType == JPEG_WASTON)//convert back to gray
//			{
//				vector<Mat> planes;
//				Tensor<float,3> colorVerNode(floatVerNode.size());
//
//				planes.push_back(floatVerNode[0]);
//				planes.push_back(cv::Mat::zeros(floatVerNode.size().height,floatVerNode.size().width,floatVerNode.type()));
//				planes.push_back(cv::Mat::zeros(floatVerNode.size().height,floatVerNode.size().width,floatVerNode.type()));
//				cv::merge(planes,colorVerNode[0]);
//				colorVerNode = colorVerNode.CvtColor<3>(CV_YCrCb2RGB);
//			//	colorVerNode.Print();
//				floatVerNode = colorVerNode.CvtColor<1>(CV_RGB2GRAY);
//
//			}


void MTC::ParseCfg(const string& cfgname)
{
  ifstream cfgstream(cfgname, ifstream::in);//read configure file
  if (!cfgstream)
    {
      cerr << "Failed to open config file: `" << cfgname << "'" << endl;
      exit(EXIT_FAILURE);
    }
  ScanFile(cfgstream);
}
void MTC::ScanFile(istream& in)
{
  do
    {
      string line;
      getline(in, line);
      ScanLine(line);
    }
  while(!!in);
}

void MTC::ScanLine(string& line)
{
  /* strip any leading whitespace */
  size_t start = line.find_first_not_of(" \t\n\r");
  if (start == string::npos)
    {
      /* blank line */
      return;
    }
  if (line[start] == '#')
    {
      /* comment line */
      return;
    }
  /* look for first whitespace or ':' after the option end */
  size_t option_end = line.find_first_of(": \t\n\r",start);
  string option = line.substr(start, option_end - start);

  /* look for ':', eat up any whitespace first */
  start = line.find_first_not_of(" \t\n\r", option_end);
  if (start == string::npos)
    {
      /* error: badly formatted line */
      return;
    }
  if (line[start] != ':')
    {
      /* error: badly formatted line */
      return;
    }

  /* look for start of value string -- eat up any leading whitespace */
  start = line.find_first_not_of(" \t\n\r", ++start);
  if (start == string::npos)
    {
      /* error: badly formatted line */
      return;
    }

  /* extract the value part, which may contain embedded spaces
       * by searching for a word at a time, until we hit a comment or end of line */
  size_t value_end = start;
  do
    {
      if (line[value_end] == '#')
        {
          /* rest of line is a comment */
          value_end--;
          break;
        }
      value_end = line.find_first_of(" \t\n\r", value_end);
      /* consume any white space, incase there is another word.
         * any trailing whitespace will be removed shortly */
      value_end = line.find_first_not_of(" \t\n\r", value_end);
    }
  while (value_end != string::npos);
  /* strip any trailing space from value*/
  value_end = line.find_last_not_of(" \t\n\r", value_end);

  string value;
  if (value_end >= start)
    {
      value = line.substr(start, value_end +1 - start);
    }
  else
    {
      /* error: no value */
      return;
    }

  start = value.find_first_of("(");
  value_end = value.find_first_of(")");
  if (start == string::npos)// fill in single augument parameters
    {
      if (option.compare("ImageName")==0)
        {
          this->Init(value);
        }
      else if (option.compare("CodeMode")==0)
        {
          if (value == "MTC")
            this->SetCodingMode(CodingMode::CODING_MODE_MTC);
          else if (value=="JPEG")
            this->SetCodingMode(CodingMode::CODING_MODE_JPEG)	;
          else if (value=="PQI")
            this->SetCodingMode(CodingMode::CODING_MODE_PQI);
          else if (value=="TPQI")
            this->SetCodingMode(CodingMode::CODING_MODE_TPQI);
          else if (value=="PTP")
            this->SetCodingMode(CodingMode::CODING_MODE_POST_TP);
          else if (value=="TPSS")
            this->SetCodingMode(CodingMode::CODING_MODE_TPSS);
          else
            CV_Error(CV_StsBadFlag,"Wrong coding mode");
        }
      else if (option.compare("QP")==0)
        this->SetJpegQFactor(boost::lexical_cast<int>(value));
      else if (option.compare("QualityThreshold")==0)
        this->SetSTSIMQualigyThrd(boost::lexical_cast<double>(value));
      else if (option.compare("SideMatchingThreshold")==0)
        this->SetMSEThrd(boost::lexical_cast<double>(value));
      else if (option.compare("QualityElevationFactor")==0)
        this->SetQFactor(boost::lexical_cast<double>(value));
      else if (option.compare("CandidateNumber")==0)
        this->SetCandidNum(boost::lexical_cast<int>(value));
      else if (option.compare("PreVarianceThreshold")==0)
        this->SetVarThrd(boost::lexical_cast<double>(value));
      else if (option.compare("InitQSize")==0)
        this->SetInitQSize(boost::lexical_cast<int>(value));
      else if (option.compare("PQIRectType")==0)
        this->SetPQIRectType(boost::lexical_cast<bool>(value));
      else if (option.compare("PQILFBlockSize")==0)
        this->SetPQILFBlockSize(boost::lexical_cast<int>(value));
      else if (option.compare("MatchingMethod")==0)
        {
          MatchingMethod temp=MatchingMethod(0);
          if (value =="MATCHING_MSE")
            temp = MatchingMethod::MATCHING_MSE;
          else if (value =="MATCHING_SAT")
            temp = MatchingMethod::MATCHING_SAT;
          else if (value =="MATCHING_VAR")
            temp = MatchingMethod::MATCHING_VAR;
          else if (value =="MATCHING_SAD")
            temp = MatchingMethod::MATCHING_SAD;
          else if (value=="MATCHING_MSE_CONSTRAINT")
            temp = MatchingMethod::MATCHING_MSE_CONSTRAINT;
          else if (value=="MATCHING_HIERARCHY")
            temp = MatchingMethod::MATCHING_HIERARCHY;
          else if (value=="MATCHING_OPENCV")
            temp = MatchingMethod::MATCHING_OPENCV;
          else if (value=="MATCHING_DIRECT")
            temp = MatchingMethod::MATCHING_DIRECT;
          else
            CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);
          this->SetMatchingMethod(temp);
        }
      else if (option.compare("LightingCorrection")==0)
        {
          LightingCorrectionType temp=LightingCorrectionType(0);
          if (value=="NO_LIGHTING_CORRECTION")
            temp = LightingCorrectionType::NO_LIGHTING_CORRECTION;
          else if (value =="HAS_LIGHTING_CORRECTION")
            temp=LightingCorrectionType::HAS_LIGHTING_CORRECTION;
          else if (value =="PQI_LF_ENCODING")
            temp=LightingCorrectionType::PQI_LF_ENCODING;
          else if (value =="PREDEF_LIGHTING")
            temp=LightingCorrectionType::PREDEF_LIGHTING;
          else if (value =="POISSON_LC")
            temp=LightingCorrectionType::POISSON_LC;
          else
            CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);

          this->SetLCType(temp);
        }
      else if (option.compare("PostBLending")==0)
        {
          PostBlendingType temp=PostBlendingType(0);
          if (value=="NO_POST_BLENDING")
            temp=PostBlendingType::NO_POST_BLENDING;
          else if (value=="POST_BLENDING_ONLINE")
            temp=PostBlendingType::POST_BLENDING_ONLINE;
          else if (value=="POST_BLENDING_OFFLINE")
            temp=PostBlendingType::POST_BLENDING_OFFLINE;
          else
            CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);
          this->SetPBType(temp);
        }
      else if (option.compare("JPEG_Modification")==0)
        {
          JPEGType temp=JPEGType(0);
          if (value=="JPEG_BASELINE")
            temp=JPEGType::JPEG_BASELINE;
          else if (value=="JPEG_WASTON")
            temp=JPEGType::JPEG_WASTON;
          else if (value=="JPEG_ADAPTIVE")
            temp =JPEGType::JPEG_ADAPTIVE;
          else if (value=="JPEG_LF_3")
            temp = JPEGType::JPEG_LF_3;
          else if (value=="JPEG_LF_6")
            temp =JPEGType::JPEG_LF_6;
          else	if (value=="JPEG_CUSTOM_TBL")
            temp =JPEGType::JPEG_CUSTOM_TBL;
          else
            CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);

          this->SetJpegQTblType(temp);
        }
      else if (option.compare("Metric_Type")==0)
        {
          MetricModifier temp=MetricModifier(0);
          if (value=="STSIM2_BASELINE")
            temp=MetricModifier::STSIM2_BASELINE;
          else if (value=="PQI_METRIC")
            temp=MetricModifier::PQI_METRIC;
          else if (value=="STSIM2_ADT_NAIVE")
            temp=MetricModifier::STSIM2_ADT_NAIVE;
          else if (value=="STSIM2_NEW_L1")
            temp=MetricModifier::STSIM2_NEW_L1;
          else if (value=="STSIM2_NEW_L2")
            temp=MetricModifier::STSIM2_NEW_L2;
          else if (value=="STSIM2_NEW_L3")
            temp=MetricModifier::STSIM2_NEW_L3;
          else if (value=="STSIM2_NEW_L4")
            temp=MetricModifier::STSIM2_NEW_L4;
          else if (value=="MSE_BASELINE")
            temp = MetricModifier::MSE_BASELINE;
          else if (value=="M_DIST")
            temp=MetricModifier::MAHALANOBIS_DIST;
          else if (value=="LRI_METRIC")
            temp=MetricModifier::LRI_METRIC;
          else if (value=="SE_MSE")
            temp=MetricModifier::SE_MSE;
          else if (value=="STSIM2_SE_MSE")
            temp=MetricModifier::STSIM2_SE_MSE;
          else if (value=="STSIM2_TUNE")
            temp=MetricModifier::STSIM2_TUNE;
          else if (value=="SVM_METRIC")
            temp=MetricModifier::SVM_METRIC;
          else if (value=="STSIM3_LSE")
            temp=MetricModifier::STSIM3_LSE;
          else
            CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);
          this->SetMetricModifier(temp);
        }
      else if (option.compare("STSIM_Pooling")==0)
        {
          FeaturePoolType temp = FeaturePoolType(0);
          if (value=="STSIM2_POOL_AVE")
            temp=FeaturePoolType::FEATURE_POOL_AVE;
          else if (value=="STSIM2_POOL_MIN")
            temp=FeaturePoolType::FEATURE_POOL_MIN;
          else
            CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);

          this->SetSTSIM2PoolType(temp);
        }
      else if (option.compare("BLENDING_METHOD")==0)
        {
          BlendingMethod temp = BlendingMethod(0);
          if (value=="SHORTEST_PATH_BLENDING")
            temp = BlendingMethod::SHORTEST_PATH_BLENDING;
          else if (value=="GRADIENT_BLENDING")
            temp = BlendingMethod::GRADIENT_BLENDING;
          else if (value=="NO_BLENDING")
            temp = BlendingMethod::NO_BLENDING;
          else
            CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);
          this->SetBlendingMethod(temp);
        }
      else if (option.compare("FootRegionSelect")==0)
        {
          this->SetFootComputeRegion(boost::lexical_cast<int>(value));
        }
      else if (option.compare("FootComputeMethod")==0)
        {
          this->SetFootComputeMethod(boost::lexical_cast<int>(value));
        }
      /* else if (option.compare("DEBUG_X")==0)
      {
         DEBUG_X = boost::lexical_cast<int>(value);
      }
      else if (option.compare("DEBUG_Y")==0)
      {
         DEBUG_Y = boost::lexical_cast<int>(value);
      }
      else if (option.compare("DEBUG_DISP_X")==0)
      {
         DEBUG_DISP_X = boost::lexical_cast<int>(value);
      }
      else if (option.compare("DEBUG_DISP_Y")==0)
      {
         DEBUG_DISP_Y = boost::lexical_cast<int>(value);
      }
      else if (option.compare("DEBUG_SIZE")==0)
      {
         DEBUG_SIZE = boost::lexical_cast<int>(value);
      }*/
      else
        {
          CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);
        }
    }
  else
    {
      vector<int> vals;
      string tempstr="";
      unsigned int idx=start+1;
      do
        {
          if (value.compare(idx,1,",")==0||value.compare(idx,1,")")==0)
            {
              vals.push_back(boost::lexical_cast<int>(tempstr));
              tempstr="";
            }
          else
            {
              tempstr+=(value[idx]);
            }
          idx++;
        }
      while ( idx < value.size()) ;

      if (option.compare("InitBlockSize")==0)
        {
          this->SetInitBlockSize(Size3(vals[0],vals[1],vals[2]));
          this->InitExt();
        }
      else if (option.compare("STSIMSubWinSize")==0)
        {
          this->SetSTSIMSubWinSize(Size3(vals[0],vals[1],vals[2]));
        }
      else if (option.compare("STSIMSubWinStep")==0)
        {
          this->SetSTSIMSubWinStep(Size3(vals[0],vals[1],vals[2]));
        }
      else if (option.compare("OverlapSize")==0)
        {
          this->SetOverlapSize(Size3(vals[0],vals[1],vals[2]));
        }
      else if (option.compare("SearchStep")==0)
        {
          this->SetSearchStep(Size3(vals[0],vals[1],vals[2]));
        }
      else
        {
          CV_Error(CV_StsBadFlag,"wrong parameter"+option+":"+value);
        }
    }

  /* store the value in option */
  //storePair(opts, true, false, option, value);
}

void MTC::SetFootComputeRegion(int region)
{
  this->footComputeRegion= region;
}
void MTC::SetFootComputeMethod(int method)
{
  this->footComputeMethod = method;
}

CompareCriteria MTC::ParseCriteria(MetricModifier metricModifer)
{
  if(metricModifer == MetricModifier::STSIM2_BASELINE ||
     metricModifer == MetricModifier::SSIM2_DCT||
     metricModifer == MetricModifier::STSIM2_ADT_NAIVE||
     metricModifer == MetricModifier::STSIM2_NEW_L1||
     metricModifer == MetricModifier::STSIM2_NEW_L2 ||
     metricModifer == MetricModifier::STSIM2_NEW_L3 ||
     metricModifer == MetricModifier::STSIM2_NEW_L4 ||
     metricModifer == MetricModifier::STSIM2_TUNE ||
     metricModifer == MetricModifier::STSIM3_LSE
     )
    return CompareCriteria::COMPARE_CRITERIA_SSIM;
  else if(metricModifer == MetricModifier::PQI_METRIC)
    return CompareCriteria::COMPARE_CRITERIA_INTERP;
  else if(metricModifer == MetricModifier::LRI_METRIC)
    return CompareCriteria::COMPARE_CRITERIA_LRI;
  else if(metricModifer == MetricModifier::SE_MSE||metricModifer==MetricModifier::MSE_BASELINE||metricModifer == MetricModifier::STSIM2_SE_MSE)
    return CompareCriteria::COMPARE_CRITERIA_MSE;
  else if(metricModifer == MetricModifier::MAHALANOBIS_DIST)
    return CompareCriteria::COMPARE_CRITERIA_MAHALANOBIS;
  else if (metricModifer == MetricModifier::SAD_BASELINE)
    return CompareCriteria::COMPARE_CRITERIA_SAD;
  else if (metricModifer ==MetricModifier::SVM_METRIC)
    return CompareCriteria::COMPARE_CRITERIA_SVM;
  else
    return CompareCriteria::COMPARE_CRITERIA_OTHER;

}

void MTC::InitExt()
{
  ensembleExt = this->ensemble.ExtendBoundary(this->initSize/2);
  for (int z=0; z< this->ensemble.size().depth; z++)
    {
      Mat temp;
      cv::copyMakeBorder(ensemble.GetFrame(z),temp,this->initSize.height/2,this->initSize.height/2, this->initSize.width/2, this->initSize.width/2, BORDER_REFLECT);
      ensembleExt.SetFrame(z,temp);
    }
  return;
}
}

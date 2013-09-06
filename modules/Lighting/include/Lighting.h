#ifndef LIGHTING_H
#define LIGHTING_H
#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
#include <TensorLite.h>
#include <ThinPlateSpline.h>
using namespace tensor;
namespace lighting {

#ifndef LIGHTING_ENCODE_TABLE
#define LIGHTING_ENCODE_TABLE
static const float base8[] = {16,11,12,14,12,10,10000}; //(0,0),(0,1),(1,0),(2,0),(2,2),(0,2),others
static const float base16[]= {16,14,14,12,13,11,10000};
static const float base32[]= {16,15,15,14,14,14,10000};

static const int light_ac_length[] = {2,2,2,3,4,5,6,7,8,9,10,11,12};//cata:0,1,2,....,14
static const int light_dc_length[] = {6,5,4,3,2,1,7,8,9,10,11,12};//no runlength, cata 0 is the EOB which should not be used, cata: 0,1,2,3,....,14

#endif


class Lighting
{
public:
  Lighting();
  //! return mat  has double data type
  EXPORTLIB Mat LSFitting(const Mat& im, int order = 2) const;

  EXPORTLIB Mat BuildLightingPlane(const Mat& im, const Mat& param, int order = 2) const;
  EXPORTLIB Mat LightingCorrection(Mat& changeFrom, const Mat& changeTo,bool saveCodeLength = true);
  EXPORTLIB Mat LightingCorrection(Mat& changeFrom, const Mat& changeTo, const Tensor<double,1>& VQCodebook);
  EXPORTLIB int GetLightingCodeLength() const;
  EXPORTLIB void SetLightingCodeLength(int l);
  //EXPORTLIB Mat ComputeTPSS(const Mat& im, double p) const;//thin plate spline smoothing
  EXPORTLIB void RecordLighting(void);
  EXPORTLIB Mat SearchCodeword(const Mat& val, const Mat& VQCodeBook);
  EXPORTLIB vector<Vec<double,1> > GetTagLighting() const;
  EXPORTLIB vector<Vec<double,1> > GetCanLighting() const;
  int codeLength;
  vector<double> lightingDCTCoeffStat;
  vector<Vec<double,1> > lightTag;
  vector<Vec<double,1> > lightCan;

};



EXPORTLIB Mat ComputeTPSS(const Mat& im, double p);
EXPORTLIB Mat ComputeSAT(const Mat& im);//compute sumed-area-table
EXPORTLIB Mat ComputeSAT(const Mat& S, const Point3i& sPos, const Point3i& ePos);

}
#endif // LIGHTING_H

#ifndef IO_H
#define IO_H
#ifdef WIN32
#define EXPORTLIB __declspec(dllexport)
#else
#define EXPORTLIB
#endif
#include <TensorLite.h>
namespace metric{
 EXPORTLIB bool readMahalanobis(vector<Mat>& covars,string filename); // 1D
 EXPORTLIB bool readMahalanobis2(vector<Mat>& covars,string filename); //2D

}
#endif // IO_H



#ifndef FILTERBANK_H
#define FILTERBANK_H
#include "TensorLite.h"

class FilterBank{
public:
  typedef double value_type;
  typedef tensor::Tensor<value_type,1> real_data_type;
  typedef tensor::Tensor<value_type,2> data_type;
  typedef const data_type&  c_data_ref;
  typedef data_type& data_ref;
  typedef const real_data_type& c_real_ref;

  virtual ~FilterBank(){}
  virtual int decompose(void)= 0;
  virtual vector<data_type>& getSpaceDomainPyr(void) = 0;
  virtual int buildFilters(int maxscale, int K) = 0;
  virtual int expand(int rh, int rw, int border,value_type value)=0;
  virtual int updateData(c_data_ref data)=0;
  virtual int getDir()=0;
  virtual int getLevel()=0;
protected:
  virtual int deleteFilters(void)=0;


};

#endif // FILTERBANK_H

/* 
   Copyright (C) 1988 Free Software Foundation
   written by Doug Lea (dl@rocky.oswego.edu)

   This file is part of the GNU C++ Library.  This library is free
   software; you can redistribute it and/or modify it under the terms of
   the GNU Library General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your
   option) any later version.	This library is distributed in the hope
   that it will be useful, but WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the GNU Library General Public License for more details.
   You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free Software
   Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef __FComplex_h__
#define __FComplex_h__ 1

#define __ATT_FComplex__

#include <iostream>
#include <cmath>
#include "precision.h"

using std::istream;
using std::ostream;
using std::ws;

class FComplex
{
#ifdef __ATT_FComplex__
 public:
#else
 protected:
#endif

  Real re;
  Real im;

 public:

  FComplex() {}
  FComplex(Real r, Real i=0) : re(r), im(i) {}
  FComplex(const FComplex& y) : re(y.re), im(y.im) {}
	
  ~FComplex() {}

  Real real() const {return re;}
  Real imag() const {return im;}

  const FComplex& operator = (const FComplex& y);
	
  const FComplex& operator += (const FComplex& y);
  const FComplex& operator += (Real y);
  const FComplex& operator -= (const FComplex& y);
  const FComplex& operator -= (Real y);
  const FComplex& operator *= (const FComplex& y);
  const FComplex& operator *= (Real y);
  const FComplex& operator /= (const FComplex& y);
  const FComplex& operator /= (Real y);
	
  void error(char* msg) const;
};

// inline members

inline const FComplex& FComplex::operator = (const FComplex& y)

{ 
  re = y.re; im = y.im; return *this; 
} 

inline const FComplex& FComplex::operator += (const FComplex& y)
{ 
  re += y.re;  im += y.im; return *this; 
}

inline const FComplex& FComplex::operator += (Real y)
{ 
  re += y; return *this; 
}

inline const FComplex& FComplex::operator -= (const FComplex& y)
{ 
  re -= y.re;  im -= y.im; return *this; 
}

inline const FComplex& FComplex::operator -= (Real y)
{ 
  re -= y; return *this; 
}

inline const FComplex& FComplex::operator *= (const FComplex& y)
{  
  Real r = re * y.re - im * y.im;
  im = re * y.im + im * y.re; 
  re = r; 
  return *this; 
}

inline const FComplex& FComplex::operator *= (Real y)
{  
  re *= y; im *= y; return *this; 
}

inline const FComplex& FComplex::operator /= (const FComplex& y)
{
  register double t1,t2,t3;
  t2=1.0/(y.re*y.re+y.im*y.im);
  t1=t2*y.re; t2 *= y.im; t3=re;
  re *= t1; re += im*t2;
  im *= t1; im -= t3*t2;
  return *this;
}

inline const FComplex& FComplex::operator /= (Real y)
{
  re /= y;
  im /= y;
  return *this;
}

//	functions

inline int operator == (const FComplex& x, const FComplex& y)
{
  return x.re == y.re && x.im == y.im;
}

inline int operator == (const FComplex& x, Real y)
{
  return x.im == 0.0 && x.re == y;
}

inline int operator != (const FComplex& x, const FComplex& y)
{
  return x.re != y.re || x.im != y.im;
}

inline int operator != (const FComplex& x, Real y)
{
  return x.im != 0.0 || x.re != y;
}

inline FComplex operator - (const FComplex& x)
{
  return FComplex(-x.re, -x.im);
}

inline FComplex conj(const FComplex& x)
{
  return FComplex(x.re, -x.im);
}

inline FComplex operator + (const FComplex& x, const FComplex& y)
{
  return FComplex(x.re+y.re, x.im+y.im);
}

inline FComplex operator + (const FComplex& x, Real y)
{
  return FComplex(x.re+y, x.im);
}

inline FComplex operator + (Real x, const FComplex& y)
{
  return FComplex(x+y.re, y.im);
}

inline FComplex operator - (const FComplex& x, const FComplex& y)
{
  return FComplex(x.re-y.re, x.im-y.im);
}

inline FComplex operator - (const FComplex& x, Real y)
{
  return FComplex(x.re-y, x.im);
}

inline FComplex operator - (Real x, const FComplex& y)
{
  return FComplex(x-y.re, -y.im);
}

inline FComplex operator * (const FComplex& x, const FComplex& y)
{
  return FComplex(x.re*y.re-x.im*y.im, x.re*y.im+x.im*y.re);
}

inline FComplex multconj(const FComplex& x, const FComplex& y)
{
  return FComplex(x.re*y.re+x.im*y.im,x.im*y.re-x.re*y.im);
}

inline FComplex operator * (const FComplex& x, Real y)
{
  return FComplex(x.re*y, x.im*y);
}

inline FComplex operator * (Real x, const FComplex& y)
{
  return FComplex(x*y.re, x*y.im);
}

inline FComplex operator / (const FComplex& x, const FComplex& y)
{
  register double t1,t2;
  t2=1.0/(y.re*y.re+y.im*y.im);
  t1=t2*y.re; t2 *= y.im;
  return FComplex(x.im*t2+x.re*t1, x.im*t1-x.re*t2);
}

inline FComplex operator / (const FComplex& x, Real y)
{
  return FComplex(x.re/y,x.im/y);
}

inline FComplex operator / (Real x, const FComplex& y)
{
  register double factor;
  factor=1.0/(y.re*y.re+y.im*y.im);
  return FComplex(x*y.re*factor,-x*y.im*factor);
}

inline Real real(const FComplex& x)
{
  return x.re;
}

inline Real imag(const FComplex& x)
{
  return x.im;
}

inline Real abs2(const FComplex& x)
{
  return x.re*x.re+x.im*x.im;
}

inline Real abs(const FComplex& x)
{
  return sqrt(abs2(x));
}

inline Real arg(const FComplex& x)
{
  return x.im != 0.0 ? atan2(x.im, x.re) : 0.0;
}

// Return the principal branch of the square root (non-negative real part).
inline FComplex sqrt(const FComplex& x)
{
  Real mag=abs(x);
  if(mag == 0.0) return FComplex(0.0,0.0);
  else if(x.re > 0) {
    Real re=sqrt(0.5*(mag+x.re));
    return FComplex(re,0.5*x.im/re);
  } else {
    Real im=sqrt(0.5*(mag-x.re));
    if(x.im < 0) im=-im;
    return FComplex(0.5*x.im/im,im);
  }
}

inline FComplex polar(Real r, Real t)
{
  return FComplex(r*cos(t), r*sin(t));
}

// FComplex exponentiation
inline FComplex pow(const FComplex& z, const FComplex& w)
{
  Real u=w.re;
  Real v=w.im;
  if(z == 0.0) return w == 0.0 ? 1.0 : 0.0;
  Real logr=0.5*log(abs2(z));
  Real th=arg(z);
  Real phi=logr*v+th*u;
  return exp(logr*u-th*v)*FComplex(cos(phi),sin(phi));
}

inline FComplex pow(const FComplex& z, Real u)
{
  if(z == 0.0) return u == 0.0 ? 1.0 : 0.0;
  Real logr=0.5*log(abs2(z));
  Real theta=u*arg(z);
  return exp(logr*u)*FComplex(cos(theta),sin(theta));
}

inline istream& operator >> (istream& s, FComplex& y)
{
  char c;
  s >> ws >> c;
  if(c == '(') {
    s >> y.re >> c;
    if(c == ',') s >> y.im >> c;
    else y.im=0.0;
  } else {
    s.putback(c);
    s >> y.re; y.im=0.0;
  }
  return s;
}

inline ostream& operator << (ostream& s, const FComplex& y)
{
  s << "(" << y.re << "," << y.im << ")";
  return s;
}

inline bool isfinite(const FComplex& z)
{
#ifdef _WIN32
  return _finite(z.re) && _finite(z.im);
#else  
  return !(std::isinf(z.re) || std::isnan(z.re) || std::isinf(z.im) || std::isnan(z.im));
#endif  
}

#endif

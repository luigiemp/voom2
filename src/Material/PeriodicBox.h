// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__PeriodicBox_h__)
#define __PeriodicBox_h__

#include "VoomMath.h"
#include "Node.h"

using namespace std;
using namespace voom;

namespace voom
{
  class PeriodicBox {
  public:
    
    PeriodicBox(double hx, double hy) { _h = hx, hy;}
    
    virtual void setSize(Vector2D & sz) { _h = sz; }
    
    virtual void setSize(double hx, double hy) { _h = hx, hy; }
    
    virtual void mapPoint(Vector2D & x) const { 
      for(int i=0; i<2; i++) {
	x(i) = fmod( x(i) , _h(i) );
	if( x(i) < 0.0 ) x(i) += _h(i);
      }
      return;
    }

    virtual void mapPointX(Vector2D & x) const {
      x(0) = fmod( x(0) , _h(0) );
      if( x(0) < 0.0 ) x(0) += _h(0);
      return;
    }

    virtual void mapPointY(Vector2D & x) const {
      x(1) = fmod( x(1) , _h(1) );
      if( x(1) < 0.0 ) x(1) += _h(1);
      return;
    }
    
    virtual void mapDistance(Vector2D & dx) const {
      for(int i=0; i<2; i++) {
	dx(i) = dx(i) - 
	  ( dx(i)>0.0 ? std::floor(dx(i)/_h(i)+0.5) : std::ceil(dx(i)/_h(i)-0.5))*_h(i);
	  }      
      return;
    }
    
    virtual bool inside(Vector2D & x) const {
      for(int i=0; i<2; i++) {    
	if( x(i) < 0.0 || x(i) > _h(i) ) return(false);
      }
      return(true);
    }
    
    virtual void stretchX(double strtchFact) {
      _h[0] = strtchFact*_h[0];
    }
    
    virtual void stretchY(double strtchFact) {
      _h[1] = strtchFact*_h[1];
    }

    virtual void stretchXY(double strtchFact) {
      _h[0] = strtchFact*_h[0];
      _h[1] = strtchFact*_h[1];
    }

    virtual const Vector2D & size() const {return _h;}

    virtual void setShear(double shear) {;}

    virtual void setShearX(double shear) {;}

    virtual void setShearY(double shear) {;}

  protected:
    Vector2D _h;
  };
};
#endif //__PeriodicBox_h__

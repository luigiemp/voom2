// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__LeesEdwards_h__)
#define __LeesEdwards_h__

#include "VoomMath.h"
#include "PeriodicBox.h"

using namespace std;
using namespace voom;

namespace voom
{
  /*! A periodic cell with parallelagram shape, obtained by simple
      shear of a rectangular cell (PeriodicBox).  Distances and images are
      computed by using the forward and inverse mappings for simple
      shear:

      x_1 = X_1 + \Lambda X_2 
      x_2 = X_2

      X_1 = x_1 - \Lambda x_2
      X_2 = x_2
    
      The idea is to apply the inverse mapping to compute a
      rectangular version of the box (X_1, X_2) from the current
      sheared version (x_1, x_2), then apply the standard formulas for
      rectangular periodic BCs to the positions in the rectangular
      domain, and finally transform those back to sheared domain.

      For distance vectors, the mapping works the same as follows:
      U = Y - X
      u = y - x

      u_1 = y_1 - x_1 
          = Y_1 + \Lamba Y_2 - X_1 - \Lambda X_2
          = U_1 + \Lambda U_2
     
      u_2 = y_2 - x_2
          = Y_2 - Y_2
	  = U_2
   */
  class LeesEdwards : public PeriodicBox {
  public:
    
    LeesEdwards(double hx, double hy, double shear=0.0, int shearD=0) : PeriodicBox(hx,hy), _shear(shear), _shearDirection(shearD) {} 
    
    void mapPoint(Vector2D & x) const { 
      Vector2D X;
      inverse(x,X);      
      for(int i=0; i<2; i++) {
	X(i) = fmod( X(i) , _h(i) );
	if( X(i) < 0.0 ) X(i) += _h(i);
      }

      mapping(X,x);

      return;
    }
    
    void mapDistance(Vector2D & dx) const {
      Vector2D dX;
      inverse(dx,dX);
      for(int i=0; i<2; i++) {
	dX(i) = dX(i) - 
	  ( dX(i)>0.0 ? std::floor(dX(i)/_h(i)+0.5) : std::ceil(dX(i)/_h(i)-0.5))*_h(i);
      }      
      mapping(dX,dx);
      return;
    }

    double shear() const {return _shear;}

    void setShear(double shear) {
      setShearX(shear);
    }
    
    void setShearX(double shear) { 
      _shear = shear;
      _shearDirection = 0;
    }

    void setShearY(double shear) {
      _shear = shear;
      _shearDirection = 1;
    }

    bool inside(Vector2D & x) const {
      Vector2D X;
      inverse(x,X);
      for(int i=0; i<2; i++) {    
	if( X(i) < 0.0 || X(i) > _h(i) ) return(false);
      }
      return(true);
    }

    inline void mapping(const Vector2D & X, Vector2D & x) const {
      if(_shearDirection == 0) {
	x(0) = X(0) + _shear*X(1);
	x(1) = X(1);
      }
      else {
	x(0) = X(0);
	x(1) = X(1) - _shear*X(0);
      }
    }

    inline void inverse(const Vector2D & x, Vector2D & X) const {
      if(_shearDirection == 0) {
	X(0) = x(0) - _shear*x(1);
	X(1) = x(1);
      }
      else {
	X(0) = x(0);
	X(1) = x(1) + _shear*x(0);
      }
    }

  protected:
    double _shear;
    int _shearDirection;
  };
};
#endif //__LeesEdwards_h__

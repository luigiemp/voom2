// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Andrew R. Missel
//                University of California Los Angeles
//                 (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

// TO DO: //
// (1) 


#if !defined(__Motor_h__)
#define __Motor_h__

#include "Node.h"
#include "Element.h"
#include "VoomMath.h"
#include "Spring.h"
#include <iostream>

using namespace tvmet;
using namespace std;
using namespace voom;

namespace voom
{

  template<int N>
  class Motor : public Element {
    
  public:
    
    typedef tvmet::Vector<double,N> VectorND;
    typedef BrownianNode<N> Node_t;
    typedef typename std::vector< Node_t* > Node_tContainer;
    typedef typename Node_tContainer::iterator Node_tIterator;
    
    Motor(const VectorND & motPos, double k) 
      : _startPos(motPos), _k(k) { 

      _d0 = 0.0;
      
      _rng.seed((unsigned int)time(0));
    }
    
    Motor(const VectorND & motPos, double k, double d0) : _startPos(motPos), _k(k), _d0(d0) {

      _rng.seed((unsigned int)time(0));
    }

    
    void setFVParams(double v0, double Fmax, double vp_Fmax, double vm_Fmax, double negV_slope) {
      _v0 = v0; 
      _Fmax = Fmax; 
      _vp_Fmax = vp_Fmax; 
      _vm_Fmax = vm_Fmax; 
      _negV_slope = negV_slope; 
    }
	
    void setDetachParams(double r0, double alpha) { 	
      _r0 = r0;
      _alpha = alpha; 
    }
	
    double getVelocity(double force) {
      double v;
      if(-force < _Fmax) {
	v = _v0 + (force/_Fmax)*(_v0 - _vp_Fmax);
      }
      else {
	v = _vm_Fmax + _negV_slope*(force + _Fmax); 
      }
      return v;
    }
    
    double getDetachProb(double force, double dt) {
      if(-force >= _Fmax) {
	return 1.0;
      }
      else {
	return 1.0 - exp(-_r0*pow(_Fmax/(_Fmax+force),_alpha)*dt);
      }
    }
    
    void setTimeStep(double dt) {
      _dt = dt;
    }
    
    VectorND getEndPoint1() {
      VectorND endPoint1;
      if(_epsi1 >= 0.0) {
	const VectorND & x1A = (*(_Fil1+_Fil1pos))->point();
	const VectorND & x1B = (*(_Fil1+_Fil1pos+1))->point();
	endPoint1 = x1A*(1.0-_epsi1)+x1B*_epsi1;
      }
      else {
	endPoint1 = _pos;
      }
      return endPoint1;
    }
    
    VectorND getEndPoint2() {
      VectorND endPoint2;
      if(_epsi2 >= 0.0) {
	const VectorND & x2A = (*(_Fil2+_Fil2pos))->point();
	const VectorND & x2B = (*(_Fil2+_Fil2pos+1))->point();
	endPoint2 = x2A*(1.0-_epsi2)+x2B*_epsi2;
      }
      else {
	endPoint2 = _pos;
      }
      return endPoint2;
    }
    
    VectorND getStartPoint() {
      return _startPos;
    }
    
    void setPosition(Node_tContainer & Fil1Cont, int Fil1pos, Node_tContainer & Fil2Cont, int Fil2pos, VectorND & intersection1, VectorND & intersection2) {
      _Fil1 = Fil1Cont.begin();
      _Fil1pos = Fil1pos;
      _Fil1size = Fil1Cont.size();
      _Fil2 = Fil2Cont.begin();
      _Fil2pos = Fil2pos;
      _Fil2size = Fil2Cont.size();
   
      _epsi1 = norm2(intersection1-(*(_Fil1+_Fil1pos))->point())/norm2((*(_Fil1+_Fil1pos+1))->point()-(*(_Fil1+_Fil1pos))->point());
      _epsi2 = norm2(intersection2-(*(_Fil2+_Fil2pos))->point())/norm2((*(_Fil2+_Fil2pos+1))->point()-(*(_Fil2+_Fil2pos))->point());
      _pos = (intersection1+intersection2)/2.0;
      _attached = true;
    }
    
    void stepMotor() {
      const VectorND & x1A = (*(_Fil1+_Fil1pos))->point();
      const VectorND & x1B = (*(_Fil1+_Fil1pos+1))->point();
      const VectorND & x2A = (*(_Fil2+_Fil2pos))->point();
      const VectorND & x2B = (*(_Fil2+_Fil2pos+1))->point();
      VectorND endPoint1; 
      endPoint1 = x1A*(1.0-_epsi1)+x1B*_epsi1;
      VectorND endPoint2;
      endPoint2 = x2A*(1.0-_epsi2)+x2B*_epsi2;
      double d  = norm2(endPoint2-endPoint1);
      double f = -_k*(d-_d0);
      double rn = _rng.random();
      if(rn < getDetachProb(f,_dt)) { // check to see if the motor should detach //
	detach();
      }
      else { // if not, step the motor and reset position, etc. (checking to make sure it hasn't stepped off) //
	bool doDetach = false;
	double len1 = norm2(x1B-x1A);
	double len2 = norm2(x2B-x2A);
	double cos_Fil1 = 0.0;
	double cos_Fil2 = 0.0;
	if(d != 0.0) {
	  cos_Fil1 = dot(endPoint1-endPoint2,x1B-x1A)/(d*len1);
	  cos_Fil2 = dot(endPoint2-endPoint1,x2B-x2A)/(d*len2);
	}
	double v1 = getVelocity(f*cos_Fil1);
	double v2 = getVelocity(f*cos_Fil2);
	double dx1 = v1*_dt;
	double dx2 = v2*_dt;
	
	_epsi1 += dx1/len1;
	endPoint1 = x1A*(1.0-_epsi1)+x1B*_epsi1;
	if(_epsi1 > 1.0) { // check to see if we have stepped over the edge of the current rod //
	  if(_Fil1pos < _Fil1size - 2) {
	    const VectorND & x1C = (*(_Fil1+_Fil1pos+2))->point();
	    double len1prime = norm2(x1C-x1B);
	    double cos_Fil1prime = dot(x1B-x1A,x1C-x1B)/(len1*len1prime);
	    double dx1prime = (_epsi1 - 1.0)*len1*cos_Fil1prime;
	    _epsi1 = dx1prime/len1prime;
	    _Fil1pos++;
	    endPoint1 = x1B*(1.0-_epsi1)+x1C*_epsi1;
	  }
	  else {
	    doDetach = true;
	  }
	}

	_epsi2 += dx2/len2;
	endPoint2 = x2A*(1.0-_epsi2)+x2B*_epsi2;
	if(_epsi2 > 1.0) { // check to see if we have stepped over the edge of the current rod //
	  if(_Fil2pos < _Fil2size - 2) {
	    const VectorND & x2C = (*(_Fil2+_Fil2pos+2))->point();
	    double len2prime = norm2(x2C-x2B);
	    double cos_Fil2prime = dot(x2B-x2A,x2C-x2B)/(len2*len2prime);
	    double dx2prime = (_epsi2 - 1.0)*len2*cos_Fil2prime;
	    _epsi2 = dx2prime/len2prime;
	    _Fil2pos++;
	    endPoint2 = x2B*(1.0-_epsi2)+x2C*_epsi2;
	  }
	  else {
	    doDetach = true;
	  }
	}

	_pos = (endPoint1+endPoint2)/2.0;
	if(doDetach == true) {
	  detach();
	}
      }
    }
    
    void detach() {
      _attached = false;
      _epsi1 = -1.0;
      _epsi2 = -1.0;
      std::cout << "detached!" << std::endl;
    }
    
    bool isAttached() {
      return _attached;
    }

    void compute(bool f0, bool f1, bool f2) {
      const VectorND & x1A = (*(_Fil1+_Fil1pos))->point();
      const VectorND & x1B = (*(_Fil1+_Fil1pos+1))->point();
      const VectorND & x2A = (*(_Fil2+_Fil2pos))->point();
      const VectorND & x2B = (*(_Fil2+_Fil2pos+1))->point();
      VectorND endPoint1;
      endPoint1 = x1A*(1.0-_epsi1)+x1B*_epsi1;
      VectorND endPoint2; 
      endPoint2 = x2A*(1.0-_epsi2)+x2B*_epsi2;
      double d  = norm2(endPoint2-endPoint1);

      if(f0) {
	_energy = 0.5*_k*sqr(d-_d0);
      }
      
      if(f1) { // assuming that the motor has already stepped, compute the new force and add this force to the proper nodes //
	for(int i=0; i<N; i++) {
	  double f = -_k*(d-_d0)*(endPoint1(i)-endPoint2(i))/d;
	  (*(_Fil1+_Fil1pos))->addForce(i, -f*(1.0-_epsi1));
          (*(_Fil1+_Fil1pos+1))->addForce(i, -f*_epsi1);
          (*(_Fil2+_Fil2pos))->addForce(i, f*(1.0-_epsi2));
	  (*(_Fil2+_Fil2pos+1))->addForce(i, f*_epsi2);
	}
      }

      return;
    }

    double stiffness() const {return _k;}

    void setStiffness(double k) { _k = k; }

  private:
    
    VectorND _pos;
    VectorND _startPos;
    double _radOfInfluence;
    Node_tIterator _Fil1;
    int _Fil1pos;
    int _Fil1size;
    Node_tIterator _Fil2;
    int _Fil2pos;
    int _Fil2size;
    double _epsi1;
    double _epsi2;
    double _k;
    double _d0;
    double _dt;
    
    bool _attached;
    
    ranlib::Uniform<double> _rng;
    
    // these static variables contain information on the force-velocity curve, which is common to all Motor objects //
    double _v0;
    double _Fmax;
    double _vp_Fmax;
    double _vm_Fmax;
    double _negV_slope;
    
    // these static variables contain information on the detachment probability //
    
    double _r0;
    double _alpha;
    
  };
};

#endif // __Motor_h__

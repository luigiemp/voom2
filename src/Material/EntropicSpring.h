// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Andrew R. Missel
//                University of California Los Angeles
//                 (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Entropic_Spring_h__)
#define __Entropic_Spring_h__

#include "Node.h"
#include "Spring.h"
#include "VoomMath.h"
#include <cmath>

namespace voom
{

  template<int N>
  class EntropicSpring : public Spring<N> {
    
  public: 

    using Spring<N>::_nodeA;
    using Spring<N>::_nodeB;
    using Spring<N>::_k;
    using Spring<N>::_d0;
    using Spring<N>::_energy;
    using Spring<N>::_baseNodes;

    typedef tvmet::Vector<double,N> VectorND;
    typedef DeformationNode<N> Node_t;


    EntropicSpring(Node_t * nodeA, Node_t * nodeB, double entropk, double Lp, double L0, double mu) : Spring<N> (nodeA,nodeB,entropk), _mu(mu), _Lp(Lp), _entropk(entropk) {
      
      _lin = false;
      
      VectorND sep;
      sep = nodeA->position() - nodeB->position();
      _L0 = norm2(sep);
      _d0 = _L0;

      _Larc = _L0*(1.0 + (_L0/(6.0*_Lp)));

      _kel = _mu/_Larc;

      double a1 = 90.0*_entropk*_Lp*(1.0+(_L0/(3.0*_Lp)))/pow(_L0,4);

      _kent = a1;

      // now find the extension and force at the elastic stiffness //

      double startForce = 1.0;
      double curStiff = stifffromForce(startForce);
      while(curStiff < _kel) {
	startForce *= 2.0;
	curStiff = stifffromForce(startForce);
      }
      double maxF = startForce;
      double minF = 0.0;
      double curF = 0.5*(minF+maxF);
      curStiff = stifffromForce(curF);
      while(fabs(maxF-minF)/(fabs(maxF)+fabs(minF)) > 1.0e-6) {
	if(curStiff < _kel) {
	  minF = curF;
	}
	else {
	  maxF = curF;
	}
	curF = 0.5*(maxF+minF);
	curStiff = stifffromForce(curF);
	
	// DEBUG //

	// std::cout << "L = " << LfromForce(curF) << ", F = " << curF << ", k = " << curStiff << std::endl;
      }

      if(_kent >= _kel) {
	_force_maxL = 0.0;
	_energy_maxL = 0.0;
	_maxL = _L0;
      }

      else {
	_maxL = LfromForce(curF);
	_force_maxL = curF;
      }

      // std::cout << "Elastic stiffness = " << _kel << ", current stiffness = " << curStiff << std::endl;
     
      // now find length at which buckling should set in //
      double maxL = _L0;
      double minL = 0.0;
      double curL = 0.5*(maxL + minL);
      curF = forcefromL(curL,-sqr(M_PI/curL)*_entropk,0.0);
      double minLstep = (maxL-minL)/1000.0;
      double stepL = minLstep < (maxL-minL)/2.0 ? minLstep : (maxL-minL)/2.0;
      double sd = stiffDeriv(curF,curL,stepL);
      curStiff = stifffromForce(curF,curL);
      
      while(fabs(maxL-minL)/(fabs(maxL)+fabs(minL)) > 1.0e-6) {
	if(sd > 0.0) {
	  maxL = curL;
	}
	else {
	  minL = curL;
	}
	curL = 0.5*(maxL+minL);
	curF = forcefromL(curL,-sqr(M_PI/curL)*_entropk,0.0);

	stepL = minLstep < (maxL-minL)/2.0 ? minLstep : (maxL-minL)/2.0;

	sd = stiffDeriv(curF,curL,stepL);
	curStiff = stifffromForce(curF,curL);
      }
      
      _medL = curL;
      _force_medL = curF;
      _kbuck = curStiff;

      // now we can compute the coefficients //

      double dlmax = _maxL - _L0;
      double dlmed = _medL - _L0;
      double Fmax = _force_maxL;
      double Fmed = _force_medL;
      double kbuck = _kbuck;
      double kel = _kel;

      if(_kent >= _kel) {
	_alph1 = 0.0;
      }

      else {
	_alph1 = (_kel-_kent)/dlmax;
      }

      if(_kent >= _kel && _kbuck < _kel) {
	_alph2 = (_kbuck-_kel)/dlmed;
      }
      else if(_kent >= _kel && _kbuck >= _kel) {
	std::cout << "Error!  Buckling stiffness greater than elastic stiffness." << std::endl;
	_alph2 = 0.0;
	_lin = true;
      }
      else {
	_alph2 = (_kbuck-_kent)/dlmed;
      }

      if(_kent >= _kel) {
	_kent = _kel;
      }
      else {
	_force_maxL = _kent*dlmax + (_alph1/2.0)*sqr(dlmax);
	_energy_maxL = (_kent/2.0)*sqr(dlmax) + (_alph1/6.0)*pow(dlmax,3);
      }      

      //  _energy_maxL = (_alph2/6.0)*pow(dlmax,3) + (_kent/_alph1)*((1.0/_alph1)*(exp(_alph1*dlmax)-1.0)-dlmax);
      
      //_energy_maxL = (_kent/_alph)*((1.0/_alph)*(exp(_alph*dlmax)-1.0)-dlmax);
      
      // std::cout << "Force at Lmax = " << Fmax << "; computed version is " << _force_maxL << std::endl;
      
      // now, medL //

      if(_kbuck > _kel) {
	_force_medL = 0.0;
	_energy_medL = 0.0;
	_medL = _L0;
      }
      else {
	_force_medL = _kent*dlmed + (_alph2/2.0)*sqr(dlmed);
	_energy_medL = (_kent/2.0)*sqr(dlmed) + (_alph2/6.0)*pow(dlmed,3);
      }      

      //_force_medL = _kbuck*dlmed + ((_kent-_kbuck)/_beta)*(exp(_beta*dlmed)-1.0);
      //_energy_medL = (_kbuck/2.0)*sqr(dlmed) + ((_kent-_kbuck)/_beta)*((1.0/_beta)*(exp(_beta*dlmed)-1.0) - dlmed);
      
      // std::cout << "Force at Lmed = " << Fmed << "; computed version is " << _force_medL << std::endl;
      
      // now determine Lmin //
      
      // for now, set Lmin to 0 //
      
      _minL = 0.0;
      
      // set force at Lmin... //
      
      _force_minL = _force_medL - kbuck*(_medL-_minL);
      
      
      // print out info for debugging //
      
      //       std::cout << "Force @ maximum stiffness = " << Fmax << ", computed version = " << _force_maxL << std::endl
      // 		<< "Force @ minimum stiffness = " << Fmed << ", computed version = " << _force_medL << std::endl;
      
      if(!_lin) {
        assert(_force_maxL >= 0.0);
        assert(_force_medL <= 0.0);
	
        assert(_maxL >= _L0);
        assert(_medL <= _L0 && _medL >= 0.0);
      }
      //checkConsistency();
      
    }

    void setLin(double mult) {
      _lin = true;
      _kent *= mult;
    }

    bool getLin() {
      return _lin;
    }

    double forcefromL(double L, double Fmin, double Fmax) {
      double Fcur = 0.5*(Fmin+Fmax);
      double lhs = forceLHS(L,Fcur);
      double rhs = forceRHS(L,Fcur);

      double meas = 2.0*fabs(Fmax-Fmin)/(fabs(Fmax)+fabs(Fmin));
      while(meas > 1.0e-6) {
	if(lhs>rhs) {
	  if(L < _L0) {
	    Fmin = Fcur;
	  }
	  else {
	    Fmax = Fcur;
	  }
	}
	else {
	  if(L < _L0) {
	    Fmax = Fcur;
	  }
	  else {
	    Fmin = Fcur;
	  }
	}
	Fcur = 0.5*(Fmax + Fmin);
	lhs = forceLHS(L,Fcur);
	rhs = forceRHS(L,Fcur);
	
	meas = 2.0*fabs(Fmax-Fmin)/(fabs(Fmax)+fabs(Fmin));
      }

      return Fcur;
    }

//     double LfromForce(double f, double Lmin, double Lmax) {
//       double Lcur = 0.5*(Lmin + Lmax);
      
//       if(f>0.0) {
// 	double lhs = forceLHS(Lcur,f);
// 	double rhs = forceRHS(Lcur,f);
// 	double meas = 2.0*fabs(lhs-rhs)/(fabs(lhs)+fabs(rhs));
// 	while(meas > 1.0e-3) {
	  
// 	}
//       }

//     }

    double LfromForce(double f) {
      double L2r;
      if(f>=0.0) {
	double minL = _L0;
	double maxL = _Larc;
	double curL = 0.5*(minL+maxL);
	double lhs = forceLHS(curL,f);
	double rhs = forceRHS(curL,f);
	double meas = 2.0*fabs(maxL-minL)/(fabs(minL)+fabs(maxL));
	while(meas > 1.0e-6) {
	  if(lhs>rhs) {
	    minL = curL;
	  }
	  else {
	    maxL = curL;
	  }
	  curL = 0.5*(maxL+minL);
	  lhs = forceLHS(curL,f);
	  rhs = forceRHS(curL,f);
	  meas = 2.0*fabs(maxL-minL)/(fabs(minL)+fabs(maxL));
	}
	
	L2r = curL;
      }
      else {
	double minL = 0.0;
	double maxL = _L0;
	double curL = 0.5*(minL+maxL);
	double lhs = forceLHS(curL,f);
	double rhs = forceRHS(curL,f);
	double meas = 2.0*fabs(maxL-minL)/(fabs(minL)+fabs(maxL));
	while(meas > 1.0e-6) {
	  if(lhs>rhs) {
	    maxL = curL;
	  }
	  else {
	    minL = curL;
	  }
	  curL = 0.5*(maxL+minL);
	  lhs = forceLHS(curL,f);
	  rhs = forceRHS(curL,f);
	  meas = 2.0*fabs(maxL-minL)/(fabs(minL)+fabs(maxL));
	}
	
	L2r = curL;
      }

      return L2r;
    }
    
    double stiffDeriv(double f) {
      double smallStep = (1.0e-4)*_force_maxL;
      double l1=LfromForce(f-(smallStep/2.0));
      double l2=LfromForce(f+(smallStep/2.0));
      double k1=stifffromForce(f-(smallStep/2.0));
      double k2=stifffromForce(f+(smallStep/2.0));

      return (k2-k1)/(l2-l1);
    }

    double stiffDeriv(double f, double L) {
      double smallStep = (1.0e-4)*_Larc;
      double L1,f1;
      double L2,f2;
      if(f < 0.0) {
	f2 = f;
	L2 = L;
	L1 = L2 - smallStep;
	f1 = forcefromL(L1,-sqr(M_PI/L1)*_entropk,f2);
      }
      else {

      }

      double s1 = stifffromForce(f1,L1);
      double s2 = stifffromForce(f2,L2);

      return (s2-s1)/(L2-L1);

    }

    double stiffDeriv(double f, double L, double stepL) {
      double L1,f1;
      double L2,f2;
      if(f < 0.0) {
	f2 = f;
	L2 = L;
	L1 = L2 - stepL;
	f1 = forcefromL(L1,-sqr(M_PI/L1)*_entropk,f2);
      }
      else {

      }

      double s1 = stifffromForce(f1,L1);
      double s2 = stifffromForce(f2,L2);

      return (s2-s1)/(L2-L1);

    }
    
    double stifffromForce(double f) {
      double stiff2r;
      double y = sqrt(fabs(f)/_entropk);
      double L = LfromForce(f);
      if(f>=0) {
	double num = f + (_entropk/(2.0*_Lp))*y*((1.0/tanh(y*L)) - y*L*sqr(1.0/sinh(y*L)));
	double denom = _Larc - L - (L/(4.0*_Lp*y))*((1.0/tanh(y*L)) - y*L*sqr(1.0/sinh(y*L)));
	stiff2r = num/denom;
      }
      else {
	double num = f + (_entropk/(2.0*_Lp))*y*((1.0/tan(y*L)) - y*L*sqr(1.0/sin(y*L)));
	double denom = _Larc - L + (L/(4.0*_Lp*y))*((1.0/tan(y*L)) - y*L*sqr(1.0/sin(y*L)));
	stiff2r = num/denom;
      }

      return stiff2r;
    }

    double stifffromForce(double f, double L) {
      double stiff2r;
      double y = sqrt(fabs(f)/_entropk);
      if(f>=0) {
	double num = f + (_entropk/(2.0*_Lp))*y*((1.0/tanh(y*L)) - y*L*sqr(1.0/sinh(y*L)));
	double denom = _Larc - L - (L/(4.0*_Lp*y))*((1.0/tanh(y*L)) - y*L*sqr(1.0/sinh(y*L)));
	stiff2r = num/denom;
      }
      else {
	double num = f + (_entropk/(2.0*_Lp))*y*((1.0/tan(y*L)) - y*L*sqr(1.0/sin(y*L)));
	double denom = _Larc - L + (L/(4.0*_Lp*y))*((1.0/tan(y*L)) - y*L*sqr(1.0/sin(y*L)));
	stiff2r = num/denom;
      }

      return stiff2r;
    }
    
    double stiff2ndDeriv(double f) {
      double smallStep = 1.0e-8;
      double f1 = f-smallStep;
      double f2 = f+smallStep;
      double L1 = LfromForce(f1);
      double L = LfromForce(f);
      double L2 = LfromForce(f2);
      
      double k2 = stifffromForce(f2);
      double k = stifffromForce(f);
      double k1 = stifffromForce(f1);

      return (2.0/(L2-L1))*(((k2-k)/(L2-L))-((k-k1)/(L-L1)));
    }

    double forceLHS(double L, double f) {
      return (_Larc-L)*f;
    }

    double forceRHS(double L, double f) {
      double prefact = _entropk/(2.0*_Lp);
      double x = sqrt(fabs(f)/_entropk)*L;
      
      if(f>=0.0) {
	return (prefact*(x*(1.0/tanh(x))-1.0));
      }

      else {
	return (prefact*(x*(1.0/tan(x))-1.0));
      }
    }

    double alphLHS(double alph) {
      double dlmax = _maxL - _L0;
      
      return alph*(_force_maxL - (dlmax/2.0)*(_kel-_kent*exp(alph*dlmax)));
    }

    double alphRHS(double alph) {
      double dlmax = _maxL - _L0;
      
      return _kent*(exp(alph*dlmax)-1.0);
    }

    double betaLHS(double beta) {
      double dlmed = _medL - _L0;

      return beta*(_kbuck*dlmed - _force_medL)/(_kent - _kbuck);
    }

    double betaRHS(double beta) {
      double dlmed = _medL - _L0;

      return (1.0 - exp(beta*dlmed));
    }

    virtual void compute(bool f0, bool f1, bool f2) {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double L = norm2(xB-xA);
      
      if(f0) {
	_energy = computeEnergy(L);
      }
      
      if(f1) {
	double fmag = computeForce(L);
	for(int i=0; i<N; i++) {	  
	  double f = fmag*(xA(i)-xB(i))/L;
	  _nodeA->addForce(i, f);
	  _nodeB->addForce(i,-f);
	}
	
      }

      return;
    }

    virtual double stiffness() {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double L = norm2(xB-xA);

      return computeStiffness(L);
    }

    virtual double stiffness(double L) {
      return computeStiffness(L);
    }

    double computeEnergy(double L) {
      double en;

//       if(L>=_L0 && L<_maxL) {
// 	// en = (_alph2/6.0)*pow(L-_L0,3) + (_kent/_alph1)*((1.0/_alph1)*(exp(_alph1*(L-_L0))-1.0)-(L-_L0));
// 	en = (_kent/_alph)*((1.0/_alph)*(exp(_alph*(L-_L0))-1.0)-(L-_L0));
//       }
      
//       else if(L>=_maxL) {
// 	en = _energy_maxL + _force_maxL*(L-_maxL) + (_kel/2.0)*sqr(L-_maxL);
//       }

//       else if(L>=_minL && L<_L0) {
// 	en = (_kbuck/2.0)*sqr(L-_L0) + ((_kent-_kbuck)/_beta)*((1.0/_beta)*(exp(_beta*(L-_L0))-1.0) - (L-_L0));
//       }

//       else {

//       }

      if(!_lin) {

	if(L>=_L0 && L<_maxL) {
	  en = (_kent/2.0)*sqr(L-_L0) + (_alph1/6.0)*pow(L-_L0,3);
	}
	
	else if(L>=_maxL) {
	  en = _energy_maxL + _force_maxL*(L-_maxL) + (_kel/2.0)*sqr(L-_maxL);
	}
	
	else if(L>=_medL && L<_L0) {
	  en = (_kent/2.0)*sqr(L-_L0) + (_alph2/6.0)*pow(L-_L0,3);
	}
	
	else if(L>=_minL && L<_medL) {
	  en = _energy_medL + _force_medL*(L-_medL) + (_kbuck/2.0)*sqr(L-_medL);
	}
	
	else {
	  
	}
      }

      else {
	en = 0.5*_kent*sqr(L-_L0);
      }

      return en;
    }

    double computeForce(double L) {
      double f;

//       if(L>=_L0 && L<_maxL) {
// 	// f = (_alph2/2.0)*sqr(L-_L0) + (_kent/_alph1)*(exp(_alph1*(L-_L0))-1.0);
// 	f = (_kent/_alph)*(exp(_alph*(L-_L0))-1.0);
//       }
      
//       else if(L>=_maxL) {
// 	f = _force_maxL + _kel*(L-_maxL);
//       }

//       else if(L>=_minL && L<_L0) {
// 	f = _kbuck*(L-_L0) + ((_kent-_kbuck)/_beta)*(exp(_beta*(L-_L0))-1.0);
//       }

//       else {
// 	f = _force_minL + _kel*(L-_minL);
//       }

      if(!_lin) {

	if(L>=_L0 && L<_maxL) {
	  f = _kent*(L-_L0) + (_alph1/2.0)*sqr(L-_L0);
	}
	
	else if(L>=_maxL) {
	  f = _force_maxL + _kel*(L-_maxL);
	}
	
	else if(L>=_medL && L<_L0) {
	  f = _kent*(L-_L0) + (_alph2/2.0)*sqr(L-_L0);
	}
	
	else if(L>=_minL && L<_medL) {
	  f = _force_medL + _kbuck*(L-_medL);
	}
	
	else {
	  
	}
      }

      else {
	f = _kent*(L-_L0);
      }

      return f;
    }

    double computeStiffness(double L) {
      double k;

//       if(L>=_L0 && L<_maxL) {
// 	// k = _alph2*(L-_L0) + _kent*exp(_alph1*(L-_L0));
// 	k = _kent*exp(_alph*(L-_L0));
//       }
      
//       else if(L>=_maxL) {
// 	k = _kel;
//       }

//       else if(L>=_minL && L<_L0) {
// 	k = _kbuck + (_kent-_kbuck)*exp(_beta*(L-_L0));
//       }

//       else {
// 	k = _kel;
//       }

      if(!_lin) {
	
	if(L>=_L0 && L<_maxL) {
	  k = _kent + _alph1*(L-_L0);
	}
	
	else if(L>=_maxL) {
	  k = _kel;
	}
	
	else if(L>=_medL && L<_L0) {
	  k = _kent + _alph2*(L-_L0);
	}
	
	else if(L>=_minL && L<_medL) {
	  k = _kbuck;
	}
	
	else {
	  
	}
	
      }

      else {
	k = _kent;
      }

      if(k>_kel) {
	assert(((k-_kel)/_kel) < 1.0e-6);
	k = _kel;
      }

      return k;
    }

    void printOutVals(std::string fileName) {
      std::ofstream flname(fileName.c_str());
      flname << "#L\tF\tk\tE" << std::endl;
      
      double stepSize = (_maxL - _medL)/10000.0;

      double startPt = 0.5*_minL;
      double endPt = 0.5*(_maxL+_Larc);

      int nPts = (int)((endPt-startPt)/stepSize);
      
      for(int n=0; n<=nPts; n++) {
	double L = startPt + n*stepSize;
	flname << L << "\t" << computeForce(L) << "\t" << computeStiffness(L) << "\t" << computeEnergy(L) << std::endl;
      }

      flname.close();
    }

    virtual double getLength() {
      return _L0;
    }

//     virtual double intitialStiff() {
//       _k = computeStiffness(_L0);
//       return _k;
//     }


    virtual bool checkConsistency() {
      double stepSize = 1.0e-8;
      double L = 0.5*(_minL + _medL);

      double ehigh = computeEnergy(L+stepSize);
      double elow = computeEnergy(L);
      double forceComp = (ehigh-elow)/stepSize;
      double forceInt = computeForce(L+(stepSize/2.0));

      std::cout << "Checking consistency: force @ L = " << L << " is" << std::endl
		<< "computed from energy: " << forceComp << std::endl
		<< "computed from function: " << forceInt << std::endl << std::endl;

      L = 0.5*(_medL + _L0);
	    
      ehigh = computeEnergy(L+stepSize);
      elow = computeEnergy(L);
      forceComp = (ehigh-elow)/stepSize;
      forceInt = computeForce(L+(stepSize/2.0));
      
      std::cout << "Checking consistency: force @ L = " << L << " is" << std::endl
		<< "computed from energy: " << forceComp << std::endl
		<< "computed from function: " << forceInt << std::endl << std::endl;

      L = 0.5*(_maxL + _L0);
	    
      ehigh = computeEnergy(L+stepSize);
      elow = computeEnergy(L);
      forceComp = (ehigh-elow)/stepSize;
      forceInt = computeForce(L+(stepSize/2.0));
      
      std::cout << "Checking consistency: force @ L = " << L << " is" << std::endl
		<< "computed from energy: " << forceComp << std::endl
		<< "computed from function: " << forceInt << std::endl << std::endl;

      L = 0.5*(_maxL+ _Larc);
	    
      ehigh = computeEnergy(L+stepSize);
      elow = computeEnergy(L);
      forceComp = (ehigh-elow)/stepSize;
      forceInt = computeForce(L+(stepSize/2.0));
      
      std::cout << "Checking consistency: force @ L = " << L << " is" << std::endl
		<< "computed from energy: " << forceComp << std::endl
		<< "computed from function: " << forceInt << std::endl << std::endl;
    }


  protected:
    
    double _L0;
    double _entropk;
    double _Lp;
    double _Larc;
    double _mu;
    double _kel; // elastic stiffness //
    double _kbuck;
    double _kent;

    bool _lin;
    
    double _maxL; // length at which elastic stiffness takes over //
    double _force_maxL; // force at maximum length //
    double _energy_maxL;

    double _medL; // length at which buckling stiffness takes over //
    double _force_medL; // force at medium length //
    double _energy_medL;

    double _minL; // minimum length //
    double _force_minL; // force at minimum length //
    double _energy_minL;

    double _alph1;
    double _alph2;
    //double _alph;
    //double _beta;

    // double coeffs[5];
    
  };
};

#endif // __EntropicSpring_h__

// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Andrew R. Missel
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__TwoBodyPotential_h__)
#define __TwoBodyPotential_h__
#endif

#include "VoomMath.h"
#include "Element.h"
#include "PeriodicBox.h"
//#include "SemiflexibleGel.h"
#include "NodeBase.h"
#include "Node.h"

using namespace std;
using namespace voom;

namespace voom
{
  class TwoBodyPotential {
  public:

    typedef BrownianNode<2> DefNode;
    typedef std::vector<DefNode*> DefNodeContainer;
    
    enum FunctionType {
      expon,
      yukawa,
      gaussian
    };

    TwoBodyPotential() {}

    TwoBodyPotential(PeriodicBox * box, FunctionType ft, double maxDist, double prefact, double elemlength, double scale, int nChrgs) : _box(box), _ft(ft), _prefact(prefact), _elemlength(elemlength), _scale(scale), _nChrgs(nChrgs) {
      double stpRx = (maxDist)/99.0;
      double stpRy = (maxDist)/99.0;
      double stpTheta = M_PI/99.0;
      
      _Rxlims[0] = -maxDist/2.0;
      _Rxlims[1] = maxDist/2.0;
      _stpRx = stpRx;
     
      _Rylims[0] = -maxDist/2.0;
      _Rylims[1] = maxDist/2.0;
      _stpRy = stpRy;
      
      _thetalims[0] = 0.0;
      _thetalims[1] = M_PI;
      _stpTheta = stpTheta;
      makeTables();
    }

    void resetTables(PeriodicBox * box, double maxDist, FunctionType ft, double prefact, double elemlength, double scale) {
      _box = box;
      _ft = ft;
      _prefact = prefact;
      _elemlength = elemlength;
      _scale = scale;
      double stpRx = (maxDist)/99.0;
      double stpRy = (maxDist)/99.0;
      double stpTheta = M_PI/99.0;
      
      _Rxlims[0] = -maxDist/2.0;
      _Rxlims[1] = maxDist/2.0;
      _stpRx = stpRx;
     
      _Rylims[0] = -maxDist/2.0;
      _Rylims[1] = maxDist/2.0;
      _stpRy = stpRy;
      
      _thetalims[0] = 0.0;
      _thetalims[1] = M_PI;
      _stpTheta = stpTheta;
      makeTables();
    }

    void checkTable(DefNodeContainer & fil1, DefNodeContainer & fil2) {
      double energy = 0.0;

      DefNode* n1A = fil1[0];
      DefNode* n1B = fil1[1];
      DefNode* n2A = fil2[0];
      DefNode* n2B = fil2[1];

      Vector2D n1;
      Vector2D n2;
      
      n1 = n1B->point() - n1A->point();
      n2 = n2B->point() - n2A->point();

      double theta1 = atan2(n1[1],n1[0]);
      double theta2 = atan2(n2[1],n2[0]);

      Vector2D r1;
      Vector2D r2;
      Vector2D R;
      r1 = .5*(n1A->point() + n1B->point());
      r2 = .5*(n2A->point() + n2B->point());
      R = r2 - r1;
      _box->mapDistance(R);

      // rotate coordinate system so that n1 points along x //
      Vector2D Rnew;
      Rnew[0] = cos(theta1)*R[0]+sin(theta1)*R[1];
      Rnew[1] = cos(theta1)*R[1]-sin(theta1)*R[0];

      Vector2D n2new;
      n2new[0] = cos(theta1)*n2[0] + sin(theta1)*n2[1];
      n2new[1] = cos(theta1)*n2[1] - sin(theta1)*n2[0];
      double theta = atan2(n2new[1],n2new[0]);
      if(theta < 0.0) theta += M_PI;

      R = Rnew;

      std::cout << "\nTable check:\nR = (" << R[0] << ", " << R[1] << "), theta = " << theta << "." << std::endl;

      int row = (int)((Rnew[0]-_Rxlims[0])/_stpRx);
      int column = (int)((Rnew[1]-_Rylims[0])/_stpRy);
      int depth = (int)((theta-_thetalims[0])/_stpTheta);
      
      double Rxd = (Rnew[0]-(_Rxlims[0]+row*_stpRx))/_stpRx;
      double Ryd = (Rnew[1]-(_Rylims[0]+column*_stpRy))/_stpRy;
      double thetad = (theta-(_thetalims[0]+depth*_stpTheta))/_stpTheta;

      double t1,t2,t3,t4;
      double Ry1,Ry2;

      // do trilinear interpolation for energy //
    
      t1 = thetad*_energyTable(row,column,depth+1)+(1.0-thetad)*_energyTable(row,column,depth);
      t2 = thetad*_energyTable(row,column+1,depth+1)+(1.0-thetad)*_energyTable(row,column+1,depth);
      t3 = thetad*_energyTable(row+1,column,depth+1)+(1.0-thetad)*_energyTable(row+1,column,depth);
      t4 = thetad*_energyTable(row+1,column+1,depth+1)+(1.0-thetad)*_energyTable(row+1,column+1,depth);
      
      Ry1 = Ryd*t2+(1.0-Ryd)*t1;
      Ry2 = Ryd*t4+(1.0-Ryd)*t3;
      
      energy = Rxd*Ry2 + (1.0-Rxd)*Ry1;
      
      std::cout << "Energy = " << energy << "; actual energy = " << integrateEnergy(R,theta) << std::endl;
      
      // do trilinear interpolation for force //
      Vector2D force;
      t1 = thetad*_forceXTable(row,column,depth+1)+(1.0-thetad)*_forceXTable(row,column,depth);
      t2 = thetad*_forceXTable(row,column+1,depth+1)+(1.0-thetad)*_forceXTable(row,column+1,depth);
      t3 = thetad*_forceXTable(row+1,column,depth+1)+(1.0-thetad)*_forceXTable(row+1,column,depth);
      t4 = thetad*_forceXTable(row+1,column+1,depth+1)+(1.0-thetad)*_forceXTable(row+1,column+1,depth);
      
      Ry1 = Ryd*t2+(1.0-Ryd)*t1;
      Ry2 = Ryd*t4+(1.0-Ryd)*t3;
      
      force[0] = Rxd*Ry2 + (1.0-Rxd)*Ry1;
      
      t1 = thetad*_forceYTable(row,column,depth+1)+(1.0-thetad)*_forceYTable(row,column,depth);
      t2 = thetad*_forceYTable(row,column+1,depth+1)+(1.0-thetad)*_forceYTable(row,column+1,depth);
      t3 = thetad*_forceYTable(row+1,column,depth+1)+(1.0-thetad)*_forceYTable(row+1,column,depth);
      t4 = thetad*_forceYTable(row+1,column+1,depth+1)+(1.0-thetad)*_forceYTable(row+1,column+1,depth);
      
      Ry1 = Ryd*t2+(1.0-Ryd)*t1;
      Ry2 = Ryd*t4+(1.0-Ryd)*t3;
      
      force[1] = Rxd*Ry2 + (1.0-Rxd)*Ry1;
      
      Vector2D actforce;
      actforce = integrateForce(R,theta);

      std::cout << "Force = (" << force[0] << ", " << force[1] << "); actual force = (" << actforce[0] << ", " << actforce[1] << ")." << std::endl;
      
      // do trilinear interpolation for torque //
      t1 = thetad*_torqueTable(row,column,depth+1)+(1.0-thetad)*_torqueTable(row,column,depth);
      t2 = thetad*_torqueTable(row,column+1,depth+1)+(1.0-thetad)*_torqueTable(row,column+1,depth);
      t3 = thetad*_torqueTable(row+1,column,depth+1)+(1.0-thetad)*_torqueTable(row+1,column,depth);
      t4 = thetad*_torqueTable(row+1,column+1,depth+1)+(1.0-thetad)*_torqueTable(row+1,column+1,depth);
      
      Ry1 = Ryd*t2+(1.0-Ryd)*t1;
      Ry2 = Ryd*t4+(1.0-Ryd)*t3;
      
      double torque = Rxd*Ry2 + (1.0-Rxd)*Ry1;
      
      std::cout << "Torque = " << torque << "; actual torque = " << integrateTorque(R,theta) << ".\n" << std::endl;
    }
      
    double compute(DefNodeContainer & fil1, DefNodeContainer & fil2, bool f0, bool f1, bool f2) {
      double energy = 0.0;
      
      DefNode* n1A = fil1[0];
      DefNode* n1B = fil1[1];
      DefNode* n2A = fil2[0];
      DefNode* n2B = fil2[1];
      
      Vector2D n1;
      Vector2D n2;
      
      n1 = n1B->point() - n1A->point();
      n2 = n2B->point() - n2A->point();
      
      double theta1 = atan2(n1[1],n1[0]);
      double theta2 = atan2(n2[1],n2[0]);
      
      Vector2D r1;
      Vector2D r2;
      Vector2D R;
      r1 = .5*(n1A->point() + n1B->point());
      r2 = .5*(n2A->point() + n2B->point());
      R = r2 - r1;
      _box->mapDistance(R);
      
      // rotate coordinate system so that n1 points along x //
      Vector2D Rnew;
      Rnew[0] = cos(theta1)*R[0]+sin(theta1)*R[1];
      Rnew[1] = cos(theta1)*R[1]-sin(theta1)*R[0];

      Vector2D n2new;
      n2new[0] = cos(theta1)*n2[0] + sin(theta1)*n2[1];
      n2new[1] = cos(theta1)*n2[1] - sin(theta1)*n2[0];
      double theta = atan2(n2new[1],n2new[0]);
      if(theta < 0.0) theta += M_PI;
      
      int row = (int)((Rnew[0]-_Rxlims[0])/_stpRx);
      int column = (int)((Rnew[1]-_Rylims[0])/_stpRy);
      int depth = (int)((theta-_thetalims[0])/_stpTheta);
      
      double Rxd = (Rnew[0]-(_Rxlims[0]+row*_stpRx))/_stpRx;
      double Ryd = (Rnew[1]-(_Rylims[0]+column*_stpRy))/_stpRy;
      double thetad = (theta-(_thetalims[0]+depth*_stpTheta))/_stpTheta;
      
      double t1,t2,t3,t4;
      double Ry1,Ry2;
      
      // do trilinear interpolation for energy //
      if(f0) {
	t1 = thetad*_energyTable(row,column,depth+1)+(1.0-thetad)*_energyTable(row,column,depth);
	t2 = thetad*_energyTable(row,column+1,depth+1)+(1.0-thetad)*_energyTable(row,column+1,depth);
	t3 = thetad*_energyTable(row+1,column,depth+1)+(1.0-thetad)*_energyTable(row+1,column,depth);
	t4 = thetad*_energyTable(row+1,column+1,depth+1)+(1.0-thetad)*_energyTable(row+1,column+1,depth);
	
	Ry1 = Ryd*t2+(1.0-Ryd)*t1;
	Ry2 = Ryd*t4+(1.0-Ryd)*t3;
	
	energy = Rxd*Ry2 + (1.0-Rxd)*Ry1;
      }

      // do trilinear interpolation for force //
      if(f1) {
	Vector2D force;
	t1 = thetad*_forceXTable(row,column,depth+1)+(1.0-thetad)*_forceXTable(row,column,depth);
	t2 = thetad*_forceXTable(row,column+1,depth+1)+(1.0-thetad)*_forceXTable(row,column+1,depth);
	t3 = thetad*_forceXTable(row+1,column,depth+1)+(1.0-thetad)*_forceXTable(row+1,column,depth);
	t4 = thetad*_forceXTable(row+1,column+1,depth+1)+(1.0-thetad)*_forceXTable(row+1,column+1,depth);
	
	Ry1 = Ryd*t2+(1.0-Ryd)*t1;
	Ry2 = Ryd*t4+(1.0-Ryd)*t3;
	
	force[0] = Rxd*Ry2 + (1.0-Rxd)*Ry1;
	
	t1 = thetad*_forceYTable(row,column,depth+1)+(1.0-thetad)*_forceYTable(row,column,depth);
	t2 = thetad*_forceYTable(row,column+1,depth+1)+(1.0-thetad)*_forceYTable(row,column+1,depth);
	t3 = thetad*_forceYTable(row+1,column,depth+1)+(1.0-thetad)*_forceYTable(row+1,column,depth);
	t4 = thetad*_forceYTable(row+1,column+1,depth+1)+(1.0-thetad)*_forceYTable(row+1,column+1,depth);
	
	Ry1 = Ryd*t2+(1.0-Ryd)*t1;
	Ry2 = Ryd*t4+(1.0-Ryd)*t3;
	
	force[1] = Rxd*Ry2 + (1.0-Rxd)*Ry1;
	// rotate force back to original coordinate system //
	Vector2D tmpForce(force);
	force[0] = cos(theta1)*tmpForce[0] - sin(theta1)*tmpForce[1];
	force[1] = cos(theta1)*tmpForce[1] + sin(theta1)*tmpForce[0];
	
	// do trilinear interpolation for torque //
	t1 = thetad*_torqueTable(row,column,depth+1)+(1.0-thetad)*_torqueTable(row,column,depth);
	t2 = thetad*_torqueTable(row,column+1,depth+1)+(1.0-thetad)*_torqueTable(row,column+1,depth);
	t3 = thetad*_torqueTable(row+1,column,depth+1)+(1.0-thetad)*_torqueTable(row+1,column,depth);
	t4 = thetad*_torqueTable(row+1,column+1,depth+1)+(1.0-thetad)*_torqueTable(row+1,column+1,depth);
	
	Ry1 = Ryd*t2+(1.0-Ryd)*t1;
	Ry2 = Ryd*t4+(1.0-Ryd)*t3;
	
	double torque = Rxd*Ry2 + (1.0-Rxd)*Ry1;
	
	// add forces to filament 1 //
	Vector2D f1A;
	f1A = force/2.0;
	f1A[0] += (torque/_elemlength)*sin(theta1);
	f1A[1] -= (torque/_elemlength)*cos(theta1);
	n1A->updateForce(f1A);
	
	Vector2D f1B;
	f1B = force/2.0;
	f1B[0] -= (torque/_elemlength)*sin(theta1);
	f1B[1] += (torque/_elemlength)*cos(theta1);
	n1B->updateForce(f1B);	
	
	// add forces to filament 2 //
	Vector2D f2A;
	f2A = -force/2.0;
	f2A[0] -= (torque/_elemlength)*sin(theta2);
	f2A[1] += (torque/_elemlength)*cos(theta2);
	n2A->updateForce(f2A);

	Vector2D f2B;
	f2B = -force/2.0;
	f2B[0] += (torque/_elemlength)*sin(theta2);
	f2B[1] -= (torque/_elemlength)*cos(theta2);
	n2B->updateForce(f2B);

      }
      
      return energy;
    }

    double potentialenergy(Vector2D & p1, Vector2D & p2) {
      if(_ft == gaussian) {
	return _prefact*exp(-sqr(norm2(p1-p2)/_scale));
      }
      else if(_ft == expon) {
	return _prefact*exp(-norm2(p1-p2)/_scale);
      }
    }

    double integrateEnergy(Vector2D & R, double theta) {
      Vector2D n1(1.0,0.0);
      Vector2D n2(cos(theta),sin(theta));
      double energy = 0.0;
      double chargeSep = _elemlength/(_nChrgs-1);
      Vector2D x1;
      Vector2D x2;
      for(int i=0; i<_nChrgs; i++) {
	x1 = ((-_elemlength/2.0)+chargeSep*i)*n1;
	for(int j=0; j<_nChrgs; j++) {
	  x2 = R + ((-_elemlength/2.0)+chargeSep*j)*n2;
	  energy += potentialenergy(x1,x2);
	}
      }
      return energy;
    }

    Vector2D integrateForce(Vector2D & R, double theta) {
      Vector2D n1(1.0,0.0);
      Vector2D n2(cos(theta),sin(theta));
      Vector2D force(0.0,0.0);
      double chargeSep = _elemlength/(_nChrgs-1);
      Vector2D x1;
      Vector2D x2;
      for(int i=0; i<_nChrgs; i++) {
	x1 = ((-_elemlength/2.0)+chargeSep*i)*n1;
	for(int j=0; j<_nChrgs; j++) {
	  x2 = R + ((-_elemlength/2.0)+chargeSep*j)*n2;
	  if(_ft==gaussian) {
	    force += (2.0/sqr(_scale))*potentialenergy(x1,x2)*(x2-x1);
	  }
	  else if(_ft==expon) {
	    Vector2D diff;
	    diff = x2 - x1;
	    double diffnorm = norm2(diff);
	    diff /= diffnorm;
	    force += (1.0/_scale)*potentialenergy(x1,x2)*diff;
	  }
	}
      }
      return force;
    }
    
    double integrateTorque(Vector2D & R, double theta) {
      Vector2D n1(1.0,0.0);
      Vector2D n2(cos(theta),sin(theta));
      double torque = 0.0;
      double chargeSep = _elemlength/(_nChrgs-1);
      Vector2D x1;
      Vector2D x2;
      double s1;
      double s2;
      for(int i=0; i<_nChrgs; i++) {
	s1 = ((-_elemlength/2.0)+chargeSep*i);
	x1 = s1*n1;
	for(int j=0; j<_nChrgs; j++) {
	  s2 = ((-_elemlength/2.0)+chargeSep*i);
	  x2 = R + s2*n2;
	  if(_ft==gaussian) {
	    torque += (2.0/sqr(_scale))*potentialenergy(x1,x2)*(s1*s2*sin(theta)+s2*(R[1]*cos(theta)-R[0]*sin(theta)));
	  }
	  else if(_ft==expon) {
	    Vector2D diff;
	    diff = x2 - x1;
	    double diffnorm = norm2(diff);
	    torque += (1.0/(_scale*diffnorm))*potentialenergy(x1,x2)*(s1*s2*sin(theta)+s2*(R[1]*cos(theta)-R[0]*sin(theta)));
	  }
	}
      }
      return torque;

    }

    void makeTables() {
      _energyTable.resize(100,100,100);
      _forceXTable.resize(100,100,100);
      _forceYTable.resize(100,100,100);
      _torqueTable.resize(100,100,100);
      
      for(int i=0; i<100; i++) {
	for(int j=0; j<100; j++) {
	  for(int k=0; k<100; k++) {
	    Vector2D R;
	    R[0] = _Rxlims[0] + i*_stpRx;
	    R[1] = _Rylims[0] + j*_stpRy;
	    double theta = _thetalims[0] + k*_stpTheta;
	    _energyTable(i,j,k) = integrateEnergy(R,theta);
	    Vector2D tmpForce;
	    tmpForce = integrateForce(R,theta);
	    _forceXTable(i,j,k) = tmpForce[0];
	    _forceYTable(i,j,k) = tmpForce[1];
	    _torqueTable(i,j,k) = integrateTorque(R,theta);
	  }
	}
	if((i+1)%10==0) std::cout << "TwoBodyPotential: Making lookup tables: " << (i+1) << " percent done." << std::endl;
      }
    }

  private:

    PeriodicBox * _box;

    FunctionType _ft;
    double _prefact;
    double _elemlength;
    double _scale;
    int _nChrgs;

    Vector2D _Rxlims;
    double _stpRx;
    Vector2D _Rylims;
    double _stpRy;
    Vector2D _thetalims;
    double _stpTheta;
    Array3D _energyTable;
    Array3D _forceXTable;
    Array3D _forceYTable;
    Array3D _torqueTable;

    int _maxSamps;

  };
}

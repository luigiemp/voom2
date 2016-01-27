// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__AngleSpring_h__)
#define __AngleSpring_h__

#include "Node.h"
#include "Element.h"

namespace voom
{

  template<int N>
  class AngleSpring : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef DeformationNode<N> Node_t;    

    AngleSpring(Node_t * nodeA, Node_t * nodeB, Node_t * nodeC, double k) 
      : _nodeA(nodeA), _nodeB(nodeB), _nodeC(nodeC), _k(k) { 
      _baseNodes.push_back(nodeA);
      _baseNodes.push_back(nodeB);
      _baseNodes.push_back(nodeC);
    }

    // angle theta = /_ABC
    // 
    //  A       C
    //    \   /
    //      B
    // E = 0.5*k*(theta^2) \approx k*(1-cos(theta))
    //
    void compute(bool f0, bool f1, bool f2) {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      const VectorND & xC = _nodeC->point();

      VectorND tBA;
      tBA = xA-xB;
      double LBA = norm2(tBA);
      tBA /= LBA;
      
      VectorND tBC;
      tBC = xC-xB;
      double LBC = norm2(tBC);
      tBC /= LBC;

      double cosABC = dot(tBA,tBC);

      if(f0) {
	_energy = _k*(1.0+cosABC);
      }

      if(f1) {	
	VectorND fA;
	fA = -cosABC * tBA;
	fA += tBC;
	fA *=  _k/LBA;

	VectorND fC;
	fC = -cosABC * tBC;
	fC += tBA;
	fC *=  _k/LBC;

	for(int i=0; i<N; i++) {
	  _nodeA->addForce( i, fA(i) );
	  _nodeB->addForce( i, -fA(i)-fC(i) );
	  _nodeC->addForce( i, fC(i)) ;
	}
      }

      return;
    }
    
    double stiffness() const {return _k;}

    void setStiffness(double k) { _k = k; }

    double meanSegmentExtension() {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      const VectorND & xC = _nodeC->point();
      double dAB = norm2(xA-xB);
      double dBC = norm2(xC-xB);
      return 0.5*(dAB+dBC);
    }
    
    double meanSegmentLength() {
      const VectorND & xA = _nodeA->position();
      const VectorND & xB = _nodeB->position();
      const VectorND & xC = _nodeC->position();
      double dAB = norm2(xA-xB);
      double dBC = norm2(xC-xB);
      return 0.5*(dAB+dBC);
    }

    double getAngleABC() {
      VectorND tAB;
      tAB = _nodeB->point() - _nodeA->point();
      double lAB = norm2(tAB);
      tAB /= lAB;
      VectorND tBC;
      tBC = _nodeC->point() - _nodeB->point();
      double lBC = norm2(tBC);
      tBC /= lBC;
      double angAB = atan2(tAB[1],tAB[0]);
      // rotate tBC //
      VectorND temptBC(tBC);
      tBC[0] = cos(angAB)*temptBC[0]+sin(angAB)*temptBC[1];
      tBC[1] = cos(angAB)*temptBC[1]-sin(angAB)*temptBC[0];
      double angABC = atan2(tBC[1],tBC[0]);
	
      return angABC;
    }

  private:
    
    Node_t * _nodeA;
    Node_t * _nodeB;
    Node_t * _nodeC;
    double _k;
    
  };
};

#endif // __AngleSpring_h__

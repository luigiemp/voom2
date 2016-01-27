// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Andrew R. Missel
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__PinchForce_h__)
#define __PinchForce_h__

#include "Node.h"
#include "Element.h"
#include "VoomMath.h"
#include <iostream>
#include "PeriodicBox.h"

using namespace tvmet;
using namespace std;
using namespace voom;

namespace voom
{

  template<int N>
  class PinchForce : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef BrownianNode<N> Node_t;

    PinchForce(Node_t * node1, Node_t * node2, double f0, PeriodicBox * box) 
      : _node1(node1), _node2(node2), _f0(f0), _box(box), _state(true) {
      
      VectorND p1 = _node1->position();
      VectorND p2 = _node2->position();

      VectorND diff;
      diff = p1-p2;
      _box->mapDistance(diff);
      _d0 = norm2(diff);

    }

    void setBox(PeriodicBox * pb) {
      _box = pb;
    }

    void compute(bool f0, bool f1, bool f2) {

      _energy = 0.0;

      if(_state) {
	const VectorND & x1 = _node1->point();
	const VectorND & x2 = _node2->point();
	//Periodic BC
	
	VectorND dx(0.0);
	dx = x2 - x1;
	_box->mapDistance(dx);
	double d = norm2(dx);
	
	if(f0) {
	  _energy = _f0*(d-_d0);
	}
	
	if(f1) {
	  for(int i=0; i<N; i++) {
	    double f = (dx(i)/d)*_f0;
	    _node1->addForce(i, -f);
	    _node2->addForce(i, f);
	  }
	  
	}
      }
      
      return;
    }
    
    double pinchF() const {return _f0;}
    
    void turnOn() {
      _state = true;
    }
    
    void turnOff() {
      _state = false;
    }

    void setPinchF(double f0) { _f0 = f0; }

    std::pair<Node_t *,Node_t *> getNodes() {
      std::pair<Node_t *, Node_t *> nds = pair<Node_t *, Node_t *>(_node1,_node2);
      return nds;
    }

//     Node_tContainer getNodes(){
//       Node_tContainer nodes;
//       nodes.push_back(_node1);
//       nodes.push_back(_node2);
//       return nodes;
//     }

//     Node_t * getNode(int a) {
//       if(a==0) {
// 	return _node1;
//       } else if(a==1) {
// 	return _node2;
//       } else {
// 	return 0;
//       }
//     }
  private:
    
    Node_t * _node1;
    Node_t * _node2;
    bool _state;
    double _f0;
    double _d0;
    PeriodicBox * _box;
    
  };
};

#endif // __PinchForce_h__

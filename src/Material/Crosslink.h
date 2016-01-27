// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Crosslink_h__)
#define __Crosslink_h__

#include "Node.h"
#include "Element.h"
#include "VoomMath.h"
#include "Spring.h"
#include <iostream>
#include "PeriodicBox.h"

using namespace tvmet;
using namespace std;
using namespace voom;

namespace voom
{

  template<int N>
  class Crosslink : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef BrownianNode<N> Node_t;
    typedef std::vector< Node_t*> Node_tContainer;
    typedef typename Node_tContainer::iterator Node_tIterator;

    Crosslink(double k, PeriodicBox * box) : _k(k), _d0(0.0), _box(box) {
      int a = 0;
      int id = a;
      BrownianNode<2>::Point X;
      X = 0.0 , 0.0;
      //X = x1A*(1-_epsi1)+x1B*_epsi1;
      NodeBase::DofIndexMap idx(2);
      idx[0]=2*a; idx[1]=2*a+1;
      _node1 = new Node_t(id,idx,X,X);
      _node1->setId(id);

      a = a+1;
      id = a;
      X = 0.0, 0.0;
      //X = x2A*(1-_epsi2)+x2B*_epsi2;
      idx[0]=2*a; idx[1]=2*a+1;
      _node2 = new Node_t(id,idx,X,X);
      _node2->setId(id);
      
    }

    Crosslink(Node_t * node1A, Node_t * node1B, Node_t * node2A, Node_t * node2B, double epsi1, double epsi2, double k, PeriodicBox * box, bool relaxed) 
      : _node1A(node1A), _node1B(node1B), _node2A(node2A), _node2B(node2B), _epsi1(epsi1), _epsi2(epsi2), _k(k), _box(box) { 
      //changed from const & to variables
      VectorND x1A = _node1A->point();
      VectorND x1B = _node1B->point();
      VectorND x2A = _node2A->point();
      VectorND x2B = _node2B->point();

      int a = 0;
      int id = a;
      BrownianNode<2>::Point X1;
      //X = 0.0 , 0.0;
      X1 = x1A*(1.0-_epsi1)+x1B*_epsi1;
      NodeBase::DofIndexMap idx(2);
      idx[0]=2*a; idx[1]=2*a+1;
      _node1 = new Node_t(id,idx,X1,X1);
      _node1->setId(id);

      a = a+1;
      id = a;
      BrownianNode<2>::Point X2;
      X2 = x2A*(1.0-_epsi2)+x2B*_epsi2;
      idx[0]=2*a; idx[1]=2*a+1;
      _node2 = new Node_t(id,idx,X2,X2);
      _node2->setId(id);
      
      if(relaxed) {
      	BrownianNode<2>::Point Xsep;
	Xsep = X2 - X1;
	_box->mapDistance(Xsep);
	_d0 = norm2(Xsep);	
      }
      else {
      	_d0=0.0;
      }
    }

    void setBox(PeriodicBox * pb) {
      _box = pb;
    }

    void setPosition(Node_t * node1A, Node_t * node1B, Node_t * node2A, Node_t * node2B, double epsi1, double epsi2) {
      _node1A = node1A; 
      _node1B = node1B; 
      _node2A = node2A;
      _node2B = node2B; 
      _epsi1 = epsi1; 
      _epsi2 = epsi2;

      VectorND x1A = _node1A->point();
      VectorND x1B = _node1B->point();
      VectorND x2A = _node2A->point();
      VectorND x2B = _node2B->point();

      for(int i=0; i<N; i++) {	  
        double p1 = x1A(i)*(1.0-_epsi1)+x1B(i)*_epsi1;
        double p2 = x2A(i)*(1.0-_epsi2)+x2B(i)*_epsi2;
        _node1->setPoint(i, p1);
	_node2->setPoint(i, p2);
      }
 
    }

    void compute(bool f0, bool f1, bool f2) {

      VectorND x1A = _node1A->point();
      VectorND x1B = _node1B->point();
      VectorND x2A = _node2A->point();
      VectorND x2B = _node2B->point();

      //Periodic BC
      
      for(int i=0; i<N; i++) {	  

        double p1 = x1A(i)*(1.0-_epsi1)+x1B(i)*_epsi1;
        double p2 = x2A(i)*(1.0-_epsi2)+x2B(i)*_epsi2;
        _node1->setPoint(i, p1);
	_node2->setPoint(i, p2);
      }
 
      const VectorND & x1 = _node1->point();
      const VectorND & x2 = _node2->point();
      VectorND dx(0.0);
      dx = x1-x2;
      _box->mapDistance(dx);
      double d  = norm2(dx);
      
      if(f0) {
	_energy = 0.5*_k*sqr(d-_d0);
      }
      
      if(f1) {	
	for(int i=0; i<N; i++) {
	  double f = -_k*(d-_d0)*(dx(i))/d;
	  _node1A->addForce(i, -f*(1.0-_epsi1));
          _node1B->addForce(i, -f*_epsi1);
          _node2A->addForce(i, f*(1.0-_epsi2));
	  _node2B->addForce(i, f*_epsi2);
	}
      }
      return;
    }
    
    double stiffness() const {return _k;}

    void setStiffness(double k) { _k = k; }

    Node_tContainer getNodes(){
      Node_tContainer nodes;
      nodes.push_back(_node1);
      nodes.push_back(_node2);
      return nodes;
    }

    Node_t * getNode(int a) {
      if(a==0) {
	return _node1;
      } else if(a==1) {
	return _node2;
      } else {
	return 0;
      }
    }
  private:
    
    Node_t * _node1;
    Node_t * _node2;
    Node_t * _node1A;
    Node_t * _node1B;
    Node_t * _node2A;
    Node_t * _node2B;
    double _epsi1;
    double _epsi2;
    double _k;
    double _d0;
    PeriodicBox * _box;
    
  };
};

#endif // __Crosslink_h__

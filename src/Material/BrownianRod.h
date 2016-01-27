
#if !defined(__BrownianRod_h__)
#define __BrownianRod_h__

#include<random/normal.h>
#include "Node.h"
#include "Element.h"
#include "VoomMath.h"

namespace voom
{

  template<int N>
  class BrownianRod : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef BrownianNode<N> Node_t;    
    typedef std::vector< Node_t *> NodeContainer;

    BrownianRod(Node_t * nodeA, Node_t * nodeB, double viscosity, double kT, double dt) 
      : _nodeA(nodeA), _nodeB(nodeB), _viscosity(viscosity), _kT(kT), _timeStep(dt) { 
      const VectorND & XA = _nodeA->position();
      const VectorND & XB = _nodeB->position();
      _d0 = norm2(XB-XA);      
      _baseNodes.push_back(_nodeA);
      _baseNodes.push_back(_nodeB);

      NodeContainer nodes;
      nodes.push_back(_nodeA);
      nodes.push_back(_nodeB);
      _nodes = nodes;
      //_nodes.push_back(_nodeA);
      //_nodes.push_back(_nodeB);


      // seed random number generator
      _rng.seed((unsigned int)time(0));

    }

    BrownianRod(Node_t * nodeA, Node_t * nodeB, double viscosity, double kT, double dt, double d) 
      : _nodeA(nodeA), _nodeB(nodeB), _viscosity(viscosity), _kT(kT), _timeStep(dt), _d0(d) {
      _baseNodes.push_back(_nodeA);
      _baseNodes.push_back(_nodeB);

      NodeContainer nodes;
      nodes.push_back(_nodeA);
      nodes.push_back(_nodeB);
      _nodes = nodes;
      //_nodes.push_back(_nodeA);
      //_nodes.push_back(_nodeB);
      
      // seed random number generator
      _rng.seed((unsigned int)time(0));

    }

    void compute(bool f0, bool f1, bool f2) {
      _energy = 0.0;

      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double l  = norm2(xB-xA);
      
      // compute drag
      _Dpara = 2.0*M_PI*l*_viscosity;
      _Dperp = 2.0*_Dpara;
      
      //Node_t::Matrix<2> M(0.0);
      Tensor2D D(0.0);
      for(int i=0; i<N; i++) D(i,i) = _Dperp;
      
      VectorND epara(0.0);
      epara = (xB-xA)/l;
      
      // HW for MO: compute unit vector to bond
      VectorND eperp(0.0);
      // If N=2
      eperp(0) = -epara(1);
      eperp(1) =  epara(0);
      
      for(int i=0; i<N; i++) {
	for(int j=0; j<N; j++) {
	  D(i,j) += (_Dpara-_Dperp)*epara(i)*epara(j);
	}
      }
      
      if(f2) {
	double fpara = sqrt( 2.0*_kT*_Dpara/_timeStep ) * _rng.random();
	double fperp = sqrt( 2.0*_kT*_Dperp/_timeStep ) * _rng.random();
      /*       double fpara = sqrt( 2.0*_kT*_timeStep*_Mpara ) * _rng.random(); */
      /*       double fperp = sqrt( 2.0*_kT*_timeStep*_Mperp ) * _rng.random(); */
      
      //Mo
      //ofstream Brown("brown.dat", ofstream::app);
      //Mo
	for(int i=0; i<N; i++) {	  
	
	  // add parallel and perpendicular Brownian forces
	  double f = 0.5* ( fpara*epara(i) + fperp*eperp(i) );
	  _nodeA->addForce(i, f);
	  _nodeB->addForce(i, f);
	  
	  //Mo
	  //Brown << f*2.0 << ' ';
	}
      }
      //Mo
      //Brown << std::endl;
      //Brown.close();
      
      
      D *= 0.5;
      _nodeA->addDrag(D);
      _nodeB->addDrag(D); 
      
      return;
    }
    
    double viscosity() const {return _viscosity;}

    void setViscosity(double v) { _viscosity = v; }

    NodeContainer getNodes() { return _nodes;}

  private:
    
    Node_t * _nodeA;
    Node_t * _nodeB;
    NodeContainer _nodes;
    double _viscosity;
    double _Dpara;
    double _Dperp;
    double _d0;
    double _timeStep;
    double _kT;

    ranlib::NormalUnit<double> _rng;
 
   
  };
};

#endif // __BrownianRod_h__

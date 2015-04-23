//-*-C++-*_
/*!
  \file LoopShellShape.h
  
  \brief Shape function implementation for Loop Shell Element
*/

#if !defined(__LoopShellShape_h__)
#define __LoopShellShape_h__

#include "Shape.h"

namespace voom{

  class LoopShellShape: public Shape {

  public:
    typedef Eigen::Matrix<double,12,Dynamic> SubdivisionMatrix;
    typedef Eigen::Vector3i CornerValences; //Valences of nodes at three corners

    //! LoopShellShape constructor fills in N, DN and DDN
    LoopShellShape(int nodes, CornerValences & V, const VectorXd & Point) {
      _nodes = nodes;
      _coords = Point;
      _Valences = V;
      _N.resize(nodes, 0.0); //12 shape functions
      _DN.resize(nodes, Vector2d::Zero()); //partial derivatives in 'v' and 'w' directions
      _DDN.resize(nodes, Matrix2d::Zero()); //partial deriatives in ['vv', 'vw'; 'wv' 'ww'] directions
      update(Point);
    }

   
    //! Update recomputes N, DN and DDN at a Point assuming it's in a regular patch
    void update(const VectorXd & Point);

    //! Get number of shape functions
    uint getShapeFunctionNum() {return _nodes; };

    //! getN returns shape function values at a given node a.
    Real getN(const uint a) {
      return _N[a];
    };
    
    //! GetDN returns shape function derivatives at node a, in direction 'v'/'w'.
    Real getDN(const uint a, const uint i) {
      return _DN[a](i);
    };
   
    //! GetDDN returns shape function second derivatives at node a
    //! DDN is a vector of Matrix2d 
    Real getDDN(const uint a, const uint i, const uint j) {
      return _DDN[a](i,j);
    };
	
    
  private:
    vector<Real> _N;        //require: vector dim = _nodes
    vector<Vector2d > _DN;  //---------- do ------------
    vector<Matrix2d> _DDN;  //---------- do ------------
    Vector2d _coords;
    CornerValences _Valences;

    int _nodes; //Regular patch: 12, Irregular patch: variable

    void _computeFunctions( const SubdivisionMatrix & S, bool needSubdivide);

    void _computeDerivatives( const SubdivisionMatrix & S, bool needSubdivide);

    void _computeSecondDerivatives( const SubdivisionMatrix & S, bool needSubdivide);

    //! compute subdivision Matrix
    void _computeSubdivisionMatrix(SubdivisionMatrix & S, 
				    bool needSubdivide );
    
    void _convertParaCoords();

    void _initialize(const int nodes, const VectorXd & paraCoords);
    
    void _initialize(const int nodes) {
      Vector2d paraCoords(1.0/3.0, 1.0/3.0);
      _initialize(nodes, paraCoords);
    };

  };
} // namespace voom

#endif

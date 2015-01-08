#include "HexShape.h"

namespace voom {

//! Update recomputes N and DN at a new Point
  void HexShape::update(const VectorXd & Point)
  {
    // Shape functions
    _N[0] = 0.125*(1-Point(0))*(1-Point(1))*(1+Point(2));
    _N[1] = 0.125*(1+Point(0))*(1-Point(1))*(1+Point(2));
    _N[2] = 0.125*(1+Point(0))*(1-Point(1))*(1-Point(2));
    _N[3] = 0.125*(1-Point(0))*(1-Point(1))*(1-Point(2));
    _N[4] = 0.125*(1-Point(0))*(1+Point(1))*(1+Point(2));
    _N[5] = 0.125*(1+Point(0))*(1+Point(1))*(1+Point(2));  
    _N[6] = 0.125*(1+Point(0))*(1+Point(1))*(1-Point(2));
    _N[7] = 0.125*(1-Point(0))*(1+Point(1))*(1-Point(2));

    // Shape functions derivatives
    _DN[0](0) = -0.125*(1-Point(1))*(1+Point(2));
    _DN[0](1) = -0.125*(1-Point(0))*(1+Point(2));
    _DN[0](2) =  0.125*(1-Point(0))*(1-Point(1));
    
    _DN[1](0) =  0.125*(1-Point(1))*(1+Point(2));
    _DN[1](1) = -0.125*(1+Point(0))*(1+Point(2));
    _DN[1](2) =  0.125*(1+Point(0))*(1-Point(1));
    
    _DN[2](0) =  0.125*(1-Point(1))*(1-Point(2)); 
    _DN[2](1) = -0.125*(1+Point(0))*(1-Point(2));
    _DN[2](2) = -0.125*(1+Point(0))*(1-Point(1));
    
    _DN[3](0) = -0.125*(1-Point(1))*(1-Point(2));
    _DN[3](1) = -0.125*(1-Point(0))*(1-Point(2));
    _DN[3](2) = -0.125*(1-Point(0))*(1-Point(1));
    
    _DN[4](0) = -0.125*(1+Point(1))*(1+Point(2));
    _DN[4](1) =  0.125*(1-Point(0))*(1+Point(2));
    _DN[4](2) =  0.125*(1-Point(0))*(1+Point(1));
    
    _DN[5](0) =  0.125*(1+Point(1))*(1+Point(2));
    _DN[5](1) =  0.125*(1+Point(0))*(1+Point(2));
    _DN[5](2) =  0.125*(1+Point(0))*(1+Point(1));   
    
    _DN[6](0) =  0.125*(1+Point(1))*(1-Point(2));  
    _DN[6](1) =  0.125*(1+Point(0))*(1-Point(2));
    _DN[6](2) = -0.125*(1+Point(0))*(1+Point(1));
    
    _DN[7](0) = -0.125*(1+Point(1))*(1-Point(2));
    _DN[7](1) =  0.125*(1-Point(0))*(1-Point(2));
    _DN[7](2) = -0.125*(1-Point(0))*(1+Point(1)); 
    
  } // update
						       
} // namespace voom

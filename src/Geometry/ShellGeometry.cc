// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.5  2005/05/23 17:43:20  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file ShellGeometry.cc

  \brief ShellGeometry is a structure which computes and contains all
  of the geometric quantities relevant to calculations on a curved surface.

*/

#include "ShellGeometry.h"

namespace voom
{
  void ShellGeometry::_updateMetrics()
  {
    //
    //  Calculate the covariant and contravariant metric tensors
    //
    //  Tensor2D m_MetricTensorInverse; a^(0,0), a^(0,1), a^(1,0) and a^(1,1)
    //  Tensor2D m_MetricTensor;        a_(0,0), a_(0,1), a_(1,0) and a_(1,1)
    _metricTensor(0,0) =  _a[0].dot(_a[0]);
    _metricTensor(0,1) =  _a[0].dot(_a[1]);
    _metricTensor(1,1) =  _a[1].dot(_a[1]);
    _metricTensor(1,0) = _metricTensor(0,1);

    _metric = _metricTensor(0,0) * _metricTensor(1,1)
      - _metricTensor(0,1) * _metricTensor(1,0);

    //assert(_metric > 0);
    if( !( _metric > 0.0 ) ) {
      std::cerr << "a0 = " << _a[0] << std::endl
		<< "a1 = " << _a[1] << std::endl
		<< "a_{ij} = " << _metricTensor << std::endl
		<< "det(a_{ij}) = "<< _metric << std::endl
		<< "Exiting." << std::endl;
      exit(0);
    }

    _metricTensorInverse(0,0) =   _metricTensor(1,1) / _metric;
    _metricTensorInverse(0,1) = - _metricTensor(0,1) / _metric;
    _metricTensorInverse(1,1) =   _metricTensor(0,0) / _metric;
    _metricTensorInverse(1,0) =   _metricTensorInverse(0,1);

    _metric = sqrt(_metric);
  }


  void ShellGeometry::_updateDual()
  {
    //
    //  Calculate m_Dual 
    //
    //  a^{\alpha}, \alpha = 1, 2;
    //
    _aDual[0] 	= _metricTensorInverse(0,0) * _a[0] 
      + _metricTensorInverse(0,1) * _a[1];

    _aDual[1] 	= _metricTensorInverse(1,0) * _a[0] 
      + _metricTensorInverse(1,1) * _a[1];
  }

  void ShellGeometry::_updateDualPartials()
  {
    //
    //
    //  Calculate derivatives of dual w.r.t curvilinear coordinates
    //
    //  a(alpha, beta) = a^{alpha}_{beta}
    //////////////////////////////////////////////////////////////////////////////
    //
    // for alpha = 0, 1 and beta = 0, 1
    //
    // a(alpha, beta) = - sum { m_Dual(mu) * 
    //				( m_Dual(alpha) dot m_BasisDeriv(mu, beta)) }
    //
    // alpha = 0, beta = 0, mu = 0, 1
    //
    _aDualPartials[0] = - (  _aDual[0] * _aDual[0].dot(_aPartials[0])
			     + _aDual[1] * _aDual[0].dot( _aPartials[2] )
			     + _d * _aDual[0].dot( _dPartials[0]) );
    //
    // alpha = 0, beta = 1, mu = 0, 1
    _aDualPartials[2] = - (  _aDual[0] *  _aDual[0].dot( _aPartials[2])
			     + _aDual[1] *  _aDual[0].dot( _aPartials[1]) 
			     + _d        *  _aDual[0].dot( _dPartials[1]  ) );
    //
    // alpah = 1, beta = 0, mu = 0, 1
    _aDualPartials[3] = - (  _aDual[0] *  _aDual[1].dot(_aPartials[0])
			     + _aDual[1] *  _aDual[1].dot( _aPartials[2])  
			     + _d        *  _aDual[1].dot( _dPartials[0]  ) );
    //
    // alpha = 0, beta = 1, mu = 0, 1
    _aDualPartials[1] = - (  _aDual[0] *  _aDual[1].dot( _aPartials[2])
			     + _aDual[1] *  _aDual[1].dot( _aPartials[1]) 
			     + _d        *  _aDual[1].dot( _dPartials[1]  ) );

  }

  void ShellGeometry::_updateDirector()
  {
    //
    //  Calculate the covariant basis vector a_(2).
    //
    _d =  _a[0].cross( _a[1]);
    _d /= _metric;
  }

  void ShellGeometry::_updateDirectorPartials()
  {
    //
    //  Caclulate the derivatives of director w.r.t curvilinear coordinates
    //   dPartials[0] and  dPartials[1]  
    //
    _dPartials[0] = - (   _aDual[0] * _d.dot( _aPartials[0])
			  + _aDual[1] * _d.dot(  _aPartials[2]) );
    _dPartials[1] = - (   _aDual[0] * _d. dot( _aPartials[2])
			  + _aDual[1] *  _d.dot( _aPartials[1]) );
  }

  void ShellGeometry::_updateGamma()
  {
    vector<Matrix2d> g_abPartials(2,Matrix2d::Zero());

    g_abPartials[0] << 2*_a[0].dot(_aPartials[0]), 
      _a[0].dot(_aPartials[2]) + _a[1].dot(_aPartials[0]),
      _a[1].dot(_aPartials[0]) + _a[0].dot(_aPartials[2]),
      2*_a[1].dot(_aPartials[2]);  

    g_abPartials[1] << 2*_a[0].dot(_aPartials[2]),
      _a[0].dot(_aPartials[1]) + _a[1].dot(_aPartials[2]),
      _a[1].dot(_aPartials[2]) + _a[0].dot(_aPartials[1]),
      2*_a[1].dot(_aPartials[1]);  


    for (int m=0; m < 2; m++){
      for (int i=0; i < 2; i++){
	for (int j=0; j < 2; j++){
	  _GammaM_ij[m](i,j) = 0.5*_metricTensorInverse(m,0)*(g_abPartials[j](0,i) + g_abPartials[i](0,j) - g_abPartials[0](i,j) ) +
	    0.5*_metricTensorInverse(m,1)*(g_abPartials[j](1,i) + g_abPartials[i](1,j) - g_abPartials[1](i,j) );
	  
	}
      }
    }
  }
	
  ShellGeometry& ShellGeometry::operator = ( const ShellGeometry& g )
  {
    _a = g._a;
    _aPartials = g._aPartials;
    _d = g._d;
    _dPartials = g._dPartials;
    _aDual = g._aDual;
    _aDualPartials = g._aDualPartials;
    _metricTensor = g._metricTensor;
    _metricTensorInverse = g._metricTensorInverse;
    _metric = g._metric;
    _metricInverse = g._metricInverse;
    _GammaM_ij = g._GammaM_ij;
    return *this;
  } // end of overloading operator =
	
  const Vector2d ShellGeometry::gklGammaI_kl(){

    //cout << "GM_ij "<< GammaM_ij[1](0,1)<< " ... " << GammaM_ij[1](1,0)<<endl;
    Vector2d ret = Vector2d::Zero();
    ret(0) = ( _metricTensorInverse*_GammaM_ij[0] ).trace();
    ret(1) = ( _metricTensorInverse*_GammaM_ij[1] ).trace();
    return ret;
  }

	
} // namespace voom

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
  \file ShellGeometry.h

  \brief ShellGeometry is a structure which computes and contains all
  of the geometric quantities relevant to calculations on a curved surface.

*/

#if !defined(__ShellGeometry_h__)
#define __ShellGeometry_h__
 
#include<vector>
#include<iostream>
#include "voom.h"

namespace voom
{

  
  /*!  Structure which computes and contains all of the geometric
    quantities relevant to calculations on a curved surface.
  */
  class ShellGeometry 
  {
  private:
    vector<Vector3d>  _a;
    vector<Vector3d>  _aPartials;
    Vector3d   	      _d;
    vector<Vector3d>  _dPartials;
    vector<Vector3d>  _aDual;
    vector<Vector3d>  _aDualPartials;
    
    Matrix2d _metricTensor;
    Matrix2d _metricTensorInverse;

    double _metric;
    double _metricInverse;
		
  public:
    // Default Constructor
    ShellGeometry() {
      _a.resize(2, Vector3d::Zero()); // (a_1, a_2)
      _aPartials.resize(3, Vector3d::Zero()); // (a_11, a_22, a_12)
      _d = Vector3d::Zero(); // d
      _dPartials.resize(2, Vector3d::Zero()); // (d_1, d_2)
      _aDual.resize(2, Vector3d::Zero()); // (a^1, a^2)
      _aDualPartials.resize(3, Vector3d::Zero()); // (a^11, a^22, a^12)
      _metricTensor << 
	0.0, 0.0,
	0.0, 0.0;
      _metricTensorInverse <<
	0.0, 0.0,
	0.0, 0.0;
      _metric = 0.0;
      _metricInverse = 0.0;
    };
    
    //! Constructor
    ShellGeometry( const vector<Vector3d> &	a, 
		   const vector<Vector3d> & aPartials ) {
      update( a, aPartials );
    }
    
    void update( const  vector<Vector3d> & 	a,
		 const  vector<Vector3d> & aPartials ) {
      // recompute everything from a and aPartials 
    
      _a = a;
      _aPartials = aPartials;

      _updateMetrics();
      _updateDual();
      _updateDirector();
      _updateDirectorPartials();
      _updateDualPartials();			
    }

    const vector<Vector3d> &   a()	   const {return _a;}
    const vector<Vector3d> & aPartials() const {return _aPartials;}
    const Vector3d &		       d()	   const {return _d;}
    const vector<Vector3d> &   dPartials() const {return _dPartials;}
    const vector<Vector3d> &   aDual()	   const {return _aDual;}
    const vector<Vector3d> & aDualPartials() const {return _aDualPartials;}
    
    const Matrix2d & metricTensor()	  const {return _metricTensor;}
    const Matrix2d & metricTensorInverse() const {return _metricTensorInverse;}

    double metric()	const  {return _metric;}
    double metricInverse() const {return _metricInverse;}
		
    //! overload operator =
    ShellGeometry& operator = ( const ShellGeometry& g );
	  
  private:
    void _updateMetrics();
    void _updateDual();
    void _updateDirector();
    void _updateDualPartials();
    void _updateDirectorPartials();

  };


} // namespace voom
#endif // __ShellGeometry_h__

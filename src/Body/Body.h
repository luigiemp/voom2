//-*-C++-*-

#ifndef __Body_h__
#define __Body_h__

#include "voom.h"
#include "Mesh.h"
#include "Result.h"

// Include files for Writing Output:
#include <boost/lexical_cast.hpp>
#include <vtkVersion.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTensor.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkDoubleArray.h>


namespace voom {
  class Body {
  private :
    Body() {}; // Default constructor is private because it should never be used by a derived class.
  
  public:
    //! Constructor where different nodes have different DoF
    Body(Mesh* myMesh, State* myState)
      : _myMesh(myMesh), _myState(myState) {}

    //! Destructor
    virtual ~Body() {};

    //! Body compute function
    virtual void compute(Result* R) = 0;

    //! Check consistency - common to all bodies
    void checkConsistency(Result* R, Real perturbationFactor, int request = 6,
			  Real h = 1e-6, Real tol = 1e-6);

    //! Return number of unique material pointers
    // virtual int getNumMat() = 0;
    
    //! GetMesh
    Mesh* getMesh() {
      return _myMesh;
    }

    //! GetState
    State* getState() {
      return _myState;
    }

    //! Write body to paraview file
    virtual void writeOutputVTK(const string OutputFile, int step);

  protected:
    Mesh*        _myMesh;
    State*       _myState;
  
  }; // Body

} // namespace voom

#endif

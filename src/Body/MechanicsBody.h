//-*-C++-*-
#ifndef __MechanicsBody_h__
#define __MechanicsBody_h__

#include "Body.h"
#include "MechanicsMaterial.h"

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
#include <vtkDoubleArray.h>

namespace voom{

  // Body
  class MechanicsBody: public Body {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh and materials. */
    MechanicsBody(Mesh* aMesh, const int NodeDoF,
		  vector<MechanicsMaterial * > Materials,
		  Result* R);

    //! Destructor
    ~MechanicsBody() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (int i = 0; i < _materials.size(); i++)
	  UNIQUEmaterials.insert(_materials[i]);

      for (set<MechanicsMaterial *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++)
		delete (*it);
    };

    //! Initialize field
    // From constant value
    void initializeField(Result* R, Real fact = 1.0) {
      const int numNodes = _myMesh->getNumberOfNodes();
      const int dim = _myMesh->getDimension();
      
      for (int i = 0; i < numNodes; i++) { // i -> node number
	for (int j = 0; j < dim; j++) {
	  R->setField(i*dim+j, _myMesh->getX(i,j)*fact);
	}
      }
    };

    int getNumMat() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (int i = 0; i < _materials.size(); i++)
	UNIQUEmaterials.insert(_materials[i]);

      return UNIQUEmaterials.size();
    };

    int getTotNumMatProp() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (int i = 0; i < _materials.size(); i++)
	UNIQUEmaterials.insert(_materials[i]);
      // Assume all materials are of the same type

      return ( UNIQUEmaterials.size() * (_materials[0]->getMaterialParameters()).size() );
    }

    vector<MechanicsMaterial * > getMaterials() {
      return _materials;
    }

    //! Compute function
    void compute(Result* R);

    //! Check consistency of gradg and Hg
    void checkDmat(Result* R, Real perturbationFactor, Real h, Real tol);

    //! Write mechanics body to paraview file
    void writeOutputVTK(const string OutputFile, int step, Result* R);



  protected:
    //! Compute Deformation Gradient
    void computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl, Result* R);

    //! Compute Green Lagrangian Strain Tensor
    void computeGreenLagrangianStrainTensor(vector<Matrix3d> & Elist, GeomElement* geomEl, Result* R);

    //! List of Material data at each QP in the model
    vector<MechanicsMaterial * > _materials;

  };

} // namespace voom

#endif

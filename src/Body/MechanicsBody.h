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
// #include <vtkTensor.h>
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

    struct InvStruct
    {
      // Finite kinematics request and result type
      InvStruct() {};
    
      vector<vector<Real > > InvValue;
      vector<string > InvName;

      Real AverageValue(int ind) {
	Real AvgInv = 0.0;
	int NumQP = InvValue[ind].size();
	for (int i = 0; i < NumQP; i++) {
	  AvgInv += InvValue[ind][i];
	}
	AvgInv /= double(NumQP);
	return AvgInv;
      }

      vector<Real > InvQP(int q) {
	int NumInv = InvValue.size();
	vector<Real > QPInvariants(NumInv, 0.0);
	for (int i = 0; i < NumInv; i++) {
	  QPInvariants[i] = InvValue[i][q];
	}
	return  QPInvariants;
      }
    
    }; // struct InvStruct



    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh and materials. */
    MechanicsBody(Mesh* myMesh, State *myState,
		  vector<MechanicsMaterial * > Materials);

    //! Destructor
    ~MechanicsBody() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (int i = 0; i < _materials.size(); i++)
	  UNIQUEmaterials.insert(_materials[i]);

      for (set<MechanicsMaterial *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++)
		delete (*it);
    };

    //! Initialize field - Only of the nodes related to this body
    // From constant value
    void initializeField(Real fact = 1.0);

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
    void writeOutputVTK(const string OutputFile, int step);
    void writeQPdataVTK(const string OutputFile, int step);



  protected:
    //! Compute Deformation Gradient
    void computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl);

    //! Compute Green Lagrangian Strain Tensor
    void computeGreenLagrangianStrainTensor(vector<Matrix3d> & Elist, GeomElement* geomEl);

    //! Compute deformation invariants
    InvStruct computeInvariants(GeomElement* geomEl, int ElNum);

    //! List of Material data at each QP in the model
    vector<MechanicsMaterial * > _materials;
  };

} // namespace voom

#endif

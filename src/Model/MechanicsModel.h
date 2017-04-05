//-*-C++-*-
#ifndef __MechanicsModel_h__
#define __MechanicsModel_h__

#include "Model.h"
#include "MechanicsMaterial.h"
#include "EigenResult.h"

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

  // Model
  class MechanicsModel: public Model {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs.
     */
    MechanicsModel(Mesh* aMesh, const uint NodeDoF,  Result* myResult,
		   vector<MechanicsMaterial * > Materials,
		   int PressureFlag = 0, Mesh* SurfaceMesh = NULL,
		   int NodalForcesFlag = 0,
		   int _resetFlag = 1,
		   int _springBCflag = 0);

		   // const vector<string > & ElMatType,
		   // const map<string, MechanicsMaterial* > & ElMaterials);

    //! Input-file-based Constructor
    // MechanicsModel(Mesh* myMesh, const string inputFile, const uint NodeDoF);

    //! Destructor
    ~MechanicsModel() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++)
	  UNIQUEmaterials.insert(_materials[i]);

      for (set<MechanicsMaterial *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++)
		delete (*it);
    };



    //! Initialize field
    // From constant value
    void initializeField(const Real value = 1.0) {
      vector<Real> field;
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();
      field.resize(numNodes*_nodeDoF);
      
      for (uint i = 0; i < numNodes; i++) // i -> node number
	for (uint j = 0; j < dim; j++)
	  field[i*dim+j] = _myMesh->getX(i,j)*value; // value = isotropic expansion/shrinking
      _myResult->initializeField(field);
    };

    //! From array
    void initializeField(const Real* value) {
      vector<Real> field;
      field.assign(value, value+field.size());
      _myResult->initializeField(field);
    };

    //! Linearized update
    void linearizedUpdate(const Real* localValues, Real fact) {
      const int nLocalDof = (_myMesh->getNumberOfNodes())*_nodeDoF;
      for(uint i = 0; i < nLocalDof; i++)
	_myResult->_field[i] += fact*localValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      // const uint dim = _myMesh->getDimension();
      // assert( id < _field.size() && dof < dim );
      _myResult->_field[id*_nodeDoF + dof] += value;
    }

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int dof, const Real value) {
      _myResult->_field[dof] += value;
    }

    void setField(uint dof, Real value) {
      _myResult->_field[dof] = value;
    }
    void setField(const Real* value) {
      _myResult->_field.assign(value, value + _myResult->_field.size());
    };

    void getField(vector<Real > & x) {
      assert(x.size() == _myResult->_field.size());
      x = _myResult->_field;
    }

    void setPrevField(vector<Real> & prevField) {
      _prevField = prevField;
    }

    void setPrevField() {
      _prevField = _myResult->_field;
    };

    void printField() {
      int i = 0;
      while (i < _myResult->_field.size()) {
	for (uint j = 0; j < _nodeDoF; j++) {
	  cout << _myResult->_field[i] << " ";
	  i++;
	}
	cout << endl;
      }
    }

    void writeField(string OutputFile, int step) {
      // Create outputFile name
      stringstream FileNameStream;
      FileNameStream << OutputFile << step << ".dat";
      ofstream out;
      out.open( (FileNameStream.str()).c_str() );

      out << _myResult->_field.size() << endl;
      for (uint i = 0; i < _myResult->_field.size(); i++) {
	out << setprecision(15) << _myResult->_field[i] << endl;
      }
      out.close();
    }

    uint getNumMat() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++)
	UNIQUEmaterials.insert(_materials[i]);

      return UNIQUEmaterials.size();
    }

    uint getTotNumMatProp() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++)
	UNIQUEmaterials.insert(_materials[i]);
      // Assume all materials are of the same type

      return ( UNIQUEmaterials.size() * (_materials[0]->getMaterialParameters()).size() );
    }

    vector<MechanicsMaterial * > getMaterials() {
      return _materials;
    }

    void setResetFlag(int ResetFlag) {
      _resetFlag = ResetFlag;
    }

    void setPressureFlag(int PressureFlag) {
      _pressureFlag = PressureFlag;
    }

    void setNodalForcesFlag(int NodalForcesFlag) {
      _nodalForcesFlag = NodalForcesFlag;
    }

    //! Write output
    // void writeOutputVTK(const string OutputFile, int step);

    //! Write VTK output for normals of pressure
    // void writePressurePolyData(string OutputFile, int step);

    //! Write VTK output for linear spring
    // void writeLinearSpringPolyData(string OutputFile, int step);

    //! Write VTK output for torsional spring
    // void writeTorsionalSpringPolyData(string OutputFile, int step);

    //! Solve the system
    void compute(Result* R = this->_myResult);

    //! Finalize Compute (Optional method which includes computations done at the end of a solve step)
    // void finalizeCompute(); // Goes to spring body

    // Apply pressure
    // void applyPressure();

    // Update pressure
    // void updatePressure(Real Pressure) {
    //  _pressure = Pressure;
    // }

    // Update nodal forces // ????????
    void updateNodalForces(vector<int > * ForcesID, vector<Real > * Forces) {
      _forcesID = ForcesID; _forces = Forces;
    }

    // Check consistency of gradg and Hg
    void checkDmat(EigenResult * R, Real perturbationFactor, Real h, Real tol);

    // Functions for applying spring BC
    // Initialize _springNodes (nodes at which spring BC are applied) and _springElements (elements connected to spring nodes)
    // void initSpringBC(const string SpNodes, Mesh* SpMesh, Real SpringK);
    // void computeNormals();
    // vector<Triplet<Real > > applySpringBC();

    // Functions for applying Torsional Spring BC
    // void initTorsionalSpringBC(const string torsionalSpringNodes, Real torsionalSpringK);
    //! Computes centroid (x,y) of all the nodes but right now assumes long axis is z.
    // void computeCentroid();
    // void computeTangents();
    // vector<Triplet<Real> > applyTorsionalSpringBC();

    Real computeRefVolume();
    Real computeCurrentVolume();

  protected:
    //! Compute Deformation Gradient
    void computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl);

    //! Compute Green Lagrangian Strain Tensor
    void computeGreenLagrangianStrainTensor(vector<Matrix3d> & Elist, GeomElement* geomEl);

    //! List of Material data at each QP in the model
    vector<MechanicsMaterial * > _materials;

    // It should not be done here - maye we should have bodies and forms models from bodies
    // int _pressureFlag;
    // Real _pressure;
    // Mesh* _surfaceMesh;

    int _nodalForcesFlag;      // ?????
    vector<int > * _forcesID;  // ?????
    vector<Real > * _forces;   // ?????

    vector<Real > _prevField;

    // int _resetFlag; // ?????

    // Spring BC
    // int _springBCflag;
    // vector<int > _spNodes;
    // Mesh* _spMesh;
    // Real _springK;
    // vector<vector<int > > _spNodesToEle;
    // vector<Vector3d > _spNormals;

    // Torsional Spring BC
    // int _torsionalSpringBCflag;
    // vector<int> _torsionalSpringNodes;
    // Vector3d _centroidLocation;
    // Real _torsionalSpringK;
    // vector<Vector3d> _spTangents;
  };

} // namespace voom

#endif

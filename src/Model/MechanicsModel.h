//-*-C++-*-
#ifndef __MechanicsModel_h__
#define __MechanicsModel_h__

#include "Model.h"
#include "MechanicsMaterial.h"
#include "EigenResult.h"

namespace voom{

  // Model Results
  class MechanicsModel: public Model {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs. 
     */
    MechanicsModel(Mesh* aMesh, vector<MechanicsMaterial * > Materials, 
		   const uint NodeDoF,
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
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();

      for (uint i = 0; i < numNodes; i++) 
	for (uint j = 0; j < dim; j++)
	  _field[i*dim+j] = _myMesh->getX(i,j)*value; // value = isotropic expansion/shrinking
    };

    //! From array
    void initializeField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    //! Linearized update
    void linearizedUpdate(const Real* localValues, Real fact) {
      const int nLocalDof = (_myMesh->getNumberOfNodes())*_nodeDoF;
      for(uint i = 0; i < nLocalDof; i++)
	_field[i] += fact*localValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      // const uint dim = _myMesh->getDimension();
      // assert( id < _field.size() && dof < dim );
      _field[id*_nodeDoF + dof] += value;
    }

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int dof, const Real value) {
      _field[dof] += value;
    }
    
    void setField(uint dof, Real value) {
      _field[dof] = value;
    }
    void setField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    void getField(vector<Real > & x) {
      assert(x.size() == _field.size());
      x = _field;
    }

    void setPrevField() {
      _prevField = _field;
    };

    void printField() {
      int i = 0;
      while (i < _field.size()) {
	for (uint j = 0; j < _nodeDoF; j++) {
	  cout << _field[i] << " ";
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
  
      out << _field.size() << endl;
      for (uint i = 0; i < _field.size(); i++) {
	out << setprecision(15) << _field[i] << endl;
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
    void writeOutputVTK(const string OutputFile, int step); 

    //! Solve the system
    void compute(Result * R);

    // Apply pressure
    void applyPressure(Result * R);

    // Update pressure
    void updatePressure(Real Pressure) {
      _pressure = Pressure;
    }

    // Update nodal forces
    void updateNodalForces(vector<int > * ForcesID, vector<Real > * Forces) {
      _forcesID = ForcesID; _forces = Forces;
    }

    // Check consistency of gradg and Hg
    void checkDmat(EigenResult * R, Real perturbationFactor, Real h, Real tol);

    // Functions for applying spring BC
    // Initialize _springNodes (nodes at which spring BC are applied) and _springElements (elements connected to spring nodes)
    void initSpringBC(const string SpNodes, Mesh* SpMesh, Real SpringK);
    void computeNormals();
    vector<Triplet<Real > > applySpringBC(Result & R);
    
    Real computeRefVolume();
    Real computeCurrentVolume();

  protected:
    //! Compute Deformation Gradient
    void computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl);

    //! Compute Green Lagrangian Strain Tensor
    void computeGreenLagrangianStrainTensor(vector<Matrix3d> & Elist, GeomElement* geomEl);

    //! List of Material data at each QP in the model
    vector<MechanicsMaterial * > _materials;

    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]
    vector<Real > _field;

    // It should not be done here - maye we should have bodies and forms models from bodies
    int _pressureFlag;
    Real _pressure;
    Mesh* _surfaceMesh;

    int _nodalForcesFlag;
    vector<int > * _forcesID;
    vector<Real > * _forces;
    
    vector<Real > _prevField;

    int _resetFlag;

    // Spring BC
    int _springBCflag;
    vector<int > _spNodes;
    Mesh* _spMesh;
    Real _springK;
    vector<vector<int > > _spNodesToEle;
    vector<Vector3d > _spNormals;
  };

} // namespace voom

#endif

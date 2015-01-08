//-*-C++-*-
#ifndef __MechanicsModel_h__
#define __MechanicsModel_h__

#include "EllipticModel.h"
#include "MechanicsMaterial.h"

namespace voom{

  // Model Results
  class MechanicsModel: public EllipticModel {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs. 
     */
    MechanicsModel(Mesh* aMesh, 
		   const vector<string > & ElMatType, 
		   const map<string, MechanicsMaterial* > & ElMaterials);		   
    
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



    //! Linearized update (local and ghost solution)
    // From solution array (Used by Solver)
    void linearizedUpdate(const Real* localValues, 
			  const Real* ghostValues) {
      const int nLocalDof = (_myMesh->getLocalDoF()).size();
      for(uint i = 0; i < nLocalDof; i++)
	_field[i] += localValues[i];
      for(uint i = 0; i < (_myMesh->getGhostDoF()).size(); i++)
      _field[i + nLocalDof] += ghostValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      const uint dim = _myMesh->getDimension();
      assert( id < _field.size() && dof < dim );
      _field[id*dim + dof] += value;
    }
    
    void setField(uint dof, Real value) {
      _field[dof] = value;
    }

    void PrintField() {
      for (uint i = 0; i < _field.size(); i++)
	cout << _field[i] << endl;
    }
    
    //! Write output
    void writeOutput(const string OutputFile, const string format = "BINARY"); 

    //! Solve the system
    void compute(EllipticResult & R);

    //! Return component of residual vector given nodal id and DoF
    /*Real getResidualComponent(EllipticResult & R, int ID, int dim) {
      assert(dim == 0);
      return 0.0; //R.getResidual(ID);
    };
    */



  protected:
    //! Compute Deformation Gradient
    void computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl);

    //! List of Material data at each element in the model
    // (need to be modified for history dependent materials, e.g. plasticity)
    vector<MechanicsMaterial * > _materials;

    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]
    vector<Real > _field;
  };

} // namespace voom

#endif

//-*-C++-*-
#ifndef __PoissonModel_h__
#define __PoissonModel_h__

#include "EllipticModel.h"
#include "DiffusionMaterial.h"

namespace voom{

  class PoissonModel: public EllipticModel {
  public:

    //! Constructor
    PoissonModel(Mesh* myMesh, 
		 const vector<string > & ElMatType, 
		 const map<string, DiffusionMaterial* > & ElMaterials);

    //! Destructor
    ~PoissonModel() {
      set<DiffusionMaterial *> UNIQUEmaterials;
      for (uint i =0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      for (set<DiffusionMaterial *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++) 
		delete (*it);
    };

    //! Initialize field
    // From constant value
    void initializeField(const Real value = 0.0) {
      _field.assign(_field.size(), value);
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
      assert( id < _field.size() && dof == 0 );
      _field[id] += value;
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
    void compute(EllipticResult& result);

    //! Return component of residual vector given nodal id and DoF
    /*Real getResidualComponent(EllipticResult & R, int ID, int dim) {
      assert(dim == 0);
      return 0.0; //R.getResidual(ID);
    };
    */

   

protected:
    // Default constructor cannot be used, Material need to be provided
    PoissonModel();
    //! List of material data at each quad point in the model
    vector<DiffusionMaterial *> _materials;

    //! Solution value at all nodes, local and ghost
    vector<Real > _field;

  };

} // namespace voom

#endif

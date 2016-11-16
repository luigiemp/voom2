//-*-C++-*-
#ifndef __FoundationModel_h__
#define __FoundationModel_h__

#include "Model.h"
#include "EigenResult.h"

namespace voom{

  // Model Results
  class FoundationModel: public Model {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs.
    */
    FoundationModel(Mesh* aMesh, vector<MechanicsMaterial * > Materials, const uint NodeDoF, Real SpringK);

    //! Destructor
    ~FoundationModel() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++)
      UNIQUEmaterials.insert(_materials[i]);

      for (set<MechanicsMaterial *>::iterator it = UNIQUEmaterials.begin();
      it != UNIQUEmaterials.end(); it++)
      delete (*it);
    };

    // **** Begin derived function implementation **** //
    //! Solve the system
    virtual void compute(Result * R);

    virtual uint getNumMat() {
      set<MechanicsMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++)
      UNIQUEmaterials.insert(_materials[i]);
      return UNIQUEmaterials.size();
    }

    //! Initialize field
    // From constant value
    virtual void initializeField(const Real value = 1.0) {
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();

      for (uint i = 0; i < numNodes; i++)
        for (uint j = 0; j < dim; j++)
          _field[i*dim+j] = _myMesh->getX(i,j)*value; // value = isotropic expansion/shrinking
    };

    //! From array
    virtual void initializeField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    //! Linearized update
    virtual void linearizedUpdate(const Real* localValues, Real fact) {
      const int nLocalDof = (_myMesh->getNumberOfNodes())*_nodeDoF;
      for(uint i = 0; i < nLocalDof; i++)
      _field[i] += fact*localValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    virtual void linearizedUpdate(const int id, const int dof, const Real value) {
      _field[id*_nodeDoF + dof] += value;
    }

    // One value at the time (Node ID, dof index, value)
    virtual void linearizedUpdate(const int dof, const Real value) {
      _field[dof] += value;
    }

    virtual void setField(uint dof, Real value) {
      _field[dof] = value;
    }
    virtual void setField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    virtual void getField(vector<Real > & x) {
      assert(x.size() == _field.size());
      x = _field;
    }

    virtual void setPrevField() {
      _prevField = _field;
    };

    virtual void printField() {
      int i = 0;
      while (i < _field.size()) {
        for (uint j = 0; j < _nodeDoF; j++) {
          cout << _field[i] << " ";
          i++;
        }
        cout << endl;
      }
    }

    virtual void writeField(string OutputFile, int step) {
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

    //! Write output
    virtual void writeOutputVTK(const string OutputFile, int step);
    // **** End derived function implementation **** //



    // Functions for applying spring BC
    // Initialize _springNodes (nodes at which spring BC are applied) and _springElements (elements connected to spring nodes)
    void initSpringBC(const string SpNodes, Mesh* SpMesh, Real SpringK);
    void computeNormals();
    vector<Triplet<Real > > applySpringBC(Result & R);

  protected:
    //! List of Material data at each QP in the model
    // vector<MechanicsMaterial * > _materials;

    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]
    vector<Real > _field;

    //! Previous field - required here for computing spring force
    //! Mutated using setPrevField()
    vector<Real > _prevField;

    int _resetFlag;

    // Spring BC
    vector<int > _spNodes;
    Real _springK;
    vector<vector<int > > _spNodesToEle;
    vector<Vector3d > _spNormals;
  };

} // namespace voom

#endif

//-*-C++-*-
#ifndef __State_h__
#define __State_h__

#include "voom.h"

namespace voom {
  //! State is a class which contains nodal positions and DoF

  class State
  {
  public:

    // Constructor --- Notice: it does not initialize _X or _Gdof vectors
    State(): _DOFcount(0) {};
    State(int NumNodes): _DOFcount(0) { 
      _X.reserve(NumNodes);
      _Gdof.reserve(NumNodes);
      _NodeDof.reserve(NumNodes);
      _phi.reserve(NumNodes);
    }; 

    // Accessor and mutators
    void insertNode(Vector3d Xnew, int dofPerNode) {
      _X.push_back(Xnew);
      for (int i = 0; i < dofPerNode; i++) {
	_phi.push_back(0.0);
      }
      _Gdof.push_back(_DOFcount);
      _NodeDof.push_back(dofPerNode);
      _DOFcount += dofPerNode;
    }

    int getDOFcount() {return _DOFcount; }
    int getGdof(int node) {return _Gdof[node]; }
    int getNodeDof(int node) {return _NodeDof[node]; }
    Vector3d getX(int node) {return _X[node]; }
    Real getX(int node, int ind) {return _X[node](ind); }
    int getXsize() {return _X.size(); }

    Real getPhi(int Gdof) {return _phi[Gdof]; };
    Real getPhi(int node, int Ldof) {
      assert(Ldof < _NodeDof[node]);
      return _phi[_Gdof[node] + Ldof]; 
    };
    void getPhi(vector<Real > & Phi) {Phi = _phi; };

    void setPhi(int Gdof, Real value) {_phi[Gdof] = value; };
    void setPhi(int node, int Ldof, Real value) {
      assert(Ldof < _NodeDof[node]);
      _phi[_Gdof[node] + Ldof] = value; 
    };
    void setPhi(const vector<Real > & Phi, Real fact = 1.0) {
      assert(_phi.size() == Phi.size());
      for (int i=0; i++; i<Phi.size()) {
	_phi[i] = Phi[i]*fact;
      }; 
    }
     


    void linearizedUpdate(int Gdof, Real value) {_phi[Gdof] += value; }; // If different linearizeupdate is needed for different DoF type, this function can be moved to Bodies
    void linearizedUpdate(int node, int Ldof, Real value) {
      assert(Ldof < _NodeDof[node]);
      _phi[_Gdof[node] + Ldof] += value; };
    void linearizedUpdate(const vector<Real > & DeltaPhi, Real fact = 1.0) {
      assert(_phi.size() == DeltaPhi.size());
      for (int i=0; i<DeltaPhi.size(); i++) {
	_phi[i] += DeltaPhi[i]*fact;
      }; 
    }
    void linearizedUpdate(const Real * values, Real fact = 1.0) {
      for (int i=0; i<_DOFcount; i++) {
	_phi[i] += values[i]*fact;
      }; 
    }
    
    void printPhi() {
      int i = 0;
      while (i < _phi.size()) {
	for (int j = 0; j < _NodeDof[j]; j++) {
	  cout << _phi[i] << " ";
	  i++;
	}
	cout << endl;
      }
    }

    virtual void writePhi(string OutputFile, int step) {
      stringstream FileNameStream;
      FileNameStream << OutputFile << step << ".dat";
      ofstream out;
      out.open( (FileNameStream.str()).c_str() );

      out << _phi.size() << endl;
      for (int i = 0; i < _phi.size(); i++) {
	out << setprecision(15) << _phi[i] << endl;
      }
      out.close();
    };

    // Protected
    vector<Vector3d > _X;
    vector<int >      _Gdof;
    vector<int >      _NodeDof;
    vector<Real >     _phi;
    int               _DOFcount;  // Keeps track of the total number of DOF in the problem
 
  }; // End of declaration of State base class

} // namespace voom

#endif




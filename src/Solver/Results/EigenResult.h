//-*-C++-*-
#ifndef __EigenResult_h__
#define __EigenResult_h__
#include "Result.h"

namespace voom {
  //! Eigen Elliptic Result is derived from the general EllipticResult class and contains Eigen storage objects. Will be used by Eigen solver class

  class EigenResult: public Result
  {

  public:

    //! Constructor and destructor
    EigenResult(int PbDoF, int NumMatProp): _numMatProp(NumMatProp)
    {
      _field = VectorXd(_pbDoF);
      // Create results structures
      _stiffness = new SparseMatrix<Real >(_pbDoF, _pbDoF);
      _residual = new VectorXd(_pbDoF);
      
      _Hg = new SparseMatrix<Real >(NumMatProp, NumMatProp);
      _Gradg = new VectorXd(NumMatProp);
    };

    ~EigenResult() {
      delete _Hg;
      delete _Gradg;
      delete _stiffness;
      delete _residual;
    };

    // Functions relating to field
    void fieldResize(int size) {_field.resize(size);};
    void initializeField(const vector<Real> field){
	assert(field.size() == _field.size());
	_field = field;
    };
    void initializeField(int ind, Real value)
    int  getFieldSize() {return _field.size();};
    // All values at once
    void linearizedUpdate(const Real* locaValues, Real fact) {
      // No checking on size of localValues, prone to seg fault! 
      for(uint i = 0; i < _field.size(); i++)
	_field[i] += fact*localValues[i]; };
    // One value at the time (Node ID, local DoF index, value)
    void linearizedUpdate(const int id, const int LocalDoF, const Real value, const int nodeDoF = 3) {
      _field[id*nodeDoF + LocalDoF] += value; };
    // One value at the time (Global DoF index, value)
    void linearizedUpdate(const int GlobalDoF, const Real value) {
      _field[GlobalDoF] += value;
    }

    void setField(uint GlobalDoF, Real value) {
      _field[GlobalDoF] = value; };
    void getField(vector<Real > & x) {
      assert(x.size() == _field.size());
      x = _field; };

    void printField(const int nodeDoF = 3) {
      int i = 0;
      while (i < _field.size()) {
	for (uint j = 0; j < nodeDoF; j++) {
	  cout << _field[i] << " ";
	  i++;
	}
	cout << endl;
      }
    };

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
    };



    // Return NumMatProp and PbDoF
    int getNumMatProp() {
      return _numMatProp;
    };
    
    //! Interface (Mutators)
    // Set functions
    void setResidual(int ind, Real value) {
      (*_residual)(ind) = value;
    };
    
    // Reset function
    void resetResidualToZero() {
      *_residual = VectorXd::Zero(_pbDoF);
    }

    void resetStiffnessToZero() {
      _stiffness->setZero();
    }

    void resetGradgToZero() {
      *_Gradg = VectorXd::Zero(_numMatProp);
    }

    void resetHgToZero() {
      _Hg->setZero();
    }



    // Add functions
    void addResidual(int ind, Real value) {
      (*_residual)(ind) += value;
    }
    
    void addStiffness(int indRow, int indCol, Real value) {
      _stiffness->coeffRef(indRow, indCol) += value;
    };

    void FinalizeGlobalStiffnessAssembly() {
      _stiffness->setFromTriplets(KtripletList.begin(), KtripletList.end());
      _stiffness->makeCompressed();
    };

    void addGradg(int ind, Real value) {
      (*_Gradg)(ind) += value;
    }

    void addHg(int indRow, int indCol, Real value) {
      _Hg->coeffRef(indRow, indCol) += value;
    }

    void FinalizeHgAssembly() {
      _Hg->setFromTriplets(HgtripletList.begin(), HgtripletList.end());
      _Hg->makeCompressed();
    };
  


    //! Interface (Accessors)
    Real getResidual(int ind) {
      return (*_residual)(ind);
    }

    Real getStiffness(int indRow, int indCol) {
      return _stiffness->coeff(indRow, indCol);
    };

    Real getGradg(int ind) {
      return (*_Gradg)(ind);
    }

    Real getHg(int indRow, int indCol) {
      return _Hg->coeff(indRow, indCol);
    };



  protected:
    int _pbDoF;
    int _numMatProp;

    VectorXd *_field;
    VectorXd *_residual;
    vector<Triplet<Real > > KtripletList, HgtripletList;
    SparseMatrix<Real > *_stiffness;
    
    SparseMatrix<Real > *_Hg;
    VectorXd *_Gradg;

  }; // EigenResult 
  
}; // namespace voom

#endif

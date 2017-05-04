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
    EigenResult(int PbDoF, int NumMatProp): Result(PbDoF), _numMatProp(NumMatProp)
    {
      _field    = VectorXd::Zero(_pbDoF);
      _residual = VectorXd::Zero(_pbDoF);
      _Gradg    = VectorXd::Zero(NumMatProp);
      
      _Hg = new SparseMatrix<Real >(NumMatProp, NumMatProp);
      _stiffness = new SparseMatrix<Real >(_pbDoF, _pbDoF);
    };


    ~EigenResult() {
      delete _Hg;
      delete _stiffness;
    };


    // Initialization function
    void initializeResults(int NumEntries) {
      KtripletList.reserve(NumEntries);
      HgtripletList.reserve(square(_numMatProp));
    };


    // NumMat
    int getNumMatProp() {
      return _numMatProp;
    };
    

    // Field
    void setField(int ind, Real value) {
      _field(ind) = value;
    };
    void linearizedUpdate(const vector<Real > & Field, Real fact) {
      assert(_field.size() == Field.size());
      _field += Field*fact; 
    };
    void linearizedUpdate(int ind, Real value) {
      _field(ind) += value;
    };
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
    Real getField(int ind) {
      return _field(ind); 
    };



    // Residual
    void resetResidualToZero() {
      _residual = VectorXd::Zero(_pbDoF);
    };
    void setResidual(int ind, Real value) {
      _residual(ind) = value;
    };
    void addResidual(int ind, Real value) {
      _residual(ind) += value;
    };
    Real getResidual(int ind) {
      return _residual(ind);
    };
    


    // Stiffness
    void resetStiffnessToZero() {
      _stiffness->setZero();
    };
    void addStiffness(int indRow, int indCol, Real value) {
      _stiffness->coeffRef(indRow, indCol) += value;
    };
    void FinalizeGlobalStiffnessAssembly() {
      _stiffness->setFromTriplets(KtripletList.begin(), KtripletList.end());
      _stiffness->makeCompressed();
      KtripletList.clear();
    };
    Real getStiffness(int indRow, int indCol) {
      return _stiffness->coeff(indRow, indCol);
    };



    // Gradg
    void resetGradgToZero() {
      _Gradg = VectorXd::Zero(_numMatProp);
    };
    void addGradg(int ind, Real value) {
      _Gradg(ind) += value;
    };
    Real getGradg(int ind) {
      return _Gradg(ind);
    };



    // Hg
    void resetHgToZero() {
      _Hg->setZero();
    };
    void addHg(int indRow, int indCol, Real value) {
      _Hg->coeffRef(indRow, indCol) += value;
    };
    void FinalizeHgAssembly() {
      _Hg->setFromTriplets(HgtripletList.begin(), HgtripletList.end());
      _Hg->makeCompressed();
      HgtripletList.clear();
    };
    Real getHg(int indRow, int indCol) {
      return _Hg->coeff(indRow, indCol);
    };



  protected:
    int _numMatProp;

    VectorXd _field;
    VectorXd _residual;
    VectorXd _Gradg;
    vector<Triplet<Real > > KtripletList, HgtripletList;

    SparseMatrix<Real > *_stiffness;    
    SparseMatrix<Real > *_Hg;

  }; // EigenResult 
  
}; // namespace voom

#endif

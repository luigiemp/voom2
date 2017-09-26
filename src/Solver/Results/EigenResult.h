//-*-C++-*-
#ifndef __EigenResult_h__
#define __EigenResult_h__
#include "Result.h"
#include "VoomMath.h"

namespace voom {
  //! Eigen Elliptic Result is derived from the general EllipticResult class and contains Eigen storage objects. Will be used by Eigen solver class

  class EigenResult: public Result
  {

  public:

    //! Constructor and destructor
    EigenResult(int PbDoF, int NumMatProp): Result(PbDoF), _numMatProp(NumMatProp)
    {
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
      _KtripletList.reserve(NumEntries);
      _HgtripletList.reserve(square(_numMatProp));
    };



    // NumMat
    int getNumMatProp() {
      return _numMatProp;
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
      _KtripletList.clear();
    };
    void addStiffness(int indRow, int indCol, Real value) {
      _KtripletList.push_back( Triplet<Real >( indRow, indCol, value ) );
    };
    void FinalizeGlobalStiffnessAssembly() {
      _stiffness->setFromTriplets(_KtripletList.begin(), _KtripletList.end());
      _stiffness->makeCompressed();
      _KtripletList.clear();
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
      _HgtripletList.clear();
    };
    void addHg(int indRow, int indCol, Real value) {
      _HgtripletList.push_back( Triplet<Real >( indRow, indCol, value ) );
      // _Hg->coeffRef(indRow, indCol) += value;
    };
    void FinalizeHgAssembly() {
      _Hg->setFromTriplets(_HgtripletList.begin(), _HgtripletList.end());
      _Hg->makeCompressed();
      _HgtripletList.clear();
    };
    Real getHg(int indRow, int indCol) {
      return _Hg->coeff(indRow, indCol);
    };



  public:
    int _numMatProp;

    VectorXd _residual;
    VectorXd _Gradg;
    vector<Triplet<Real > > _KtripletList, _HgtripletList;

    SparseMatrix<Real > *_stiffness;    
    SparseMatrix<Real > *_Hg;

  }; // EigenResult 
  
}; // namespace voom

#endif

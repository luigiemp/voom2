//-*-C++-*-
#ifndef __EPResult_h__
#define __EPResult_h__
#include "Model.h"

namespace voom{
  /*!
    Results for Reaction diffusion equation relating to cardiac
    electrophysiology. Defines a basic interface to access/update the 
    result values
   */

  struct EPResult {
    // Default constructor
    EPResult(): _request(0) {};

    // Virtual function to set epresults members
    virtual void setRequest(int request) { _request = request; }

    // Pure Virtual members
    // Manipulate Diffusion solve terms
    virtual void addDiffusionLHSMatrix(const int localRow, const int localCol,
				       const Real value) = 0;
    virtual void addDiffusionRHSMatrix(const int localRow, const int localCol,
				       const Real value) = 0;    

    // Manipulate Ionic solve terms
    virtual void addIonicLHSMatrix(const int localRow, const int localCol,
				   const Real value) = 0;
    virtual void addIonicRHSMatrix(const int localRow, const int localCol,
				   const Real value) = 0;    

    // Manipulate Ionic Current
    virtual void addIonicCurrent(const int localIndex, const Real value) = 0;

    // Finalize assembly of matrices
    virtual void finalizeMatrixAssembly() = 0;

  private:
    int _request;
  };
}
#endif

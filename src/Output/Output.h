//-*-C++-*-
/*!
  \file Output.h
  \brief A base output class which creates output files
 */
#ifndef __Output_h__
#define __Output_h__

#include "voom.h"

namespace voom {

  class Output {
  public:
    Output(string outputFileName)
      : _outputFileName(outputFileName)
    {
      // AVSP: Todo: Fill in Constructor
    }

    //! Destructor
    virtual ~Output() {};

    enum DataType {POINT_DATA, CELL_DATA};

    //! Function to insert a scalar
    virtual void insertScalar(DataType type, string scalarName, vector<double>& scalarArray) = 0;

    //! Function to insert a vector
    virtual void insertVector(DataType type, string vectorName, vector<vector <double> >& vectorArray) = 0;

    //! Function to insert tensor
    virtual void insertTensor(DataType type, string tensorName, vector<Matrix3d>& tensorArray) = 0;

  private :
    Output() {}; // Default constructor is private because it should never be used by a derived class.

  protected:
    string _outputFileName;  
  }; // Output

} // namespace voom

#endif

//-*-C++-*-
/*!
  \file Output.h
  \brief A output class which creates VTK type output
 */
#ifndef __VTKOutput_h__
#define __VTKOutput_h__

#include "output.h"
#include "GeomElement.h"

// Todo: 07.27.2016 - Refactor code to use const refs everywhere.

namespace voom {

  class VTKOutput : public Output {
  public:

    enum vtkDataTypes{STRUCTUREDGRID,
		      UNSTRUCTUREDGRID,
                      POLYDATA};

    //! Constructor
    VTKOutput(string outputFileName, vtkDataTypes type);

    //! Destructor
    ~VTKOutput() {};

    //! Function to insert a scalar
    void insertScalar(DataType type, string scalarName, vector<double>& scalarArray);

    //! Function to insert a vector
    void insertVector(DataType type, string vectorName, vector<vector <double> >& vectorArray);

    //! Function to insert tensor
    void insertTensor(DataType type, string tensorName, vector<Matrix3d>& tensorArray);

    //! Function to insert nodal positions
    void insertNodalPositions(const vector<Vector3d>& nodalPositions);

    //! Function to insert connectivity
    void insertConnectivity(const vector<GeomElement *>& elements);

    //! Write output
    void writeOutput()
    {
      _writer.SetFileName(_outputFileName.c_str());
      _writer.SetInput(_dataSet);
      _writer.Write();
    }

  private :
    VTKOutput() {}; // Default constructor is private because it should never be used by a derived class.

  protected:
   vtkSmartPointer<vtkDataWriter> _writer;
   vtkDataSet*                    _dataSet;

  }; // VTKOutput

} // namespace voom

#endif

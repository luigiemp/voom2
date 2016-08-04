//-*-C++-*-
/*!
  \file VTKOutput.h
  \brief A output class which creates VTK type output
 */
#ifndef __VTKOutput_h__
#define __VTKOutput_h__

#include "Output.h"
#include "GeomElement.h"

#include <vtkSmartPointer.h>
#include <vtkStructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkDataWriter.h>
#include <vtkDataSet.h>

// Todo: 07.27.2016 - Refactor code to use const refs everywhere.

namespace voom {

  class VTKOutput : public Output {
  public:

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
      _writer->SetFileName(_outputFileName.c_str());
      _writer->SetInput(_dataSet);
      _writer->Write();
    }

  protected:
   vtkSmartPointer<vtkDataWriter> _writer;
   vtkSmartPointer<vtkDataSet>    _dataSet;

  }; // VTKOutput

} // namespace voom

#endif

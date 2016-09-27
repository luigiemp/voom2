#include "VTKOutput.h"

namespace voom {

  VTKOutput::VTKOutput(string outputFileName, vtkDataTypes type):
    Output(outputFileName) {
  
    switch (type) {
      case STRUCTUREDGRID:
        _writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
        _dataSet = vtkSmartPointer<vtkStructuredGrid>::New();
        break;
      case UNSTRUCTUREDGRID:
        _writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        _dataSet = vtkSmartPointer<vtkUnstructuredGrid>::New();
        break;
      case POLYDATA:
        _writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        _dataSet = vtkSmartPointer<vtkPolyData>::New();
        break;
      default:
        cout << "** ERROR: Unknown VTK Dataset type!" << endl;
        exit(EXIT_FAILURE);
    }
  }

  void VTKOutput::insertScalar(DataType type, string scalarName, vector<double>& scalarArray) {
    vtkSmartPointer<vtkDoubleArray> vtkScalarArray = vtkSmartPointer<vtkDoubleArray>::New();
    vtkScalarArray->SetName(scalarName.c_str());
    vtkScalarArray->SetNumberOfComponents(1);
    vtkScalarArray->SetNumberOfTuples(scalarArray.size());
    
    for (int array_iter = 0; array_iter < scalarArray.size(); array_iter++)
      vtkScalarArray->SetTuple1(array_iter, scalarArray[array_iter]);
  }

  void VTKOutput::insertVector(DataType type, string vectorName, vector<vector <double> >& vectorArray) {
    vtkSmartPointer<vtkDoubleArray> vtkVectorArray = vtkSmartPointer<vtkDoubleArray>::New();
    vtkVectorArray->SetName(vectorName.c_str());
    vtkVectorArray->SetNumberOfComponents(vectorArray[0].size());
    vtkVectorArray->SetNumberOfTuples(vectorArray.size());

    for (int array_iter = 0; array_iter < vectorArray.size(); array_iter++) {
      for (int comp_iter = 0; comp_iter < vectorArray[0].size(); comp_iter++)
        vtkVectorArray->SetComponent(array_iter, comp_iter, vectorArray[array_iter][comp_iter]);
    }
  }

  void VTKOutput::insertTensor(DataType type, string tensorName, vector<Matrix3d>& tensorArray) {
    vtkSmartPointer<vtkDoubleArray> vtkVectorArray = vtkSmartPointer<vtkDoubleArray>::New();
    vtkVectorArray->SetName(tensorName.c_str());
    vtkVectorArray->SetNumberOfComponents(9);
    vtkVectorArray->SetNumberOfTuples(tensorArray.size());

    for (int array_iter = 0; array_iter < tensorArray.size(); array_iter++) {
      for (int row_iter = 0; row_iter < 3; row_iter++) {
        for (int col_iter = 0; col_iter < 3; col_iter++) {
          vtkVectorArray->SetComponent(array_iter, row_iter*3 + col_iter, (tensorArray[array_iter])(row_iter, col_iter));
        }
      }
    }
  }
}

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
}

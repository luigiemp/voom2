#include "voom.h"
#include "Output.h"
#include "VTKOutput.h"

using namespace std;
using namespace voom;

int main(int argc, char** argv) {

  Output* output = new VTKOutput("LuigiFileName.vtk", UNSTRUCTUREDGRID);

}

//-*-C++-*-
/*!
  \file Model.h
  \brief A base model class creates basic entities needs for solving in parallel
  and stores generic boundary conditions. 
 */
#ifndef __Model_h__
#define __Model_h__

#include "voom.h"
#include "Mesh.h"

namespace voom {
  class Model {
  private :
    Model() {}; // Default constructor is private because it should never be used by a derived class.
  
  public:
    //! Constructor using an input file ... to be changed after implementation of input class - no longer current
    // Model(Mesh* myMesh, const string inputFile, const uint NodeDoF);

    //! Constructor to be used
    Model(Mesh* aMesh, const uint NodeDoF)
      : _myMesh(aMesh), _nodeDoF(NodeDoF) {
      // WSK: evetually we want to copy the mesh object that is
      // pointed to.  To do this we need to implement a virtual copy
      // method to be used like this:
      //
    }

    //! Destructor
    virtual ~Model() {};

    //! Get DoF per Node
    uint getDoFperNode() {
      return _nodeDoF;
    }
    
    //! GetMesh
    Mesh * getMesh() {
      return _myMesh;
    }

    //! Initialize field
    // From constant value
    virtual void initializeField(const Real value = 0.0) = 0;
 
    //! From array
    virtual void initializeField(const Real* value) = 0;

    //! Linearized update (local and ghost solution)
    // From solution array (Used by Solver)
    // virtual void linearizedUpdate(const Real* locaValues, 
    // 				  const Real* ghostValues) = 0;
    virtual void linearizedUpdate(const Real* locaValues) = 0;

    // One value at the time (Node ID, DoF index, value)
    virtual void linearizedUpdate(const int id, const int dof, const Real value) = 0;

    // Set solution field
    virtual void setField(uint dof, Real value) = 0;

    // Print field
    virtual void PrintField() = 0;

    //! Output generation
    /*!
      Output format is
      nLocalNodes   Dof_per_Node   Name_of_Dof1 Name_of_Dof2 ....
      nElements     Output_per_Element  Name_output1 Name_output2 ...
      nSteps

      First all nodal outputs are written
      Dof_Value forall nodes
      Output_Value forall elements
      
      This is repeated nSteps times
    */
    virtual void writeOutput(const string OutputFile, const string format = "BINARY") = 0;

    // Initialize output file
    void initializeOutput(const int DOF_PER_NODE   , const char* DOF_NAMES,
			  const int OUTPUT_PER_ELEM, const char* OUTPUT_NAMES,
			  const int NSTEPSconst,
			  const string OutputFile, const string format = "BINARY");

    /*! 
      Gateway to Model class. 
      \param mesh Pointer to parent Mesh class
      \param inputFile Input v2i file from which loading information is 
      parsed to determine the type of analysis to be performed.
    */
    // static Model* New(Mesh* mesh, const string inputFile, const uint NodeDoF);

    


  protected:
    Mesh*        _myMesh;

    //! DoF per node
    uint         _nodeDoF;

    // //! Parse Input deck to get list of Material Names
    // void _getMaterialNamesList(vector<string > & names);
  
  }; // Model

} // namespace voom

#endif

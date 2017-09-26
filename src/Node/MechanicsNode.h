//-*-C++-*-

#ifndef __MechanicsNode_h__
#define __MechanicsNode_h__
#include "BaseNode.h"

namespace voom {
  class MechanicsNode : public BaseNode {
  public:
    //! Constructor
    MechanicsNode(int ID, vector<int > DOF, Vector3d X): 
      BaseNode(ID, DOF, X), _field(X) {}
    MechanicsNode(int ID, vector<int > DOF, Vector3d X, Vector3d Field): 
      BaseNode(ID, DOF, X), _field(Field) {}

    //! Accessors and Mutators
    Vector3d getField() { return _field; }
    void setField(Vector3d Field) { _field = Field; }
    void linearizeUpdate(Vector3d DeltaField) { _field += DeltaField; }
    
  protected:
    //! Reference position
    Vector3d _field;
  };
}
#endif

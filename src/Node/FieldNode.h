//-*-C++-*-
/*!
  \file FieldNode.h
  \brief Field Node derived from base node. Field node needs position to 
  define it. Hence the reference position of the node is known.
*/

#ifndef __FieldNode_h__
#define __FieldNode_h__
#include "BaseNode.h"

namespace voom {
  class FieldNode : public BaseNode {
  public:
    //! Constructor
    FieldNode(unsigned int id, unsigned int dof, VectorXd& X): 
      BaseNode(id, dof), _X(X) {};

    //! Default Constructor
    FieldNode(unsigned int id, unsigned int dof):BaseNode(id, dof), 
						 _X(VectorXd::Zero(dof)) {};

    //! Get reference position
    const VectorXd& getX() { return _X; }

    //! Set reference position
    void setX(VectorXd& X) { 
      assert(_X.size() == X.size() );
      _X = X;
    }
    
  protected:
    //! Reference position
    VectorXd     _X;
  private:
    //! Contructor with no arguments cannot be used by users.
    FieldNode() {};
  };
}
#endif

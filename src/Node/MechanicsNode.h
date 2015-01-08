//-*-C++-*-
/*!
  \file MechanicsNode.h
  \brief Mechanics node derived from Field Node. In addition to reference
  position we also need displacement at the node. In case of small displacement 
  problems it is better to store displacement instead of current position. 
  This is purely from numerical standpoint.
*/

#ifndef __MechanicsNode_h__
#define __MechanicsNode_h__

#include "FieldNode.h"

namespace voom{
  class MechanicsNode : public FieldNode {
  public:
    //! Constructor
    MechanicsNode(unsigned int id, unsigned int dof, VectorXd& X):
      FieldNode(id, dof, X), _u(VectorXd::Zero(X.size()) ) {};
    
    //! Get current position
    const VectorXd& getu() { return _u; }

    //! Set displacement
    void setu(const VectorXd& u);

    //! Set particular displacement component
    void setu(unsigned int id, Real value);

    //! Get property associated with a node
    vector<string> getPropertyList();

    /*!
      Given a propertyName as input the Real data is filled with values
      and returned to the user. The integer length denotes the number of
      values stored in the Real array.
      \param propertyName String name which defines which property needs to 
      be returned
      \param data Real array where results are stored. Memory is 
      allocated in this routine. Whatever is stores in data is cleared 
      initially and values are stored in it.
      \param length Number of entires stored in the array data.
    */
    virtual void getProperty(string propertyName, Real** data, uint& length);

    //! Linear Update
    virtual void linearizedUpdate(const VectorXd& data);

    //! Get particular value in the property list.
    virtual void linearizedUpdate(unsigned int dof, Real data);
    
  protected:
    //! Displacement 
    VectorXd     _u;
  };
}
#endif

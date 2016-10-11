//-*-C++-*-
/*!
  \file BaseNode.h
  \brief General node field class. It stores the nodal id, the number of 
  degrees of the node, vector<string> which stores list of properties
  which are stores in the node and Eigen::VectorXd which holds the numerical
  values corresponding to the property.
 */
#ifndef _BaseNode_
#define _BaseNode_

#include "voom.h"

namespace voom
{
  //! A node base class holds id, dof, property values
  class BaseNode  {
  public:
    //! Default constructor
    BaseNode() : _id(0), _dof(0) {};
    //! Construct from id
    BaseNode(unsigned int id) : _id(id), _dof(0) {};
    //! Construct from id and dof
    BaseNode(unsigned int id, unsigned int dof): _id(id), _dof(dof) {};
    //! accessors and mutators
    const unsigned int & getId() const {return _id;}

    //! Interfaces
    int getDoF() { return _dof; }

    /*!
      Update to be used by solvers. Updates all the internal property for the
      node.
    */
    virtual void linearizedUpdate(const VectorXd& data) = 0;

    //! Get particular value in the property list.
    virtual void linearizedUpdate(unsigned int dof, Real data) = 0;

    //! Output interface
    virtual vector<string> getPropertyList() = 0;

    //! Get data associated with a property.
    virtual void getProperty(string propertyName, Real** data, uint& length) {};
    
  protected:
    //! Nodal ID
    unsigned int   _id; 

    //! Degree of freedom associated with the node
    unsigned int   _dof;
  }; // class BaseNode
} // namespace voom

#endif // _BaseNode_

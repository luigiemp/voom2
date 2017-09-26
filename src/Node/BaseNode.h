//-*-C++-*-
/*!
  \file BaseNode.h
  \brief General node field class. It stores the:
  - node id;
  - the degrees of freedom of the node;
  - the position of the node;
  - the field values at the nodes in Eigen::VectorXd
 */
#ifndef _BaseNode_
#define _BaseNode_

#include "voom.h"

namespace voom
{
  //! A node base class holds id, dof, property values
  class BaseNode  {
  private:
    BaseNode() {}; // No node should be declared without id and corresponding dof
  public:
    //! Constructor from id, and dof
    BaseNode(int ID, vector<int > DOF) : _id(ID), _dof(DOF) {};
    //! Constructor from id, dof, and position
    BaseNode(int ID, vector<int > DOF, Vector3d X) : _id(ID), _dof(DOF), _X(X) {};
   
    //! accessors
    int getId() {return _id;}
    const vector<int > & getDoF() {return _dof;}
    const Vector3d & getX() {return _X;};

    //! Set reference position
    void setX(Vector3d & X) { 
      assert(_X.size() == X.size() );
      _X = X;
    }
    
  protected:
    //! Nodal ID
    int _id; 
    //! Node global DoF
    vector<int > _dof;
    //! Node reference position
    Vector3d _X;
    
  }; // class BaseNode
} // namespace voom

#endif // _BaseNode_

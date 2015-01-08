#include "MechanicsNode.h"

namespace voom{
  //! Set displacement
  void MechanicsNode::setu(const VectorXd& u) {
    assert( _u.size() == u.size() );
    _u = u;
  }
  
  //! Set particular displacement component
  void MechanicsNode::setu(unsigned int id, Real value) {
    assert( id < _u.size() );
    _u[id] = value;
  }
  
  //! Get property associated with a node
  vector<string> MechanicsNode::getPropertyList(){
    vector<string> property(2); 
    property[0] = "POSITION";
    property[1] = "DISPLACEMENT";
    return property;
  }

  void MechanicsNode::getProperty(string propertyName,Real** data,uint& length){
    delete *data; *data = NULL;
    Real * value;
    if ( propertyName == "POSITION") {
      value = new Real[_X.size()];
      for(unsigned int i = 0; i < _X.size(); i++) value[i] = _X(i);
      length = _X.size();
    } else if ( propertyName == "DISPLACEMENT" ) {
      value = new Real[_X.size()];
      for(unsigned int i = 0; i < _u.size(); i++) value[i] = _u(i);
      length = _u.size();
    } else {
      cerr << "** ERROR: Unknown property name " << propertyName 
	   << "provided\n";
	length = 0;
    }
    *data = value;
  }

  // Linearized Update - update displacement
  void MechanicsNode::linearizedUpdate(const VectorXd& data) {
    assert(data.size() == _u.size() );
    _u += data;
  }

  //! Get particular value in the property list.
  void MechanicsNode::linearizedUpdate(unsigned int dof, Real data) {
    assert( dof < _u.size() );
    _u(dof) += data;
  }
} // namespace

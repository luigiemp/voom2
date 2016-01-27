// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Dirichlet_h__)
#define __Dirichlet_h__

#include <blitz/array.h>
#include <vector>
#include "NodeBase.h"
#include "Node.h"
#include "Constraint.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"

namespace voom
{

  //! Predictor-corrector Constraint class for a Dirichlet boundary condition
  template<int N>
  class PeriodicTie : public Constraint
  {
  public:

    typedef DeformationNode<N> NodeType;
    typedef typename NodeType::Point ForceVN;

    //! Constructor
    PeriodicTie( NodeType *  master, NodeType * slave, PeriodicBox * box, bool relaxLength=false )
      : _master(master), 
	_force_master(0.0), _force_slaves(0.0),
	_box(box)
    {
      assert(_master!=0);
      assert(_box!=0);

      _slaves.push_back(slave);

      if(!relaxLength) {
        _restLength = 0.0;
      }
      else {
	ForceVN x(0.0);
	x = _master->point() - slave->point();
	_box->mapDistance(x);
	_restLength = norm2(x);
	if(_restLength < 1.0e-3) _restLength = 0.0;
      }
    }

    PeriodicTie( NodeType *  master, PeriodicBox * box )
      : _master(master), 
	_force_master(0.0), _force_slaves(0.0),
	_box(box)
    {
      assert(_master!=0);
      assert(_box!=0);

    }

    
    ~PeriodicTie() {}

    virtual void addSlave(NodeType * newSlave) { _slaves.push_back(newSlave); }
 
    //! Adjust DOF to meet constraint; assume the initial 
    virtual void predict();

    //! Correct out of balance force conjugate to DOF
    virtual void correct();

    //void zeroOutSlaves(bool f1, bool f2) {
     // for(typename vector<NodeType*>::iterator sls=_slaves.begin(); sls!=_slaves.end(); sls++) {
        //if(f1) for(int df=0; df<(*sls)->dof(); df++) (*sls)->setForce(df,0.0);
	//if(f2) for(int df=0; df<(*sls)->dof(); df++) (*sls)->setStiffness(df,0.0);
      //}   
    //}

    NodeType * masterNode() { return _master; }

    std::vector<NodeType *> & slaveNodes() { return _slaves; }

    const double getForceMaster(int i) const {
      assert( 0 <= i && i < N );
      return _force_master(i);
    }

    const double getForceSlave(int i) const {
      assert( 0 <= i && i < N );
      return _force_slaves(i);
    }

    const ForceVN & forceMaster() const { return _force_master;}

    const ForceVN & forceSlave() const { return _force_slaves;}

    void setBox(PeriodicBox * pb) {
      _box = pb;
    }

  private:

    //! The master node to which the slave will be tied
    NodeType * _master;

    //! The slave nodes to be tied to the master
    std::vector<NodeType *> _slaves;

    // pre-corrected force on the master
    ForceVN _force_master;

    // pre-corrected force on the slaves
    ForceVN _force_slaves;

    // Periodic BC object
    PeriodicBox * _box;

    double _restLength;
  };
  
  template<int N>
  void PeriodicTie<N>::predict() {
    // set slave points to master point //
    for(typename vector<NodeType*>::iterator sls=_slaves.begin(); sls!=_slaves.end(); sls++) {
      // distance vector between two nodes
      ForceVN x(0.0);
      x = _master->point() - (*sls)->point();

      // mod the distance vector respecting periodic BCs
      _box->mapDistance(x); 

      // move the slave along the modded distance vector
      //double xnorm = norm2(x);
      //if(_restLength > 1.0e-6) {
      //  x *= (1.0 - (_restLength/xnorm));
      //}     
      x += (*sls)->point(); 
      (*sls)->setPoint(x);
      for(int df=0; df<(*sls)->dof(); df++) (*sls)->setForce(df,0.0);
      for(int df=0; df<(*sls)->dof(); df++) (*sls)->setStiffness(df,0.0);
    }
    // check to make sure we've moved correctly //
    //if(_restLength > 1.0e-6 && xnorm > 1.0e-6) {
    //  x = _master->point() - _slave->point();
    //  _box->mapDistance(x);
    //  xnorm = norm2(x);
    //  if(abs(xnorm-_restLength) > 1.0e-6) std::cout << "Error: slave point not moved correctly." << std::endl;
    //}

    return;
  }

  template<int N>
  void PeriodicTie<N>::correct() {
    // store pre-corrected forces
    _force_master = _master->force();
    _force_slaves = 0.0;
    for(typename vector<NodeType*>::iterator sls=_slaves.begin(); sls!=_slaves.end(); sls++) {
      ForceVN oldForce;
      oldForce = (*sls)->force();
      _force_slaves += oldForce;
      oldForce *= -1.0;
      (*sls)->updateForce(oldForce);
    }
    

    //ForceVN x(0.0);
    //x = _master->point() - _slave->point();
    //_box->mapDistance(x);
    //double xnorm = norm2(x);

    //if(_restLength > 1.0e-6) {
    //  x /= xnorm;
    //  _force_master = x*dot(_force_master,x);
    //  _force_slave = x*dot(_force_slave,x);
    //}

    // assemble: sum forces together
    _master->updateForce(_force_slaves);
    ForceVN newForce;
    newForce = _force_slaves + _force_master;
    for(typename vector<NodeType*>::iterator sls=_slaves.begin(); sls!=_slaves.end(); sls++) {
      (*sls)->updateForce(newForce); 
    }


    return;
  }


} // namespace voom

#endif // __PeriodicTie_h__

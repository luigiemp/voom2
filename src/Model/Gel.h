// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file SemiflexibleGel.h

  \brief SemiflexibleGel is a concrete class derived from Body, implementing
  the concept of a collection of cross-linked semiflexible polymers (i.e., beams)

*/

#if !defined(__Gel_h__)
#define __Gel_h__


#include<blitz/array.h>
#include<random/exponential.h>
#include<vector>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <string>
#include <fstream>
#include "Body.h"

#include "Spring.h"
#include "EntropicSpring.h"
#include "AngleSpring.h"
#include "BrownianRod.h"

#include "Constraint.h"
#include "voom.h"
#include "Node.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "Crosslink.h"
#include "Motor.h"
#include "Grid.h"
#include "TwoBodyPotential.h"
#include "IntersectionFinder.h"
#include "PeriodicTie.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
#include "PinchForce.h"
#include "NematicProbTable.h"
#include "SemiflexibleInput.h"


//#include "AffinityElement.h"
//#include "AffinityMeasure.h"

#include "ViscousRegularizer.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{


  /*!  
    Concrete class for a semiflexible gel
  */
  template<int N>
  class SemiflexibleGel : public Body
  {
  public:
    
    // typedefs
    typedef BrownianNode<N> DefNode;
    typedef typename std::vector< DefNode* > DefNodeContainer;
    typedef typename DefNodeContainer::iterator DefNodeIterator;
    typedef typename DefNodeContainer::const_iterator ConstDefNodeIterator;

    typedef DeformationNode<N> BaseDefNode;

    typedef Spring<N> Bond;
    typedef std::vector< Bond* > BondContainer;
    typedef typename BondContainer::iterator BondIterator;
    typedef typename BondContainer::const_iterator ConstBondIterator;
		
    typedef AngleSpring<N> Angle;
    typedef std::vector< Angle* > AngleContainer;
    typedef typename AngleContainer::iterator AngleIterator;
    typedef typename AngleContainer::const_iterator ConstAngleIterator;

    typedef BrownianRod<N> Rod;
    typedef std::vector< Rod* > RodContainer;
    typedef typename RodContainer::iterator RodIterator;
    typedef typename RodContainer::const_iterator ConstRodIterator;
		
    typedef std::vector< Constraint* > ConstraintContainer;
    typedef typename ConstraintContainer::iterator ConstraintIterator;
    typedef typename ConstraintContainer::const_iterator ConstConstraintIterator;
	
    //Mo	
    typedef Crosslink<N> Clink;
    typedef std::vector< Clink* > CrosslinkContainer;
    typedef typename CrosslinkContainer::iterator CrosslinkIterator;
    typedef typename CrosslinkContainer::const_iterator ConstCrosslinkIterator;
  
    typedef Motor<N> MolMot;
    typedef std::vector< MolMot* > MotorContainer;
    typedef typename MotorContainer::iterator MotorIterator;
    typedef typename MotorContainer::const_iterator ConstMotorIterator;

    typedef PinchForce<N> Pinch;
    typedef std::vector< Pinch* > PinchContainer;
    typedef typename PinchContainer::iterator PinchIterator;
    typedef typename PinchContainer::const_iterator ConstPinchIterator;
    
    typedef typename std::set<DefNode*> PinchNodeSet;
    typedef typename PinchNodeSet::iterator PinchNodeIterator;
		
    typedef tvmet::Vector<double,N> VectorND;
    typedef tvmet::Matrix<double,N,N> TensorND;

    typedef typename std::map< DefNode*, DefNode* > CrosslinkNodeMap;
    typedef typename CrosslinkNodeMap::iterator CLNMiter;

    typedef typename std::map< double, int > CrosslinkDistFreq;
    typedef typename CrosslinkDistFreq::iterator ClDFIter;
    typedef std::pair< int, int > intPair;
    typedef typename std::multimap< DefNode*, int >::iterator SNiter;

    typedef std::pair< double, double > doublePair;
    typedef typename std::vector< doublePair > doublePairContainer;
    
    typedef std::pair< double, doublePair > doublePairWErrors;
    typedef typename std::vector< doublePairWErrors > doublePairWErrorsContainer;

    typedef typename std::map< std::string, std::string > PropertyList;
    typedef typename PropertyList::iterator PropertyIterator;
    typedef typename PropertyList::const_iterator ConstPropertyIterator;

    typedef typename std::vector< std::pair<VectorND,TensorND> > StrainField;
    typedef typename StrainField::iterator StrainFieldIterator;
	
    struct Filament {

      Filament(const DefNodeContainer & n, double kappa, double mu, double viscosity, double kT, double dt, double minLength);

      ~Filament();

      const VectorND & point();

      DefNodeContainer nodes;
      BondContainer    bonds;
      AngleContainer   angles;
      RodContainer     rods;
      std::vector< double > clinks;

      VectorND pt;
    };

    typedef std::vector< Filament* > FilamentContainer;
    typedef typename FilamentContainer::iterator FilamentIterator;
    typedef typename FilamentContainer::const_iterator ConstFilamentIterator;
    
    typedef Grid<Filament,Filament,N> FilGrid;
    typedef Grid<DefNode,BaseDefNode,N> NodeGrid;
//     typedef Grid<AffinityElement,AffinityElement,N> AffElementGrid;
    typedef Grid<Bond,Bond,N> BondGrid;

    struct TempCrosslink {
      //VectorND location;
      int baseFil;
      std::map< int, VectorND > otherFils;
      bool active;
      DefNode* clNode;
      PeriodicTie<N> * ptie;
    };
 
    struct TempFilament {
      VectorND start;
      VectorND end;
      std::vector< std::pair<VectorND,VectorND> > filSegs;
      std::multimap<double,TempCrosslink *> crossFils;
    };

    struct TempBox {
      std::vector<TempCrosslink *> boxCLs;
    };

    typedef typename std::vector<TempFilament *> TempFilamentContainer;

    struct Segment {
      int filID;
      int nodeID1,nodeID2;

      double length;

      double critStrain;

      double buckleStrain;
      
      double bendE,stretchE;
      
      //std::vector<Segment*> neighbors;

      VectorND pos;
      VectorND pt;

      const VectorND & position() { return pos; }

      const VectorND & point() { return pt; }
      
    };

    struct buckleComp {
      bool operator() (const Segment* a, const Segment* b) {
	if(a->critStrain < b->critStrain) return true;
	else return false;
      }

    };

    //! Default Constructor
    SemiflexibleGel() {
      _output=paraview;
    }

    //! A constructor that uses the SemiflexibleInput class to bring in parameters
    SemiflexibleGel(DefNodeContainer & dNodes, SemiflexibleInput * pInp );

    //! Constructor that implements adaptive meshing
    SemiflexibleGel(DefNodeContainer & dNodes, PeriodicBox * box, double filDens, double filLength, const string & bondType, bool cutOffEnds, double minLength, const PropertyList & properties);

    //! virtual destructor
    virtual ~SemiflexibleGel() { 
      for(int i=0; i<_filaments.size(); i++) delete(_filaments[i]);
    }

    //! reset nodal positions to initial positions, reset energy to 0
    void resetGel() {
      for(FilamentIterator f=_filaments.begin(); f!=_filaments.end(); f++) {
	for(DefNodeIterator n=f->nodes.begin(); n!=f->nodes.end(); n++) {
	  (*n)->setPoint((*n)->position());
	}
      }
      _energy = 0.0;
    }

    void updateNodalPoints(std::string & strainedGelFN);
    
    void removePrestress();
    
    int removeCrosslinks(TempFilamentContainer& tmpFils, double minLength);
    
    int collapseCrosslinks(TempFilamentContainer & tmpFils, double minLength);

    //! Do mechanics on Filaments
    void compute( bool f0, bool f1, bool f2 );

    void printParaview(const std::string fileName) const;
    
    void storeGel(std::string fileName);

    void storeSparseGel(std::string fileName, TempFilamentContainer & tmpFils);

    void printforAbaqus(std::string fileName, double mu, double lB, double minLength, int method);

    void addFilament(const DefNodeContainer & n, double kappa, double mu, double viscosity, double kT, double dt, double minLength) {
      Filament * f = new Filament(n,kappa,mu,viscosity,kT,dt,minLength);
      _filaments.push_back(f);

    }

    const FilamentContainer & filaments() const { return _filaments; }

    const Filament * filament(int a) const { return _filaments[a]; }

    Filament * filament(int a) {return _filaments[a]; }

    void affineShearX(double shear) {      
      for(FilamentIterator fi=_filaments.begin(); fi!=_filaments.end(); fi++) {
	for(DefNodeIterator dni=(*fi)->nodes.begin(); dni!=(*fi)->nodes.end(); dni++) {
	  VectorND nodePos;
	  nodePos = (*dni)->position();
	  nodePos[0] += shear*(nodePos[1]-(_box->size()[1]/2.0));
	  (*dni)->setPoint(nodePos);
	}
      }
    }

    void printConnectionData(std::string & connFileName);

    void printStiffenedSegments(std::string & stiffSegFile, double crit);
    
    void stiffenRandomSegments(double frac, double mult);

    //void reportStiffened

    void moveCLNodes(Filament * f);

    void setGrid(FilGrid * g) { _grid = g; }

    void addTwoBodyPotential(TwoBodyPotential * tbp) { _tbp.push_back(tbp); }

    void addCrosslink( Clink * c ) { _crosslinks.push_back( c ); }

    const CrosslinkContainer & crosslinks() const { return _crosslinks; }

    const ConstraintContainer & constraints() const { return _constraints; }

    void attachCrosslink(Clink * cl) { // eventually, find nearest filaments and attempt to attach; for now, just take two filaments and attach at intersection point //
      attachCrosslink(cl, *_filaments.begin(),*(_filaments.begin()+1));
    }
  
    void attachCrosslink(Clink * cl, const Filament * f1, const Filament * f2);

    bool attachCrosslink(Filament * f1, Filament * f2, double kcl, double relax);

    void addMotor(MolMot * mot) { _motors.push_back(mot); }

    void addMotor(VectorND & p, double k, double d0) {
      MolMot * mot = new Motor<N>(p,k,d0);
      _motors.push_back(mot);
    }

    void addMotor(VectorND & p, double k) {
      MolMot * mot = new Motor<N>(p,k);
      _motors.push_back(mot);
    }

    void attachMotor(MolMot * mot) { // eventually, find nearest filaments and attempt to attach; for now, just take two filaments and attach at intersection point //
      attachMotor(mot, *_filaments.begin(),*(_filaments.begin()+1));
    }

    void attachMotor(MolMot * mot, const Filament * f1, const Filament * f2);

    void attachMotor(const Filament * f1, const Filament * f2);

    const MotorContainer & motors(){ return _motors; }

    void addConstraint( Constraint * c ) { _constraints.push_back( c ); }

    void addPinches(double pinchDensity, double a, double tol, double f0, bool sameFilament);
    
    void addPinches(int nPinches, double a, double tol, double f0, bool sameFilament);

    void addPinch(DefNode * n1, DefNode * n2, double f0);

    void addPinch(double a, double tol, double f0);

    Pinch* pinch(int i) { return _pinches[i]; }
    
    std::vector<Pinch*> & pinches() { return _pinches; }

    //void addPinch(double f0, bool springy, DefNodeContainer & dNodes, double kBond, double kAngle, double visc, double kT, double dt, double kcl);

    void turnOffPinches() {
      for(PinchIterator pi = _pinches.begin(); pi!=_pinches.end(); pi++) {
	(*pi)->turnOff();
      }
    }
    
    void turnOnPinches() {
      for(PinchIterator pi = _pinches.begin(); pi!=_pinches.end(); pi++) {
	(*pi)->turnOn();
      }
    }

    void turnOnPinches(double f0) {
      for(PinchIterator pi = _pinches.begin(); pi!=_pinches.end(); pi++) {
	(*pi)->turnOn();
	(*pi)->setPinchF(f0);
      }
    }

    void setBox( PeriodicBox * box ) { _box = box; }

    PeriodicBox * box() { return _box; }

    double getMeanCLsep() { return _meanCLsep; }
    
    double getNematicOP() { return _nematicOP; }

    VectorND & getNemDirector() { return _nemDirector; }

    double getMeanFilLen() { return _meanFilLen; }

    double crosslinkenergy();

    double filenergy() { return bendingenergy() + stretchingenergy(); }

    double bendingenergy();

    double stretchingenergy();

    double motorenergy();
    
    double pinchenergy();

    double parallelenergy();
    
    double perpenergy();

    CrosslinkDistFreq & getCrossDistro() {
      CrosslinkDistFreq & cdf = _crossDistFreqs;
      return cdf;
    }

    std::map< double, int > & getLengthDistro() {
      std::map< double, int > & ld = _filLenFreqs;
      return ld;
    }

    std::map< double, int > & getNematicDistro() {
      std::map< double, int > & nd = _nematicFreqs;
      return nd;
    }

    std::set<DefNode *> & pinchNodes() { return _pinchNodes; }

    std::map< doublePair, doublePair > getAngularEnergyDistro();
    
    std::multimap< double, std::vector<double> > getDensityEnergyDistro(double scale);

    void cutOffEndsandCCD(double kcl, DefNodeContainer & dNodes);
    
    void computeCrossDistro(double kcl);

    void computeNematicDistro(double nemAngle);

    void computeFilLenDistro();

    void computeNematicOP();

    std::vector< std::pair<double,double> > computeNemCorrelations(double minSep, double step, double tol);

    void printAngles(std::string & angleFile);

    bool isSlave(DefNode* node) {
      if(_crossNodeMap.find(node) == _crossNodeMap.end()) return false;
      else if(_crossNodeMap[node] == node) return false;
      else return true;
    }

    bool isMaster(DefNode* node) {
      if(_crossNodeMap.find(node) == _crossNodeMap.end()) return false;
      else if(_crossNodeMap[node] != node) return false;
      else return true;
    }

    bool areLinked(DefNode* node1, DefNode* node2) {
      if(_crossNodeMap.find(node1) == _crossNodeMap.end() || _crossNodeMap.find(node2) == _crossNodeMap.end()) return false;
      else {
	DefNode* n1mast = _crossNodeMap.find(node1)->second;
	DefNode* n2mast = _crossNodeMap.find(node2)->second;
	if(node1 == n2mast || node2 == n1mast) return true;
	else return false;
      }
    }
    
    bool checkCrosslinks();

    void checkParallelForces();

    void printCLFilDist();
    
    void printInitialBends();

    void printBigBends();

    doublePairContainer affineMeasurement(double minLength, double stepSize, double maxLength, double shear, std::string measureType);

    doublePairWErrorsContainer affineMeasurementHeadLevine(double minLength, double stepSize, double maxLength, double shear, double largestFil, bool getAngularDist);

    doublePairWErrorsContainer affineMeasurementHeadLevineInterpolated(double minLength, double stepSize, double maxLength, double shear, double largestFil);

    void affineBoxesMeasurement(VectorND & boxsize, double pairDist, double shear, double largestFil, std::string fileName, bool doCorrTest);

    doublePairContainer energyCorrelationFunction(double boxs, double maxlen);

    void computeNonaffinityLengthDensityCorrelation(double boxsize, double maxdist, double maxFL, double shear);

    doublePairContainer computeCorrelationFunction(std::map<DefNode *, double> & dataPts, double minSep, double maxSep, double step, double tol);

    doublePairContainer computeCrossCorrelationFunction(std::map< DefNode *, std::pair<double,double> > & dataPts, double minSep, double maxSep, double step, double tol);

    double computeCorrelation(std::vector<double> & dat1, std::vector<double> & dat2);

    void computeCrossCorrelations(double len, double shear, std::string & fileName);

    double computeBucklingEnergy(double shear, double kap, double mu);

    void computeBucklingMap(VectorND & boxsize, double shear, double bendfrac, std::string fileName, bool doCorrTest);

    void buckleOPCalc(double gridSize, std::map<double,std::string> strainedGelFiles, bool doCorrs);

    void cooperativeBuckleMeasure(double gridSize, double bendfrac, double l_B, std::map<double,std::string> strainedGelFiles, bool doCorrs);

    void setViscReg(double visc) {
      std::vector<NodeBase*> allNodes;
      for(FilamentIterator fi=_filaments.begin(); fi!=_filaments.end(); fi++) {
	for(DefNodeIterator ni=(*fi)->nodes.begin(); ni!=(*fi)->nodes.end(); ni++) {
	  allNodes.push_back(*ni);
	}
      }
      
      if(_viscReg!=0) delete _viscReg;
      
      _viscReg = new ViscousRegularizer(allNodes,visc);
    }

    ViscousRegularizer* viscReg() { return _viscReg; }

    void makeFinalPosMap(std::vector<VectorND> & finalPos) {
      _finalPosMap.clear();
      int nNodesTotal = 0;
      int nFils = _filaments.size();
      for(int fn=0; fn<nFils; fn++) {
	int nNodes = filament(fn)->nodes.size();
	for(int nn=0; nn<nNodes; nn++) {
	  nNodesTotal++;
	}
      }

      assert(nNodesTotal == finalPos.size());

      nNodesTotal = 0;
      for(int fn=0; fn<nFils; fn++) {
	int nNodes = filament(fn)->nodes.size();
	for(int nn=0; nn<nNodes; nn++) {
	  _finalPosMap.insert(pair<DefNode*,VectorND>(filament(fn)->nodes[nn],finalPos[nNodesTotal]));
	  nNodesTotal++;
	}
      }     
      
    }

    std::map<DefNode*,VectorND> & finalPosMap() { return _finalPosMap; }
    

  private:

    //! Filaments
    FilamentContainer 	_filaments;		

    //! Crosslinks
    CrosslinkContainer 	_crosslinks;		

    //! Constraints
    ConstraintContainer _constraints;

    //! Motors
    MotorContainer _motors;

    PinchContainer _pinches;

    PeriodicBox * _box;

    CrosslinkNodeMap _crossNodeMap;

    PinchNodeSet _pinchNodes;

    std::map< DefNode*, int > _nSlavesMap;

    CrosslinkDistFreq _crossDistFreqs;

    std::set<DefNode*> _crosslinkNodes;

    std::map< double, int > _filLenFreqs;

    std::map< double, int > _nematicFreqs;

    double _meanCLsep;
    
    double _nematicOP;
    
    VectorND _nemDirector;

    double _meanFilLen;

    std::vector<TwoBodyPotential*> _tbp;
    
    FilGrid * _grid;

    ViscousRegularizer* _viscReg;

    std::map<DefNode*, VectorND> _finalPosMap;

  };  
} // namespace voom

#include "Gel.icc"

#endif // __Gel_h__

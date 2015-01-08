// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2010 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*!
  \file voom.h
  \brief Inlcudes and typedef Real that should be available to all VOOM classes.
 */

#ifndef _VOOM_H_
#define _VOOM_H_

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>
#include <string>

#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


namespace voom {
  // Split a string. Used by text parsing routines
  inline vector<string> splitString(string line, string delim) {
    char str[ line.size() + 1]; strcpy( str, line.c_str() );
    char* pch = strtok (str, delim.c_str());
    string word; vector<string> words;
    while (pch != NULL) {
      words.push_back( string(pch) );
      pch = strtok(NULL, delim.c_str());
    }
    return words;
  }

  // Change double to float for single precision
  typedef double Real;

  //! Enumeration for requesting computed results 
  enum ComputeRequest {
    NOTHING   	= 0,
    ENERGY    	= 1,
    FORCE     	= 2,
    STIFFNESS 	= 4
  };

}; // namespace voom

#endif // _VOOM_H_

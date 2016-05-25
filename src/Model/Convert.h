//-*-C++-*-
#ifndef __Convert_h__
#define __Convert_h__

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <typeinfo>
#include <cstdlib>

namespace voom{

  class Convert {

  public:
    
    template <typename T>
      static std::string T_to_string(T const &val) 
      {
	std::ostringstream ostr;
	ostr << val;

	return ostr.str();
      }
		
    template <typename T>
      static T string_to_T(std::string const &val) 
      {
	std::istringstream istr(val);
	T returnVal;
	if (!(istr >> returnVal))
	  cout << "Not a valid " + (std::string)typeid(T).name() + " received!\n";
	return returnVal;
      }

  };

}

#endif

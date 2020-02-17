/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_PARAM_DESCRIPTION_H
#define CS_PARAM_DESCRIPTION_H

#include <cstdlib>
#include <cstdio>
#include <string>

namespace cseis_system {

 static int const TYPE_NONE = -1;

/**
* General description of Cseis parameter definition
* Stores information such as name and a short help description for Cseis user parameter definitions
*
* @author Bjorn Olofsson
* @date   2007
*/
class csParamDescription {
public:
  csParamDescription( char const* name, char const* desc, char const* descExtra, int type ) :
    myName(name), myDesc(desc), myDescExtra(descExtra), myType(type) {}
  csParamDescription(csParamDescription const& obj) : myName(obj.myName), myDesc(obj.myDesc), myDescExtra(obj.myDescExtra), myType(obj.myType) {}
  ~csParamDescription() {}
  /**
   * @return Name of parameter
   */
  inline char const* name() const { return myName; }
  /**
   * @return Free-form text description of parameter
   */
  inline char const* desc() const { return myDesc; }
  /**
   * @return Additional text description of parameter providijng more details
   */
  inline char const* descExtra() const { return myDescExtra; }
  /**
   * @return Parameter type: 
   */
  inline int type() const { return myType; }
  /**
   * Dump parameter to standard output
   */
  void dump() const { std::printf(" Name: '%s', Desc: '%s', DescExtra: '%s', Type: %d\n", myName, myDesc, myDescExtra, myType); }

private:
  csParamDescription();
  csParamDescription& operator=( csParamDescription const& obj );
  char const* myName;
  char const* myDesc;
  char const* myDescExtra;
  int         myType;
};

}
#endif


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_HEADER_INFO_H
#define CS_HEADER_INFO_H

#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {

/**
* Header info.
*
* Defines one trace header
*
* @author Bjorn Olofsson
* @date 2007
*/
class csHeaderInfo {
public:
  csHeaderInfo() :
    type(TYPE_UNKNOWN), nElements(0), name(""), description("") {}
  csHeaderInfo( type_t theType, std::string const& theName, std::string const& theDesc ) :
    type(theType), nElements(1), name(theName), description(theDesc) {}
  csHeaderInfo( type_t theType, std::string const& theName, std::string const& theDesc, int theElements ) :
    type(theType), nElements(theElements), name(theName), description(theDesc) {}
  ~csHeaderInfo() {}

  void set( csHeaderInfo const& info ) {
    type = info.type;
    name = info.name;
    description = info.description;
  }
  void set( type_t theType, std::string const& theName, std::string const& theDesc ) {
    type = theType;
    name = theName;
    description = theDesc;
  }
  void set( type_t theType, std::string const& theName, std::string const& theDesc, int theElements ) {
    set( theType, theName, theDesc );
    nElements = theElements;
  }

  type_t type;
  int    nElements;
  std::string name;
  std::string description;
private:
  csHeaderInfo( csHeaderInfo const& obj ) {
    name = obj.name;
    type = obj.type;
    description = obj.description;
    nElements = obj.nElements;
  }

  /*  bool operator==( const csHeaderInfo& obj ) const;
      bool operator!=( const csHeaderInfo& obj ) const;
      bool operator>=( const csHeaderInfo& obj ) const;
      bool operator<=( const csHeaderInfo& obj ) const;
      bool operator>( const csHeaderInfo& obj ) const;
      bool operator<( const csHeaderInfo& obj ) const;
  */
};

} // namespace

#endif


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGY_HEADER_INFO_H
#define CS_SEGY_HEADER_INFO_H

#include "geolib_defines.h"
#include <string>

namespace cseis_geolib {

/**
 * SEG-Y trace header definition
 * @author Bjorn Olofsson
 * @date 2006
 */
class csSegyHeaderInfo {
public:
 csSegyHeaderInfo( int theByteLoc, int theByteSize, cseis_geolib::type_t theInType, cseis_geolib::type_t theOutType, std::string const& theName, std::string const& theDesc ) :
  byteLoc(theByteLoc),
    byteSize(theByteSize),
    inType(theInType),
    outType(theOutType),
    name(theName),
    description(theDesc) {}
 csSegyHeaderInfo( csSegyHeaderInfo const& obj ) :
    byteLoc(obj.byteLoc),
    byteSize(obj.byteSize),
    inType(obj.inType),
    outType(obj.outType),
    name(obj.name),
    description(obj.description) {}
public:
  int byteLoc;
  /// Number of bytes
  int byteSize;
  /// Type of SEGY header in data file
  cseis_geolib::type_t inType;
  /// Type of SEGY header after being read in
  cseis_geolib::type_t outType;
  std::string const name;
  std::string const description;
private:
  //  csSegyHeaderInfo( csSegyHeaderInfo const& obj );
};

} // end namespace
#endif


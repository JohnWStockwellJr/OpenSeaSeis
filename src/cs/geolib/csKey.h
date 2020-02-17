/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_KEY_H
#define CS_KEY_H

#include <string>

namespace cseis_geolib {

/**
* Helper class containing one string value.
*
* This class was created to be used in conjunction with csCollection, as a representation of one
* string value which can perform a comparison using '==' operator, required by csCollection for the
* use of the 'contain' method
*
* @author Bjorn Olofsson
* @date 2006
*/
class csKey {
public:
  csKey() : value("") {}
  csKey( std::string theValue ) : value( theValue ) {}
  csKey( char const* theValue ) : value( theValue ) {}
  csKey& operator=( csKey const& obj ) {
    value = obj.getValue();
    return *this;
  }
  bool operator==( csKey const& obj ) {
    return( !value.compare( obj.getValue() ) );
  }
  std::string getValue() {
    return value;
  }
  std::string const getValue() const {
    return value;
  }
private:
  std::string value;
};

} // namespace

#endif


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_EXCEPTION_H
#define CS_EXCEPTION_H

#include <string>

namespace cseis_geolib {

/**
 * Simple exception class
 *
 * @author Bjorn Olofsson
 * @date 2005
 */
class csException {
public:
  csException();
  csException( std::string const& message );
  csException( int dummy, char const* text );
  csException( char const* text, ... );
  csException( csException const& obj );
  ~csException();
  const char* getMessage();
protected:
  std::string myMessage;
};

} // namespace

#endif


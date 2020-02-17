/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_I_READER_H
#define CS_I_READER_H

#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {

  class csFlexHeader;

/**
 * Pure virtual reader class = Java-style 'interface'.
 * Defines interface for seismic data reader.
 * To be implemented by specific seisimc reader classes.
 */
class csIReader {
 public:
  csIReader();
  virtual ~csIReader();

  virtual int numTraces() const = 0;
  virtual int getCurrentTraceIndex() const = 0;
  virtual bool moveToTrace( int traceIndex ) = 0;
  virtual bool setHeaderToPeek( std::string const& headerName ) = 0;
  virtual bool setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType ) = 0;
  virtual bool peekHeaderValue( cseis_geolib::csFlexHeader* value, int traceIndex ) = 0;
  //  virtual bool peekHeaderValue( cseis_geolib::csFlexNumber* value, int traceIndex = -1 ) = 0;
};

} // end namespace

#endif

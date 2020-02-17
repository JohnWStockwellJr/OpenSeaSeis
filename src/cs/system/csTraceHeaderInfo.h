/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_HEADER_INFO_H
#define CS_TRACE_HEADER_INFO_H

#include <string>

namespace cseis_system {

/**
* Trace header info.
* Defines one trace header.
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTraceHeaderInfo {
public:
  csTraceHeaderInfo( char theType, std::string const& theName, std::string const& theDesc ) :
    type(theType), nElements(1), name(theName), description(theDesc) {}
  csTraceHeaderInfo( char theType, int theElements, std::string const& theName, std::string const& theDesc ) :
    type(theType), nElements((short int)theElements), name(theName), description(theDesc) {}
  ~csTraceHeaderInfo() {}
  char  type;
  short nElements;
  std::string name;
  std::string description;
private:
  csTraceHeaderInfo(csTraceHeaderInfo const& obj);
};

} // namespace

#endif


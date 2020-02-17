/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_HEADER_INFO_POOL_H
#define CS_TRACE_HEADER_INFO_POOL_H

#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {
  template<typename T> class csVector;
}

namespace cseis_system {

class csTraceHeaderInfo;

/**
* Trace header info pool
* This pool keeps all trace header info objects that are in use.
* Info objects are never deleted until the program terminates.
*
* @author Bjorn Olofsson
* @date   2007
*/
class csTraceHeaderInfoPool {
public:
  csTraceHeaderInfoPool();
  ~csTraceHeaderInfoPool();
  csTraceHeaderInfo const* createTraceHeaderInfo( cseis_geolib::type_t type, std::string const& name, std::string const& description, int nElements = 1 );
//  void dump();
private:
  cseis_geolib::csVector<csTraceHeaderInfo*>* myInfoList;
};

} // namespace
#endif


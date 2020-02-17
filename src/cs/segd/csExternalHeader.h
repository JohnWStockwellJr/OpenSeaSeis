/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_EXTERNAL_HEADER_H
#define CS_EXTERNAL_HEADER_H

#include "csSegdHeader.h"

namespace cseis_segd {

class csNavHeader;
class csGCS90Header;
class csNavInterface;

class csExternalHeader : public csBaseExternalHeader {
public:
  csExternalHeader( int numBytes, int navSystemID, int navInterfaceID );
  virtual ~csExternalHeader();
  virtual void extractHeaders( byte const* buffer );
  virtual void extractHeaders( byte const* buffer, int totalNumBytes );
  virtual void dump( std::ostream& os );
  virtual void setCommonRecordHeaders( commonRecordHeaderStruct& comRecHdr );
 private:
  csNavHeader*    myNavHeader;
  csGCS90Header*  myGCS90Header;
  csNavInterface* myNavInterface;
};

} // end namespace
#endif


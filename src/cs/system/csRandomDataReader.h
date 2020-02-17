/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_RANDOM_READER_H
#define CS_RANDOM_READER_H

#include <cstdio>
#include <string>
#include "geolib_defines.h"
#include "csIReader.h"

namespace cseis_geolib {
  class csFlexHeader;
  class csSegyReader;
  class csIReader;
}
namespace cseis_io {
  class csRSFReader;
}
namespace cseis_system {
  class csTraceHeaderDef;
  class csSuperHeader;
  class csSeismicReader;
  class csTraceHeader;
}

namespace cseis_system {

/**
 * Random line reader
 * - Implements csIReader interface
 * - Encapsulates reader objects of several seismic formats with one interface
 */
class csRandomDataReader : public cseis_geolib::csIReader {
public:

  static const int FORMAT_CSEIS = 11;
  static const int FORMAT_SEGY  = 12;
  static const int FORMAT_RSF   = 13;
public:
  csRandomDataReader( std::string filename, bool enableRandomAccess, int numTracesToBuffer );
  ~csRandomDataReader();

  bool readFileHeader( cseis_system::csSuperHeader* shdr, cseis_system::csTraceHeaderDef* hdef, std::FILE* stream = NULL );

  virtual int numTraces() const;
  virtual int getCurrentTraceIndex() const;
  virtual bool moveToTrace( int traceIndex );
  virtual bool setHeaderToPeek( std::string const& headerName );
  virtual bool setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType );
  virtual bool peekHeaderValue( cseis_geolib::csFlexHeader* value, int traceIndex );

  bool readTrace( float* samples, int numSamples, cseis_system::csTraceHeader* trcHdr = NULL );
  void setTraceHeaderValues( cseis_system::csTraceHeader* trcHdr );

  int fileFormat() const;
  //  void closeFile();

private:
  void setSEGYTraceHeaders( cseis_system::csTraceHeaderDef* hdef );
  void setRSFTraceHeaders( cseis_system::csTraceHeaderDef* hdef );
  void readSegyTraceHeaders( cseis_system::csTraceHeader* trcHdr );
  void readRSFTraceHeaders( cseis_system::csTraceHeader* trcHdr );

  cseis_system::csSeismicReader* mySeismicReader;
  cseis_geolib::csSegyReader* mySegyReader;
  cseis_io::csRSFReader* myRSFReader;
  cseis_geolib::csIReader* myIReader;
  
  int* myHdrIndex;
  cseis_geolib::type_t* myHdrType;

  char* myHdrValueBlock;
  int myHdrValueBlockSize;

  int myFileFormat;
};

} // end namespace
#endif

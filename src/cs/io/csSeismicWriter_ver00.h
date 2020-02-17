/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEISMIC_WRITER_VER_H
#define CS_SEIMSIC_WRITER_VER_H

#include <cstdio>
#include <string>

namespace cseis_geolib {
  class csHeaderInfo;
}

namespace cseis_io {

class csSeismicIOConfig;

/**
 * Seismic file writer, Cseis format
 * Base class for version dependent seismic writer
 *
 * @author Bjorn Olofsson
 * @date 2006
 */
class csSeismicWriter_ver {
 public:
  static int const VERSION_SEISMIC_WRITER = 00;
  static int const DEFAULT_BUFFERED_BYTES = 10000000;
 public:
  csSeismicWriter_ver( std::string filename );
  csSeismicWriter_ver( std::string filename, int numTracesBuffer, bool overwrite = true );
  ~csSeismicWriter_ver();
  bool writeFileHeader( csSeismicIOConfig const* config );
  bool writeTrace( float* samples, char const* hdrValueBlock );
  void close();
public:
  short myVersionMinor;
  short myVersionMajor;
  char myVersionText[4];
  std::string myFileName;
  FILE* myFile;
  /// Size in bytes of trace header value block = All trace values of one trace
  int myByteSizeHdrValueBlock;
  /// Size in bytes of data samples of one trace
  int myByteSizeSamples;
  
private:
  void open( bool overwrite );
  void initialize();
  bool writeCurrentDataBuffer();
  void resizeDataBuffer( int newSize );
    
  char* myTempBuffer;
  int myTempByteSize;
  int myByteLoc;
  void checkBufferSize( int sizeAdd );
  void appendChar( char value );
  void appendInt( int value );
  void appendFloat( float value );
  void appendString( char const* value, int size );

  char* myDataBuffer;
  int   myDataBufferSize;
  int   myCurrentDataBufferSize;
  int   myNumBufferTraces;
};

} // end namespace
#endif


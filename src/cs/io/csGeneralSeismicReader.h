/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_GENERAL_SEISMIC_READER_H
#define CS_GENERAL_SEISMIC_READER_H

#include <cstdio>
#include <string>
#include "geolib_defines.h"
#include "csIReader.h"

namespace cseis_geolib {
  class csFlexHeader;
  class csFlexNumber;
  class csIOSelection;
}

namespace cseis_io {
  class csSeismicReader_ver;
  class csSeismicIOConfig;

/**
 * General seismic file reader, Cseis format
 *
 * @author Bjorn Olofsson
 * @date 2006
 */
class csGeneralSeismicReader : public cseis_geolib::csIReader {
public:
  csGeneralSeismicReader( std::string filename, bool enableRandomAccess, int numTracesToBuffer );
  ~csGeneralSeismicReader();

  bool readFileHeader();
  bool moveToTrace( int traceIndex, int numTracesToRead );
  bool moveToTrace( int traceIndex );
  float const* readTraceReturnPointer();
  bool readTrace( float* samples );
   
  void closeFile();

  int hdrIntValue( int hdrIndex ) const;
  float hdrFloatValue( int hdrIndex ) const;
  double hdrDoubleValue( int hdrIndex ) const;
  csInt64_t hdrInt64Value( int hdrIndex ) const;
  std::string hdrStringValue( int hdrIndex ) const;

  int numTraceHeaders() const;
  int headerIndex( std::string const& headerName ) const;
  std::string headerName( int hdrIndex ) const;
  std::string headerDesc( int hdrIndex ) const;
  cseis_geolib::type_t headerType( int hdrIndex ) const;
  int headerElements( int hdrIndex ) const;
  csSeismicIOConfig const* getSeismicConfig() const { return myConfig; } 
 
  int numTraces() const;
 
  bool setHeaderToPeek( std::string const& headerName );
  bool setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType );

  bool peekHeaderValue( cseis_geolib::csFlexHeader* hdrValue, int traceIndex = -1 );

  bool setSelection( std::string const& hdrValueSelectionText, std::string const& headerName, int sortOrder, int sortMethod );
  void setSelectionStep1( std::string const& hdrValueSelectionText, std::string const& headerName, int sortOrder, int sortMethod );
  bool setSelectionStep2( int& numTracesToRead );
  bool setSelectionStep3();

  /**
   * @return Number of selected traces
   */
  int getNumSelectedTraces() const;
  /**
   * @param traceIndex  Index of consecutively selected trace
   * @return Value of selected trace
   */
  cseis_geolib::csFlexNumber const* getSelectedValue( int traceIndex ) const;
  int getSelectedIndex( int traceIndex ) const;
  int getCurrentTraceIndex() const;

private:
  void setByteLocation();
  bool performIOSelection();

  csSeismicReader_ver* myReader;
  csSeismicIOConfig* myConfig;

  bool    myIsFileHeaderRead;
  int*    myByteLocation;

  float*  myTraceBuffer;
  
  char*   myHdrValueBlock;          // byte header buffer
  int*    myIntPtr;                 // header buffer as int
  float*  myFloatPtr;               // header buffer as float
  
  /// Header to check. 
  int  myHdrCheckByteOffset;
  cseis_geolib::type_t myHdrCheckType;
  int   myHdrCheckByteSize;
  char* myHdrCheckBuffer;
  cseis_geolib::csIOSelection* myIOSelection;
};

} // end namespace
#endif


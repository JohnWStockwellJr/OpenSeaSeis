/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGD_TRACE_HEADER_H
#define CS_SEGD_TRACE_HEADER_H

#include "geolib_defines.h"
#include "csSegdDefines.h"
#include "csHeaderInfo.h"

namespace cseis_geolib {

template<typename T> class csVector;
class csFlexHeader;

}
//


namespace cseis_segd {

/**
 * Class that holds SEG-D header values
 *
 * @author Bjorn Olofsson
 * @date 2011
 */
class csSegdHdrValues {

public:
  csSegdHdrValues();
  ~csSegdHdrValues();
  int numHeaders() const;
  cseis_geolib::csHeaderInfo const* info( int hdrIndex ) const;

  char const* headerDesc( int hdrIndex ) const;
  char const* headerName( int hdrIndex ) const;
  int headerIndex( std::string name ) const;
  cseis_geolib::type_t headerType( int hdrIndex ) const;
    
  int intValue( int hdrIndex ) const;
  float floatValue( int hdrIndex ) const;
  double doubleValue( int hdrIndex ) const;
  std::string stringValue( int hdrIndex ) const;
  
  void setIntValue( int hdrIndex, int value );
  void setFloatValue( int hdrIndex, float value );
  void setDoubleValue( int hdrIndex, double value );
  void setStringValue( int hdrIndex, std::string const& value );

  /**
   * Dump header values
   */
  void dump( FILE* fout );

  void setFileHeaderValues( commonFileHeaderStruct const& comFileHdr );
  void setRecordHeaderValues( commonRecordHeaderStruct const& comRecHdr );
  void setTraceHeaderValues( commonTraceHeaderStruct const& comTrcHdr );

protected:
  cseis_geolib::csFlexHeader* getHandleHdrValues();

private:
  void init();

  cseis_geolib::csFlexHeader* myHdrValues;
  std::string* myHdrNames;
  std::string* myHdrDesc;
  int myNumHeaders;
  int myTraceCounter;
  //--------------------------------------------------
  // File headers
  //
  /*
    /// SEGD file format code
    int hdrId_formatCode;
    /// SEGD revision number
    int hdrId_revisionNum[2];
    /// Manufacturer code
    int hdrId_manufactCode;
    /// Manufacturer serial number
    int hdrId_manufactSerialNum;
    /// Source line name
    char hdrId_srcLineName[16];
    /// Record length [ms]  (= time of last sample..?)
    int hdrId_recordLength_ms;
    /// Sample interval [us]
    int hdrId_sampleInt_us;
    /// Number of samples
    int hdrId_numSamples;
    /// Sample bit size
    int hdrId_sampleBitSize;
    /// Total number of channels
    int hdrId_totalNumChan;
    /// Number of seismic channels
    int hdrId_numSeismicChan; 
    /// Number of auxiliary channels
    int hdrId_numAuxChan;
  */


  /// SEGD file format code
  int hdrId_formatCode;
  /// SEGD revision number
  int hdrId_revisionNum;
  /// Manufacturer code
  int hdrId_manufactCode;
  /// Manufacturer serial number
  int hdrId_manufactSerialNum;
  /// Source line name
  int hdrId_srcLineName;
  /// Record length [ms]  (= time of last sample..?)
  int hdrId_recordLength_ms;
  /// Sample interval [us]
  int hdrId_sampleInt_us;
  /// Number of samples
  int hdrId_numSamples;
  /// Sample bit size
  int hdrId_sampleBitSize;
  /// Total number of channels
  int hdrId_totalNumChan;
  /// Number of seismic channels
  int hdrId_numSeismicChan; 
  /// Number of auxiliary channels
  int hdrId_numAuxChan;

  //--------------------------------------------------
  // Record headers
  //
  /// Shot record field file number (FFID)
  int hdrId_fileNum;
  /// Shot time
  int hdrId_year;
  int hdrId_julianDay;
  int hdrId_hour;
  int hdrId_minute;
  int hdrId_second;
  int hdrId_microSecond;

  /// Shot point number
  int hdrId_shotNum;
  /// Source easting [m]
  int hdrId_srcEasting;
  /// Source northing [m]
  int hdrId_srcNorthing;
  /// Source elevation [m]
  int hdrId_srcElev;

  //---------------------------------------------------
  // Trace headers
  //
  /// Channel number
  int hdrId_chanNum;
  /// Channel type identifier
  int hdrId_chanTypeID;
  /// Trace edit flag...
  int hdrId_traceEdit;
  /// Receiver line number
  int hdrId_rcvLineNumber;
  /// Receiver point number
  int hdrId_rcvPointNumber;
  /// Receiver point index
  int hdrId_rcvPointIndex;
  /// Receiver easting [m]
  int hdrId_rcvEasting;
  /// Receiver northing [m]
  int hdrId_rcvNorthing;
  /// Receiver elevation [m]
  int hdrId_rcvElevation;
  /// Sensor serial number
  int hdrId_serialNumber;
  /// Sensor ID
  int hdrId_sensor;
  /// Inclinometer 'ICV' values
  int hdrId_incl_i;
  int hdrId_incl_c;
  int hdrId_incl_v;

  /// Sequential trace number
  int hdrId_trcno;
  /// Sequential record/file number
  int hdrId_ffid;

  /// Time of first sample [s]
  int hdrId_time_samp1;
  /// Time of first sample, microsecond fraction [us]
  int hdrId_time_samp1_us;

};

} // end namespace
#endif


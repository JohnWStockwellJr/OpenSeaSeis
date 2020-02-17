/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegdHdrValues.h"
#include "csFlexHeader.h"
#include "csException.h"
#include "csStandardHeaders.h"
#include "csGeolibUtils.h"
#include <cstdlib>
#include <cstring>
#include <string>

/*
#include "csSegdHdrMap.h"
#include "csSegdHeaderInfo.h"
#include "csSegdHeader.h"
#include "csFlexNumber.h"
#include "csByteConversions.h"
*/

using namespace cseis_geolib;
using namespace cseis_segd;

csSegdHdrValues::csSegdHdrValues() {
  myNumHeaders = 0;
  myHdrValues  = NULL;
  myHdrNames   = NULL;
  myHdrDesc    = NULL;

  init();
}

csSegdHdrValues::~csSegdHdrValues() {
  if( myHdrValues != NULL ) {
    delete [] myHdrValues;
    myHdrValues = NULL;
  }
  if( myHdrNames != NULL ) {
    delete [] myHdrNames;
    myHdrNames = NULL;
  }
  if( myHdrDesc != NULL ) {
    delete [] myHdrDesc;
    myHdrDesc = NULL;
  }
}

void csSegdHdrValues::init() {
  int id = 0;

  hdrId_fileNum   = id++;
  hdrId_year      = id++;
  hdrId_julianDay = id++;
  hdrId_hour      = id++;
  hdrId_minute    = id++;
  hdrId_second    = id++;
  hdrId_microSecond = id++;

  hdrId_shotNum     = id++;
  hdrId_srcEasting  = id++;
  hdrId_srcNorthing = id++;
  hdrId_srcElev     = id++;

  //---------------------------------------------------
  // Trace headers
  //
  hdrId_chanNum      = id++;
  hdrId_chanTypeID   = id++;
  hdrId_traceEdit    = id++;
  hdrId_rcvLineNumber = id++;
  hdrId_rcvPointNumber = id++;
  hdrId_rcvPointIndex = id++;
  hdrId_rcvEasting   = id++;
  hdrId_rcvNorthing  = id++;
  hdrId_rcvElevation = id++;
  hdrId_serialNumber = id++;
  hdrId_sensor = id++;
  hdrId_incl_i = id++;
  hdrId_incl_c = id++;
  hdrId_incl_v = id++;

  hdrId_trcno  = id++;

  // File headers
  hdrId_formatCode      = id++;
  hdrId_revisionNum     = id++;
  hdrId_manufactCode    = id++;
  hdrId_manufactSerialNum = id++;
  hdrId_srcLineName     = id++;
  hdrId_recordLength_ms = id++;
  hdrId_sampleInt_us    = id++;
  hdrId_numSamples      = id++;
  hdrId_sampleBitSize   = id++;
  hdrId_totalNumChan    = id++;
  hdrId_numSeismicChan  = id++; 
  hdrId_numAuxChan      = id++;

  // Absolute time
  hdrId_time_samp1    = id++;
  hdrId_time_samp1_us = id++;

  myNumHeaders = id;
  myHdrValues = new cseis_geolib::csFlexHeader[myNumHeaders];
  for( int ihdr = 0; ihdr < myNumHeaders; ihdr++ ) {
    myHdrValues[ihdr] = 0;
  }
  myHdrNames = new std::string[myNumHeaders];
  myHdrDesc  = new std::string[myNumHeaders];

  myHdrNames[hdrId_fileNum]   = cseis_geolib::HDR_FFID.name;
  myHdrNames[hdrId_year]      = cseis_geolib::HDR_TIME_YEAR.name;
  myHdrNames[hdrId_julianDay] = cseis_geolib::HDR_TIME_DAY.name;
  myHdrNames[hdrId_hour]      = cseis_geolib::HDR_TIME_HOUR.name;
  myHdrNames[hdrId_minute]    = cseis_geolib::HDR_TIME_MIN.name;
  myHdrNames[hdrId_second]    = cseis_geolib::HDR_TIME_SEC.name;
  myHdrNames[hdrId_microSecond] = cseis_geolib::HDR_TIME_MSEC.name;

  myHdrNames[hdrId_shotNum]     = cseis_geolib::HDR_SOURCE.name;
  myHdrNames[hdrId_srcEasting]  = cseis_geolib::HDR_SOU_X.name;
  myHdrNames[hdrId_srcNorthing] = cseis_geolib::HDR_SOU_Y.name;
  myHdrNames[hdrId_srcElev]     = cseis_geolib::HDR_SOU_ELEV.name;

  myHdrNames[hdrId_chanNum]     = cseis_geolib::HDR_CHAN.name;
  myHdrNames[hdrId_chanTypeID]  = cseis_geolib::HDR_TRC_TYPE.name;
  myHdrNames[hdrId_traceEdit]   = cseis_geolib::HDR_TRC_EDIT.name;

  myHdrNames[hdrId_rcvLineNumber]  = cseis_geolib::HDR_REC_LINE.name;
  myHdrNames[hdrId_rcvPointNumber] = cseis_geolib::HDR_RCV.name;
  myHdrNames[hdrId_rcvPointIndex]  = cseis_geolib::HDR_REC_INDEX.name;
  myHdrNames[hdrId_rcvElevation]   = cseis_geolib::HDR_REC_ELEV.name;
  myHdrNames[hdrId_rcvEasting]     = cseis_geolib::HDR_REC_X.name;
  myHdrNames[hdrId_rcvNorthing]    = cseis_geolib::HDR_REC_Y.name;
  myHdrNames[hdrId_serialNumber]   = cseis_geolib::HDR_SERIAL.name;
  myHdrNames[hdrId_sensor] = cseis_geolib::HDR_SENSOR.name;
  myHdrNames[hdrId_incl_i] = cseis_geolib::HDR_INCL_I.name;
  myHdrNames[hdrId_incl_c] = cseis_geolib::HDR_INCL_C.name;
  myHdrNames[hdrId_incl_v] = cseis_geolib::HDR_INCL_V.name;

  myHdrNames[hdrId_trcno]   = cseis_geolib::HDR_TRCNO.name;

  myHdrNames[hdrId_formatCode]      = "SEG-D format code";
  myHdrNames[hdrId_revisionNum]     = "SEG-D revision code";
  myHdrNames[hdrId_manufactCode]    = "SEG-D manufacturer code";
  myHdrNames[hdrId_manufactSerialNum] = "SEG-D manufacturer serial number";
  myHdrNames[hdrId_srcLineName]     = "Source line name";
  myHdrNames[hdrId_recordLength_ms] = "SEG-D Record length [ms]";
  myHdrNames[hdrId_sampleInt_us]    = cseis_geolib::HDR_SAMPINT_US.name;
  myHdrNames[hdrId_numSamples]      = cseis_geolib::HDR_NSAMP.name;
  myHdrNames[hdrId_sampleBitSize]   = "SEG-D sample bit size";
  myHdrNames[hdrId_totalNumChan]    = "SEG-D total number of channels";
  myHdrNames[hdrId_numSeismicChan]  = "SEG-D number of seismic channels"; 
  myHdrNames[hdrId_numAuxChan]      = "SEG-D number of auxiliary channels";

  myHdrNames[hdrId_time_samp1]      = "time_samp1"; //"Time of first sample [s]";
  myHdrNames[hdrId_time_samp1_us]   = "time_samp1_us"; //Time of first sample, us fraction [us]";

  /*
  trcHdr->setIntValue(vars->hdrId_time_nano, vars->segdReader->extendedHdr()->nanoSeconds());
  trcHdr->setIntValue(vars->hdrId_sline, 0 );
  trcHdr->setIntValue(vars->hdrId_seq, 0 );
  trcHdr->setFloatValue(vars->hdrId_sou_z, 0 );
  trcHdr->setFloatValue(vars->hdrId_rec_wdep, 0 );
  trcHdr->setIntValue( vars->hdrId_fileno, vars->currentFile+1 );
  */

  // Initialise header types by setting to default values
  commonFileHeaderStruct comFileHdr;
  char revNumStr[4];
  sprintf(revNumStr, "%1d.%-1d", comFileHdr.revisionNum[0], comFileHdr.revisionNum[1]);
  revNumStr[3] = '\0';

  commonRecordHeaderStruct comRecHdr;
  commonTraceHeaderStruct comTrcHdr;
  setFileHeaderValues( comFileHdr );
  setRecordHeaderValues( comRecHdr );
  setTraceHeaderValues( comTrcHdr );
  myTraceCounter = 0;
}

void csSegdHdrValues::setFileHeaderValues( commonFileHeaderStruct const& comFileHdr ) {
  myHdrValues[hdrId_formatCode]      .setIntValue(comFileHdr.formatCode);
  char revNumStr[4];
  sprintf(revNumStr, "%1d.%-1d", comFileHdr.revisionNum[0], comFileHdr.revisionNum[1]);
  revNumStr[3] = '\0';
  std::string revNumStr2(revNumStr);
  myHdrValues[hdrId_revisionNum]     .setStringValue(revNumStr2);
  myHdrValues[hdrId_manufactCode]    .setIntValue(comFileHdr.manufactCode);
  myHdrValues[hdrId_manufactSerialNum].setIntValue(comFileHdr.manufactSerialNum);
  std::string srcLineNameStr(comFileHdr.srcLineName);
  myHdrValues[hdrId_srcLineName]     .setStringValue(srcLineNameStr);
  myHdrValues[hdrId_recordLength_ms] .setIntValue(comFileHdr.recordLength_ms);
  myHdrValues[hdrId_sampleInt_us]    .setIntValue(comFileHdr.sampleInt_us);
  myHdrValues[hdrId_numSamples]      .setIntValue(comFileHdr.numSamples);
  myHdrValues[hdrId_sampleBitSize]   .setIntValue(comFileHdr.sampleBitSize);
  myHdrValues[hdrId_totalNumChan]    .setIntValue(comFileHdr.totalNumChan);
  myHdrValues[hdrId_numSeismicChan]  .setIntValue(comFileHdr.numSeismicChan);
  myHdrValues[hdrId_numAuxChan]      .setIntValue(comFileHdr.numAuxChan);
}
void csSegdHdrValues::setRecordHeaderValues( commonRecordHeaderStruct const& comRecHdr ) {
  int year = comRecHdr.shotTime.year;
  if( year < 1900 ) {
    if( year > 30 ) {
      year += 1900;
    }
    else {
      year += 2000;
    }
  }

  myHdrValues[hdrId_fileNum]    .setIntValue(comRecHdr.fileNum);
  myHdrValues[hdrId_year]       .setIntValue(year);
  myHdrValues[hdrId_julianDay]  .setIntValue(comRecHdr.shotTime.julianDay);
  myHdrValues[hdrId_hour]       .setIntValue(comRecHdr.shotTime.hour);
  myHdrValues[hdrId_minute]     .setIntValue(comRecHdr.shotTime.minute);
  myHdrValues[hdrId_second]     .setIntValue(comRecHdr.shotTime.second);
  myHdrValues[hdrId_microSecond].setIntValue(comRecHdr.shotTime.microSecond);

  myHdrValues[hdrId_shotNum]    .setIntValue(comRecHdr.shotNum);
  myHdrValues[hdrId_srcEasting] .setDoubleValue(comRecHdr.srcEasting);
  myHdrValues[hdrId_srcNorthing].setDoubleValue(comRecHdr.srcNorthing);
  myHdrValues[hdrId_srcElev]    .setFloatValue(comRecHdr.srcElev);

  csInt64_t startTime = csGeolibUtils::date2UNIXmsec( year,
                                                      comRecHdr.shotTime.julianDay,
                                                      comRecHdr.shotTime.hour,
                                                      comRecHdr.shotTime.minute,
                                                      comRecHdr.shotTime.second,
                                                      comRecHdr.shotTime.microSecond ); //extendedHdr()->nanoSeconds()/1000 );
  myHdrValues[hdrId_time_samp1]      .setIntValue( (int)(startTime/1000LL) );
  myHdrValues[hdrId_time_samp1_us]   .setIntValue( (int)(startTime%1000LL)*1000 );
}
void csSegdHdrValues::setTraceHeaderValues( commonTraceHeaderStruct const& comTrcHdr ) {
  myHdrValues[hdrId_chanNum]       .setIntValue(comTrcHdr.chanNum);
  myHdrValues[hdrId_chanTypeID]    .setIntValue(comTrcHdr.chanTypeID);
  myHdrValues[hdrId_traceEdit]     .setIntValue(comTrcHdr.traceEdit);
  myHdrValues[hdrId_rcvLineNumber] .setIntValue(comTrcHdr.rcvLineNumber);
  myHdrValues[hdrId_rcvPointNumber].setIntValue(comTrcHdr.rcvPointNumber);
  myHdrValues[hdrId_rcvPointIndex] .setIntValue(comTrcHdr.rcvPointIndex);
  myHdrValues[hdrId_rcvEasting]    .setDoubleValue(comTrcHdr.rcvEasting);
  myHdrValues[hdrId_rcvNorthing]   .setDoubleValue(comTrcHdr.rcvNorthing);
  myHdrValues[hdrId_rcvElevation]  .setFloatValue(comTrcHdr.rcvElevation);
  myHdrValues[hdrId_serialNumber]  .setIntValue(comTrcHdr.serialNumber);
  myHdrValues[hdrId_sensor].setIntValue(comTrcHdr.sensor);
  myHdrValues[hdrId_incl_i].setIntValue(comTrcHdr.incl_i);
  myHdrValues[hdrId_incl_c].setIntValue(comTrcHdr.incl_c);
  myHdrValues[hdrId_incl_v].setIntValue(comTrcHdr.incl_v);

  myHdrValues[hdrId_trcno].setIntValue(++myTraceCounter);
}

void csSegdHdrValues::dump( FILE* fout ) {
  fprintf(stdout,"csSegdHdrValues::dump()\n");
  for( int ihdr = 0; ihdr < myNumHeaders; ihdr++ ) {
    fprintf(stdout,"#%3d:  %s\n", ihdr+1, myHdrValues[ihdr].toString().c_str());
  }
}


int csSegdHdrValues::numHeaders() const {
  return myNumHeaders;
}

int csSegdHdrValues::intValue( int hdrIndex ) const {
  if( hdrIndex < numHeaders() ) {
    return myHdrValues[hdrIndex].intValue();
  }
  else {
    throw( csException("hdrIntValue: Wrong header index passed") );
  }
}
float csSegdHdrValues::floatValue( int hdrIndex ) const {
  if( hdrIndex < numHeaders() ) {
    return myHdrValues[hdrIndex].floatValue();
  }
  else {
    throw( csException("hdrFloatValue: Wrong header index passed") );
  }
}
double csSegdHdrValues::doubleValue( int hdrIndex ) const {
  if( hdrIndex < numHeaders() ) {
    return myHdrValues[hdrIndex].doubleValue();
  }
  else {
    throw( csException("hdrDoubleValue: Wrong header index passed") );
  }
}

std::string csSegdHdrValues::stringValue( int hdrIndex ) const {
  return myHdrValues[hdrIndex].stringValue();
}

cseis_geolib::csFlexHeader* csSegdHdrValues::getHandleHdrValues() {
  return myHdrValues;
}


//-----------------------------------------------------------------------------------------
char const* csSegdHdrValues::headerDesc( int hdrIndex ) const {
  if( hdrIndex < numHeaders() ) {
    return myHdrDesc[hdrIndex].c_str();
  }
  else {
    throw( csException("hdrDesc: Wrong header index passed") );
  }
}
//-----------------------------------------------------------------------------------------
int csSegdHdrValues::headerIndex( std::string name ) const {
  for( int ihdr = 0; ihdr < myNumHeaders; ihdr++ ) {
    if( !name.compare(myHdrNames[ihdr]) ) {
      return ihdr;
    }
  }
  return -1;
}
//-----------------------------------------------------------------------------------------
char const* csSegdHdrValues::headerName( int hdrIndex ) const {
  if( hdrIndex < numHeaders() ) {
    return myHdrNames[hdrIndex].c_str();
  }
  else {
    throw( csException("hdrName: Wrong header index passed") );
  }
}
//-----------------------------------------------------------------------------------------
type_t csSegdHdrValues::headerType( int hdrIndex ) const {
  if( hdrIndex < numHeaders() ) {
    return myHdrValues[hdrIndex].type();
  }
  else {
    throw( csException("hdrType: Wrong header index passed") );
  }
}

void csSegdHdrValues::setIntValue( int hdrIndex, int value ) {
  if( hdrIndex < numHeaders() ) {
    myHdrValues[hdrIndex].setIntValue( value );
  }
  else { 
    throw( csException("setIntValue: Wrong header index passed") );
  }
}
void csSegdHdrValues::setFloatValue( int hdrIndex, float value ) {
  if( hdrIndex < numHeaders() ) {
    myHdrValues[hdrIndex].setFloatValue( value );
  }
  else {
    throw( csException("setFloatValue: Wrong header index passed") );
  }
}
void csSegdHdrValues::setDoubleValue( int hdrIndex, double value ) {
  if( hdrIndex < numHeaders() ) {
    myHdrValues[hdrIndex].setDoubleValue( value );
  }
  else {
    throw( csException("setDoubleValue: Wrong header index passed") );
  }
}
void csSegdHdrValues::setStringValue( int hdrIndex, std::string const& value ) {
  throw( csException("setStringValue: Not implemented yet") );
}


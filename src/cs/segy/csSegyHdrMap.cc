/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegyHdrMap.h"
#include "csSegyDefines.h"
#include "csSegyHeaderInfo.h"
#include "csVector.h"
#include "csGeolibUtils.h"
#include "csFlexNumber.h"
#include "csFlexHeader.h"
#include "csStandardHeaders.h"
#include "geolib_string_utils.h"
#include <fstream>
#include <cmath>

using namespace cseis_geolib;

csSegyHdrMap::csSegyHdrMap( int defaultMap, bool initScalars ) {
  myDefaultMap = defaultMap;
  init( initScalars );
}
csSegyHdrMap::csSegyHdrMap( int defaultMap, bool initScalars, std::string filename_hdrmap ) {
  myDefaultMap = defaultMap;
  init( initScalars );
  readHdrMapExternalFile( filename_hdrmap );
}
csSegyHdrMap::csSegyHdrMap( csSegyHdrMap const* hdrMap ) {
  myDefaultMap    = hdrMap->myDefaultMap;
  myHdrID_scalar_elev  = hdrMap->myHdrID_scalar_elev;
  myHdrID_scalar_coord = hdrMap->myHdrID_scalar_coord;
  myHdrID_scalar_stat  = hdrMap->myHdrID_scalar_stat;

  myHdrList        = new csVector<csSegyHeaderInfo const*>();
  myCoordHdrIDList = new csVector<int>();
  myElevHdrIDList  = new csVector<int>();
  myStatHdrIDList  = new csVector<int>();
  myIsInitScalars = false;

  for( int i = 0; i < hdrMap->myHdrList->size(); i++ ) {
    myHdrList->insertEnd( new csSegyHeaderInfo( *(hdrMap->header(i)) ) );
  }
  /*
  for( int i = 0; i < hdrMap->myCoordHdrIDList->size(); i++ ) {
    myCoordHdrIDList->insertEnd( hdrMap->myCoordHdrIDList->at(i) );
  }
  for( int i = 0; i < hdrMap->myElevHdrIDList->size(); i++ ) {
    myElevHdrIDList->insertEnd( hdrMap->myElevHdrIDList->at(i) );
  }
  */
}

csSegyHdrMap::csSegyHdrMap() {
  myDefaultMap = SEGY_STANDARD;
  init( true );
}

csSegyHdrMap::~csSegyHdrMap() {
  if( myHdrList != NULL ) {
    for( int i = 0; i < myHdrList->size(); i++ ) {
      delete myHdrList->at( i );
    }
    delete myHdrList;
  }
  if( myCoordHdrIDList != NULL ) {
    delete myCoordHdrIDList;
    myCoordHdrIDList= NULL;
  }
  if( myStatHdrIDList != NULL ) {
    delete myStatHdrIDList;
    myStatHdrIDList= NULL;
  }
  if( myElevHdrIDList != NULL ) {
    delete myElevHdrIDList;
    myElevHdrIDList= NULL;
  }
}

csSegyHeaderInfo const* csSegyHdrMap::header( int hdrIndex ) const {
  if( hdrIndex < myHdrList->size() ) {
    return myHdrList->at( hdrIndex );
  }
  return NULL;  // throw csException("Error");
}
csSegyHeaderInfo const* csSegyHdrMap::header( std::string hdrName ) const {
  int hdrIndex = -1;
  if( getHdrIndex( hdrName, &hdrIndex ) ) {
    return myHdrList->at( hdrIndex );
  }
  else {
    return NULL;  // throw csException("Error");
  }
}

int csSegyHdrMap::numHeaders() const {
  return myHdrList->size();
}

//--------------------------------------------------------------------------------
//
//
void csSegyHdrMap::initScalars() {
  myIsInitScalars = true;

  myCoordHdrIDList->clear();
  myElevHdrIDList->clear();
  myStatHdrIDList->clear();

  // Prepare headers that need to be scaled automatically:
  if( !getHdrIndex( HDR_SCALAR_ELEV.name, &myHdrID_scalar_elev ) ) {
    myHdrID_scalar_elev = -1;
    //    throw( csException("Elevation scalar not defined in SEGY trace header mapping. Program bug..") );
  }
  if( !getHdrIndex( HDR_SCALAR_COORD.name, &myHdrID_scalar_coord ) ) {
    myHdrID_scalar_coord = -1;
    //    throw( csException("Coordinate scalar not defined in SEGY trace header mapping. Program bug..") );
  }
  if( !getHdrIndex( HDR_SCALAR_STAT.name, &myHdrID_scalar_stat ) ) {
    myHdrID_scalar_stat = -1;
    //    throw( csException("Statics scalar not defined in SEGY trace header mapping. Program bug..") );
  }
  if( myDefaultMap == SEGY_SU_ONLY ) {
    getHdrIndex( "scalel", &myHdrID_scalar_elev );
    getHdrIndex( "scalco", &myHdrID_scalar_coord );
  }

  int hdrIndex;
  if( getHdrIndex( HDR_REC_X.name, &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_REC_Y.name, &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_SOU_X.name, &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_SOU_Y.name, &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_BIN_X.name, &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_BIN_Y.name, &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );

  if( getHdrIndex( HDR_REC_ELEV.name, &hdrIndex ) )  myElevHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_REC_DATUM.name, &hdrIndex ) ) myElevHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_REC_WDEP.name, &hdrIndex ) )  myElevHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_SOU_Z.name, &hdrIndex ) )     myElevHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_SOU_ELEV.name, &hdrIndex ) )  myElevHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_SOU_DATUM.name, &hdrIndex ) ) myElevHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_SOU_WDEP.name, &hdrIndex ) )  myElevHdrIDList->insertEnd( hdrIndex );  

  if( getHdrIndex( HDR_STAT_REC.name, &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_STAT_SOU.name, &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_STAT_TOT.name, &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_DELAY_TIME.name, &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_MUTE_START.name, &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( HDR_MUTE_END.name, &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( "lag_time_a", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( "lag_time_b", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( "uphole_time_sou", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
  if( getHdrIndex( "uphole_time_rec", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );

  if( myDefaultMap == csSegyHdrMap::SEGY_SU_BOTH || myDefaultMap == csSegyHdrMap::SEGY_SU_ONLY ) {
    if( getHdrIndex( "gx", &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "gy", &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "sx", &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "sy", &hdrIndex ) ) myCoordHdrIDList->insertEnd( hdrIndex );

    if( getHdrIndex( "gelev", &hdrIndex ) )  myElevHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "gdel", &hdrIndex ) ) myElevHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "gwdep", &hdrIndex ) )  myElevHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "sdepth", &hdrIndex ) )     myElevHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "selev", &hdrIndex ) )  myElevHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "sdel", &hdrIndex ) ) myElevHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "swdep", &hdrIndex ) )  myElevHdrIDList->insertEnd( hdrIndex );

    /*
    if( getHdrIndex( "sut", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "gut", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "gstat", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "sstat", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "tstat", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "delrt", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "laga", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "lagb", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "muts", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    if( getHdrIndex( "mute", &hdrIndex ) ) myStatHdrIDList->insertEnd( hdrIndex );
    */
  }
}

//--------------------------------------------------------------------------------
//
//
bool csSegyHdrMap::getHdrIndex( std::string const& name, int* hdrIndex ) const {
  int hdrIndexByteLoc = 0;
  return getHdrIndex( name, hdrIndex, 0, &hdrIndexByteLoc );
}
bool csSegyHdrMap::getHdrIndex( std::string const& name, int* hdrIndex, int byteLocIn, int* hdrIndexByteLoc ) const {
  *hdrIndexByteLoc = 0;
  for( int ihdr = 0; ihdr < myHdrList->size(); ihdr++ ) {
    if( !myHdrList->at(ihdr)->name.compare( name ) ) {
      *hdrIndex = ihdr;
      return true;
    }
    if( myHdrList->at(ihdr)->byteLoc < byteLocIn ) *hdrIndexByteLoc = ihdr+1;
  }
  *hdrIndex = -1;
  return false;
}

//-----------------------------------------------------------------------------------------
bool csSegyHdrMap::addHeader( int byteLoc, int byteSize, type_t inType, csHeaderInfo const& hdr ) {
  return addHeader( hdr.name, byteLoc, byteSize, inType, hdr.type, hdr.description );
}

//-----------------------------------------------------------------------------------------
bool csSegyHdrMap::addHeader( csSegyHeaderInfo const& hdr ) {
  return addHeader( hdr.name, hdr.byteLoc, hdr.byteSize, hdr.inType, hdr.outType, hdr.description );
}
//-----------------------------------------------------------------------------------------
bool csSegyHdrMap::addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theInType ) {
  return addHeader( theName, theByteLoc, theByteSize, theInType, theInType, "" );
}
//-----------------------------------------------------------------------------------------
bool csSegyHdrMap::addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theInType, std::string const& theDesc ) {
  return addHeader( theName, theByteLoc, theByteSize, theInType, theInType, theDesc );
}
//-----------------------------------------------------------------------------------------
bool csSegyHdrMap::addHeader( std::string const& theName, int theByteLoc, int theByteSize, type_t theInType, type_t theOutType, std::string const& theDesc ) {
  if( myIsInitScalars ) {
    throw( csException("csSegyHdrMap::addHeader: Accessing method to set new header AFTER initializing headers. This is a program bug in the calling method") );
  }
  if( theByteLoc > cseis_segy::SIZE_TRCHDR ) {
    throw csException("Specified byte location is larger than actual SEGY header block");
  }
  else if( theByteLoc < 0 ) {
    throw csException("Specified byte location is smaller than zero");
  }
  else if( theByteLoc+theByteSize > cseis_segy::SIZE_TRCHDR ) {
    throw csException("Specified byte location (+ byte size) is larger than SEGY trace header block");
  }

  int hdrIndex;
  int hdrIndexByteLoc;
  if( !getHdrIndex( theName, &hdrIndex, theByteLoc, &hdrIndexByteLoc ) ) {
    myHdrList->insertEnd( new csSegyHeaderInfo( theByteLoc, theByteSize, theInType, theOutType, theName, theDesc ) );
    // myHdrList->insert( new csSegyHeaderInfo( theByteLoc, theByteSize, theInType, theOutType, theName, theDesc ), hdrIndexByteLoc );
    return true;
  }
  else if( myReplaceExistingHeader ) {
    removeHeader( hdrIndex );
    if( !getHdrIndex( theName, &hdrIndex, theByteLoc, &hdrIndexByteLoc ) ) {
      myHdrList->insertEnd( new csSegyHeaderInfo( theByteLoc, theByteSize, theInType, theOutType, theName, theDesc ) );
      //  myHdrList->insert( new csSegyHeaderInfo( theByteLoc, theByteSize, theInType, theOutType, theName, theDesc ), hdrIndexByteLoc );
      return true;
    }
    else {
      return false;
    }
  }
  else {
    return false;
  }
}
bool csSegyHdrMap::contains( std::string hdrName, int* hdrIndex ) const {
  return getHdrIndex( hdrName, hdrIndex );
}
int csSegyHdrMap::headerIndex( std::string hdrName ) const {
  int hdrIndex = -1;
  getHdrIndex( hdrName, &hdrIndex );
  return hdrIndex;
}
bool csSegyHdrMap::removeHeader( int hdrIndex ) {
  if( hdrIndex < myHdrList->size() ) {
    myHdrList->remove( hdrIndex );
    return true;
  }
  else {
    return false;
  }
}
bool csSegyHdrMap::removeHeader( std::string theName ) {
  int hdrIndex;
  if( getHdrIndex( theName, &hdrIndex ) ) {
    myHdrList->remove( hdrIndex );
    return true;
  }
  return false;
}
void csSegyHdrMap::removeAll() {
  for( int i = 0; i < myHdrList->size(); i++ ) {
    delete myHdrList->at( i );
  }
  myHdrList->clear();
}

// Special methods for application of coordinate scalars
int csSegyHdrMap::numCoordHeaders() const {
  return myCoordHdrIDList->size();
}
int csSegyHdrMap::numElevHeaders() const {
  return myElevHdrIDList->size();
}
int csSegyHdrMap::numStatHeaders() const {
  return myStatHdrIDList->size();
}
int csSegyHdrMap::getCoordScalarHeaderIndex() const {
  return myHdrID_scalar_coord;
}
int csSegyHdrMap::getElevScalarHeaderIndex() const {
  return myHdrID_scalar_elev;
}
int csSegyHdrMap::getStatScalarHeaderIndex() const {
  return myHdrID_scalar_stat;
}
int csSegyHdrMap::getCoordHeaderIndex( int index ) const {
  return myCoordHdrIDList->at(index);
}
int csSegyHdrMap::getElevHeaderIndex( int index ) const {
  return myElevHdrIDList->at(index);
}
int csSegyHdrMap::getStatHeaderIndex( int index ) const {
  return myStatHdrIDList->at(index);
}

//--------------------------------------------------------------------------------
//
//
void csSegyHdrMap::applyCoordinateScalar( cseis_geolib::csFlexHeader* hdrValues ) const {
  if( myHdrID_scalar_elev < 0 || myHdrID_scalar_coord < 0 ) throw( csException("csSegyHdrMap::applyCoordinateScalar: Inconsistent call to method. This method is not implemented in the best possible way. Refactoring required. Program bug.") );
  // 1) Convert scalar values into signed shorts
  //  hdrValues[myHdrID_scalar_elev].setIntValue( (int)(short)hdrValues[myHdrID_scalar_elev].intValue() );
  //  hdrValues[myHdrID_scalar_coord].setIntValue( (int)(short)hdrValues[myHdrID_scalar_coord].intValue() );
  //  hdrValues[myHdrID_scalar_stat].setIntValue( (int)(short)hdrValues[myHdrID_scalar_stat].intValue() );

  double scalarElev  = (double)hdrValues[myHdrID_scalar_elev].intValue();
  double scalarCoord = (double)hdrValues[myHdrID_scalar_coord].intValue();
  double scalarStat  = (double)hdrValues[myHdrID_scalar_stat].intValue();

//  fprintf(stdout,"csSegyHdrMap:read: Scalars: %f %f\n", scalarElev, scalarCoord );
  
  if( scalarElev  < 0 ) scalarElev  = -1/scalarElev;
  if( scalarCoord < 0 ) scalarCoord = -1/scalarCoord;
  if( scalarStat  < 0 ) scalarStat = -1/scalarStat;

  if( scalarElev != 0.0 ) {
    for( int i = 0; i < myElevHdrIDList->size(); i++ ) {
      int headerIndex = myElevHdrIDList->at(i);
      hdrValues[headerIndex].setFloatValue( (float)( (double)(hdrValues[headerIndex].intValue()) * scalarElev ) );
      //  fprintf(stdout,"Scale elev %1d %d %f\n", i, hdrValues[headerIndex].intValue(), hdrValues[headerIndex].floatValue() );
    }
  }
  if( scalarCoord != 0.0 ) {
    for( int i = 0; i < myCoordHdrIDList->size(); i++ ) {
      int headerIndex = myCoordHdrIDList->at(i);
      hdrValues[headerIndex].setDoubleValue( (double)(hdrValues[headerIndex].intValue()) * scalarCoord );
      //      fprintf(stdout,"Scale header %1d %d %f   %s  %d %d\n", i, hdrValues[headerIndex].intValue(), hdrValues[headerIndex].doubleValue(), myHdrList->at(headerIndex)->name.c_str(),
      //         myHdrList->at(headerIndex)->byteLoc, myHdrList->at(headerIndex)->byteSize );
    }
  }
  if( scalarStat != 0.0 ) {
    for( int i = 0; i < myStatHdrIDList->size(); i++ ) {
      int headerIndex = myStatHdrIDList->at(i);
      hdrValues[headerIndex].setDoubleValue( (double)(hdrValues[headerIndex].intValue()) * scalarStat );
      //  fprintf(stdout,"Scale header %1d %d %f\n", i, hdrValues[headerIndex].intValue(), hdrValues[headerIndex].doubleValue() );
    }
  }
}

void csSegyHdrMap::applyCoordinateScalarWriting( cseis_geolib::csFlexHeader* hdrValues ) const {
  if( myHdrID_scalar_elev < 0 || myHdrID_scalar_coord < 0 ) throw( csException("csSegyHdrMap::applyCoordinateScalar: Inconsistent call to method. This method is not implemented in the best possible way. Refactoring required. Program bug.") );
  double scalarElev  = (double)hdrValues[myHdrID_scalar_elev].intValue();
  double scalarCoord = (double)hdrValues[myHdrID_scalar_coord].intValue();
  double scalarStat  = (double)hdrValues[myHdrID_scalar_stat].intValue();
  // Swap scalar polarity on output
  hdrValues[myHdrID_scalar_elev].setIntValue( -scalarElev );
  hdrValues[myHdrID_scalar_coord].setIntValue( -scalarCoord );
  hdrValues[myHdrID_scalar_stat].setIntValue( -scalarStat );
  if( scalarElev  < 0 ) scalarElev  = -1/scalarElev;
  if( scalarCoord < 0 ) scalarCoord = -1/scalarCoord;
  if( scalarStat < 0 ) scalarStat = -1/scalarStat;

  // fprintf(stdout,"csSegyHdrMap:write: Scalars: %f %f\n", scalarElev, scalarCoord );

  if( scalarElev != 0.0 ) {
    for( int i = 0; i < myElevHdrIDList->size(); i++ ) {
      int headerIndex = myElevHdrIDList->at(i);
      hdrValues[headerIndex].setIntValue( (int)round(hdrValues[headerIndex].floatValue()*scalarElev) );
      //fprintf(stdout,"Scale elev %1d %d %f\n", i, hdrValues[headerIndex].intValue(), hdrValues[headerIndex].floatValue() );
    }
  }
  if( scalarCoord != 0.0 ) {
    for( int i = 0; i < myCoordHdrIDList->size(); i++ ) {
      int headerIndex = myCoordHdrIDList->at(i);
      hdrValues[headerIndex].setIntValue( (int)round(hdrValues[headerIndex].doubleValue()*scalarCoord) );
      //fprintf(stdout,"Scale header %1d %d %f\n", i, hdrValues[headerIndex].intValue(), hdrValues[headerIndex].doubleValue() );
    }
  }
  if( scalarStat != 0.0 ) {
    for( int i = 0; i < myStatHdrIDList->size(); i++ ) {
      int headerIndex = myStatHdrIDList->at(i);
      hdrValues[headerIndex].setIntValue( (int)round(hdrValues[headerIndex].doubleValue()*scalarStat) );
      //fprintf(stdout,"Scale header %1d %d %f\n", i, hdrValues[headerIndex].intValue(), hdrValues[headerIndex].doubleValue() );
    }
  }
}

//-------------------------------------------------------------------
void csSegyHdrMap::init( bool initScalars_in ) {
  myIsInitScalars = false;
  myReplaceExistingHeader = false;
  myHdrList = NULL;
  myCoordHdrIDList = NULL;
  myElevHdrIDList  = NULL;
  myStatHdrIDList  = NULL;

  myHdrList = new csVector<csSegyHeaderInfo const*>();
  myCoordHdrIDList = new csVector<int>();
  myElevHdrIDList  = new csVector<int>();
  myStatHdrIDList  = new csVector<int>();

  if( myDefaultMap != SEGY_PSEGY && myDefaultMap != csSegyHdrMap::SEGY_SU_ONLY && myDefaultMap != csSegyHdrMap::NONE ) {
    addHeader( "trcseq_line", 0, 4, TYPE_INT, "Trace sequence number within line" );
    addHeader( 4, 4, TYPE_INT, HDR_TRCNO );
    addHeader( 8, 4, TYPE_INT, HDR_FFID );
    addHeader( 12, 4, TYPE_INT, HDR_CHAN );
    addHeader( 16, 4, TYPE_INT, HDR_SOURCE );
    addHeader( 20, 4, TYPE_INT, HDR_CMP );
    addHeader( 24, 4, TYPE_INT, HDR_CHAN_ENS );
    addHeader( 28, 2, TYPE_INT, HDR_TRC_TYPE );
    addHeader( 30, 2, TYPE_INT, HDR_FOLD_VERT );
    addHeader( 32, 2, TYPE_INT, HDR_FOLD );
    addHeader( 34, 2, TYPE_INT, HDR_DATA_TYPE );
    addHeader( 36, 4, TYPE_INT, HDR_OFFSET );
    addHeader( 40, 4, TYPE_INT, HDR_REC_ELEV );
    addHeader( 44, 4, TYPE_INT, HDR_SOU_ELEV );
    addHeader( 48, 4, TYPE_INT, HDR_SOU_Z );
    addHeader( 52, 4, TYPE_INT, HDR_REC_DATUM );
    addHeader( 56, 4, TYPE_INT, HDR_SOU_DATUM );
    addHeader( 60, 4, TYPE_INT, HDR_SOU_WDEP );
    addHeader( 64, 4, TYPE_INT, HDR_REC_WDEP );
    addHeader( 68, 2, TYPE_INT, HDR_SCALAR_ELEV );
    addHeader( 70, 2, TYPE_INT, HDR_SCALAR_COORD );
    addHeader( 72, 4, TYPE_INT, HDR_SOU_X );
    addHeader( 76, 4, TYPE_INT, HDR_SOU_Y );
    addHeader( 80, 4, TYPE_INT, HDR_REC_X);
    addHeader( 84, 4, TYPE_INT, HDR_REC_Y);
    addHeader( 88, 2, TYPE_INT, HDR_UNIT_COORD);

    addHeader(  98, 2, TYPE_INT, HDR_STAT_SOU );
    addHeader( 100, 2, TYPE_INT, HDR_STAT_REC );
    addHeader( 102, 2, TYPE_INT, HDR_STAT_TOT );
    addHeader( 108, 2, TYPE_INT, HDR_DELAY_TIME );

    addHeader( 110, 2, TYPE_INT, HDR_MUTE_START );
    addHeader( 112, 2, TYPE_INT, HDR_MUTE_END);
    addHeader( 114, 2, TYPE_USHORT, HDR_NSAMP );
    addHeader( 116, 2, TYPE_USHORT, HDR_SAMPINT_US );

    addHeader( 156, 2, TYPE_INT, HDR_TIME_YEAR);
    addHeader( 158, 2, TYPE_INT, HDR_TIME_DAY);
    addHeader( 160, 2, TYPE_INT, HDR_TIME_HOUR);
    addHeader( 162, 2, TYPE_INT, HDR_TIME_MIN);
    addHeader( 164, 2, TYPE_INT, HDR_TIME_SEC);
    addHeader( 166, 2, TYPE_INT, HDR_TIME_CODE);

    addHeader( 214, 2, TYPE_INT, HDR_SCALAR_STAT );
  }

  if( myDefaultMap == csSegyHdrMap::NONE ) {
    addHeader( 156, 2, TYPE_INT, HDR_TIME_YEAR);
    addHeader( 158, 2, TYPE_INT, HDR_TIME_DAY);
    addHeader( 160, 2, TYPE_INT, HDR_TIME_HOUR);
    addHeader( 162, 2, TYPE_INT, HDR_TIME_MIN);
    addHeader( 164, 2, TYPE_INT, HDR_TIME_SEC);
    addHeader( 166, 2, TYPE_INT, HDR_TIME_CODE);
  }
  else if( myDefaultMap == SEGY_STANDARD ) {
    addHeader( "vel_weathering", 90, 2, TYPE_INT, TYPE_FLOAT, "Weathering velocity");
    addHeader( "vel_subweathering", 92, 2, TYPE_INT, TYPE_FLOAT, "Subweathering velocity");
    addHeader( "uphole_time_sou", 94, 2, TYPE_INT, TYPE_FLOAT, "Uphole time at source [ms]");
    addHeader( "uphole_time_rec",96, 2, TYPE_INT, TYPE_FLOAT, "Uphole time at group [ms]");

    addHeader( "lag_time_a", 104, 2, TYPE_INT, TYPE_FLOAT, "Lag time A" );
    addHeader( "lag_time_b", 106, 2, TYPE_INT, TYPE_FLOAT, "Lag time B" );

    addHeader( 118, 2, TYPE_INT, HDR_GAIN_TYPE );
    addHeader( 120, 2, TYPE_INT, HDR_GAIN );
    addHeader( "gain_init", 122, 2, TYPE_INT, "Instrument early or initial gain [dB]" );
    addHeader( "correlated", 124, 2, TYPE_INT, "Correlated" );
    addHeader( "sweep_freq_start", 126, 2, TYPE_INT, "Sweep frequency at start [Hz]" );
    addHeader( "sweep_freq_end", 128, 2, TYPE_INT, "Sweep frequency at end [Hz]" );
    addHeader( "sweep_len", 130, 2, TYPE_INT, "Sweep length [ms]" );
    addHeader( "sweep_type", 132, 2, TYPE_INT, "Sweep type" );
    addHeader( "sweep_taper_start", 134, 2, TYPE_INT, "Sweep trace taper length at start [ms]" );
    addHeader( "sweep_taper_end", 136, 2, TYPE_INT, "Sweep trace taper length at end [ms]" );
    addHeader( "taper_type", 138, 2, TYPE_INT, "Taper type" );
    addHeader( "filt_alias_freq", 140, 2, TYPE_INT, "Alias filter frequency [Hz]" );
    addHeader( "filt_alias_db", 142, 2, TYPE_INT, "Alias filter slope [dB/oct]" );
    addHeader( "filt_notch_freq", 144, 2, TYPE_INT, "Notch filter frequency [Hz]" );
    addHeader( "filt_notch_db", 146, 2, TYPE_INT, "Notch filter slope [dB/oct]" );
    addHeader( 148, 2, TYPE_INT, HDR_FILT_LOW_FREQ ); // Corrected 14/11/4
    addHeader( 150, 2, TYPE_INT, HDR_FILT_HIGH_FREQ );// Corrected 14/11/4
    addHeader( 152, 2, TYPE_INT, HDR_FILT_LOW_DB );
    addHeader( 154, 2, TYPE_INT, HDR_FILT_HIGH_DB );

    addHeader( "trc_weight_factor", 168, 2, TYPE_INT, "Trace weighting factor");
    addHeader( "rcv_roll_switch", 170, 2, TYPE_INT, "Geophone group number of roll switch position one");
    addHeader( "rcv_first_trc", 172, 2, TYPE_INT, "Geophone group number of first trace within original field record");
    addHeader( "rcv_last_trc", 174, 2, TYPE_INT, "Geophone group number of last trace within original field record");
    addHeader( "gap_size", 176, 2, TYPE_INT, "Gap size");
    addHeader( "over_travel", 178, 2, TYPE_INT, "Over travel associated with taper at beginning or end of line");

    addHeader( 180, 4, TYPE_INT, HDR_BIN_X);
    addHeader( 184, 4, TYPE_INT, HDR_BIN_Y);
    addHeader( 188, 4, TYPE_INT, HDR_ROW);
    addHeader( 192, 4, TYPE_INT, HDR_COL);

    addHeader( "sp_num", 196, 4, TYPE_INT, "Shotpoint number (SEGY)");
    addHeader( "scalar_sp", 200, 2, TYPE_INT, "Scalar to be applied to shotpoint number");
    addHeader( "trc_val_measure", 202, 2, TYPE_INT, "Trace value measurement unit");
    addHeader( "trans_const", 204, 6, csSegyHdrMap::SEGY_HDR_TYPE_6BYTE, TYPE_DOUBLE, "Transduction constant");
    addHeader( "trans_unit", 210, 2, TYPE_INT, "Transduction units");
    addHeader( "device_ident", 212, 2, TYPE_INT, "Device/trace identifier");

    addHeader( "sou_type", 216, 2, TYPE_INT, "Source type/orientation");
    addHeader( "sou_direction", 218, 6, csSegyHdrMap::SEGY_HDR_TYPE_6BYTE, TYPE_DOUBLE, "Source energy direction with respect to the source orientation");
//    addHeader( "sou_measure", 224, 5, ???, TYPE_DOUBLE, "Source measurement");
    addHeader( "sou_measure_unit", 230, 2, TYPE_INT, "Source measurement unit");
  }
  else if( myDefaultMap == csSegyHdrMap::SEGY_SU || myDefaultMap == csSegyHdrMap::SEGY_SU_ONLY || myDefaultMap == csSegyHdrMap::SEGY_SU_BOTH ) {

    if( myDefaultMap == csSegyHdrMap::SEGY_SU_BOTH || myDefaultMap == csSegyHdrMap::SEGY_SU_ONLY ) {
      // Trace headers that are already defined with different names in Seaseis
      addHeader( "tracl", 0, 4, TYPE_INT, "Trace sequence number within line" );
      addHeader( "tracr", 4, 4, TYPE_INT,"Trace sequence number within file" );
      addHeader( "fldr", 8, 4, TYPE_INT, "Original field record number" );
      addHeader( "tracf", 12, 4, TYPE_INT, "Trace number within the original field record" );
      addHeader( "ep", 16, 4, TYPE_INT, "Energy source point number" );
      addHeader( "cdp", 20, 4, TYPE_INT, "Ensemble number" );
      addHeader( "cdpt",  24, 4, TYPE_INT, "Trace number within the ensemble" );
      addHeader( "trid",  28, 2, TYPE_INT, "Trace identification code" );
      addHeader( "nvs", 30, 2, TYPE_INT, "Number of vertically summed traces yielding this trace" );
      addHeader( "nhs", 32, 2, TYPE_INT, "Number of horizontally summed traces yielding this trace" );
      addHeader( "duse", 34, 2, TYPE_INT, "Data use" );

      addHeader( "gelev", 40, 4, TYPE_INT, TYPE_FLOAT, "Receiver group elevation" );
      addHeader( "selev", 44, 4, TYPE_INT, TYPE_FLOAT, "Surface elevation at source" );
      addHeader( "sdepth", 48, 4, TYPE_INT, TYPE_FLOAT, "Source depth below surface" );
      addHeader( "gdel", 52, 4, TYPE_INT, TYPE_FLOAT, "Datum elevation at receiver group" );
      addHeader( "sdel", 56, 4, TYPE_INT, TYPE_FLOAT, "Datum elevation at source" );
      addHeader( "swdep", 60, 4, TYPE_INT, TYPE_FLOAT, "Water depth at source" );
      addHeader( "gwdep", 64, 4, TYPE_INT, TYPE_FLOAT, "Water depth at group" );
      addHeader( "scalel", 68, 2, TYPE_INT, "Scalar to be applied to all elevations and depths" );
      addHeader( "scalco", 70, 2, TYPE_INT, "Scalar to be applied to all coordinates" );
      addHeader( "sx", 72, 4, TYPE_INT, TYPE_DOUBLE, "Source coordinate - X" );
      addHeader( "sy", 76, 4, TYPE_INT, TYPE_DOUBLE, "Source coordinate - Y" );
      addHeader( "gx", 80, 4, TYPE_INT, TYPE_DOUBLE, "Group coordinate - X" );
      addHeader( "gy", 84, 4, TYPE_INT, TYPE_DOUBLE, "Group coordinate - Y" );
      addHeader( "counit", 88, 2, TYPE_INT, "Coordinate units" );

      addHeader( "sstat",  98, 2, TYPE_INT, TYPE_FLOAT, "Source static correction [ms]" );
      addHeader( "gstat", 100, 2, TYPE_INT, TYPE_FLOAT, "Group static correction [ms]" );
      addHeader( "tstat", 102, 2, TYPE_INT, TYPE_FLOAT, "Total static correction [ms]" );

      addHeader( "delrt", 108, 2, TYPE_INT, TYPE_FLOAT, "Delay recording time" );
      addHeader( "muts", 110, 2, TYPE_INT, TYPE_FLOAT, "Mute start time [ms]" );
      addHeader( "mute", 112, 2, TYPE_INT, TYPE_FLOAT, "Mute end time [ms]" );
      addHeader( "ns", 114, 2, TYPE_USHORT, TYPE_INT, "Number of samples" );
      addHeader( "dt", 116, 2, TYPE_USHORT, TYPE_INT, "Sample interval [us]" );
      addHeader( "ns", 114, 2, TYPE_INT, "Number of samples" );
      addHeader( "dt", 116, 2, TYPE_INT, "Sample interval [us]" );

      addHeader( "year", 156, 2, TYPE_INT, "Year");
      addHeader( "day", 158, 2, TYPE_INT, "Day of year");
      addHeader( "hour", 160, 2, TYPE_INT, "Hour");
      addHeader( "minute", 162, 2, TYPE_INT, "Minute");
      addHeader( "sec", 164, 2, TYPE_INT, "Second");
      addHeader( "timbas", 166, 2, TYPE_INT, "Time bases code");

    }
    addHeader( "wevel", 90, 2, TYPE_INT, TYPE_FLOAT, "Weathering velocity");
    addHeader( "swevel", 92, 2, TYPE_INT, TYPE_FLOAT, "Subweathering velocity");
    addHeader( "sut", 94, 2, TYPE_INT, TYPE_FLOAT, "Uphole time at source [ms]");
    addHeader( "gut",96, 2, TYPE_INT, TYPE_FLOAT, "Uphole time at group [ms]");

    addHeader( "laga", 104, 2, TYPE_INT, TYPE_FLOAT, "Lag time A" );
    addHeader( "lagb", 106, 2, TYPE_INT, TYPE_FLOAT, "Lag time B" );

    addHeader( "gain", 118, 2, TYPE_INT, "Gain type" );
    addHeader( "igc", 120, 2, TYPE_INT, "Instrument gain constant [dB]" );
    addHeader( "igi", 122, 2, TYPE_INT, "Instrument early or initial gain [dB]" );
    addHeader( "corr", 124, 2, TYPE_INT, "Correlated" );
    addHeader( "sfs", 126, 2, TYPE_INT, "Sweep frequency at start [Hz]" );
    addHeader( "sfe", 128, 2, TYPE_INT, "Sweep frequency at end [Hz]" );
    addHeader( "slen", 130, 2, TYPE_INT, "Sweep length [ms]" );
    addHeader( "styp", 132, 2, TYPE_INT, "Sweep type" );
    addHeader( "stas", 134, 2, TYPE_INT, "Sweep trace taper length at start [ms]" );
    addHeader( "stae", 136, 2, TYPE_INT, "Sweep trace taper length at end [ms]" );
    addHeader( "tatyp", 138, 2, TYPE_INT, "Taper type" );
    addHeader( "afilf", 140, 2, TYPE_INT, "Alias filter frequency [Hz]" );
    addHeader( "afils", 142, 2, TYPE_INT, "Alias filter slope [dB/oct]" );
    addHeader( "nofilf", 144, 2, TYPE_INT, "Notch filter frequency [Hz]" );
    addHeader( "nofils", 146, 2, TYPE_INT, "Notch filter slope [dB/oct]" );
    addHeader( "lcf", 148, 2, TYPE_INT, "Low-cut frequency [Hz]" );
    addHeader( "hcf", 150, 2, TYPE_INT, "High-cut frequency [Hz]" );
    addHeader( "lcs", 152, 2, TYPE_INT, "Low-cut slope [dB/oct]" );
    addHeader( "hcs", 154, 2, TYPE_INT, "High-cut slope [dB/oct]" );

    addHeader( "trwf", 168, 2, TYPE_INT, "Trace weighting factor");
    addHeader( "grnors", 170, 2, TYPE_INT, "Geophone group number of roll switch position one");
    addHeader( "grnofr", 172, 2, TYPE_INT, "Geophone group number of first trace within original field record");
    addHeader( "grnlof", 174, 2, TYPE_INT, "Geophone group number of last trace within original field record");
    addHeader( "gaps", 176, 2, TYPE_INT, "Gap size");
    addHeader( "otrav", 178, 2, TYPE_INT, "Over travel associated with taper at beginning or end of line");

    addHeader( "d1", 180, 4, TYPE_FLOAT, "Sample spacing for non-seismic data");
    addHeader( "f1", 184, 4, TYPE_FLOAT, "First sample location for non-seismic data");
    addHeader( "d2", 188, 4, TYPE_FLOAT, "Sample spacing between traces");
    addHeader( "f2", 192, 4, TYPE_FLOAT, "First trace location");

    addHeader( "ungpow", 196, 4, TYPE_INT, "Negative power used for dynamic range compression");
    addHeader( "unscale", 200, 4, TYPE_INT, "Reciprocal of scaling factor to normalize range");
    addHeader( "ntr", 204, 2, TYPE_INT, "Number of traces");
    addHeader( "mark", 208, 2, TYPE_INT, "Mark selected traces");

    // Add standard SEGY trace headers that have the same name in SU as in Seaseis
    if( myDefaultMap == csSegyHdrMap::SEGY_SU_ONLY ) {
      addHeader( 36, 4, TYPE_INT, HDR_OFFSET );
      // addHeader( "", 214, 2, TYPE_INT, SCALAR_STAT );
    }
  }
  else if( myDefaultMap == SEGY_OBC ) {
    addHeader( 180, 4, TYPE_INT, HDR_SOU_LINE);
    addHeader( 184, 4, TYPE_INT, HDR_REC_LINE);
    addHeader( 188, 2, TYPE_INT, HDR_GUN_SEQ);
    addHeader( 190, 2, TYPE_INT, HDR_SENSOR);
    addHeader( 192, 4, TYPE_INT, HDR_RCV);
    addHeader( 196, 4, TYPE_INT, HDR_SERIAL);
    addHeader( 200, 4, TYPE_INT, HDR_SEQ);
    addHeader( 204, 4, TYPE_INT, HDR_INCL_I);
    addHeader( 208, 4, TYPE_INT, HDR_INCL_C);
    addHeader( 212, 4, TYPE_INT, HDR_INCL_V);
    addHeader( 216, 4, TYPE_INT, HDR_ORIENT_I);
    addHeader( 220, 4, TYPE_INT, HDR_ORIENT_C);
    addHeader( 224, 4, TYPE_INT, HDR_ORIENT_V);
    addHeader( 228, 4, TYPE_INT, HDR_CBL_AZIM);
  }
  else if( myDefaultMap == SEGY_SEND ) {
    addHeader( 190, 2, TYPE_INT, HDR_SENSOR);
    addHeader( 192, 4, TYPE_INT, HDR_RCV);
    addHeader( 200, 2, TYPE_INT, HDR_SEQ);
    addHeader( 212, 2, TYPE_INT, HDR_SERIAL);
  }
  else if( myDefaultMap == SEGY_ARMSS ) {
    addHeader( "armss_timing_word", 90, 2, TYPE_INT, "Time word (us between shot and first sample)" );
    addHeader( "armss_time0_s", 92, 4, TYPE_INT, "Time of first sample [s]" );
    addHeader( "armss_time0_us", 96, 4, TYPE_INT, "Time of first sample [us]" );
    addHeader( "armss_shottime_s", 100, 4, TYPE_INT, "Shot time [s]" );
    addHeader( "armss_shottime_us", 104, 4, TYPE_INT, "Shot time [us]" );
    addHeader( "armss_mes_gain1", 120, 2, TYPE_INT, "MES gain 1" );
    addHeader( "armss_mes_gain2", 122, 2, TYPE_INT, "MES gain 2" );

    addHeader( 140, 2, TYPE_INT, HDR_FILT_HIGH_FREQ);
    addHeader( 142, 2, TYPE_INT, HDR_FILT_HIGH_DB);
    addHeader( 148, 2, TYPE_INT, HDR_FILT_LOW_FREQ);
    addHeader( 152, 2, TYPE_INT, HDR_FILT_LOW_DB);

    addHeader( 190, 2, TYPE_INT, HDR_SENSOR);
    addHeader( 192, 4, TYPE_INT, HDR_RCV);
    addHeader( 212, 2, TYPE_INT, HDR_SERIAL);
    addHeader( 200, 2, TYPE_INT, HDR_SEQ);
    addHeader( 216, 2, TYPE_INT, HDR_DC);
    addHeader( 230, 2, TYPE_INT, HDR_AN_ROLL);
    addHeader( 232, 2, TYPE_INT, HDR_AN_TILT);

    addHeader( "armss_ain1", 234, 2, TYPE_INT, "Accelerometer 1 (ARMSS ain1)" );
    addHeader( "armss_ain2", 236, 2, TYPE_INT, "Accelerometer 2 (ARMSS ain2)" );
    addHeader( "armss_ain3", 238, 2, TYPE_INT, "Accelerometer 3 (ARMSS ain3)" );
  }
  else if( myDefaultMap == SEGY_PSEGY ) {

    addHeader( 0, 4, TYPE_INT, HDR_TRCNO );
    addHeader( 8, 4, TYPE_INT, HDR_FFID );
    addHeader( 12, 4, TYPE_INT, HDR_CHAN );

    addHeader( 28, 2, TYPE_INT, HDR_TRC_TYPE );

    addHeader( 68, 2, TYPE_INT, HDR_SCALAR_ELEV );
    addHeader( 70, 2, TYPE_INT, HDR_SCALAR_COORD );

    addHeader( 72, 4, TYPE_INT, HDR_SOU_X );
    addHeader( 76, 4, TYPE_INT, HDR_SOU_Y );
    addHeader( 80, 4, TYPE_INT, HDR_REC_X);
    addHeader( 84, 4, TYPE_INT, HDR_REC_Y);
    addHeader( 88, 2, TYPE_INT, HDR_UNIT_COORD);

    addHeader( 156, 2, TYPE_INT, HDR_TIME_YEAR);
    addHeader( 158, 2, TYPE_INT, HDR_TIME_DAY);
    addHeader( 160, 2, TYPE_INT, HDR_TIME_HOUR);
    addHeader( 162, 2, TYPE_INT, HDR_TIME_MIN);
    addHeader( 164, 2, TYPE_INT, HDR_TIME_SEC);
    addHeader( 166, 2, TYPE_INT, HDR_TIME_CODE);


    addHeader( 114, 2, TYPE_USHORT, HDR_NSAMP );
    addHeader( 116, 2, TYPE_USHORT, HDR_SAMPINT_US );

    addHeader( "ps_gain",       120, 2, TYPE_INT, "P-SEGY gain amplifier" );

    addHeader( "ps_st_name",    180, 5, TYPE_STRING, "P-SEGY station name code" );
    addHeader( "ps_sens_serial",186, 7, TYPE_STRING, "P-SEGY sensor serial code" );
    addHeader( "ps_chan_name",  194, 4, TYPE_STRING, "P-SEGY channel name code" );

    addHeader( "ps_stat_high",  198, 2, TYPE_INT, "P-SEGY static shift [ms] (high 2 bytes)" );
    addHeader( "ps_sampint_us", 200, 4, TYPE_INT, "P-SEGY sample interval [us] (4 bytes)" );
    addHeader( "ps_format",     204, 2, TYPE_INT, "P-SEGY data format (0=16bit, 1=32bit))" );
    addHeader( "ps_time_samp1_ms", 206, 2, TYPE_INT, "P-SEGY time of first sample [ms]" );
    addHeader( "ps_trig_year",  208, 4, TYPE_INT, "P-SEGY trigger time year" );
    addHeader( "ps_trig_day",   210, 4, TYPE_INT, "P-SEGY trigger time julian day" );
    addHeader( "ps_trig_hour",  212, 4, TYPE_INT, "P-SEGY trigger time hours" );
    addHeader( "ps_trig_min",   214, 4, TYPE_INT, "P-SEGY trigger time minutes" );
    addHeader( "ps_trig_sec",   216, 4, TYPE_INT, "P-SEGY trigger time seconds" );
    addHeader( "ps_trig_msec",  218, 4, TYPE_INT, "P-SEGY trigger time milliseconds" );
    addHeader( "ps_scale_fac",  220, 4, TYPE_FLOAT, "P-SEGY scale factor" );
    addHeader( "ps_inst_serial",224, 2, TYPE_INT, "P-SEGY Instrument serial number" );
    addHeader( "ps_nsamp",      228, 4, TYPE_INT, "P-SEGY number of samples (4 bytes)" );
    addHeader( "ps_max_val",    232, 4, TYPE_INT, "P-SEGY Maximum value in counts" );
    addHeader( "ps_min_val",    236, 4, TYPE_INT, "P-SEGY Minimum value in counts" );
  }
  else if( myDefaultMap == SEGY_NODE_OLD ) {
    addHeader( "segy_orig_rec_heading",     30, 2, TYPE_INT, "" );
    addHeader( "segy_orig_tiltx",  94, 2, TYPE_INT, "" );
    addHeader( "segy_orig_tilty",  96, 2, TYPE_INT, "" );
    addHeader( "segy_sub_samp_t0", 104, 2, TYPE_INT, "" );
    addHeader( "segy_comp",        212, 2, TYPE_INT, "Component (1:Z,2:X,3:Y,4:P)" );
    addHeader( "segy_node_id",     214, 2, TYPE_INT, "" );
    addHeader( "segy_station",     216, 2, TYPE_INT, "" );
    addHeader( "segy_geo_heading", 218, 2, TYPE_INT, "" );
    addHeader( "segy_tiltx",       220, 2, TYPE_INT, "" );
    addHeader( "segy_tilty",       222, 2, TYPE_INT, "" );
    addHeader( "segy_sr_azim",     224, 2, TYPE_INT, "" );
    addHeader( "segy_lineno",      228, 4, TYPE_INT, "" );
  }
  else if( myDefaultMap == SEGY_NODE ) {
    addHeader( 118, 2, TYPE_INT, HDR_SEQ );
    addHeader( 120, 2, TYPE_INT, HDR_SAIL_LINE );
    addHeader( 122, 2, TYPE_INT, HDR_SOU_LINE );
    addHeader( 124, 2, TYPE_INT, HDR_GUN_SEQ );
    addHeader( 126, 2, TYPE_INT, HDR_SOU_INDEX );

    addHeader( 128, 4, TYPE_INT, HDR_RCV );
    addHeader( 132, 2, TYPE_INT, HDR_REC_LINE );
    addHeader( "rec_point", 134, 2, TYPE_INT, "Receiver point number" );
    addHeader( 136, 2, TYPE_INT, HDR_REC_INDEX );
    addHeader( 138, 2, TYPE_INT, HDR_NODE );

    addHeader( 168, 4, TYPE_INT, HDR_TIME_USEC );

    addHeader( 172, 2, TYPE_INT, HDR_AN_AZIM );
    addHeader( 174, 2, TYPE_INT, HDR_AN_TILTX );
    addHeader( 176, 2, TYPE_INT, HDR_AN_TILTY );

    addHeader( "an_azim_data", 180, 2, TYPE_INT, TYPE_FLOAT, "Azimuth angle from North (data derived) [deg]" );
    addHeader( "an_tiltx_data", 182, 2, TYPE_INT, TYPE_FLOAT, "Tilt angle, X direction (data derived) [deg]" );
    addHeader( "an_tilty_data", 184, 2, TYPE_INT, TYPE_FLOAT, "Tilt angle, Y direction (data derived) [deg]" );

    addHeader( 186, 2, TYPE_INT, HDR_SR_AZIM );
    addHeader( "stat_res_us", 188, 2, TYPE_INT, "Residual sub-sample static shift [us]" );
    addHeader( 190, 2, TYPE_INT, HDR_SENSOR );
    addHeader( "orient_flag", 192, 2, TYPE_INT, "Re-orientation applied: Yes = 1, No = 0" );

    addHeader( 194, 2, TYPE_INT, HDR_TRC_EDIT );

    addHeader( "sou_x2", 196, 4, TYPE_INT, TYPE_DOUBLE, "Secondary source X position [m]" );
    addHeader( "sou_y2", 200, 4, TYPE_INT, TYPE_DOUBLE, "Secondary source Y position [m]" );
    addHeader( "rec_x2", 204, 4, TYPE_INT, TYPE_DOUBLE, "Secondary receiver X position [m]" );
    addHeader( "rec_y2", 208, 4, TYPE_INT, TYPE_DOUBLE, "Secondary receiver Y position [m]" );

    addHeader( 212, 2, TYPE_INT, HDR_HEADING );

    addHeader( "dc", 236, 4, TYPE_INT, TYPE_FLOAT, "DC bias" );
  }
  else if( myDefaultMap == SEGY_DAS ) {
    removeHeader( HDR_CMP.name );
    removeHeader( HDR_CHAN_ENS.name );
    removeHeader( HDR_STAT_SOU.name );
    removeHeader( HDR_MUTE_START.name );
    removeHeader( HDR_MUTE_END.name );

    addHeader( "das_rec_elev_kb", 20, 4, TYPE_INT, TYPE_FLOAT, "DAS - Kelly Bushing receiver elevation" );
    addHeader( "das_sou_elev_kb", 24, 4, TYPE_INT, TYPE_FLOAT, "DAS - Kelly Bushing source elevation" );
    addHeader( "das_stage", 98, 2, TYPE_SHORT, TYPE_INT, "DAS - Perforation stage number" );

    addHeader( "das_perf", 100, 2, TYPE_SHORT, TYPE_INT, "DAS - Perforation cluster number" );

    addHeader( "das_shot_tvd", 110, 4, TYPE_INT, TYPE_FLOAT, "DAS - True vertical depth of shot below KB" );
    addHeader( 168, 4, TYPE_INT, HDR_TIME_USEC);
    addHeader( 174, 4, TYPE_INT, HDR_SEQ);
    addHeader( "das_acqseq", 178, 4, TYPE_INT, TYPE_INT, "DAS - Acquisition sequence number" );

    addHeader( 180, 4, TYPE_INT, HDR_REC_LINE);
    addHeader( "rec_point", 184, 4, TYPE_INT, "Receiver point number" );
    addHeader( "das_rec_point_scale", 188, 4, TYPE_INT, "Receiver point scalar" );
    addHeader( 192, 4, TYPE_INT, HDR_SOU_LINE);
    addHeader( "sp_num", 196, 4, TYPE_INT, "Shotpoint number (SEGY)");
    addHeader( "das_sp_num_scale", 200, 2, TYPE_INT, "Shotpoint scalar" );
    addHeader( "trace_unit", 202, 2, TYPE_INT, "Trace unit: -1(Other) 0(unknown) 1(Pa) 2(V) 3(mV) 4(A) 5/m 6(m/s) 7(m/s2) 8(N) 9(W)" );
    addHeader( "das_rec_tvd", 204, 4, TYPE_INT, TYPE_FLOAT, "DAS - True vertical depth of receiver below KB" );
    addHeader( "das_rec_tvdss", 208, 4, TYPE_INT, TYPE_FLOAT, "DAS - True vertical depth of receiver below MS" );
    addHeader( "device_ident", 212, 2, TYPE_INT, "Device/trace identifier");
    addHeader( "das_ping_rate", 232, 2, TYPE_INT, "Interrogator unit ping rate [kHz]");
    addHeader( 234, 4, TYPE_INT, HDR_REC_Z);
    addHeader( "das_gauge_len", 238, 2, TYPE_INT, "Interrogator unit gauge length [ns]");
  }
  //  if( myDefaultMap != csSegyHdrMap::NONE && initScalars_in ) {
  if( initScalars_in ) {
    initScalars();  
  }
}
void csSegyHdrMap::dump( FILE* stream ) const {
  //  fprintf(stream," SEGY trace header map\n", info->info->name.c_str(), info->description.c_str() );
  fprintf(stream,"BYTE TYPE     OUT_TYPE NAME                     DESCRIPTION                     # BYTESIZE\n" );
  for( int ihdr = 0; ihdr < numHeaders(); ihdr++ ) {
    csSegyHeaderInfo const* info = header( ihdr );
    std::string text("\"" + info->description + "\"");
    type_t type = info->inType;
    if( type == TYPE_INT && info->byteSize == 2 ) {
      type = TYPE_SHORT;
    }
    fprintf(stream,"%4d %-8s %-8s %-24s %-30s  # %-4d\n", info->byteLoc+1, csGeolibUtils::typeText(type),
            csGeolibUtils::typeText(info->outType),
            info->name.c_str(), text.c_str(), info->byteSize );
  }
}
//--------------------------------------------------------------------------------
//
void csSegyHdrMap::readHdrMapExternalFile( std::string filename ) {
  std::ifstream* fileAscii = new std::ifstream();
  fileAscii->open( filename.c_str(), std::ios::in );
  if( fileAscii->fail() ) {
    delete fileAscii;
    throw( csException("Could not open file '%s'\n", filename.c_str()) );
  }

  myReplaceExistingHeader = true; // Replace existing headers with same name

  char buffer[1024];
  int counterLines = 0;
  while( !fileAscii->getline( buffer, 1024 ).eof() ) {
    if( fileAscii->fail() ) {
      delete fileAscii;
      throw( csException("Unknown error occurred when reading line #%d, input file '%s'", counterLines+1, filename.c_str()) );
    }
    counterLines++;

    cseis_geolib::csVector<std::string> tokenList;
    tokenize( buffer, tokenList );

    if( tokenList.size() == 0 ) continue;
    if( tokenList.size() < 3 ) {
      delete fileAscii;
      throw( csException("Too few columns in line #%d, input file '%s':\n'%s'. Need at least 3.", counterLines+1, filename.c_str(), buffer) );
    }
    int byteLoc  = atoi(tokenList.at(0).c_str()) - 1; // -1: Convert to C style index
    std::string typeNameIn  = tokenList.at(1);
    std::string typeNameOut = tokenList.at(2);
    std::string name        = tokenList.at(3);
    std::string desc        = "";
    if( tokenList.size() > 4 ) desc = tokenList.at(4);
    if( byteLoc < 0 || byteLoc >= 240 ) {  // Internal byte location is between 0-239
      throw( csException("Inconsistent byte location in line #%d, input file '%s':\n'%s'. Must be number between 1-240",
                         counterLines+1, filename.c_str(), buffer) );
    }

    if( !addHeader( byteLoc, typeNameIn, typeNameOut, name, desc ) ) {
      delete fileAscii;
      throw( csException("Cannot add trace header '%s'. Header already exists.\n", name.c_str()) );
    }
  }
  myReplaceExistingHeader = false;
}

bool csSegyHdrMap::addHeader( int byteLoc, std::string const& typeNameIn, std::string const& typeNameOut, std::string const& name, std::string const& desc ) {
  int byteSize = 0;
  type_t typeIn  = TYPE_UNKNOWN;
  type_t typeOut = TYPE_UNKNOWN;

  if( !typeNameIn.compare("int") ) {
    typeIn   = TYPE_INT;
    byteSize = 4;
  }
  else if( !typeNameIn.compare("short") ) {
    typeIn   = TYPE_SHORT;
    byteSize = 2;
  }
  else if( !typeNameIn.compare("ushort") ) {
    typeIn   = TYPE_USHORT;
    byteSize = 2;
  }
  else if( !typeNameIn.compare("float") ) {
    typeIn   = TYPE_FLOAT;
    byteSize = 4;
  }
  else if( !typeNameIn.compare("float_ibm") ) {
    typeIn   = TYPE_FLOAT_IBM;
    byteSize = 4;
  }
  else if( !typeNameIn.compare("double") ) {
    typeIn   = TYPE_DOUBLE;
    byteSize = 8;
  }
  else if( !typeNameIn.compare("2") ) {
    typeIn   = TYPE_SHORT;
    byteSize = 2;
  }
  else if( !typeNameIn.compare("4") ) {
    typeIn   = TYPE_INT;
    byteSize = 4;
  }
  else {
    throw csException("Unknown or unsupported trace header type: %s", typeNameIn.c_str());
  }

  if( !typeNameOut.compare( "int" ) ) {
    typeOut = TYPE_INT;
  }
  else if( !typeNameOut.compare( "float" ) ) {
    typeOut = TYPE_FLOAT;
  }
  else if( !typeNameOut.compare( "double" ) ) {
    typeOut = TYPE_DOUBLE;
  }
  else if( !typeNameOut.compare( "int64" ) ) {
    typeOut = TYPE_INT64;
  }
  else if( !typeNameOut.compare( "string" ) ) {
    typeOut = TYPE_STRING;
    throw csException("String headers are currently not supported");
  }
  else {
    throw csException("Unknown or unsupported trace header type: %s", typeNameOut.c_str());
  }

  if( csStandardHeaders::isStandardHeader( name ) ) {
    csHeaderInfo const* stdInfo = csStandardHeaders::get( name );
    if( typeOut != stdInfo->type ) {
      throw csException("SeaSeis standard trace header with name '%s' exists, but with different type: %s  (user specified type: %s))",
                        name.c_str(), csGeolibUtils::typeText(stdInfo->type), csGeolibUtils::typeText(typeOut) );
    }
    return addHeader( name.c_str(), byteLoc, byteSize, typeIn, stdInfo->type, stdInfo->description.c_str() );
  }
  else {
    return addHeader( name.c_str(), byteLoc, byteSize, typeIn, typeOut, desc.c_str() );
  }
}



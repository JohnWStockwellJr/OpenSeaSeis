/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csException.h"
#include "csVector.h"
#include "csSort.h"
#include "csStandardHeaders.h"
#include <string>

using namespace cseis_geolib;

//----------------------------------------------------------
//
bool csStandardHeaders::isStandardHeader( csHeaderInfo const* info_in ) {
  csHeaderInfo const* info = csStandardHeaders::get( info_in->name );
  if( info == NULL ) {
    return false;
  }  
  else if( info_in == info ) {  // Check if pointers point to same header info object
    return true;
  }
  else if( info->type == info_in->type && !info->description.compare(info_in->description) ) {
    return true;
  }
  else {
    return false;
  }
}
//----------------------------------------------------------
//
bool csStandardHeaders::isStandardHeader( std::string const& name ) {
  if( csStandardHeaders::get( name ) != NULL ) return true;
  return false;
}
//----------------------------------------------------------
//
csHeaderInfo const* csStandardHeaders::get( std::string const& name ) {
  static csVector<csHeaderInfo const*> hdrInfoList(40);
  if( hdrInfoList.size() == 0 ) {
    csStandardHeaders::getAll( &hdrInfoList );
  }
  csHeaderInfo info( TYPE_UNKNOWN, name, "" );
  int pos = csStandardHeaders::findPosition( info.name, &hdrInfoList );
  if( pos >= 0 && pos < hdrInfoList.size() ) {
    return hdrInfoList.at(pos);
  }    
  else {
    return NULL;
  }
}
//----------------------------------------------------------
//
void csStandardHeaders::getAll( cseis_geolib::csVector<csHeaderInfo const*>* hdrList ) {
  if( hdrList == NULL ) throw csException("csStandardHeaders::getStandardHeaders: Passed NULL pointer. This is a program bug in the calling function");
  hdrList->clear();

  cseis_geolib::csVector<csHeaderInfo const*> hdrListTMP(40);

  hdrListTMP.insert( &HDR_AOFFSET );
  hdrListTMP.insert( &HDR_AN_AZIM );
  hdrListTMP.insert( &HDR_AN_INCI );
  hdrListTMP.insert( &HDR_AN_ROLL );
  hdrListTMP.insert( &HDR_AN_TILT );
  hdrListTMP.insert( &HDR_AN_TILTX );
  hdrListTMP.insert( &HDR_AN_TILTY );

  hdrListTMP.insert( &HDR_BIN_X );
  hdrListTMP.insert( &HDR_BIN_Y );
  
  hdrListTMP.insert( &HDR_CBL_AZIM );
  hdrListTMP.insert( &HDR_CCP );
  hdrListTMP.insert( &HDR_CHAN );
  hdrListTMP.insert( &HDR_CHAN_ENS );
  hdrListTMP.insert( &HDR_CMP );
  hdrListTMP.insert( &HDR_CMP_X );
  hdrListTMP.insert( &HDR_CMP_Y );
  hdrListTMP.insert( &HDR_CMP_Z );
  hdrListTMP.insert( &HDR_COL );

  hdrListTMP.insert( &HDR_DATA_TYPE );
  hdrListTMP.insert( &HDR_DC );
  hdrListTMP.insert( &HDR_DELAY_TIME );

  hdrListTMP.insert( &HDR_FFID );
  hdrListTMP.insert( &HDR_FILENO );
  hdrListTMP.insert( &HDR_FILT_HIGH_FREQ );
  hdrListTMP.insert( &HDR_FILT_HIGH_DB );
  hdrListTMP.insert( &HDR_FILT_LOW_FREQ );
  hdrListTMP.insert( &HDR_FILT_LOW_DB );
  hdrListTMP.insert( &HDR_FOLD );
  hdrListTMP.insert( &HDR_FOLD_VERT );
  hdrListTMP.insert( &HDR_GAIN );
  hdrListTMP.insert( &HDR_GAIN_TYPE );
  hdrListTMP.insert( &HDR_GUN_SEQ );
  hdrListTMP.insert( &HDR_HEADING );
  hdrListTMP.insert( &HDR_INCL_I );
  hdrListTMP.insert( &HDR_INCL_C );
  hdrListTMP.insert( &HDR_INCL_V );
  hdrListTMP.insert( &HDR_MUTE_START );
  hdrListTMP.insert( &HDR_MUTE_END );
  hdrListTMP.insert( &HDR_NSAMP );
  hdrListTMP.insert( &HDR_NODE );
  hdrListTMP.insert( &HDR_OFFSET );
  hdrListTMP.insert( &HDR_ORIENT_I );
  hdrListTMP.insert( &HDR_ORIENT_C );
  hdrListTMP.insert( &HDR_ORIENT_V );

  hdrListTMP.insert( &HDR_RCV );
  hdrListTMP.insert( &HDR_REC_INDEX );
  hdrListTMP.insert( &HDR_REC_X );
  hdrListTMP.insert( &HDR_REC_Y );
  hdrListTMP.insert( &HDR_REC_Z );
  hdrListTMP.insert( &HDR_REC_DATUM );
  hdrListTMP.insert( &HDR_REC_ELEV );
  hdrListTMP.insert( &HDR_REC_WDEP );
  hdrListTMP.insert( &HDR_REC_LINE );
  hdrListTMP.insert( &HDR_ROW );

  hdrListTMP.insert( &HDR_SAMPINT_US );
  hdrListTMP.insert( &HDR_SCALAR_ELEV );
  hdrListTMP.insert( &HDR_SCALAR_COORD );
  hdrListTMP.insert( &HDR_SCALAR_STAT );
  hdrListTMP.insert( &HDR_SENSOR );
  hdrListTMP.insert( &HDR_SEQ );
  hdrListTMP.insert( &HDR_SERIAL );
  hdrListTMP.insert( &HDR_SAIL_LINE );
  hdrListTMP.insert( &HDR_SOU_LINE );
  hdrListTMP.insert( &HDR_SOURCE );
  hdrListTMP.insert( &HDR_SOU_INDEX );
  hdrListTMP.insert( &HDR_SOU_X );

  hdrListTMP.insert( &HDR_SOU_Y );
  hdrListTMP.insert( &HDR_SOU_Z );
  hdrListTMP.insert( &HDR_SOU_DATUM );
  hdrListTMP.insert( &HDR_SOU_ELEV );
  hdrListTMP.insert( &HDR_SOU_WDEP );
  hdrListTMP.insert( &HDR_SR_AZIM );

  hdrListTMP.insert( &HDR_STAT_REC );
  hdrListTMP.insert( &HDR_STAT_SOU );
  hdrListTMP.insert( &HDR_STAT_TOT );
  hdrListTMP.insert( &HDR_STAT_RES );

  hdrListTMP.insert( &HDR_TIME_CODE );
  hdrListTMP.insert( &HDR_TIME_YEAR );
  hdrListTMP.insert( &HDR_TIME_DAY );
  hdrListTMP.insert( &HDR_TIME_HOUR );
  hdrListTMP.insert( &HDR_TIME_MIN );
  hdrListTMP.insert( &HDR_TIME_SEC );
  hdrListTMP.insert( &HDR_TIME_MSEC );
  hdrListTMP.insert( &HDR_TIME_USEC );
  hdrListTMP.insert( &HDR_TIME_NANO );
  hdrListTMP.insert( &HDR_TIME_SAMP1 );
  hdrListTMP.insert( &HDR_TIME_SAMP1_US );
  hdrListTMP.insert( &HDR_TAPENO );
  hdrListTMP.insert( &HDR_TRC_TYPE );
  hdrListTMP.insert( &HDR_TRCNO );
  hdrListTMP.insert( &HDR_TRC_EDIT );
  hdrListTMP.insert( &HDR_UNIT_COORD );
  hdrListTMP.insert( &HDR_WDEP );

  int nHdrs = hdrListTMP.size();
  int* idSorted = new int[nHdrs];
  std::string* nameSorted = new std::string[nHdrs];
  for( int ihdr = 0; ihdr < nHdrs; ihdr++ ) {
    idSorted[ihdr]   = ihdr;
    nameSorted[ihdr] =  hdrListTMP.at(ihdr)->name;
  }
  csSort<std::string>().simpleSortIndex( nameSorted, nHdrs, idSorted );
  for( int ihdr = 0; ihdr < nHdrs; ihdr++ ) {
    hdrList->insertEnd( hdrListTMP.at(idSorted[ihdr]) );
  }
  delete [] idSorted;
  delete [] nameSorted;
}
//----------------------------------------------------------
//
int csStandardHeaders::findPosition( std::string const& name, cseis_geolib::csVector<cseis_geolib::csHeaderInfo const*> const* infoList ) {
  int indexLeft  = 0;
  int indexRight = infoList->size() - 1;

  if( infoList->at(indexLeft)->name > name ) {
    return -1;
  }
  else if( infoList->at(indexRight)->name < name ) {
    return infoList->size();
  }
  else if( infoList->at(indexRight)->name == name ) {
    return indexRight;
  }
  else if( infoList->at(indexLeft)->name == name ) {
    return indexLeft;
  }

  while( (indexRight-indexLeft) > 1 ) {
    int atIndex = (indexLeft+indexRight)/2;
    if( infoList->at(atIndex)->name < name ) {
      indexLeft = atIndex;
    }
    else if( infoList->at(atIndex)->name > name ) {
      indexRight = atIndex;
    }
    else {
      return atIndex;
    }
  }
  return -1;
}


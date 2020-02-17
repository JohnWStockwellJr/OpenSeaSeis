/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csStandardSegdHeader.h"
#include "csStandardHeaders.h"
#include "csVector.h"
#include "csHeaderInfo.h"
#include "csException.h"

using namespace cseis_geolib;
using namespace cseis_segd;

void csStandardSegdHeader::setStandardHeaders( int segdHeaderMap, cseis_geolib::csVector<cseis_geolib::csHeaderInfo const*>* stdHdrList ) {
  if( segdHeaderMap == HDR_MAPPING_STANDARD || segdHeaderMap == HDR_MAPPING_OBC || segdHeaderMap == HDR_MAPPING_TOWED ) {
    stdHdrList->insertEnd( &HDR_TRCNO );  // Trace seq num within line
    stdHdrList->insertEnd( &HDR_FFID );   // Original field record number
    stdHdrList->insertEnd( &HDR_CHAN );  // Trace num within original field record
    stdHdrList->insertEnd( &HDR_SOURCE );   // Source point num
    stdHdrList->insertEnd( &HDR_CMP );   // Ensemble num (CMP...)
    stdHdrList->insertEnd( &HDR_CHAN_ENS ); // Trace num within ensemble
    stdHdrList->insertEnd( &HDR_TRC_TYPE ); // Trace identification code / trace type
//    stdHdrList->insertEnd( new standardHeaderStruct(30,2,HDR_ID_FOLD_VERT ); // Num vertically summed traces
    stdHdrList->insertEnd( &HDR_FOLD );      // Num horizontally summed fold ('normal' fold)
    stdHdrList->insertEnd( &HDR_DATA_TYPE ); // Data use/data type (1:production, 2:test)
    stdHdrList->insertEnd( &HDR_OFFSET );    // Offset
    stdHdrList->insertEnd( &HDR_REC_ELEV );  // Rcv elevation
    stdHdrList->insertEnd( &HDR_SOU_ELEV );  // Surface elevation at source
    stdHdrList->insertEnd( &HDR_SOU_Z );     // Source depth below surface
    stdHdrList->insertEnd( &HDR_REC_DATUM );  // Datum elevation at rcv
    stdHdrList->insertEnd( &HDR_SOU_DATUM );  // Datum elevation at src
    stdHdrList->insertEnd( &HDR_REC_WDEP );  // Water depth at rcv
    stdHdrList->insertEnd( &HDR_SOU_WDEP );  // Water depth at src
    stdHdrList->insertEnd( &HDR_SCALAR_ELEV ); // Elevation scalar
    stdHdrList->insertEnd( &HDR_SCALAR_COORD );// Coordinates scalar

    stdHdrList->insertEnd( &HDR_SOU_X );
    stdHdrList->insertEnd( &HDR_SOU_Y );
    stdHdrList->insertEnd( &HDR_REC_X );
    stdHdrList->insertEnd( &HDR_REC_Y );
    stdHdrList->insertEnd( &HDR_UNIT_COORD );

    stdHdrList->insertEnd( &HDR_MUTE_START );
    stdHdrList->insertEnd( &HDR_MUTE_END );
    stdHdrList->insertEnd( &HDR_NSAMP );
    stdHdrList->insertEnd( &HDR_SAMPINT_US );

    stdHdrList->insertEnd( &HDR_TIME_YEAR );
    stdHdrList->insertEnd( &HDR_TIME_DAY );
    stdHdrList->insertEnd( &HDR_TIME_HOUR );
    stdHdrList->insertEnd( &HDR_TIME_MIN );
    stdHdrList->insertEnd( &HDR_TIME_SEC );
    stdHdrList->insertEnd( &HDR_TIME_NANO );
    stdHdrList->insertEnd( &HDR_TIME_CODE );
  }
  else {
    throw( csException("csStandardSegdHeader::setStandardHeaders: Standard header set not supported (%d)", segdHeaderMap) );
  }

  if( segdHeaderMap == HDR_MAPPING_OBC ) {
    stdHdrList->insertEnd( &HDR_SOU_LINE );
    stdHdrList->insertEnd( &HDR_REC_LINE );
    stdHdrList->insertEnd( &HDR_GUN_SEQ );
    stdHdrList->insertEnd( &HDR_SENSOR );
    stdHdrList->insertEnd( &HDR_RCV );
    stdHdrList->insertEnd( &HDR_SERIAL );
    stdHdrList->insertEnd( &HDR_SEQ );
    stdHdrList->insertEnd( &HDR_INCL_I );
    stdHdrList->insertEnd( &HDR_INCL_C );
    stdHdrList->insertEnd( &HDR_INCL_V );
    stdHdrList->insertEnd( &HDR_ORIENT_I );
    stdHdrList->insertEnd( &HDR_ORIENT_C );
    stdHdrList->insertEnd( &HDR_ORIENT_V );
    stdHdrList->insertEnd( &HDR_CBL_AZIM );
  }
}


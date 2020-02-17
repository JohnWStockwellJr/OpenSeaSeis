/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csNavHeader.h"
#include "csSegdFunctions.h"
#include "csException.h"
#include <ostream>
#include <cstring>

using namespace cseis_segd;
using namespace cseis_geolib;
using std::endl;

csNavHeader::csNavHeader( int navHeaderID ) {
  myNavHeaderID = navHeaderID;
  myStructPtr = NULL;
  mySize      = 0;
  switch( myNavHeaderID ) {
    case NAV_HEADER_PSI:
      myStructPtr = new NavHeader_PSI_Struct();
      mySize = 122;  // May also be 115 bytes..?
      break;
    case NAV_HEADER_LABO:
      mySize = 166;
      myStructPtr = new NavHeader_LABO_Struct();
      break;
    case NAV_HEADER_HYDRONAV_VER6:
      mySize = 176 + 16 + 16;  // 16 header description block + 16 bytes GCS90 tag header (should be put somewhere else...)
      myStructPtr = new NavHeader_HYDRONAV_ver6_Struct();
      break;
    case NAV_HEADER_HYDRONAV_VER1:
      mySize = 176;
      break;
    case NAV_HEADER_NONE:
      break;
    default:
      throw( csException("Unknown nav header ID. Cannot generate navigation header object") );
  }
}
csNavHeader::~csNavHeader() {
  if( myStructPtr != NULL ) {
    switch( myNavHeaderID ) {
      case NAV_HEADER_PSI:
        delete (NavHeader_PSI_Struct*)myStructPtr;
        break;
      case NAV_HEADER_LABO:
        delete (NavHeader_PSI_Struct*)myStructPtr;
        break;
      case NAV_HEADER_HYDRONAV_VER6:
        delete (NavHeader_HYDRONAV_ver6_Struct*)myStructPtr;
        break;
      case NAV_HEADER_HYDRONAV_VER1:
        delete (NavHeader_HYDRONAV_ver1_Struct*)myStructPtr;
    }
    myStructPtr = NULL;
  }
}
void csNavHeader::extractHeaders( byte const* buffer, int numBytes ) {
  if( numBytes < mySize ) return;

  if( myNavHeaderID == NAV_HEADER_LABO || myNavHeaderID == NAV_HEADER_PSI ) {
    NavHeader_PSI_Struct* psi;
    if( myNavHeaderID == NAV_HEADER_LABO ) {
      psi = &( reinterpret_cast<NavHeader_LABO_Struct*>(myStructPtr)->psi );
    }
    else {
      psi = reinterpret_cast<NavHeader_PSI_Struct*>(myStructPtr);
    }
    psi->headerID[2]  = '\0';
    psi->programRevision[4] = '\0';
    psi->shotTime[21] = '\0';
    psi->timeRef[3]   = '\0';
    psi->lineName[16] = '\0';
  
    memcpy( psi->headerID, &buffer[0], 2 );
    psi->length = string2int( &buffer[2], 4 );
    memcpy( psi->programRevision, &buffer[6], 4 );
    psi->lineStatus = string2int( &buffer[10], 2 );
    memcpy( psi->shotTime, &buffer[12], 21 );  // 2+2+2+1+6+4+2+2
    memcpy( psi->timeRef, &buffer[33], 3 );
    psi->shotNumber = string2int( &buffer[36], 6 );
    memcpy( psi->lineName, &buffer[42], 16 );
    psi->masterLatitude  = string2double( &buffer[58], 11 );
    psi->masterLongitude = string2double( &buffer[69], 11 );
    psi->waterDepth      = string2float( &buffer[80], 6 );
    psi->sourceLatitude  = string2double( &buffer[86], 11 );
    psi->sourceLongitude = string2double( &buffer[97], 11 );
    psi->masterGyro      = string2float( &buffer[108], 5 );
    psi->masterCMG       = string2float( &buffer[113], 5 );
    psi->masterSpeed     = string2float( &buffer[118], 4 );
    
    if( myNavHeaderID == NAV_HEADER_LABO ) {
      NavHeader_LABO_Struct* labo;
      labo = reinterpret_cast<NavHeader_LABO_Struct*>(myStructPtr);
      labo->masterEasting   = string2double( &buffer[122], 11 );
      labo->masterNorthing  = string2double( &buffer[133], 11 );
      labo->sourceDeltaEasting  = string2double( &buffer[144], 7 );
      labo->sourceDeltaNorthing = string2double( &buffer[151], 7 );
      labo->lineBearing     = string2float( &buffer[158], 5 );
      labo->julianDay       = string2int( &buffer[163], 3 );
    }
  }
  else if( myNavHeaderID == NAV_HEADER_HYDRONAV_VER6 ) {
    NavHeader_HYDRONAV_ver6_Struct* hydro = reinterpret_cast<NavHeader_HYDRONAV_ver6_Struct*>(myStructPtr);
    hdrDesc_HYDRONAV_Struct* hdrDesc = &hydro->hdrDesc;
    tagDesc_HYDRONAV_Struct* tagDesc = &hydro->tagDesc;
    
    hdrDesc->version    = UINT8(&buffer[0]);
    hdrDesc->revision   = UINT8(&buffer[1]);
    hdrDesc->length     = UINT16(&buffer[2]);
    hdrDesc->numTagBlocks = UINT8(&buffer[4]);
    hdrDesc->quality    = UINT8(&buffer[5]);

    int bytePos = 16;
    memcpy( tagDesc->desc, &buffer[bytePos+0], 8 );
    tagDesc->desc[8] = '\0';
    tagDesc->version    = UINT8(&buffer[bytePos+8]);
    tagDesc->revision   = UINT8(&buffer[bytePos+9]);
    tagDesc->length     = UINT16(&buffer[bytePos+10]);
    tagDesc->dateFormat = UINT8(&buffer[bytePos+12]);

/*
    tag_HYDRONAV_Struct* tag;
    if( myNavHeaderID == NAV_HEADER_HYDRONAV ) {
      tag = &( reinterpret_cast<NavHeader_HYDRONAV_Struct*>(myStructPtr)->tag );
    }
    else {
      tag = &( reinterpret_cast<NavHeader_HYDRONAV716_Struct*>(myStructPtr)->tag );
    }
    memcpy( tag->desc, &buffer[0], 8 );
    tag->desc[8] = '\0';
    tag->version  = buffer[8];
    tag->revision = buffer[9];
    memcpy( tag->len, &buffer[10], 2 );
    tag->len[2] = '\0';
    tag->dateFormat = buffer[12];
    tag->null[0]     = buffer[13];
    tag->null[1]     = buffer[14];
    tag->eot         = buffer[15];
*/
//    NavHeader_HYDRONAV_Struct* hydro = reinterpret_cast<NavHeader_HYDRONAV_Struct*>(myStructPtr);
    memcpy( hydro->header, &buffer[bytePos+16], 8 );
    hydro->header[8] = '\0';
    memcpy( hydro->time, &buffer[bytePos+24], 15 );
    hydro->time[15]  = '\0';
    hydro->GPSclock  = string2int( &buffer[bytePos+40], 1 );
    memcpy( hydro->vesselId, &buffer[bytePos+42], 8 );
    hydro->vesselId[8] = '\0';

    hydro->stationNum  = string2int( &buffer[bytePos+51], 7 );
    hydro->FSID        = string2int( &buffer[bytePos+59], 10 );

    hydro->northing    = string2double( &buffer[bytePos+70], 11 );
    hydro->easting     = string2double( &buffer[bytePos+83], 11 );

    hydro->vesselSpeed = string2float( &buffer[bytePos+96], 6 );
    hydro->bearing     = string2float( &buffer[bytePos+105], 6 );
    hydro->waterDepth  = string2float( &buffer[bytePos+114], 5 );
    memcpy( hydro->date, &buffer[bytePos+121], 10 );
    hydro->date[10] = '\0';

    hydro->distOffline = string2int( &buffer[bytePos+132], 5 );
    hydro->distInline  = string2int( &buffer[bytePos+139], 8 );
    hydro->lineStatus  = string2int( &buffer[bytePos+149], 2 );

    hydro->GPSHeight = string2float( &buffer[bytePos+151], 5 );
    memcpy( hydro->lineName, &buffer[bytePos+158], 16 );
    hydro->lineName[16] = '\0';
    hydro->eod = buffer[bytePos+174];
    hydro->eot = buffer[bytePos+175];

    bytePos = 16 + 176;
    tagDesc = &hydro->tagDescGCS90;
    memcpy( tagDesc->desc, &buffer[bytePos+0], 8 );
    tagDesc->desc[8] = '\0';
    tagDesc->version    = UINT8(&buffer[bytePos+8]);
    tagDesc->revision   = UINT8(&buffer[bytePos+9]);
    tagDesc->length     = UINT16(&buffer[bytePos+10]);
    tagDesc->dateFormat = UINT8(&buffer[bytePos+12]);
  }
}

void csNavHeader::dump( std::ostream& os ) {
  if( myNavHeaderID == NAV_HEADER_LABO || myNavHeaderID == NAV_HEADER_PSI ) {
    NavHeader_PSI_Struct* psi;
    if( myNavHeaderID == NAV_HEADER_LABO ) {
      psi = &( reinterpret_cast<NavHeader_LABO_Struct*>(myStructPtr)->psi );
    }
    else {
      psi = reinterpret_cast<NavHeader_PSI_Struct*>(myStructPtr);
    }
    os <<
      "headerID             : " << psi->headerID << endl <<
      "length of header(excl):" << psi->length << endl <<
      "programRevision      : " << psi->programRevision << endl <<
      "lineStatus           : " << psi->lineStatus << endl <<
      "shotTime             : " << psi->shotTime << endl <<
      "timeRef              : " << psi->timeRef << endl <<
      "shotNumber           : " << psi->shotNumber << endl <<
      "lineName             : " << psi->lineName << endl <<
      "masterLatitude       : " << psi->masterLatitude << endl <<
      "masterLongitude      : " << psi->masterLongitude << endl <<
      "waterDepth           : " << psi->waterDepth << endl <<
      "sourceLatitude       : " << psi->sourceLatitude << endl <<
      "sourceLongitude      : " << psi->sourceLongitude << endl <<
      "masterGyro           : " << psi->masterGyro << endl <<
      "masterCMG            : " << psi->masterCMG << endl <<
      "masterSpeed          : " << psi->masterSpeed << endl;
    if( myNavHeaderID == NAV_HEADER_LABO ) {
      NavHeader_LABO_Struct* labo = reinterpret_cast<NavHeader_LABO_Struct*>(myStructPtr);
      os <<
        "masterEasting        : " << labo->masterEasting << endl <<
        "masterNorthing       : " << labo->masterNorthing << endl <<
        "sourceDeltaEasting   : " << labo->sourceDeltaEasting << endl <<
        "sourceDeltaNorthing  : " << labo->sourceDeltaNorthing << endl <<
        "lineBearing          : " << labo->lineBearing << endl <<
        "julianDay            : " << labo->julianDay << endl;
    }
  }
  else if( myNavHeaderID == NAV_HEADER_HYDRONAV_VER6 ) {
    NavHeader_HYDRONAV_ver6_Struct* hydro = reinterpret_cast<NavHeader_HYDRONAV_ver6_Struct*>(myStructPtr);
    hdrDesc_HYDRONAV_Struct* hdrDesc = &hydro->hdrDesc;
    tagDesc_HYDRONAV_Struct* tagDesc = &hydro->tagDesc;

    os <<
      "Version            : " << hdrDesc->version << endl <<
      "Revision           : " << hdrDesc->revision << endl <<
      "Length             : " << hdrDesc->length << endl <<
      "Num tag blocks     : " << hdrDesc->numTagBlocks << endl <<
      "Quality            : " << hdrDesc->quality << endl;

    os <<
      "Tag description    : " << tagDesc->desc << endl <<
      "Version            : " << tagDesc->version << endl <<
      "Revision           : " << tagDesc->revision << endl <<
      "Length             : " << tagDesc->length << endl <<
      "Date format        : " << tagDesc->dateFormat << endl;

      os <<
        "Header      : " << hydro->header << endl <<
        "Time        : " << hydro->time << endl <<
        "GPS clock   : " << hydro->GPSclock << endl <<
        "Vessel ID   : " << hydro->vesselId << endl <<
        "Station num : " << hydro->stationNum << endl <<
        "FSID        : " << hydro->FSID << endl <<
        "Northing    : " << hydro->northing << endl <<
        "Easting     : " << hydro->easting << endl <<
        "Vessel speed: " << hydro->vesselSpeed << endl <<
        "Bearing     : " << hydro->bearing << endl <<
        "Water depth : " << hydro->waterDepth << endl <<
        "Date        : " << hydro->date << endl <<
        "Dist Offline: " << hydro->distOffline << endl <<
        "Dist line   : " << hydro->distInline << endl <<
        "Line status : " << hydro->lineStatus << endl <<        
        "GSP height  : " << hydro->GPSHeight << endl <<
        "Line name   : " << hydro->lineName << endl <<
        "EOD         : " << hydro->eod << endl <<
        "EOT         : " << hydro->eot << endl;
    tagDesc = &hydro->tagDescGCS90;
    os << endl <<
      "Tag description    : " << tagDesc->desc << endl <<
      "Version            : " << tagDesc->version << endl <<
      "Revision           : " << tagDesc->revision << endl <<
      "Length             : " << tagDesc->length << endl <<
      "Date format        : " << tagDesc->dateFormat << endl;

  }
  else {
    os << "..nav header not recognized." << endl;
    return;
  }
}

void csNavHeader::setCommonRecordHeaders( commonRecordHeaderStruct& comRecHdr ) {
  if( myNavHeaderID == NAV_HEADER_HYDRONAV_VER6 ) {
    NavHeader_HYDRONAV_ver6_Struct* hydro = reinterpret_cast<NavHeader_HYDRONAV_ver6_Struct*>(myStructPtr);
//    hdrDesc_HYDRONAV_Struct* hdrDesc = &hydro->hdrDesc;
//    tagDesc_HYDRONAV_Struct* tagDesc = &hydro->tagDesc;

      
    // hydro->vesselSpeed = string2float( &buffer[bytePos+96], 6 );
    // hydro->bearing     = string2float( &buffer[bytePos+105], 6 );
    //    hydro->waterDepth
    comRecHdr.shotNum    = hydro->stationNum;
    comRecHdr.srcEasting = hydro->easting;
    comRecHdr.srcNorthing = hydro->northing;
    //    comRecHdr.srcElev    = hydro->;
  }
}



/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csNavInterface.h"
#include "csSegdFunctions.h"
#include <ostream>
#include <iostream>
#include <cstring>

using namespace cseis_segd;
using namespace std;

csNavInterface::csNavInterface( int interfaceID ) {
  myInterfaceID = interfaceID;

  blockCompass1 = NULL;
  blockBird1    = NULL;
  blockBird2    = NULL;

  myHeaderSeemsOK = true;
}
csNavInterface::~csNavInterface() {
  if( blockCompass1 != NULL ) {
    delete [] blockCompass1;
    blockCompass1 = NULL;
  }
  if( blockBird1 != NULL ) {
    delete [] blockBird1;
    blockBird1 = NULL;
  }
  if( blockBird2 != NULL ) {
    delete [] blockBird2;
    blockBird2 = NULL;
  }
}

void csNavInterface::extractHeaders( byte const* buffer, int numBytes ) {
  if( numBytes < 4 ) return;
  int byteLoc = 0;
  if( myInterfaceID == CM_DIGI_COMP_A ) {
    time[8] = '\0';
    memcpy( time, &buffer[0], 8 );
    eventNumber  = string2int( &buffer[8], 5 );
    numCompasses = string2int( &buffer[13], 2 );
    blockCompass1 = new compass1Struct[numCompasses];
    byteLoc = 15;
    for( int i = 0; i < numCompasses; i++ ) {
      blockCompass1[i].status = buffer[byteLoc];
      blockCompass1[i].unitNumber = (short)string2int( &buffer[byteLoc+1], 2 );
      blockCompass1[i].heading    = (float)(string2int( &buffer[byteLoc+3], 4 ))*0.1;
      byteLoc += 7;
    }
    byteLoc += 2; // ???
    numBirds = string2int( &buffer[byteLoc], 2 );
    blockBird1 = new bird1Struct[numBirds];
    byteLoc += 2;
    for( int i = 0; i < numBirds; i++ ) {
      blockBird1[i].birdStatus = buffer[byteLoc];
      blockBird1[i].tempStatus = buffer[byteLoc+1];
      blockBird1[i].unitNumber = (short)string2int( &buffer[byteLoc+2], 2 );
      blockBird1[i].depth      = ((float)(string2int( &buffer[byteLoc+4], 4 )))*0.01;
      blockBird1[i].finAngle   = ((float)(string2int( &buffer[byteLoc+8], 4 ))*0.1) - 18.5;  // degrees
      blockBird1[i].temperature= ((float)(string2int( &buffer[byteLoc+12], 4 )))*0.1; // deg
      byteLoc += 16;
    }
  }
  else if( myInterfaceID == CM_DIGI_PSIB || myInterfaceID == CM_DIGI_TS ) {
    syncCharacters[4] = '\0';
    memcpy( syncCharacters, &buffer[0], 4 );
    byteCountAll = UINT16(&buffer[4]);
    if( numBytes < byteCountAll ) return;
    
    blockBird2 = new bird2Struct();
  
    int byteLoc = 6;
    int byteCount = UINT16_swap(&buffer[byteLoc]);  // Swap LSB and MSB first
    numBirds = ( byteCount - 4 ) / 5;  // 4bytes = 2b byteCount + 1b messageType + 1b checkByte,  5bytes = remainder
  // Heading
    int nBytes = extractBlock( HEADING, &blockBird2->blockHeading, &buffer[byteLoc] );
    byteLoc += 3+nBytes;
  
  // Depth
    nBytes = extractBlock( DEPTH, &blockBird2->blockDepth, &buffer[byteLoc] );
    byteLoc += 3+nBytes;
  
  // FinAngle
    nBytes = extractBlock( FINANGLE, &blockBird2->blockFinAngle, &buffer[byteLoc] );
    byteLoc += 3+nBytes;
  }
  mySize = byteLoc;
}

int csNavInterface::extractBlock( int blockType, bird2HeaderStruct* headerBlock, byte const* buffer ) {
  headerBlock->byteCount = UINT16_swap(&buffer[0]);  // Swap LSB and MSB first
  headerBlock->messageType = buffer[2];
  bird2DataStruct* block = new bird2DataStruct[numBirds];
  headerBlock->block = block;

  int byteLoc = 3;
  char tmp;
  for( int ibird = 0; ibird < numBirds; ibird++ ) {
    block[ibird].unitNumber = buffer[byteLoc++];
    if( block[ibird].unitNumber > 64 && block[ibird].unitNumber < 127 ) {
      block[ibird].unitNumber -= 64;
    }
    else {
      block[ibird].unitNumber -= 160;
    }
    block[ibird].cableNumber  = buffer[byteLoc++] + 1;
    block[ibird].vesselNumber = buffer[byteLoc++] + 1;
    block[ibird].valueDecimal = buffer[byteLoc++];
    switch( blockType ) {
      case HEADING:
        block[ibird].value        = UINT8(&buffer[byteLoc]);       // Bits 0-3
        block[ibird].status       = (buffer[byteLoc] & 0x8) >> 3;  // Bit 4
        break;
      case DEPTH:
        block[ibird].status       = (buffer[byteLoc] & 0x1);    // Bit 7
        if( block[ibird].status == 0 ) { // feet
          tmp = (buffer[byteLoc] & 0x10) >> 4;  // Bit 3
          block[ibird].value     = tmp; //((int)tmp == 0) ? 1 : -1;
        }
        else { // meter
          tmp = buffer[byteLoc] & 0xd0;
          block[ibird].value      = UINT8( &tmp );  // Bits 0-2
          tmp = (buffer[byteLoc] & 010) >> 4;
          block[ibird].value      *= ((int)tmp == 0) ? 1 : -1;
        }
        break;
      case FINANGLE:
        tmp = (buffer[byteLoc] & 0x10) >> 4;
        block[ibird].value        = (buffer[byteLoc] & 0x10) >> 4; // Bit 3, sign
        block[ibird].status       = (buffer[byteLoc] & 0x8) >> 3;      // Bit 4, status
        break;
    }
    byteLoc += 1;
  }
  headerBlock->checkByte = buffer[byteLoc++];

  return byteLoc-3;
}

void csNavInterface::dump( std::ostream& os ) {
  if( myInterfaceID == CM_DIGI_COMP_A ) {
    os <<
      "Time             : " << time << endl <<
      "Event number     : " << eventNumber << endl <<
      " # Compasses     : " << numCompasses << endl;
    for( int i = 0; i < numCompasses; i++ ) {
      os << 
        "Compass #" << (i+1) << endl <<
        " Status(C:good,c:bad): " << blockCompass1[i].status << endl <<
        " Unit number         : " << blockCompass1[i].unitNumber << endl <<
        " Heading [deg]       : " << blockCompass1[i].heading << endl;
    }
    os << endl << " # Birds     : " << numBirds << endl;
    for( int i = 0; i < numBirds; i++ ) {
      os << 
        "Bird #" << (i+1) << endl <<
        " Status(B:good,b:bad): " << blockBird1[i].birdStatus << endl <<
        " Status(T:good,t:bad): " << blockBird1[i].tempStatus << endl <<
        " Unit number         : " << blockBird1[i].unitNumber << endl <<
        " Depth [m]           : " << blockBird1[i].depth << endl <<
        " Fin Angle [deg]     : " << blockBird1[i].finAngle << endl <<
        " Temperature [deg]   : " << blockBird1[i].temperature << endl;
    }
  }
  else if( myInterfaceID == CM_DIGI_PSIB ||
           myInterfaceID == CM_DIGI_TS ) {
      os <<
        "Sync characters      : " << syncCharacters << endl <<
        "Byte count (all)     : " << byteCountAll << endl <<
        " *** COMPASS Heading block ******************************" << endl;
      csNavInterface::bird2HeaderStruct const* b = &blockBird2->blockHeading;
      os << " Message type        : " << b->messageType << endl <<  // should be displayed hexadecimal
        " Byte count          : " << b->byteCount << endl;
      for( int ibird = 0; ibird < numBirds; ibird++ ) {
        os <<
          " Bird #  " << (ibird+1) << " : " << endl <<
          "  Unit number          : " << b->block[ibird].unitNumber << endl <<
          "  Cable number         : " << b->block[ibird].cableNumber << endl <<
          "  Vessel number        : " << b->block[ibird].vesselNumber << endl <<
          "  Heading [0.1deg]     : " << b->block[ibird].valueDecimal << endl <<
          "  Heading [deg]        : " << b->block[ibird].value << endl <<
          "  Status               : " << b->block[ibird].status << endl;
      }
      os << "Check byte           : " << b->checkByte << endl;

      b = &blockBird2->blockDepth;
      os <<
        " *** COMPASS Bird Depth block ************************" << endl <<
        " Message type        : " << b->messageType << endl <<  // should be displayed hexadecimal
        " Byte count          : " << b->byteCount << endl;
      for( int ibird = 0; ibird < numBirds; ibird++ ) {
        os <<
          " Bird #  " << (ibird+1) << " : " << endl <<
          "  Unit number          : " << b->block[ibird].unitNumber << endl <<
          "  Cable number         : " << b->block[ibird].cableNumber << endl <<
          "  Vessel number        : " << b->block[ibird].vesselNumber << endl <<
          "  Depth, [ft]or[0.1m]  : " << b->block[ibird].valueDecimal << endl <<
          "  ft:Sign, m:Depth     : " << b->block[ibird].value << endl <<
          "  Unit (0=ft, 1=m)     : " << b->block[ibird].status << endl;
      }
      os << "Check byte           : " << b->checkByte << endl;
      b = &blockBird2->blockFinAngle;
      os <<
        " *** COMPASS Fin Angle block ************************" << endl <<
        " Message type        : " << b->messageType << endl <<  // should be displayed hexadecimal
        " Byte count          : " << b->byteCount << endl;
      for( int ibird = 0; ibird < numBirds; ibird++ ) {
        os <<
          " Bird #  " << (ibird+1) << " : " << endl <<
          "  Unit number          : " << b->block[ibird].unitNumber << endl <<
          "  Cable number         : " << b->block[ibird].cableNumber << endl <<
          "  Vessel number        : " << b->block[ibird].vesselNumber << endl <<
          "  FinAngle [0.1deg]    : " << b->block[ibird].valueDecimal << endl <<
          "  Sign (0:pos, 1:neg)  : " << b->block[ibird].value << endl <<
          "  Status (0:ok,1:bad)  : " << b->block[ibird].status << endl;
      }
      os << "Check byte           : " << b->checkByte << endl;
  }
  else {
    os << "Unknown nav interface type..." << endl;
  }
}



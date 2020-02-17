/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_NAV_INTERFACE_H
#define CS_NAV_INTERFACE_H

#include "csSegdDefines.h"
#include <ostream>

namespace cseis_segd {

/**
 * SEGD Navigation interface
 * @author Bjorn Olofsson
 * @date 2006
 */
class csNavInterface {
public:
  csNavInterface( int interfaceID );
  ~csNavInterface();

  int numCompasses;
  int numBirds;

  // CM_DIGI_COMP_A :
  char time[8+1];
  int  eventNumber;
  struct compass1Struct {
    char  status;
    short unitNumber;
    float heading;
  };
  struct bird1Struct {
    char  birdStatus;
    char  tempStatus;  // Temperature
    short unitNumber;
    float depth;
    float finAngle;
    float temperature;
  };
  compass1Struct* blockCompass1;
  bird1Struct*    blockBird1;

  //----------------------------------------------------
  // CM_DIGI_PSIB, CM_DIGI_PS
  //
  char  syncCharacters[4+1];
  short byteCountAll;

  struct bird2DataStruct {
    short unitNumber;
    short cableNumber;
    short vesselNumber;
    short valueDecimal;
    short value;
    short status;
  };
  struct bird2HeaderStruct {
    short byteCount;
    short messageType;
    short checkByte;
    bird2DataStruct* block;
  };
  struct bird2Struct {
    bird2HeaderStruct blockHeading;
    bird2HeaderStruct blockDepth;
    bird2HeaderStruct blockFinAngle;
  };
  bird2Struct* blockBird2;

  //-------------------------------------
  //
  void extractHeaders( byte const* buffer, int numBytes );
  void dump( std::ostream& os );
  int size() const { return mySize; }
private:
  bool myHeaderSeemsOK;
  csNavInterface( csNavInterface const& obj );
  int extractBlock( int blockType, bird2HeaderStruct* headerBlock, byte const* buffer );
  int myInterfaceID; // e.g. CM_DIGI_COMP_A
  int mySize;

  static int const HEADING  = 44;
  static int const DEPTH    = 45;
  static int const FINANGLE = 46;
};

} // namespace
#endif


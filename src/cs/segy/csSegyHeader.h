/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGY_HEADER_H
#define CS_SEGY_HEADER_H

namespace cseis_geolib {
  
/**
* Segy Header
* @author Bjorn Olofsson
* @date 2006
*/
class csSegyHeader {
public:
  static int const SIZE_CHARHDR;// = 3200;
  static int const SIZE_TRCHDR;//  = 240;
  static int const SIZE_BINHDR;//  = 400;

  static int const AUTO;// = -1;
  static int const DATA_FORMAT_IEEE;//  = 5;
  static int const DATA_FORMAT_IBM;//   = 1;
  static int const DATA_FORMAT_INT32;// = 2;
  static int const DATA_FORMAT_INT16;// = 3;

  static int const BYTE_LOC_SCALAR_ELEV;//  = 68;
  static int const BYTE_LOC_SCALAR_COORD;// = 70;
  static int const BYTE_LOC_SCALAR_STAT;//  = 214;

};  // END class

}  // END namespace

#endif


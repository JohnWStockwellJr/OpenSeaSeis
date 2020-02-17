/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSegyHeader.h"

 int const cseis_geolib::csSegyHeader::SIZE_CHARHDR = 3200;
 int const cseis_geolib::csSegyHeader::SIZE_TRCHDR  = 240;
 int const cseis_geolib::csSegyHeader::SIZE_BINHDR  = 400;

 int const cseis_geolib::csSegyHeader::AUTO = -1;
 int const cseis_geolib::csSegyHeader::DATA_FORMAT_IEEE  = 5;
 int const cseis_geolib::csSegyHeader::DATA_FORMAT_IBM   = 1;
 int const cseis_geolib::csSegyHeader::DATA_FORMAT_INT32 = 2;
 int const cseis_geolib::csSegyHeader::DATA_FORMAT_INT16 = 3;

 int const cseis_geolib::csSegyHeader::BYTE_LOC_SCALAR_ELEV  = 68;
 int const cseis_geolib::csSegyHeader::BYTE_LOC_SCALAR_COORD = 70;
 int const cseis_geolib::csSegyHeader::BYTE_LOC_SCALAR_STAT  = 214;

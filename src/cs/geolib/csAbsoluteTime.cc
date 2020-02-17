/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csAbsoluteTime.h"
#include <cstring>

using namespace cseis_geolib;

csAbsoluteTime::csAbsoluteTime() {
  myTime_days = 0;
  myTime_sec  = 0;
  myTime_us   = 0;
}

csAbsoluteTime::csAbsoluteTime( byte_t const* buffer) {
  set( buffer );
}
void csAbsoluteTime::set( byte_t const* buffer ) {
  memcpy( &myTime_days, &buffer[0], 4 );
  memcpy( &myTime_sec,  &buffer[4], 4 );
  memcpy( &myTime_us,   &buffer[8], 4 );
}
csAbsoluteTime::csAbsoluteTime( int days, int sec, int us ) {
  set( days, sec, us );
}

csAbsoluteTime::~csAbsoluteTime() {
}

void csAbsoluteTime::get( int& days, int& sec, int& us ) const {
  days = myTime_days;
  sec  = myTime_sec;
  us   = myTime_us;
}
void csAbsoluteTime::set( int days, int sec, int us ) {
  myTime_days = days;
  myTime_sec  = sec;
  myTime_us   = us;
}

void csAbsoluteTime::set( csAbsoluteTime const* obj ) {
  set( obj->myTime_days, obj->myTime_sec, obj->myTime_us );
}

std::string csAbsoluteTime::getDate() const {
  std::string text = "DATE CONVERSION NOT IMPLEMENTED YET";
  return text;
}

int csAbsoluteTime::writeBytes( byte_t* buffer ) const {
  memcpy( &buffer[0], &myTime_days, 4 );
  memcpy( &buffer[4], &myTime_sec, 4 );
  memcpy( &buffer[8], &myTime_us, 4 );
  return 12;
}


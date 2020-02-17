/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TIME_H
#define CS_TIME_H

#include "geolib_defines.h"
#include "csGeolibUtils.h"
#include "geolib_math.h"
#include <cmath>

namespace cseis_geolib {

 struct csTime_t {
   public:
   csTime_t() { s = 0; us = 0; }
   csTime_t( int second, int usecond ) { s = second; us = usecond; }
   csTime_t( csTime_t const& t ) { s = t.s; us = t.us; }
   void add( csTime_t const& t ) { add(t.s, t.us); }
   void subtract( csTime_t const& t ) { add(-t.s, -t.us); }
   void add( int t_s, int t_us ) { s += t_s; us += t_us; int add = SIGN(us)*(int)( abs(us)/1000000 ); s += add; us -= add*1000000; if(us<0) { s -= 1; us += 1000000; } }
//   void add( csTime_t const& t ) { s += t.s; us += t.us; if(us>0) { int add = -us += 1000000; s -= 1; } int add = (int)(t.us/1000000); s += add; us -= add*1000000; }
//   void add( int time_s, int time_us ) { s += time_s; us += time_us; int add = (int)(time_us/1000000); s += add; us -= add*1000000; }
   bool operator<( const csTime_t& obj ) const { return( s < obj.s || (s == obj.s && us < obj.us) ); }
   bool operator>( const csTime_t& obj ) const { return( obj.s < s || (obj.s == s && obj.us < us) ); }
   csTime_t& operator=( csTime_t const& obj ) { s = obj.s; us = obj.us; return *this; }
   int s;
   int us;
   double time_ms() const { return( (double)s*1000.0 + (double)us/1000.0 ); }
 };

 struct csDate_t {
   public:
   int year;
   int julianDay;
   int hour;
   int min;
   int sec;
   int usec;
   char text[18];
   csDate_t() { init(); }
   csDate_t( int unixSec ) {
     init(); 
     csGeolibUtils::UNIXsec2date( year, julianDay, hour, min, sec, unixSec );
   }
   int unixTime() {
     return csGeolibUtils::date2UNIXsec( year, julianDay, hour, min, sec );
   }
   void init() {
     text[17] = '\0';
     year = 1970;
     julianDay = 1;
     hour = 0;
     min  = 0;
     sec  = 0;
     usec = 0;
   }
    char const* getString() {
     std::sprintf(text,"%04d%03d%02d%02d%02d.%06d", year, julianDay, hour, min, sec, usec);
     return text;
   }
 };


} // end namespace
#endif


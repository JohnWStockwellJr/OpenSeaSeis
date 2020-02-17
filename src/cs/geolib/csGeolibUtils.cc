/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csGeolibUtils.h"
#include "geolib_defines.h"
#include <cmath>
#include <ctime>
#include <string>

using namespace cseis_geolib;

cseis_geolib::type_t csGeolibUtils::text2Type( std::string text ) {
  if( !text.compare("int") ) {
    return cseis_geolib::TYPE_INT;
  }
  else if( !text.compare("float") ) {
    return cseis_geolib::TYPE_FLOAT;
  }
  else if( !text.compare("float_ibm") ) {
    return cseis_geolib::TYPE_FLOAT_IBM;
  }
  else if( !text.compare("double") ) {
    return cseis_geolib::TYPE_DOUBLE;
  }
  else if( !text.compare("int64") ) {
    return cseis_geolib::TYPE_INT64;
  }
  else if( !text.compare("short") ) {
    return cseis_geolib::TYPE_SHORT;
  }
  else if( !text.compare("ushort") ) {
    return cseis_geolib::TYPE_USHORT;
  }
  else if( !text.compare("char") ) {
    return cseis_geolib::TYPE_CHAR;
  }
  else if( !text.compare("string") ) {
    return cseis_geolib::TYPE_STRING;
  }
  else {
    return cseis_geolib::TYPE_UNKNOWN;
  }
}

csInt64_t csGeolibUtils::date2UNIXmsec( int year, int month, int day_of_month, int hour, int min, int sec, int msec ) {
  int julianDay = csGeolibUtils::date2JulianDay( year, month, day_of_month );
  return csGeolibUtils::date2UNIXmsec( year, julianDay, hour, min, sec, msec );
}
int csGeolibUtils::date2UNIXsec( int year, int month, int day_of_month, int hour, int min, int sec ) {
  int julianDay = csGeolibUtils::date2JulianDay( year, month, day_of_month );
  return csGeolibUtils::date2UNIXsec( year, julianDay, hour, min, sec );
}
int csGeolibUtils::date2JulianDay( int year, int month, int day ) {
  // Compute jdn1 = Julian day of 1-January-YEAR
  int month0 = 1;
  int day0   = 1;
  int a = int((14-month0)/12);
  int y = year + 4800 - a;
  int m = month0 + 12*a - 3;
  int jdn1 = day0 + int((153*m+2)/5) + 365*y + int(y/4) - int(y/100) + int(y/400); // - 32045;

  // Compute jdn2 = Julian day of DAY-MON-YEAR
  a = int((14-month)/12);
  y = year + 4800 - a;
  m = month + 12*a - 3;
  int jdn2 = day + int((153*m+2)/5) + 365*y + int(y/4) - int(y/100) + int(y/400); // - 32045;

  return( jdn2-jdn1+1 );
}

int csGeolibUtils::date2UNIXsec( int year, int day, int hour, int min, int sec ) {
  int shortyear;
  int a4, a100, a400;
  int intervening_leap_days;
  int days;
  //  double matlabTo1970sec = (double)MATLAB_SEC_01JAN1970;
  shortyear = year - 1900;

  a4 = (shortyear >> 2) + 475 - ! (shortyear & 3);
  a100 = a4 / 25 - (a4 % 25 < 0);
  a400 = a100 >> 2;
  intervening_leap_days = (a4 - 492) - (a100 - 19) + (a400 - 4);

  days = (365 * (shortyear - 70) + intervening_leap_days + (day - 1));

  return( 60 * (60 * (24 * days + hour) + min) + sec );
}

csInt64_t csGeolibUtils::date2UNIXmsec( int year, int day, int hour, int min, int sec, int msec ) {
  int unixSec = csGeolibUtils::date2UNIXsec( year, day, hour, min, sec );
  return csInt64_t( 1000LL * (csInt64_t)(unixSec) + (csInt64_t)msec );
}

void csGeolibUtils::UNIXsec2date( int& year, int& julianDay, int& hour, int& min, int& sec, int unixSec ) {
  time_t tsec = (time_t)unixSec;
  struct tm* tm;
  tm   = gmtime( &tsec );
  year = int(tm->tm_year+1900);
  julianDay  = tm->tm_yday + 1;  // tm_yday starts at 0
  hour = tm->tm_hour;
  min  = tm->tm_min;
  sec  = tm->tm_sec;
}

void csGeolibUtils::UNIXsec2date( int& year, int& month, int& day, int& hour, int& min, int& sec, int unixSec ) {
  time_t tsec = (time_t)unixSec;
  struct tm* tm;
  tm   = gmtime( &tsec );
  year = int(tm->tm_year+1900);
  month = tm->tm_mon + 1;  // tm month starts at 0
  day  = tm->tm_mday;
  hour = tm->tm_hour;
  min  = tm->tm_min;
  sec  = tm->tm_sec;
}
void csGeolibUtils::julianDay2date( int& month, int& day, int year, int julianDay ) {
  month = 0;
  day   = 0;
  fprintf(stderr,"julianDay2date: NOT IMPLEMENTED YET\n");
}

void csGeolibUtils::printDate( FILE* stream, int unixSec ) {
  int year;
  int month;
  int day;
  int hour;
  int min;
  int sec;
  csGeolibUtils::UNIXsec2date( year, month, day, hour, min, sec, unixSec );

  fprintf(stream,"%4d-%02d-%02d %02d:%02d:%02d", year,month,day, hour,min,sec);
}
std::string csGeolibUtils::UNIXsec2dateString( int unixSec ) {
  int year;
  int month;
  int day;
  int hour;
  int min;
  int sec;
  csGeolibUtils::UNIXsec2date( year, month, day, hour, min, sec, unixSec );

  char text[20];
  sprintf(text,"%4d-%02d-%02d %02d:%02d:%02d", year,month,day, hour,min,sec);
  text[19] = '\0';
  return text;
}


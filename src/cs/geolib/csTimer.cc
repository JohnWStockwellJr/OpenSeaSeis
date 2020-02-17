/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csTimer.h"
#include <cstdio>

using namespace cseis_geolib;

csTimer::csTimer() {
}

#ifdef PLATFORM_WINDOWS
void csTimer::start() {
  ftime( &tStart );
}
 
double csTimer::getElapsedTime() {
  ftime( &tEnd );
  return timeb_subtract( &tStart, &tEnd );
}

double csTimer::timeb_subtract( struct timeb* t1, struct timeb* t2 ) {
  double time_ms  = (double)((t2->time - t1->time) * 1000);
  int msec_extra;
  if( t2->time > t1->time ) {
    msec_extra = t2->millitm - t1->millitm;
  }
  else {
    msec_extra = t2->millitm - t1->millitm;
  }
  time_ms += (double)msec_extra;

  return( time_ms/1000.0 );
}

#endif

#ifdef PLATFORM_SGI
void csTimer::start() {
  gettimeofday( &tStart );
}

double csTimer::getElapsedTime() {
  gettimeofday( &tEnd );
  return timeval_subtract( &tStart, &tEnd );
}

/// Subtract two timeval objects. Return time difference in seconds
double csTimer::timeval_subtract( struct timeval* t1, struct timeval* t2 ) {
  if( t2->tv_usec < t1->tv_usec ) {
    int nsec = (t1->tv_usec - t2->tv_usec) / 1000000 + 1;
    t1->tv_usec -= 1000000 * nsec;
    t1->tv_sec += nsec;
  }
  if( t2->tv_usec - t1->tv_usec > 1000000 ) {
    int nsec = (t2->tv_usec - t1->tv_usec) / 1000000;
    t1->tv_usec += 1000000 * nsec;
    t1->tv_sec -= nsec;
  }
  
  /* Compute the time. tv_usec is certainly positive. */
  double time_us  = (double)((t2->tv_sec - t1->tv_sec) * 1000000);
  time_us += (double)(t2->tv_usec - t1->tv_usec);
  return( time_us/1000000.0 );
}
#endif

#if defined(PLATFORM_LINUX) || defined(PLATFORM_APPLE)
void csTimer::start() {
  ftime( &tStart );
}

double csTimer::getElapsedTime() {
  ftime( &tEnd );
  return timeb_subtract( &tStart, &tEnd );
}

double csTimer::timeb_subtract( struct timeb* t1, struct timeb* t2 ) {
  double time_ms  = (double)((t2->time - t1->time) * 1000);
  int msec_extra;
  if( t2->time > t1->time ) {
    msec_extra = t2->millitm - t1->millitm;
  }
  else {
    msec_extra = t2->millitm - t1->millitm;
  }
  time_ms += (double)msec_extra;

  return( time_ms/1000.0 );
}

#endif



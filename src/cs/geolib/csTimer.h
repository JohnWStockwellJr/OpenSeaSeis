/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TIMER_H
#define CS_TIMER_H

#include "geolib_platform_dependent.h"
#include <ctime>
#include <sys/time.h>

//#ifdef PLATFORM_WINDOWS
#include <sys/timeb.h>
//#endif
//#ifdef PLATFORM_LINUX
//#include <sys/timeb.h>
//#endif
//#ifdef PLATFORM_SGI
//#include <sys/times.h>
//#endif

namespace cseis_geolib {

/**
* Simple Timer
* Computes elapsed CPU time, deals with all the nasty platform dependencies
* How to use:
*   Call start() to start timer
*   Call getElapsedTime() to get snapshot of elapsed time
*/
class csTimer {
public:
  csTimer();
  ~csTimer() {}
  /// Start timer (reset timer to 0)
  void start();
  /// @return elapsed time in seconds
  double getElapsedTime();
#ifdef PLATFORM_WINDOWS
/// Subtract two timeval objects. Return time difference in seconds
  double timeb_subtract( struct timeb* tStart, struct timeb* tEnd );
private:
  struct timeb tStart;
  struct timeb tEnd;
#endif
#if defined(PLATFORM_LINUX) || defined(PLATFORM_APPLE)
/// Subtract two timeval objects. Return time difference in seconds
  double timeb_subtract( struct timeb* tStart, struct timeb* tEnd );
private:
  struct timeb tStart;
  struct timeb tEnd;
#endif
#ifdef PLATFORM_SGI
/// Subtract two timeval objects. Return time difference in seconds
  double timeval_subtract( struct timeval* tStart, struct timeval* tEnd );
private:
  struct timeval tStart;
  struct timeval tEnd;
#endif
};

} // namespace
#endif


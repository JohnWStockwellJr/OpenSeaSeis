/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_ABSOLUTE_TIME_H
#define CS_ABSOLUTE_TIME_H

#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {

  /**
   * Absolute time class
   *
   * Contains absolute time
   *
   * @author Bjorn Olofsson
   * @date 2008
   */
class csAbsoluteTime {
 public:
  static int const BYTE_SIZE = 12;

 public:
  csAbsoluteTime();
  csAbsoluteTime( byte_t const* buffer);
  csAbsoluteTime( int days, int sec, int us );
  ~csAbsoluteTime();

  void set( csAbsoluteTime const* obj );
  void set( byte_t const* buffer );
  void set( int days, int sec, int us );
  void get( int& days, int& sec, int& us ) const;
  std::string getDate() const;

  int writeBytes( byte_t* buffer ) const;

  /// @return packed String object containing abolute time for disk storage purposes
  //std::string serialize();

 private:
  csAbsoluteTime( csAbsoluteTime const& obj );

  int myTime_days; // Full days since 01-Jan-1970
  int myTime_sec;  // Fraction: Full seconds
  int myTime_us;   // Fraction: Microseconds

};


} // end namespace


#endif


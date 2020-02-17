/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_POINT_3D_H
#define CS_POINT_3D_H

#include <cstdio>
#include <string>

namespace cseis_geolib {

class csPoint3D {
public:
  csPoint3D();
  csPoint3D( double x_in, double y_in, double z_in );
  csPoint3D( csPoint3D const& p );
  ~csPoint3D();

  void set( double x_in, double y_in, double z_in );
  double dotProduct( csPoint3D const& p ) const;
  csPoint3D crossProduct( csPoint3D const& p ) const;
  /**
   * @return Length of 3D vector (assuming 3D vector)
   */
  double vectorLength() const;
  /**
   * @param vec: Second vector
   * @return Angle between two 3D vectors
   */
  double vectorAngle( csPoint3D const& vec ) const;
  /**
   * @param vec: Second vector
   * @return Angle between two 3D vectors
   */
  double vectorAngle_rad( csPoint3D const& vec ) const;
  /**
   * @param point: Second point
   * @return Distance between two 3D points
   */
  double distance( csPoint3D const& point ) const;
  void rotateAzim( double azim_rad, csPoint3D const& pOrigin );

  void dump( std::FILE* stream ) const;

  csPoint3D operator+( csPoint3D const& p ) const;
  csPoint3D operator-( csPoint3D const& p ) const;
  csPoint3D operator*( double scalar ) const;
  csPoint3D operator/( double value ) const;
  csPoint3D& operator*=( double scalar );
  csPoint3D& operator/=( double value );
  csPoint3D& operator+=( double value );
  csPoint3D& operator-=( double value );
  csPoint3D& operator=( csPoint3D const & p );

  std::string toString() const;
 public:
  double x;
  double y;
  double z;
};

} // END namespace
#endif

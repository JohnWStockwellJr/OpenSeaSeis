/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_POINT_2D_H
#define CS_POINT_2D_H

#include <vector>
#include <cstdio>

namespace cseis_geolib {

class csPoint2D {
public:
  csPoint2D();
  csPoint2D( double x_in, double z_in );
  csPoint2D( csPoint2D const& p );
  void rotateAxis( double angleRot );
  bool isPointInPolygon( cseis_geolib::csPoint2D const* polygon, int numPoints ) const;
  bool isPointInPolygon( std::vector<cseis_geolib::csPoint2D> const& polygon ) const;
  bool isPointInPolygon( std::vector<cseis_geolib::csPoint2D> const* polygon ) const;
  double dotProduct( csPoint2D const& p ) const;
  double vectorLength() const;
  double vectorAngle( csPoint2D const& vec ) const;
  double distance( csPoint2D const& p ) const;

  csPoint2D operator+( csPoint2D const& p ) const;
  csPoint2D operator-( csPoint2D const& p ) const;

  void dump( FILE* stream ) const;

 public:
  double x;
  double z;
};

} // END namespace
#endif

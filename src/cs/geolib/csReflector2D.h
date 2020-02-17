/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_REFLECTOR_2D_H
#define CS_REFLECTOR_2D_H

#include "csPoint2D.h"
#include "csPoint3D.h"

namespace cseis_geolib {

class csReflector2D {
public:
  csReflector2D( double x0, double z0, double angle_deg );
  csReflector2D( cseis_geolib::csPoint2D pOrig, double angle_deg );
  csReflector2D( cseis_geolib::csPoint2D orig, cseis_geolib::csPoint2D direction );
  double compute( double x ) const;
  void set( cseis_geolib::csPoint2D orig, cseis_geolib::csPoint2D direction );
  cseis_geolib::csPoint2D direction() const;
  cseis_geolib::csPoint2D intersection( csReflector2D const& line ) const;

private:
  void compute_internal();
  cseis_geolib::csPoint2D myOrig;
  cseis_geolib::csPoint2D myDirection;
  double myAngle_deg;
  double myTanAngle;
};

} // END namespace
#endif

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csReflector2D.h"
#include "csPoint2D.h"
#include "csPlane3D.h"
#include <cmath>

using namespace cseis_geolib;

csReflector2D::csReflector2D( double x0, double z0, double angle_deg ) {
  myOrig.x = x0;
  myOrig.z = z0;
  myAngle_deg = angle_deg;
  compute_internal();
}
csReflector2D::csReflector2D( cseis_geolib::csPoint2D pOrig, double angle_deg ) {
  myOrig.x = pOrig.x;
  myOrig.z = pOrig.z;
  myAngle_deg = angle_deg;
  compute_internal();
}
csReflector2D::csReflector2D( cseis_geolib::csPoint2D orig, cseis_geolib::csPoint2D direction ) {
  myOrig = orig;
  myDirection = direction;
  myAngle_deg = atan2( myDirection.z, myDirection.x ) * 180.0 / M_PI;
  myTanAngle = tan(myAngle_deg * M_PI / 180.0);
}
double csReflector2D::compute( double x ) const {
  return( myOrig.z + (x - myOrig.x) * myTanAngle );
}
void csReflector2D::set( cseis_geolib::csPoint2D orig, cseis_geolib::csPoint2D direction ) {
  myOrig = orig;
  myDirection = direction;
}
cseis_geolib::csPoint2D csReflector2D::direction() const {
  return myDirection;
}
csPoint2D csReflector2D::intersection( csReflector2D const& line ) const {
  csPoint2D point( 0, 0 );
  double upper = ( line.myOrig.z - myOrig.z ) * line.myDirection.x + ( myOrig.x - line.myOrig.x ) * line.myDirection.z;
  double lower = line.myDirection.x * myDirection.z - line.myDirection.z * myDirection.x;
  if( fabs(lower) < 1.0e-12 ) return point;
  double lamda = upper / lower;
  point.x = myOrig.x + lamda * myDirection.x;
  point.z = myOrig.z + lamda * myDirection.z;
  return point;
}
void csReflector2D::compute_internal() {
  myTanAngle = tan(myAngle_deg * M_PI / 180.0);
  double z2 = compute( myOrig.x+100 );
  myDirection.x = 100;
  myDirection.z = z2 - myOrig.z;
}

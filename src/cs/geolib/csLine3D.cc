/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include <cmath>
#include "csPoint3D.h"
#include "csLine3D.h"

using namespace cseis_geolib;
using namespace std;

csLine3D::csLine3D() : origin(), direction(), point2(), myOffset(0), myLength(0), myAzimuth(0) {
}
csLine3D::csLine3D( cseis_geolib::csPoint3D const& originPoint, cseis_geolib::csPoint3D const& directionVector ) {
  origin    = originPoint;
  direction = directionVector;
  point2    = origin + direction;
  compute_internal();
}
csLine3D::~csLine3D() {
}
void csLine3D::setFromPoints( cseis_geolib::csPoint3D const& p1, cseis_geolib::csPoint3D const& p2 ) {
  origin  = p1;
  point2  = p2;
  direction = point2 - origin;
  compute_internal();
}
cseis_geolib::csPoint3D csLine3D::pointAt( double lamda ) const {
  return cseis_geolib::csPoint3D( origin.x + lamda*direction.x, origin.y + lamda*direction.y, origin.z + lamda*direction.z );
}
csPoint3D csLine3D::computePointFromX( double x ) const {
  csPoint3D p( x, 0, 0 );
  double lamda = abs(direction.x) > 1.0e-12 ? ( x - origin.x ) / direction.x : 0.0;
  p.y = origin.y + lamda * direction.y;
  p.z = origin.z + lamda * direction.z;
  return p;
}
csPoint3D csLine3D::computePointFromY( double y ) const {
  csPoint3D p( 0, y, 0 );
  double lamda = abs(direction.y) > 1.0e-12 ? ( y - origin.y ) / direction.y : 0.0;
  p.x = origin.x + lamda * direction.x;
  p.z = origin.z + lamda * direction.z;
  return p;
}
csPoint3D csLine3D::computePointFromZ( double z ) const {
  csPoint3D p( 0, 0, z );
  double lamda = abs(direction.z) > 1.0e-12 ? ( z - origin.z ) / direction.z : 0.0;
  p.x = origin.x + lamda * direction.x;
  p.y = origin.y + lamda * direction.y;
  return p;
}
double csLine3D::distanceToPoint( csPoint3D const& point ) const {
  double lenDirVec = direction.vectorLength();
  if( lenDirVec < 1.0e-5 ) return 0.0;
  csPoint3D tmp1 = origin - point;
  csPoint3D tmp2 = direction.crossProduct( tmp1 );
  double distance = tmp2.vectorLength() / lenDirVec;
  return distance;
}
double csLine3D::distanceLineToPoint( csPoint3D const& linePoint1, csPoint3D const& linePoint2, csPoint3D const& point ) {
  csPoint3D lineDirection = linePoint2 - linePoint1;
  double lenDirVec = lineDirection.vectorLength();
  if( lenDirVec < 1.0e-5 ) return 0.0;
  csPoint3D tmp1 = linePoint1 - point;
  csPoint3D tmp2 = lineDirection.crossProduct( tmp1 );
  double distance = tmp2.vectorLength() / lenDirVec;
  return distance;
}
void csLine3D::compute_internal() {
  myAzimuth = atan2( direction.x, direction.y );
  myOffset  = sqrt( direction.x*direction.x + direction.y*direction.y );
  myLength  = sqrt( myOffset*myOffset + direction.z*direction.z );
}
/*csPoint3D csLine3D::intersection( csLine3D const& line ) const {
  csPoint3D point( 0, 0, 0 );
  csPoint3D pOrigDiff = line.origin - origin;
  double upper = pOrigDiff.dotProduct( line.direction );
  double lower = line.direction.x * direction.z - line.direction.z * direction.x;
  if( abs(lower) < 1.0e-12 ) return
}
*/

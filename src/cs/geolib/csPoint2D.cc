/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csPoint2D.h"
#include <cmath>

using namespace cseis_geolib;
using namespace std;

csPoint2D::csPoint2D() {
  x = 0;
  z = 0;
}
csPoint2D::csPoint2D( double x_in, double z_in ) {
  x = x_in;
  z = z_in;
}
csPoint2D::csPoint2D( csPoint2D const& p ) {
  x = p.x;
  z = p.z;
}
void csPoint2D::rotateAxis( double angleRot ) {
  double angleRot_an = angleRot * M_PI / 180.0;
  double sina = sin(angleRot_an);
  double cosa = cos(angleRot_an);
  double x_rot =  x*cosa + z*sina;
  double z_rot = -x*sina + z*cosa;
  x = x_rot;
  z = z_rot;
} 
bool csPoint2D::isPointInPolygon( csPoint2D const* polygon, int numPoints ) const {
  bool test = false;
  for( int i = 0, j = numPoints-1; i < numPoints; j = i++ ) {
    if( ((polygon[i].z > z) != ( polygon[j].z > z)) &&
        (x < (polygon[j].x-polygon[i].x) * (z-polygon[i].z) / (polygon[j].z - polygon[i].z) + polygon[i].x) )
      {
        test = !test;
      }
  }
  return test;
}
bool csPoint2D::isPointInPolygon( std::vector<csPoint2D> const& polygon ) const {
  return isPointInPolygon( &polygon );
}
bool csPoint2D::isPointInPolygon( std::vector<csPoint2D> const* polygon ) const {
  bool test = false;
  int numPoints = polygon->size();
  for( int i = 0, j = numPoints-1; i < numPoints; j = i++ ) {
    if( ((polygon->at(i).z > z) != ( polygon->at(j).z > z)) &&
        (x < (polygon->at(j).x-polygon->at(i).x) * (z-polygon->at(i).z) / (polygon->at(j).z - polygon->at(i).z) + polygon->at(i).x) )
      {
        test = !test;
      }
  }
  return test;
}
double csPoint2D::dotProduct( csPoint2D const& p ) const {
  return( p.x * x  +  p.z * z );
}
double csPoint2D::vectorLength() const {
  return( sqrt( x*x + z*z ) );
}
double csPoint2D::vectorAngle( csPoint2D const& vec ) const {
  double upper = dotProduct( vec );
  double lower = vectorLength() * vec.vectorLength();
  if( abs(lower) > 1.0e-9 ) {
    return( acos( upper / lower ) * 180.0 / M_PI );
  }
  else {
    return 0.0;
  }
}
void csPoint2D::dump( FILE* stream ) const {
  fprintf(stream,"%.4f %.4f\n", x, z );
}
csPoint2D csPoint2D::operator+( csPoint2D const& p ) const {
  return( csPoint2D( x+p.x, z+p.z ) );
}
csPoint2D csPoint2D::operator-( csPoint2D const& p ) const {
  return( csPoint2D( x-p.x, z-p.z ) );
}
double csPoint2D::distance( csPoint2D const& p ) const {
  return( sqrt( pow(x-p.x,2) + pow(z-p.z,2) ) );
}

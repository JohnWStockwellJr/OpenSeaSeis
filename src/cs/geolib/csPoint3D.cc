/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csPoint3D.h"
#include <cmath>
#include <cstdio>

using namespace cseis_geolib;
using namespace std;

csPoint3D::csPoint3D() {
  x = 0;
  y = 0;
  z = 0;
}
csPoint3D::csPoint3D( double x_in, double y_in, double z_in ) {
  x = x_in;
  y = y_in;
  z = z_in;
}
csPoint3D::csPoint3D( csPoint3D const& p ) {
  x = p.x;
  y = p.y;
  z = p.z;
}
csPoint3D::~csPoint3D() {
} 

void csPoint3D::set( double x_in, double y_in, double z_in ) {
  x = x_in;
  y = y_in;
  z = z_in;
}
csPoint3D csPoint3D::operator+( csPoint3D const& p ) const {
  return( csPoint3D( x+p.x, y+p.y, z+p.z ) );
}
csPoint3D csPoint3D::operator-( csPoint3D const& p ) const {
  return( csPoint3D( x-p.x, y-p.y, z-p.z ) );
}
csPoint3D csPoint3D::operator*( double scalar ) const {
  return( csPoint3D( x*scalar, y*scalar, z*scalar ) );
}
csPoint3D csPoint3D::operator/( double value ) const {
  return( csPoint3D( x/value, y/value, z/value ) );
}
csPoint3D& csPoint3D::operator*=( double scalar ) {
  x *= scalar;
  y *= scalar;
  z *= scalar;
  return( *this );
}
csPoint3D& csPoint3D::operator/=( double value ) {
  x /= value;
  y /= value;
  z /= value;
  return( *this );
}
csPoint3D& csPoint3D::operator+=( double value ) {
  x += value;
  y += value;
  z += value;
  return( *this );
}
csPoint3D& csPoint3D::operator-=( double value ) {
  x -= value;
  y -= value;
  z -= value;
  return( *this );
}
csPoint3D& csPoint3D::operator=( csPoint3D const & p ) {
  x = p.x;
  y = p.y;
  z = p.z;
  return( *this );
}
double csPoint3D::distance( csPoint3D const& point ) const {
  return( sqrt( pow(x-point.x,2) + pow(y-point.y,2) + pow(z-point.z,2) ) );
}
// Rotate around origin
void csPoint3D::rotateAzim( double azim_rad, csPoint3D const& pOrigin ) {
  double dx = x - pOrigin.x;
  double dy = y - pOrigin.y;
  x = dx * cos( azim_rad ) - dy * sin( azim_rad ) + pOrigin.x;
  y = dy * sin( azim_rad ) + dy * cos( azim_rad ) + pOrigin.y;
}
double csPoint3D::dotProduct( csPoint3D const& p ) const {
  return( p.x * x  +  p.y * y  +  p.z * z );
}
csPoint3D csPoint3D::crossProduct( csPoint3D const& p ) const {
  // 'this' x 'p'
  csPoint3D result( y*p.z - z*p.y, z*p.x - x*p.z, x*p.y - y*p.x );
  return( result );
}
double csPoint3D::vectorLength() const {
  return( sqrt( x*x + y*y + z*z ) );
}
double csPoint3D::vectorAngle( csPoint3D const& vec ) const {
  double upper = dotProduct( vec );
  double lower = vectorLength() * vec.vectorLength();
  if( abs(lower) > 1.0e-9 ) {
    return( acos( upper / lower ) * 180.0 / M_PI );
  }
  else {
    return 0.0;
  }
}
double csPoint3D::vectorAngle_rad( csPoint3D const& vec ) const {
  double upper = dotProduct( vec );
  double lower = vectorLength() * vec.vectorLength();
  if( abs(lower) > 1.0e-9 ) {
    return( acos( upper / lower ) );
  }
  else {
    return 0.0;
  }
}
void csPoint3D::dump( FILE* stream ) const {
  fprintf(stream,"%.4f %.4f %.4f\n", x, y, z );
}
std::string csPoint3D::toString() const {
  char text[40];
  sprintf(text,"%10.6e %10.6e %10.6e", x, y, z );
  return std::string( text );
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

//********************************************************************************
//
// Solves for isotropic CCP position using Tessmer & Behle's (1988) formula
//
//********************************************************************************
//
//
// (i) vpvs : vp/vs ratio
// (i) offset: src-rcv offset [m]
// (i) depth : Reflector depth [m]
// (o) offset_pside: src-ccp offset [m]
#include <cmath>

namespace cseis_geolib {

void ccp_offset_iso( double vpvs, double offset, double depth, double* offset_pside ) {

  double k, h, p, q, d, y_tilde, y;
  double a_right, b_right;
  double ab_left;
  double m1, m2, m3, m4;

  k = ( vpvs*vpvs + 1.0 ) / ( vpvs*vpvs - 1.0 );
  h = ( 0.5 * (offset/depth) );
  h = h * h;
  p = (1.0/9.0) * ( -4.0*h*h - 2.0*h - 0.25  );
  q = (1.0/54.0) * ( -16*h*h*h - 12.0*h*h - 3.0*h + 13.5*(1.0 - k*k)*h - 0.25 );
  d = q*q + p*p*p;
  y_tilde = pow(( -q + sqrt(d) ),(1.0/3.0)) + pow(( -q - sqrt(d) ),(1.0/3.0));
  y = y_tilde - (h/3.0) + (1.0/6.0);

  a_right = 0.5*(h-y-0.5) - ( k*sqrt(h) )/( sqrt(2.0*(y+h) - 1.0) );
  b_right = 0.5*(h-y-0.5) + ( k*sqrt(h) )/( sqrt(2.0*(y+h) - 1.0) );
  ab_left = 0.5*sqrt( 2.0*(y+h)-1.0 );

  m1 = 0.0;
  m2 = 0.0;
  m3 = 0.0;
  m4 = 0.0;
  *offset_pside = -999.9;
  if( a_right >= 0.0 ) {
    m1 = -ab_left + sqrt(a_right);
    m2 = -ab_left - sqrt(a_right);
    if( m1*depth <= 0.5*offset ) {
      *offset_pside = 0.5*offset + m1*depth;
    }
    else if( m2*depth <= 0.5*offset ) {
      *offset_pside = 0.5*offset + m2*depth;
    }
  }
  if( b_right >= 0.0 ) {
    m3 =  ab_left + sqrt(b_right);
    m4 =  ab_left - sqrt(b_right);
    if( m3*depth <= 0.5*offset ) {
      *offset_pside = 0.5*offset + m3*depth;
    }
    else if( m4*depth <= 0.5*offset ) {
      *offset_pside = 0.5*offset + m4*depth;
    }
  }
}

} // namespace


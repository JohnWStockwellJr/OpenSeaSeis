/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <cstdio>
#include "csInterpolation2D.h"
#include "csPoint2D.h"

using namespace std;
using namespace cseis_geolib;

csInterpolation2D::csInterpolation2D( int method ) {
  myFactor1 = 0;
  myFactor2 = 0;
  myFactor3 = 0;
  myFactor4 = 0;
  myMethod = method;
// myAlpha = 0;
// myBeta = 0;
}
csInterpolation2D::~csInterpolation2D() {
}

bool csInterpolation2D::prepareBilinear( cseis_geolib::csPoint2D const& p1, cseis_geolib::csPoint2D const& p2, cseis_geolib::csPoint2D const& pNew ) {
  return( prepare( p1, p2, p1, p2, pNew ) );
}

bool csInterpolation2D::prepare( cseis_geolib::csPoint2D const& p1,
                                 cseis_geolib::csPoint2D const& p2,
                                 cseis_geolib::csPoint2D const& p3,
                                 cseis_geolib::csPoint2D const& p4,
                                 cseis_geolib::csPoint2D const& pNew ) {
  if( myMethod == METHOD_BILINEAR ) {
    // Assume that p1 and p2 define the far corners of a rectangular area
    double tmp0 = ( p2.z - p1.z )*( p2.x - p1.z );
    if( abs(tmp0) > 1.0e-14 ) {
      tmp0 = 1.0/tmp0;
    }
    else {
      tmp0 = 1.0;
    }
    double tmp_row1 = p2.x - pNew.x;
    double tmp_row2 = pNew.x - p1.x;
    double tmp_col1 = p2.z - pNew.z;
    double tmp_col2 = pNew.z - p1.z;
    myFactor1 = tmp_col1 * tmp_row1 * tmp0;
    myFactor2 = tmp_col2 * tmp_row1 * tmp0;
    myFactor3 = tmp_col1 * tmp_row2 * tmp0;
    myFactor4 = tmp_col2 * tmp_row2 * tmp0;
    return true;
  }
  //
  // Points must be sorted into zig-zag line !!!
  // 
  double a = -p1.x + p3.x;
  double b = -p1.x + p2.x;
  double c = p1.x - p2.x - p3.x + p4.x;
  double d = pNew.x - p1.x;
  double e = -p1.z + p3.z;
  double f = -p1.z + p2.z;
  double g = p1.z - p2.z - p3.z + p4.z;
  double h = pNew.z - p1.z;

  double be = b*e;
  double af = a*f;
  double dg = d*g;
  double ch = c*h;

  double ce = c*e;
  double ag = a*g;

  double temp1 = (2*ce - 2*ag);
  double temp2 = (2*c*f - 2*b*g);
  if( temp1 == 0 || temp2 == 0 ) {
    myFactor1 = 1.0;
    myFactor2 = 0.0;
    myFactor3 = 0.0;
    myFactor4 = 0.0;
    // Linear interpolation;
    /*    double dist1 = sqrt( pow(p1.x-pNew.x, 2.0) + pow(p1.z-pNew.z,2.0) );
    double dist2 = sqrt( pow(p2.x-pNew.x, 2.0) + pow(p2.z-pNew.z,2.0) );
    double totalDist = (dist1 + dist2);
    myFactor1 = dist2/totalDist;
    myFactor2 = dist1/totalDist;
    myFactor3 = 0; //dist2/totalDist;
    myFactor4 = 0; //dist1/totalDist; */
    return false;
  }

  double temp1_sqrt = -4 * (ce - ag) * (d*f - b*h) + pow((be - af + dg - ch),2);
  if( temp1_sqrt < 0 ) temp1_sqrt = 0;

  // Solution #1
  double alpha1 = -(be - af + dg - ch + sqrt( temp1_sqrt ))/temp1;
  double beta1  =  (be - af - dg + ch + sqrt( temp1_sqrt ))/temp2;

  // Solution #2
  double alpha2 =   (-be + af - dg + ch + sqrt( temp1_sqrt )) / temp1;
  double beta2  = -((-be + af + dg - ch + sqrt( temp1_sqrt )) / temp2);

  double alphaFinal;
  double betaFinal;
  if( alpha1 >= 0 && alpha1 <= 1.0 && beta1 >= 0 && beta1 <= 1.0 ) {
    alphaFinal = alpha1;
    betaFinal  = beta1;
  }
  else if( alpha2 >= 0 && alpha2 <= 1.0 && beta2 >= 0 && beta2 <= 1.0 ) {
    alphaFinal = alpha2;
    betaFinal  = beta2;
  }
  else {
    double t1 = pow( alpha1*alpha1 + beta1*beta1, 2.0 );
    double t2 = pow( alpha2*alpha2 + beta2*beta2, 2.0 );
    if( t1 < t2 ) {
      alphaFinal = alpha1;
      betaFinal  = beta1;
    }
    else {
      alphaFinal = alpha2;
      betaFinal  = beta2;
    }
  }
  myFactor1 = (1 - alphaFinal) * (1 - betaFinal);
  myFactor2 = (1 - alphaFinal) * betaFinal;
  myFactor3 = alphaFinal * (1 - betaFinal);
  myFactor4 = alphaFinal * betaFinal;
  //  fprintf(stderr,"alpha/beta: %f %f   (%f %f  %f %f)\n", alphaFinal, betaFinal, alpha1, beta1, alpha2, beta2);
  return true;
}

double csInterpolation2D::compute( double val1, double val2, double val3, double val4 ) {
  return( myFactor1 * val1 + myFactor2 * val2 + myFactor3 * val3 + myFactor4 * val4 );
}
//double csInterpolation2D::compute( double val1, double val2, double val3, double val4 ) {
//  return( (1 - myAlpha) * ((1 - myBeta) * val1 + myBeta * val2) + myAlpha * ((1 - myBeta) * val3 + myBeta * val4) );
//}

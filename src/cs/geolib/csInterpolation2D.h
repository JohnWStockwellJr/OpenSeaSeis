/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_INTERPOLATION_2D_H
#define CS_INTERPOLATION_2D_H

namespace cseis_geolib {

class csPoint2D;
/**
 * Provides interpolation methods
 */
class csInterpolation2D {
 public:
  static const int METHOD_BILINEAR = 40;
  static const int METHOD_GENERAL  = 41;

 public:
  csInterpolation2D( int method );
  ~csInterpolation2D();
  bool prepare( cseis_geolib::csPoint2D const& p1,
                cseis_geolib::csPoint2D const& p2,
                cseis_geolib::csPoint2D const& p3,
                cseis_geolib::csPoint2D const& p4,
                cseis_geolib::csPoint2D const& pNew );
  bool prepareBilinear( cseis_geolib::csPoint2D const& p1, cseis_geolib::csPoint2D const& p2, cseis_geolib::csPoint2D const& pNew );
  double compute( double val1, double val2, double val3, double val4 );

 private:
  csInterpolation2D();
  int myMethod;
  double myFactor1;
  double myFactor2;
  double myFactor3;
  double myFactor4;
  //  double myAlpha;
  //  double myBeta;
};

} // END namespace

#endif


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_PROFILE_H
#define CS_PROFILE_H

#include <cstdlib>
#include "csPoint3D.h"
#include "csLine3D.h"

namespace cseis_geolib {

class csProfile : csLine3D {
 public:
  static int const EXTRAPOLATION_CONSTANT = 1;
  static int const EXTRAPOLATION_LINEAR   = 2;

 public:
  csProfile();
  csProfile( cseis_geolib::csPoint3D const& p1, cseis_geolib::csPoint3D const& p2 );
  csProfile( cseis_geolib::csPoint3D const& p1, cseis_geolib::csPoint3D const& p2, int extrapolationMethod );
  ~csProfile();

  virtual cseis_geolib::csPoint3D pointAt( double lamda ) const;
  double depthAt( double lamda ) const;
  double offsetAt( double lamda ) const;

  double step() const { return myStep; }
  int numSteps() const { return myNumSteps; }
  double offset() const { return myOffset; }

  int initDepths( double step );
  void initDepths( int numSteps, double const* depths = NULL );
  void setDepth( int index, double depth );

 protected:
  void init( int extrapolationMethod);
  double* myDepths;
  double myStep;
  int myNumSteps;
  int myExtrapolationMethod;
};


} // end namespace

#endif


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "geolib_math.h"
#include <cmath>

using namespace cseis_geolib;

/**
* Rotate XYZ vectors using ICV matrix. This only affects roll & tilt
*/
void rotate_xyz_icv( float* xTrace,
         float* yTrace,
         float* zTrace,
         float* weight_i,
         float* weight_c,
         float* weight_v,
         int nSamples,
         bool apply )
{
  float xValue;
  float yValue;
  float zValue;
  
  if( !apply ) {
    // Apply transposed rotation matrix
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i]*weight_i[0] + yTrace[i]*weight_i[1] + zTrace[i]*weight_i[2];
      yValue = xTrace[i]*weight_c[0] + yTrace[i]*weight_c[1] + zTrace[i]*weight_c[2];
      zValue = xTrace[i]*weight_v[0] + yTrace[i]*weight_v[1] + zTrace[i]*weight_v[2];

      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
  else {  // if APPLY
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i]*weight_i[0] + yTrace[i]*weight_c[0] + zTrace[i]*weight_v[0];
      yValue = xTrace[i]*weight_i[1] + yTrace[i]*weight_c[1] + zTrace[i]*weight_v[1];
      zValue = xTrace[i]*weight_i[2] + yTrace[i]*weight_c[2] + zTrace[i]*weight_v[2];

      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
}


/***********************************************************************
 *
 *
 */

/**
 * Rotate XYZ vectors in horizontal plane by given azimuth
 * Mode 'apply':  Apply clockwise rotation to coordinate system, assuming 2D coordinate system such as X/East, Y/North
 * Mode 'remove': Apply ANTI-clockwise rotation to coordinate system
*/
void rotate_xy_azim( float* xTrace,
         float* yTrace,
         float azim,
         int nSamples,
         bool apply )
{
  float xValue;
  float yValue;

  if( !apply ) {
    // Apply transposed rotation matrix
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * cos(azim) + yTrace[i] * sin(azim);
      yValue = -xTrace[i] * sin(azim) + yTrace[i] * cos(azim);
      xTrace[i] = xValue;
      yTrace[i] = yValue;
    }
  }
  else {  // if APPLY
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * cos(azim) - yTrace[i] * sin(azim);
      yValue = xTrace[i] * sin(azim) + yTrace[i] * cos(azim);
      xTrace[i] = xValue;
      yTrace[i] = yValue;
    }
  }
}

/***********************************************************************
 *
 *
 */

/**
* Rotate XYZ vectors by given roll & tilt angle
*/
void rotate_xyz_roll_tilt( float* xTrace,
         float* yTrace,
         float* zTrace,
         float roll,
         float tilt,
         int nSamples,
         bool apply )
{
  float xValue;
  float yValue;
  float zValue;
  float i3, i4, i5, c3, c4, c5, v3, v4, v5;

  i3 = cos(tilt);
  i4 = sin(tilt)*cos(roll);
  i5 = sin(tilt)*sin(roll);

  c3 = 0.0;
  c4 = sin(roll);
  c5 = -cos(roll);

  v3 = -sin(tilt);
  v4 = cos(tilt)*cos(roll);
  v5 = cos(tilt)*sin(roll);

  if( !apply ) {
    // Apply transposed rotation matrix
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * i3 + yTrace[i] * i4 + zTrace[i] * i5;
      yValue = xTrace[i] * c3 + yTrace[i] * c4 + zTrace[i] * c5;
      zValue = xTrace[i] * v3 + yTrace[i] * v4 + zTrace[i] * v5;
      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
  else {  // if APPLY
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * i3 + yTrace[i] * c3 + zTrace[i] * v3;
      yValue = xTrace[i] * i4 + yTrace[i] * c4 + zTrace[i] * v4;
      zValue = xTrace[i] * i5 + yTrace[i] * c5 + zTrace[i] * v5;
      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
}


//******************************************************************************
/**
* Rotate XYZ vectors by given tilt X & tilt Y angle
*/
void rotate_xyz_tiltxy( float* xTrace,
                        float* yTrace,
                        float* zTrace,
                        float tiltx,
                        float tilty,
                        int nSamples,
                        bool apply )
{
  double xValue;
  double yValue;
  double zValue;
  double i3, i4, i5, c3, c4, c5, v3, v4, v5;

  tiltx = -tiltx;
  tilty = -tilty;

  double sin_tiltx = sin(tiltx);
  double cos_tiltx = cos(tiltx);
  if( std::abs(cos_tiltx) < 1e-30 ) {
    cos_tiltx = 1e-30;
  }
  double sinth = sin(tilty)/cos_tiltx;
  double tmp   = cos_tiltx * cos_tiltx - sin(tilty) * sin(tilty);
  if( tmp < 0 ) tmp = 0;   // Avoid potential sqrt(-1)
  double costh = sqrt( tmp ) / cos_tiltx;

  i3 = cos_tiltx;
  i4 = -sin_tiltx*sinth;
  i5 = -sin_tiltx*costh;

  c3 = 0.0;
  c4 = costh;
  c5 = -sinth;

  v3 = sin_tiltx;
  v4 = cos_tiltx*sinth;
  v5 = cos_tiltx*costh;

  if( !apply ) {
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * i3 + yTrace[i] * i4 + zTrace[i] * i5;
      yValue = xTrace[i] * c3 + yTrace[i] * c4 + zTrace[i] * c5;
      zValue = xTrace[i] * v3 + yTrace[i] * v4 + zTrace[i] * v5;
      xTrace[i] = (float)xValue;
      yTrace[i] = (float)yValue;
      zTrace[i] = (float)zValue;
    }
  }
  else {  // if APPLY
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * i3 + yTrace[i] * c3 + zTrace[i] * v3;
      yValue = xTrace[i] * i4 + yTrace[i] * c4 + zTrace[i] * v4;
      zValue = xTrace[i] * i5 + yTrace[i] * c5 + zTrace[i] * v5;
      xTrace[i] = (float)xValue;
      yTrace[i] = (float)yValue;
      zTrace[i] = (float)zValue;
    }
  }

}

//******************************************************************************
//
//

void rotate_galperin( float* xTrace, float* yTrace, float* zTrace, int nSamples, bool apply )
{
  // float a = 35.2644; // [deg]
  // float b = 45.0;    // [deg]
  static float const R11 = 0.81649648;  // cosa
  static float const R13 = 0.57735042;  // sina
  static float const R21 = -0.40824839; // -sinb*sina
  static float const R22 = 0.70710678;  // cosb
  static float const R23 = 0.57735020;  // sinb*cosa
  static float const R31 = -0.40824839; // -cosb*sina
  static float const R32 = -0.70710678; // -sinb
  static float const R33 = 0.57735020;  // cosb*cosa

  float xValue;
  float yValue;
  float zValue;

  if( !apply ) {
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * R11 + yTrace[i] * R21 + zTrace[i] * R31;
      yValue =                   yTrace[i] * R22 + zTrace[i] * R32;
      zValue = xTrace[i] * R13 + yTrace[i] * R23 + zTrace[i] * R33;
      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
  // Apply: rotate coordinate system from x(hor) y(hor) z(vert/down) to Galperin xyz
  // All components after rotation point downwards, new x component has same azimuth as original x component
  else {  // APPLY
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * R11                   + zTrace[i] * R13;
      yValue = xTrace[i] * R21 + yTrace[i] * R22 + zTrace[i] * R23;
      zValue = xTrace[i] * R31 + yTrace[i] * R32 + zTrace[i] * R33;
      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
}

// Method that showcases how to rotate from XYZ to Galperin configuration
void rotate_to_galperin( float rx, float ry, float rz, float& sx, float& sy, float& sz )
{
  float a = 35.2644;
  float b = 45.0;

  float a_rad = a * M_PI / 180.0;
  float b_rad = b * M_PI / 180.0;

  float sina = sin(a_rad);
  float cosa = cos(a_rad);
  float sinb = sin(b_rad);
  float cosb = cos(b_rad);

  sx =       cosa*rx +                sina*rz;
  sy = -sinb*sina*rx + cosb*ry + sinb*cosa*rz;
  sz = -cosb*sina*rx - sinb*ry + cosb*cosa*rz;
}

void rotate_general( float* xTrace, float* yTrace, float* zTrace, int nSamples, bool apply,
  float R11,
  float R12,
  float R13,
  float R21,
  float R22,
  float R23,
  float R31,
  float R32,
  float R33 )
{
  float xValue;
  float yValue;
  float zValue;

  if( !apply ) {
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * R11 + yTrace[i] * R21 + zTrace[i] * R31;
      yValue = xTrace[i] * R12 + yTrace[i] * R22 + zTrace[i] * R32;
      zValue = xTrace[i] * R13 + yTrace[i] * R23 + zTrace[i] * R33;
      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
  // Apply: rotate coordinate system from x(hor) y(hor) z(vert/down) to Galperin xyz
  // All components after rotation point downwards, new x component has same azimuth as original x component
  else {  // APPLY
    for( int i = 0; i < nSamples; i++ ) {
      xValue = xTrace[i] * R11 + yTrace[i] * R12 + zTrace[i] * R13;
      yValue = xTrace[i] * R21 + yTrace[i] * R22 + zTrace[i] * R23;
      zValue = xTrace[i] * R31 + yTrace[i] * R32 + zTrace[i] * R33;
      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }

}


/**
 * Rotate coordinate system in horizontal plane by given rotation angle
 * Apply clockwise rotation to coordinate system, assuming 2D coordinate system such as X/East, Y/North
 *
 * @param xTrace
 * @param yTrace
 * @param angle_rot  Rotation angle in radians
 * @param nSamples   Number of samples
*/
void rotate_xy( float* xTrace,
                float* yTrace,
                float angle_rot,
                int nSamples )
{
  float xValue;
  float yValue;

  float sina = sin( angle_rot );
  float cosa = cos( angle_rot );

  for( int i = 0; i < nSamples; i++ ) {
    xValue =  xTrace[i] * cosa - yTrace[i] * sina;
    yValue =  xTrace[i] * sina + yTrace[i] * cosa;
    xTrace[i] = xValue;
    yTrace[i] = yValue;
  }
}


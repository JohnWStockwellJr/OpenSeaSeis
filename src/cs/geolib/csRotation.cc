/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csRotation.h"
#include <cmath>
#include <cstdio>

using namespace cseis_geolib;

//
// 15/01/08  Changed values (better accuracy)
//
// For Trilobit:
// R =  ( -0.81649658        0.40824829        0.40824829  )
//      (  0                -0.70710678        0.70710678 )
//      (  0.57735027        0.57735027        0.57735027  )
//
//  = 1 / [6] ( -2    1   1  )
//            (  0  -[3] [3] )
//            ( [2]  [2] [2] )
// [x] : Square root of x
//

// float a = 35.2644; // [deg]
// float b = 45.0;    // [deg]
double const csRotation::GALPERIN_R11 =  0.81649658; // cosa
double const csRotation::GALPERIN_R12 =  0.0;        // 0
double const csRotation::GALPERIN_R13 =  0.57735027; // sina
double const csRotation::GALPERIN_R21 = -0.40824829; // -sinb*sina
double const csRotation::GALPERIN_R22 =  0.70710678; // cosb
double const csRotation::GALPERIN_R23 =  0.57735027; // sinb*cosa
double const csRotation::GALPERIN_R31 = -0.40824829; // -cosb*sina
double const csRotation::GALPERIN_R32 = -0.70710678; // -sinb
double const csRotation::GALPERIN_R33 =  0.57735027; // cosb*cosa

csRotation::csRotation() {
  myNumSamples = 0;
  myMode    = ROT_MODE_APPLY;
  myRotType = ROT_TYPE_UNKNOWN;
}
csRotation::csRotation( int nSamples, int mode ) {
  myNumSamples = nSamples;
  myMode       = mode;
  myRotType    = ROT_TYPE_UNKNOWN;
}
csRotation::~csRotation() {
}
void csRotation::setMode( int rotationMode ) {
  myMode    = rotationMode;
}
void csRotation::setRotationType( int rotType ) {
  myRotType = rotType;
  if( myRotType == ROT_TYPE_GALPERIN ) {
    set_galperin();
  }
  else if( myRotType == ROT_TYPE_TRILOBIT ) {
    set_trilobit();
  }
}
void csRotation::setNumSamples( int nSamples ) {
  myNumSamples = nSamples;
}
//--------------------------------------------------------------------------
//
void csRotation::rotate3d_weights( float* xTrace, float* yTrace, float* zTrace,
  float* weight_i, float* weight_c, float* weight_v )
{
  R11 = weight_i[0];
  R12 = weight_c[0];
  R13 = weight_v[0];
  R21 = weight_i[1];
  R22 = weight_c[1];
  R23 = weight_v[1];
  R31 = weight_i[2];
  R32 = weight_c[2];
  R33 = weight_v[2];

  myRotType = ROT_TYPE_UNKNOWN;

  rotate3d( xTrace, yTrace, zTrace );
}

void csRotation::rotate2d_azim( float* xTrace, float* yTrace, float azim )
{
  R11 = cos(azim);
  R21 = sin(azim);
  R12 = -R21;
  R22 = R11;

  myRotType = ROT_TYPE_UNKNOWN;

  rotate2d_xy( xTrace, yTrace );
}

void csRotation::rotate3d_roll_tilt( float* xTrace, float* yTrace, float* zTrace, float roll, float tilt )
{
  R11 = cos(tilt);
  R21 = sin(tilt)*cos(roll);
  R31 = sin(tilt)*sin(roll);

  R12 = 0.0;
  R22 = sin(roll);
  R32 = -cos(roll);

  R13 = -sin(tilt);
  R23 = cos(tilt)*cos(roll);
  R33 = cos(tilt)*sin(roll);

  myRotType = ROT_TYPE_UNKNOWN;

  rotate3d( xTrace, yTrace, zTrace );
}

void csRotation::rotate3d_galperin( float* xTrace, float* yTrace, float* zTrace )
{
  if( myRotType != ROT_TYPE_GALPERIN ) {
    set_galperin();
    myRotType = ROT_TYPE_GALPERIN;
  }
  rotate3d( xTrace, yTrace, zTrace );
}

void csRotation::rotate3d_armss( float* xTrace, float* yTrace, float* zTrace, float armss_roll, float armss_tilt ) {
  double sinr = sin(armss_roll);
  double cosr = cos(armss_roll);
  double sint = sin(armss_tilt);
  double cost = cos(armss_tilt);
/*
  // Galperin and azimuth first, then roll & tilt. NOTE: This will only work for mode REMOVE!
  rotate3d_armss_raw( xTrace, yTrace, zTrace );
  R11 = cost;
  R12 = -sinr*sint;
  R13 = cosr*sint;
  R21 = 0.0;
  R22 = cosr;
  R23 = sinr;
  R31 = -sint;
  R32 = -sinr*cost;
  R33 = cosr*cost;

*/
  // All together
  double a = 162.0;
  double a_rad = a * M_PI / 180.0;
  double sina = sin(a_rad);
  double cosa = cos(a_rad);

  double a11 = GALPERIN_R11*cosa + GALPERIN_R12*sina;
  double a12 = GALPERIN_R11*(-sina) + GALPERIN_R12*cosa;
  double a13 = GALPERIN_R13;
  double a21 = GALPERIN_R21*cosa + GALPERIN_R22*sina;
  double a22 = GALPERIN_R21*(-sina) + GALPERIN_R22*cosa;
  double a23 = GALPERIN_R23;
  double a31 = GALPERIN_R31*cosa + GALPERIN_R32*sina;
  double a32 = GALPERIN_R31*(-sina) + GALPERIN_R32*cosa;
  double a33 = GALPERIN_R33;

  double b11 = cost;
  double b12 = 0.0;
  double b13 = sint;
  double b21 = -sinr*sint;
  double b22 = cosr;
  double b23 = sinr*cost;
  double b31 = -cosr*sint;
  double b32 = -sinr;
  double b33 = cosr*cost;

  R11 = a11*b11 + a12*b21 + a13*b31;
  R12 = a11*b12 + a12*b22 + a13*b32;
  R13 = a11*b13 + a12*b23 + a13*b33;

  R21 = a21*b11 + a22*b21 + a23*b31;
  R22 = a21*b12 + a22*b22 + a23*b32;
  R23 = a21*b13 + a22*b23 + a23*b33;

  R31 = a31*b11 + a32*b21 + a33*b31;
  R32 = a31*b12 + a32*b22 + a33*b32;
  R33 = a31*b13 + a32*b23 + a33*b33;

  myRotType = ROT_TYPE_UNKNOWN;

  rotate3d( xTrace, yTrace, zTrace );

/*
  // TEMP: Roll only
  R11 = 1.0;
  R12 = 0.0;
  R13 = 0.0;
  R21 = 0.0;
  R22 = cosr;
  R23 = sinr;
  R31 = 0.0;
  R32 = -sinr;
  R33 = cosr;

  rotate3d( xTrace, yTrace, zTrace );

  // TEMP: Tilt only
  R11 = cost;
  R12 = 0.0;
  R13 = sint;
  R21 = 0.0;
  R22 = 1.0;
  R23 = 0.0;
  R31 = -sint;
  R32 = 0.0;
  R33 = cost;
  rotate3d( xTrace, yTrace, zTrace );
*/
}
void csRotation::rotate3d_armss( float* xTrace, float* yTrace, float* zTrace )
{
  // ARMSS rotation, APPLY:
  // 1) Rotation in XY plane by 162deg, anti-clock-wise rotation of coordinate system
  // 2) Galperin rotation
  double a = 162.0;
  double a_rad = a * M_PI / 180.0;
  double sina = sin(a_rad);
  double cosa = cos(a_rad);

  // Forward rotation / Apply:
  // R = R_galperin * (  cos(a) sin(a) )
  //                  ( -sin(a) cos(a) )

  R11 = GALPERIN_R11*cosa + GALPERIN_R12*sina;
  R12 = GALPERIN_R11*(-sina) + GALPERIN_R12*cosa;
  R13 = GALPERIN_R13;
  R21 = GALPERIN_R21*cosa + GALPERIN_R22*sina;
  R22 = GALPERIN_R21*(-sina) + GALPERIN_R22*cosa;
  R23 = GALPERIN_R23;
  R31 = GALPERIN_R31*cosa + GALPERIN_R32*sina;
  R32 = GALPERIN_R31*(-sina) + GALPERIN_R32*cosa;
  R33 = GALPERIN_R33;

  myRotType = ROT_TYPE_UNKNOWN;

  rotate3d( xTrace, yTrace, zTrace );
}

void csRotation::rotate3d_trilobit( float* xTrace, float* yTrace, float* zTrace ) {
  rotate3d( xTrace, yTrace, zTrace );
}

//**************************************************************************
//
//
//**************************************************************************

void csRotation::rotate3d( float* xTrace, float* yTrace, float* zTrace )
{
  float xValue;
  float yValue;
  float zValue;

  if( myMode == ROT_MODE_REMOVE ) {
    for( int i = 0; i < myNumSamples; i++ ) {
      xValue = (float)( xTrace[i] * R11 + yTrace[i] * R21 + zTrace[i] * R31 );
      yValue = (float)( xTrace[i] * R12 + yTrace[i] * R22 + zTrace[i] * R32 );
      zValue = (float)( xTrace[i] * R13 + yTrace[i] * R23 + zTrace[i] * R33 );
      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
  // Apply: rotate coordinate system from x(hor) y(hor) z(vert/down) to Galperin xyz
  // All components after rotation point downwards, new x component has same azimuth as original x component
  else {  // ROT_MODE_APPLY
    for( int i = 0; i < myNumSamples; i++ ) {
      xValue = (float)( xTrace[i] * R11 + yTrace[i] * R12 + zTrace[i] * R13 );
      yValue = (float)( xTrace[i] * R21 + yTrace[i] * R22 + zTrace[i] * R23 );
      zValue = (float)( xTrace[i] * R31 + yTrace[i] * R32 + zTrace[i] * R33 );
      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
}
void csRotation::rotate2d_xy( float* xTrace, float* yTrace )
{
  float xValue;
  float yValue;

  if( myMode == ROT_MODE_REMOVE ) {
    for( int i = 0; i < myNumSamples; i++ ) {
      xValue = (float)( xTrace[i] * R11 + yTrace[i] * R21 );
      yValue = (float)( xTrace[i] * R12 + yTrace[i] * R22 );
      xTrace[i] = xValue;
      yTrace[i] = yValue;
    }
  }
  else {  // if APPLY
    for( int i = 0; i < myNumSamples; i++ ) {
      xValue = (float)( xTrace[i] * R11 + yTrace[i] * R12 );
      yValue = (float)( xTrace[i] * R21 + yTrace[i] * R22 );
      xTrace[i] = xValue;
      yTrace[i] = yValue;
    }
  }
}

void csRotation::set_galperin() {
  R11 = GALPERIN_R11;
  R12 = GALPERIN_R12;
  R13 = GALPERIN_R13;
  R21 = GALPERIN_R21;
  R22 = GALPERIN_R22;
  R23 = GALPERIN_R23;
  R31 = GALPERIN_R31;
  R32 = GALPERIN_R32;
  R33 = GALPERIN_R33;
}
void csRotation::set_trilobit() {
  R11 = -GALPERIN_R11;
  R12 = -GALPERIN_R12;
  R13 = GALPERIN_R13;
  R21 = -GALPERIN_R21;
  R22 = -GALPERIN_R22;
  R23 = GALPERIN_R23;
  R31 = -GALPERIN_R31;
  R32 = -GALPERIN_R32;
  R33 = GALPERIN_R33;
}

/**
* Rotate XYZ vectors using ICV matrix. This only affects roll & tilt
*/
void csRotation::rotate_xyz_icv( float* xTrace,
         float* yTrace,
         float* zTrace,
         float* weight_i,
         float* weight_c,
         float* weight_v )
{
  float xValue;
  float yValue;
  float zValue;

  //  fprintf(stdout,"modeIN : %d, xyzTrace[0]: %f %f %f\n", mode, xTrace[0], yTrace[0], zTrace[0] );

  if( myMode == csRotation::ROT_MODE_REMOVE ) {
    // Apply transposed rotation matrix
    for( int i = 0; i < myNumSamples; i++ ) {
      xValue = xTrace[i]*weight_i[0] + yTrace[i]*weight_i[1] + zTrace[i]*weight_i[2];
      yValue = xTrace[i]*weight_c[0] + yTrace[i]*weight_c[1] + zTrace[i]*weight_c[2];
      zValue = xTrace[i]*weight_v[0] + yTrace[i]*weight_v[1] + zTrace[i]*weight_v[2];

      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
  else {  // if APPLY
    for( int i = 0; i < myNumSamples; i++ ) {
      xValue = xTrace[i]*weight_i[0] + yTrace[i]*weight_c[0] + zTrace[i]*weight_v[0];
      yValue = xTrace[i]*weight_i[1] + yTrace[i]*weight_c[1] + zTrace[i]*weight_v[1];
      zValue = xTrace[i]*weight_i[2] + yTrace[i]*weight_c[2] + zTrace[i]*weight_v[2];

      xTrace[i] = xValue;
      yTrace[i] = yValue;
      zTrace[i] = zValue;
    }
  }
}


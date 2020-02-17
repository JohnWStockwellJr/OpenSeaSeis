/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_ROTATION_H
#define CS_ROTATION_H

namespace cseis_geolib {

  /**
   * Perform XYZ rotation.
   *
   * - 2D, 3D...
   * @author Bjorn Olofsson
   * @date 2008
   */
class csRotation {
public:
  static int const ROT_MODE_APPLY    = 0;
  static int const ROT_MODE_REMOVE   = 1;

  static int const ROT_TYPE_AZIM     = 2;
  static int const ROT_TYPE_GALPERIN = 3;
  static int const ROT_TYPE_ARMSS    = 4;
  static int const ROT_TYPE_TRILOBIT = 5;
  static int const ROT_TYPE_UNKNOWN  = -1;

  // float a = 35.2644; // [deg]
  // float b = 45.0;    // [deg]
  static double const GALPERIN_R11;// = 0.81649648;  // cosa
  static double const GALPERIN_R12;// = 0.0;         // 0
  static double const GALPERIN_R13;// = 0.57735042;  // sina
  static double const GALPERIN_R21;// = -0.40824839; // -sinb*sina
  static double const GALPERIN_R22;// = 0.70710678;  // cosb
  static double const GALPERIN_R23;// = 0.57735020;  // sinb*cosa
  static double const GALPERIN_R31;// = -0.40824839; // -cosb*sina
  static double const GALPERIN_R32;// = -0.70710678; // -sinb
  static double const GALPERIN_R33;// = 0.57735020;  // cosb*cosa

public:
  csRotation();
  csRotation( int nSamples, int mode );
  ~csRotation();
  void setMode( int rotationMode );
  void setRotationType( int rotType );
  void setNumSamples( int nSamples );

  void rotate3d_weights( float* xTrace, float* yTrace, float* zTrace,
    float* weight_i, float* weight_c, float* weight_v );
  void rotate2d_azim( float* xTrace, float* yTrace, float azim );
  void rotate3d_roll_tilt( float* xTrace, float* yTrace, float* zTrace, float roll, float tilt );

  void rotate3d_galperin( float* xTrace, float* yTrace, float* zTrace );
  void rotate3d_armss( float* xTrace, float* yTrace, float* zTrace );
  void rotate3d_armss( float* xTrace, float* yTrace, float* zTrace, float armss_roll, float armss_tilt );
  void rotate3d_trilobit( float* xTrace, float* yTrace, float* zTrace );

  void rotate_xyz_icv( float* xTrace, float* yTrace, float* zTrace,
    float* weight_i, float* weight_c, float* weight_v );

private:
  void set_galperin();
  void set_armss();
  void set_trilobit();
  void rotate3d( float* xTrace, float* yTrace, float* zTrace );
  void rotate2d_xy( float* xTrace, float* yTrace );

  int myNumSamples;
  int myMode;
  int myRotType;

// 3D rotation matrix
  double R11;
  double R12;
  double R13;
  double R21;
  double R22;
  double R23;
  double R31;
  double R32;
  double R33;
};

} // namespace
#endif


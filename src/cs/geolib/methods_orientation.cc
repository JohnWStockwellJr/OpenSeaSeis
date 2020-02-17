/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include "geolib_math.h"
#include <cstdio>
#include <string>

namespace cseis_geolib {

/**
* Compute receiver orientation azimuth
* from recorded polarisation vector and source-receiver direction vector
*
* @return azimuth [rad]
*/
void orientation_compute_azimuth(
  float const* vec_polarisation,  // Unit polarisation vector in as-laid coordinate system
  double rcvx,     // Receiver X position [m]
  double rcvy,     // Receiver Y position [m]
  double srcx,     // Source X position [m]
  double srcy,     // Source Y position [m]
  float* angle_azim // Computed azimuth [rad]
)
{
  float dx = (float)(rcvx - srcx);
  float dy = (float)(rcvy - srcy);
  float angle_sr_azim = atan2( dx, dy );
  if( angle_sr_azim < 0.0 ) {
    angle_sr_azim += M_PI*2.0;
  }

  float vecx = cos( angle_sr_azim );  // Northing
  float vecy = sin( angle_sr_azim );  // Easting

  *angle_azim = atan2( vec_polarisation[0]*vecy - vec_polarisation[1]*vecx,
                       vec_polarisation[0]*vecx + vec_polarisation[1]*vecy );
}

/**
* Compute orientation angles (Roll and Tilt) from direct arrival polarisation (vec_polarisation)
* and source-receiver directin vector (vec_direction)
*
* To use incidence angle, set v1v2_ratio = 1.0
* To use transmission angle, set v1v2_ratio != 1.0
*/
void orientation_compute_roll_tilt (
  float const* vec_polarisation, // Unit polarisation vector in as-laid coordinate system
  double rcvx,     // Receiver X position [m]
  double rcvy,     // Receiver Y position [m]
  double rcvz,     // Receiver Z position [m]
  double srcx,     // Source X position [m]
  double srcy,     // Source Y position [m]
  double srcz,     // Source Z position [m]
  float v1v2_ratio,         // Velocity ratio at seabed
  int   isOBCMode,          // true(1) if source is above receiver, false(0) if source is beneath receiver
  float angle_azim,         // Receiver azimuth [rad]
  float* angle_tilt,        // Tilt angle [rad]
  float* angle_roll )       // Roll angle [rad]
{
  float dx = (float)(rcvx - srcx);
  float dy = (float)(rcvy - srcy);
  float dz = (float)(rcvz - srcz);

  float daoffset = (float)std::max( sqrt(dx*dx + dy*dy + dz*dz), 0.0001 );
  float angle_incidence = acos(dz/daoffset);  // Always positive, assumes dz > 0
  float angle_transmitted = sin(angle_incidence) / v1v2_ratio;
  if( fabs(angle_transmitted) <= 1.0 ) {
    angle_transmitted = asin( angle_transmitted );  // Always positive
  }
  else {
    angle_transmitted = M_PI_2;
  }

  float angle_sr_azim = atan2( dx, dy );
  if( angle_sr_azim < 0.0 ) {
    angle_sr_azim += M_PI*2.0;
  }

  // Compute unit source vector
  float vec_direction[3];
  vec_direction[0] = sin( angle_transmitted ) * cos( angle_sr_azim );  // Northing
  vec_direction[1] = sin( angle_transmitted ) * sin( angle_sr_azim );  // Easting
  vec_direction[2] = cos( angle_transmitted );  // Z, positive downwards

  float b_i_azim, b_c_azim, b_v_azim;  // Source vector, rotated by receiver azimuth
  float angle_epsilon, angle_kappa;
  float b_v_tilt_azim;  // Source vector, rotated by azimuth and tilt angle
  float b_vec_length;
  float tmp_tilt, tmp_roll;
  float temp;

  // Rotate source vector into inline/xline/vertical coordinate system
  b_i_azim =  cos(angle_azim) * vec_direction[0] + sin(angle_azim) * vec_direction[1];
  b_c_azim = -sin(angle_azim) * vec_direction[0] + cos(angle_azim) * vec_direction[1];
  b_v_azim = vec_direction[2];
  
  b_vec_length = sqrt( b_i_azim*b_i_azim + b_v_azim*b_v_azim );    // Length of unit source vector in inline/vertical plane
  angle_kappa = atan2( b_i_azim, b_v_azim );  // Incidence angle in inline/vertical plane. e [-90,90]
  temp = vec_polarisation[0] / b_vec_length;
  if( fabs(temp) <= 1.0 ) {
    angle_epsilon = acos( temp );  // Angle between source direction and cable
  }
  else {
    angle_epsilon = 0.0;
  }

  if( isOBCMode ) {
    tmp_tilt = angle_kappa - M_PI_2 + angle_epsilon;
  }
  else {
    tmp_tilt = angle_kappa - M_PI_2 - angle_epsilon;
    if( tmp_tilt > M_PI_2 ) {
      tmp_tilt -= M_PI*2.0;
    }
  }
  // Apply tilt rotation to source vector
  b_v_tilt_azim = sin(tmp_tilt) * b_i_azim + cos(tmp_tilt) * b_v_azim;

  // Roll angle (inversion of rotation equation)
  float top    = b_v_tilt_azim*vec_polarisation[2] + b_c_azim*vec_polarisation[1];
  float bottom = b_v_tilt_azim*vec_polarisation[1] - b_c_azim*vec_polarisation[2];
  tmp_roll  = atan2( top, bottom );
  // Wrap to 0-360deg
  if( tmp_roll < 0.0 ) {
    tmp_roll += 2.0*M_PI;
  }
  else if( tmp_roll >= 2.0*M_PI ) {
    tmp_roll -= 2.0*M_PI;
  }

  *angle_roll = tmp_roll;
  *angle_tilt = tmp_tilt;

  //  fprintf(stdout,"sr_azim: %9.3f  %9.3f  %9.4f %9.4f %9.4f %9.4f %9.4f\n", angle_sr_azim, angle_transmitted, angle_kappa, angle_epsilon,  b_v_tilt_azim, top, bottom );
}

} // namespace


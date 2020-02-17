/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_LINE_3D_H
#define CS_LINE_3D_H

#include <cstdio>
#include "csPoint3D.h"

namespace cseis_geolib {

/**
 * 3D Line representation.
 * Line is defined by origin point and direction vector
 */
class csLine3D {
public:
  csLine3D();
  csLine3D( cseis_geolib::csPoint3D const& originPoint, cseis_geolib::csPoint3D const& directionVector );
  virtual ~csLine3D();
  virtual void setFromPoints( cseis_geolib::csPoint3D const& p1, cseis_geolib::csPoint3D const& p2 );
  virtual cseis_geolib::csPoint3D pointAt( double lamda ) const;
  virtual csPoint3D computePointFromX( double x ) const;
  virtual csPoint3D computePointFromY( double y ) const;
  virtual csPoint3D computePointFromZ( double z ) const;
  /**
   * Compute distance of this line to a point in 3D space
   * @param: point: Point from which distance to line shall be computed 
   */
  virtual double distanceToPoint( csPoint3D const& point ) const;
  virtual cseis_geolib::csPoint3D p1() const { return origin; }
  virtual cseis_geolib::csPoint3D p2() const { return point2; }
  virtual cseis_geolib::csPoint3D delta() const { return direction; }
  /// Offset between end points (= horizontal distance)
  virtual double offset() const { return myOffset; }
  /// Distance between end points (=3D distance)
  virtual double length() const { return myLength; }
  /// Azimuth of line
  virtual double azimuth() const { return myAzimuth; }
  /**
   * Compute minimum distance between a line and a point in 3D space
   * @param linePoint1: First point defining the line
   * @param linePoint2: Second point defining the line
   * @param point:      Point from which distance to line shall be computed 
   */
  static double distanceLineToPoint( csPoint3D const& linePoint1, csPoint3D const& linePoint2, csPoint3D const& point );
  //  csPoint3D intersection( csLine3D const& line ) const;
public:
  csPoint3D origin;
  csPoint3D direction;
  csPoint3D point2;
protected:
  virtual void compute_internal();
  double myOffset;
  double myLength;
  double myAzimuth;
};

} // END namespace
#endif

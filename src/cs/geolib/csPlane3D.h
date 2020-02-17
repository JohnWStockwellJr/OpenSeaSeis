/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_PLANE_3D_H
#define CS_PLANE_3D_H

#include "csPoint3D.h"
#include "csPoint2D.h"
#include "csLine3D.h"
#include <cstdio>
 
namespace cseis_geolib {

template <typename T> class csVector;

struct rayPathParam {
  double distanceSrcCMP;  // Ray path distance from source to receiver
  double distanceRcvCMP;  // Ray path distance from receiver to CMP
  double offRatioSrcCMP;  // Normalized horizontal offset from source to CMP. 0: At source position, 1: At receiver position, 0.5: Mid-way between...
  cseis_geolib::csPoint3D cmp;
};

class csProjection2D {
 public:
  csProjection2D( cseis_geolib::csPoint3D const& xvec, cseis_geolib::csPoint3D const& yvec );
  cseis_geolib::csPoint2D apply( cseis_geolib::csPoint3D const& point ) const;
 private:
  cseis_geolib::csPoint3D myXVec;
  cseis_geolib::csPoint3D myYVec;
};

class csPlane3D {
public:
  csPlane3D();
  /**
   * @param normalVector Normal vector on plane
   * @param origin       Point on plane
   */
  csPlane3D( csPoint3D const& normalVector, csPoint3D const& origin );
  csPlane3D( csPoint3D const& normalVector, double d_in );
  csPlane3D( csPoint3D const& normalVector, csPoint3D const& origin, double d_in );
  ~csPlane3D();

  void set( csPoint3D const& normalVector, double d_in );
  void set( csPoint3D const& normalVector, csPoint3D const& origin, double d_in );
  void setOrigin( double x, double y, double z );
  
  void setFromPoints( csPoint3D p1, csPoint3D p2, csPoint3D p3 );
  bool fitToPoints( cseis_geolib::csVector<cseis_geolib::csPoint3D*> const* pList, std::FILE* errStream );
  bool fitToPoints( double const* xval, double const* yval, double const* zval, int npoints, FILE* errStream );
  
  double computeZ( double x, double y ) const;
  void computeZ( csPoint3D& p ) const;
  double distanceToPoint( csPoint3D const& point ) const;
  double dipAngle() const;
  double azimuthAngle() const;
  csPoint3D normalVector() const;
  csPoint3D origin() const;
  double getA() const;
  double getB() const;
  double getC() const;
  double getD() const;
  cseis_geolib::csLine3D intersection( csPlane3D const& plane2 ) const;
  cseis_geolib::csPoint3D intersection( csLine3D const& line ) const;
  /**
   * Project point onto plane
   * @param point: Input point
   * @return Projected point on plane
   */
  cseis_geolib::csPoint3D projectPoint( csPoint3D const& point ) const;
  /**
   * Project point onto 2D-equivalent plane
   */
  cseis_geolib::csPoint2D projectPoint2D( cseis_geolib::csPoint3D const& pointOnPlane );
  /**
   * Make sure the plane's normal vector is normalized to 1
   */
  void normalizeNormalVector();
  /**
   * Compute three random points on plane
   */
  void compute3points( cseis_geolib::csPoint3D& p1, cseis_geolib::csPoint3D& p2, cseis_geolib::csPoint3D& p3 ) const;
  
  /**
   * Compute reflection ray path (distance) for one source-receiver pair
   */
  rayPathParam computeReflectionPoint( cseis_geolib::csPoint3D const& pSrc, cseis_geolib::csPoint3D const& pRcv, std::FILE* debugStream = NULL ) const;

 private:
  csPlane3D( csPlane3D const& plane );
  bool prepareFitToPoints_internal( int npoints, std::FILE* errStream );
  bool fitToPoints_internal( int npoints, std::FILE* errStream );
  void computeDipAzim_from_normalVector();
  void prepareProjection2D();

  csPoint3D myNormalVector;
  csPoint3D myOrigPoint;
  double myD;
  
  double  myDip_deg;
  double  myAzimuth_deg;
  
  double** myMatU_temp;
  double*  myVecB_temp;

  csProjection2D* myProjection2D;
};

} // END namespace
#endif

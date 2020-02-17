/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csPlane3D.h"
#include "csReflector2D.h"
#include "csPoint3D.h"
#include "csPoint2D.h"
#include "csVector.h"
#include "geolib_methods.h"
#include <cmath>

using namespace cseis_geolib;
using namespace std;

csProjection2D::csProjection2D( cseis_geolib::csPoint3D const& xvec, cseis_geolib::csPoint3D const& yvec ) : myXVec( xvec ), myYVec( yvec ) {
}
cseis_geolib::csPoint2D csProjection2D::apply( cseis_geolib::csPoint3D const& point ) const {
  cseis_geolib::csPoint2D p2d( point.dotProduct(myXVec), point.dotProduct(myYVec) );
  return p2d;
}

csPlane3D::csPlane3D() {
  myOrigPoint.set( 0, 0, 0 );
  myNormalVector.set( 0, 0, 0 );
  myD = 0;

  myMatU_temp = NULL;
  myVecB_temp = NULL;
  myProjection2D = NULL;
}
csPlane3D::csPlane3D( csPoint3D const& normalVector, csPoint3D const& origin ) {
  myNormalVector = normalVector;
  myOrigPoint    = origin;
  myD = -( myNormalVector.x * myOrigPoint.x + myNormalVector.y * myOrigPoint.y + myNormalVector.z * myOrigPoint.z );

  myMatU_temp = NULL;
  myVecB_temp = NULL;
  myProjection2D = NULL;
}
csPlane3D::csPlane3D( csPoint3D const& normalVector, double d_in ) {
  myNormalVector = normalVector;
  myOrigPoint.set( 0, 0, 0 );
  myD = d_in;

  myMatU_temp = NULL;
  myVecB_temp = NULL;
  myProjection2D = NULL;
}
csPlane3D::csPlane3D( csPoint3D const& normalVector, csPoint3D const& origin, double d_in ) {
  myNormalVector = normalVector;
  myOrigPoint    = origin;
  myD = d_in;

  myMatU_temp = NULL;
  myVecB_temp = NULL;
  myProjection2D = NULL;
}
csPlane3D::~csPlane3D() {
}
void csPlane3D::set( csPoint3D const& normalVector, csPoint3D const& origin, double d_in ) {
  myNormalVector = normalVector;
  myOrigPoint = origin;
  myD = d_in;
}
void csPlane3D::set( csPoint3D const& normalVector, double d_in ) {
  myNormalVector = normalVector;
  myD = d_in;
}
void csPlane3D::setOrigin( double x0, double y0, double z0 ) {
  myOrigPoint.x = x0;
  myOrigPoint.y = y0;
  myOrigPoint.z = z0;
}

double csPlane3D::computeZ( double x, double y ) const {
  if( abs(myNormalVector.z) < 1.0e-12 ) return myD;
  double temp = x * myNormalVector.x  +  y * myNormalVector.y + myD;
  double z = (-temp/myNormalVector.z);
  return z;
}
void csPlane3D::computeZ( csPoint3D& p ) const {
  p.z = computeZ( p.x, p.y );
}
double csPlane3D::distanceToPoint( csPoint3D const& point ) const {
  // a) Compute vector connecting given point and plane origin point: vec  = point - myNormalVector
  // b) Project vector onto normal vector by using dot product:       dot  = myNormalVector * vec
  // c) Normalize by length of normal vector:                         dist = dot / myNormalVector.vectorLentgh() 
  return( abs(myNormalVector.dotProduct( myOrigPoint ) - myNormalVector.dotProduct( point )) / myNormalVector.vectorLength() );
}
//double csPlane3D::computeFromXZ( double x, double z ) const {
//  double z = (x-myOrigPoint.x) * myCoefficients[0]  +  (y-myOrigPoint.y) * myCoefficients[1]  +  myCoefficients[2];
//  return z;
//}
double csPlane3D::dipAngle() const {
  return myDip_deg;
}
double csPlane3D::azimuthAngle() const {
  return myAzimuth_deg;
}
csPoint3D csPlane3D::normalVector() const {
  return myNormalVector;
}
csPoint3D csPlane3D::origin() const {
  return myOrigPoint;
}
double csPlane3D::getA() const {
  return myNormalVector.x;
}
double csPlane3D::getB() const {
  return myNormalVector.y;
}
double csPlane3D::getC() const {
  return myNormalVector.z;
}
double csPlane3D::getD() const {
  return myD;
}

void csPlane3D::computeDipAzim_from_normalVector() {
  csPoint3D vecZ( 0, 0, -1 );
  myDip_deg = acos(  vecZ.dotProduct( myNormalVector )  ) * 180.0 / M_PI; // No need to normalize since both vectors have length 1.0
  myAzimuth_deg = fmod( 360.0 + atan2( myNormalVector.x, myNormalVector.y ) * 180.0 / M_PI, 360.0 );
}
void csPlane3D::setFromPoints( csPoint3D p1, csPoint3D p2, csPoint3D p3 ) {
  myOrigPoint = p1;
  csPoint3D vec12( p2.x-p1.x, p2.y-p1.y, p2.z-p1.z );
  csPoint3D vec13( p3.x-p1.x, p3.y-p1.y, p3.z-p1.z );

  myNormalVector    = vec12.crossProduct( vec13 );
  double normScalar = myNormalVector.vectorLength();
  if( myNormalVector.z > 0 ) normScalar *= -1.0; // Make sure normal vector sits 'on top' of reflector, i.e. it points upwards for sub-horizontal planes
  myNormalVector   /= normScalar;

  myD = -( myNormalVector.x * p1.x + myNormalVector.y * p1.y + myNormalVector.z * p1.z );

  // Compute angles:
  computeDipAzim_from_normalVector();
  //  csPoint3D vecZ( 0, 0, -1 );
  //  myDip_deg = acos(  vecZ.dotProduct( myNormalVector )  ) * 180.0 / M_PI; // No need to normalize since both vectors have length 1.0
  //  myAzimuth_deg = fmod( 360.0 + atan2( myNormalVector.x, myNormalVector.y ) * 180.0 / M_PI, 360.0 );
}
bool csPlane3D::fitToPoints( double const* xval, double const* yval, double const* zval, int npoints, FILE* errStream ) {
  bool ret = prepareFitToPoints_internal( npoints, errStream );
  if( !ret ) return ret;

  for( int ip = 0; ip < npoints; ip++ ) {
    myMatU_temp[ip][0] = xval[ip]; // A
    myMatU_temp[ip][1] = yval[ip]; // B
    myMatU_temp[ip][2] = 1.0; // D
    //    myMatU_temp[ip][3] = 1.0; // D
    myVecB_temp[ip]    = -zval[ip]; // C = 1
  }

  return fitToPoints_internal( npoints, errStream );
}
bool csPlane3D::fitToPoints( cseis_geolib::csVector<cseis_geolib::csPoint3D*> const* pList, FILE* errStream ) {
  int npoints = pList->size();
  bool ret = prepareFitToPoints_internal( npoints, errStream );
  if( !ret ) return ret;

  for( int ip = 0; ip < npoints; ip++ ) {
    csPoint3D* p = pList->at(ip);
    myMatU_temp[ip][0] = p->x;
    myMatU_temp[ip][1] = p->y;
    myMatU_temp[ip][2] = 1.0;
    //    myMatU_temp[ip][3] = 1.0;
    myVecB_temp[ip]    = -p->z;
  }
  return fitToPoints_internal( npoints, errStream );
}

bool csPlane3D::prepareFitToPoints_internal( int npoints, FILE* errStream ) {
  int nCols = 3;
  if( npoints < nCols ) {
    if( errStream != NULL ) fprintf(errStream,"csPlane3D::computeFromPoints: Too few points provided (=%d). Need at least %d points.\n", npoints, nCols);
    return false;
  }
  myMatU_temp = new double*[npoints];
  for( int i = 0; i < npoints; i++ ) {
    myMatU_temp[i] = new double[nCols];;
  }
  myVecB_temp = new double[npoints];
  return true;
}

bool csPlane3D::fitToPoints_internal( int npoints, FILE* errStream ) {
  int nCols = 3;

  double** mat_v = new double*[nCols];
  for( int i = 0; i < nCols; i++ ) {
    mat_v[i] = new double[nCols];;
  }
  double* vec_w    = new double[npoints];
  double* vec_work = new double[npoints];

  int ret = svd_decomposition( myMatU_temp, npoints, nCols, vec_w, mat_v, vec_work );
  if( ret == 0 ) {
    if( errStream != NULL ) fprintf(errStream,"csPlane3D::fitToPoints: Unknown error occurred in SVD\n");
  }
  else {
    double coefficients[3];
    svd_linsolve( myMatU_temp, vec_w, mat_v, npoints, nCols, myVecB_temp, coefficients, vec_work );
    myNormalVector.x = coefficients[0];
    myNormalVector.y = coefficients[1];
    myNormalVector.z = 1.0; //coefficients[2];
    myD = coefficients[2];
    // Normal vector on 3D plane = ( A B C ), constant D
    // dz/dx = A = normVec.x
    // dz/dy = B = normVec.y
    // dz/dz = C = normVec.z
    // z0    = D = myD

    double length = myNormalVector.vectorLength();
    myNormalVector /= -length; // Normalize normal vector and flip it so that it point upwards, at least for sub-horizontal planes
    myD /= -length; // D must be normalized in the same way

    computeDipAzim_from_normalVector();
    //    fprintf(stderr,"Normal vector: %f %f %f   dip/azim: %f %f\n", myNormalVector.x, myNormalVector.y, myNormalVector.z, myDip_deg, myAzimuth_deg );
  }

  for( int i = 0; i < npoints; i++ ) {
    delete [] myMatU_temp[i];
  }
  for( int i = 0; i < nCols; i++ ) {
    delete [] mat_v[i];
  }
  delete [] myMatU_temp;
  delete [] mat_v;
  delete [] vec_w;
  delete [] vec_work;
  delete [] myVecB_temp;

  if( ret == 0 ) return false;
  return true;
}

cseis_geolib::csLine3D csPlane3D::intersection( csPlane3D const& plane2 ) const {
  // Not tested yet!
  cseis_geolib::csLine3D line;
  line.direction = myNormalVector.crossProduct( plane2.myNormalVector );
  double dotProduct = line.direction.dotProduct( line.direction );

  csPoint3D u1( myNormalVector * plane2.myD );
  csPoint3D u2( plane2.myNormalVector * (-myD));

  csPoint3D pOrig( u1 + u2 );
  pOrig = pOrig.crossProduct( line.direction );
  pOrig /= dotProduct;
  line.origin = pOrig;

  return line;
}
cseis_geolib::csPoint3D csPlane3D::intersection( csLine3D const& line ) const {
  cseis_geolib::csPoint3D point( 0, 0, 0 );

  double upper = -( myNormalVector.dotProduct( line.origin ) + myD );
  double lower = myNormalVector.dotProduct( line.direction );
  if( abs(lower) < 1.0e-12 ) return point;
  double lamda = upper / lower;
  point.x = line.origin.x + lamda * line.direction.x;
  point.y = line.origin.y + lamda * line.direction.y;
  point.z = line.origin.z + lamda * line.direction.z;
  return point;
}
cseis_geolib::csPoint3D csPlane3D::projectPoint( csPoint3D const& point ) const {
  csPoint3D vec1 = myOrigPoint - point;
  double lamda = vec1.dotProduct( myNormalVector );
  csPoint3D vec2 = myNormalVector * lamda;
  csPoint3D pointProjected = point + vec2;
  return pointProjected;
}

void csPlane3D::prepareProjection2D() {
  normalizeNormalVector();
  csPoint3D p1;
  csPoint3D p2;
  csPoint3D p3;
  compute3points( p1, p2, p3 );
  csPoint3D xvec   = p2 - myOrigPoint;
  double vecLength = xvec.vectorLength();
  if( vecLength < 1.0e-5 ) {
    xvec  = p3 - myOrigPoint;
    vecLength = xvec.vectorLength();
  }
  if( vecLength > 1.0e-9 ) {
    xvec /= vecLength;
  }
  csPoint3D yvec = myNormalVector.crossProduct( xvec );
  if( myProjection2D != NULL ) delete myProjection2D;
  myProjection2D = new csProjection2D( xvec, yvec );
}

/**
 * Compute plane 2D projection point
 */
cseis_geolib::csPoint2D csPlane3D::projectPoint2D( cseis_geolib::csPoint3D const& pointOnPlane ) {
  if( myProjection2D == NULL ) prepareProjection2D();
  return( myProjection2D->apply( pointOnPlane ) );
}

void csPlane3D::normalizeNormalVector() {
  double length = myNormalVector.vectorLength();
  if( length != 0.0 ) {
    myNormalVector /= length;
    myD /= length; // D must be normalized in the same way
  }
}

void csPlane3D::compute3points( csPoint3D& p1, csPoint3D& p2, csPoint3D& p3 ) const {
  csPoint3D point(0,0,0);
  p1 = myOrigPoint;
  if( abs(myNormalVector.y) < 1.0e-12 || abs(myNormalVector.z) < 1.0e-12 ) {
    point.x = 1.0;
  }
  else {
    point.y = 1.0;
  }
  p2 = myNormalVector.crossProduct( point );
  p3 = myNormalVector.crossProduct( p2 );
  p2 *= 100;
  p3 *= 100;
  p2 = p2 + myOrigPoint;
  p3 = p3 + myOrigPoint;
}
//--------------------------------------------------------------------------------
//
//
rayPathParam csPlane3D::computeReflectionPoint( cseis_geolib::csPoint3D const& pSrc, cseis_geolib::csPoint3D const& pRcv, std::FILE* debugStream ) const {
  //rayPathParam computeRayPath( cseis_geolib::csPlane3D const& planeRefl, cseis_geolib::csPoint3D const& pSrc, cseis_geolib::csPoint3D const& pRcv, std::FILE* rayStream, std::FILE* stream, bool debug ) {
  double nominalDist = 1000.0;

  // 1) Compute plane that passes through source & receiver position and is normal to the plane reflector
  csPoint3D vecNormRefl = normalVector();

  csPoint3D tmp1( ( pSrc + pRcv ) * 0.5 ); // Point at mid position between source & receiver, at the surface
  csPoint3D tmp2( vecNormRefl * nominalDist );  // Elongate the reflector plane's normal vector by nominal distance
  csPoint3D p3( tmp1 + tmp2 );  // Compute third point of ray path plane by moving mid point down along reflector normal by a nominal distance

  csPlane3D planeRayPath;
  planeRayPath.setFromPoints( pSrc, pRcv, p3 );
  csPoint3D vecNormRay = planeRayPath.normalVector();

  // 2) Compute distance from source and receiver points to reflector plane
  //  - Compute reflector dip in ray path 2D plane (source-receiver-reflector plane)
  double distReflSrc = distanceToPoint( pSrc );
  double distReflRcv = distanceToPoint( pRcv );

  double distReflDiff = abs( distReflSrc - distReflRcv );
  double distSR       = pSrc.distance( pRcv );
  double dipAngleRefl2D_rad = asin( distReflDiff / distSR );
  double dipAngleRefl2D_deg = dipAngleRefl2D_rad * 180.0 / M_PI;

  // 3) Construct 2D geometry: Source, receiver, dipping reflector
  // Point on reflector:
  double minDist2Refl = std::min( distReflSrc, distReflRcv );
  double maxDist2Refl = std::max( distReflSrc, distReflRcv );
  csPoint2D p2d_src( 0, 0 );
  csPoint2D p2d_rcv( distSR, 0 );
  csPoint2D p2d_orig( p2d_src.x + sin(dipAngleRefl2D_rad) * maxDist2Refl, p2d_src.z + cos(dipAngleRefl2D_rad) * maxDist2Refl );
  csReflector2D refl2d( p2d_orig, -dipAngleRefl2D_deg );  // Using the negative dip angble makes sure that the reflector is deepest below the 2D 'source' and dipping up to the 2D 'receiver'

  // 4) Prepare iteration: Define left-most and right-most points on reflector to test
  //   - Compute normal vector on 2D reflector
  double xmin = 0;
  double xmax = distSR + sin(dipAngleRefl2D_rad) * minDist2Refl;
  csPoint2D p2d_left( xmin, refl2d.compute( xmin ) );
  csPoint2D p2d_right( xmax, refl2d.compute( xmax ) );
  csPoint2D p2d_mid;
  csPoint2D vecNormRefl2d( refl2d.direction().z, -refl2d.direction().x );

  if( debugStream ) {
    fprintf(debugStream,"%f %f %f RAY PATH POINT3\n", p3.x, p3.y, p3.z);
    fprintf(debugStream,"REFL_PLANE_ABCD      %12.6f %12.6f %12.6f %12.6f\n", getA(), getB(), getC(), getD() );
    fprintf(debugStream,"RAYPATH_PLANE_ABCD   %12.6f %12.6f %12.6f %12.6f\n", planeRayPath.getA(), planeRayPath.getB(), planeRayPath.getC(), planeRayPath.getD() );
    fprintf(debugStream,"REFL_PLANE_NORMAL    %12.6f %12.6f %12.6f\n", vecNormRefl.x, vecNormRefl.y, vecNormRefl.z);
    fprintf(debugStream,"RAYPATH_PLANE_NORMAL %12.6f %12.6f %12.6f\n", vecNormRay.x, vecNormRay.y, vecNormRay.z);
    fprintf(debugStream,"MISC_Distance_from_source_receiver_points  %.3f  %.3f   SR distance: %.2f\n", distReflSrc, distReflRcv, distSR );
    fprintf(debugStream,"MISC_Distance_difference  %.3f  Dip angle: %.3f\n", distReflDiff, dipAngleRefl2D_deg );
    fprintf(debugStream,"MISC_Source_Receiver_locations  %.2f %.2f   lineOrig: %.2f %.2f\n", p2d_src.x, p2d_rcv.x, p2d_orig.x, p2d_orig.z );

    int num = 20;
    double dx = (xmax - xmin ) / (num-1);
    for( int i = 0; i < num; i++ ) {
      double x = i*dx + xmin;
      double y = refl2d.compute( x );
      fprintf(debugStream,"REFL_PLANE2D  %.2f %.2f\n", x, y);
    }
  }

  double angleDiff = 999;
  double angleThreshold = 0.001;
  int maxIter = 100;
  int counter = 0;

  // 5) Iteration: Test mid point: Does it have same incidence(src) & reflection(rcv) angle?
  //  - Stop iteration when difference between the two angles falls below threshold
  do {
    p2d_mid.x = 0.5 * ( p2d_left.x + p2d_right.x );
    p2d_mid.z = refl2d.compute( p2d_mid.x );

    csPoint2D vecLeft( p2d_src - p2d_mid );
    csPoint2D vecRight( p2d_rcv - p2d_mid );

    double angleLeft  = vecNormRefl2d.vectorAngle( vecLeft );
    double angleRight = vecNormRefl2d.vectorAngle( vecRight );

    angleDiff = angleLeft - angleRight;

    if( angleDiff > 0 ) {
      p2d_right = p2d_mid;
    }
    else {
      p2d_left = p2d_mid;
    }
    counter += 1;
    //    fprintf(stream,"Iteration: #%d  AngleDiff: %f    Midpoint: %f %f   Angles: %.4f %.4f\n", counter, angleDiff, p2d_mid.x, p2d_mid.z, angleLeft, angleRight);
  } while( abs(angleDiff) > angleThreshold || counter > maxIter );

  //csPoint2D pfake( p2d_src.x, 500 );
  // Reflector2D line2d_sr( pfake, p2d_rcv - p2d_src );
  csReflector2D line2d_sr( p2d_src, p2d_rcv - p2d_src );
  csReflector2D line2d_cmp( p2d_mid, vecNormRefl2d );
  csPoint2D p2d_cmp_atSurface = line2d_sr.intersection( line2d_cmp );
  double src2cmp_offsetRatio = abs(p2d_cmp_atSurface.x - p2d_src.x) / p2d_src.distance( p2d_rcv );
  // Remember that for equivalent 2D geometry, source & receiver positions may have been swapped
  // In the equivalent 2D geometry, the reflector is deepest below the 2D 'source' and dipping up to the 2D 'receiver'
  if( distReflSrc < distReflRcv ) {
    src2cmp_offsetRatio = 1.0 - src2cmp_offsetRatio;
  }
  // Compute CMP position on 3D reflector by projecting surface CMP position down to 3D reflector by going along the plane's normal vector
  csLine3D lineCMP;
  lineCMP.origin.x = pSrc.x +  src2cmp_offsetRatio * (pRcv.x - pSrc.x);
  lineCMP.origin.y = pSrc.y +  src2cmp_offsetRatio * (pRcv.y - pSrc.y);
  lineCMP.origin.z = pSrc.z +  src2cmp_offsetRatio * (pRcv.z - pSrc.z);
  lineCMP.direction = normalVector();

  rayPathParam param;
  param.distanceSrcCMP   = p2d_mid.distance( p2d_src );  // Ray path distance from source to receiver
  param.distanceRcvCMP   = p2d_mid.distance( p2d_rcv );  // Ray path distance from receiver to CMP
  //  param.distanceTotal    = param.distanceSrcCMP + param.distanceRcvCMP;  // Total ray path distance
  param.offRatioSrcCMP   = src2cmp_offsetRatio; // Normalized horizontal offset from source to CMP. 0: At source position, 1: At receiver position, 0.5: Mid-way between...
  param.cmp              = intersection( lineCMP );

  if( debugStream ) {
    fprintf(debugStream,"MISC_LINECMP_ORIG  %.2f %.2f %.2f\n",  lineCMP.origin.x,  lineCMP.origin.y,  lineCMP.origin.z);

    fprintf(debugStream,"RAYPATH2D_SR   %.2f %.2f\n", p2d_src.x, p2d_src.z);
    fprintf(debugStream,"RAYPATH2D_SR   %.2f %.2f\n", p2d_rcv.x, p2d_rcv.z);
    fprintf(debugStream,"RAYPATH2D_SR\n");

    fprintf(debugStream,"RAYPATH2D_SCR  %.2f %.2f\n", p2d_src.x, p2d_src.z);
    fprintf(debugStream,"RAYPATH2D_SCR  %.2f %.2f\n", p2d_mid.x, p2d_mid.z);
    fprintf(debugStream,"RAYPATH2D_SCR  %.2f %.2f\n", p2d_rcv.x, p2d_rcv.z);
    fprintf(debugStream,"RAYPATH2D\n");

    fprintf(debugStream,"RAYPATH2D_NORMAL %.2f %.2f\n", p2d_mid.x, p2d_mid.z);
    fprintf(debugStream,"RAYPATH2D_NORMAL %.2f %.2f\n", p2d_cmp_atSurface.x, p2d_cmp_atSurface.z);
    fprintf(debugStream,"RAYPATH2D_NORMAL\n");
  }

  return param;
}

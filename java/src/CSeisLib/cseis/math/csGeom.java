/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.math;

import cseis.general.csLinePosition;
import java.awt.geom.*;

/**
 * Provides geometrical methods
 * @author Bjorn Olofsson
 */
public class csGeom {
  public IntersectionProperties lineIntersection( csLinePosition lp1, csLinePosition lp2 ) {
    double bx = lp1.x2 - lp1.x1;
    double by = lp1.y2 - lp1.y1;
    double dx = lp2.x2 - lp2.x1;
    double dy = lp2.y2 - lp2.y1;
    
    double b_dot_d_perp = bx*dy - by*dx;
    
    if(b_dot_d_perp == 0) return null;
    
    double cx = lp2.x1-lp1.x1;
    double cy = lp2.x1-lp1.y1;
    
    IntersectionProperties prop = new csGeom.IntersectionProperties();
    prop.l1 = (cx*dy - cy*dx) / b_dot_d_perp;
    prop.x  = lp1.x1 + prop.l1*bx;
    prop.y  = lp1.y1 + prop.l1*by;
    if( Math.abs(dx) > Math.abs(dy) ) {
      prop.l2 = ( prop.x - lp2.x1 ) / dx;
    }
    else {
      prop.l2 = ( prop.y - lp2.y1 ) / dy;
    }
    return prop;
  }
  /**
   * Compute area of polygon
   * @param p Polygon corner points (non-crossing line segements)
   * @return Size of polygon area, in squared unit of input points
   */
  public static double polygonAreaSize( Point2D.Double[] p ) {
    int numPoints = p.length;
    double area = 0.0;
    for( int i = 0; i < numPoints-1; i++ ){
      area += p[i].x*p[i+1].y - p[i+1].x*p[i].y;
    }
    // Add last cross product for the case where polygon is not closed:
    // (If the last polygon point is the same as the first one, this cross product will be zero)
    area += p[numPoints-1].x*p[0].y - p[0].x*p[numPoints-1].y;
    return Math.abs(area/2.0);
  }
  public static boolean isPointInPolygon( Point2D.Double[] polygon, Point2D.Double point ) {
    int numPolPoints = polygon.length;
    boolean test = false;
    for( int i = 0, j = numPolPoints-1; i < numPolPoints; j = i++) {
      if( ((polygon[i].y > point.y) != ( polygon[j].y > point.y)) &&
        (point.x < (polygon[j].x-polygon[i].x) * (point.y-polygon[i].y) / (polygon[j].y - polygon[i].y) + polygon[i].x) )
      {
        test = !test;
      }
    }
    return test;
  }

  public class IntersectionProperties {
    public double l1;
    public double l2;
    public double x;
    public double y;
  }
  
  public static void main( String[] args ) {
    double xp[] = { 1, 7, 7, 7, 4, 1, -1, -1 };
    double yp[] = { 4, 1, 4, 7, 5, 7,  5,  2 };
    Point2D.Double[] p = new Point2D.Double[xp.length];
    for( int i = 0; i < xp.length; i++ ) {
      p[i] = new Point2D.Double( xp[i], yp[i] );
    }    
    double area = csGeom.polygonAreaSize( p );
    System.out.println("Area = " + area );
  }
}



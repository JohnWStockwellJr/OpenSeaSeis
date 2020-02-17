/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import java.awt.geom.Point2D;

/**
 * 2D line representation
 * @author Bjorn Olofsson
 */
public class csLine2D {
  private double grad_x;
  private double grad_y;
  private Point2D.Double myP1;
  private Point2D.Double myP2;
  private double myLength;
  
  public csLine2D() {
    myP1 = new Point2D.Double(0,0);
    myP2 = new Point2D.Double(0,0);
    init();
  }
  public csLine2D( csLine2D line_in ) {
    this( line_in.myP1, line_in.myP2 );
  }
  public csLine2D( Point2D.Double p1_in, double grad_x_in, double grad_y_in ) {
    myP1 = p1_in;
    grad_x = grad_x_in;
    grad_y = grad_y_in;
    myP2 = pos( 1.0 );
    myLength = Math.sqrt( grad_x*grad_x + grad_y*grad_y );
  }
  public csLine2D( Point2D.Double p1_in, Point2D.Double p2_in ) {
    myP1 = p1_in;
    myP2 = p2_in;
    init();
  }
  public void setPoints( Point2D.Double p1_in, Point2D.Double p2_in ) {
    myP1 = p1_in;
    myP2 = p2_in;
    init();
  }
  public void setLength( double length ) {
    if( myLength != 0.0 ) {
      double scalar = length/myLength;
      grad_x *= scalar;
      grad_y *= scalar;
      myP2 = new Point2D.Double( myP1.x + grad_x, myP1.y + grad_y );
    }
    else {
//      System.out.println("ERROR");
    }
    myLength = length;
  }
  public Point2D.Double p1() { return myP1; }
  public Point2D.Double p2() { return myP2; }
  public double grad_x() { return grad_x; }
  public double grad_y() { return grad_y; }

  // normDist: Normalised distance, between 0 and 1
  public Point2D.Double pos( double normDist ) {
    return new Point2D.Double( myP1.x + normDist*grad_x, myP1.y + normDist*grad_y );
  }

  private void init() {
    grad_x = myP2.x - myP1.x;
    grad_y = myP2.y - myP1.y;
    myLength = Math.sqrt( grad_x*grad_x + grad_y*grad_y );
  }
  
  /**
   * Compute intersection point between two lines
   * 
   * @param line2           (i) Line to intersect with
   * @param intersect       (o) Intersection point
   * @param angle_tolerance (i) Do not compute intersection if angle between lines is smaller than specified tolerance
   * @return false if no intersection point could be found (--> lines are sub-parallel)
   */
  public boolean intersection( csLine2D line2, Point2D.Double intersect, double angle_tolerance ) {
//    double length1 = myLength;
    double an_rad1 = ( 2*Math.PI + Math.atan2( grad_x, grad_y ) ) % Math.PI;

    double grad_x2 = line2.grad_x;
    double grad_y2 = line2.grad_y;
//    double length2 = Math.sqrt( grad_x2*grad_x2 + grad_y2*grad_y2 );
    double an_rad2 = ( 2*Math.PI + Math.atan2( grad_x2, grad_y2 ) ) % Math.PI;
    
    if( Math.abs( an_rad1 - an_rad2 )*180/Math.PI < angle_tolerance ) {
      return false;
    }

    double tmp1 = grad_x*grad_y2 - grad_y*grad_x2;
    if( Math.abs( tmp1 ) == 0.0 ) return false;

    double lamda2 = (grad_y*( line2.myP1.x - myP1.x ) + grad_x*( myP1.y - line2.myP1.y ) ) / tmp1;
    double lamda1;
    if( Math.abs(grad_x) > Math.abs(grad_y) ) {
      lamda1 = ( ( line2.myP1.x - myP1.x ) + lamda2*grad_x2 ) / grad_x;
    }
    else if( Math.abs(grad_y) != 0.0 ) {
      lamda1 = ( ( line2.myP1.y - myP1.y ) + lamda2*grad_y2 ) / grad_y;
    }
    else {
      return false;
    }
    
    intersect.x = myP1.x + grad_x * lamda1;
    intersect.y = myP1.y + grad_y * lamda1;

    return true;
  }
  /**
   * Compute distance between this line and point p
   * @param p (i) Point
   */
  public double distance( Point2D.Double p ) {
    Point2D.Double intersection = new Point2D.Double();
    return distance( p, intersection );
  }
  /**
   * Compute distance between this line and point p
   * @param p (i) Point
   * @param intersection (o) Intersection point of normal vector on line that passes through point p
   * @return
   */
  public double distance( Point2D.Double p, Point2D.Double intersection ) {
    double offset = ( (p.x - myP1.x) * (myP2.x - myP1.x) + (p.y - myP1.y) * (myP2.y - myP1.y) ) / ( myLength * myLength );
    
//    if( returnOn && (offset < 0.0 || offset > 1.0) ) {
//      return 0;   // closest point does not fall within the line segment
//    }

    intersection.x = myP1.x + offset * (myP2.x - myP1.x);
    intersection.y = myP1.y + offset * (myP2.y - myP1.y);
    double dx = intersection.x - p.x;
    double dy = intersection.y - p.y;
    double distance = Math.sqrt( dx*dx + dy*dy );

    return distance;
  }

}



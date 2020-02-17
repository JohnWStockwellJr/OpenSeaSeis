/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import java.awt.geom.Point2D;

/**
 * 2D line segment representation.<br>
 * The line segment is represented by its two corner points.
 * @author Bjorn Olofsson
 */
public class csLinePosition {
  public double x1;
  public double y1;
  public double x2;
  public double y2;
  public double length;  // derived from xy's

  public csLinePosition() {
    x1 = 0;
    x2   = 0;
    y1 = 0;
    x2   = 0;
    length = 0;
  }
  public csLinePosition( Point2D.Double p1, Point2D.Double p2 ) {
    x1 = p1.x;
    y1 = p1.y;
    x2   = p2.x;
    y2   = p2.y;
    double dx = x2 - x1;
    double dy = y2 - y1;
    length = Math.sqrt( dx*dx + dy*dy );
  }
  public csLinePosition( csLinePosition p ) {
    x1 = p.x1;
    x2   = p.x2;
    y1 = p.y1;
    y2   = p.y2;
    length = p.length;
  }
  public void computeLength() {
    double dx = x2 - x1;
    double dy = y2 - y1;
    length = Math.sqrt( dx*dx + dy*dy );
  }
  /**
   * Compute angle between two lines
   * @param line2
   * @return angle
   */
  public double angleTo( csLinePosition line2 ) {
    double grad_x1 = x2 - x1;
    double grad_y1 = y2 - y1;
    double an_rad1 = ( 2*Math.PI + Math.atan2( grad_x1, grad_y1 ) ) % Math.PI;

    double grad_x2 = line2.x2 - line2.x1;
    double grad_y2 = line2.y2 - line2.y1;
    double an_rad2 = ( 2*Math.PI + Math.atan2( grad_x2, grad_y2 ) ) % Math.PI;
    
    return( Math.abs( an_rad1 - an_rad2 )*180/Math.PI );
  }
  // Assumes that length has been computed beforehand
  public boolean intersection( csLinePosition p, Point2D.Double intersect, double dist_tolerance1, double dist_tolerance2, double angle_tolerance ) {
    double grad_x1 = x2 - x1;
    double grad_y1 = y2 - y1;
    double length1 = length;
    double an_rad1 = ( 2*Math.PI + Math.atan2( grad_x1, grad_y1 ) ) % Math.PI;

    double grad_x2 = p.x2 - p.x1;
    double grad_y2 = p.y2 - p.y1;
    double length2 = Math.sqrt( grad_x2*grad_x2 + grad_y2*grad_y2 );
    double an_rad2 = ( 2*Math.PI + Math.atan2( grad_x2, grad_y2 ) ) % Math.PI;
    
    if( Math.abs( an_rad1 - an_rad2 )*180/Math.PI < angle_tolerance ) {
      return false;
    }
    if( Math.abs( grad_x1 ) == 0.0 ) return false;
    double tmp1 = grad_y1 / grad_x1;
    double tmp2 = grad_y2 - grad_x2*tmp1;
    if( Math.abs( tmp2 ) == 0.0 ) return false;

    double lamda2 = (tmp1*( p.x1 - x1 ) + ( y1 - p.y1 ) ) / tmp2;
    double lamda1 = ( ( p.x1 - x1 ) + lamda2*grad_x2 ) / grad_x1;
    
    // Check tolerances
    if( lamda1 < 0 || lamda1 > 1.0 ) {
      double delta = lamda1 < 0 ? Math.abs( lamda1 ) : lamda1-1.0;
      if( delta*length1 > dist_tolerance1 ) return false;
    }
    if( lamda2 < 0 || lamda2 > 1.0 ) {
      double delta = lamda2 < 0 ? Math.abs( lamda2 ) : lamda2-1.0;
      if( delta*length2 > dist_tolerance2 ) return false;
    }
    intersect.x = x1 + grad_x1 * lamda1;
    intersect.y = y1 + grad_y1 * lamda1;
    
    return true;
  }
}



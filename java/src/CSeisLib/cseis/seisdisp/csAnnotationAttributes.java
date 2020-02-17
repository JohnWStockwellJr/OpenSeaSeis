/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

/**
 * Annotation attributes
 * @author Bjorn Olofsson
 */
public class csAnnotationAttributes {
  public boolean omitRepeating;
  public boolean showSequential;
  public boolean fixedTraceLabelStep;
  public int traceLabelStep;
  
  public csAnnotationAttributes() {
    omitRepeating       = true;
    showSequential      = true;
    fixedTraceLabelStep = false;
    traceLabelStep      = 10;
  }
}


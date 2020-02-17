/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

/**
 * Class defining a graphical curve.
 * @author Bjorn Olofsson
 */
public class csCurve {
  private static int staticCurveCounter = 1;
  public float[] dataX;
  public float[] dataY;
  public csCurveAttributes attr;
  public int curveID;

  public csCurve( float[] valuesX_in, float[] valuesY_in, csCurveAttributes attr_in ) {
    this(valuesX_in,valuesY_in,attr_in,csCurve.staticCurveCounter);
  }
  public csCurve( float[] valuesX_in, float[] valuesY_in, csCurveAttributes attr_in, int curveID_in ) {
    attr    = attr_in;
    dataX   = valuesX_in;
    dataY   = valuesY_in;
    curveID = curveID_in;
    csCurve.staticCurveCounter += 1;
  }

}



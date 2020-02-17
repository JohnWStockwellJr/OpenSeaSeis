/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.math;

/**
 * Provides static methods for mathematical data operations on pairs of data points (2D curves).
 * @author Bjorn Olofsson
 */
public class csDataOperation {
  private int myNumValuesOrig;
  private float[] myXValuesOrig;
  private float[] myYValuesOrig;

  private float[] myXValues;
  private float[] myYValues;
  
  public csDataOperation( float[] xValues, float[] yValues ) {
    myNumValuesOrig = xValues.length;
    myXValuesOrig = xValues;
    myYValuesOrig = yValues;
    myXValues = xValues;
    myYValues = yValues;
  }
  /**
   * Interpolate data points to regular grid, using mean values
   * @param increment Increment of new data points (X dimension)
   */
  public void grid_mean( float increment ) {
    float xValueStart = myXValuesOrig[0] - ( myXValuesOrig[0] % increment );
    int nValuesNew = (int)( (myXValuesOrig[myNumValuesOrig-1]-xValueStart)/increment + 0.5);
    if( nValuesNew >= myNumValuesOrig ) {
      myXValues = myXValuesOrig;
      myYValues = myYValuesOrig;
      return;
    }
    float[] xValuesNew = new float[nValuesNew];
    float[] yValuesNew = new float[nValuesNew];

    int indexOld   = 0;
    int indexStart = indexOld;
    for( int i = 0; i < nValuesNew; i++ ) {
      float xValue     = xValueStart + (float)i * increment;
      float xValueHalf = xValue + increment*0.5f;
      
      while( indexOld < myNumValuesOrig-1 && myXValuesOrig[indexOld] < xValueHalf ) {
        indexOld += 1;
      }

      float yValue = 0.0f;
      if( indexOld != indexStart ) {
        for( int index = indexStart; index < indexOld; index++ ) {
          yValue += myYValuesOrig[index];
        }
      }
      else {
        yValue = myYValuesOrig[indexOld];
      }
      
      xValuesNew[i] = xValue;
      yValuesNew[i] = yValue;
      indexStart = indexOld;
    }
    myXValues = xValuesNew;
    myYValues = yValuesNew;
  }
  public int numValues() {
    return myXValues.length;
  }
  public float[] xValues() {
    return myXValues;
  }
  public float[] yValues() {
    return myYValues;
  }
}



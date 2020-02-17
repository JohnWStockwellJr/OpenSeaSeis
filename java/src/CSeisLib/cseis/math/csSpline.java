/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.math;

public class csSpline {
  /**
   * Prepares 1D spline interpolation
   *
   * @param x float[]   Coordinates (need to be increasing monotonically, minimum 3 values)
   * @param y float[]   Input data values
   * @param yp1 float   First derivative at left edge  (if unsure, use 0.99e30f)
   * @param ypn float   First derivative at right edge (if unsure, use 0.99e30f)
   * @param y2 float[]  Calculated first derivatives
   */
  public static void spline( float x[], float y[], float yp1, float ypn, float y2[] ) {
    int n = x.length;
    if( n < 3 ) return;
    float[] uArray = new float[n-1];
    if( yp1 >= 0.99e30f ) {
      y2[0] = uArray[0] = 0.0f;
    }
    else {
      y2[0] = -0.5f;
      uArray[0] = ( 3.0f / ( x[1] - x[0] ) ) * ( ( y[1] - y[0] ) / ( x[1] - x[0] ) - yp1 );
    }
    for( int i = 1; i <= n - 2; i++ ) {
      float sig = ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] );
      float p = sig * y2[i-1] + 2.0f;
      y2[i] = ( sig - 1.0f ) / p;
      uArray[i] = ( y[i+1] - y[i] ) / ( x[i+1] - x[i] ) - ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
      uArray[i] = ( 6.0f * uArray[i] / ( x[i+1] - x[i-1] ) - sig * uArray[i-1] ) / p;
    }

    float qn, un;
    if( ypn >= 0.99e30f ) {
      qn = un = 0.0f;
    }
    else {
      qn = 0.5f;
      un = ( 3.0f / ( x[n-1] - x[n-2] ) ) * ( ypn - ( y[n-1] - y[n-2] ) / ( x[n-1] - x[n-2] ) );
    }
    y2[n-1] = ( un - qn * uArray[n-2] ) / ( qn * y2[n-2] + 1.0f );

    for( int k = n - 2; k >= 0; k-- ) {
      y2[k] = y2[k] * y2[k + 1] + uArray[k];
    }
  }
//-----------------------------------------------------------
  /**
   * 1D spline interpolation
   *
   * @param xa float[]   Coordinates (need to be increasing monotonically)
   * @param ya float[]   Input data values
   * @param y2a float[]  First derivatives as calculated by {@link spline}
   * @param x float      Coordinate of interpolated value
   * @return float       Interpolated value
   */
  public static float splint(float xa[], float ya[], float y2a[], float x) {
    int n = xa.length;
    float result;
    int klo,khi,k;
    float h,b,a;

    klo = 0;
    khi = n-1;
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (xa[k] > x) khi=k;
      else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) {
      new Exception("Badly sorted input array: Values are not increasing.").printStackTrace();
      return 0.0f;
    }
    a = (xa[khi]-x)/h;
    b = (x-xa[klo])/h;
    result = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0f;
    return result;
  }
  /**
   * 1D spline interpolation
   *
   * @param xa float[]   Coordinates (need to be increasing monotonically)
   * @param x_in float   Coordinate of interpolated value
   * @param attr         Further attributes   
   * @return float       Interpolated value
   */
  public static void splint_computeHelpAttr( float xa[], float[] x_in, csSplineHelpAttr[] attr ) {
    int numCoordinates = x_in.length;
    int numSeedPoints = xa.length;
    int klo,khi,k;
    float h,b,a;

    for( int icoord = 0; icoord < numCoordinates; icoord++ ) {
      float x1 = x_in[icoord];
      csSplineHelpAttr attr1 = attr[icoord];
      klo = 0;
      khi = numSeedPoints-1;
      while( khi-klo > 1 ) {
        k = (khi+klo) >> 1;
        if( xa[k] > x1 ) khi = k;
        else klo = k;
      }
      h = xa[khi] - xa[klo];
      a = (xa[khi]-x1)/h;
      b = (x1-xa[klo])/h;
      
      attr1.a1 = a;
      attr1.b1 = b;
      attr1.a3 = a*a*a-a;
      attr1.b3 = b*b*b-b;
      attr1.h2 = h*h/6.0f;
      attr1.khi = khi;
      attr1.klo = klo;
    }
//    result = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0f;
//    result = a1*ya[klo] + b1*ya[khi] + ( a3*y2a[klo] + b3*y2a[khi] ) * h2;
  }
//----------------------------------------------------------------
  /**
   * Prepares 2D spline interpolation
   *
   * @param x2a float[]   Dimension 2 coordinates (need to be increasing monotonically)
   * @param ya float[][]  Input data array
   * @param y2a float[][] Second derivates for dimension 2
   */
  public static void splie2(float x2a[], float[][] ya, float[][] y2a) {
    int m = ya.length;
    for( int j = 0; j < m; j++ ) {
      spline( x2a, ya[j], 1.0e30f, 1.0e30f, y2a[j] );
    }
  }
//----------------------------------------------------------------
  /**
   * 2D spline interpolation
   *
   * @param x1a float[]   Dimension 1 coordinates (need to be increasing monotonically)
   * @param x2a float[]   Dimension 2 coordinates (need to be increasing monotonically)
   * @param ya float[][]  Input data array
   * @param y2a float[][] Second derivates for dimension 2, as calculated by method {@link splie2}
   * @param x1 float      Coordinate (dim 1) of interpolate value
   * @param x2 float      Coordinate (dim 2) of interpolate value
   * @return float        Interpolated value at (x1,x2)
   */
  public static float splin2(float x1a[], float x2a[], float[][] ya, float[][] y2a,
      float x1, float x2)
  {
    int m = ya.length;

    float[] ytmp = new float[m];
    float[] yytmp = new float[m];

    for( int j = 0; j < m; j++ ) {
      yytmp[j] = splint( x2a, ya[j], y2a[j], x2 );
    }
    spline( x1a, yytmp, 1.0e30f, 1.0e30f, ytmp );
    float result = splint( x1a, yytmp, ytmp, x1 );
    return result;
  }
  /**
   * 2D spline interpolation, Preparation for MORE EFFICIENCY
   *
   * @param x1a float[]   Dimension 1 coordinates (need to be increasing monotonically)
   * @param x2a float[]   Dimension 2 coordinates (need to be increasing monotonically)
   * @param ya float[][]  Input data array
   * @param y2a float[][] Second derivates for dimension 2, as calculated by method {@link splie2}
   * @param x2 float      Coordinate (dim 2) of interpolate value
   * @param yytmp float[] Help field
   * @param ytmp float[]  Help field
   */
  public static void splin2_prepare( float x1a[], float x2a[], float[][] ya, float[][] y2a, float x2, float[] yytmp, float[] ytmp )
  {
    int m = ya.length;

    for( int j = 0; j < m; j++ ) {
      yytmp[j] = splint( x2a, ya[j], y2a[j], x2 );
    }
    spline( x1a, yytmp, 1.0e30f, 1.0e30f, ytmp );
  }
//----------------------------------------------------------------
  public static void main( String[] args ) {
/*
    float[] xValues = { 1, 11, 21, 31, 41, 51, 61, 71, 81, 91, 101 };
    float[] yValues = { 0,  0,  1,  1,  1,  0,  0,  1,  0,  1,   0 };
    float[] yDerivativeValues = new float[11];
    float yFirstDerivative = 0.0f;
    float yLastDerivative = 0.0f;
    Spline.spline(xValues,yValues,yFirstDerivative,yLastDerivative,yDerivativeValues);

    float[] yInterpolated = new float[101];
    for( int i = 1; i <= 101; i++ ) {
      yInterpolated[i-1] = Spline.splint(xValues,yValues,yDerivativeValues,(float)i);
      System.out.println( (i-1) + " " + yInterpolated[i-1]);
    }
*/

    float[] x1Values = { 1, 11, 21, 31, 41 };
    float[] x2Values = { 1, 11, 21, 31, 41 };
    float[][] yValues = {
        {0,  0,  0,  0,  0},
        {0,  0,  0,  1,  0},
        {0,  0,  1,  0,  0},
        {0,  1,  0,  0,  0},
        {0,  0,  0,  0,  0}
    };
    float[][] yDerivative2Values = new float[5][5];
    csSpline.splie2( x2Values, yValues, yDerivative2Values );

    float[][] yInterpolated = new float[41][41];
    for( int x1 = 1; x1 <= 41; x1++ ) {
//      System.out.print( x1 + " " );
      for( int x2 = 1; x2 <= 41; x2++ ) {
        yInterpolated[x1-1][x2-1] =
            csSpline.splin2( x1Values, x2Values, yValues, yDerivative2Values, (float)x1, (float)x2 );
        System.out.print( yInterpolated[x1-1][x2-1] + " " );
      }
      System.out.print( "\n" );
    }

  }
}



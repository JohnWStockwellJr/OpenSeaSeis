/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.image.*;
import java.awt.Color;
import cseis.math.csSpline;
import cseis.seis.csISeismicTraceBuffer;

/**
 * Renderer for seismic display.<br>
 * This class contains methods performing 2D rendering of seismic traces.<br>
 * <br>
 * Sub-classes are 32bit and 8bit realisations of this class, where color pixels are represented
 * by either 32bit integer values, and 8bit standard indexed bytes, respectively.
 * @author Bjorn Olofsson
 */
public abstract class csSeisViewRenderer implements csISeisViewListener {
  /// Full height of seismic view area
  protected int myImageHeight;
  /// Full width of seismic view area
  protected int myImageWidth;
  /// Seismic display settings
  protected csSeisDispSettings mySettings;
  protected csSeisView mySeisView;

  protected int myNumSamples;
  protected int myNumTraces;
  protected double mySampleInt;
  /// Derivatives of seismic amplitudes, required for spline interpolation
  protected float[][] myVA_yDerivative2Values = null;
  /// Flag indicating that derivatives are still being computed
  protected boolean myIsComputingDerivatives;
  /// Reference to seismic traces currently buffered in memory, containing all seismic trace samples
  protected csISeismicTraceBuffer myTraceBuffer;

  /// Constant scalar to apply to each trace. This is usually 1.
  /// When using 'trace' scaling, the trace scalar becomes the inverse of the mean or max of each individual trace.
  protected float[] myTraceScalar;
  /// The 'wiggle central amplitude' is usually 0; this is where the wiggle 'zeor' line is plotted.
  /// When using 'range' scaling, the wiggle central amplitude becomes the mean value of the specified min/max range values.
  protected float myWiggleCentreAmplitude;
  /// General display scalar. This local field is for convenience. The display scalar is also stored in the display settings (mySettings)
  protected float myDispScalar;
  /// Trace 'step'. Only every myStepTrace'th trace is drawn to the screen.
  /// This is to enhance performance when the display is squeezed below one pixel/trace.
  protected int myStepTrace = 1;
  
  protected float myBIAS = 0;

  public csSeisViewRenderer( csSeisView seisview, csSeisDispSettings ds ) {
    mySeisView = seisview;
    myTraceScalar    = null;
    
    myWiggleCentreAmplitude = 0.0f;
    myDispScalar = 1.0f;

    myVA_yDerivative2Values = null;
    myIsComputingDerivatives = false;
    mySettings = ds;
    changedSettings( ds );
  }
  public void resetTraceScaling( float dispScalar, float wiggleCentreAmplitude, float[] traceScalar ) {
    myDispScalar = dispScalar;
    myWiggleCentreAmplitude = wiggleCentreAmplitude;
    myTraceScalar = traceScalar;
  }
  public void reset( csISeismicTraceBuffer traceBuffer, double sampleInt ) {
    myTraceBuffer = traceBuffer;
    myNumSamples = myTraceBuffer.numSamples();
    myNumTraces  = myTraceBuffer.numTraces();
    mySampleInt   = sampleInt;

    myVA_yDerivative2Values = null;
    if( mySettings.isVIDisplay || mySettings.wiggleType == csSeisDispSettings.WIGGLE_TYPE_CUBIC ) {
      resetSplineFields();
    }
  }
  //-----------------------------------------------------
  protected void resetSplineFields( boolean doReset ) {
    if( myVA_yDerivative2Values == null && doReset ) resetSplineFields();
  }
  protected void resetSplineFields() {
    myIsComputingDerivatives = true;
    if( myVA_yDerivative2Values == null ) {
      myVA_yDerivative2Values = new float[myNumTraces][myNumSamples];
    }
    float[] va_ypos = new float[myNumSamples];
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      va_ypos[isamp] = isamp;
    }
    for( int traceIndex = 0; traceIndex < myNumTraces; traceIndex++ ) {
      float[] samples = myTraceBuffer.samples(traceIndex);
      csSpline.spline( va_ypos, samples, 1.0e30f, 1.0e30f, myVA_yDerivative2Values[traceIndex] );
    }
    myIsComputingDerivatives = false;
  }
  @Override
  public void traceBufferChanged( csISeismicTraceBuffer traceBuffer ) {
    // Nothing to be done
  }
  @Override
  public void changedSettings( csSeisDispSettings ds ) {
    boolean resetSpline = ( myVA_yDerivative2Values == null &&
        (ds.isVIDisplay != mySettings.isVIDisplay || ds.wiggleType != mySettings.wiggleType) &&
        (ds.isVIDisplay || ds.wiggleType == csSeisDispSettings.WIGGLE_TYPE_CUBIC) );
    
    mySettings = ds;
      
    if( resetSpline ) {
      resetSplineFields();
    }
  }
  public void mouseExited() {}
  protected abstract boolean resetScreenBuffer( int width, int height, DataBuffer dataBuffer );
  protected abstract void copyBufferHorz( int width, int height, int xfrom, int xto );
  protected abstract void copyBufferVert( int width, int height, int yfrom, int yto );
  
  /**
   * Repaint method, step 2
   * Clear or repaint background with VA display, then repaint wiggles on top
   * Finally, repaint time lines
   * 
   * @param yMin minimum y view (pixel) value to repaint
   * @param yMax maximum y view (pixel) value to repaint
   * @param xMin minimum x view (pixel) value to repaint
   * @param xMax maximum x view (pixel) value to repaint
   */
  protected abstract void repaintStep2( int yMin, int yMax, int xMin, int xMax );
  /**
   * Repaint time lines over given screen area
   * 
   * @param yMin minimum/maximum pixel values to repaint
   * @param yMax
   * @param xMin
   * @param xMax
   */
  protected abstract void repaintLines( int yMin, int yMax, int xMin, int xMax );
  /**
   * Repaint wiggle display, for given trace range and y pixel range
   * Basic concepts:
   *  a) Set pixel integer array. Do not use Java shapes, they're much too slow. Do not use bitmap, this is also too slow.
   *  b) Paint wiggles
   * 
   * @param yMin
   * @param yMax
   * @param minTrace
   * @param maxTrace
   * @param colorTrace
   */
  protected void repaintWiggleTrace( int yMin, int yMax, int minTrace, int maxTrace, Color colorTrace ) {
    repaintWiggleTrace( yMin, yMax, minTrace, maxTrace, colorTrace, true );
  }
  protected abstract void repaintWiggleTrace( int yMin, int yMax, int minTrace, int maxTrace, Color colorTrace, boolean doFill );
  /**
   * Repaint seismic display with discrete colours (nearest sample value tp screen pixel), without any interpolation
   * Variable density plot
   * 
   * @param yMin Minimum/maximum x pixel values to repaint
   * @param yMax
   * @param xMin
   * @param xMax
   */
  protected void repaintVADiscrete( int yMin, int yMax, int xMin, int xMax ) {
    repaintVADiscrete( yMin, yMax, xMin, xMax, false );
  }
  protected abstract void repaintVADiscrete( int yMin, int yMax, int xMin, int xMax, boolean doHighlight );
  /**
   * Repaint seismic display using linear or 1D spline interpolation, without interpolation in the horizontal direction
   * Variable density plot
   * 
   * @param yMin
   * @param yMax
   * @param xMin
   * @param xMax
   */
  protected abstract void repaintVAVertical( int yMin, int yMax, int xMin, int xMax );
  /**
   * Repaint seismic display using 2D spline interpolation
   * Variable density plot
   * 
   * @param yMin
   * @param yMax
   * @param xMin
   * @param xMax
   */
  protected abstract void repaintVA2DSpline( int yMin, int yMax, int xMin, int xMax );
  public void horzScrollChanged(int scrollValue) {
  }
  public void vertScrollChanged(int scrollValue) {
  }
  public void sizeChanged( java.awt.Dimension size ) {
  }
}



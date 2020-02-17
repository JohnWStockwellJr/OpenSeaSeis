/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.image.*;

import java.awt.Color;
import java.awt.image.DataBufferInt;
import cseis.math.csSpline;

/**
 * 32bit renderer for seismic display.<br>
 * Color pixels are represented by 32bit integer values.
 * @author Bjorn Olofsson
 */
public class csSeisViewRenderer32bit extends csSeisViewRenderer {
  private int[] myBuffer;
  private static int COLOR_BLACK_32BIT = 0;

  public csSeisViewRenderer32bit( csSeisView seisview, csSeisDispSettings ds ) {
    super( seisview, ds );
  }
  protected boolean resetScreenBuffer( int width, int height, DataBuffer dataBuffer ) {
    myImageWidth  = width;
    myImageHeight = height;
    myBuffer = ((DataBufferInt)dataBuffer).getData();
    return true;
  }
  protected void copyBufferHorz( int width, int height, int xfrom, int xto ) {
    if( xto < xfrom ) {
      for( int iy = 0; iy < myImageHeight; iy++ ) {
        int indexRow = iy*myImageWidth;
        for( int ix = 0; ix < width; ix++ ) {
          myBuffer[ix+xto+indexRow] = myBuffer[ix+xfrom+indexRow];
        }
      }
    }
    else {
      for( int iy = 0; iy < myImageHeight; iy++ ) {
        int indexRow = iy*myImageWidth;
        for( int ix = width-1; ix >= 0; ix-- ) {
          myBuffer[ix+xto+indexRow] = myBuffer[ix+xfrom+indexRow];
        }
      }
    }    
  }
  protected void copyBufferVert( int width, int height, int yfrom, int yto ) {
    if( yto < yfrom ) {
      for( int ix = 0; ix < myImageWidth; ix++ ) {
        for( int iy = 0; iy < height; iy++ ) {
          myBuffer[ix+myImageWidth*(iy+yto)] = myBuffer[ix+myImageWidth*(iy+yfrom)];
        }
      }
    }
    else {
      for( int ix = 0; ix < myImageWidth; ix++ ) {
        for( int iy = height-1; iy >= 0; iy-- ) {
          myBuffer[ix+myImageWidth*(iy+yto)] = myBuffer[ix+myImageWidth*(iy+yfrom)];
        }
      }
    }
  }

  
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
  protected void repaintStep2( int yMin, int yMax, int xMin, int xMax ) {
    if( myIsComputingDerivatives ) {
      return;
    }
    // System.out.println("Paint step 2  xmin/xmax  ymin/ymax " + xMin + " " + xMax + "   " + yMin + " " + yMax);
    // The following two statements catch the case when the painted area is smaller than the visible area
    // This may happen when the seisview window is smaller than the seispane window(=viewport)
    if( xMax > myImageWidth-1 )  xMax = myImageWidth-1;
    if( yMax > myImageHeight-1 ) yMax = myImageHeight-1;
    
    myStepTrace = 1;
    if( mySettings.zoomHorz < 1.0 ) {
      myStepTrace = (int)(1.0 / mySettings.zoomHorz);
    }

    int colorBackgroundRGB = Color.white.getRGB();
    // (1) Zero all pixels in given rectangle
    if( !mySettings.isVIDisplay ) {
      for( int iy = yMin; iy <= yMax; iy++ ) {
        int row = iy*myImageWidth;
        for( int ix = xMin; ix <= xMax; ix++ ) {
          myBuffer[ix+row] = colorBackgroundRGB;
        }
      }
    }
    // (1)  Alternatively to zeroing , paint VA display as background.
    else {
      // Only zero out left/right margin:
      for( int ix = 0; ix <= mySeisView.myMarginLeftRight-mySeisView.myViewPositionHorz; ix++ ) {
        for( int iy = yMin; iy <= yMax; iy++ ) {
          int row = iy*myImageWidth;
          myBuffer[ix+row] = colorBackgroundRGB;  // Exception in thread "AWT-EventQueue-0" java.lang.ArrayIndexOutOfBoundsException: 3408
        }
      }
      float pos = mySeisView.xModel2View( 0.0f, mySeisView.myNumTraces-1 );
      for( int ix = Math.max((int)pos, xMin); ix <= xMax; ix++ ) {
        for( int iy = yMin; iy <= yMax; iy++ ) {
          int row = iy*myImageWidth;
          myBuffer[ix+row] = colorBackgroundRGB;
        }
      }
      int xMinVA = Math.max( xMin, mySeisView.myMarginLeftRight-mySeisView.myViewPositionHorz );  // This is not correct yet!
      int xMaxVA = Math.min( xMax, (int)pos );
      switch( mySettings.viType ) {
        case csSeisDispSettings.VA_TYPE_DISCRETE:
          repaintVADiscrete( yMin, yMax, xMin, xMax );
          break;
        case csSeisDispSettings.VA_TYPE_VERTICAL:
          repaintVAVertical( yMin, yMax, xMin, xMax );
          break;
        case csSeisDispSettings.VA_TYPE_2DSPLINE:
          repaintVA2DSpline( yMin, yMax, xMinVA, xMaxVA );
          break;
      }
    }
    // (2) Paint wiggle display
    if( mySettings.showWiggle || mySettings.isNegFill || mySettings.isPosFill ) {
      // Determine min/max trace that contribute to the given rectangle
      int minTrace = (int)mySeisView.xView2Trace( xMin ) - (int)(mySettings.traceClip+0.5f) - 1;
      int maxTrace = (int)mySeisView.xView2Trace( xMax ) + (int)(mySettings.traceClip+0.5f) + 1;
      minTrace = Math.min( Math.max( minTrace, 0 ), myNumTraces-1 );
      maxTrace = Math.min( Math.max( maxTrace, 0 ), myNumTraces-1 );
  
      repaintWiggleTrace( yMin, yMax, minTrace, maxTrace, Color.black );
    }

    // (3) Paint time lines (or frequency or other domain...)
    if( mySettings.showTimeLines ) repaintLines( yMin, yMax, xMin, xMax );    
  }
  /**
   * Repaint time lines over given screen area
   * 
   * @param yMin minimum/maximum pixel values to repaint
   * @param yMax
   * @param xMin
   * @param xMax
   */
  protected void repaintLines( int yMin, int yMax, int xMin, int xMax ) {
    int sampMinorInc = (int)( mySettings.timeLineMinorInc / mySampleInt + 0.5f );
    int sampMajorInc = (int)( mySettings.timeLineMajorInc / mySampleInt + 0.5f );
    if( sampMinorInc <= 0 ) sampMinorInc = 1;
    if( sampMajorInc <= 0 ) sampMajorInc = 1;
    
    float sampleIndex = mySeisView.yView2Model( yMin );
    int minSamp = (int)sampleIndex;
    sampleIndex = mySeisView.yView2Model( yMax );
    int maxSamp = (int)sampleIndex + 1;
    
    int firstSamp = ( (int)( minSamp / sampMinorInc - 0.1f ) + 1 ) * sampMinorInc;
    
    for( int isamp = firstSamp; isamp <= maxSamp; isamp += sampMinorInc ) {
      int ypix = (int)mySeisView.yModel2View( isamp );
      if( ypix < 0 ) continue;
      int row = ypix*myImageWidth;
      if( ypix <= yMax ) {
        for( int ix = xMin; ix <= xMax; ix++ ) {
          myBuffer[ ix + row ] = COLOR_BLACK_32BIT;
        }
      }
      if( isamp % sampMajorInc == 0 ) {
        if( ypix < -1 || ypix > yMax-1 ) continue;
        row += myImageWidth;
        for( int ix = xMin; ix <= xMax; ix++ ) {
          myBuffer[ ix + row ] = COLOR_BLACK_32BIT;
        }
      }
    }
  }
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
  protected void repaintWiggleTrace( int yMin, int yMax, int minTrace, int maxTrace, Color colorTrace, boolean doFill ) {
    int colorTrace32bit = colorTrace.getRGB();
    //    long time1 = System.currentTimeMillis();   // Comment out to compute time of wiggle paint operation
//    for( int ii = 0; ii < 10; ii++ ) {

    // Trace step: Do not paint every single trace when display is squeezed to less than 0.5 pixels per trace.
    // Also, artificially boost samples amplitudes by trace step factor so that wiggles extend over more traces when display is squeezed.
    minTrace = (int)(minTrace/myStepTrace + 0.5f) * myStepTrace;
    maxTrace = (int)(maxTrace/myStepTrace + 0.5f) * myStepTrace;

    float traceClipPixels = mySeisView.xModel2View( myStepTrace*mySettings.traceClip ) - mySeisView.xModel2View( 0 );  // Boost trace clip by stepTrace for squeezed display
    int xAddum = 2000;  // Add 2000 to avoid xModel2View to become negative. This leads to problems when converting between float and integer

    if( yMin > 0 ) yMin -= 1;

    // Main loop: Repaint trace after trace. Always start with right-hand trace. This will bring positive wiggle to the front.
    for( int traceIndex = maxTrace; traceIndex >= minTrace; traceIndex -= myStepTrace ) {
      float[] samples = myTraceBuffer.samples(traceIndex);
      int xViewZeroCentre = (int)( mySeisView.xModel2View( myWiggleCentreAmplitude, traceIndex ) + xAddum + 0.5f );
      int xViewZero = (int)( mySeisView.xModel2View( traceIndex ) + xAddum + 0.5f );

      // Set some variables upfront, for 'previous' sample index
      int isamp1 = -1;
      int isamp2 = -1;
      float value1 = 0.0f;
      float value2 = 0.0f;
      float dx = 0.0f;

      int isampPrevious = -1;
      float sampleIndex = 0;
      float value   = 0;
      float valNext = 0;
      float valPrev = 0;
      float oneSixth = 1.0f/6.0f;

      // Loop from yMin-2: First two passes are only for setting 'valuePrev' and 'value'
      for( int yView = yMin-2; yView <= yMax; yView++ ) {
        valPrev = value;
        value   = valNext;

        sampleIndex = mySeisView.yView2Model( yView+1 );
        if( sampleIndex >= myNumSamples ) sampleIndex = myNumSamples-1.001f;
        if( sampleIndex <= 0 ) sampleIndex = 0.001f;
        isamp1 = (int)sampleIndex;
        if( isamp1 != isampPrevious && isamp1+1 < myNumSamples ) {
          isamp2 = ( isamp1+1 <= myNumSamples-1 ) ? isamp1+1 : isamp1;
          value1 = myStepTrace*samples[isamp1] + myBIAS;   // Boost sample values by stepTrace when display is squeezed a lot
          value2 = myStepTrace*samples[isamp2] + myBIAS;
          dx = value2 - value1;
        }
        isampPrevious = isamp1;

        if( mySettings.wiggleType == csSeisDispSettings.WIGGLE_TYPE_LINEAR ) {
          valNext = (sampleIndex-isamp1)*dx + value1;
        }
        else {  // mySettings.wiggleType == csSeisDispSettings.WIGGLE_TYPE_CUBIC ) {
          float a = isamp2-sampleIndex;
          float b = 1.0f - a;
          valNext = a*value1 + b*value2 + ( a*( a*a - 1 ) * myVA_yDerivative2Values[traceIndex][isamp1] + 
              b*( b*b - 1 ) * myVA_yDerivative2Values[traceIndex][isamp2] ) * oneSixth;
        }
        if( yView < yMin ) continue;  // First two passes are only for setting 'valuePrev' and 'value'

        int xView = (int)( mySeisView.xModel2View( mySettings.polarity*value, traceIndex ) + xAddum + 0.5f );

        int xStart = xView;
        int xEnd   = xView;

        // xView is the view x coordinate of the sample value at the current sample position
        // Determine which line shall be plotted to make a smooth line between adjacent samples
        int xMidPrev = (int)( mySeisView.xModel2View( mySettings.polarity*(valPrev+value)/2, traceIndex ) + xAddum );
        int xMidNext = (int)( mySeisView.xModel2View( mySettings.polarity*(valNext+value)/2, traceIndex ) + xAddum );
        if( xView > xMidPrev ) {
          if( xMidNext > xView ) {  // Wiggle curve is INCREASING both before and after current pixel
            xStart = xMidPrev + 1;
            xEnd   = xMidNext;
          }
          else {  // Wiggle curve exhibits MAXIMUM at current pixel
            xEnd   = xView;
            xStart = Math.min( xEnd, Math.min( xMidPrev+1, xMidNext+1 ) );
          }
        }
        else {
          if( xMidNext < xView ) {  // Wiggle curve is DECREASING both before and after current pixel
            xStart = xMidNext + 1;
            xEnd   = xMidPrev;
          }
          else {  // Wiggle curve exhibits MINIMUM at current pixel
            xStart = xView;
            xEnd   = Math.max( Math.max( xMidPrev, xMidNext ), xStart );
          }
        }

        if( mySettings.doTraceClipping ) {
          if( xViewZero - xStart > traceClipPixels ) {
            xStart = xViewZero - (int)traceClipPixels;
            xEnd = Math.max( xStart, xEnd );
          }
          if( xEnd - xViewZero > traceClipPixels ) {
            xEnd   = xViewZero + (int)traceClipPixels;
            xStart = Math.min( xStart, xEnd );
          }
        }

        // Fill wiggle with color
        if( (mySettings.isPosFill || mySettings.isNegFill) && doFill) {
          int colorVar = mySettings.wiggleColorMap.getColorRGB( mySettings.polarity*value );
          // Only fill if wiggle actually extends out from zero line (xStart > xViewZeroCentre)
          if( mySettings.isPosFill && xStart > xViewZeroCentre ) {
            int color = mySettings.wiggleColorPos.getRGB();
            if( mySettings.isVariableColor ) color = colorVar;
            for( int xpix = xViewZeroCentre-xAddum; xpix < xStart-xAddum; xpix++ ) {
              if( xpix < 0 || xpix >= myImageWidth ) continue;
              myBuffer[xpix + myImageWidth * yView ] = color;
            }
          }
          else if( mySettings.isNegFill && xStart < xViewZeroCentre ) {
            int color = mySettings.wiggleColorNeg.getRGB();
            if( mySettings.isVariableColor ) color = colorVar;
            for( int xpix = myStepTrace*xStart-xAddum+1; xpix <= xViewZeroCentre-xAddum; xpix++ ) {
              if( xpix < 0 || xpix >= myImageWidth ) continue;
              myBuffer[xpix + myImageWidth * yView ] = color;
            }
          }
        } // END: Fill wiggle with color
        // Paint wiggle trace
        if( mySettings.showWiggle ) {
          for( int xpix = xStart-xAddum; xpix <= xEnd-xAddum; xpix++ ) {
            if( xpix < 0 || xpix >= myImageWidth ) continue;
            myBuffer[xpix + myImageWidth*yView] = colorTrace32bit;
          }
        }
      }  // END Loop over all yView
      
      // Plot zero line, but only if zero line is not further away from trace than the trace clip value.
      // Note: If 'range' min/max values are specified, zero line may not coincide with actual trace zero amplitude.
      if( mySettings.showZeroLines && Math.abs( xViewZeroCentre-xViewZeroCentre ) < traceClipPixels ) {
        int xpix = xViewZeroCentre-xAddum;
        if( xpix >= 0 && xpix < myImageWidth ) {
          for( int yView = yMin; yView <= yMax; yView++ ) {
            myBuffer[xpix + myImageWidth*yView] = colorTrace32bit;
          }
        }
      }
    }
    
//    } // END ii
//    long time2 = System.currentTimeMillis();
//    System.out.println("NEW Elapsed time: " + (time2-time1) + "\n" );
    //      System.out.println("Number of paint operations " + this.myZoomHorz + " " + this.myZoomVert + "  --  " + "  min/max sample " + minSamp + " " + maxSamp );
  }

  protected void repaintVADiscrete( int yMin, int yMax, int xMin, int xMax, boolean doHighlight ) {
    // Determine min/max trace that contribute to the given rectangle
    int minTrace = (int)mySeisView.xView2Trace( xMin ) - (int)(mySettings.traceClip+0.5f) - 1;
    int maxTrace = (int)mySeisView.xView2Trace( xMax ) + (int)(mySettings.traceClip+0.5f) + 1;
    minTrace = Math.min( Math.max( minTrace, 0 ), myNumTraces-1 );
    maxTrace = Math.min( Math.max( maxTrace, 0 ), myNumTraces-1 );

    int xAddum = 2000;
    for( int traceIndex = maxTrace; traceIndex >= minTrace; traceIndex-- ) {
      float[] samples = myTraceBuffer.samples(traceIndex);
      int xmin = (int)( mySeisView.xModel2View( traceIndex-myStepTrace*0.5f ) + xAddum ) - xAddum;
      int xmax = (int)( mySeisView.xModel2View( traceIndex+myStepTrace*0.5f ) + xAddum ) - xAddum;
      if( xmin < xMin ) xmin = xMin;
      if( xmax > xMax ) xmax = xMax;
      for( int yView = yMin; yView <= yMax; yView++ ) {
        int sampleIndex = (int)( mySeisView.yView2Model( yView ) + 0.5f );
        if( sampleIndex < 0 ) sampleIndex = 0;
        if( sampleIndex > myNumSamples-1 ) sampleIndex = myNumSamples-1;
          // SCALE_TYPE_TRACE is not trivial. The following is a quick fix to make it work:
          if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE ) {
            mySettings.viColorMap.setScalar( myDispScalar*myTraceScalar[traceIndex] );
          }
        int colorRGB = mySettings.viColorMap.getColorRGB( mySettings.polarity*samples[sampleIndex] );
        if( doHighlight ) {
          colorRGB ^= 0x00aaaaff;
        }
        int row = myImageWidth*yView;
        for( int xpix = xmin; xpix <= xmax; xpix++ ) {
          myBuffer[xpix + row] = colorRGB;
        }
      }
    }
//    System.out.println("Given xmin/max: " + xMin + " " + xMax + "    ACTUAL xmin/xmax: " + minmin+ " " + maxmax + "  ---- trace: " + traceMinMin + " " + traceMaxMax );
  }
  /**
   * Repaint seismic display using linear or 1D spline interpolation, without interpolation in the horizontal direction
   * Variable density plot
   * 
   * @param yMin
   * @param yMax
   * @param xMin
   * @param xMax
   */
  protected void repaintVAVertical( int yMin, int yMax, int xMin, int xMax ) {
    // Determine min/max trace that contribute to the given rectangle
    int minTrace = (int)mySeisView.xView2Trace( xMin ) - (int)(mySettings.traceClip+0.5f) - 1;
    int maxTrace = (int)mySeisView.xView2Trace( xMax ) + (int)(mySettings.traceClip+0.5f) + 1;
    minTrace = Math.min( Math.max( minTrace, 0 ), myNumTraces-1 );
    maxTrace = Math.min( Math.max( maxTrace, 0 ), myNumTraces-1 );

    int xAddum = 2000;
    for( int traceIndex = maxTrace; traceIndex >= minTrace; traceIndex-- ) {
      float[] samples = myTraceBuffer.samples(traceIndex);
      int xmin = (int)( mySeisView.xModel2View( traceIndex-myStepTrace*0.5f ) + xAddum ) - xAddum;
      int xmax = (int)( mySeisView.xModel2View( traceIndex+myStepTrace*0.5f ) + xAddum ) - xAddum;
      if( xmin < xMin ) xmin = xMin;
      if( xmax > xMax ) xmax = xMax;
      
      for( int yView = yMin; yView <= yMax; yView++ ) {
        float sampleIndex = mySeisView.yView2Model( yView );
        if( sampleIndex >= myNumSamples-1 ) sampleIndex = myNumSamples-1.001f;
        if( sampleIndex <= 0 ) sampleIndex = 0.001f;
        int isamp1 = (int)sampleIndex;
        int isamp2 = isamp1 + 1;
        float value1 = samples[isamp1];
        float value2 = samples[isamp2];
        float dx = value2 - value1;
        float value;
        if( mySettings.wiggleType == csSeisDispSettings.WIGGLE_TYPE_LINEAR ) {
          value = (sampleIndex-isamp1)*dx + value1;
        }
        else {  // CUBIC
          float a = isamp2-sampleIndex;
          float b = 1.0f - a;
          value = a*value1 + b*value2 + ( a*( a*a - 1 ) * myVA_yDerivative2Values[traceIndex][isamp1] + 
              b*( b*b - 1 ) * myVA_yDerivative2Values[traceIndex][isamp2] ) * 1.0f/6.0f;
        }
        // SCALE_TYPE_TRACE is not trivial. The following is a quick fix to make it work:
        if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE ) {
          mySettings.viColorMap.setScalar( myDispScalar*myTraceScalar[traceIndex] );
        }
        int colorRGB = mySettings.viColorMap.getColorRGB( mySettings.polarity*value );
        int row = myImageWidth*yView;
        for( int xpix = xmin; xpix <= xmax; xpix++ ) {
          myBuffer[xpix + row] = colorRGB;
        }
      }
    }

  }
  /**
   * Repaint seismic display using 2D spline interpolation
   * Variable density plot
   * 
   * @param yMin
   * @param yMax
   * @param xMin
   * @param xMax
   */
  protected void repaintVA2DSpline( int yMin, int yMax, int xMin, int xMax ) {
    if( myNumTraces < 2 ) {
      repaintVADiscrete( yMin, yMax, xMin, xMax );
      return;
    }
    // Include 5 additional traces on each side to get correct result for spline interpolation
    // These limits are only required for computing interpolated spline values using csSpline.spline().
    // Repainting is only done within the given smaller limits
    int minTrace = (int)(mySeisView.xView2Trace( xMin ) - 5.2f);
    int maxTrace = Math.round(mySeisView.xView2Trace( xMax ) + 5.2f);

    minTrace = Math.min( Math.max( minTrace, 0 ), myNumTraces-1 );
    maxTrace = Math.min( Math.max( maxTrace, 0 ), myNumTraces-1 );
    if( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE ) {
      mySettings.viColorMap.setScalar( myDispScalar );
    }
    minTrace = Math.max( Math.round(minTrace/myStepTrace) * myStepTrace, 0);
    maxTrace = Math.min(Math.round(maxTrace/myStepTrace + 0.5f) * myStepTrace, myNumTraces-1);
    // Compute x view position for each trace that shall be painted, preset arrays for 2D interpolation
    int nTraces = (maxTrace-minTrace)/myStepTrace + 1;
    // Make sure there are least two traces to be repainted. One is too few for the 2D interpolation
    if( nTraces == 1 ) {
      if( minTrace == 0 ) {
        maxTrace = minTrace + myStepTrace;
        if( maxTrace >= myNumTraces ) {
          maxTrace = myNumTraces-1;
          minTrace = maxTrace - myStepTrace;
        }
      }
      else {
        minTrace = maxTrace - myStepTrace;
        if( minTrace < 0 ) {
          minTrace = 0;
          maxTrace = minTrace + myStepTrace;
        }
      }
      nTraces = 2;
    }
    float[] xViewAtTrace  = new float[nTraces];
    float[] seqTraceIndex = new float[nTraces];
    for( int traceIndex = maxTrace; traceIndex >= minTrace; traceIndex -= myStepTrace ) {
      int traceIndexReduced = (traceIndex - minTrace) / myStepTrace;
      float xViewZero = mySeisView.xModel2View( traceIndex );
      xViewAtTrace[traceIndexReduced]  = xViewZero;
      seqTraceIndex[traceIndexReduced] = traceIndex;
    }

    float[] valueAtTrace = new float[nTraces];   // Array holding computed values for all traces at one yView position
    float[] va2D_yDerivative2Values = new float[nTraces]; // Second derivatives for 2D spline interpolation
    // Slice through repainted area one by one row from smallest yView to largest yView
    for( int yView = yMin; yView <= yMax; yView++ ) {
      float sampleIndex = mySeisView.yView2Model( yView );
      if( sampleIndex >= myNumSamples-1 ) sampleIndex = myNumSamples-1.001f;
      if( sampleIndex <= 0 ) sampleIndex = 0.001f;
      int isamp1 = (int)sampleIndex;
      int isamp2 = isamp1 + 1;

      // At current yView position, compute 1D interpolated spline value for all traces
      for( int traceIndex = maxTrace; traceIndex >= minTrace; traceIndex -= myStepTrace ) {
        int traceIndexReduced = ( traceIndex - minTrace ) / myStepTrace;
        float[] samples = myTraceBuffer.samples(traceIndex);
//        float xViewZero = xViewAtTrace[traceIndexReduced];

        // Quick fix  to make it work...
        float value1 = samples[isamp1];
        float value2 = samples[isamp2];
        float a = isamp2-sampleIndex;
        float b = 1.0f - a;
        valueAtTrace[traceIndexReduced] = a*value1 + b*value2 + ( a*( a*a - 1 ) * myVA_yDerivative2Values[traceIndex][isamp1] + 
            b*( b*b - 1 ) * myVA_yDerivative2Values[traceIndex][isamp2] ) * 1.0f/6.0f;
      }
      // Compute spline derivatives for second dimension, across traces
      csSpline.spline( seqTraceIndex, valueAtTrace, 1.0e30f, 1.0e30f, va2D_yDerivative2Values );

      // At current yView position, compute 2D spline values for all x pixel positions
      int row = myImageWidth*yView;

      for( int xpix = xMin; xpix <= xMax; xpix++ ) {
        float traceIndexReduced = ( mySeisView.xView2Trace( xpix ) - minTrace ) / (float)myStepTrace;
//        float traceIndexReduced = xView2Trace( xpix ) - minTrace;
        if( traceIndexReduced >= nTraces-1 ) traceIndexReduced = nTraces-1.001f;
        if( traceIndexReduced <= 0 ) traceIndexReduced = 0.001f;
        int trace1 = (int)traceIndexReduced;
        int trace2 = trace1+1;
        float value1 = valueAtTrace[trace1];
        float value2 = valueAtTrace[trace2];
        float a = trace2-traceIndexReduced;
        float b = 1.0f - a;
        float value = a*value1 + b*value2 + ( a*( a*a - 1 ) * va2D_yDerivative2Values[trace1] + 
            b*( b*b - 1 ) * va2D_yDerivative2Values[trace2] ) * 1.0f/6.0f;
        int colorRGB = mySettings.viColorMap.getColorRGB( mySettings.polarity*value );
        myBuffer[xpix + row] = colorRGB;
      }
    }
  }
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.processing;

import cseis.math.csSpline;
import cseis.seis.csDataBuffer;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csSeismicData;
import java.awt.BorderLayout;
import java.awt.GridLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Processing module: Interpolation
 * @author 2013 Felipe Punto
 */
public class csProcessingInterpolation implements csIProcessing {
  public static final String NAME = "Interpolation";

  private int myInterpolationFactor;
//  private final float mySampleInt;
  private final JTextField myTextFactor;

  private final JPanel myParameterPanel;
  private boolean myIsActive;
  
  public csProcessingInterpolation( int interpolationFactor, float sampleInt ) {
//    mySampleInt = sampleInt;
    
    myIsActive = true;
    
    myInterpolationFactor = interpolationFactor;
    myParameterPanel = new JPanel( new GridLayout(1,1) );
    JPanel panel1 = new JPanel( new BorderLayout() );
    JLabel labelFactor = new JLabel("Interpolation factor (2,3,4): ");
    myTextFactor = new JTextField( "" + interpolationFactor );

    myTextFactor.setToolTipText("<html><i>Specify 2 to interpolate 2:1 traces</i></html>");
    panel1.add( labelFactor, BorderLayout.WEST );
    panel1.add( myTextFactor, BorderLayout.CENTER );
    myParameterPanel.add( panel1 );
  }
  @Override
  public boolean isActive() {
    return myIsActive;
  }
  @Override
  public void setActive( boolean doSet ) {
    myIsActive = doSet;
  }
  @Override
  public String retrieveParameters() {
    try {
      myInterpolationFactor = Integer.parseInt( myTextFactor.getText() );
      if( myInterpolationFactor < 2 || myInterpolationFactor > 4 ) myInterpolationFactor = 2;
    }
    catch( NumberFormatException e ) {
      return "Invalid entry for window length";
    }
    return null;
  }
  @Override
  public JPanel getParameterPanel() {
    return myParameterPanel;
  }
  @Override
  public String getName() {
    return csProcessingInterpolation.NAME;
  }
  @Override
  public void apply( csISeismicTraceBuffer traceBufferIn, csDataBuffer dataBufferOut ) {
    int numSamples = traceBufferIn.numSamples();
    int numTraces  = traceBufferIn.numTraces();
//    int numTracesOut = 2*numTraces - 1;
    int numTracesOut = numTraces;
    for( int traceIndex = 0; traceIndex < numTracesOut; traceIndex++ ) {
//    for( int traceIndex = 0; traceIndex < numTraces; traceIndex++ ) {
      float[] samplesOut = new float[numSamples];
//      for( int isampOut = 0; isampOut < numSamples; isampOut++ ) {
//        samplesOut[isampOut] = isampOut;
//      }
      dataBufferOut.addDataTrace( new csSeismicData(samplesOut) );
    }

    
    /// Derivatives of seismic amplitudes, required for spline interpolation
    float[][] va_yDerivative2Values = new float[numTraces][numSamples];
    float[] va_ypos = new float[numSamples];
    for( int isamp = 0; isamp < numSamples; isamp++ ) {
      va_ypos[isamp] = isamp;
    }
    for( int traceIndex = 0; traceIndex < numTraces; traceIndex++ ) {
      float[] samples = traceBufferIn.samples(traceIndex);
      csSpline.spline( va_ypos, samples, 1.0e30f, 1.0e30f, va_yDerivative2Values[traceIndex] );
    }

    int minTraceIn = 0;
    int maxTraceIn = numTraces-1;
    int stepTrace = 1;
    float[] va_xpos = new float[numTraces];
    for( int traceIndex = maxTraceIn; traceIndex >= minTraceIn; traceIndex -= stepTrace ) {
      int traceIndexReduced = (traceIndex - minTraceIn) / stepTrace;
      va_xpos[traceIndexReduced] = traceIndex;
    }

    float[] valueAtTrace = new float[numTraces];   // Array holding computed values for all traces at one yView position
    float[] va2D_yDerivative2Values = new float[numTraces]; // Second derivatives for 2D spline interpolation
//    int numSamplesOut = 2*numSamples - 1;

    // Slice through repainted area one by one row from smallest yView to largest yView
    for( int isampOut = 0; isampOut < numSamples-1; isampOut++ ) {
      float sampleIndexOut = isampOut * 1.0f;
      int isamp1 = isampOut; //(int)sampleIndexOut;
      int isamp2 = isamp1 + 1;
//      System.out.println("Sample " + isampOut);

      // At current sample position, compute 1D interpolated spline value for all traces
      for( int traceIndex = maxTraceIn; traceIndex >= minTraceIn; traceIndex -= stepTrace ) {
        int traceIndexReduced = ( traceIndex - minTraceIn ) / stepTrace;
        float[] samples = traceBufferIn.samples(traceIndex);
        float value1 = samples[isamp1];
        float value2 = samples[isamp2];
        float a = isamp2-sampleIndexOut;
        float b = 1.0f - a;
        valueAtTrace[traceIndexReduced] = a*value1 + b*value2 + ( a*( a*a - 1 ) * va_yDerivative2Values[traceIndex][isamp1] + 
            b*( b*b - 1 ) * va_yDerivative2Values[traceIndex][isamp2] ) * 1.0f/6.0f;
      }
      // Compute spline derivatives for second dimension, across traces
      csSpline.spline( va_xpos, valueAtTrace, 1.0e30f, 1.0e30f, va2D_yDerivative2Values );

      // Interpolate new traces
      for( int itrcOut = 0; itrcOut < numTracesOut; itrcOut++ ) {
        float traceOut = itrcOut * 0.5f;
        int traceIn1 = (int)traceOut;
        int traceIn2 = traceIn1 + 1;
        if( traceIn2 >= numTraces ) {
          traceIn2 -= 1;
          traceIn1 -= 1;
        }
//        System.out.println("Trace " + traceOut + "  " + trace2 + " / " + numTraces);
        float value1 = traceBufferIn.samples(traceIn1)[isamp1];
        float value2 = traceBufferIn.samples(traceIn2)[isamp1];
        float a = traceIn2-traceOut;
        float b = 1.0f - a;
        float value = a*value1 + b*value2 + ( a*( a*a - 1 ) * va2D_yDerivative2Values[traceIn1] + 
            b*( b*b - 1 ) * va2D_yDerivative2Values[traceIn2] ) * 1.0f/6.0f;
        dataBufferOut.getDataTrace(itrcOut).samples()[isamp1] = value;
      }
    } // END for isampOut
    System.out.println("Interpolation complete");
  } // END: apply()

}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.processing;

import cseis.seis.csDataBuffer;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csSeismicData;
import java.awt.BorderLayout;
import java.awt.GridLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Processing module: AGC
 * @author 2013 Felipe Punto
 */
public class csProcessingAGC implements csIProcessing {
  public static final String NAME = "AGC";

  private float mySampleInt;
  private int myWindowHalfLengthInSamples;
  private int myWindowIncInSamples;
  private JTextField myTextWindow;
  private JTextField myTextStep;

  private JPanel myParameterPanel;
  private boolean myIsActive;
  
  public csProcessingAGC( float fullWindowLengthInUnits, float windowIncInUnits, float sampleInt ) {
    mySampleInt = sampleInt;
    
    myWindowHalfLengthInSamples = 500;
    myWindowIncInSamples = 0;
    myIsActive = true;
    
    myParameterPanel = new JPanel( new GridLayout(2,1) );
    JPanel panel1 = new JPanel( new BorderLayout() );
    JPanel panel2 = new JPanel( new BorderLayout() );
    JLabel labelWindow    = new JLabel("Window length: ");
    myTextWindow = new JTextField( "" + fullWindowLengthInUnits );
    JLabel labelStep      = new JLabel("Window increment: ");
    myTextStep   = new JTextField( "" + windowIncInUnits );
    labelWindow.setToolTipText("Window length in units of trace");
    myTextWindow.setToolTipText("Window length in units of trace");
    labelStep.setToolTipText("<html>Increment in units of trace.<br><i>Specify 0 for sliding window</i></html>");
    myTextStep.setToolTipText("<html>Increment in units of trace.<br><i>Specify 0 for sliding window</i></html>");
    panel1.add( labelWindow, BorderLayout.WEST );
    panel1.add( myTextWindow, BorderLayout.CENTER );
    panel2.add( labelStep, BorderLayout.WEST );
    panel2.add( myTextStep, BorderLayout.CENTER );
    myParameterPanel.add( panel1 );
    myParameterPanel.add( panel2 );
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
      float fullWindowLengthInUnits = Float.parseFloat( myTextWindow.getText() );
      myWindowHalfLengthInSamples = (int)Math.round( fullWindowLengthInUnits / (2.0f*mySampleInt) );
      if( myWindowHalfLengthInSamples < 1 ) myWindowHalfLengthInSamples = 1;
    }
    catch( NumberFormatException e ) {
      return "Invalid entry for window length";
    }
    try {
      float windowIncInUnits = Float.parseFloat( myTextStep.getText() );
      myWindowIncInSamples = (int)Math.round( windowIncInUnits / (2.0f*mySampleInt) );
      if( myWindowIncInSamples < 1 ) myWindowIncInSamples = 0;
    }
    catch( NumberFormatException e ) {
      return "Invalid entry for window increment";
    }
    return null;
  }
  @Override
  public JPanel getParameterPanel() {
    return myParameterPanel;
  }
  @Override
  public String getName() {
    return csProcessingAGC.NAME;
  }
  @Override
  public void apply( csISeismicTraceBuffer traceBufferIn, csDataBuffer dataBufferOut ) {
    int numSamples = traceBufferIn.numSamples();
    for( int itrc = 0; itrc < traceBufferIn.numTraces(); itrc++ ) {
      float[] samples = traceBufferIn.samples( itrc );
      float[] samplesOut = new float[numSamples];
      for( int isampOut = 0; isampOut < numSamples; isampOut++ ) {
        int samp1 = Math.max( 0, isampOut-myWindowHalfLengthInSamples );
        int samp2 = Math.min( numSamples-1, isampOut+myWindowHalfLengthInSamples );
        int nsamp = samp2-samp1+1;
        if( nsamp < 1 ) continue;
        float rms = 0.0f;
        for( int isampIn = samp1; isampIn <= samp2; isampIn++ ) {
          rms += samples[isampIn] * samples[isampIn];
        }
        rms = (float)Math.sqrt( rms / (float)nsamp );
        if( rms > 1.0e-30f ) {
          samplesOut[isampOut] = samples[isampOut] / rms;
        }
      }
      dataBufferOut.addDataTrace( new csSeismicData(samplesOut) );
    }
  }

}


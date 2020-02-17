/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.processing;

import cseis.jni.csNativeFilter;
import cseis.seis.csDataBuffer;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csSeismicData;
import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import java.text.DecimalFormat;
import javax.swing.Box;
import javax.swing.JCheckBox;

/**
 * Processing module: Bandpass Filter
 * @author Felipe Punto
 */
public class csProcessingFilter implements csIProcessing {
  public static final String NAME = "BPFilter";

  private float mySampleInt;
  private int myNumSamples;

  private final JCheckBox myBoxIsLowPass;
  private final JCheckBox myBoxIsHighPass;
  
  private final JTextField myTextLPFreq;
  private final JTextField myTextLPSlope;
  private final JTextField myTextHPFreq;
  private final JTextField myTextHPSlope;

  private final JPanel myParameterPanel;
  private boolean myIsActive;
  private csNativeFilter myNativeFilter;
  
  public csProcessingFilter( float highPassFreq, float highPassSlope, float lowPassFreq, float lowPassSlope, int numSamples, float sampleInt ) {
    mySampleInt = sampleInt;
    myNumSamples = numSamples;

    myBoxIsLowPass  = new JCheckBox("Low-pass filter", true);
    myBoxIsHighPass = new JCheckBox("High-pass filter", true);
    
    myIsActive = true;
    DecimalFormat format = new DecimalFormat(".00");
    myParameterPanel = new JPanel( new BorderLayout() );
    JPanel panelInner = new JPanel( new GridLayout(3,3) );
    JLabel labelFreq  = new JLabel("Frequency [Hz]");
    JLabel labelSlope = new JLabel("Slope [dB/oct]");
    
    myTextLPFreq = new JTextField( "" + format.format(lowPassFreq) );
    myTextLPSlope = new JTextField( "" + format.format(lowPassSlope) );
    myTextHPFreq = new JTextField( "" + format.format(highPassFreq) );
    myTextHPSlope = new JTextField( "" + format.format(highPassSlope) );

    myTextHPFreq.setToolTipText("High-pass -3dB filter point [Hz]");
    myTextLPFreq.setToolTipText("Low-pass -3dB filter point [Hz]");
    myTextHPSlope.setToolTipText("High-pass slope [dB/oct]");
    myTextLPSlope.setToolTipText("Low-pass slope [dB/oct]");
    
    panelInner.add( Box.createHorizontalGlue() );
    panelInner.add( labelFreq );
    panelInner.add( labelSlope );

    panelInner.add( myBoxIsLowPass );
    panelInner.add( myTextLPFreq );
    panelInner.add( myTextLPSlope );

    panelInner.add( myBoxIsHighPass );
    panelInner.add( myTextHPFreq );
    panelInner.add( myTextHPSlope );
    myParameterPanel.add( panelInner );

    myBoxIsLowPass.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        boolean isSelected = ((JCheckBox)e.getSource()).isSelected();
        myTextLPFreq.setEnabled( isSelected );
        myTextLPSlope.setEnabled( isSelected );
      }
    });
    myBoxIsHighPass.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        boolean isSelected = ((JCheckBox)e.getSource()).isSelected();
        myTextHPFreq.setEnabled( isSelected );
        myTextHPSlope.setEnabled( isSelected );
      }
    });
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
    float lowPassFreq;
    float lowPassSlope;
    float highPassFreq;
    float highPassSlope;
    try {
      lowPassFreq = Float.parseFloat( myTextLPFreq.getText() );
      if( lowPassFreq < 0 ) lowPassFreq = 0;
    }
    catch( NumberFormatException e ) {
      return "Invalid entry for low-cut frequency";
    }
    try {
      lowPassSlope = Float.parseFloat( myTextLPSlope.getText() );
      if( lowPassSlope <= 0 ) lowPassSlope = 0.01f;
    }
    catch( NumberFormatException e ) {
      return "Invalid entry for low-cut slope";
    }
    try {
      highPassFreq = Float.parseFloat( myTextHPFreq.getText() );
      if( highPassFreq < 0 ) highPassFreq = 0;
//      if( myHighCutFreq > my ) myHighCutFreq = 0;
    }
    catch( NumberFormatException e ) {
      return "Invalid entry for high-cut frequency";
    }
    try {
      highPassSlope = Float.parseFloat( myTextHPSlope.getText() );
      if( highPassSlope < 0 ) highPassSlope = 0.01f;
    }
    catch( NumberFormatException e ) {
      return "Invalid entry for high-cut slope";
    }
    if( !myBoxIsLowPass.isSelected() ) lowPassFreq = -1;  // Indicate that no low pass filter shall be applied
    if( !myBoxIsHighPass.isSelected() ) highPassFreq = -1;  // Indicate that no high pass filter shall be applied
    myNativeFilter = new csNativeFilter( myNumSamples, mySampleInt, lowPassFreq, lowPassSlope, highPassFreq, highPassSlope );
    return null;
  }
  @Override
  public JPanel getParameterPanel() {
    return myParameterPanel;
  }
  @Override
  public String getName() {
    return csProcessingFilter.NAME;
  }
  @Override
  public void apply( csISeismicTraceBuffer traceBufferIn, csDataBuffer dataBufferOut ) {
    for( int itrc = 0; itrc < traceBufferIn.numTraces(); itrc++ ) {
      float[] samplesIn  = traceBufferIn.samples( itrc );
      float[] samplesOut = new float[myNumSamples];
      myNativeFilter.performFilter( samplesIn, samplesOut );
      dataBufferOut.addDataTrace( new csSeismicData(samplesOut) );
    }
  }

}


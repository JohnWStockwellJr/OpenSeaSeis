/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.processing;

import cseis.seis.csDataBuffer;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csSeismicData;
import java.awt.BorderLayout;
import java.awt.Font;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * Processing module: Test (only for testing purposes, not a real module)
 * @author 2013 Felipe Punto
 */
public class csProcessingTest implements csIProcessing {
  public static final String NAME = "Test";

  private JPanel myParameterPanel;
  private boolean myIsActive;
  
  public csProcessingTest() {
    myIsActive = false;
    myParameterPanel = new JPanel( new BorderLayout() );
    JLabel label = new JLabel("No parameters to set up", JLabel.CENTER);
    Font font = label.getFont();
    label.setFont( new Font(font.getName(), Font.ITALIC, font.getSize()) );
    myParameterPanel.add( label );
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
    return null;
  }
  @Override
  public JPanel getParameterPanel() {
    return myParameterPanel;
  }
  @Override
  public String getName() {
    return csProcessingTest.NAME;
  }
  @Override
  public void apply( csISeismicTraceBuffer traceBufferIn, csDataBuffer dataBufferOut ) {
    int numSamples = traceBufferIn.numSamples();
    for( int itrc = 0; itrc < traceBufferIn.numTraces(); itrc++ ) {
      float[] samples = traceBufferIn.samples( itrc );
      float[] samplesOut = new float[numSamples];
      for( int isampOut = 0; isampOut < numSamples; isampOut++ ) {
        samplesOut[isampOut] = (float)Math.abs( samples[isampOut] );
      }
      dataBufferOut.addDataTrace( new csSeismicData(samplesOut) );
    }
  }

}


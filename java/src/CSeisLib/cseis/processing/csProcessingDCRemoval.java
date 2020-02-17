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
 * Processing module: DC removal
 * @author 2013 Felipe Punto
 */
public class csProcessingDCRemoval implements csIProcessing {
  public static final String NAME = "DC removal";
  private boolean myIsActive;
  private JPanel myParameterPanel;
  public csProcessingDCRemoval() {
    myIsActive = false;
    myParameterPanel = new JPanel( new BorderLayout() );
    JLabel label = new JLabel("No parameters to set up", JLabel.CENTER);
    Font font = label.getFont();
    label.setFont( new Font(font.getName(), Font.ITALIC, font.getSize()) );
    myParameterPanel.add( label );
  }
  @Override
  public void apply( csISeismicTraceBuffer traceBufferIn, csDataBuffer dataBufferOut ) {
    int numSamples = traceBufferIn.numSamples();
    for( int itrc = 0; itrc < traceBufferIn.numTraces(); itrc++ ) {
      float[] samples = traceBufferIn.samples( itrc );
      float[] samplesOut = new float[numSamples];
      float dcValue = 0.0f;
      for( int isamp = 0; isamp < numSamples; isamp++ ) {
        dcValue += samples[isamp];
      }
      dcValue /= (float)numSamples;
      for( int isamp = 0; isamp < numSamples; isamp++ ) {
        samplesOut[isamp] = samples[isamp] - dcValue;
      }
      dataBufferOut.addDataTrace( new csSeismicData(samplesOut) );
    }
  }

  @Override
  public boolean isActive() {
    return myIsActive;
  }

  @Override
  public void setActive(boolean doSet) {
    myIsActive = doSet;
  }

  @Override
  public String getName() {
    return csProcessingDCRemoval.NAME;
  }

  @Override
  public JPanel getParameterPanel() {
    return myParameterPanel;
  }

  @Override
  public String retrieveParameters() {
    return null;
  }
}


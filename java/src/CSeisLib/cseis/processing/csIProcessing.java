/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.processing;

import cseis.seis.csDataBuffer;
import cseis.seis.csISeismicTraceBuffer;
import javax.swing.JPanel;

/**
 * Processing interface
 * @author 2013 Felipe Punto
 */
public interface csIProcessing {
  public boolean isActive();
  public void setActive( boolean doSet );
  public String getName();
  public JPanel getParameterPanel();
  public void apply( csISeismicTraceBuffer traceInputBuffer, csDataBuffer dataOutputBuffer );
  public String retrieveParameters();
}


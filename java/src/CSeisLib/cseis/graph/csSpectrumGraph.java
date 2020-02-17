/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

import javax.swing.JPanel;
import javax.swing.JLabel;
import java.awt.BorderLayout;
import java.awt.Dimension;
import cseis.jni.csFFTObject;

/**
 * Specialised 2D graph for displaying of amplitude spectra.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csSpectrumGraph extends JPanel implements csIGraph2DListener {
  private csGraph2D myGraph;
  private JLabel myStatusBar;

  public csSpectrumGraph() {
    this("Amplitude spectrum", "Frequency [Hz]", "Amplitude" );
  }
  public csSpectrumGraph( String title, String axisXLabel, String axisYLabel ) {
    super( new BorderLayout() );
    myGraph = new csGraph2D();
    myGraph.setPreferredSize( new Dimension(800,600) );
    myGraph.addGraph2DListener(this);
    myStatusBar = new JLabel( "---", JLabel.LEFT );
    csGraphAttributes graphAttr = myGraph.getGraphAttributes();
    graphAttr.title  = title;
    graphAttr.xLabel = axisXLabel;
    graphAttr.yLabel = axisYLabel;
//    myGraphAttr.autoScaleAxes = false;
//    myGraphAttr.axisScaleX = csGraphAttributes.AXIS_SCALE_LINEAR;
//    myGraphAttr.axisScaleY = csGraphAttributes.AXIS_SCALE_DB_NEG;
//    myDialogSettings.changedSettings( myGraphAttr );
    
    float[] values1 = new float[1];
    float[] values2 = new float[1];
    values1[0] = 0;
    values2[0] = 0;
    csCurveAttributes attr = new csCurveAttributes();
    attr.pointType = csCurveAttributes.POINT_TYPE_NONE;
    attr.lineColor = java.awt.Color.black;
    attr.lineSize = 1;
    attr.fillColor = new java.awt.Color(100,255,100);

    myGraph.addCurve( 0, values1, values2, attr );

    add( myGraph, BorderLayout.CENTER );
    add( myStatusBar, BorderLayout.SOUTH );
  }
  public void update( csFFTObject fftObject ) {
    int nValues = fftObject.amplitude.length;
    float[] valuesX = new float[nValues];
    float[] valuesY = new float[nValues];

    for( int ivalue = 0; ivalue < nValues; ivalue++ ) {
      valuesX[ivalue] = (float)(fftObject.freqInc*ivalue);
      valuesY[ivalue] = (float)fftObject.amplitude[ivalue];
    }

    myGraph.updateCurve( 0, valuesX, valuesY, true );
  }
  @Override
  public void graph2DValues(float xModel, float yModel) {
    myStatusBar.setText("Frequency: " + xModel + " Hz   Amplitude: " + yModel);
  }
  @Override
  public void graphChanged() {
  }
}



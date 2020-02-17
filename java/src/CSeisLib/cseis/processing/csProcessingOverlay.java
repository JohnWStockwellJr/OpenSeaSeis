/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.processing;

import cseis.seisdisp.csISeisOverlay;
import cseis.seisdisp.csSeisView;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

/**
 * Processing overlay.<br>
 * This class defines a label which shows the processing steps applied to the current data set
 * plotted on top of a seismic view (csSeisView).
 * This is done by implementing the standard interface csISeisOverlay.
 * @author 2013 Felipe Punto
 */
public class csProcessingOverlay implements csISeisOverlay {
  private String myName;
  public csProcessingOverlay( String name ) {
    myName = name;
  }
  public void setName( String name ) {
    myName = name;
  }
  public void addProcess( String name ) {
    if( myName == null || myName.length() == 0 ) {
      myName = name;
    }
    else {
      myName += " + " + name;
    }
  }
  public String getName() {
    return myName;
  }
  @Override
  public void draw( csSeisView seisview, Graphics2D g ) {
    int margin = 1;
    Color colorText = Color.WHITE;
    Color colorBackground = Color.darkGray;

    FontMetrics fm = g.getFontMetrics();
    Rectangle2D rect = fm.getStringBounds( myName, g );

    g.setColor( colorBackground );
    g.fillRect( 0, 0, (int)rect.getWidth()+margin*2, (int)rect.getHeight()+margin*2 );

    g.setColor( colorText );
    g.drawString( myName, margin, fm.getAscent()+margin );
  }
  
}


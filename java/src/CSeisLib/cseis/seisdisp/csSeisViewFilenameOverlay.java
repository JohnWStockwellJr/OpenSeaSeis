/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;

/**
 * Filename overlay.<br>
 * This class defines a label which shows the file name of the current data set
 * plotted on top of a seismic view (csSeisView).
 * This is done by implementing the standard interface csISeisOverlay.
 * @author 2013 Felipe Punto
 */
public class csSeisViewFilenameOverlay implements csISeisOverlay {
  private String myFilename;
  public csSeisViewFilenameOverlay( String filename ) {
    myFilename = filename;
  }
  public void setFilename( String filename ) {
    myFilename = filename;
  }
  public String getFilename() {
    return myFilename;
  }
  @Override
  public void draw( csSeisView seisview, Graphics2D g ) {
    Rectangle rectVisible = seisview.getVisibleRect();
    int margin = 1;
    Color colorText = Color.black;
    Color colorBackground = new Color( 255, 255, 255, 200 );

    FontMetrics fm = g.getFontMetrics();
    Rectangle2D rect = fm.getStringBounds( myFilename, g );
    int x = rectVisible.width - fm.stringWidth( myFilename ) - margin*2;
    
    g.setColor( colorBackground );
    g.fillRect( x, 0, (int)rect.getWidth()+margin*2, (int)rect.getHeight()+margin*2 );

    g.setColor( colorText );
    g.drawString( myFilename, x+margin, fm.getAscent()+margin );
  }
  
}


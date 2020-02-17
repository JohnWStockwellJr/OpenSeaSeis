/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.Icon;

/**
 * Square color icon, defined by a square size and a color.<br>
 * This icon can for example be used as the icon in a JLabel.
 * @author Bjorn Olofsson
 */
public class csSquareIcon implements Icon {
  private int mySize;
  private Color myColor;
  private csColorMap myColorMap;
  private boolean myIsSingleColor;

  public csSquareIcon( int size, Color color ) {
    myIsSingleColor = true;
    myColor  = color;
    mySize  = size;
  }
  public csSquareIcon( int size, csColorMap colorMap ) {
    myIsSingleColor = false;
    myColorMap  = colorMap;
    mySize = size;
  }
  public Color getColor() {
    return myColor;
  }
  public void setColor( Color color ) {
    myColor = color;
  }
  public csColorMap getColorMap() {
    return myColorMap;
  }
  @Override
  public void paintIcon( Component c, Graphics g, int x, int y ) {
    Graphics2D g2d = (Graphics2D) g.create();

    g2d.setColor( Color.black );
    g2d.drawRect( x, y, mySize-1, mySize-1 );

    if( myIsSingleColor ) {
      g2d.setColor(myColor);
      g2d.fillRect(x+1, y+1, mySize-2, mySize-2 );
    }
    else {
      float minValue = (float)myColorMap.getMinValue();
      float maxValue = (float)myColorMap.getMaxValue();
      float valueRange = maxValue - minValue;
      for( int i = 0; i < mySize-2; i++ ) {
        float value = maxValue - valueRange * (float)i / (float)(mySize-3);
        g2d.setColor( new Color(myColorMap.getColorRGB(value)) );
        g2d.drawLine( x+1, y+i+1, x+mySize-2, y+i+1 );
      }
    }

    g2d.dispose();
  }
  
  @Override
  public int getIconWidth() {
    return mySize;
  }
  
  @Override
  public int getIconHeight() {
    return mySize;
  }
}



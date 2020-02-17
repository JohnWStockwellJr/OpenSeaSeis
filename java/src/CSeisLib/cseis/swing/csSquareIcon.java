/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import javax.swing.Icon;

/**
 * A square button with a specific color.
 * @author Bjorn Olofsson
 */
public class csSquareIcon implements Icon {
  private int mySize;
  private Color myColor;

  public csSquareIcon( Color color, int size ) {
    myColor = color;
    mySize  = size;
  }
  public int getIconHeight() {
    return mySize;
  }
  public int getIconWidth() {
    return mySize;
  }
  public void paintIcon(Component c, Graphics g, int x, int y) {
    Graphics2D g2 = (Graphics2D)g;
    g2.setColor(myColor);
    g2.fill( new Rectangle2D.Float(x,y,(float)mySize,(float)mySize) );
  }
}



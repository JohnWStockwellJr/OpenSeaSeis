/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import javax.swing.Icon;
import javax.swing.JFrame;
import javax.swing.JRadioButton;

/**
 * Radio button displaying a colored circle.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csColorRadioButton extends JRadioButton {
  private RoundIcon myRoundIcon;

  public csColorRadioButton() {
    this("");
  }
  public csColorRadioButton( String text ) {
    this( text, Color.black );
  }
  public csColorRadioButton( Color color ) {
    this( "", color );
  }
  public csColorRadioButton( String text, Color color ) {
    super( text );
    myRoundIcon = new RoundIcon( color );
    setIcon( myRoundIcon );
    repaint();
  }
  public void changeColor( Color color ) {
    myRoundIcon.changeColor( color );
    setIcon( myRoundIcon );
    repaint();
  }
  class RoundIcon implements Icon {
    private Color myColor;
    public RoundIcon( Color color ) {
      myColor = color;
    }
    public void changeColor( Color color ) {
      myColor = color;
    }
    @Override
    public void paintIcon( Component comp, Graphics g, int x, int y ) {
      Graphics2D g2 = (Graphics2D)g;
      g2.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
      g2.setColor( myColor );
      g2.fillOval( x, y, getIconWidth(), getIconHeight() );
    }
    @Override
    public int getIconWidth() {
      return 10;
    }
    @Override
    public int getIconHeight() {
      return 10;
    }
  }
  public static void main( String[] args ) {
    JFrame frame = new JFrame();
    frame.getContentPane().add( new csColorRadioButton(), java.awt.BorderLayout.NORTH );
    frame.getContentPane().add( new csColorRadioButton( Color.yellow ), java.awt.BorderLayout.CENTER );
    frame.getContentPane().add( new csColorRadioButton( Color.red ), java.awt.BorderLayout.SOUTH );
    frame.setSize( new Dimension( 400, 400) );
    frame.setLocation( 400, 200 );
    frame.setVisible( true );
  }
}


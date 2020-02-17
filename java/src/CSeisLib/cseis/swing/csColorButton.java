/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.util.ArrayList;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.Color;
import java.awt.Component;

/**
 * Square colored button, extending JButton
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csColorButton extends JButton {
  public static final int APPEARANCE_CENTER = 0;
  public static final int APPEARANCE_RIGHT  = 1;
  public static final int APPEARANCE_ABOVE  = 2;
  private Color myColor;
  private Component myParentComponent;
  private ArrayList<csColorChangeListener> myListeners;
  private int myIdentNumber;
  private boolean myShowDialog;

  public csColorButton(Component parent, Color color) {
    this( parent, color, 0 );
  }
  public csColorButton(Component parent, Color color, int id ) {
    this( parent, color, id, true );
  }
  public csColorButton(Component parent, Color color, int id, boolean showChooserDialog ) {
    super();
    myIdentNumber = id;
    myColor = color;
    myParentComponent = parent;
    myShowDialog = showChooserDialog;
    setToolTipText("Select color");
    setMargin( new Insets(0,0,0,0) );
    updateIcon();

    setHorizontalAlignment (SwingConstants.CENTER);
    setHorizontalTextPosition (SwingConstants.CENTER);
    setVerticalAlignment (SwingConstants.CENTER);
    setVerticalTextPosition (SwingConstants.BOTTOM);
    setIconTextGap (0);
    
    if( myShowDialog ) {
      addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          csColorButton b = (csColorButton)e.getSource();
          Color color = JColorChooser.showDialog(b, "Select color", myColor);
          if( color != null && color != myColor) {
            b.setColor(color);
            myColor = color;
            updateIcon();
          }
          if( myListeners != null ) {
            fireColorChangeEvent( myColor );
          }
          myParentComponent.requestFocus();
        }
      });
    }
  }
  public int getID() {
    return myIdentNumber;
  }
  public void setID( int id ) {
    myIdentNumber = id;
  }
  public void setScreenAppearence( int appearenceFlag ) {
    
  }
  private void updateIcon() {
    setIcon(new csSquareIcon(myColor, 16));
  }
  public void addColorChangeListener( csColorChangeListener listener ) {
    if( myListeners == null ) {
      myListeners = new ArrayList<csColorChangeListener>(1);
    }
    if( !myListeners.contains(listener) ) {
      myListeners.add(listener);
    }
  }
  public void removeColorChangeListener( csColorChangeListener listener ) {
    if( myListeners == null ) {
      myListeners = new ArrayList<csColorChangeListener>(1);
    }
    if( myListeners.contains(listener) ) {
      myListeners.remove(listener);
    }
  }
  public void setColor(Color color) {
    myColor = color;
    updateIcon();
  }
  public Color getColor() {
    return myColor;
  }
  private void fireColorChangeEvent( Color color ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      ((csColorChangeListener)myListeners.get(i)).colorChanged( this, color );
    }
  }
}



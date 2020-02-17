/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.awt.Component;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JFrame;

/**
 * JFrame holding one csDockPane when undocked.<br>
 * 
 * When a csDockPane is undocked, a new csDockFrame is opened with a csDockPane inside.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csDockFrame extends JFrame implements csIDockPaneListener {
  private csDockPane myPane;
  private boolean myIsEOF;
  public csDockFrame( csDockPane pane, Component component ) {
    super( pane.attr().title );
    myPane = pane;
    myIsEOF = false;
    getContentPane().add(pane);
    pack();
    setLocationRelativeTo(component);
    myPane.addPaneListener(this);
    addWindowListener(new WindowAdapter() {
      @Override
      public void windowClosing(WindowEvent e) {
        if( !myIsEOF ) {
          closeFrame();
          myPane.changeDockingState();
          myPane.hidePane();
          myIsEOF = true;
        }
      }
    });
  }
  private void closeFrame() {
    myPane.removePaneListener(this);
    myIsEOF = true;
    setVisible(false);
    dispose();
  }
  @Override
  public void hidePane( csDockPane pane ) {
    closeFrame();
  }
  @Override
  public void dockPane( csDockPane pane, boolean dock ) {
    if( dock ) closeFrame();
  }
  @Override
  public void scrollbars(boolean showScrollbars) {
  }
  @Override
  public void closePane(csDockPane pane) {
  }
  @Override
  public void maximizePane(csDockPane pane) {
  }
}


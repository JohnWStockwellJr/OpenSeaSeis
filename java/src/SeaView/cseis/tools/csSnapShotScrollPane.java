/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.tools;

import javax.swing.*;
import java.awt.*;
import java.awt.image.*;

/**
 * Scroll pane for snapshot thumbnails.<br>
 * Contains csSnapShotPanel.
 * @author 2011 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csSnapShotScrollPane extends JScrollPane {
  private int myTotalWidth  = 120;
  private int myTotalHeight = 120;
  JScrollPane myScrollPane;
  private csSnapShotPanel myPanel;

  public csSnapShotScrollPane() {
    super(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    int scrollBarWidth = super.getPreferredSize().width;
    setPreferredSize(new Dimension(myTotalWidth+scrollBarWidth,myTotalHeight));

    myPanel = new csSnapShotPanel();
    super.setViewportView(myPanel);
    addMouseListener(myPanel);
    addMouseMotionListener(myPanel);
  }
  protected csSnapShotPanel getPanel() {
    return myPanel;
  }
  public void addImage( BufferedImage image, String infoText ) {
    myPanel.addImage(image,infoText);
  }
  public void removeImage( int index ) {
    myPanel.removeImage(index);
  }
  public void addSnapShotListener( csISnapShotListener listener ) {
    myPanel.addSnapShotListener(listener);
  }
}



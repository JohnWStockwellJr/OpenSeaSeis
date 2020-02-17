/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.tools;

import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * Frame that captures display snapshots created by user.<br>
 * Snapshots can be browsed, moved and deleted.
 * @author 2011 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csSnapShotFrame extends JDialog implements csISnapShotListener, MouseListener, KeyListener {
  private csSnapShotScrollPane mySnapShotScrollPane;
  private JPanel myMainPanel;
  private int myMainImageIndex;
  private ArrayList<BufferedImage> mySnapShotImageList;
  private ArrayList<String> mySnapShotInfoList;
  boolean myIsMouseInWindow = false;

  public csSnapShotFrame( JFrame parentFrame, ImageIcon appIcon ) {
    super(parentFrame,"Snapshot Viewer");
    setModal(false);

    mySnapShotImageList = new ArrayList<BufferedImage>();
    mySnapShotInfoList  = new ArrayList<String>();
    myMainImageIndex = -1;

    mySnapShotScrollPane = new csSnapShotScrollPane();
    mySnapShotScrollPane.addSnapShotListener( this );
    if( appIcon != null ) setIconImage( appIcon.getImage() );

    addMouseListener(this);
    KeyboardFocusManager manager = KeyboardFocusManager.getCurrentKeyboardFocusManager();
    manager.addKeyEventDispatcher(new SnapShotKeyDispatcher(this));

    myMainPanel = new JPanel() {
      public void paintComponent( Graphics g ) {
        if( myMainImageIndex >= 0 && myMainImageIndex < mySnapShotImageList.size() ) {
          g.drawImage(mySnapShotImageList.get(myMainImageIndex), 0, 0, null);
        }
      }
    };
    getContentPane().add(myMainPanel, BorderLayout.CENTER);
    getContentPane().add(mySnapShotScrollPane,BorderLayout.EAST);

    pack();
    setSize( new Dimension(700,800) );
  }
  public void addImage( BufferedImage image, String infoText ) {
    mySnapShotImageList.add(image);
    mySnapShotInfoList.add(infoText);
    mySnapShotScrollPane.addImage(image, infoText);
    refreshFrameInfoText();
  }
  public void refreshFrameInfoText() {
    if( myMainImageIndex >= 0 && mySnapShotInfoList.size() > 0 ) {
      setTitle( "Image #" + (myMainImageIndex+1) + ": " + mySnapShotInfoList.get(myMainImageIndex).toString());
    }
    else {
      setTitle( "Snapshot Viewer");
    }
  }
  @Override
  public void moveSnapShot(int fromItemIndex, int toItemIndex) {
    if( fromItemIndex != toItemIndex ) {
      toItemIndex = Math.min( Math.max(toItemIndex, 0), mySnapShotImageList.size() );
      if( fromItemIndex >= 0 && fromItemIndex < mySnapShotImageList.size() ) {
        if( toItemIndex > fromItemIndex ) toItemIndex -= 1;
        BufferedImage img = mySnapShotImageList.remove(fromItemIndex);
        mySnapShotImageList.add(toItemIndex,img);
        String infoText = mySnapShotInfoList.remove(fromItemIndex);
        mySnapShotInfoList.add(toItemIndex,infoText);
      }
      refreshFrameInfoText();
    }
  }
  @Override
  public void selectSnapShot(int itemIndex, boolean forceSelection) {
    int indexSave = myMainImageIndex;
    if( itemIndex >= 0 && itemIndex < mySnapShotImageList.size() ) {
      myMainImageIndex = itemIndex;
    }
    else {
      myMainImageIndex = -1;
    }
    if( indexSave != myMainImageIndex || forceSelection ) {
      repaint();
      myMainPanel.invalidate();
      myMainPanel.repaint();
    }
    refreshFrameInfoText();
  }

  @Override
  public void deleteSnapShot(int itemIndex) {
    if( itemIndex >= 0 && itemIndex < mySnapShotImageList.size() ) {
      mySnapShotImageList.remove(itemIndex);
      mySnapShotInfoList.remove(itemIndex);
      mySnapShotScrollPane.removeImage(itemIndex);
      myMainImageIndex = Math.min( itemIndex, mySnapShotImageList.size()-1 );
      refreshFrameInfoText();
    }
  }
  @Override
  public void keyPressed(KeyEvent event) {
    if( myIsMouseInWindow || mySnapShotScrollPane.getPanel().isMouseInWindow() ) {
      mySnapShotScrollPane.getPanel().keyPressed(event);
    }
  }
  @Override
  public void keyTyped(KeyEvent e) {
  }
  @Override
  public void keyReleased(KeyEvent e) {
  }

  @Override
  public void mouseClicked(MouseEvent e) {
  }
  @Override
  public void mousePressed(MouseEvent e) {
  }
  @Override
  public void mouseReleased(MouseEvent e) {
  }
  @Override
  public void mouseEntered(MouseEvent e) {
    myIsMouseInWindow = true;
  }
  @Override
  public void mouseExited(MouseEvent e) {
    myIsMouseInWindow = false;
  }
  //***********************************************************************************
  //***********************************************************************************
  //***********************************************************************************
  private class SnapShotKeyDispatcher implements KeyEventDispatcher {
    private KeyListener myListener;
    public SnapShotKeyDispatcher( KeyListener listener ) {
      super();
      myListener = listener;
    }
    @Override
    public boolean dispatchKeyEvent(KeyEvent e) {
      if( e.getID() == KeyEvent.KEY_PRESSED) {
        myListener.keyPressed( e );
      }
      else if( e.getID() == KeyEvent.KEY_RELEASED ) {
        myListener.keyReleased( e );
      }
      else if( e.getID() == KeyEvent.KEY_TYPED ) {
        myListener.keyTyped( e );
      }
      return false;
    }
  }
}



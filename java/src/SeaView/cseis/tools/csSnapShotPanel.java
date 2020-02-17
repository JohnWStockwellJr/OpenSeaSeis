/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


/*
 * SeaSeis - Seismic processing system for seabed (OBS) data.
 * CSeisLib - Seismic display library
 * SeaView  - Seismic viewer prototype
 * Copyright (c) 2012, Bjorn Olofsson.
 * All rights reserved.
 */


package cseis.tools;

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;

/**
 * Panel containing snapshot thumbnails
 * Thumbnails can be selected, dragged and deleted
 */
@SuppressWarnings("serial")
public class csSnapShotPanel extends JPanel implements Scrollable,
        MouseListener, MouseMotionListener, KeyListener {
  private int myMaxUnitIncrement = 1;
  private int myTotalWidth  = 120;
  private int myTotalHeight = 120;
  private int myMargin = 5;
  private int myImageWidth  = myTotalWidth - 2*myMargin;
  private int myImageHeight = myTotalHeight - 2*myMargin;
  private ArrayList<Image> mySnapShotThumbList;
  private boolean myIsMouseInWindow = false;
  private int myMousePosY;
  private int myMousePressedPosY;
  private int myMouseDraggedIndex;
  private int myMousePressedIndex;
  private ArrayList<csISnapShotListener> myListeners;
  private int mySelectedIndex;
  private boolean myIsMouseDragged = false;

  public csSnapShotPanel() {
    myListeners = new ArrayList<csISnapShotListener>();
    mySnapShotThumbList = new ArrayList<Image>();
    mySelectedIndex = -1;
    setOpaque(true);
    setBackground(Color.white);
    myMaxUnitIncrement = 1;
    setAutoscrolls(true); // Enable drag events
    addMouseMotionListener(this);
    addMouseListener(this);
  }
  public void addSnapShotListener( csISnapShotListener listener ) {
    myListeners.add(listener);
  }
  public void paintComponent(Graphics g) {
    Graphics2D g2 = (Graphics2D)g;
    Rectangle rectVisible = getVisibleRect();
    g.setColor(Color.white);
    g.clearRect(rectVisible.x, rectVisible.y, rectVisible.width, rectVisible.height);
    g2.setStroke( new BasicStroke(1.0f));

    if( mySelectedIndex >= 0 && mySelectedIndex < mySnapShotThumbList.size() ) {
      int ypos = myTotalHeight * mySelectedIndex;
      g2.setColor(Color.orange);
      g2.setStroke( new BasicStroke(3.0f));
      g2.drawRect(myMargin/2, ypos+myMargin/2, myTotalWidth-myMargin, myTotalHeight-myMargin );
      g2.setStroke( new BasicStroke(1.0f));
    }
    for( int i = 0; i < mySnapShotThumbList.size(); i++ ) {
      int ypos = myTotalHeight * i;
      g.setColor(Color.black);
      if( !myIsMouseDragged || i != myMousePressedIndex ) {
        g.drawImage(mySnapShotThumbList.get(i), myMargin, ypos+myMargin, null);
        g.drawRect(myMargin, ypos+myMargin, myTotalWidth-2*myMargin, myTotalHeight-2*myMargin );
        if( !myIsMouseDragged && myIsMouseInWindow && myMousePosY >= ypos && myMousePosY < ypos+myTotalHeight ) {
          g.setColor(Color.red);
          g.drawRect(myMargin/2, ypos+myMargin/2, myTotalWidth-myMargin, myTotalHeight-myMargin );
        }
        g.drawString("Image #"+(i+1), myMargin, ypos+myTotalHeight-myMargin);
      }
    }
    if( myIsMouseDragged && myMousePressedIndex >= 0 && myMousePressedIndex < mySnapShotThumbList.size()) {
      int mousePos = myMousePosY-myMousePressedPosY;
      myMouseDraggedIndex = (int)((mousePos+myTotalHeight)/myTotalHeight);
      int yposNew  = myMouseDraggedIndex*myTotalHeight;
      g.setColor(Color.orange);
      g.drawRect(myMargin/2, yposNew-2, myTotalWidth-myMargin, 5 );

      // Draw image with 50% transparency
      g.setColor(Color.black);
      Image inputImage = mySnapShotThumbList.get(myMousePressedIndex);
      BufferedImage tempImage = new BufferedImage(inputImage.getWidth(null), inputImage.getHeight(null), BufferedImage.TYPE_INT_ARGB);
      Graphics gb = tempImage.getGraphics();
      gb.drawImage(inputImage, 0, 0, null);
      // Create a rescale filter op that makes the image 50% opaque
      float[] scales = { 1f, 1f, 1f, 0.5f };
      float[] offsets = new float[4];
      RescaleOp rop = new RescaleOp(scales, offsets, null);
      // Draw the image, applying the filter
      g2.drawImage(tempImage, rop, myMargin, mousePos+myMargin);

      g.drawRect(myMargin, mousePos+myMargin, myTotalWidth-2*myMargin, myTotalHeight-2*myMargin );
      g.drawString("Image #"+(myMouseDraggedIndex+1), myMargin, mousePos+myTotalHeight-myMargin);
    }
  }
  public void addImage( BufferedImage image, String infoText ) {
    Image thumb = image.getScaledInstance(myImageWidth, myImageHeight, BufferedImage.SCALE_SMOOTH);
    mySnapShotThumbList.add(thumb);
    mySelectedIndex = mySnapShotThumbList.size()-1;
    invalidate();
    repaint();
    fireSelectSnapShotEvent( mySnapShotThumbList.size()-1, false );
  }
  public void removeImage( int index ) {
    mySnapShotThumbList.remove(index);
    if( mySelectedIndex > index ) {
      mySelectedIndex -= index;
    }
    else {
      mySelectedIndex = Math.min(index, mySnapShotThumbList.size()-1);
    }
    invalidate();
    repaint();
    fireSelectSnapShotEvent( mySelectedIndex, true );
  }
  public boolean isMouseInWindow() {
    return myIsMouseInWindow;
  }
  public void mouseMoved(MouseEvent e) {
    myMousePosY = e.getY();
    repaint();
  }
  public void mouseDragged(MouseEvent e) {
    myMousePosY = e.getY();
    myIsMouseDragged = true;
//      Rectangle r = new Rectangle(e.getX(), e.getY(), 1, 1);
//      scrollRectToVisible(r);
    repaint();
  }
  public void mouseClicked(MouseEvent e) {
    if( SwingUtilities.isLeftMouseButton(e) ) {
      myMousePosY = e.getY();
      myIsMouseDragged = false;
      mySelectedIndex = (int)(myMousePosY/myTotalHeight);
      fireSelectSnapShotEvent(mySelectedIndex,false);
      repaint();
    }
  }

  @Override
  public void mousePressed(MouseEvent e) {
    if( SwingUtilities.isLeftMouseButton(e) ) {
      myMousePosY = e.getY();
      myMousePressedIndex = (int)(myMousePosY/myTotalHeight);
      myMousePressedPosY = myMousePosY - myMousePressedIndex*myTotalHeight;
    }
  }

  @Override
  public void mouseReleased(MouseEvent e) {
    if( myIsMouseDragged && myMousePressedIndex >= 0 && myMousePressedIndex != myMouseDraggedIndex ) {
      myMouseDraggedIndex = Math.min( Math.max(myMouseDraggedIndex, 0), mySnapShotThumbList.size() );
      int dyMoved    = myMouseDraggedIndex - myMousePressedIndex;
      int dyDragged  = myMouseDraggedIndex - mySelectedIndex;
      int dyPressed  = myMousePressedIndex - mySelectedIndex;
      if( dyPressed == 0 ) {
        if( dyMoved > 0 ) mySelectedIndex = myMouseDraggedIndex-1;
        else mySelectedIndex = myMouseDraggedIndex;
      }
      else if(dyPressed < 0) mySelectedIndex -= 1;
      else if(dyDragged <= 0) mySelectedIndex += 1;
      fireMoveSnapShotEvent( myMousePressedIndex, myMouseDraggedIndex );
      fireSelectSnapShotEvent( mySelectedIndex, true );

      if( myMouseDraggedIndex > myMousePressedIndex ) myMouseDraggedIndex -= 1;
      Image img = mySnapShotThumbList.remove(myMousePressedIndex);
      mySnapShotThumbList.add(myMouseDraggedIndex,img);
      repaint();
    }
    myIsMouseDragged = false;
    myMousePosY = e.getY();
    repaint();
  }
  @Override
  public void mouseEntered(MouseEvent e) {
    myIsMouseInWindow = true;
    myMousePosY = e.getY();
    repaint();
  }
  @Override
  public void mouseExited(MouseEvent e) {
    myIsMouseInWindow = false;
    repaint();
  }
  @Override
  public void keyPressed(KeyEvent event) {
    int code = event.getKeyCode();
    int index = mySelectedIndex;
    boolean isMove = false;
    switch( code ) {
    case KeyEvent.VK_UP:
    case KeyEvent.VK_LEFT:
      index -= 1;
      isMove = true;
      break;
    case KeyEvent.VK_DOWN:
    case KeyEvent.VK_RIGHT:
      index += 1;
      isMove = true;
      break;
    }
    if( isMove ) {
      if( index < 0 ) index = mySnapShotThumbList.size()-1;
      else if( index >= mySnapShotThumbList.size() ) index = 0;
      if( index >= 0 ) {
        mySelectedIndex = index;
        fireSelectSnapShotEvent(mySelectedIndex,false);
        repaint();
      }
    }
    if( code == KeyEvent.VK_DELETE ) {
      fireDeleteSnapShotEvent(mySelectedIndex);
    }
  }
  @Override
  public void keyTyped(KeyEvent e) {
  }
  @Override
  public void keyReleased(KeyEvent e) {
  }

  //----------------------------------------------------
  // Implement Scrollable interface
  //
  public Dimension getPreferredSize() {
    return new Dimension(myTotalWidth, myTotalHeight*mySnapShotThumbList.size());
  }
  public Dimension getPreferredScrollableViewportSize() {
    return getPreferredSize();
  }
  public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
    int currentPosition = 0;
    if( orientation == SwingConstants.HORIZONTAL ) {
      currentPosition = visibleRect.x;
    }
    else {
      currentPosition = visibleRect.y;
    }
    //Return the number of pixels between currentPosition
    //and the nearest tick mark in the indicated direction.
    if( direction < 0 ) {
      int newPosition = currentPosition - (currentPosition / myMaxUnitIncrement) * myMaxUnitIncrement;
      return (newPosition == 0) ? myMaxUnitIncrement : newPosition;
    }
    else {
      return ((currentPosition / myMaxUnitIncrement) + 1) * myMaxUnitIncrement - currentPosition;
    }
  }
  public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
    if (orientation == SwingConstants.HORIZONTAL) {
      return visibleRect.width - myMaxUnitIncrement;
    }
    else {
      return visibleRect.height - myMaxUnitIncrement;
    }
  }
  public boolean getScrollableTracksViewportWidth() {
    return false;
  }
  public boolean getScrollableTracksViewportHeight() {
    return false;
  }
  public void setMaxUnitIncrement(int pixels) {
    myMaxUnitIncrement = pixels;
  }
  //----------------------------------------------------
  // Event handling
  //
  private void fireSelectSnapShotEvent( int itemIndex, boolean forceSelection ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).selectSnapShot(itemIndex,forceSelection);
    }
  }
  private void fireDeleteSnapShotEvent( int itemIndex ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).deleteSnapShot(itemIndex);
    }
  }
  private void fireMoveSnapShotEvent( int fromItemIndex, int toItemIndex ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).moveSnapShot(fromItemIndex,toItemIndex);
    }
  }
}



/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.xcseis;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import javax.swing.JPanel;

public class csFlowView extends JPanel implements MouseListener, MouseMotionListener, KeyListener {

  public static int TEMP_COUNTER = 0;
  public int myTempID;
  private static final int DEFAULT_VERTICAL_SIZE = 30;
  private static final int DEFAULT_MODULE_WIDTH   = 140;
  private static final int DEFAULT_MODULE_HEIGHT  = 24;
  private static final int MODULE_BORDER  = 3;
  private static final int MODULE_INDENT  = 30;
  private static final int DEFAULT_VERTICAL_MARGIN = 10;
  private static final int DEFAULT_FONT_SIZE = 13;

  private double myModelY2ViewConst;
  private double myModelY2ViewScalar;
  private ArrayList<Integer> myModulePosInFlow;
  private ArrayList<String> myModuleNames;
  private ArrayList<String> myModuleComments;
  private csFlowPanel myPanel;
  private int myCurrentModuleIndex = -1;
  private int myTotalTextSize;
  private ArrayList<csIFlowViewListener> myListeners;
  private FlowViewKeyDispatcher myKeyDispatcher = null;

  public csFlowView( csFlowPanel panel, String text ) {
    myTempID = csFlowView.TEMP_COUNTER++;
    myListeners = new ArrayList<csIFlowViewListener>();
    myPanel = panel;
    myModulePosInFlow = new ArrayList<Integer>();
    myModuleNames     = new ArrayList<String>();
    myModuleComments  = new ArrayList<String>();
    myModelY2ViewConst  = 1.0;
    myModelY2ViewScalar = 1.0;
    myTotalTextSize = text.length();

    boolean prevLineIsComment = false;
    int prevLinePos = 0;
    int counter = 0;
    while( counter < text.length() ) {
      while( counter < text.length() && text.charAt(counter) == ' ' ) {
        counter += 1;
      }
      if( counter == text.length() ) break;
      if( text.charAt(counter) == '$' ) {
        if( prevLineIsComment ) {
          myModulePosInFlow.add( prevLinePos );
        }
        else {
          myModulePosInFlow.add( counter );
        }
        int start = counter+1;
        while( counter < text.length() && text.charAt(counter) != ' ' && text.charAt(counter) != '\n' ) {
          counter += 1;
        }
        myModuleNames.add(text.substring(start,counter));
        prevLineIsComment = false;
      }
//      else if( text.charAt(counter) == '#' ) {
//        int start = counter+1;
//        while( counter < text.length() && text.charAt(counter) != ' ' && text.charAt(counter) != '\n' ) {
//          counter += 1;
//        }
//        myModuleNames.add(text.substring(start,counter));
//      }
      else if( text.charAt(counter) == '\n' ) {
        counter += 1;
      }
      else {
        prevLineIsComment = false;
        if( text.charAt(counter) == '#' ) {
          if( !prevLineIsComment ) {
            prevLineIsComment = true;
            prevLinePos = counter;
          }
        }
        // Forward to next line
        while( counter < text.length() && text.charAt(counter) != '\n' ) {
          counter += 1;
        }
      }
    } // END while
    int numModules = myModuleNames.size();
    int ySize = numModules * DEFAULT_VERTICAL_SIZE + 2*DEFAULT_VERTICAL_MARGIN;
    int xSize = 500;
    setPreferredSize(new Dimension(xSize,ySize));
    setMinimumSize(new Dimension(xSize,ySize));
    addMouseListener(this);
    addMouseMotionListener(this);
//    addKeyListener(this);

    computeModel2View();
  }
  public int addModule( String moduleName, String moduleText ) {
    int addNumLetters = moduleText.length();
    myCurrentModuleIndex += 1;
    if( myCurrentModuleIndex > myModuleNames.size() ) myCurrentModuleIndex = myModuleNames.size();
    int currentPos = myTotalTextSize;
    if( myCurrentModuleIndex < myModuleNames.size() ) currentPos = myModulePosInFlow.get(myCurrentModuleIndex);
    myModulePosInFlow.add( myCurrentModuleIndex, currentPos );
    for( int i = myCurrentModuleIndex+1; i < myModulePosInFlow.size(); i++ ) {
      myModulePosInFlow.set(i, myModulePosInFlow.get(i) + addNumLetters );
    }
    myModuleNames.add( myCurrentModuleIndex, moduleName);
    myTotalTextSize += addNumLetters;
    repaint();
    return currentPos;
  }
  private void deleteModule( int moduleIndex ) {
/*    int numPositionsToDelete = myModulePosInFlow.get(moduleIndex);
    if( moduleIndex < myModuleNames.size()-1 ) {
      numPositionsToDelete = myModulePosInFlow.get(moduleIndex+1) - myModulePosInFlow.get(moduleIndex);
      for( int imodule = moduleIndex; imodule < myModuleNames.size(); imodule++ ) {
        
      }
    }
    else {
      numPositionsToDelete = myTotalTextSize - myModulePosInFlow.get(moduleIndex);
    }
  */
  }
  public void addFlowViewListener( csIFlowViewListener listener ) {
    myListeners.add(listener);
  }
  public void removeFlowViewListener( csIFlowViewListener listener ) {
    myListeners.remove(listener);
  }
  private void computeModel2View() {
    myModelY2ViewConst  = (DEFAULT_VERTICAL_SIZE-DEFAULT_MODULE_HEIGHT)/2 + DEFAULT_VERTICAL_MARGIN;
    myModelY2ViewScalar = DEFAULT_VERTICAL_SIZE;
//      int ypos = i*DEFAULT_VERTICAL_SIZE + (DEFAULT_VERTICAL_SIZE-DEFAULT_MODULE_HEIGHT)/2 + DEFAULT_VERTICAL_MARGIN;
//      ypos = (i+1)*DEFAULT_VERTICAL_SIZE - (DEFAULT_VERTICAL_SIZE-DEFAULT_MODULE_HEIGHT)/2 + DEFAULT_VERTICAL_MARGIN;
//      g2.drawString(text, xpos+DEFAULT_MODULE_WIDTH/2-width/2, ypos-MODULE_BORDER*2);
  }
  public float modelY2View( float yModel ) {
    return( (float)(myModelY2ViewScalar * yModel + myModelY2ViewConst ) );
  }
  public float viewY2Model( float yView ) {
    return( (float)((yView - myModelY2ViewConst) / myModelY2ViewScalar ) );
  }
  public void paintComponent( Graphics g ) {
    Graphics2D g2 = (Graphics2D)g;
    paintAll( g2 );
  }
  private void paintAll( Graphics2D g2 ) {
    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
    g2.setFont(new Font("SansSerif",Font.BOLD,DEFAULT_FONT_SIZE));
    int value = myPanel.getVerticalValue();
    int maxValue = myPanel.getMaxVerticalValue();
    g2.setColor(Color.white);
    g2.fill(getVisibleRect());
    g2.setColor(Color.blue);

    int xpos = 10;
    FontMetrics metrics = g2.getFontMetrics();
    for( int i = 0; i < myModuleNames.size(); i++ ) {
      String text = myModuleNames.get(i);
      int width = metrics.stringWidth(text);
      if( text.compareTo("ELSEIF") == 0 || text.compareTo("ELSE") == 0 || text.compareTo("ENDIF") == 0 || text.compareTo("ENDSPLIT") == 0 ) {
        xpos -= MODULE_INDENT;
      }

      int ypos = i*DEFAULT_VERTICAL_SIZE + (DEFAULT_VERTICAL_SIZE-DEFAULT_MODULE_HEIGHT)/2 + DEFAULT_VERTICAL_MARGIN;
      g2.setColor(Color.gray);
      g2.fillRect(xpos, ypos, DEFAULT_MODULE_WIDTH, DEFAULT_MODULE_HEIGHT );
      g2.setColor(new Color(154,220,154));
      g2.fillRect(xpos+MODULE_BORDER, ypos+MODULE_BORDER, DEFAULT_MODULE_WIDTH-2*MODULE_BORDER, DEFAULT_MODULE_HEIGHT-2*MODULE_BORDER );
      g2.setColor(Color.white);
      g2.drawRect(xpos+MODULE_BORDER, ypos+MODULE_BORDER, DEFAULT_MODULE_WIDTH-2*MODULE_BORDER, DEFAULT_MODULE_HEIGHT-2*MODULE_BORDER );
      g2.setColor(Color.black);
      g2.drawRect(xpos, ypos, DEFAULT_MODULE_WIDTH, DEFAULT_MODULE_HEIGHT );

      if( myCurrentModuleIndex == i ) {
        g2.setColor(Color.black);
        g2.setStroke(new BasicStroke(3.0f));
        g2.drawRect(xpos-1, ypos-1, DEFAULT_MODULE_WIDTH+2, DEFAULT_MODULE_HEIGHT+2 );
        g2.setStroke(new BasicStroke(1.0f));
      }

      g2.setColor(Color.black);
      ypos = (i+1)*DEFAULT_VERTICAL_SIZE - (DEFAULT_VERTICAL_SIZE-DEFAULT_MODULE_HEIGHT)/2 + DEFAULT_VERTICAL_MARGIN; 
      g2.drawString(text, xpos+DEFAULT_MODULE_WIDTH/2-width/2, ypos-MODULE_BORDER*2);

      if( text.compareTo("IF") == 0 || text.compareTo("ELSEIF") == 0 || text.compareTo("ELSE") == 0 || text.compareTo("SPLIT") == 0 ) {
        xpos += MODULE_INDENT;
      }
    }
  }

  public void mouseClicked(MouseEvent e) {
    float ymodel = viewY2Model(e.getY());
    int ymodelInt = (int)ymodel;
    if( myCurrentModuleIndex != ymodelInt ) {
      myCurrentModuleIndex = ymodelInt;
      if( myCurrentModuleIndex >= 0 && myCurrentModuleIndex < myModuleNames.size() ) {
        fireModuleSelectEvent( myCurrentModuleIndex );
      }
      repaint();
    }
  }
  public void mousePressed(MouseEvent e) {
  }
  public void mouseReleased(MouseEvent e) {
  }
  public void mouseEntered(MouseEvent e) {
  }
  public void mouseExited(MouseEvent e) {
  }
  public void mouseDragged(MouseEvent e) {
  }
  public void mouseMoved(MouseEvent e) {
  }
  private void fireModuleSelectEvent( int moduleIndex ) {
    String moduleName = myModuleNames.get(moduleIndex);
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).selectModuleInFlow( moduleIndex, moduleName );
    }
  }
  private void fireModuleDeleteEvent( int moduleIndex ) {
    String moduleName = myModuleNames.get(moduleIndex);
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).deleteModuleInFlow( moduleIndex, moduleName );
    }
  }
  public void grabKeys() {
    KeyboardFocusManager manager = KeyboardFocusManager.getCurrentKeyboardFocusManager();
    myKeyDispatcher = new FlowViewKeyDispatcher(this);
    manager.addKeyEventDispatcher( myKeyDispatcher );
  }
  public void releaseKeys() {
    KeyboardFocusManager manager = KeyboardFocusManager.getCurrentKeyboardFocusManager();
    manager.removeKeyEventDispatcher( myKeyDispatcher );
  }
  public void keyTyped(KeyEvent e) {
  }
  public void keyPressed(KeyEvent e) {
    if( e.getKeyCode() == KeyEvent.VK_DOWN ) {
      if( myCurrentModuleIndex < myModuleNames.size()-1 ) {
        myCurrentModuleIndex += 1;
        fireModuleSelectEvent( myCurrentModuleIndex );
        repaint();
      }
    }
    else if( e.getKeyCode() == KeyEvent.VK_UP ) {
      if( myCurrentModuleIndex > 0 ) {
        myCurrentModuleIndex -= 1;
        fireModuleSelectEvent( myCurrentModuleIndex );
        repaint();
      }
    }
    else if( e.getKeyCode() == KeyEvent.VK_DELETE ) {
      if( myCurrentModuleIndex < 0 || myCurrentModuleIndex >= myModuleNames.size() ) return;
      if( myCurrentModuleIndex == myModuleNames.size()-1 ) {
        deleteModule( myCurrentModuleIndex );
        fireModuleDeleteEvent( myCurrentModuleIndex );
      }
      else {
        fireModuleDeleteEvent( myCurrentModuleIndex );
      }
      repaint();
    }
  }
  public void keyReleased(KeyEvent e) {
  }
  
  //***********************************************************************************
  //***********************************************************************************
  //***********************************************************************************
  private class FlowViewKeyDispatcher implements KeyEventDispatcher {
    private KeyListener myListener;
    public FlowViewKeyDispatcher( KeyListener listener ) {
      super();
      myListener = listener;
    }
    @Override
    public boolean dispatchKeyEvent(KeyEvent e) {
      System.out.println("Dispatcher key event " + e.getKeyChar());
      if( e.getID() == KeyEvent.KEY_PRESSED) {
        myListener.keyPressed( e );
      }
//        else if( e.getID() == KeyEvent.KEY_RELEASED ) {
//          myListener.keyReleased( e );
//        }
//        else if( e.getID() == KeyEvent.KEY_TYPED ) {
//          myListener.keyTyped( e );
//        }
      return false;
    }
  }
}



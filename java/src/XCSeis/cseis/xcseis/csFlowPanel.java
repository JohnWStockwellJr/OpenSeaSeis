/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.xcseis;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.BadLocationException;
import javax.swing.text.StyledDocument;

public class csFlowPanel extends JPanel implements MouseListener, DocumentListener {
  private JTextPane myTextPaneFlow;
  private csFlowView myFlowView = null;
  private JScrollPane myPaneFlow;
  private StyledDocument myFlowModelDoc;
  private boolean myIsFlowView = false;
  private XCSeis myXCSeis;
  private JPopupMenu myPopupMenu;
  private JMenuItem myMenuOpenSeaView;
  private JMenuItem myMenuDummy;
  private String myCurrentFilename;
  private String myFlowPathName;
  private boolean myIsUpdated;
  
  public csFlowPanel( XCSeis xcseis, String flowPathName ) {
    super( new BorderLayout() );
    myIsUpdated = false;
    myXCSeis = xcseis;
    myFlowPathName = flowPathName;
    myPopupMenu = new JPopupMenu();
    myMenuOpenSeaView = new JMenuItem("Open data set with SeaView");
    myMenuOpenSeaView.setEnabled(false);
    myPopupMenu.add( myMenuOpenSeaView );
    myMenuOpenSeaView.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        myXCSeis.openInSeaView(myCurrentFilename);
      }
    });
    myFlowModelDoc = new csFlowDocument();
    myTextPaneFlow  = new JTextPane();
    myTextPaneFlow.setFont( myXCSeis.getSystemFont() );
    myTextPaneFlow.setEditable( true );
    myTextPaneFlow.setStyledDocument(myFlowModelDoc);
    myPaneFlow       = new JScrollPane( myTextPaneFlow, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );
    add( myPaneFlow,  BorderLayout.CENTER );
    addMouseListener(this);
    myTextPaneFlow.addMouseListener(this);
  }
  public String flowPathName() {
    return myFlowPathName;
  }
  // !CHANGE!
  // listener gets mixed up if several flows are open. Each flowview object is linked to same listener..??
  public void setFlowView( boolean doSet ) {
    if( myIsFlowView != doSet ) {
      myIsFlowView = doSet;
      if( myIsFlowView ) {
        myFlowView = new csFlowView(this,myTextPaneFlow.getText());
        myFlowView.addFlowViewListener(myXCSeis);
        myPaneFlow.remove( myTextPaneFlow );
        myPaneFlow.setViewportView( myFlowView );
        myFlowView.grabKeys();
      }
      else {
        myPaneFlow.remove( myFlowView );
        myPaneFlow.setViewportView( myTextPaneFlow );
        myFlowView.removeFlowViewListener(myXCSeis);
        myFlowView.releaseKeys();
        myFlowView = null;
      }
      myPaneFlow.invalidate();
      myPaneFlow.repaint();
    }
  }
  public int getVerticalValue() {
    return myPaneFlow.getVerticalScrollBar().getValue();
  }
  public int getMaxVerticalValue() {
    return myPaneFlow.getVerticalScrollBar().getMaximum();
  }
  public void updateModel( StringBuffer strBuffer ) {
    try {
      csFlowDocument doc = new csFlowDocument();
      doc.insertString( 0, strBuffer.toString(), null );
      myFlowModelDoc.removeDocumentListener(this);
      myFlowModelDoc = doc;
      myTextPaneFlow.setStyledDocument(myFlowModelDoc);
      myIsUpdated = false;
      SwingUtilities.invokeLater( new Runnable() {
        public void run() {
          myFlowModelDoc.addDocumentListener(csFlowPanel.this);
        }
      });
    }
    catch( BadLocationException e ) {
      e.printStackTrace();
      return;
    }
  }
  private void showPopupMenu( int xpos, int ypos ) {
    int posCaret = myTextPaneFlow.getCaretPosition();
    
    String text = myTextPaneFlow.getText();
    int pos1 = posCaret;
    int pos2 = posCaret;
    while( pos1 > 0 && text.charAt(pos1) != ' ' && text.charAt(pos1) != '"' && text.charAt(pos1) != '\n' ) {
      pos1 -= 1;
    }
    pos1 += 1;
    while( pos2 < text.length() && text.charAt(pos2) != ' ' && text.charAt(pos2) != '"' && text.charAt(pos2) != '\n' ) {
      pos2 += 1;
    }
    if( pos2 < pos1 || pos2 < 0 ) return;
    myTextPaneFlow.select(pos1, pos2); // Highlight selected text
    String selectedText = text.substring(pos1, pos2).trim();
    // a) Check whether this is a module name
    // b) check whether this is a valid file name
    if( selectedText.charAt(0) == '$' ) {
      JPopupMenu popupMenu = new JPopupMenu();
      JMenuItem menu = new JMenuItem("View module help");
      final String moduleName = selectedText.substring(1);
      menu.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          myXCSeis.showModuleHelp(moduleName);
        }
      });
      popupMenu.add( menu );
      xpos -= myPaneFlow.getHorizontalScrollBar().getValue();
      ypos -= myPaneFlow.getVerticalScrollBar().getValue();
      popupMenu.show( this, xpos, ypos );
      return;
    }
    myCurrentFilename = text.substring(pos1, pos2).trim();
    
    if ( !new java.io.File(myCurrentFilename).exists() ) {
      boolean success = false;
      String flowPath = myXCSeis.pathActiveFlow();
      int index = flowPath.lastIndexOf( java.io.File.separatorChar );
      if( index >= 0 ) {
        String flowRootPath = flowPath.substring(0,index);
        String newName = flowRootPath + java.io.File.separatorChar + myCurrentFilename;
        success = new java.io.File(newName).exists();
        if( success ) myCurrentFilename = newName;
      }
      if( !success ) {
        myXCSeis.updateStatus("File does not exist or is not readable:\n" + myCurrentFilename);
//        JOptionPane.showMessageDialog(this, "File does not exist or is not readable:\n" + myCurrentFilename,
//                "File not found", JOptionPane.ERROR_MESSAGE);
        return;
      }
    }
//    myMenuOpenSeaView.setText("Open " + myCurrentFilename + " with SeaView");
    myMenuOpenSeaView.setText("Open file with SeaView");
    myMenuOpenSeaView.setEnabled(true);
    xpos -= myPaneFlow.getHorizontalScrollBar().getValue();
    ypos -= myPaneFlow.getVerticalScrollBar().getValue();
    myPopupMenu.show( this, xpos, ypos );
  }
  public void addModule( String moduleName, String moduleText ) {
    int posCaret = myTextPaneFlow.getCaretPosition();
    if( myIsFlowView ) {
//      myFlowView.
//      myFlowView = new csFlowView(this,myTextPaneFlow.getText());
      posCaret = myFlowView.addModule( moduleName, moduleText );
    }
//    else {
      try {
        myFlowModelDoc.insertString( posCaret, moduleText, null );
      } catch (BadLocationException ex) {
        JOptionPane.showMessageDialog( this, "Cannot add module to flow./nWhy? ..not sure..", "Error", JOptionPane.ERROR_MESSAGE);
      }
//    }
  }
  public String getText() {
    return myTextPaneFlow.getText();
  }
  public void updateFont( Font font ) {
    myTextPaneFlow.setFont(font);
  }
  public void mouseClicked(MouseEvent e) {
    if( SwingUtilities.isRightMouseButton(e) ) {
      myTextPaneFlow.setCaretPosition(myTextPaneFlow.viewToModel(e.getPoint())); 
      showPopupMenu(e.getX(),e.getY());
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
  public void insertUpdate(DocumentEvent e) {
    doUpdate(true);
  }
  public void removeUpdate(DocumentEvent e) {
    doUpdate(true);
  }
  public void changedUpdate(DocumentEvent e) {
    doUpdate(true);
  }
  private void doUpdate( boolean doUpdate ) {
    if( !myIsUpdated && doUpdate ) {
      myXCSeis.updateFlowTab( this );
      myIsUpdated = doUpdate;
    }
  }
  public void resetUpdated() {
    myIsUpdated = false;
  }
  public boolean isUpdated() {
    return myIsUpdated;
  }

}



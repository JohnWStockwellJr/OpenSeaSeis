/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.xcseis;

import javax.swing.*;
import javax.swing.text.*;
import javax.swing.text.html.*;
import cseis.jni.csNativeModuleHelp;
import cseis.seaview.csSeaViewActions;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

public class csModuleListingPanel extends JScrollPane {
  private final JTextPane myTextPane;
  private csIModuleListingListener myListener;
  private ModuleSelectionPopupMenu myPopupMenu;
  private int myTextLength;
  private String myLinkedModuleName;
  private boolean myIsHyperlinkActive;

  public csModuleListingPanel() {
    super( JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );

    myIsHyperlinkActive = false;
    myLinkedModuleName = "";
    myTextPane = new JTextPane();
    myTextPane.setEditable(false);
    myTextPane.setEditorKit(new HTMLEditorKit());
    csNativeModuleHelp help = new csNativeModuleHelp();
    myTextPane.setText("<html>" + help.moduleHtmlListing() + "</html>");
    MutableAttributeSet attrs = myTextPane.getInputAttributes();
    int size = StyleConstants.getFontSize(attrs);
    StyleConstants.setFontSize(attrs, size - 1);
    StyledDocument doc = myTextPane.getStyledDocument();
    doc.setCharacterAttributes(0, doc.getLength() + 1, attrs, false);
    myTextPane.setCaretPosition(0);
    myTextLength = myTextPane.getText().length();

    myTextPane.addMouseListener( new MouseAdapter() {
      @Override
      public void mousePressed(MouseEvent e) {
        if(SwingUtilities.isRightMouseButton(e)) {
          if( myIsHyperlinkActive ) {
            myPopupMenu = new ModuleSelectionPopupMenu();
            myPopupMenu.show( myTextPane, e.getX(), e.getY() );
          }
        }
      }
    });
    myTextPane.addHyperlinkListener(new HyperlinkListener() {
      public void hyperlinkUpdate( HyperlinkEvent hle ) {
        myLinkedModuleName = hle.getDescription().replace("#", "").trim();
        if( hle.getEventType().equals(HyperlinkEvent.EventType.ACTIVATED ) ) {
          myIsHyperlinkActive = true;
          if( myListener != null ) myListener.selectModule(myLinkedModuleName);
        }
        else if( hle.getEventType().equals( HyperlinkEvent.EventType.EXITED ) ) {
          myIsHyperlinkActive = false;
        }
        else if( hle.getEventType().equals( HyperlinkEvent.EventType.ENTERED ) ) {
          myIsHyperlinkActive = true;
        }
      }
    });
    super.setViewportView( myTextPane );
  }
  public void addModuleListingListener( csIModuleListingListener listener ) {
    myListener = listener;
  }

  //-----------------------------------------------------------------------------
  class ModuleSelectionPopupMenu extends JPopupMenu {
    ModuleSelectionPopupMenu() {
      super();
      JMenuItem itemAddModule = new JMenuItem( "Add module " + myLinkedModuleName + " to flow" );
      add( itemAddModule );
      itemAddModule.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed( ActionEvent e ) {
          if( myListener != null ) {
            myListener.addModule( myLinkedModuleName );
          }
        }
      });
      
    }
  }
}

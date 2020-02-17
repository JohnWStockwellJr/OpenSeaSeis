/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/**
 * Recent file menu.<br>
 * A recent file menu keeps path names of files which have recently been accessed.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csRecentFileMenu extends JMenu {
  /**
   * Maximum width in pixels of text in menu items. File names are shortened to fit this width
   */
  public int MAX_WIDTH_PIXELS = 450;
  /**
   * Maximum number of files allowed in recent file menu
   */
  private int myMaxNumFiles;
  /**
   * List of file menu listeners
   */
  private Vector<csFileMenuListener> myListeners;
  /**
   * true if menu shall be disabled automatically when it contains no items.
   * false if menu shall always stay enabled.
   */
  private boolean myAutoDisable = true;
//---------------------------------------------------------------
  public csRecentFileMenu( String name ) {
    this(name,0);
  }
  public csRecentFileMenu( String name, int maxNumFiles ) {
    super( name );
    myMaxNumFiles = maxNumFiles;
    myListeners = new Vector<csFileMenuListener>();
    setEnabled(false);
  }
  public void addFile( String filenameFullPath ) {
    if( super.getItemCount() == 0 ) {
      super.setEnabled(true);
    }
//    String fileName = file.getAbsolutePath();
    boolean found = false;
    int itemPos = -1;
    int nItems = super.getItemCount();
    for( int i = 0; i < nItems; i++ ) {
      if( super.getItem(i) instanceof RecentFileItem ) {
        if( ((RecentFileItem)getItem(i)).filenameFullPath.compareTo( filenameFullPath ) == 0 ) {
          found = true;
          itemPos = i;
          break;
        }
      }
      else {
        // Throw exception
      }
    }
    if( !found ) {
      if( super.getItemCount() == myMaxNumFiles ) {
        super.remove( myMaxNumFiles-1 );
      }
    }
    else if( itemPos != 0 ) {
      super.remove( itemPos );
    }
    else {
      return;
    }
    RecentFileItem item = new RecentFileItem( filenameFullPath );
    super.add( item, 0 );
    updateMnemonics( 0 );
  }
  public void removeFile( String filenameFullPath ) {
    int nItems = super.getItemCount();
    for( int i = 0; i < nItems; i++ ) {
      if( super.getItem( i ) instanceof RecentFileItem ) {
        if( ( ( RecentFileItem )getItem( i ) ).filenameFullPath.compareTo( filenameFullPath ) == 0 ) {
          super.remove( i );
          updateMnemonics(i);
          break;
        }
      }
    }
    if( super.getItemCount() == 0 && myAutoDisable ) {
      super.setEnabled(false);
    }
  }
  public java.util.List<String> getFileList() {
    int nItems = getItemCount();
    Vector<String> list = new Vector<String>(nItems);
    for( int i = 0; i < nItems; i++ ) {
      if( getItem(i) instanceof RecentFileItem ) {
        list.add(((RecentFileItem)getItem( i )).filenameFullPath );
      }
    }
    return list;
  }
  public int getMaxNumFiles() {
    return myMaxNumFiles;
  }
  public void setMaxNumFiles( int maxNumFiles ) {
    myMaxNumFiles = maxNumFiles;
  }
  public void setAutoDisable( boolean doSetAutoDisable ) {
    if( myAutoDisable != doSetAutoDisable && getItemCount() == 0) {
      if( doSetAutoDisable ) {
        setEnabled(false);
      }
      else {
        setEnabled(true);
      }
    }
    myAutoDisable = doSetAutoDisable;
  }
  private void updateMnemonics( int startPos ) {
    int nItems = super.getItemCount();
    for( int i = startPos; i < nItems; i++ ) {
      if( super.getItem( i ) instanceof RecentFileItem ) {
        RecentFileItem item = ( RecentFileItem )getItem( i );
        item.setText( (i+1) + " " + item.filenameShort );
        item.setMnemonic( getPositionMnemonic(i) );
        item.position = i;
      }
    }
  }
  private int getPositionMnemonic( int pos ) {
    int mnemonic = -1;
    if( pos < 9 ) {
      switch( pos ) {
        case 0:
          mnemonic = KeyEvent.VK_1;
          break;
        case 1:
          mnemonic = KeyEvent.VK_2;
          break;
        case 2:
          mnemonic = KeyEvent.VK_3;
          break;
        case 3:
          mnemonic = KeyEvent.VK_4;
          break;
        case 4:
          mnemonic = KeyEvent.VK_5;
          break;
        case 5:
          mnemonic = KeyEvent.VK_6;
          break;
        case 6:
          mnemonic = KeyEvent.VK_7;
          break;
        case 7:
          mnemonic = KeyEvent.VK_8;
          break;
        case 8:
          mnemonic = KeyEvent.VK_9;
          break;
      }
    }
    return mnemonic;
  }
  private void moveToFirstPosition( RecentFileItem item ) {
    int nItems = super.getItemCount();
    for( int i = 0; i < nItems; i++ ) {
      if( getItem( i ).equals( item ) ) {
        if( item.position != 0 ) {
          super.remove( i );
          super.add( item, 0 );
          updateMnemonics( 0 );
        }
        break;
      }
    }
  }
  public void addFileMenuListener( csFileMenuListener listener ) {
    if( !myListeners.contains(listener) ) {
      myListeners.add(listener);
    }
  }
  public void removeFileMenuListener( csFileMenuListener listener ) {
    if( myListeners.contains(listener) ) {
      myListeners.remove(listener);
    }
  }
  public void fireFileSelectedEvent( csFileMenuEvent event ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      ((csFileMenuListener)myListeners.get(i)).fileSelected( event );
    }
  }
//------------------------------------------------------
  class RecentFileItem extends JMenuItem {
    String filenameFullPath;
    String filenameShort;
    int position = 0;
    RecentFileItem( String filenameFullPath_in ) {
      super("TEMP_NAME");
      this.filenameFullPath = filenameFullPath_in;

      FontMetrics metrics = getFontMetrics(getFont());
      int size = metrics.stringWidth(filenameFullPath);
      int charsToChop = (int)( ((double)(size-MAX_WIDTH_PIXELS))/(double)size * filenameFullPath.length() );
      if( charsToChop > 0 ) {
        filenameShort = "..." + filenameFullPath.substring( charsToChop );
      }
      else {
        filenameShort = filenameFullPath;
      }
      setToolTipText( filenameFullPath );
      addActionListener(new ActionListener() {
        public void actionPerformed( ActionEvent e ) {
          if( position != 0 ) moveToFirstPosition( (RecentFileItem)e.getSource() );
          fireFileSelectedEvent(new csFileMenuEvent(e.getSource(),new File(filenameFullPath)) );
        }
      });
    }
  }
}



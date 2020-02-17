/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToggleButton;
import javax.swing.border.BevelBorder;

/**
 * Panel which can be docked and undocked from a csDockPaneManager.<br>
 * A csDockPane provides a top ribbon containing buttons for the various dock options.<br>
 * 
 * To create a new csDockPane, provide a JPanel -the content of the csDockPane- and the csDockManager.
 * Next, add the csDockPane to the csDockPaneManager.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csDockPane extends JPanel {
  public static final Color COLOR_ACTIVE = new Color(255, 255, 124);
  public static final Color COLOR_INACTIVE = new JLabel("123").getBackground();

  private JScrollPane myScrollPaneTop;
  private JPanel myTopBar;
  private JLabel myLabelTitle;
  private JButton myButtonClose;
  private JButton myButtonMaximize;
  private JButton myButtonHide;
  private JButton myButtonDocking;
  private JToggleButton myButtonShowScroll;
  private JToggleButton myButtonSync;
  
  private csDockPaneManager myManager;
  private JPanel myMainPane;
  private csDockPaneAttr myAttr;
  private int myDockingState;
  
  private ArrayList<csIDockPaneListener> myListeners;
  private ArrayList<csIDockPaneSyncListener> mySyncListeners;

  public csDockPane( csDockPaneManager manager, JPanel mainPane, String title, csDockPaneButtonSelection bselect ) {
    this( manager, mainPane, title, bselect, title );
  }  
  public csDockPane( csDockPaneManager manager, JPanel mainPane, String title, csDockPaneButtonSelection bselect,
          String titleDescription ) {
    super( new BorderLayout() );
    myManager  = manager;
    myMainPane = mainPane;
    myAttr = new csDockPaneAttr( title );
    myAttr.titleDescription = titleDescription;
    myDockingState = csDockPaneManager.STATE_DOCKED;
    myListeners = new ArrayList<csIDockPaneListener>();
    mySyncListeners = new ArrayList<csIDockPaneSyncListener>();

    setBorder( BorderFactory.createBevelBorder(BevelBorder.RAISED) );

    myLabelTitle = new JLabel( title, JLabel.CENTER );
    myLabelTitle.setToolTipText(titleDescription);
    myLabelTitle.setPreferredSize( new Dimension(0,myLabelTitle.getPreferredSize().height) );
    
    myLabelTitle.addMouseListener( new MouseListener() {
      @Override
      public void mouseClicked(MouseEvent e) {
        // setActivePane already implicitely called when mouse button is pressed
        // myManager.setActivePane( csDockPane.this );
      }
      @Override
      public void mousePressed(MouseEvent e) {
        myManager.pressPane( csDockPane.this, e.getLocationOnScreen() );
      }
      @Override
      public void mouseReleased(MouseEvent e) {
        // releasePane may not be needed here since it's already done in mouseExited() below.
        // ..no harm if called more than once
        myManager.releasePane( csDockPane.this );
      }
      @Override
      public void mouseEntered(MouseEvent e) {
      }
      @Override
      public void mouseExited(MouseEvent e) {
        myManager.releasePane( csDockPane.this );
      }
    });
    myLabelTitle.addMouseMotionListener( new MouseMotionListener() {
      @Override
      public void mouseDragged(MouseEvent e) {
        myManager.dragPane( csDockPane.this, e.getLocationOnScreen() );
      }
      @Override
      public void mouseMoved(MouseEvent e) {
        //  Nothing
      }
    });

    myButtonClose = new JButton( csDockPaneActions.getIcon(csDockPaneActions.CloseAction) );
    myButtonMaximize = new JButton( csDockPaneActions.getIcon(csDockPaneActions.MaximizeAction) );
    myButtonHide = new JButton( csDockPaneActions.getIcon(csDockPaneActions.HideAction) );
    myButtonShowScroll = new JToggleButton( csDockPaneActions.getIcon(csDockPaneActions.ScrollbarAction) );
    myButtonSync = new JToggleButton( csDockPaneActions.getIcon(csDockPaneActions.SyncAction) );
    myButtonSync.setSelected(false);
    myButtonDocking = new JButton( csDockPaneActions.getIcon(csDockPaneActions.UndockAction) );

    myButtonClose.setToolTipText( csDockPaneActions.ACTION_DESC[csDockPaneActions.CloseAction] );
    myButtonHide.setToolTipText( csDockPaneActions.ACTION_DESC[csDockPaneActions.HideAction] );
    myButtonShowScroll.setToolTipText( csDockPaneActions.ACTION_DESC[csDockPaneActions.ScrollbarAction] );
    myButtonSync.setToolTipText( csDockPaneActions.ACTION_DESC[csDockPaneActions.SyncAction] );
    myButtonDocking.setToolTipText( csDockPaneActions.ACTION_DESC[csDockPaneActions.UndockAction] );
    myButtonMaximize.setToolTipText( csDockPaneActions.ACTION_DESC[csDockPaneActions.MaximizeAction] );

    Dimension preferredSize = new Dimension( 24, 24 );
    myButtonClose.setPreferredSize(preferredSize);
    myButtonMaximize.setPreferredSize(preferredSize);
    myButtonHide.setPreferredSize(preferredSize);
    myButtonShowScroll.setPreferredSize(preferredSize);
    myButtonSync.setPreferredSize(preferredSize);
    myButtonDocking.setPreferredSize(preferredSize);

    int xp = 0;
    JPanel panelButtons = new JPanel( new GridBagLayout());
    if( bselect.sync ) {
      panelButtons.add( myButtonSync, new GridBagConstraints(
          xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
    if( bselect.scrollbars ) {
      panelButtons.add( myButtonShowScroll, new GridBagConstraints(
          xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
    if( bselect.maximize ) {
      panelButtons.add( myButtonMaximize, new GridBagConstraints(
          xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
    if( bselect.dock ) {
      panelButtons.add( myButtonDocking, new GridBagConstraints(
          xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
    if( bselect.hide ) {
      panelButtons.add( myButtonHide, new GridBagConstraints(
          xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
    if( bselect.close ) {
      panelButtons.add( myButtonClose, new GridBagConstraints(
          xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
    
    myTopBar = new JPanel(new BorderLayout());
    myTopBar.setBorder(BorderFactory.createLineBorder(Color.gray));
    if( myAttr.isActive ) {
      myTopBar.setBackground(COLOR_ACTIVE);
    }
    else {
      myTopBar.setBackground(COLOR_INACTIVE);
    }
    myTopBar.add( myLabelTitle, BorderLayout.CENTER );
    myTopBar.add( panelButtons, BorderLayout.EAST );

    // Place top bar in viewport of scroll pane to avoid layout problems when squeezed
    myScrollPaneTop = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_NEVER, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    myScrollPaneTop.getViewport().add(myTopBar);
    myScrollPaneTop.setBorder(BorderFactory.createEmptyBorder());

    add(myScrollPaneTop,BorderLayout.NORTH);
    add(myMainPane,BorderLayout.CENTER);
    
    myButtonClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        closePane();
      }
    });
    myButtonDocking.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        changeDockingState();
        if( myButtonDocking.getToolTipText().compareTo(csDockPaneActions.ACTION_DESC[csDockPaneActions.DockAction]) == 0 ) {
          myButtonDocking.setIcon( csDockPaneActions.getIcon(csDockPaneActions.UndockAction) );
          myButtonDocking.setToolTipText(csDockPaneActions.ACTION_DESC[csDockPaneActions.UndockAction]);
        }
        else {
          myButtonDocking.setIcon( csDockPaneActions.getIcon(csDockPaneActions.DockAction) );
          myButtonDocking.setToolTipText(csDockPaneActions.ACTION_DESC[csDockPaneActions.DockAction]);
        }
      }
    });
    myButtonMaximize.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        fireEventMaximizePane();
      }
    });
    myButtonHide.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        hidePane();
      }
    });
    myButtonShowScroll.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        setScrollVisible( myButtonShowScroll.isSelected() );
      }
    });
    myButtonSync.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if( myButtonSync.isSelected() ) {
          myButtonSync.setIcon( csDockPaneActions.getIcon(csDockPaneActions.UnsyncAction) );
          myButtonSync.setToolTipText(csDockPaneActions.ACTION_DESC[csDockPaneActions.UnsyncAction]);
          setSync( true );
        }
        else {
          myButtonSync.setIcon( csDockPaneActions.getIcon(csDockPaneActions.SyncAction) );
          myButtonSync.setToolTipText(csDockPaneActions.ACTION_DESC[csDockPaneActions.SyncAction]);
          setSync( false );
        }
      }
    });
  }
  //-------------------------------------------------------------
  public Component getTabComponent() {
    int height = 14;
    JPanel panel = new JPanel( new BorderLayout() );
    panel.setOpaque( false );
    String text = attr().title + "  ";
    JLabel label = new JLabel( text, JLabel.LEFT );
    label.setOpaque( false );
    label.setBackground( new Color(0, 0, 0, 255) );
    label.setFocusable( true );
    
    int labelWidth = label.getFontMetrics( label.getFont() ).stringWidth( text );
    if( labelWidth > myManager.getMaxTabWidth() ) label.setPreferredSize( new Dimension( myManager.getMaxTabWidth(), height ) );

    JButton button = new JButton( csDockPaneActions.getIcon(csDockPaneActions.CloseTabAction) );
    button.setToolTipText( csDockPaneActions.ACTION_DESC[csDockPaneActions.CloseTabAction] );
    button.setPreferredSize( new Dimension(height+2,height) );
    button.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        closePane();
      }
    });
    button.setFocusable( false );
    button.setRolloverEnabled( true );
    panel.add( label, BorderLayout.CENTER );
    panel.add( button, BorderLayout.EAST );
    return panel;
  }
  public int getState() {
    return myDockingState;
  }
  public void setState( int dockingState ) {
    myDockingState = dockingState;
  }
  public csDockPaneAttr attr() {
    return new csDockPaneAttr( myAttr );
  }
  public void setAttr( csDockPaneAttr attr ) {
    myAttr = new csDockPaneAttr( attr );
    setTitle( myAttr.title, myAttr.titleDescription );
  }
  public void setTitle( String title, String titleDescription ) {
    myAttr.title = title;
    myAttr.titleDescription = titleDescription;
    myLabelTitle.setText(title);
    myLabelTitle.setToolTipText(titleDescription);
  }
  protected void setActive( boolean active ) {
    if( myAttr.isActive ) myTopBar.setBackground(COLOR_INACTIVE);
    myAttr.isActive = active;
    if( myAttr.isActive ) myTopBar.setBackground(COLOR_ACTIVE);
  }
  public JPanel getMainPanel() {
    return myMainPane;
  }
  public boolean closePane() {
/*
    int option = JOptionPane.showConfirmDialog( this,
      "Are you sure you want to close data set\n" +
            "'"+myAttr.title+"' ?",
      "Close panel",
      JOptionPane.YES_NO_OPTION );
    if( option != JOptionPane.YES_OPTION ) {
      return false;
    }
  */  
//    fireEventHidePane(); // Tell pane listeners to hide pane
    fireEventClosePane(); // Tell pane listeners to close pane
    return true;
  }
  public void showPane() {
    setVisible(true);
  }
  public void hidePane() {
    fireEventHidePane();
    setVisible(false);
    if( !myAttr.isDocked ) {
      changeDockingState();
    }
  }
  public void setSync( boolean isSync ) {
    myAttr.isSync = isSync;
    fireEventSync( isSync );
  }
  public void setTopBarVisible( boolean isTitleBarVisible ) {
    myScrollPaneTop.setVisible(isTitleBarVisible);
    revalidate();
  }
  public void setScrollVisible( boolean isScrollVisible ) {
    myAttr.isScrollVisible = isScrollVisible;
    fireEventScrollbars( isScrollVisible );
  }
  public void changeDockingState() {
    fireEventDockPane( !myAttr.isDocked );
    myAttr.isDocked = !myAttr.isDocked;
  }
  public void addSyncListener( csIDockPaneSyncListener listener ) {
    mySyncListeners.add(listener);
  }
  public void removeSyncListener( csIDockPaneSyncListener listener ) {
    mySyncListeners.remove(listener);
  }
  public void addPaneListener( csIDockPaneListener listener ) {
    myListeners.add(listener);
  }
  public void removePaneListener( csIDockPaneListener listener ) {
    myListeners.remove(listener);
  }
  private void fireEventHidePane() {
    for( int i = 0 ; i < myListeners.size(); i++ ) {
      myListeners.get(i).hidePane( this );
    }
  }
  private void fireEventClosePane() {
    for( int i = 0 ; i < myListeners.size(); i++ ) {
      myListeners.get(i).closePane( this );
    }
  }
  private void fireEventMaximizePane() {
    for( int i = 0 ; i < myListeners.size(); i++ ) {
      myListeners.get(i).maximizePane( this );
    }
  }
  private void fireEventDockPane( boolean dock ) {
    for( int i = 0 ; i < myListeners.size(); i++ ) {
      myListeners.get(i).dockPane( this, dock );
    }
  }
  private void fireEventScrollbars( boolean showScrollbars ) {
    for( int i = 0 ; i < myListeners.size(); i++ ) {
      myListeners.get(i).scrollbars( showScrollbars );
    }
  }
  private void fireEventSync( boolean isSync ) {
    for( int i = 0 ; i < mySyncListeners.size(); i++ ) {
      mySyncListeners.get(i).syncStateChanged( this, isSync );
    }
  }
}


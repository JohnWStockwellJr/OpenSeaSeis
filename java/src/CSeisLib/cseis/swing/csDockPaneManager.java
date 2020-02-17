/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.awt.event.ActionEvent;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * JPanel containing and managing a collection of csDockPane objects.<br>
 * A csDockPaneManager contains and manages a collection of csDockPane objects.
 * It defines several layout options, provides drag and sync functionality etc.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csDockPaneManager extends JPanel implements csIDockPaneListener {
  public static final int LAYOUT_TABS       = 0;
  public static final int LAYOUT_ONE_ROW    = 1;
  public static final int LAYOUT_ONE_COL    = 2;
  public static final int LAYOUT_TWO_ROWS   = 3;
  public static final int LAYOUT_TWO_COLS   = 4;
  public static final int LAYOUT_THREE_ROWS = 5;
  public static final int LAYOUT_THREE_COLS = 6;
  
  public static final String LAYOUT_TEXT[] = {
    "Tabbed panes",
    "One row",
    "One column",
    "Two rows",
    "Two columns",
    "Three rows",
    "Three columns",
  };
  
  public static final int STATE_DOCKED   = 11;
  public static final int STATE_UNDOCKED = 12;
  public static final int STATE_HIDDEN   = 13;

  /// Layout manager for displaying Swing panels (=dock panes)
  private GridLayout myGridLayout;
  /// Layout option. For valid layouts, see above (LAYOUT_ONE_ROW...)
  private int myLayoutOption;
  /// Panes contained in manager
  private ArrayList<csDockPane> myPanes;
  /// Number of visible panes (all panes which are not 'hidden')
  private int myNumPanesVisible;
  /// Index of currently active pane. -1 if no pane is active.
  private int myActivePaneIndex;

  private JTabbedPane myTabbedPane;
  private boolean myIsNewPaneInProgress;
  
  // Pane dragging
  private Point myDragPoint;
  private Point myDragPointOrigin;
  private int myDragPaneIndex;
  private boolean myIsDragging = false;
  private BufferedImage myDragImage = null;
  private ArrayList<csIDockPaneSelectionListener> myListeners;
  private boolean myDoActivateNewPanes;
  private boolean myShowActivePaneFeedback;

  private int myMaxTabWidth;
  
  public csDockPaneManager() {
    this( false );
  }
  public csDockPaneManager( boolean doActivateNewPanes ) {
    myNumPanesVisible = 0;
    myActivePaneIndex = -1;
    myDoActivateNewPanes = doActivateNewPanes;
    myLayoutOption = csDockPaneManager.LAYOUT_ONE_ROW;
    myGridLayout = new GridLayout();
    myTabbedPane = null;
    myIsNewPaneInProgress = false;
    myShowActivePaneFeedback = true;
    myMaxTabWidth = 9999;
    
    myPanes = new ArrayList<csDockPane>();
    myListeners  = new ArrayList<csIDockPaneSelectionListener>();
    setLayout( myGridLayout );
  }
  public void setMaxTabWidth( int maxTabWidth ) {
    myMaxTabWidth = maxTabWidth > 50 ? maxTabWidth : 50;
    if( myTabbedPane != null ) {
      int numTabs = myTabbedPane.getTabCount();
      if( myPanes.size() < numTabs ) numTabs = myPanes.size();
      for( int i = 0; i < numTabs; i++ ) {
        myTabbedPane.setTabComponentAt( i, myPanes.get(i).getTabComponent() );
      }
    }
  }
  public int getMaxTabWidth() {
    return myMaxTabWidth;
  }
  // Show graphical feedback when active pane is selected (=highlight panel for a short time)
  public void setActivePaneFeedback( boolean showFeedback ) {
    myShowActivePaneFeedback = showFeedback;
  }
  public static int getNumLayoutOptions() {
    return LAYOUT_TEXT.length;
  }
  public static ImageIcon getIcon( int layoutOption ) {
    switch( layoutOption ) {
      case csDockPaneManager.LAYOUT_TABS:
        return cseis.resources.csResources.getIcon("icon_tabs.png");
      case csDockPaneManager.LAYOUT_ONE_ROW:
        return cseis.resources.csResources.getIcon("icon_1row.png");
      case csDockPaneManager.LAYOUT_ONE_COL:
        return cseis.resources.csResources.getIcon("icon_1col.png");
      case csDockPaneManager.LAYOUT_TWO_ROWS:
        return cseis.resources.csResources.getIcon("icon_2rows.png");
      case csDockPaneManager.LAYOUT_TWO_COLS:
        return cseis.resources.csResources.getIcon("icon_2cols.png");
      case csDockPaneManager.LAYOUT_THREE_ROWS:
        return cseis.resources.csResources.getIcon("icon_3rows.png");
      case csDockPaneManager.LAYOUT_THREE_COLS:
        return cseis.resources.csResources.getIcon("icon_3cols.png");
      default:
        return null;
    }
  }
  public static String getDescription( int layoutOption ) {
    switch( layoutOption ) {
      case csDockPaneManager.LAYOUT_TABS:
        return "Use tabbed panes";
      case csDockPaneManager.LAYOUT_ONE_ROW:
        return "Arrange panes in one row";
      case csDockPaneManager.LAYOUT_ONE_COL:
        return "Arrange panes in one column";
      case csDockPaneManager.LAYOUT_TWO_ROWS:
        return "Arrange panes in two rows";
      case csDockPaneManager.LAYOUT_TWO_COLS:
        return "Arrange panes in two columns";
      case csDockPaneManager.LAYOUT_THREE_ROWS:
        return "Arrange panes in three rows";
      case csDockPaneManager.LAYOUT_THREE_COLS:
        return "Arrange panes in three columns";
      default:
        return null;
    }
  }
  public void addDockPaneSelectionListener( csIDockPaneSelectionListener listener ) {
    myListeners.add( listener );
  }
  public void removeDockPaneSelectionListener( csIDockPaneSelectionListener listener ) {
    myListeners.remove( listener );
  }
  public csDockPane addPanel( JPanel panel, String title ) {
    return addPanel( panel, title, new csDockPaneButtonSelection() );
  }
  public csDockPane addPanel( JPanel panel, String title, csDockPaneButtonSelection bselect ) {
    return addPanel( panel, title, bselect, title );
  }
  public csDockPane addPanel( JPanel panel, String title, csDockPaneButtonSelection bselect, String titleDescription ) {
    csDockPane dockPane = new csDockPane( this, panel, title, bselect, titleDescription );
    addDockPane( dockPane, false );
    return dockPane;
  }
  public csDockPane addDockPane( csDockPane pane ) {
    addDockPane( pane, true );
    return pane;
  }
  public void addDockPane( csDockPane pane, boolean refresh ) {
    int paneIndex = getPaneIndex( pane );
    if( paneIndex >= 0 ) return;  // Pane already exists
    pane.addPaneListener(this);
    myPanes.add( pane );  // Add pane to list of panes
    myIsNewPaneInProgress = true;
    if( pane.isVisible() ) {
      myNumPanesVisible += 1;
      pane.setState( csDockPaneManager.STATE_DOCKED );
      addToMainPanel( pane );  // Add pane to manager (show pane window)
      if( myDoActivateNewPanes ) setActivePane( pane );
    }
    else {
      pane.setState( csDockPaneManager.STATE_HIDDEN );
    }
    resetWindowLayout( myLayoutOption, refresh );
    myIsNewPaneInProgress = false;
  }
  private void addBackPane( csDockPane pane, int paneIndex ) {
    int counterDockedBeforePane = 0;
    // Count docked panes in front of this pane, place pane at end
    // For example, if all panes in front of this pane are hidden or undocked, place this pane in first place
    for( int i = 0; i < paneIndex; i++ ) {
      if( myPanes.get(i).getState() == STATE_DOCKED ) {
        counterDockedBeforePane += 1;
      }
    }
    pane.setState( csDockPaneManager.STATE_DOCKED );
    pane.setVisible(true);
    addToMainPanel( pane, counterDockedBeforePane );
    myNumPanesVisible += 1;
    reset();
    if( myActivePaneIndex < 0 || myActivePaneIndex >= myPanes.size() ) {
      if( myDoActivateNewPanes ) setActivePane( myPanes.get(paneIndex) );
    }
  }
  private void addToMainPanel( csDockPane pane ) {
    addToMainPanel( pane, getComponentCount() );
  }
  private void addToMainPanel( csDockPane pane, int atIndex ) {
    if( myLayoutOption != LAYOUT_TABS ) {
      add( pane, atIndex );
    }
    else {
      insertTab( pane, atIndex-1 );
    }
  }
  private void insertTab( csDockPane pane, int atIndex ) {
    myTabbedPane.insertTab( "", null, pane, "", atIndex );
    myTabbedPane.setTabComponentAt( atIndex, pane.getTabComponent() );
    myTabbedPane.setToolTipTextAt( atIndex, pane.attr().titleDescription );
  }
  public void resetActiveTitle( String title, String titleDescription ) {
    csDockPane activePane = getActiveDockPane();
    if( activePane != null ) {
      activePane.setTitle( title, titleDescription );
      if( myLayoutOption == csDockPaneManager.LAYOUT_TABS ) {
        int index = myTabbedPane.getSelectedIndex();
        myTabbedPane.setTabComponentAt( index, activePane.getTabComponent() );
      }
    }
  }
  //----------------------------------------------------------------
  //
  public void setPaneVisible( JPanel pane, boolean doSetVisible ) {
    int paneIndex = getPaneIndex( pane );
    setPaneVisible( paneIndex, doSetVisible, true );
  }
  private void setPaneVisible( int paneIndex, boolean doSetVisible, boolean refresh ) {
    if( paneIndex >= 0 && paneIndex < myPanes.size() ) {
      csDockPane pane = myPanes.get(paneIndex);
      int stateCurrent = pane.getState();
      if( stateCurrent == csDockPaneManager.STATE_HIDDEN ) {
        if( doSetVisible ) {
          addBackPane(pane,paneIndex);
        }
      }
      else if( !doSetVisible ) {
        myNumPanesVisible -= 1;
        pane.setState( csDockPaneManager.STATE_HIDDEN );
        if( stateCurrent == STATE_DOCKED ) {
          remove(pane);
        }
        else {
          pane.hidePane();
        }
      }
      pane.setVisible( doSetVisible );
      if( refresh ) reset();
    } // END if paneIndex in range
  }
  private void showHideAll( boolean doShowAll ) {
    if( doShowAll ) {
      for( int i = 0; i < myPanes.size(); i++ ) {
        csDockPane pane = myPanes.get(i);
        if( pane.getState() == STATE_HIDDEN ) {
          addBackPane(pane,i);
        }
      }
    }
    else {
      for( int i = 0; i < myPanes.size(); i++ ) {
        csDockPane pane = myPanes.get(i);
        if( pane.getState() != STATE_HIDDEN ) {
          pane.hidePane();
          if( pane.getState() == STATE_DOCKED ) {
            remove(pane);
          }
          pane.setState(STATE_HIDDEN);
          pane.setVisible(false);
        }
      }
    }
    if( doShowAll ) myNumPanesVisible = myPanes.size();
    else myNumPanesVisible = 0;
//    fireHideStateEvent( myNumPanesVisible != myPanes.size() );
    reset();
  }
  public int getState( JPanel mainPanel ) {
    int paneIndex = getMainPanelIndex( mainPanel );
    if( paneIndex >= 0 && paneIndex < myPanes.size() ) {
      return myPanes.get(paneIndex).getState();
    }
    return -1;
  }
  public void changeDocking( JPanel pane, boolean refresh ) {
    int paneIndex = getPaneIndex( pane );
    if( paneIndex >= 0 && paneIndex < myPanes.size() ) {
      changeDocking( paneIndex, true );
    }
  }
  private void changeDocking( int paneIndex, boolean refresh ) {
    csDockPane pane = myPanes.get(paneIndex);
    int stateCurrent = pane.getState();
    if( stateCurrent == STATE_UNDOCKED ) {
      pane.setState(STATE_DOCKED);
      addBackPane(pane,paneIndex);
    }
    else {
      remove(pane);
      csDockFrame frame = new csDockFrame(pane,this);
      frame.setVisible(true);
      pane.setState(STATE_UNDOCKED);
      myNumPanesVisible -= 1;
    }
//    fireHideStateEvent( myNumPanesVisible != myPanes.size() );
    if( refresh ) reset();
  }
  public void reset() {
    resetWindowLayout( myLayoutOption, true );     
  }
  public int maxLayoutOption() {
    if( myNumPanesVisible < 2 ) return LAYOUT_ONE_COL;
    else if( myNumPanesVisible < 3 ) return LAYOUT_TWO_ROWS;
    else return LAYOUT_THREE_COLS;
  }
  public int getWindowLayout() {
    return myLayoutOption;
  }
  private void resetTabbedPane( int newLayoutOption ) {
    csDockPane activePane = getActiveDockPane();
    if( newLayoutOption == LAYOUT_TABS ) {
      myTabbedPane = new JTabbedPane( JTabbedPane.TOP, JTabbedPane.SCROLL_TAB_LAYOUT );
      myTabbedPane.addChangeListener( new ChangeListener() {
        @Override
        public void stateChanged( ChangeEvent e ) {
          if( !myIsNewPaneInProgress ) {
            csDockPane pane = (csDockPane)myTabbedPane.getSelectedComponent();
            if( pane != null && !pane.attr().isActive ) setActivePane( pane );
          }
        }
      });
      for( int i = 0; i < myPanes.size(); i++ ) {
        csDockPane pane = myPanes.get(i);
        pane.setState( STATE_DOCKED );
        pane.setTopBarVisible(false);
        insertTab( pane, i );
      }
      removeAll();
      if( activePane != null ) setActivePane( activePane, newLayoutOption );
      add( myTabbedPane );
      myNumPanesVisible = myPanes.size();
      revalidate();
    }
    else {
      removeAll();
      for( int i = 0; i < myPanes.size(); i++ ) {
        csDockPane pane = myPanes.get(i);
        pane.setState( STATE_DOCKED );
        pane.setTopBarVisible(true);
        add( pane );
      }
      myTabbedPane = null;
      myNumPanesVisible = myPanes.size();
      if( activePane != null ) setActivePane( activePane, newLayoutOption );
      revalidate();
    }
  }
  public int resetWindowLayout( int layoutOption, boolean refresh ) {
    int numRows = 1;
    int numCols = 1;
    if( layoutOption == LAYOUT_TABS || myLayoutOption == LAYOUT_TABS ) {
      resetTabbedPane( layoutOption );
    }
    myLayoutOption = layoutOption;
    if( myPanes.isEmpty() ) {
      // Nothing
    }
    else if( myLayoutOption == LAYOUT_TABS ) {
      if( myTabbedPane == null ) {
        myTabbedPane = new JTabbedPane( JTabbedPane.TOP, JTabbedPane.SCROLL_TAB_LAYOUT );
        removeAll();
        add( myTabbedPane );
      }
      numRows = 1;
      numCols = 1;
    }
    else if( myLayoutOption == LAYOUT_ONE_ROW ) {
      numRows = 1;
      numCols = myNumPanesVisible / numRows;
    }
    else if( myLayoutOption == LAYOUT_ONE_COL ) {
      numCols = 1;
      numRows = myNumPanesVisible / numCols;
    }
    else if( myLayoutOption == LAYOUT_TWO_ROWS ) {
      numRows = Math.max(1,Math.min(myNumPanesVisible,2));
      numCols = myNumPanesVisible / numRows;
      if( numRows * numCols < myNumPanesVisible ) numCols += 1;
    }
    else if( myLayoutOption == LAYOUT_TWO_COLS ) {
      numCols = Math.max(1,Math.min(myNumPanesVisible,2));
      numRows = myNumPanesVisible / numCols;
      if( numRows * numCols < myNumPanesVisible ) numRows += 1;
    }
    else if( myLayoutOption == LAYOUT_THREE_ROWS ) {
      numRows = Math.max(1,Math.min(myNumPanesVisible,3));
      numCols = myNumPanesVisible / numRows;
      if( numRows * numCols < myNumPanesVisible ) {
        numCols += 1;
        if( (numRows-1) * numCols >= myNumPanesVisible ) numRows -= 1;
      }
    }
    else if( myLayoutOption == LAYOUT_THREE_COLS ) {
      numCols = Math.max(1,Math.min(myNumPanesVisible,3));
      numRows = myNumPanesVisible / numCols;
      if( numRows * numCols < myNumPanesVisible ) numRows += 1;
    }
    else {
      numRows = 1;
      numCols = myNumPanesVisible / numRows;
    }
    myGridLayout.setRows(numRows);
    myGridLayout.setColumns(numCols);
    if( refresh ) {
      revalidate();
      repaint();
    }
    return myLayoutOption;
  }
  public int numAllPanes() {
    return myPanes.size();
  }
  public int numVisiblePanes() {
    return myNumPanesVisible;
  }
  public JPanel createSelectionPanel() {
    JPanel panel = new JPanel(new GridBagLayout());
    JCheckBox checkAll = new JCheckBox("Select all",false);
    final PaneCheckBox[] buttons = new PaneCheckBox[myPanes.size()];
    panel.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Panel selection"),
        BorderFactory.createEmptyBorder(2, 2, 2, 2) ) );
    int yp = 0;
    panel.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.5, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panel.add( checkAll, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 0 ), 0, 0 ) );
    for( int i = 0; i < myPanes.size(); i++ ) {
      csDockPane pane = myPanes.get(i);
      buttons[i] = new PaneCheckBox(pane);
      panel.add( buttons[i], new GridBagConstraints(
          0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      buttons[i].addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          PaneCheckBox button = (PaneCheckBox)e.getSource();
          csDockPane pane = (csDockPane)button.pane;
          if( button.isSelected() ) {
            pane.showPane();
            addBackPane( pane, getPaneIndex(pane) );
          }
          else {
            pane.hidePane();
            setPaneVisible( pane, false );
          }
//          fireHideStateEvent( myNumPanesVisible != myPanes.size() );
        }
      });
    }
    panel.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.5, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    checkAll.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        boolean isSelected = ((JCheckBox)e.getSource()).isSelected();
        showHideAll( isSelected );
        for( int i = 0; i < myPanes.size(); i++ ) {
          buttons[i].setSelected(isSelected);
        }
      }
    });

    return panel;
  }
  private int getPaneIndex( JPanel pane_in ) {
    for( int i = 0; i < myPanes.size(); i++ ) {
      csDockPane pane = myPanes.get(i);
      if( pane.equals(pane_in) ) {
        return i;
      }
    }
    return -1;
  }
  private int getMainPanelIndex( JPanel mainPanel_in ) {
    for( int i = 0; i < myPanes.size(); i++ ) {
      csDockPane pane = myPanes.get(i);
      if( pane.getMainPanel().equals(mainPanel_in) ) {
        return i;
      }
    }
    return -1;
  }
  public csDockPane getDockPane( int index ) {
    if( index < 0 || index >= myPanes.size() ) {
      return null;
    }
    return myPanes.get(index);
  }
  public csDockPane getActiveDockPane() {
    if( myActivePaneIndex < 0 || myActivePaneIndex >= myPanes.size() ) {
      myActivePaneIndex = -1;
      return null;
    }
    return myPanes.get(myActivePaneIndex);
  }
  public JPanel getActivePanel() {
    if( myActivePaneIndex < 0 || myActivePaneIndex >= myPanes.size() ) {
      myActivePaneIndex = -1;
      return null;
    }
    return myPanes.get(myActivePaneIndex).getMainPanel();
  }
  public void setActivePane( csDockPane pane ) {
    setActivePane( pane, myLayoutOption );
    fireDockPaneSelectionEvent( pane );
  }
  public void setActivePane( csDockPane pane, int layoutOption ) {
    int paneIndex = getPaneIndex(pane);
    if( paneIndex < 0 ) return;
    for( int i = 0; i < myPanes.size(); i++ ) {
      if( i != paneIndex && myPanes.get(i).attr().isActive ) {
        myPanes.get(i).setActive(false);
      }
    }
    pane.setActive(true);
    myActivePaneIndex = paneIndex;
    if( layoutOption == LAYOUT_TABS ) myTabbedPane.setSelectedComponent(pane);
  }
  public void resetActivePane() {
    myActivePaneIndex = -1;
    // Set first visible pane as active
    for( int i = 0; i < myPanes.size(); i++ ) {
      if( myPanes.get(i).isVisible() ) {
        setActivePane( myPanes.get(i), myLayoutOption );
        fireDockPaneSelectionEvent( myPanes.get(i) );
        return;
      }
    }
    fireDockPaneSelectionEvent( null );
  }
  public void pressPane( csDockPane pane, Point point ) {
    int paneIndex = getPaneIndex(pane);
    if( paneIndex < 0 ) return;
    setActivePane( pane );
    if( pane.getState() == STATE_UNDOCKED ) return;
    myIsDragging = true;
    myDragPoint = point;
    myDragPointOrigin = point;
    myDragPaneIndex = paneIndex;
    if( myShowActivePaneFeedback ) {
      int width  = pane.getWidth();
      int height = pane.getHeight();
      myDragImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
      Graphics2D g2 = myDragImage.createGraphics();
      pane.paint(g2);
      float factor = 255/5; // Brighten by 20%
      for( int x = 0; x < width; x++ ) {
        for( int y = 0; y < height; y++ ) {
          int rgb = myDragImage.getRGB(x, y);
          int red = (rgb >> 16) & 0xFF;
          int green = (rgb >> 8) & 0xFF;
          int blue = rgb & 0xFF;
          int redNew   = (int)Math.min(255,red + factor);
          int greenNew = (int)Math.min(255,green + factor);
          int blueNew  = (int)Math.min(255,blue + factor);
          Color colorNew = new Color(redNew, greenNew, blueNew);
          myDragImage.setRGB(x, y, colorNew.getRGB());
        }
      }
      repaint();
    }
  }
  public void releasePane( csDockPane pane ) {
    if( !myIsDragging ) return;
    int paneIndex = getPaneIndex(pane);
    if( paneIndex < 0 ) return;
    myIsDragging = false;
    Dimension totalManagerSize = getSize();
    int numCols = Math.max( myGridLayout.getColumns(), 1 );
    int numRows = Math.max( myGridLayout.getRows(), 1 );
    int paneWidth  = totalManagerSize.width / numCols;
    int paneHeight = totalManagerSize.height / numRows;
    Point pOnScreen = getLocationOnScreen();
    int colCurrent = (int)( (myDragPoint.x-pOnScreen.x) / paneWidth );
    int rowCurrent = (int)( (myDragPoint.y-pOnScreen.y) / paneHeight );
    int colOrigin = (int)( (myDragPointOrigin.x-pOnScreen.x) / paneWidth );
    int rowOrigin = (int)( (myDragPointOrigin.y-pOnScreen.y) / paneHeight );
    if( colOrigin != colCurrent || rowOrigin != rowCurrent ) {
      int paneIndexNewOnScreen = Math.min(rowCurrent * numCols + colCurrent, myNumPanesVisible-1);
      csDockPane paneOrig = (csDockPane)getComponent(paneIndexNewOnScreen);
      int paneIndexNew = getPaneIndex(paneOrig);
      myPanes.remove(paneIndex);
      myPanes.add(paneIndexNew, pane);
      removeAll();
      for( int i = 0; i < myPanes.size(); i++ ) {
        csDockPane paneCurrent = myPanes.get(i);
        if( paneCurrent.isVisible() && paneCurrent.attr().isDocked ) {
          addToMainPanel( paneCurrent );
        }
      }
      if( myActivePaneIndex == paneIndex ) { // Should always be the case
        myActivePaneIndex = paneIndexNew;
      }
      reset();
    }
    repaint();
  }
  public void dragPane( csDockPane pane, Point point ) {
    int paneIndex = getPaneIndex(pane);
    if( paneIndex < 0 ) return;
    myIsDragging = true;
    myDragPoint = point;
    myDragPaneIndex = paneIndex;
    repaint();
  }
  @Override
  public void paint( Graphics g ) {
    super.paint(g);
    if( myIsDragging && myPanes.get(myDragPaneIndex).getState() == STATE_DOCKED ) {
      Graphics2D g2 = (Graphics2D)g;
      Point pOnScreen = getLocationOnScreen();
      drawDragImage( g2, pOnScreen );
    }
  }
  protected void drawDragImage( Graphics2D g2,Point pOnScreen ) {
    int numCols = Math.max( myGridLayout.getColumns(), 1 );
    int numRows = Math.max( myGridLayout.getRows(), 1 );
    Dimension totalManagerSize = getSize();
    int paneWidth  = totalManagerSize.width / numCols;
    int paneHeight = totalManagerSize.height / numRows;
    int colCurrent = (int)( (myDragPoint.x-pOnScreen.x) / paneWidth );
    int rowCurrent = (int)( (myDragPoint.y-pOnScreen.y) / paneHeight );
    int colOrigin = (int)( (myDragPointOrigin.x-pOnScreen.x) / paneWidth );
    int rowOrigin = (int)( (myDragPointOrigin.y-pOnScreen.y) / paneHeight );
    Stroke stroke = g2.getStroke();
    Color colorSave = g2.getColor();
    if( colOrigin != colCurrent || rowOrigin != rowCurrent ) {
      int xPaneOrigin = colOrigin*paneWidth;
      int yPaneOrigin = rowOrigin*paneHeight;
      g2.setColor( new JLabel().getBackground() );
      g2.fillRect(xPaneOrigin, yPaneOrigin, paneWidth, paneHeight);
    }
    int xPaneCurrent = colCurrent*paneWidth;
    int yPaneCurrent = rowCurrent*paneHeight;

    g2.drawImage(myDragImage, xPaneCurrent, yPaneCurrent, null);
    g2.setStroke( new BasicStroke(3.0f) );
    g2.setColor( Color.black );
    g2.drawRect(xPaneCurrent, yPaneCurrent, paneWidth, paneHeight);
    g2.setStroke( stroke );
    g2.setColor(colorSave);
  }
  //------------------------------------------------------------------------
  // Implementation of csIDockPaneListener interface:
  @Override
  public void hidePane(csDockPane pane) {
    setPaneVisible( pane, false );
    int paneIndex = getPaneIndex( pane );
    if( paneIndex == myActivePaneIndex ) {
      pane.setActive(false);
      resetActivePane();
    }
  }
  @Override
  public void closePane( csDockPane pane ) {
    int paneIndex = getPaneIndex( pane );
    if( paneIndex >= 0 ) {
      int stateCurrent = pane.getState();
      if( stateCurrent != STATE_HIDDEN ) {
        myNumPanesVisible -= 1;
        if( stateCurrent == STATE_DOCKED ) {
          remove(pane); // Remove pane from view
        }
      }
      myPanes.remove(paneIndex); // Completely remove pane from manager
      if( paneIndex == myActivePaneIndex ) {
        resetActivePane();
      }
      else if( paneIndex < myActivePaneIndex ) {
        myActivePaneIndex -= 1;
      }
      resetWindowLayout( myLayoutOption, true );
      fireDockPaneClosedEvent( pane );
    }
  }
  @Override
  public void dockPane( csDockPane pane, boolean dock ) {
    changeDocking( pane, true );
  }
  @Override
  public void scrollbars(boolean showScrollbars) {
  }
  @Override
  public void maximizePane( csDockPane paneIn ) {
    int paneIndex = getPaneIndex( paneIn );
    if( paneIndex >= 0 ) {
      setPaneVisible( paneIndex, true, false );
      for( int i = 0; i < myPanes.size(); i++ ) {
        if( i == paneIndex ) continue;
        csDockPane pane = myPanes.get(i);
        if( pane.getState() != STATE_HIDDEN ) {
          pane.hidePane();
          if( pane.getState() == STATE_DOCKED ) {
            remove(pane);
          }
          pane.setState(STATE_HIDDEN);
          pane.setVisible(false);
        }
      }
      myNumPanesVisible = 1;
      reset();
    }
  }
  private void fireDockPaneSelectionEvent( csDockPane pane ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).dockPaneSelected( pane );
    }
  }
  private void fireDockPaneClosedEvent( csDockPane pane ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).dockPaneClosed( pane );
    }
  }
  //*******************************************************************
  private class PaneCheckBox extends JCheckBox {
    PaneCheckBox( csDockPane pane_in ) {
      super( pane_in.attr().title, pane_in.isVisible() );
      pane = pane_in;
    }
    public JPanel pane;
  }
  //*******************************************************************
  // TEST functions
  //
  public void setNumPanes( int numPanes ) {
    if( numPanes <= myPanes.size() ) {
      removeAll();
      for( int i = 0; i < numPanes; i++ ) {
        addToMainPanel( myPanes.get(i) );
      }
      myNumPanesVisible = numPanes;
      resetWindowLayout( myLayoutOption, true );
    }
  }
	    
}


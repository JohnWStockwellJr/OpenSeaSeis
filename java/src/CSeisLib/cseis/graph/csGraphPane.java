/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

import javax.swing.*;
import java.awt.event.*;
import java.awt.*;

/**
 * Graph pane provides scrollbars and other functionality for csGraph2D<br>
 * @author Bjorn Olofsson
 */
// ...this class is a work in progress
@SuppressWarnings("serial")
public class csGraphPane extends JPanel implements csIGraph2DListener {
  private JScrollBar myScrollBarVert;
  private JScrollBar myScrollBarHorz;
  private int myPreviousScrollValueVert;
  private int myPreviousScrollValueHorz;

  private boolean myStopScrolling = false;
  private csGraph2D myGraph;
  private JPanel myPanelViewPort = null;

  private Dimension myGraphSize = null;
  private JLabel myStatusBar;

  public csGraphPane( csGraph2D graph ) {
    this( graph, true, true );
  }
  public csGraphPane( csGraph2D graph, boolean verticalScrollBar, boolean horizontalScrollBar ) {
    super( new BorderLayout() );
    myGraph = graph;
    graph.addGraph2DListener(this);

    myPreviousScrollValueVert = 0;
    myPreviousScrollValueHorz = 0;

    if( horizontalScrollBar ) {
      myScrollBarHorz = new JScrollBar( JScrollBar.HORIZONTAL );
      myScrollBarHorz.setValues( 0, 100, 0, 0 );
      myScrollBarHorz.setBlockIncrement(100);
      myScrollBarHorz.setUnitIncrement(4);
      myScrollBarHorz.setVisibleAmount( 100 );
      myScrollBarHorz.addAdjustmentListener( new AdjustmentListener() {
        public void adjustmentValueChanged( AdjustmentEvent e ) {
          if( myStopScrolling ) return;
          int valueNew = e.getValue();
          if( myPreviousScrollValueHorz != valueNew ) {
            myPreviousScrollValueHorz = valueNew;
            myGraph.resetViewPositionHorz( valueNew );
          }
        }
      });
    }
    if( verticalScrollBar ) {
      myScrollBarVert = new JScrollBar( JScrollBar.VERTICAL );
      myScrollBarVert.setValues( 0, 100, 0, 0 );
      myScrollBarVert.setBlockIncrement(100);
      myScrollBarVert.setUnitIncrement(4);
      myScrollBarVert.setVisibleAmount( 100 );
      myScrollBarVert.addAdjustmentListener( new AdjustmentListener() {
        public void adjustmentValueChanged( AdjustmentEvent e ) {
          if( myStopScrolling ) return;
          int valueNew = e.getValue();
          if( myPreviousScrollValueVert != valueNew ) {
            myPreviousScrollValueVert = valueNew;
            myGraph.resetViewPositionVert( valueNew );
          }
        }
      });
    }
    myStatusBar = new JLabel("----", JLabel.LEFT );

    // GridBagLayout necessary to place view into upper left corner of scroll pane's view port.
    myPanelViewPort = new JPanel( new GridBagLayout() );
    myPanelViewPort.add( graph, new GridBagConstraints(
        0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myPanelViewPort.add( Box.createHorizontalGlue(), new GridBagConstraints(
        1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myPanelViewPort.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, 1, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelScrollPane = new JPanel( new BorderLayout() );
    //---------------------------------------------
    if( verticalScrollBar ) {
      panelScrollPane.add( myScrollBarVert, BorderLayout.EAST );
    }
    if( horizontalScrollBar ) {
      panelScrollPane.add( myScrollBarHorz, BorderLayout.SOUTH );
    }
    panelScrollPane.add( myPanelViewPort, BorderLayout.CENTER );

    add( panelScrollPane, BorderLayout.CENTER );
    add( myStatusBar, BorderLayout.SOUTH );

    setMinimumSize( new Dimension(0,0) );


    addComponentListener( new ComponentAdapter() {
      public void componentResized(ComponentEvent e ) {
        setNewSize();
        myGraph.resetViewPositionVert( myScrollBarVert.getValue() );
        myGraph.resetViewPositionHorz( myScrollBarHorz.getValue() );
      }
    });
  }
  public void setNewSize() {
    Dimension sizeViewPortPanel = myPanelViewPort.getSize();
    int newHeight = 0;
    int newWidth  = 0;
    Dimension sizeGraph = myGraph.getPreferredSize();
    if( myGraph.isFixedXDim() ) {
      newWidth = sizeGraph.width;
    }
    else {
      newWidth = sizeViewPortPanel.width;
    }
    if( myGraph.isFixedYDim() ) {
      newHeight = sizeGraph.height;
    }
    else {
      newHeight = sizeViewPortPanel.height;
    }
    myGraphSize = new Dimension(newWidth,newHeight);
    if( myGraphSize.height != sizeGraph.height || myGraphSize.width != sizeGraph.width ) {
      myGraph.setNewVariableDim( newWidth, newHeight );
    }

    int maxScrollValueVert = Math.max( myGraphSize.height - sizeViewPortPanel.height, 0 );
    int scrollValueVert    = myGraph.getViewPositionVert();
    myPreviousScrollValueVert = -1;

    int maxScrollValueHorz = Math.max( myGraphSize.width  - sizeViewPortPanel.width, 0 );
    int scrollValueHorz    = myGraph.getViewPositionHorz();
    myPreviousScrollValueHorz = -1;

    // Complex coding to avoid refreshing of scroll bars and seismic while setting new maximum values etc. ...not desirable.
    myStopScrolling = true;
    myScrollBarVert.setBlockIncrement( (myGraphSize.height/2)+1 );
    myScrollBarVert.setMaximum( maxScrollValueVert );
    myScrollBarHorz.setBlockIncrement( (myGraphSize.width/2)+1 );
    myScrollBarHorz.setMaximum( maxScrollValueHorz );
    myStopScrolling = false;

    myScrollBarVert.setValue( scrollValueVert <= maxScrollValueVert ? scrollValueVert : maxScrollValueVert );
    myScrollBarHorz.setValue( scrollValueHorz <= maxScrollValueHorz ? scrollValueHorz : maxScrollValueHorz );
  }

  @Override
  public void graph2DValues(float xModel, float yModel) {
    myStatusBar.setText("X value: " + xModel + "  Y value: " + yModel);
  }
  public void graphChanged() {
  }
}



/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.*;

import java.util.ArrayList;

import cseis.seis.csHeaderDef;
import cseis.general.csStandard;

/**
 * Dialog window where trace headers are selected for trace annotation.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csAnnotationDialog extends JDialog {
  private static final int DEFAULT_WIDTH  = 400;
  private static final int DEFAULT_HEIGHT = 850;
  
  private csSeisPane mySeisPane;
  /// Trace header names & descriptions
  private ArrayList<csHeaderDef> myTraceHeaders;
  private JList<csHeaderDef> myListHeaders;
  /// List model containing trace headers, except those that shall be annotated
  private DefaultListModel<csHeaderDef> myListModelHeaders;
  private JList<csHeaderDef> myListAnn;
  /// List model containing trace headers that shall be annotated
  private DefaultListModel<csHeaderDef> myListModelAnn;
  private ArrayList<csHeaderDef> myTempAddHeaders;
  private ArrayList<csHeaderDef> myTempRemoveHeaders;

  private JCheckBox myBoxShowSeqTraceNum;
  /// Omit plotting repeated header values
  private JCheckBox myBoxOmitRepeatedValues;

  private JCheckBox myBoxTraceLabelStep;
  private JTextField myTextTraceLabelStep;

  private JButton myButtonAddHeader;
  private JButton myButtonRemoveHeader;
  
  private JButton myButtonOK;
  private JButton myButtonApply;
  private JButton myButtonCancel;

  public csAnnotationDialog( JFrame frame, csSeisPane seisPane, csHeaderDef[] headers ) {
    super( frame, "Set trace annotation");
    super.setModal( false );

    myTempAddHeaders = new ArrayList<csHeaderDef>();
    myTempRemoveHeaders = new ArrayList<csHeaderDef>();
    mySeisPane = seisPane;
    myBoxShowSeqTraceNum = new JCheckBox( "Show sequential trace number", true );
    myBoxShowSeqTraceNum.setToolTipText( "Select this box to plot a sequential trace number of every trace" );
    myBoxOmitRepeatedValues   = new JCheckBox( "Omit repeating header values", true );
    myBoxOmitRepeatedValues.setToolTipText( "Select this box to plot header values only when they change. Deselect to plot values for every trace" );

    myBoxTraceLabelStep = new JCheckBox( "Fixed trace label step", false );
    myBoxTraceLabelStep.setToolTipText("Select to set a fixed trace step for trace label");
    myTextTraceLabelStep = new JTextField("10");
    myTextTraceLabelStep.setEnabled( myBoxTraceLabelStep.isSelected() );

    myListModelHeaders = new DefaultListModel<csHeaderDef>();
    myListHeaders = new JList<csHeaderDef>( myListModelHeaders );
    myListHeaders.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    myListHeaders.setCellRenderer( new AnnotationCellRenderer() );
    myListModelAnn = new DefaultListModel<csHeaderDef>();
    myListAnn = new JList<csHeaderDef>( myListModelAnn );

    updateHeaderNames( headers );

    myListAnn.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    myListAnn.setCellRenderer( new AnnotationCellRenderer() );

    myButtonAddHeader = new JButton("+");
    myButtonRemoveHeader = new JButton("-");
    Dimension size = myButtonAddHeader.getPreferredSize();
    int max = Math.max( size.height, size.width );
    myButtonAddHeader.setPreferredSize(new Dimension(max,max));
    myButtonAddHeader.setMinimumSize(new Dimension(max,max));
    myButtonRemoveHeader.setPreferredSize(new Dimension(max,max));
    myButtonRemoveHeader.setMinimumSize(new Dimension(max,max));
    
    myButtonOK = new JButton("OK");
    myButtonApply = new JButton("Apply");
    myButtonCancel = new JButton("Cancel");
    this.getRootPane().setDefaultButton(myButtonApply);
    
    JPanel panelHeaders = new JPanel(new BorderLayout());
    panelHeaders.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Trace headers"),
        csStandard.INNER_EMPTY_BORDER ) );
    panelHeaders.add( new JScrollPane(myListHeaders), BorderLayout.CENTER );

    JPanel panelButtons = new JPanel(new GridBagLayout());
    panelButtons.setBorder( BorderFactory.createEmptyBorder(3, 3, 3, 3) );
    int xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.5, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonAddHeader, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonRemoveHeader, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 11, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.5, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelAnn = new JPanel(new BorderLayout());
    panelAnn.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Annotation headers"),
        csStandard.INNER_EMPTY_BORDER ) );
    panelAnn.add( new JScrollPane(myListAnn), BorderLayout.CENTER );

    JPanel panelApplyButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelApplyButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelApplyButtons.add( myButtonOK, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelApplyButtons.add( myButtonApply, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelApplyButtons.add( myButtonCancel, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelApplyButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelUpper = new JPanel(new GridBagLayout());
    JPanel panelLower = new JPanel(new GridBagLayout());
    int yp = 0;
    panelUpper.add( panelHeaders, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelUpper.add( panelButtons, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    yp = 0;
    panelLower.add( panelAnn, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelLower.add( myBoxShowSeqTraceNum, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelLower.add( myBoxOmitRepeatedValues, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelLower.add( myBoxTraceLabelStep, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.VERTICAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelLower.add( myTextTraceLabelStep, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 5, 0, 10 ), 0, 0 ) );

    JSplitPane splitPane = new JSplitPane( JSplitPane.VERTICAL_SPLIT, panelUpper, panelLower );
    splitPane.setDividerSize( 3 );
    
    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );
    panelAll.add( splitPane, BorderLayout.CENTER );
    panelAll.add( panelApplyButtons, BorderLayout.SOUTH );

    this.getContentPane().add( panelAll, BorderLayout.CENTER );
    setSize( DEFAULT_WIDTH, DEFAULT_HEIGHT );
    pack();
    splitPane.setDividerLocation(0.70);

    myButtonOK.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        apply();
        dispose();
      }
    });
    myButtonApply.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    myButtonCancel.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        cancel();
      }
    });

    myButtonAddHeader.addActionListener( new ActionListener() {
      public void actionPerformed( ActionEvent e ) {
        addAnnotationHeaders();
      }
    } );
    myButtonRemoveHeader.addActionListener( new ActionListener() {
      public void actionPerformed( ActionEvent e ) {
        removeAnnotationHeaders();
      }
    } );
    myListHeaders.addMouseListener( new MouseAdapter() {
      public void mouseClicked( MouseEvent e ) {
        if( e.getClickCount() == 2 ) {
//          int index = myListHeaders.locationToIndex(e.getPoint());
          addAnnotationHeaders();
        }
      }
    });
    myListAnn.addMouseListener( new MouseAdapter() {
      public void mouseClicked( MouseEvent e ) {
        if( e.getClickCount() == 2 ) {
//          int index = myListAnn.locationToIndex(e.getPoint());
          removeAnnotationHeaders();
        }
      }
    });
    myBoxTraceLabelStep.addActionListener( new ActionListener() {
      public void actionPerformed( ActionEvent e ) {
        myTextTraceLabelStep.setEnabled( myBoxTraceLabelStep.isSelected() );
      }
    } );
    
    addWindowListener(new WindowAdapter() {
      public void windowClosing(WindowEvent e) {
        cancel();
      }
    });
  }
  /**
   * Update trace headers
   * @param headers Array of trace header objects
   */
  public void updateHeaderNames( csHeaderDef[] headers ) {
    if( myTempAddHeaders.isEmpty() || myTempRemoveHeaders.isEmpty() ) {
      cancel();
    }

    int numHeaders = headers.length;
    myTraceHeaders = new ArrayList<csHeaderDef>(numHeaders);

    for( int ih = 0; ih < numHeaders; ih++ ) {
      myTraceHeaders.add( new csHeaderDef( headers[ih] ) );
    }
    myListModelHeaders.removeAllElements();
    java.util.Collections.sort(myTraceHeaders);
    for( int ih = 0; ih < numHeaders; ih++ ) {
      myListModelHeaders.addElement( myTraceHeaders.get(ih) );
    }

    // Special care for currently annotated headers after update with new seismic data set.
    // If annotated trace headers exist, check if new seismic data set contains same trace headers:
    // If yes, keep this trace header annotated. If not, remove it from annotated list
    int numHeadersAnnotated = myListModelAnn.size();
    if( numHeadersAnnotated > 0 ) {
      for( int ih = 0; ih < numHeadersAnnotated; ih++ ) {
        csHeaderDef header = (csHeaderDef)myListModelAnn.get(ih);
        if( myTraceHeaders.contains(header) ) {
          myListModelHeaders.removeElement( header );
        }
        else {
          myTempRemoveHeaders.add( header );
        }
      }
      for( int i = 0; i < myTempRemoveHeaders.size(); i++ ) {
        myListModelAnn.removeElement( myTempRemoveHeaders.get(i) );
      }
      myTempRemoveHeaders.clear();
      apply();
    }
  }
  private void addAnnotationHeaders() {
    if( !myListHeaders.isSelectionEmpty() ) {
      java.util.List<csHeaderDef> objList = myListHeaders.getSelectedValuesList();
      for( int i = 0; i < objList.size(); i++ ) {
        csHeaderDef newHeader = objList.get(i);
        addAnnotationHeader( newHeader );
        if( myTempRemoveHeaders.contains( newHeader ) ) {
          myTempRemoveHeaders.remove( newHeader );
        }
        else {
          myTempAddHeaders.add( newHeader );
        }
      }
    }
    else {
//      JOptionPane.showMessageDialog(osAnnotationDialog.this,
//          "To add a trace header for annotation, first select\n" +
//          "the trace header from the list above.",
//          null, JOptionPane.INFORMATION_MESSAGE);
    }
  }
  public void addAnnotationHeader( String headerName ) {
    for( int ihdr = 0; ihdr < myTraceHeaders.size(); ihdr++ ) {
      csHeaderDef header = myTraceHeaders.get( ihdr );
      if( header.name.compareTo( headerName ) == 0 ) {
        if( !myListModelAnn.contains( header ) ) {
          addAnnotationHeader( header );
          apply();
        }
        return;
      }
    }
  }
  private void removeAnnotationHeaders() {
    if( !myListAnn.isSelectionEmpty() ) {
      int firstIndex = myListAnn.getSelectedIndex();
      java.util.List<csHeaderDef> objList = myListAnn.getSelectedValuesList();
      for( csHeaderDef objList1 : objList ) {
        csHeaderDef newHeader = objList1;
        removeAnnotationHeader( newHeader );
        if( myTempAddHeaders.contains( newHeader ) ) {
          myTempAddHeaders.remove( newHeader );
        }
        else {
          myTempRemoveHeaders.add( newHeader );
        }
      }
      if( myListModelAnn.size() != 0 ) {
        myListAnn.setSelectedIndex( Math.max( 0, firstIndex-1 ) );
      }

    }
    else {
//      JOptionPane.showMessageDialog(osAnnotationDialog.this,
//          "To remove an annotated trace header, first select\n" +
//          "the trace header from the list below.",
//          null, JOptionPane.INFORMATION_MESSAGE);
    }
  }
  private void addAnnotationHeader( csHeaderDef header ) {
    myListModelHeaders.removeElement(header);
    int index = myListModelAnn.getSize();
    myListModelAnn.add( index, header );
    myListAnn.setSelectedIndex( index );
    myTempAddHeaders.add( header );
  }
  private void removeAnnotationHeader( csHeaderDef header ) {
    myListModelAnn.removeElement(header);
    int size = myListModelHeaders.getSize();
    int index = size;
    for( int i = 0; i < size; i++ ) {
      if( header.name.compareTo( ((csHeaderDef)myListModelHeaders.getElementAt(i)).name ) <= 0  ) {
        index = i;
        break;
      }
    }
    myListModelHeaders.add( index, header );
    myListHeaders.setSelectedIndex( index );
  }
  public void changeSettings( csAnnotationAttributes attr ) {
    myBoxOmitRepeatedValues.setSelected( attr.omitRepeating );
    myBoxShowSeqTraceNum.setSelected( attr.showSequential );
    myBoxTraceLabelStep.setSelected( attr.fixedTraceLabelStep );
    if( attr.traceLabelStep > 0 ) myTextTraceLabelStep.setText("" + attr.traceLabelStep );
  }
  public void getSettings( csAnnotationAttributes attr ) {
    attr.omitRepeating = myBoxOmitRepeatedValues.isSelected();
    attr.showSequential = myBoxShowSeqTraceNum.isSelected();
    attr.fixedTraceLabelStep = myBoxTraceLabelStep.isSelected();
    try {
      attr.traceLabelStep = Integer.parseInt(myTextTraceLabelStep.getText());
    }
    catch( NumberFormatException e ) {
      // Nothing
    }
  }
  /**
   * Apply changes: Add/remove annotated trace headers
   */
  private void apply() {
    myTempAddHeaders.clear();
    myTempRemoveHeaders.clear();
    int numHeaders = myListModelAnn.size();
    csHeaderDef[] headers = new csHeaderDef[ numHeaders ];
    for( int i = 0; i < numHeaders; i++ ) {
      headers[i] = (csHeaderDef)myListModelAnn.get(i);
    }
    csAnnotationAttributes attr = new csAnnotationAttributes();
    getSettings( attr );
    mySeisPane.updateAnnotationTraceHeaders( headers, attr );
  }
  private void cancel() {
    for( int i = 0; i < myTempAddHeaders.size(); i++ ) {
      removeAnnotationHeader( myTempAddHeaders.get(i) );
    }
    for( int i = 0; i < myTempRemoveHeaders.size(); i++ ) {
      addAnnotationHeader( myTempRemoveHeaders.get(i) );
    }

    myTempAddHeaders.clear();
    myTempRemoveHeaders.clear();    
    dispose();
  }
  /**
   * Cell renderer for trace header annotation table
   * 
   * @author Bjorn Olofsson
   * @date   4 Feb 2008
   */
  public class AnnotationCellRenderer extends DefaultListCellRenderer {
    private static final int CELL_WIDTH = DEFAULT_WIDTH;
    private static final int NAME_WIDTH = (int)(CELL_WIDTH/3.5);
    private static final int DESC_WIDTH = CELL_WIDTH - NAME_WIDTH;
//    private boolean myIsInitialised = false;
//    private int myMaxWidth;
//    private int myWidth10Spaces;
    
    public Component getListCellRendererComponent(
      JList<?> list,
      Object value,            // value to display
      int index,               // cell index
      boolean isSelected,      // is the cell selected
      boolean cellHasFocus)    // the list and the cell have the focus
    {
//      Component comp = super.getListCellRendererComponent( list, value, index, isSelected, cellHasFocus);
      csHeaderDef header = (csHeaderDef)value;
/*
      JLabel label = (JLabel)comp;
      Graphics g = label.getGraphics();
      if( g == null ) {
        setText(" OOOOPS ");
        return this;
      }
      FontMetrics metrics = g.getFontMetrics( g.getFont() );
      if( !myIsInitialised ) {
        myMaxWidth = 0;
        for( int ih = 0; ih < myTraceHeaders.size(); ih++ ) {
          int labelWidth  = metrics.stringWidth( myTraceHeaders.get(ih).name );
          if( labelWidth > myMaxWidth ) myMaxWidth = labelWidth;
        }
        String tenSpaces = "          ";
        myWidth10Spaces = metrics.stringWidth( tenSpaces );
        myMaxWidth += myWidth10Spaces;
        myIsInitialised = true;
      }

      String text = header.name;
      int labelWidth  = metrics.stringWidth( text );
      int numAddSpaces = (int)(( 10.0*(double)(myMaxWidth-labelWidth) )/(double)myWidth10Spaces + 0.5);
      for( int i = 0; i < numAddSpaces; i++ ) {
        text += " ";
      }
      text += header.desc;
*/
//      setPreferredSize( new Dimension(CELL_WIDTH,getPreferredSize().height) );
      
      String text = "<html><table width=\"" + CELL_WIDTH + "\"><tr><td width=\"" +
                    NAME_WIDTH + "\" align=\"left\">" + header.name + "</td><td width=\"" +
                    DESC_WIDTH + "\" align=\"left\">" + header.desc + "</td></tr></table></html>";
      setText(text);
      
      return this;
    }
  }
}



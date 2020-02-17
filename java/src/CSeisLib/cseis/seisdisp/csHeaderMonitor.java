/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import javax.swing.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import cseis.seis.*;

/**
 * Trace header monitor.<br>
 * Displays all trace headers and header values of current seismic trace.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csHeaderMonitor extends JDialog {
  /// Header description
  private String[] myHeaderDesc; 
  private JTable myHeaderTable;
  private DefaultTableModel myTableModel;
  /// Conversion of model index to view index. Trace headers are unsorted in trace buffer, and sorted by name in table view
  private int[] myModel2ViewIndexArray;
  private int myNumHeaders;
  
  public csHeaderMonitor( JFrame frame, csHeaderDef[] headers ) {
    this( frame, headers, "Trace header monitor");
  }
  public csHeaderMonitor( JFrame frame, csHeaderDef[] headers, String title ) {
    super( frame, title );

    String[] columnNames = new String[2];
    columnNames[0] = "Header name";
    columnNames[1] = "Header value";
    myTableModel = new HeaderTableModel( columnNames, headers.length );
    myHeaderTable = new HeaderTable( myTableModel );
    updateHeaderNames( headers );
    pack();
    setSize(300,400);
  }
  /**
   * Update header names
   * @param headersIN Array of trace header objects (name, description..)
   */
  public void updateHeaderNames( csHeaderDef[] headersIN ) {
    myNumHeaders = headersIN.length;
    myHeaderDesc = new String[myNumHeaders];
    myTableModel.setRowCount( myNumHeaders );
    
    // Sort input headers (without changing input array), and set index array that converts from
    // 'model' row index to view row index.
    csHeaderDef[] headersTMP = new csHeaderDef[myNumHeaders];
    myModel2ViewIndexArray = new int[myNumHeaders];
    for( int ih = 0; ih < myNumHeaders; ih++ ) {
      headersTMP[ih] = new csHeaderDef( headersIN[ih] );
      headersTMP[ih].index = ih;
    }
    java.util.List<csHeaderDef> list = java.util.Arrays.asList(headersTMP);
    java.util.Collections.sort( list );
    for( int ih = 0; ih < myNumHeaders; ih++ ) {
      myModel2ViewIndexArray[list.get(ih).index] = ih;
    }
    for( int i = 0; i < myNumHeaders; i++ ) {
      myHeaderTable.setValueAt( headersIN[i], myModel2ViewIndexArray[i], 0 );
      myHeaderTable.setValueAt( new String(""), myModel2ViewIndexArray[i], 1 );
      myHeaderDesc[myModel2ViewIndexArray[i]] = headersIN[i].desc;
    }
    getContentPane().removeAll();    
    getContentPane().add( new JScrollPane(myHeaderTable), BorderLayout.CENTER );
    if( isVisible() ) setVisible(true);  // Ensure graphics are correctly updated
  }
  /**
   * Update header values
   * @param headerValues Array of header values
   */
  public void updateValues( csHeader[] headerValues ) {
    if( myNumHeaders != headerValues.length ) {
      return;
    }
    for( int i = 0; i < myNumHeaders; i++ ) {
      myHeaderTable.setValueAt( headerValues[i], myModel2ViewIndexArray[i], 1 );
    }
  }
  /**
   * Header table
   *
   */
  public class HeaderTable extends JTable {
    HeaderTable( TableModel model ) {
      super( model );
    }
    public String getToolTipText(MouseEvent e) {
      String tip = null;
      java.awt.Point p = e.getPoint();
      int rowIndex = rowAtPoint(p);
      int colIndex = columnAtPoint(p);
      int realColumnIndex = convertColumnIndexToModel(colIndex);

      if( realColumnIndex == 0 && (rowIndex >= 0 && rowIndex < myHeaderDesc.length) ) {
        tip = myHeaderDesc[rowIndex];  // + getValueAt(rowIndex, colIndex);
      }
//      else if( realColumnIndex == 1 && (rowIndex >= 0 && rowIndex < myHeaderDesc.length) ) {
//        tip = "Value = " + super.getValueAt( rowIndex, colIndex );  // + getValueAt(rowIndex, colIndex);
//      }
      else {
        tip = super.getToolTipText(e);
      }
      return tip;
    }
  }
  /**
   * Header table
   *
   */
  public class HeaderTableModel extends DefaultTableModel {
    HeaderTableModel( Object[] columnNames, int numColumns ) {
      super( columnNames, numColumns );
    }
    @Override
    public boolean isCellEditable(int row, int col) {
      return false;
    }
  }
}



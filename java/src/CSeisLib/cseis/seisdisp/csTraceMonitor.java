/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import javax.swing.*;
import javax.swing.table.*;

import java.awt.*;
import java.awt.event.*;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

import cseis.seis.*;
import cseis.swing.csFileChooser;

/**
 * Trace value monitor.<br>
 * Displays all sample values of current seismic trace.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csTraceMonitor extends JDialog implements csISampleInfoListener {
  public static final int NUM_COLUMNS = 3;
  private JTable myTraceTable;
  private DefaultTableModel myTableModel;
  private JLabel myLabel;
  private JTextField myTextTrace;
  private JCheckBox myBoxLockTrace;
  private JCheckBox myBoxLockScroll;
  private csISeismicTraceBuffer myTraceBuffer;
  private JButton myButtonSave;
  private int myCurrentTraceIndex;
  private JPanel myMainPanel;
  private DecimalFormat floatFmt1 = new DecimalFormat("0.0");
  private float mySampleInt;
  private int myFirstTraceIndex;
  private JScrollPane myScrollPane;
  
  public csTraceMonitor( JFrame frame, csISeismicTraceBuffer traceBuffer, float sampleInt ) {
    this( frame, traceBuffer, sampleInt, "Trace sample monitor" );
  }
  public csTraceMonitor( JFrame frame, csISeismicTraceBuffer traceBuffer, float sampleInt, String title ) {
    super( frame, title );
    myCurrentTraceIndex = 0;
    myFirstTraceIndex   = 0;
    
    myLabel = new JLabel(" Trace: ");
    myTextTrace = new JTextField("0");
    myBoxLockTrace  = new JCheckBox( "Lock trace", false );
    myBoxLockScroll = new JCheckBox( "Lock scroll", false );
    myButtonSave = new JButton("Save");
    myMainPanel = new JPanel( new GridBagLayout() );
    myTextTrace.setPreferredSize( new Dimension(10,myTextTrace.getPreferredSize().height) );


    int xp = 0;
    myMainPanel.add( myLabel, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myMainPanel.add( myTextTrace, new GridBagConstraints(
        xp++, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myMainPanel.add( myBoxLockTrace, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myMainPanel.add( myBoxLockScroll, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    myMainPanel.add( myButtonSave, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//    panel.add( Box.createHorizontalGlue(), new GridBagConstraints(
//        xp++, 0, 1, 1, 0.5, 1.0, GridBagConstraints.CENTER,
//        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    String[] columnNames = new String[NUM_COLUMNS];
    columnNames[0] = "Sample index";
    columnNames[1] = "Time [ms]";
    columnNames[2] = "Sample value";
    myTableModel = new DefaultTableModel( columnNames, traceBuffer.numSamples() );
    myTraceTable = new JTable( myTableModel );

    updateBuffer( traceBuffer, 0, sampleInt );

    pack();
    setSize(400,400);
//    setPreferredSize(new Dimension(300,400));
    
    myButtonSave.addActionListener( new ActionListener() {
      public void actionPerformed( ActionEvent e ) {
        csFileChooser fc = new csFileChooser( "Save sample values to ASCII file" );
        if( fc.selectFileForWriting(csTraceMonitor.this) ) {
          String filename = fc.getSelectedFile().getAbsolutePath();
          BufferedWriter writer = null;
          try {
            writer = new BufferedWriter( new FileWriter( filename ) );
          }
          catch( FileNotFoundException exc ) {
            return;
          }
          catch( IOException exc ) {
            return;
          }
          String text = "";
          float[] samples = myTraceBuffer.samples(myCurrentTraceIndex);
          for( int i = 0; i < myTraceBuffer.numSamples(); i++ ) {
            text += " " + (i+1+myFirstTraceIndex) + " " + floatFmt1.format(mySampleInt*i) + " " + samples[i] + "\n";
          }
          try {
            writer.write( text );
            writer.close();
          }
          catch( IOException exc ) {
            return;
          }
        }
      }
    });
    myTextTrace.addActionListener( new ActionListener() {
      public void actionPerformed( ActionEvent e ) {
        try {
          int traceNumber = Integer.parseInt( myTextTrace.getText() );  // trace number starting at 1 for first trace
          update( traceNumber-1-myFirstTraceIndex, -1 );  // -1 to get to trace index, starting at 0
        }
        catch( NumberFormatException exc ) {
        }
      }
    });
  }
  public void updateBuffer( csISeismicTraceBuffer traceBuffer ) {
    updateBuffer( traceBuffer, myFirstTraceIndex, mySampleInt );
  }
  public void updateBuffer( csISeismicTraceBuffer traceBuffer, int firstTraceIndex, float sampleInt ) {
    mySampleInt   = sampleInt;
    myTraceBuffer = traceBuffer;
    myFirstTraceIndex = firstTraceIndex;
    int numSamples = myTraceBuffer.numSamples();
    myTableModel.setRowCount( numSamples );

    myCurrentTraceIndex = 0;
    for( int i = 0; i < myTraceBuffer.numSamples(); i++ ) {
      myTraceTable.setValueAt( i+1, i, 0 );
      myTraceTable.setValueAt( floatFmt1.format(sampleInt*i), i, 1 );
      myTraceTable.setValueAt( new String(""), i, 2 );
    }

    myScrollPane = new JScrollPane(myTraceTable);
    getContentPane().removeAll();    
    getContentPane().add( myMainPanel, BorderLayout.NORTH );
    getContentPane().add( myScrollPane, BorderLayout.CENTER );
    if( isVisible() ) setVisible(true);  // Ensure graphics are correctly updated
  }
  public void mouseMoved( Object source, csSampleInfo info ) {
    if( myBoxLockTrace.isSelected() ) return;
    update( info.trace, info.sample );
  }
  public void mouseClicked( Object source, csSampleInfo info ) {
  }
  private void update( int traceIndex, int sampleIndex ) {
    if( traceIndex >= 0 && traceIndex < myTraceBuffer.numTraces() && isVisible() ) {
      float[] samples = myTraceBuffer.samples( traceIndex );
      for( int i = 0; i < myTraceBuffer.numSamples(); i++ ) {
        myTraceTable.setValueAt( samples[i], i, 2 );
      }
      myTextTrace.setText("" + (traceIndex+1+myFirstTraceIndex));
      myCurrentTraceIndex = traceIndex;
      int newScrollValue = 0;
      int newRow = 0;
      if( sampleIndex >= 0 && sampleIndex < myTraceTable.getRowCount() ) {
        int max = myScrollPane.getVerticalScrollBar().getMaximum();
        newScrollValue = (int)( max * ( (float)sampleIndex / (float)myTraceTable.getRowCount() ) );
        newRow = sampleIndex;
      }
      myTraceTable.setRowSelectionInterval( newRow, newRow );
      if( !myBoxLockScroll.isSelected() ) {
        myScrollPane.getVerticalScrollBar().setValue( newScrollValue );
      }
    }
  }
}



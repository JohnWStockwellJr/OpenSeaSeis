/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.segy;

import java.awt.*;
import javax.swing.*;
import javax.swing.table.*;
import cseis.jni.csSegyTrcHeaderDefinition;
import java.util.ArrayList;

/**
 * View of SEGY header contents.
 * @author 2009 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csSegyHeaderView extends JDialog {
  public csSegyHeaderView( JFrame parentFrame, String title, ArrayList<csSegyTrcHeaderDefinition> hdrDef ) {
    this( parentFrame, title, hdrDef, 0 );
  }
  public csSegyHeaderView( JFrame parentFrame, String title, ArrayList<csSegyTrcHeaderDefinition> hdrDef, int tabIndex ) {
    super( parentFrame, title );

    JTabbedPane tabPane = new JTabbedPane();
    tabPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
    int numTabs = hdrDef.size();
    ArrayList<String> columnNames = new ArrayList<String>();
    columnNames.add( "Name" );
    columnNames.add( "Description" );
    columnNames.add( "Bytes" );
    
    for( int itab = 0; itab < numTabs; itab++ ) {
      JTable table = new JTable( new DefaultTableModel( columnNames.toArray(), hdrDef.get(itab).numHeaders() ) {
        public boolean isCellEditable( int row, int column ) { return false; }
      });
      table.getColumnModel().getColumn(0).setPreferredWidth( 60 );
      table.getColumnModel().getColumn(1).setPreferredWidth( 140 );
      table.getColumnModel().getColumn(2).setPreferredWidth( 12 );
      int numHeaders = hdrDef.get(itab).numHeaders();
      for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
        int byteLoc1 = hdrDef.get(itab).byteLoc(ihdr)+1;
        int byteLoc2 = byteLoc1 + hdrDef.get(itab).byteSize(ihdr) - 1;
        table.getModel().setValueAt( hdrDef.get(itab).name(ihdr), ihdr, 0 );
        table.getModel().setValueAt( hdrDef.get(itab).desc(ihdr), ihdr, 1 );
        table.getModel().setValueAt( byteLoc1+"-"+byteLoc2, ihdr, 2 );
      }
      table.doLayout();
      tabPane.addTab( hdrDef.get(itab).title(), new JScrollPane(table) );
    }

    tabPane.setSelectedIndex(tabIndex);
    JPanel panelText = new JPanel( new BorderLayout() );
    panelText.add( tabPane );
    getContentPane().add( panelText ,BorderLayout.CENTER);
    pack();
    setLocationRelativeTo(parentFrame);
  }
  public csSegyHeaderView( JFrame parentFrame, String title, String[] tabTitle, String[] textToShow ) {
    super( parentFrame, title );

    JTabbedPane tabPane = new JTabbedPane();
    
    int numTabs = textToShow.length;
    for( int i = 0; i < numTabs; i++ ) {
      
      JTextArea textArea = new JTextArea();
      textArea.setWrapStyleWord( false );
      textArea.setLineWrap( false );
      textArea.setEditable( false );
      textArea.setBackground( Color.white );
      textArea.setText( textToShow[i] );
      tabPane.addTab( tabTitle[i], textArea );
    }
    
    JPanel panelText = new JPanel( new BorderLayout() );
    panelText.add( tabPane );
    JScrollPane pane = new JScrollPane( panelText, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );

    getContentPane().add(pane,BorderLayout.CENTER);
    pack();
    setLocationRelativeTo(parentFrame);
  }
  public csSegyHeaderView( JFrame parentFrame, String title, String textToShow ) {
    super( parentFrame, title );

    JTextArea textArea = new JTextArea();
    textArea.setWrapStyleWord( false );
    textArea.setLineWrap( false );
    textArea.setEditable( false );
    textArea.setText( textToShow );
    textArea.setBackground( Color.white );

    JPanel panelText = new JPanel( new BorderLayout() );
    panelText.add( textArea );
    JScrollPane pane = new JScrollPane( panelText, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );

    getContentPane().add(pane,BorderLayout.CENTER);
    pack();
    setLocationRelativeTo(parentFrame);
  }
}



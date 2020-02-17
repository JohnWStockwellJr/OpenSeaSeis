/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;

import cseis.general.csColorMap;
import cseis.general.csCustomColorMap;
import cseis.seisdisp.csColorMapListItem;
import java.awt.event.ActionEvent;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.event.ListDataEvent;
import javax.swing.event.ListDataListener;

/**
 * Test color map generation.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class TestColorMapGenerator extends JFrame implements csIColorMapGeneratorListener, ListDataListener {
  private csColorMapComboBox myCombo;
  private csManagerCustomMapList myList;
  private csCustomColorMapModel myCustomColorMapModel;
  private csColorMapListModel myColorMapModel;
  private JMenuBar myMenuBar;
  private JMenu myMenuFile;
  private JMenuItem myMenuOpen;
  
  public TestColorMapGenerator() {

    myCustomColorMapModel = new csCustomColorMapModel( csColorMap.COLOR_MAP_TYPE_32BIT );
    myColorMapModel = new csColorMapListModel();
    myColorMapModel.addItem( csColorMapListItem.createStandardItem(new csCustomColorMap(Color.yellow)) );
    myCombo = new csColorMapComboBox( myCustomColorMapModel );
    myList = new csManagerCustomMapList( myCustomColorMapModel );

    myMenuBar = new JMenuBar();
    myMenuFile = new JMenu("File");
    myMenuOpen = new JMenuItem("Open...");

    myMenuFile.add( myMenuOpen );
    myMenuBar.add( myMenuFile );
    setJMenuBar( myMenuBar );
    setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );

    getContentPane().add( myCombo, BorderLayout.WEST );
    getContentPane().add( myList, BorderLayout.EAST );
    pack();
    setLocation( 500, 400 );
    
    myMenuOpen.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        openColorMapGenerator();
      }
    });
  }
  public void openColorMapGenerator() {
    final csColorMapGenerator generator = new csColorMapGenerator(
            this, new java.util.ArrayList<csCustomColorMap>(), null );
    generator.addListener( myCustomColorMapModel );
/*    generator.addListDataListener( new ListDataListener() {
      @Override
      public void intervalAdded(ListDataEvent e) {
        e.getIndex0();
        e.getIndex1();
      }
      @Override
      public void intervalRemoved(ListDataEvent e) {
      }
      @Override
      public void contentsChanged(ListDataEvent e) {
      }
    });
 * 
 */
    generator.setVisible(true);
  }
  
  public static void main( String[] args ) {
    TestColorMapGenerator t = new TestColorMapGenerator();
    t.setVisible(true);
  }

  @Override
  public void applyColorMap( csCustomColorMap cmap ) {
  }

  @Override
  public void intervalAdded(ListDataEvent e) {
    e.getIndex0();
    e.getIndex1();
  }
  @Override
  public void intervalRemoved(ListDataEvent e) {
  }
  @Override
  public void contentsChanged(ListDataEvent e) {
  }

  @Override
  public void updateColorMaps( java.util.List<csCustomColorMap> list ) {
    ArrayList<csColorMapListItem> itemList = new ArrayList<csColorMapListItem>(list.size());
    for( int i = 0; i < list.size(); i++ ) {
      itemList.add( csColorMapListItem.createStandardItem( new csCustomColorMap(list.get(i)) ) );
    }
    myColorMapModel.update( itemList );
  }

}



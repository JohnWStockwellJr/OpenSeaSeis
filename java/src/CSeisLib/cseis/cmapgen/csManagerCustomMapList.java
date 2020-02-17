/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;

import cseis.general.csColorMap;
import cseis.general.csCustomColorMap;
import cseis.seisdisp.csColorMapListItem;
import cseis.seisdisp.csComboColorMapRenderer;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListDataListener;

/**
 * Manages custom color maps for color map generator application.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csManagerCustomMapList extends JList<csColorMapListItem> {
  private csColorMapListModel myColManagerListModel;
  
  public csManagerCustomMapList() {
    this( new csColorMapListModel() );
  }
  public csManagerCustomMapList( csColorMapListModel model ) {
    myColManagerListModel = model;
    setModel(myColManagerListModel);
    setCellRenderer( new csComboColorMapRenderer() );
    setSelectionMode( ListSelectionModel.SINGLE_SELECTION );
  }
  public void addListDataListener( ListDataListener listener ) {
    myColManagerListModel.addListDataListener( listener );    
  }
  public void removeListDataListener( ListDataListener listener ) {
    myColManagerListModel.removeListDataListener( listener );    
  }
  public int getNumMaps() {
    return myColManagerListModel.getSize();
  }
  public csCustomColorMap getMap( int index ) {
    return (csCustomColorMap)myColManagerListModel.getElementAt(index).map;
  }
  public csCustomColorMap getSelectedMap() {
    if( isSelectionEmpty() ) return null;
    return (csCustomColorMap)((csColorMapListItem)getSelectedValue()).map;
  }
  public int getSelectedMapIndex() {
    if( isSelectionEmpty() ) return -1;
    return getSelectedIndex();
  }
  public void removeSelectedMap() {
    if( isSelectionEmpty() ) return;
    int itemIndex = super.getSelectedIndex();
    if( itemIndex >= 0 && itemIndex < myColManagerListModel.getSize() ) {
      myColManagerListModel.removeItem( itemIndex );
    }
  }
  public void addMap( csCustomColorMap cmap ) {
    csCustomColorMap newColorMap = new csCustomColorMap(cmap);
    int index = getColorMapIndex( cmap );
    if( index >= 0 ) {
      int option = JOptionPane.showConfirmDialog( this,
                "Another color map with the name\n'"+ cmap.toString() +"'\nalready exists.\nOK to overwrite?",
                "Create new color map", JOptionPane.YES_NO_OPTION );
      if( option != JOptionPane.YES_OPTION ) return;
      myColManagerListModel.setItem( index, csColorMapListItem.createStandardItem( newColorMap ) );
      return;
    }
    else if( defaultMapExists(cmap.toString()) ) {
      JOptionPane.showMessageDialog( this,
              "A default color map with the name\n'"+ cmap.toString() +"'\nalready exists.\nSpecify a different name.",
              "Create new color map", JOptionPane.INFORMATION_MESSAGE );
      return;
    }
    myColManagerListModel.addItem( csColorMapListItem.createStandardItem( newColorMap ) );
  }
  private int getColorMapIndex( csCustomColorMap cmap ) {
    String name = cmap.toString();
    int numItems = getModel().getSize();
    for( int i = 0; i < numItems; i++ ) {
      csColorMapListItem item = (csColorMapListItem)getModel().getElementAt(i);
      if( name.compareTo( item.map.toString()) == 0 ) {
        return i;
      }
    }
    return -1;
  }
  public boolean updateColorMap( int cmapIndex, csCustomColorMap cmap ) {
    if( cmapIndex < 0 || cmapIndex >= myColManagerListModel.getSize() ) return false;
    int index = getColorMapIndex( cmap );
    if( index >= 0 && index != cmapIndex ) {
      JOptionPane.showMessageDialog( this,
              "Another color map with the name\n'"+ cmap.toString() +"'\nalready exists.\nSpecify a different name.",
              "Update color map", JOptionPane.INFORMATION_MESSAGE );
      return false;
    }
    else if( defaultMapExists(cmap.toString()) ) {
      JOptionPane.showMessageDialog( this,
              "A default color map with the name\n'"+ cmap.toString() +"'\nalready exists.\nSpecify a different name.",
              "Update color map", JOptionPane.INFORMATION_MESSAGE );
      return false;
    }
    csCustomColorMap newColorMap = new csCustomColorMap(cmap);
    myColManagerListModel.setItem( cmapIndex, csColorMapListItem.createStandardItem( newColorMap ) );
    return true;
  }
//  public void updateColorMaps( java.util.List<csCustomColorMap> list ) {
//    ArrayList<csColorMapListItem> itemList = new ArrayList<csColorMapListItem>(list.size());
//    for( int i = 0; i < list.size(); i++ ) {
//      itemList.add( csColorMapListItem.createStandardItem( new csCustomColorMap(list.get(i)) ) );
//    }
//    myColManagerListModel.update( itemList );
//  }
  private boolean defaultMapExists( String name ) {
    for( int i = 0; i < csColorMap.NUM_DEFAULT_MAPS; i++ ) {
      if( csColorMap.TEXT[i].compareTo(name) == 0 ) {
        return true;
      }
    }
    return false;
  }
}


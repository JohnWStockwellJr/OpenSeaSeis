/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;
  
import cseis.general.csColorMap;
import cseis.general.csCustomColorMap;
import cseis.seisdisp.csColorMapListItem;
import java.util.ArrayList;
import javax.swing.*;

/**
 * Implementation of AbstractListModel to enable listing color maps in Java GUI objects.<br>
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csColorMapListModel extends AbstractListModel<csColorMapListItem> {
  private ArrayList<csColorMapListItem> myItemList;

  public csColorMapListModel() {
    this( new ArrayList<csCustomColorMap>() );
  }
  public csColorMapListModel( java.util.List<csCustomColorMap> list ) {
    myItemList = new ArrayList<csColorMapListItem>( list.size() );
    for( csCustomColorMap list1 : list ) {
      myItemList.add(csColorMapListItem.createStandardItem(list1));
    }
  }
  public void updateColorMapType( int colorMapType ) {
    for( int i = 0; i < csColorMap.NUM_DEFAULT_MAPS; i++ ) {
      myItemList.get( i ).map.setColorMapType( colorMapType );
    }
  }
  public int indexOf( Object item ) {
    if( item instanceof csColorMapListItem ) {
      return myItemList.indexOf( (csColorMapListItem)item );
    }
    else {
      return -1;
    }
  }
  protected boolean contains( Object item ) {
    if( item instanceof csColorMapListItem ) {
      return myItemList.contains( (csColorMapListItem)item );
    }
    else {
      return false;
    }
  }
  @Override
  public int getSize() {
    return myItemList.size();
  }
  @Override
  public csColorMapListItem getElementAt( int index ) {
    return( (index >= 0 && index < myItemList.size()) ? myItemList.get(index) : null );
  }
  public boolean isEmpty() {
    return( myItemList.isEmpty() );
  }
//----------------------------------------------
  public void addItem( csColorMapListItem item ) {
    myItemList.add( item );
    int newIndex = myItemList.size()-1;
    fireIntervalAdded( this, newIndex, newIndex );
  }
  public void setItem( int atIndex, csColorMapListItem item ) {
    myItemList.set( atIndex, item );
    fireContentsChanged( this, atIndex, atIndex );
  }
  public void removeItem( csColorMapListItem item ) {
    int numItems = myItemList.size();
    for( int i = 0; i < numItems; i++ ) {
      if( myItemList.get(i).map.equals( item.map) ) {
        myItemList.remove(i);
        fireIntervalRemoved( this, i, i );
        return;
      }
    }
  }
  public void removeItem( int itemIndex ) {
    if( itemIndex < 0 || itemIndex >= myItemList.size() ) return;
    myItemList.remove( itemIndex );
    fireIntervalRemoved( this, itemIndex, itemIndex );
  }
  public void update( java.util.List<csColorMapListItem> list ) {
    myItemList.clear();
    for( int i = 0; i < list.size(); i++ ) {
      myItemList.add( list.get(i) );
    }
    fireContentsChanged( this, 0, getSize()-1 );
  }
}


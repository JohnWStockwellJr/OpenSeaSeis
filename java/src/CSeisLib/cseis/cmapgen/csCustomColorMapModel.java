/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;

import cseis.general.csColorMap;
import cseis.general.csCustomColorMap;
import cseis.seisdisp.csColorMapListItem;
import java.util.ArrayList;

/**
 * Extends csColorMapListModel for custom color maps.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csCustomColorMapModel extends csColorMapListModel implements csIColorMapGeneratorListener {
  private ArrayList<csColorMapListItem> myDefaultMapList;
  
  public csCustomColorMapModel( int colorMapType ) {
    myDefaultMapList = new ArrayList<csColorMapListItem>( csColorMap.NUM_DEFAULT_MAPS );
    for( int i = 0; i < csColorMap.NUM_DEFAULT_MAPS; i++ ) {
      csColorMap cmap = new csColorMap( i, colorMapType );
      myDefaultMapList.add( csColorMapListItem.createStandardItem( cmap ) );
    }
    update( myDefaultMapList );
  }
  @Override
  public void updateColorMapType( int colorMapType ) {
    super.updateColorMapType(colorMapType);
    for( int i = 0; i < csColorMap.NUM_DEFAULT_MAPS; i++ ) {
      myDefaultMapList.get( i ).map.setColorMapType( colorMapType );
    }
  }
  @Override
  public void applyColorMap( csCustomColorMap cmap ) {
    // Nothing to be done here
  }

  @Override
  public void updateColorMaps( java.util.List<csCustomColorMap> list ) {
    if( list == null ) return;
    ArrayList<csColorMapListItem> itemList = new ArrayList<csColorMapListItem>(list.size());
    for( int i = 0; i < myDefaultMapList.size(); i++ ) {
      itemList.add( myDefaultMapList.get(i) );
    }
    for( int i = 0; i < list.size(); i++ ) {
      itemList.add( csColorMapListItem.createStandardItem( new csCustomColorMap(list.get(i)) ) );
    }
    update( itemList );
  }
}


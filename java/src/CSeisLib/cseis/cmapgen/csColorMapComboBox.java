/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;

import cseis.seisdisp.csComboColorMapRenderer;
import javax.swing.event.*;
import javax.swing.*;

/**
 * JComboBox used to list color maps.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csColorMapComboBox extends JComboBox<Object> {

  public csColorMapComboBox( csColorMapListModel listModel ) {
    setModel( new csColorMapComboModel(listModel) );
    SwingUtilities.invokeLater( new Runnable() {
      @Override
      public void run() {
        setRenderer( new csComboColorMapRenderer() );
      }
    });
  }
//----------------------------------------------------------
  public class csColorMapComboModel implements ComboBoxModel<Object>, ListDataListener {
    private Object mySelectedItem = null;
    private int mySelectedIndex   = -1;
    private csColorMapListModel myListModel;

    public csColorMapComboModel( csColorMapListModel listModel ) {
      myListModel = listModel;
      listModel.addListDataListener(this);
    }
    @Override
    public void intervalAdded( ListDataEvent event ) {
      int nItems = event.getIndex1() - event.getIndex0() + 1;

      if( mySelectedIndex >= event.getIndex0() ) {
        mySelectedIndex += nItems;
        setSelectedIndex( mySelectedIndex );
      }
    }
    @Override
    public void intervalRemoved( ListDataEvent event ) {
      int nItems = event.getIndex1() - event.getIndex0() + 1;

      if( mySelectedIndex > event.getIndex1() ) {
        mySelectedIndex -= nItems;
        setSelectedIndex( mySelectedIndex );
      }
      else if( mySelectedIndex >= event.getIndex0() ) {
        mySelectedIndex = -1;
        setSelectedIndex( mySelectedIndex );
      }
    }
    @Override
    public void contentsChanged( ListDataEvent event ) {
      if( mySelectedIndex >= myListModel.getSize() ) {
        mySelectedIndex = -1;
      }
      setSelectedIndex( mySelectedIndex );
    }
    public int getSelectedIndex() {
      return mySelectedIndex;
    }
    public void setSelectedIndex( int index ) {
      mySelectedIndex = index;
      mySelectedItem  = myListModel.getElementAt( mySelectedIndex );
    }
    @Override
    public Object getSelectedItem() {
      return mySelectedItem;
    }
    @Override
    public void setSelectedItem( Object item ) {
      if( item == null || !myListModel.contains(item) ) {
        mySelectedItem = null;
        mySelectedIndex = -1;
      }
      else {
        mySelectedItem = item;
        mySelectedIndex = myListModel.indexOf(mySelectedItem);
      }
    }
    //------------------------------------------------------------
    @Override
    public void addListDataListener( ListDataListener listener ) {
      myListModel.addListDataListener(listener);
    }
    @Override
    public Object getElementAt(int index) {
      return myListModel.getElementAt(index);
    }
    @Override
    public int getSize() {
      return myListModel.getSize();
    }
    @Override
    public void removeListDataListener(ListDataListener l) {
      myListModel.removeListDataListener(l);
    }

  }
}


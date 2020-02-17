/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.Component;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JLabel;
import javax.swing.JList;

/**
 * Renderer for JComboBox containing color maps.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csComboColorMapRenderer extends DefaultListCellRenderer {
  public csComboColorMapRenderer() {
//      setOpaque(true);
  }
  @Override
  public Component getListCellRendererComponent(
      JList<?> list,
      Object value,
      int index,
      boolean isSelected,
      boolean cellHasFocus)
  {
    Component comp   = super.getListCellRendererComponent( list, value, index, isSelected, cellHasFocus);
    JLabel label     = (JLabel)comp;
    csColorMapListItem item   = (csColorMapListItem)value;
    if( value == null ) return this;
//    label.setPreferredSize( comp.getPreferredSize() );
    label.setText( item.toString() );
    label.setIcon( item.icon );
    label.setToolTipText( item.toString() );
    return this;
  }
}


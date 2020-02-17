/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;
import javax.swing.JPopupMenu;

/**
 * Popup menu used by csGraph2D objects.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csGraphPopupMenu extends JPopupMenu {
  protected csShowGraphSettingsAction myShowGraphSettingsAction;
  protected csGraph2D myGraph;
  
  csGraphPopupMenu( csGraph2D graph ) {
    super();
    myGraph = graph;
    myShowGraphSettingsAction = new csShowGraphSettingsAction();
    
    add( myShowGraphSettingsAction );
  }
  public csShowGraphSettingsAction getGraphSettingsAction() {
    return myShowGraphSettingsAction;
  }
  class csShowGraphSettingsAction extends AbstractAction {
    public csShowGraphSettingsAction() {
      super( "Graph settings" );
      putValue( AbstractAction.SHORT_DESCRIPTION, "Graph settings" );
    }
    public void actionPerformed( ActionEvent e ) {
      myGraph.showDialog();
    }
  }
}


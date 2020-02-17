/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.awt.event.*;
import javax.swing.*;

@SuppressWarnings("serial")
public class csSemblanceView extends csSeisView {

  public csSemblanceView( JFrame parentFrame ) {
    super( parentFrame );
    myPopupMenu = new csSemblanceSeisViewPopupMenu( this );

//    mySeisViewBackground.addSeisViewListener( new csISeisViewListener() {
//      @Override
//      public void changedSettings(csSeisDispSettings settings) {
//        repaint();
//      }
//      @Override
//      public void vertScrollChanged(int scrollValue) {
//      }
//      @Override
//      public void horzScrollChanged(int scrollValue) {
//      }
//      @Override
//      public void sizeChanged(Dimension size) {
//      }
//      @Override
//      public void traceBufferChanged(csISeismicTraceBuffer traceBuffer) {
//      }
//    });
//  }
  }
  //==============================================================================
  //
  public class csSemblanceSeisViewPopupMenu extends csSeisViewPopupMenu {
    csSemblanceSeisViewPopupMenu( csSemblanceView seisView ) {
      super( seisView );
      JMenuItem itemBackgroundSettings = new JMenuItem("Display settings (background)...");
      addSeparator();
      add( myShowDispSettingsAction );
      add( itemBackgroundSettings );
      itemBackgroundSettings.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed( ActionEvent e ) {
//          csSemblanceView.this.showBackgroundSettingsDialog();
        }
      });
    }
  }

}

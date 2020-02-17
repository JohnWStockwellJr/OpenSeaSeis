/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.seisdisp.csMouseModes;
import cseis.seisdisp.csSeisPane;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;

/**
 * Toolbar located on left side of SeaView application.
 * @author 2011 Bjorn Olofsson
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csSeaViewToolBarLeft extends JToolBar {
  private final JButton myButtonZoomInVert;
  private final JButton myButtonZoomInHorz;
  private final JButton myButtonZoomOutVert;
  private final JButton myButtonZoomOutHorz;
  private final JButton myButtonZoomInBoth;
  private final JButton myButtonZoomOutBoth;
  private final JButton myButtonFitToScreen;

  private JToggleButton myButtonKillTraceMode;
  private JToggleButton myButtonSpectrum;
  private JToggleButton myButtonRubberZoom;
  private JToggleButton myButtonPanMode;
  private SeaView mySeaView;

  public csSeaViewToolBarLeft( SeaView seaview ) {
    super(JToolBar.VERTICAL);
    mySeaView = seaview;

    myButtonZoomOutBoth = new JButton( csSeaViewActions.getIcon(csSeaViewActions.ZoomOutAction) );
    myButtonZoomOutHorz = new JButton( csSeaViewActions.getIcon(csSeaViewActions.ZoomOutHorzAction) );
    myButtonZoomOutVert = new JButton( csSeaViewActions.getIcon(csSeaViewActions.ZoomOutVertAction) );
    myButtonZoomInBoth  = new JButton( csSeaViewActions.getIcon(csSeaViewActions.ZoomInAction) );
    myButtonZoomInHorz  = new JButton( csSeaViewActions.getIcon(csSeaViewActions.ZoomInHorzAction) );
    myButtonZoomInVert  = new JButton( csSeaViewActions.getIcon(csSeaViewActions.ZoomInVertAction) );
    myButtonSpectrum    = new JToggleButton( csSeaViewActions.getIcon(csSeaViewActions.SpectrumAction) );
    myButtonRubberZoom  = new JToggleButton( csSeaViewActions.getIcon(csSeaViewActions.RubberBandZoomAction) );
    myButtonPanMode     = new JToggleButton( csSeaViewActions.getIcon(csSeaViewActions.PanModeAction) );
    myButtonKillTraceMode   = new JToggleButton( csSeaViewActions.getIcon(csSeaViewActions.KillTraceAction) );
    myButtonFitToScreen     = new JButton( csSeaViewActions.getIcon(csSeaViewActions.FitToScreenAction) );

    myButtonZoomOutBoth.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ZoomOutAction] );
    myButtonZoomOutHorz.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ZoomOutHorzAction] );
    myButtonZoomOutVert.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ZoomOutVertAction] );
    myButtonZoomInBoth .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ZoomInAction] );
    myButtonZoomInHorz .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ZoomInHorzAction] );
    myButtonZoomInVert .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ZoomInVertAction] );
    myButtonSpectrum   .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.SpectrumAction] );
    myButtonRubberZoom .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.RubberBandZoomAction] );
    myButtonPanMode    .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.PanModeAction] );
    myButtonKillTraceMode.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.KillTraceAction] );
    myButtonFitToScreen    .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.FitToScreenAction] );
    
    add(myButtonRubberZoom);
    add(myButtonPanMode);
    add(myButtonSpectrum);
    add(myButtonKillTraceMode);
    addSeparator();
    add(myButtonZoomInBoth);
    add(myButtonZoomOutBoth);
    addSeparator();
    add(myButtonZoomInHorz);
    add(myButtonZoomOutHorz);
    addSeparator();
    add(myButtonZoomInVert);
    add(myButtonZoomOutVert);
    addSeparator();
    add(myButtonFitToScreen);

    myButtonRubberZoom.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.setMouseMode( csMouseModes.ZOOM_MODE, myButtonRubberZoom );
      }
    });
    myButtonPanMode.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.setMouseMode( csMouseModes.PAN_MODE, myButtonPanMode );
      }
    });
    myButtonSpectrum.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.setMouseMode( csMouseModes.SPECTRUM_MODE, myButtonSpectrum );
      }
    });
    myButtonKillTraceMode.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.setMouseMode( csMouseModes.KILL_MODE, myButtonKillTraceMode );
      }
    });
    myButtonZoomInBoth.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.zoom( csSeisPane.ZOOM_IN, csSeisPane.ZOOM_BOTH );
      }
    });
    myButtonZoomOutBoth.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.zoom( csSeisPane.ZOOM_OUT, csSeisPane.ZOOM_BOTH );
      }
    });
    myButtonZoomInHorz.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.zoom( csSeisPane.ZOOM_IN, csSeisPane.ZOOM_HORZ );

      }
    });
    myButtonZoomOutHorz.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.zoom( csSeisPane.ZOOM_OUT, csSeisPane.ZOOM_HORZ );

      }
    });
    myButtonZoomInVert.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.zoom( csSeisPane.ZOOM_IN, csSeisPane.ZOOM_VERT );
      }
    });
    myButtonZoomOutVert.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.zoom( csSeisPane.ZOOM_OUT, csSeisPane.ZOOM_VERT );
      }
    });
    myButtonFitToScreen.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.fitToScreen();
      }
    });
  }
  public int getMouseMode() {
    if( myButtonRubberZoom.isSelected() ) {
      return csMouseModes.ZOOM_MODE;
    }
    else if( myButtonPanMode.isSelected() ) {
      return csMouseModes.PAN_MODE;
    }
    else if( myButtonSpectrum.isSelected() ) {
      return csMouseModes.SPECTRUM_MODE;
    }
    else if( myButtonKillTraceMode.isSelected() ) {
      return csMouseModes.KILL_MODE;
    }
    return csMouseModes.NO_MODE;
  }
  public void resetToggleButtons() {
    myButtonRubberZoom.setSelected(false);
    myButtonSpectrum.setSelected(false);
    myButtonKillTraceMode.setSelected(false);
    myButtonPanMode.setSelected(false);
  }
  public void setToolbarEnabled( boolean set ) {
    myButtonRubberZoom.setEnabled(set);
    myButtonZoomOutBoth.setEnabled(set);
    myButtonZoomOutHorz.setEnabled(set);
    myButtonZoomOutVert.setEnabled(set);
    myButtonZoomInBoth.setEnabled(set);
    myButtonZoomInHorz.setEnabled(set);
    myButtonZoomInVert.setEnabled(set);
    myButtonFitToScreen.setEnabled(set);

    myButtonKillTraceMode.setEnabled(set);
    myButtonSpectrum.setEnabled(set);
    myButtonPanMode.setEnabled(set);
  }
}


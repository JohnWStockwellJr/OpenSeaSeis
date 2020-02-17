/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.csITimerTaskDecrementListener;
import cseis.general.csITimerTaskIncrementListener;
import cseis.general.csTimerTaskDecrement;
import cseis.general.csTimerTaskIncrement;
import cseis.swing.csDockPaneActions;
import cseis.swing.csDockPaneManager;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Timer;
import javax.swing.Box;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;

/**
 * Toolbar located on left side of SeaView application.
 * @author 2011 Bjorn Olofsson
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csSeaViewToolBarTop extends JToolBar {
  private JButton myButtonForwardSeismic;
  private JButton myButtonBackwardSeismic;
  private JButton myButtonBeginSeismic;
  private JButton myButtonEndSeismic;
  private JButton myButtonSnapShot;
  private JButton myButtonSnapShotPane;
  private JButton myButtonIncScaling;
  private JButton myButtonDecScaling;
  private float myCurrentScalar;
  private JTextField myTextScalar;
  private JCheckBox myBoxOpenInNewPane;
  private JCheckBox myBoxLockScalar;
  private JComboBox<Integer> myComboLayout;
  private JButton myButtonSelectPanes;
  private JToggleButton myButtonSyncSettings;
  private ArrayList<csISeaViewToolBarListener> myListeners;

  private SeaView mySeaView;
  
  private csTimerTaskDecrement myDecTimerTask = null;
  private csTimerTaskIncrement myIncTimerTask = null;
  
  public csSeaViewToolBarTop( SeaView seaview, csSeaViewProperties properties ) {
    super(JToolBar.HORIZONTAL);

    mySeaView = seaview;
    
    myListeners = new ArrayList<csISeaViewToolBarListener>();
    myBoxOpenInNewPane = new JCheckBox("New", true);
    myBoxOpenInNewPane.setToolTipText("Open data sets in new pane. De-select to open in active pane");
    myBoxLockScalar = new JCheckBox("Lock", false);
    myBoxLockScalar.setToolTipText("Lock scalar when opening new or when moving through data set");
    myCurrentScalar = 1.0f;
    myTextScalar = new JTextField("" + myCurrentScalar);
    int preferredHeight = myTextScalar.getPreferredSize().height;
    myTextScalar.setPreferredSize( new Dimension(150,preferredHeight) );
    myTextScalar.setMaximumSize( new Dimension(150,preferredHeight) );
    
    myButtonIncScaling = new JButton( csSeaViewActions.getIcon(csSeaViewActions.IncreaseScalingAction) );
    myButtonDecScaling = new JButton( csSeaViewActions.getIcon(csSeaViewActions.DecreaseScalingAction) );
    myButtonIncScaling.setMargin(new Insets(0,0,1,1));
    myButtonDecScaling.setMargin(new Insets(0,0,1,1));
    myButtonIncScaling.setPreferredSize(new Dimension(preferredHeight,preferredHeight));
    myButtonDecScaling.setPreferredSize(new Dimension(preferredHeight,preferredHeight));

    myButtonForwardSeismic  = new JButton( csSeaViewActions.getIcon(csSeaViewActions.ForwardSeismicAction) );
    myButtonBackwardSeismic = new JButton( csSeaViewActions.getIcon(csSeaViewActions.BackwardSeismicAction) );
    myButtonBeginSeismic    = new JButton( csSeaViewActions.getIcon(csSeaViewActions.BeginSeismicAction) );
    myButtonEndSeismic      = new JButton( csSeaViewActions.getIcon(csSeaViewActions.EndSeismicAction) );
    myButtonSnapShot        = new JButton( csSeaViewActions.getIcon(csSeaViewActions.SnapShotAction) );
    myButtonSnapShotPane    = new JButton( csSeaViewActions.getIcon(csSeaViewActions.SnapShotPaneAction) );
    myButtonSelectPanes     = new JButton( csSeaViewActions.getIcon(csSeaViewActions.SelectPanesAction) );
    myButtonSyncSettings    = new JToggleButton( csDockPaneActions.getIcon(csDockPaneActions.SyncAction) );
    myButtonSyncSettings.setSelected(false);
    
    myButtonIncScaling.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.IncreaseScalingAction] );
    myButtonDecScaling.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.DecreaseScalingAction] );
    myButtonForwardSeismic .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.ForwardSeismicAction] );
    myButtonBackwardSeismic.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.BackwardSeismicAction] );
    myButtonBeginSeismic   .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.BeginSeismicAction] );
    myButtonEndSeismic     .setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.EndSeismicAction] );
    myButtonSnapShot.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.SnapShotAction] );
    myButtonSnapShotPane.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.SnapShotPaneAction] );
    myButtonSelectPanes.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.SelectPanesAction] );
    myTextScalar.setToolTipText("Display scalar");
    myButtonSyncSettings.setToolTipText("<html>Toggle to synchronize ALL settings across sync'ed displays,<br>including display settings and zom level</html>");
    
    Integer[] comboItems = new Integer[csDockPaneManager.getNumLayoutOptions()];
    for( int i = 0; i < csDockPaneManager.getNumLayoutOptions(); i++ ) {
      comboItems[i] = i;
    }
    myComboLayout = new JComboBox<Integer>( comboItems );
    myComboLayout.setToolTipText("Select window layout");
    int height = myComboLayout.getPreferredSize().height;
    // +2 to give more room for Windows look and feel. Otherwise, icons are clipped. Not the best solution..
    myComboLayout.setPreferredSize( new Dimension(2*height+2,height+2) );
    myComboLayout.setMaximumSize( new Dimension(2*height+2,height+2) );
    myComboLayout.setRenderer( new ComboLayoutRenderer() );
    if( properties.windowLayout >= 0 && properties.windowLayout < csDockPaneManager.getNumLayoutOptions() ) {
      myComboLayout.setSelectedIndex( properties.windowLayout );
      if( properties.windowLayout == csDockPaneManager.LAYOUT_TABS ) {
        myButtonSelectPanes.setEnabled(false);
      }
    }

    add(myBoxOpenInNewPane);
    addSeparator();
    add(myBoxLockScalar);
    add(myTextScalar);
    add(myButtonIncScaling);
    add(myButtonDecScaling);
    addSeparator();
    add(myButtonSyncSettings);
    addSeparator();
    add( myButtonBackwardSeismic );
    add( myButtonForwardSeismic );
    addSeparator();
    add( myButtonBeginSeismic );
    add( myButtonEndSeismic );
    addSeparator();
    add( myButtonSnapShot );
    add( myButtonSnapShotPane );
    addSeparator();
    add(myComboLayout);
    add(myButtonSelectPanes);
    add(Box.createHorizontalGlue());
    
    myComboLayout.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        Integer layoutOption = (Integer)myComboLayout.getSelectedItem();
        if( layoutOption != null ) {
          mySeaView.resetWindowsLayout( layoutOption, true );
          if( layoutOption == csDockPaneManager.LAYOUT_TABS ) {
            myButtonSelectPanes.setEnabled(false);
          }
          else {
            myButtonSelectPanes.setEnabled(true);
          }
        }
      }
    });
    myBoxLockScalar.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.setScalarLocked( myBoxLockScalar.isSelected() );
      }
    });
    myBoxOpenInNewPane.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.setOpenInNewPanel( myBoxOpenInNewPane.isSelected() );
      }
    });
    myButtonSelectPanes.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.showSelectWindowsDialog();
      }
    });
    myTextScalar.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        try {
          myCurrentScalar = Float.parseFloat( myTextScalar.getText() );
          fireUpdateScalarEvent( myCurrentScalar );
//          mySeisPane.setScalar( scalar );
        }
        catch( NumberFormatException exc ) {
        }
      }
    });
    myButtonSyncSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.updateSyncDispSettings( myButtonSyncSettings.isSelected() );
      }
    });
    myButtonDecScaling.addMouseListener( new MouseAdapter() {
      @Override
      public void mousePressed( MouseEvent e ) {
        decScaling();
        myDecTimerTask = new csTimerTaskDecrement( new csITimerTaskDecrementListener() {
          @Override
          public void decrement() {
            decScaling();
          }
        });
        Timer timer = new Timer();
        timer.scheduleAtFixedRate(myDecTimerTask, 500, 200);
      }
      @Override
      public void mouseReleased( MouseEvent e ) {
        myDecTimerTask.cancel();
        myDecTimerTask = null;
      }
    });
    myButtonIncScaling.addMouseListener( new MouseAdapter() {
      @Override
      public void mousePressed( MouseEvent e ) {
        incScaling();
        myIncTimerTask = new csTimerTaskIncrement( new csITimerTaskIncrementListener() {
          @Override
          public void increment() {
            incScaling();
          }
        });
        Timer timer = new Timer();
        timer.scheduleAtFixedRate(myIncTimerTask, 500, 200);
      }
      @Override
      public void mouseReleased( MouseEvent e ) {
        myIncTimerTask.cancel();
        myIncTimerTask = null;
      }
    });
    myButtonBackwardSeismic.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.seismicBackward();
      }
    });
    myButtonForwardSeismic.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.seismicForward();
      }
    });
    myButtonBeginSeismic.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.seismicBegin();
      }
    });
    myButtonEndSeismic.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.seismicEnd();
      }
    });
    myButtonSnapShot.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
       mySeaView.createSnapShot(false);
      }
    });
    myButtonSnapShotPane.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        mySeaView.createSnapShot(true);
      }
    });
  }
  public void addToolBarListener( csISeaViewToolBarListener listener ) {
    myListeners.add(listener);
  }
  public void removeToolBarListener( csISeaViewToolBarListener listener ) {
    myListeners.remove(listener);
  }
  private void fireUpdateScalarEvent( float scalar ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).updateDisplayScalar(scalar);
    }
  }
  public void setMoveButtonsEnabled( boolean atBegin, boolean atEnd ) {
    myButtonBackwardSeismic.setEnabled(!atBegin);
    myButtonBeginSeismic.setEnabled(!atBegin);
    myButtonForwardSeismic.setEnabled(!atEnd);
    myButtonEndSeismic.setEnabled(!atEnd);
  }
  public void setToolbarEnabled( boolean set ) {
    myButtonSnapShot.setEnabled(set);
    myButtonSnapShotPane.setEnabled(set);
  }
  public void incScaling() {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).incDisplayScalar();
    }
  }
  public void decScaling() {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).decDisplayScalar();
    }
  }

  public void updateScalar( float newScalar, boolean isScalar ) {
    if( isScalar ) {
      myTextScalar.setText( "" + newScalar );
    }
    myTextScalar.setEnabled(isScalar);
  }
  
  public void updateSelectPanes( boolean hiddenPanesExist ) {
    if( hiddenPanesExist ) {
      myButtonSelectPanes.setIcon( csSeaViewActions.getIcon(csSeaViewActions.SelectPanesAlertAction) );
      myButtonSelectPanes.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.SelectPanesAlertAction] );
    }
    else {
      myButtonSelectPanes.setIcon( csSeaViewActions.getIcon(csSeaViewActions.SelectPanesAction) );
      myButtonSelectPanes.setToolTipText( csSeaViewActions.ACTION_DESC[csSeaViewActions.SelectPanesAction] );
    }
  }
  //*****************************************************************************
  //*****************************************************************************
  //
  private class ComboLayoutRenderer extends DefaultListCellRenderer {
    public ComboLayoutRenderer() {
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
      Integer layoutOption = (Integer)value;
      if( value == null ) return this;
      label.setIcon( csDockPaneManager.getIcon( layoutOption ) );
      label.setText("");
      setToolTipText( csDockPaneManager.getDescription( layoutOption ) );
      return this;
    }
  }
}


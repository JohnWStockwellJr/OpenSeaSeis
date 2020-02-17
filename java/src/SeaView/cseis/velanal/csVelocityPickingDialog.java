/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.velanal;

import java.awt.*;
import javax.swing.*;
import cseis.general.csStandard;
import cseis.seisdisp.csHorizonAttr;
import cseis.seisdisp.csPickOverlay;
import cseis.swing.csColorButton;
import cseis.swing.csColorChangeListener;
import cseis.velocity.csVelFunction;
import cseis.velocity.csVelPickOverlay;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.DecimalFormat;
import java.util.ArrayList;

/**
 * Dialog window providing velocity picking functionality in seismic view.<br>
 * This class is almost an exact copy of class csPickingDialog. These two should ideally be combined in the future. 
 * <br>
 */
@SuppressWarnings("serial")
public class csVelocityPickingDialog extends JDialog {
  private final static Color[] COLORS = {
    Color.magenta,
    Color.orange,
    Color.cyan,
    Color.lightGray,
    Color.darkGray,
    Color.yellow,
    Color.green
  };
  private static int ITEM_COUNTER = 1;
  private final ArrayList<PickItem> myPickItemList;
  private final csVelPickOverlay myVelPickOverlay;

  private int myActiveIndex;
  private JPanel myPanelItems;
  private final JButton myButtonAdd;
  private final JButton myButtonRemove;
  private final JButton myButtonLoad;
  private final JButton myButtonSave;

  private final JButton myButtonClose;
  private final JButton myButtonHide;
  private final JButton myButtonReset;
  
  private final csIVelPickDialogListener myListener;
  
  public csVelocityPickingDialog( JFrame frame, csVelPickOverlay velPickOverlay, String name, csIVelPickDialogListener listener ) {
    super(frame,"Velocity picking for: " + name );
    myVelPickOverlay = velPickOverlay;
    myListener = listener;

    myActiveIndex = -1;

    myButtonAdd  = new JButton("Add");
    myButtonAdd.setToolTipText("Add pick object");
    myButtonRemove = new JButton("Remove");
    myButtonRemove.setToolTipText("Remove active pick object");
    myButtonLoad  = new JButton("Load...");
    myButtonLoad.setToolTipText("Load velocity picks/function from ASCII file");
    myButtonSave  = new JButton("Save...");
    myButtonSave.setToolTipText("Save selected velocity picks to external ASCII file");
    myButtonClose = new JButton("Close");
    myButtonClose.setToolTipText("Close dialog, end velocity picking mode");
    myButtonHide = new JButton("Hide");
    myButtonHide.setToolTipText("Hide dialog, stay in velocity picking mode");
    myButtonReset = new JButton("Reset");
    myButtonReset.setToolTipText("<html>Reset velocity picking (=re-open setup dialog)<br><b>This will delete all velocity picks!</b></html>");
    myButtonRemove.setEnabled(false);
    myButtonSave.setEnabled(false);
    
    JLabel labelHelp1 = new JLabel("Left/middle mouse: Add/delete pick");
    int fontSize = labelHelp1.getFont().getSize();
    labelHelp1.setFont( new Font(Font.SANS_SERIF, Font.ITALIC, fontSize-1) );
    myPickItemList = new ArrayList<PickItem>();
    myPanelItems = new JPanel( new GridBagLayout() );

    JPanel panelActions = new JPanel(new GridBagLayout());
    panelActions.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Actions"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelItems = new JPanel(new GridBagLayout());
    panelItems.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Velocity pick objects"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    myPanelItems.add( new JLabel(""), new GridBagConstraints(
        0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    
    panelItems.add( myPanelItems, new GridBagConstraints(
        0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    int xp = 0;
    panelActions.add( myButtonAdd, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelActions.add( myButtonRemove, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelActions.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelActions.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelActions.add( myButtonLoad, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelActions.add( myButtonSave, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );

    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 1, 1, 1, 0.1, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonReset, new GridBagConstraints(
        xp++, 1, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 1, 1, 1, 0.8, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonHide, new GridBagConstraints(
        xp++, 1, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 11 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 1, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 1, 1, 1, 0.1, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( labelHelp1, new GridBagConstraints(
        0, 0, xp, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );
    panelAll.add(panelActions,BorderLayout.NORTH);
    panelAll.add(panelItems,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);
    
    myButtonAdd.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        addItem();
        refresh();
      }
    });
    myButtonRemove.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        removeItem();
        refresh();
      }
    });
    myButtonSave.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        saveVelField();
      }
    });
    myButtonLoad.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        loadVelField();
      }
    });
    myButtonClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        if( myListener != null ) myListener.velPickDialogClosing();
        dispose();
      }
    });
    myButtonHide.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        dispose();
      }
    });
    myButtonReset.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        int option = JOptionPane.showConfirmDialog( csVelocityPickingDialog.this,
            "This will remove all current velocity picks!\n" +
            "Are you sure you want to reset the velocity analysis?", "Confirm reset", JOptionPane.YES_NO_OPTION );
        if( option != JOptionPane.YES_OPTION ) return;
        if( myListener != null ) {
          myListener.velPickDialogReset();
        }
        dispose();
      }
    });
    
    getContentPane().add(panelAll);
    setSize( new Dimension(400,280) );
    setLocationRelativeTo(frame);
  }
  public void refresh() {
    myPanelItems.removeAll();
    int xp = 0;
    int counter = 0;
    JLabel labelEdit  = new JLabel("Active");
    JLabel labelName  = new JLabel("Name");
    JLabel labelColor = new JLabel("Color");
    myPanelItems.add( labelEdit, new GridBagConstraints(
      xp++, counter, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 5, 0, 5, 0 ), 0, 0 ) );
    myPanelItems.add( labelName, new GridBagConstraints(
      xp++, counter, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 5, 5, 5, 0 ), 0, 0 ) );
    myPanelItems.add( labelColor, new GridBagConstraints(
      xp++, counter++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.NONE, new Insets( 5, 11, 5, 11 ), 0, 0 ) );
    for( int item = 0; item < myPickItemList.size(); item++ ) {
      PickItem pickItem = myPickItemList.get(item);
      xp = 0;
      myPanelItems.add( pickItem.checkActive, new GridBagConstraints(
        xp++, counter, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 0 ), 0, 0 ) );
      myPanelItems.add( pickItem.textName, new GridBagConstraints(
        xp++, counter, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 5, 0 ), 0, 0 ) );
      myPanelItems.add( pickItem.colorButton, new GridBagConstraints(
        xp++, counter++, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 11, 5, 11 ), 0, 0 ) );
    }
    myPanelItems.add( Box.createVerticalGlue(), new GridBagConstraints(
      0, counter++, 3, 1, 1.0, 1.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    revalidate();
    repaint();
  }
  //--------------------------------------------------------
  private void addItem() {
    addItem( null, null );
  }
  private void addItem( java.util.List<csVelFunction> velFuncList, String name ) {
    int numItems = myPickItemList.size();
    int colorIndex = numItems % COLORS.length;
    PickItem pickItem = new PickItem( COLORS[colorIndex] );
    myPickItemList.add( pickItem );
    if( name == null ) {
      name = pickItem.textName.getText();
    }
    else {
      pickItem.textName.setText(name);
    }
    csHorizonAttr attr = new csHorizonAttr( pickItem.id, name, pickItem.colorButton.getColor(), pickItem.comboType.getSelectedIndex() );
    if( velFuncList == null ) {
      myVelPickOverlay.addVelField(attr);
      myVelPickOverlay.refresh();
    }
    else {
      myVelPickOverlay.addVelField( attr, velFuncList );
      if( myVelPickOverlay.getNumVelFields() == 1 ) myVelPickOverlay.setActiveObject( attr.getID() );
      myVelPickOverlay.refresh();
    }
    setActiveIndex( numItems );
    updateButtons();
  }
  private void removeItem() {
    if( myActiveIndex < 0 ) return;
    int velFieldID = myPickItemList.get(myActiveIndex).id;
    myVelPickOverlay.removeVelField( velFieldID );
    myPickItemList.remove( myActiveIndex );
    setActiveIndex( Math.min( myActiveIndex, myPickItemList.size()-1 ) );
    updateButtons();
    if( myActiveIndex >= 0 ) myVelPickOverlay.setActiveObject( myPickItemList.get(myActiveIndex).id );
    myVelPickOverlay.refresh();
  }
  private void setActiveItem( PickItem pickItem ) {
    int index = getPickItemIndex( pickItem );
    if( index < 0 ) return;
    setActiveIndex( index );
  }
  private void updatePickAttr( PickItem pickItem ) {
    int index = getPickItemIndex( pickItem );
    if( index < 0 ) return;
    csHorizonAttr attr = new csHorizonAttr( pickItem.id, pickItem.textName.getText(), pickItem.colorButton.getColor(), pickItem.comboType.getSelectedIndex() );
    myVelPickOverlay.updateDisplayAttr( attr );
  }
  private int getPickItemIndex( PickItem pickItem ) {
    for( int i = 0; i < myPickItemList.size(); i++ ) {
      PickItem item = myPickItemList.get(i);
      if( pickItem.equals(item) ) {
        return i;
      }
    }
    return -1;
  }
  private void setActiveIndex( int index ) {
    if( myActiveIndex < myPickItemList.size() && myActiveIndex >= 0 ) {
      myPickItemList.get(myActiveIndex).setActive(false);
    }
    myActiveIndex = index;
    if( myActiveIndex >= 0 ) {
      myPickItemList.get(myActiveIndex).setActive(true);
      myVelPickOverlay.setActiveObject( myPickItemList.get(myActiveIndex).id );
      myVelPickOverlay.refresh();
    }
  }
  private void saveVelField() {
    myListener.savePicks( myVelPickOverlay.getVelField( myPickItemList.get(myActiveIndex).id ).values(), myPickItemList.get(myActiveIndex).textName.getText() );
  }
  private void loadVelField() {
    myListener.loadPicks();
  }
  public void load( java.util.List<csVelFunction> velFuncList, String name ) {
    addItem( velFuncList, name );
    refresh();
//    myListener.updatePicks();
  }
  private void updateButtons() {
    boolean doEnable = ( myPickItemList.size() > 0 && myActiveIndex >= 0 );
    myButtonRemove.setEnabled(doEnable);
    myButtonSave.setEnabled(doEnable);
  }
  //-------------------------------------------------------------
  public class PickItem {
    private final int   id;
    final JTextField    textName;
    final JRadioButton  checkActive;
    final csColorButton colorButton;
    final JComboBox<String>     comboType;
    PickItem() {
      this( Color.black );
    }
    PickItem( Color color ) {
      id = ITEM_COUNTER++;
      textName = new JTextField( "vel" + id );
      checkActive = new JRadioButton();
      checkActive.setSelected(false);
      colorButton = new csColorButton( csVelocityPickingDialog.this, color );
      comboType   = new JComboBox<String>();
      comboType.setModel( new DefaultComboBoxModel<String>(csPickOverlay.MODE_TEXT_FIELDS) );
      comboType.setSelectedIndex(0);
      checkActive.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          setActiveItem( PickItem.this );
        }
      });
      comboType.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          updatePickAttr( PickItem.this );
        }
      });
      colorButton.addColorChangeListener( new csColorChangeListener() {
        @Override
        public void colorChanged(Object obj, Color color) {
          updatePickAttr( PickItem.this );
        }
      });
    }
    boolean isActive() {
      return checkActive.isSelected();
    }
    void setActive( boolean doSet ) {
      checkActive.setSelected( doSet );
    }
    int getID() {
      return id;
    }
  }
}



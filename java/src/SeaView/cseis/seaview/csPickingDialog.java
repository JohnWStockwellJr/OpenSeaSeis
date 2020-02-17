/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import java.awt.*;
import javax.swing.*;
import cseis.general.csStandard;
import cseis.seis.csHeaderDef;
import cseis.seisdisp.csHorizonAttr;
import cseis.seisdisp.csPickOverlay;
import cseis.swing.csColorButton;
import cseis.swing.csColorChangeListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Dialog window providing picking functionality in seismic view.<br>
 * <br>
 * @author 2011 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csPickingDialog extends JDialog {
  private static int ITEM_COUNTER = 0;
  
  private final static Color[] COLORS = {
    Color.blue,
    Color.yellow,
    Color.red,
    Color.black,
    Color.green,
    Color.orange,
    Color.white
  };
  private JFileChooser myFileChooser;
  private final ArrayList<PickItem> myPickItemList;
  private final csPickOverlay myPickOverlay;

  private int myActiveIndex;
  private JPanel myPanelItems;
  private final JButton myButtonClose;
  private final JButton myButtonHide;
  private final JButton myButtonAdd;
  private final JButton myButtonRemove;

  private final JButton myButtonLoad;
  private final JButton myButtonHeader;
  private final JButton myButtonSave;
  private final csIPickDialogListener myListener;
  
  private final csHeaderDef[] myHeaderDef;
  
  public csPickingDialog( JFrame frame, csPickOverlay pickOverlay, String name, csIPickDialogListener listener, csHeaderDef[] headerDef ) {
    super(frame,"Horizon picking for: " + name );
    myPickOverlay = pickOverlay;
    myFileChooser = null;
    myListener = listener;
    myHeaderDef = headerDef;

    myActiveIndex = -1;

    myButtonAdd  = new JButton("Add");
    myButtonAdd.setToolTipText("Add pick object");
    myButtonRemove = new JButton("Remove");
    myButtonRemove.setToolTipText("Remove active pick object");
    myButtonHeader = new JButton("Header...");
    myButtonHeader.setToolTipText("Load time picks from trace header");
    myButtonLoad  = new JButton("Load...");
    myButtonLoad.setToolTipText("Load time picks from ASCII file");
    myButtonSave  = new JButton("Save...");
    myButtonSave.setToolTipText("Save selected time picks to external ASCII file");
    myButtonClose = new JButton("Close");
    myButtonClose.setToolTipText("Close picking dialog, end picking mode");
    myButtonHide = new JButton("Hide");
    myButtonHide.setToolTipText("Hide picking dialog, stay in picking mode");
    myButtonRemove.setEnabled(false);
    myButtonSave.setEnabled(false);
    
    JLabel labelHelp1 = new JLabel("Left/middle mouse: Add/delete pick. Press down CTRL: Interpolate");
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
        BorderFactory.createTitledBorder("Time pick objects"),
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
    panelActions.add( myButtonHeader, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelActions.add( myButtonLoad, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelActions.add( myButtonSave, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );

    getRootPane().setDefaultButton(myButtonClose);

    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 1, 1, 1, 0.1, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonHide, new GridBagConstraints(
        xp++, 1, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 1, 1, 1, 0.8, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 1, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 1, 1, 1, 0.1, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( labelHelp1, new GridBagConstraints(
        0, 0, xp, 1, 1.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//    panelButtons.add( labelHelp2, new GridBagConstraints(
//        0, 1, xp, 1, 1.0, 0.0, GridBagConstraints.CENTER,
//        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );
    panelAll.add(panelActions,BorderLayout.NORTH);
    panelAll.add(panelItems,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);

    addItem();
    refresh();
    
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
        save();
      }
    });
    myButtonLoad.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        load();
      }
    });
    myButtonHeader.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        setFromHeader();
      }
    });
    myButtonClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        if( myListener != null ) myListener.closePickDialog();
        myPickOverlay.deactivate();
        dispose();
      }
    });
    myButtonHide.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        dispose();
      }
    });
    
    getContentPane().add(panelAll);
    setSize( new Dimension(500,300) );
    setLocationRelativeTo(frame);
  }
  private void refresh() {
    myPanelItems.removeAll();
    for( int item = 0; item < myPickItemList.size(); item++ ) {
      PickItem pickItem = myPickItemList.get(item);
      int xp = 0;
      myPanelItems.add( pickItem.checkActive, new GridBagConstraints(
        xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 0 ), 0, 0 ) );
      myPanelItems.add( pickItem.textName, new GridBagConstraints(
        xp++, item, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 5, 0 ), 0, 0 ) );
      myPanelItems.add( pickItem.colorButton, new GridBagConstraints(
        xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 11, 5, 11 ), 0, 0 ) );
      myPanelItems.add( pickItem.comboType, new GridBagConstraints(
        xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 5, 0 ), 0, 0 ) );
    }
    myPanelItems.add( Box.createVerticalGlue(), new GridBagConstraints(
      0, myPickItemList.size(), 4, 1, 1.0, 1.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    revalidate();
    repaint();
  }
  //--------------------------------------------------------
  private void addItem() {
    int numItems = myPickItemList.size();
    int colorIndex = numItems % COLORS.length;
    PickItem pickItem = new PickItem( COLORS[colorIndex] );
    myPickItemList.add( pickItem );
    csHorizonAttr attr = new csHorizonAttr( pickItem.getID(), pickItem.textName.getText(), pickItem.colorButton.getColor(), pickItem.comboType.getSelectedIndex() );
    myPickOverlay.addHorizon( attr );
    setActiveIndex( numItems );
    updateButtons();
  }
  private void removeItem() {
    if( myActiveIndex < 0 ) return;
    myPickOverlay.removeHorizon( myPickItemList.get(myActiveIndex).id );
    myPickItemList.remove( myActiveIndex );
    setActiveIndex( Math.min( myActiveIndex, myPickItemList.size()-1 ) );
    updateButtons();
    if( myListener != null ) myListener.updatePicks();
  }
  private void setActiveItem( PickItem pickItem ) {
    int index = getPickItemIndex( pickItem );
    if( index < 0 ) return;
    setActiveIndex( index );
  }
  private void updatePickAttr( PickItem pickItem ) {
    int index = getPickItemIndex( pickItem );
    if( index < 0 ) return;
    csHorizonAttr attr = new csHorizonAttr( pickItem.getID(), pickItem.textName.getText(), pickItem.colorButton.getColor(), pickItem.comboType.getSelectedIndex() );
    myPickOverlay.updateDisplayAttr( attr );
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
      myPickOverlay.setActiveObject( myPickItemList.get(myActiveIndex).id );
    }
  }
  private void save() {
//    HeaderSelectDialog dialog = new HeaderSelectDialog( this, "Select header to save with picks:", false, new IHeaderSelectListener() {
//      @Override
//      public void selectHeader( csHeaderDef hdef ) {
//        float[] picks = myPickOverlay.getPicks( myPickItemList.get(myActiveIndex).id );
        myListener.savePicks( myPickOverlay.getHorizon( myPickItemList.get(myActiveIndex).id ), myPickItemList.get(myActiveIndex).textName.getText() );
//      }
//    });
//    dialog.setVisible(true);
  }
  private void load() {
    if( myFileChooser == null ) myFileChooser = new JFileChooser();
    int option = myFileChooser.showOpenDialog( this );
    if( option == JFileChooser.APPROVE_OPTION ) {
      File file = myFileChooser.getSelectedFile();
      double sampleInt = myPickOverlay.getSampleInt();
      try {
        BufferedReader reader = new BufferedReader( new FileReader(file) );
        String line;
        ArrayList<Integer> traceNumList = new ArrayList<Integer>();
        ArrayList<Float> sampleIndexList = new ArrayList<Float>();
        while( (line = reader.readLine()) != null ) {
          String[] tokens = line.trim().split(" ");
          if( tokens.length >= 2 ) {
            int traceNumber = Integer.parseInt( tokens[0] );
            double time = Double.parseDouble( tokens[1] );
            traceNumList.add( traceNumber );
            sampleIndexList.add( (float)(time/sampleInt) );
          }
        }
        reader.close();
        String name = file.getAbsolutePath();
        int separatorIndex = name.lastIndexOf( java.io.File.separatorChar );
        if( separatorIndex < 0 ) separatorIndex = 0;
        int pointIndex = name.lastIndexOf('.');
        if( pointIndex >= 0 ) name = name.substring( separatorIndex+1, pointIndex );
        if( myPickOverlay.loadActivePicks( traceNumList, sampleIndexList, name ) ) {
          if( myActiveIndex >= 0 ) myPickItemList.get(myActiveIndex).textName.setText(name);
          myListener.updatePicks();
        }
        else {
          JOptionPane.showMessageDialog(this, "Unknown error (bug?) occurred when trying to set time picks", "Error", JOptionPane.ERROR_MESSAGE);
        }
        JOptionPane.showMessageDialog(this,
          "Loaded horizon picks from file\n'" + file.getAbsolutePath() +
            "'\nMin/max trace in file: " + traceNumList.get(0) + "/" + traceNumList.get(traceNumList.size()-1),
          "Info", JOptionPane.INFORMATION_MESSAGE);
      }
      catch (IOException ex) {
        JOptionPane.showMessageDialog(this, "Error occurred when loading time picks:\n" + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
      }
      
    }
  }
  private void setFromHeader() {
    HeaderSelectDialog dialog = new HeaderSelectDialog( this, "Set picks from header:", true, new IHeaderSelectListener() {
      @Override
      public void selectHeader(csHeaderDef hdef) {
        myListener.setPicksFromHeader( hdef );
      }
    });
    dialog.setVisible(true);
  }
  private void updateButtons() {
    boolean doEnable = ( myPickItemList.size() > 0 && myActiveIndex >= 0 );
    myButtonRemove.setEnabled(doEnable);
    myButtonSave.setEnabled(doEnable);
  }
  //-----------------------------------------------------------------------------
  //
  public class PickItem {
    private final int   id;
    final JTextField    textName;
    final JRadioButton  checkActive;
    final csColorButton colorButton;
    final JComboBox<String>  comboType;
    PickItem() {
      this( Color.black );
    }
    PickItem( Color color ) {
      id = ITEM_COUNTER++;
      textName = new JTextField( "pick" + id );
      checkActive = new JRadioButton();
      checkActive.setSelected(false);
      colorButton = new csColorButton( csPickingDialog.this, color );
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
  //-------------------------------------------------------------
  interface IHeaderSelectListener {
    public void selectHeader( csHeaderDef hdef );
  }
  class HeaderSelectDialog extends JDialog {
    JComboBox<csHeaderDef> headerBox;
    JButton buttonApply;
    JButton buttonClose;
    IHeaderSelectListener myHeaderSelectionListener;
    boolean myIsApplyButton;

    HeaderSelectDialog( JDialog parent, String title, boolean isApplyButton, IHeaderSelectListener listener ) {
      super( parent, "Header selection", true );
      myIsApplyButton = isApplyButton;
      myHeaderSelectionListener = listener;
      headerBox = new JComboBox<csHeaderDef>( myHeaderDef );
      buttonApply = new JButton("Apply");
      buttonClose = new JButton("Close");
      
      JPanel panel1 = new JPanel(new GridBagLayout());
      panel1.setBorder( BorderFactory.createCompoundBorder(
          BorderFactory.createTitledBorder(""),
          csStandard.INNER_EMPTY_BORDER ) );

      panel1.add( headerBox, new GridBagConstraints(
          0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
          GridBagConstraints.HORIZONTAL, new Insets( 11, 5, 11, 5 ), 0, 0 ) );
      panel1.add( Box.createVerticalBox(), new GridBagConstraints(
          0, 1, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
          GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

      getRootPane().setDefaultButton(buttonClose);

      JPanel panelButtons = new JPanel( new GridBagLayout() );
      int xp = 0;
      panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
          xp++, 0, 1, 1, 0.9, 0.0, GridBagConstraints.CENTER,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      if( myIsApplyButton ) {
        panelButtons.add( buttonApply, new GridBagConstraints(
            xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
            GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
      }
      panelButtons.add( buttonClose, new GridBagConstraints(
          xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
          GridBagConstraints.HORIZONTAL, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
      panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
          xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.CENTER,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

      JPanel panelAll = new JPanel(new BorderLayout());
      panelAll.setBorder( csStandard.DIALOG_BORDER );
      panelAll.add( new JLabel(title,JLabel.CENTER),BorderLayout.NORTH);
      panelAll.add(panel1,BorderLayout.CENTER);
      panelAll.add(panelButtons,BorderLayout.SOUTH);
    
      getContentPane().add(panelAll);
      pack();
      setLocationRelativeTo( parent );

      buttonApply.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed( ActionEvent e ) {
          int index = headerBox.getSelectedIndex();
          if( index < 0 ) return;
          csHeaderDef hdrDef = myHeaderDef[ index ];
          myHeaderSelectionListener.selectHeader(hdrDef);
        }
      });
      buttonClose.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed( ActionEvent e ) {
          dispose();
          if( !myIsApplyButton ) {
            int index = headerBox.getSelectedIndex();
            csHeaderDef hdrDef = null;
            if( index >= 0 ) hdrDef = myHeaderDef[ index ];
            myHeaderSelectionListener.selectHeader(hdrDef);
          }
        }
      });
    }
  }
}



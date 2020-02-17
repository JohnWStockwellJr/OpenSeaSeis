/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.seisdisp.csHeaderOverlay;
import cseis.general.csStandard;
import cseis.seis.csHeaderDef;
import cseis.seisdisp.csSeisView;
import cseis.swing.csColorButton;
import cseis.swing.csColorChangeListener;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.Hashtable;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * Dialog window providing functionality to overlay trace header values on top of seismic view.<br>
 * <br>
 * Up to myNumItems (8) trace headers may be selected to be displayed as graphs on top of a seismic view.
 * The horizontal axis is defined by the trace number. The vertical axis is defined by the vertical domain.
 * Each graph point is placed at its trace and the vertical location corresponding to its trace header value.
 * Values may be shifted and scaled to fit them into the min/max range of the vertical domain.
 * @author 2010 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csOverlayDialog extends JDialog {
  private static final int TEXT_MIN_WIDTH = 50;
  private static final int MAX_TRACE_INTERVAL = 20;

  private csSeisView mySeisView;

  private JButton myButtonApply;
  private JButton myButtonClose;

  private int myNumItems;
  private JCheckBox[] myCheckItem;
  private JCheckBox myBoxShowLabels;
  private csColorButton[] myColorButton;
  private JTextField[] myTextAddValue;
  private JTextField[] myTextMultValue;
  private csHeaderOverlay[] myHeaderOverlay;
//  private JComboBox[] myComboHdrNames;
//  private JComboBox[] myComboType;
  private ArrayList<JComboBox<csHeaderDef>> myComboHdrNames;
  private ArrayList<JComboBox<String>> myComboType;

  csHeaderDef[] myTraceHeaderDef;

  private JTextField myTextSize;
  private JTextField myTextLineWidth;
  private JSlider mySliderTraceInterval;
  private JSlider mySliderTransparency;


  public csOverlayDialog( JFrame parentFrame, csSeisView seisview, csHeaderDef[] traceHeaders ) {
    super( parentFrame, "Trace header overlay" );
    setModal(false);

    mySeisView = seisview;
    myTraceHeaderDef = null;

    csHeaderOverlay.Attribute attr = new csHeaderOverlay().getAttributes();
    myTextSize      = new JTextField("" + attr.size);
    int height = myTextSize.getPreferredSize().height;
    myTextSize.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextLineWidth = new JTextField("" + attr.lineWidth);
    myTextLineWidth.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );

    myNumItems  = 8;
    myCheckItem     = new JCheckBox[myNumItems];
    myColorButton   = new csColorButton[myNumItems];
    myComboHdrNames = new ArrayList<JComboBox<csHeaderDef>>(myNumItems);
    myComboType     = new ArrayList<JComboBox<String>>(myNumItems);
//    myComboHdrNames = new JComboBox[myNumItems];
//    myComboType     = new JComboBox[myNumItems];
    myTextAddValue  = new JTextField[myNumItems];
    myTextMultValue = new JTextField[myNumItems];
    myHeaderOverlay = new csHeaderOverlay[myNumItems];
    myBoxShowLabels = new JCheckBox("Show labels");
    myBoxShowLabels.setSelected(true);
    myBoxShowLabels.setToolTipText("Plot trace header names next to overlays");

    ArrayList<Color> colorList = new ArrayList<Color>();
    colorList.add(Color.blue);
    colorList.add(Color.yellow);
    colorList.add(Color.red);
    colorList.add(Color.black);
    colorList.add(Color.green);
    colorList.add(Color.orange);
    colorList.add(Color.white);

    for( int item = 0; item < myNumItems; item++ ) {
      int colorIndex = item % colorList.size();
      Color color = colorList.get(colorIndex);
      myHeaderOverlay[item] = new csHeaderOverlay();
      myHeaderOverlay[item].setColor(color);

      myCheckItem[item]     = new JCheckBox();
      myCheckItem[item].setSelected(false);
      myComboHdrNames.add( new JComboBox<csHeaderDef>() );
      myTextAddValue[item]  = new JTextField("" + attr.addValue);
      myTextMultValue[item] = new JTextField("" + attr.multValue);
      myComboType.add( new JComboBox<String>() );
      myComboType.get(item).setModel( new DefaultComboBoxModel<String>(csHeaderOverlay.TYPE_TEXT_FIELDS) );
      myComboType.get(item).setSelectedIndex(0);

      myComboHdrNames.get(item).setPreferredSize( new Dimension(120,20) );
      myTextAddValue[item].setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
      myTextMultValue[item].setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );

      myColorButton[item]   = new csColorButton(this,color,item);
      myColorButton[item].addColorChangeListener( new csColorChangeListener() {
        public void colorChanged( Object obj, Color color ) {
          int item = ((csColorButton)obj).getID();
          myHeaderOverlay[item].setColor(myColorButton[item].getColor());
          mySeisView.repaint();
        }
      });
      myTextAddValue[item].addActionListener( new ActionListener() {
        public void actionPerformed( ActionEvent e ) {
          JTextField field = (JTextField)e.getSource();
          int itemCurrent = -1;
          for( int item = 0; item < myNumItems; item++ ) {
            if( myTextAddValue[item].equals(field) ) {
              itemCurrent = item;
              break;
            }
          }
          if( itemCurrent < 0 ) return;
          applyItem( itemCurrent, true );
        }
      });
      myTextMultValue[item].addActionListener( new ActionListener() {
        public void actionPerformed( ActionEvent e ) {
          JTextField field = (JTextField)e.getSource();
          int itemCurrent = -1;
          for( int item = 0; item < myNumItems; item++ ) {
            if( myTextMultValue[item].equals(field) ) {
              itemCurrent = item;
              break;
            }
          }
          if( itemCurrent < 0 ) return;
          applyItem( itemCurrent, true );
        }
      });
      myCheckItem[item].addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
          JCheckBox box = (JCheckBox)e.getSource();
          int itemCurrent = -1;
          for( int item = 0; item < myNumItems; item++ ) {
            if( myCheckItem[item].equals(box) ) {
              itemCurrent = item;
              break;
            }
          }
          if( itemCurrent < 0 ) return;
          applyItem(itemCurrent,false);
          mySeisView.repaint();
/*
          if( e.getStateChange() == ItemEvent.SELECTED ) {
            mySeisView.removeOverlay(myHeaderOverlay[itemCurrent]);
          }
          else {
            mySeisView.addOverlay(myHeaderOverlay[itemCurrent]);
          }
 */
        }
      });
      myComboHdrNames.get(item).addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
          JComboBox<?> box = (JComboBox<?>)e.getSource();
          int itemCurrent = -1;
          for( int item = 0; item < myNumItems; item++ ) {
            if( myComboHdrNames.get(item).equals(box) ) {
              itemCurrent = item;
              break;
            }
          }
          if( itemCurrent < 0 ) return;
          if( e.getStateChange() == ItemEvent.SELECTED ) {
            applyItem(itemCurrent,true);
          }
        }
      });
      myComboType.get(item).addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
          JComboBox<?> box = (JComboBox<?>)e.getSource();
          int itemCurrent = -1;
          for( int item = 0; item < myNumItems; item++ ) {
            if( myComboType.get(item).equals(box) ) {
              itemCurrent = item;
              break;
            }
          }
          if( itemCurrent < 0 ) return;
          if( e.getStateChange() == ItemEvent.SELECTED ) {
            myHeaderOverlay[itemCurrent].setType( myComboType.get(itemCurrent).getSelectedIndex() );
            mySeisView.repaint();
          }
        }
      });

    //=============================================================================
    }
    myBoxShowLabels.addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
          apply();
          mySeisView.repaint();
        }
      });

    JPanel panelSelect = createSelectionPanel();
    int xp = 0;

//-----------------
// Trace interval panel
    JPanel panelInterval = new JPanel(new GridBagLayout());
    panelInterval.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Trace interval"),
        csStandard.INNER_EMPTY_BORDER ) );
    JPanel panelTrans = new JPanel(new GridBagLayout());
    panelTrans.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Transparency"),
        csStandard.INNER_EMPTY_BORDER ) );
    JPanel panelOther = new JPanel(new GridBagLayout());
    panelOther.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Other settings"),
        csStandard.INNER_EMPTY_BORDER ) );

    mySliderTraceInterval = new JSlider( JSlider.HORIZONTAL, 1, csOverlayDialog.MAX_TRACE_INTERVAL, 1 );
    mySliderTraceInterval.setToolTipText("Trace interval between travel time curve points");
    Hashtable<Integer,JLabel> labelsSliderInterval = new Hashtable<Integer,JLabel>();
    labelsSliderInterval.put( 1, new JLabel("1") );
    for( int i = 5; i <= csOverlayDialog.MAX_TRACE_INTERVAL; i += 5 ) {
      labelsSliderInterval.put( i, new JLabel(""+i) );
    }
    mySliderTraceInterval.setLabelTable( labelsSliderInterval );
    mySliderTraceInterval.setMinorTickSpacing(1);
    mySliderTraceInterval.setMajorTickSpacing(5);
    mySliderTraceInterval.setPaintTicks(true);
    mySliderTraceInterval.setPaintLabels(true);
    mySliderTraceInterval.setSnapToTicks(true);

    mySliderTransparency = new JSlider( JSlider.HORIZONTAL, 50, 255, 255 );
    mySliderTransparency.setToolTipText("Transparency of trace header symbol/line");
    mySliderTransparency.setMajorTickSpacing( 255 - 50 );
    mySliderTransparency.setPaintTicks(true);
    mySliderTransparency.setToolTipText("Brightness/opacity of symbols/lines");

    Hashtable<Integer,JLabel> labelsSlider = new Hashtable<Integer,JLabel>();
    labelsSlider.put( 50, new JLabel("Dim") );
    labelsSlider.put( 255, new JLabel("Bright") );
    mySliderTransparency.setLabelTable( labelsSlider );
    mySliderTransparency.setPaintLabels(true);
    mySliderTransparency.setSnapToTicks(false);

    panelInterval.add( mySliderTraceInterval, new GridBagConstraints(
        0, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelTrans.add( mySliderTransparency, new GridBagConstraints(
        0, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelOther.add( new JLabel("Symbol size:"), new GridBagConstraints(
        0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 5 ), 0, 0 ) );
    panelOther.add( myTextSize, new GridBagConstraints(
        1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelOther.add( new JLabel("Line width:"), new GridBagConstraints(
        0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 15, 10, 0, 5 ), 0, 0 ) );
    panelOther.add( myTextLineWidth, new GridBagConstraints(
        1, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 15, 0, 0, 0 ), 0, 0 ) );
    panelOther.add( myBoxShowLabels, new GridBagConstraints(
        0, 2, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 15, 10, 0, 0 ), 0, 0 ) );
    panelOther.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, 3, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    //  ----------------------------
    //  Button panel
    myButtonApply   = new JButton("Apply");
    myButtonClose = new JButton("Close");
    this.getRootPane().setDefaultButton(myButtonApply);

    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.4, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonApply, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.4, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelCenter = new JPanel( new GridBagLayout() );

    panelCenter.add( panelInterval, new GridBagConstraints(
        0, 0, 1, 1, 0.9, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelCenter.add( panelTrans, new GridBagConstraints(
        0, 1, 1, 1, 0.9, 0.1, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelCenter.add( panelOther, new GridBagConstraints(
        1, 0, 1, 2, 0.1, 0.1, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelCenter.add( Box.createVerticalGlue(), new GridBagConstraints(
        1, 2, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     JPanel panelAll = new JPanel(new BorderLayout());
     panelAll.setBorder( csStandard.DIALOG_BORDER );

     JPanel panelUpper = new JPanel(new BorderLayout());
     panelUpper.add(panelSelect,BorderLayout.NORTH);
     panelUpper.add(panelCenter,BorderLayout.SOUTH);


     panelAll.add(panelUpper,BorderLayout.NORTH);
     panelAll.add(panelButtons,BorderLayout.SOUTH);
     getContentPane().add(panelAll);

     setTraceHeaders( traceHeaders );
     pack();
     setLocationRelativeTo( parentFrame );

    //------------------------------------------------------------
    //
    myButtonApply.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    myButtonClose.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        cancel();
      }
    });
    myTextSize.addActionListener( new ActionListener() {
      public void actionPerformed( ActionEvent e ) {
        apply();
      }
    });
    myTextLineWidth.addActionListener( new ActionListener() {
      public void actionPerformed( ActionEvent e ) {
        apply();
      }
    });
    mySliderTraceInterval.addChangeListener( new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        apply();
      }
    });
    mySliderTransparency.addChangeListener( new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        apply();
      }
    });

    this.addWindowListener(new WindowAdapter() {
      public void windowClosing(WindowEvent e) {
        cancel();
      }
    });
  }
  private JPanel createSelectionPanel() {
    JPanel panelSelect = new JPanel(new GridBagLayout());
    panelSelect.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Trace header selection"),
        csStandard.INNER_EMPTY_BORDER ) );

    for( int item = 0; item < myNumItems; item++ ) {
      int xp = 0;
      panelSelect.add( myCheckItem[item], new GridBagConstraints(
          xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
      panelSelect.add( myComboHdrNames.get(item), new GridBagConstraints(
          xp++, item, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      panelSelect.add( new JLabel(" * "), new GridBagConstraints(
          xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
      panelSelect.add( myTextMultValue[item], new GridBagConstraints(
          xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      panelSelect.add( new JLabel(" + "), new GridBagConstraints(
          xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
      panelSelect.add( myTextAddValue[item], new GridBagConstraints(
          xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      panelSelect.add( myColorButton[item], new GridBagConstraints(
          xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
      panelSelect.add( myComboType.get(item), new GridBagConstraints(
          xp++, item, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
    return panelSelect;
  }

  private void cancel() {
    dispose();
  }
  //----------------------------------------------------------------------
  //
  private boolean apply() {
    int trcInterval = mySliderTraceInterval.getValue();
    if( trcInterval < 1 ) trcInterval = 1;
    if( trcInterval > csOverlayDialog.MAX_TRACE_INTERVAL ) trcInterval = csOverlayDialog.MAX_TRACE_INTERVAL;
    int transparency = mySliderTransparency.getValue();
    if( transparency < 50 ) transparency = 50;
    if( transparency > 255 ) transparency = 255;
    boolean showLabel = myBoxShowLabels.isSelected();
    int symbolSize  = 0;
    float lineWidth = 0;
    try {
      symbolSize = Integer.parseInt(myTextSize.getText());
      lineWidth  = Float.parseFloat(myTextLineWidth.getText());
    }
    catch( NumberFormatException exc ) {
      JOptionPane.showMessageDialog(this,
        "Inconsistent entry in text field",
        "Error message",
        JOptionPane.ERROR_MESSAGE);
      return false;
    }
    for( int item = 0; item < myNumItems; item++ ) {
      if( applyItem( item, false ) ) {
        myHeaderOverlay[item].setSymbolSize(symbolSize);
        myHeaderOverlay[item].setLineWidth(lineWidth);
        myHeaderOverlay[item].setTraceInterval(trcInterval);
        myHeaderOverlay[item].setTransparency(transparency);
        myHeaderOverlay[item].setShowLabel(showLabel);
      }
    }
    mySeisView.repaint();
    return true;
  }
  //--------------------------------------------------------------------------
  //
  private boolean applyItem( int item, boolean refresh ) {
    boolean isSelected = myCheckItem[item].isSelected();
    if( !isSelected ) {
      mySeisView.removeOverlay(myHeaderOverlay[item]);
      return false;
    }
    else {
      mySeisView.addOverlay(myHeaderOverlay[item]);
    }
    if( myComboHdrNames.get(item).getSelectedIndex() < 0 ) {
//      JOptionPane.showMessageDialog( this,
//        "No valid trace header selected in drop-down box #" + (item+1) + ".\n",
//        "Error message",
//        JOptionPane.ERROR_MESSAGE );
      return false;
    }
    try {
      float multValue = Float.parseFloat(myTextMultValue[item].getText());
      float addValue  = Float.parseFloat(myTextAddValue[item].getText());
      myHeaderOverlay[item].setLayoutAttributes(multValue,addValue);
    }
    catch( NumberFormatException exc ) {
      JOptionPane.showMessageDialog(this,
        "Inconsistent entry in one of the text fields",
        "Error message",
        JOptionPane.ERROR_MESSAGE);
      return false;
    }
    String hdrName = myComboHdrNames.get(item).getSelectedItem().toString();
    for( int ihdr = 0; ihdr < myTraceHeaderDef.length; ihdr++ ) {
      if( myTraceHeaderDef[ihdr].name.compareTo( hdrName ) == 0 ) {
        myHeaderOverlay[item].setHeader( hdrName, ihdr );
        break;
      }
    }
    if( refresh ) mySeisView.repaint();
    return true;
  }
  //-------------------------------------------------------------------
  //
  public void setTraceHeaders( csHeaderDef[] traceHeaders ) {
    if( traceHeaders == null ) {
      return;
    }
    if( myTraceHeaderDef == traceHeaders ) return;
//    if( myComboHdrNames[0].getModel() != null && myComboHdrNames[0].getModel().getSize() == traceHeaders.length ) return;
    myTraceHeaderDef = traceHeaders;
    csHeaderDef[] newHeaders = new csHeaderDef[traceHeaders.length];
    for( int ihdr = 0; ihdr < traceHeaders.length; ihdr++ ) {
      newHeaders[ihdr] = new csHeaderDef( traceHeaders[ihdr] );
    }
    java.util.Arrays.sort( newHeaders );
    for( int item = 0; item < myNumItems; item++ ) {
      if( myComboHdrNames.get(item) != null && myComboHdrNames.get(item).getSelectedIndex() >= 0 ) {
        String hdrName = myComboHdrNames.get(item).getSelectedItem().toString();
        myComboHdrNames.get(item).setModel( new DefaultComboBoxModel<csHeaderDef>(newHeaders) );
        myComboHdrNames.get(item).setSelectedIndex(-1);
        for( int ihdr = 0; ihdr < myTraceHeaderDef.length; ihdr++ ) {
          if( myTraceHeaderDef[ihdr].name.compareTo( hdrName ) == 0 ) {
            myHeaderOverlay[item].setHeader(hdrName, ihdr);
            for( int ihdr2 = 0; ihdr2 < newHeaders.length; ihdr2++ ) {
              if( newHeaders[ihdr2].name.compareTo( hdrName ) == 0 ) {
                myComboHdrNames.get(item).setSelectedIndex( ihdr2 );
                break;
              }
            }
            break;
          }
        } // END: for ihdr
      }
      else {
        myComboHdrNames.get(item).setModel( new DefaultComboBoxModel<csHeaderDef>(newHeaders) );
        myComboHdrNames.get(item).setSelectedIndex(-1);
      }
    } // END: for item
  }

}



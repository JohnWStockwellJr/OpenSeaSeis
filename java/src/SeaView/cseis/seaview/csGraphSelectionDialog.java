/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.seisdisp.csHeaderOverlay;
import cseis.general.csStandard;
import cseis.graph.csCurve;
import cseis.graph.csCurveAttributes;
import cseis.graph.csGraph2D;
import cseis.graph.csGraphAttributes;
import cseis.seis.csHeaderDef;
import cseis.swing.csColorButton;
import cseis.swing.csColorChangeListener;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import javax.swing.*;

/**
 * Dialog window providing functionality to display trace header graph on top of seismic view.<br>
 * <br>
 * 
 * @author 2011 Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csGraphSelectionDialog extends JDialog {
  private static final int TEXT_MIN_WIDTH = 40;

  private csSeisPaneBundle mySeisPaneBundle;
  private csGraph2D myGraph;

  private JButton myButtonApply;
  private JButton myButtonClose;
  private JButton myButtonHide;
  
  private int myNumItems;
  private JCheckBox[] myCheckItem;
  private csColorButton[] myColorButton;
  private ArrayList<JComboBox<csHeaderDef>> myComboHdrNames;
  private ArrayList<JComboBox<String>> myComboType;
  //  private JComboBox<String>[] myComboHdrNames;
 // private JComboBox<String>[] myComboType;
  private JTextField[] myTextSize;
  /// Unique identifier for each curve
  private int[] myCurveIDs;

  csHeaderDef[] myTraceHeaderDef;
  private boolean myIsUpdating;

  public csGraphSelectionDialog( JFrame parentFrame, csSeisPaneBundle seisPaneBundle, csGraph2D graph, csHeaderDef[] traceHeaders ) {
    super( parentFrame, "Trace header graph" );
    setModal(false);

    mySeisPaneBundle  = seisPaneBundle;
    myGraph    = graph;
    myTraceHeaderDef = null;
    myIsUpdating = false;

    myNumItems = 8;
    myCheckItem     = new JCheckBox[myNumItems];
    myColorButton   = new csColorButton[myNumItems];
    myComboHdrNames = new ArrayList<JComboBox<csHeaderDef>>(myNumItems);
    myComboType     = new ArrayList<JComboBox<String>>(myNumItems);
//    myComboHdrNames = new JComboBox[myNumItems];
//    myComboType     = new JComboBox[myNumItems];
    myTextSize      = new JTextField[myNumItems];
    myCurveIDs      = new int[myNumItems];

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

      myCheckItem[item]     = new JCheckBox();
      myCheckItem[item].setSelected(false);
      myComboHdrNames.add( new JComboBox<csHeaderDef>() );
      myComboType.add( new JComboBox<String>() );
      myComboType.get(item).setModel( new DefaultComboBoxModel<String>(csHeaderOverlay.TYPE_TEXT_FIELDS) );
      myComboType.get(item).setSelectedIndex(0);
      myTextSize[item] = new JTextField("2");
      myTextSize[item].setPreferredSize( new Dimension( TEXT_MIN_WIDTH, myTextSize[item].getPreferredSize().height ) );
      myComboHdrNames.get(item).setPreferredSize( new Dimension(120,20) );
      myCurveIDs[item] = -1;

      myColorButton[item]   = new csColorButton(this,color,item);
      myColorButton[item].addColorChangeListener( new csColorChangeListener() {
        public void colorChanged( Object obj, Color color ) {
          int item = ((csColorButton)obj).getID();
          csCurveAttributes attr = myGraph.getCurveAttributes(myCurveIDs[item]);
          if( attr != null ) {
            attr.lineColor  = color;
            attr.pointColor = color;
            myGraph.repaint();
          }
        }
      });
      myTextSize[item].addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          JTextField field = (JTextField)e.getSource();
          int itemCurrent = -1;
          for( int item = 0; item < myNumItems; item++ ) {
            if( myTextSize[item].equals(field) ) {
              itemCurrent = item;
              break;
            }
          }
          if( itemCurrent < 0 ) return;
          csCurveAttributes attr = myGraph.getCurveAttributes(myCurveIDs[itemCurrent]);
          if( attr != null ) {
            try {
              attr.lineSize = Integer.parseInt( myTextSize[itemCurrent].getText() );
              myGraph.repaint();
            }
            catch( Exception e2 ) {
              // Nothing
            }
          }
        }
      });
      myCheckItem[item].addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
          if( myIsUpdating ) return;
          JCheckBox box = (JCheckBox)e.getSource();
          int itemCurrent = -1;
          for( int item = 0; item < myNumItems; item++ ) {
            if( myCheckItem[item].equals(box) ) {
              itemCurrent = item;
              break;
            }
          }
          if( itemCurrent < 0 ) return;
          boolean isSelected = myCheckItem[itemCurrent].isSelected();
          if( !isSelected ) {
            removeCurve(itemCurrent);
          }
          else if( myComboHdrNames.get(itemCurrent).getSelectedIndex() < 0 ) {
            return;
          }
          else {
            addCurve(itemCurrent);
          }
        }
      });
      myComboHdrNames.get(item).addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
          if( myGraph == null ) return;
          if( myIsUpdating ) return;
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
            boolean isSelected = myCheckItem[itemCurrent].isSelected();
            if( !isSelected ) {
              return;
            }
            else {
              addCurve(itemCurrent);
            }
          }
          else if(e.getStateChange() == ItemEvent.DESELECTED) {
            removeCurve(itemCurrent);
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
            boolean isSelected = myCheckItem[itemCurrent].isSelected();
            if( isSelected && myComboHdrNames.get(itemCurrent).getSelectedIndex() < 0 ) {
//              myHeaderOverlay[itemCurrent].setType( myComboType[itemCurrent].getSelectedIndex() );
            }
          }
        }
      });
    }

    JPanel panelSelect = createSelectionPanel();
    int xp = 0;

    //  ----------------------------
    //  Button panel
    myButtonHide    = new JButton("Hide");
    myButtonHide.setToolTipText("Hide trace header graph");
    myButtonApply   = new JButton("Apply");
    myButtonApply.setToolTipText("Apply changes");
    myButtonClose = new JButton("Close");
    myButtonClose.setToolTipText("Close dialog window");
    this.getRootPane().setDefaultButton(myButtonApply);

    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonHide, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonApply, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.25, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     JPanel panelAll = new JPanel(new BorderLayout());
     panelAll.setBorder( csStandard.DIALOG_BORDER );
     
     panelAll.add(panelSelect,BorderLayout.NORTH);
     panelAll.add(panelButtons,BorderLayout.SOUTH);
     getContentPane().add(panelAll);

     setTraceHeaders( traceHeaders );
     pack();
     setLocationRelativeTo( parentFrame );

    //------------------------------------------------------------
    //
    myButtonHide.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        mySeisPaneBundle.showGraphPanel( false );
      }
    });
    myButtonApply.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    myButtonClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        cancel();
      }
    });

    this.addWindowListener(new WindowAdapter() {
      @Override
      public void windowClosing(WindowEvent e) {
        cancel();
      }
    });
  }
  private void removeCurve( int item ) {
    myGraph.removeCurve(myCurveIDs[item]);
    myCurveIDs[item] = -1;
    updateYLabel();
  }
  public void removeAllCurves() {
    for( int item = 0; item < myNumItems; item++ ) {
      myGraph.removeCurve(myCurveIDs[item]);
      myCurveIDs[item] = -1;
    }
    updateYLabel();
  }
  private void addCurve( int item ) {
    String hdrName = myComboHdrNames.get(item).getSelectedItem().toString();
    csCurveAttributes attr = new csCurveAttributes(hdrName);
    attr.filledType = csCurveAttributes.FILLED_TYPE_NONE;
    attr.pointType  = csCurveAttributes.POINT_TYPE_NONE;
    attr.lineColor  = myColorButton[item].getColor();
    attr.pointColor = myColorButton[item].getColor();
    try {
      attr.lineSize = Integer.parseInt( myTextSize[item].getText() );
    }
    catch( Exception e ) {
      // Nothing
    }
    csCurve curve   = mySeisPaneBundle.createCurve(hdrName,attr);
    if( curve == null ) return;
    myCurveIDs[item] = curve.curveID;
    myGraph.addCurve(curve);
    updateYLabel();
  }
//  private void updateCurve( int item ) {
//    String hdrName = myComboHdrNames.get(item).getSelectedItem().toString();
//    csCurveAttributes attr = new csCurveAttributes(hdrName);
//    csCurve curve   = mySeisPaneBundle.createCurve(hdrName,attr);
//    if( curve == null ) return;
//    myGraph.updateCurve(item, curve.dataX, curve.dataY, true);
//  }
  private void updateYLabel() {
    csGraphAttributes graphAttr = myGraph.getGraphAttributes();
    graphAttr.yLabel = "";
    java.util.Iterator<csCurve> iter = myGraph.getCurves().iterator();
    while( iter.hasNext() ) {
      String name = iter.next().attr.name;
      if( graphAttr.yLabel.length() > 0 ) graphAttr.yLabel += ", ";
      graphAttr.yLabel += name;
    }
    myGraph.repaintAll();
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
          xp++, item, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      panelSelect.add( myColorButton[item], new GridBagConstraints(
          xp++, item, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.NONE, new Insets( 0, 5, 0, 5 ), 0, 0 ) );
      panelSelect.add( myTextSize[item], new GridBagConstraints(
          xp++, item, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//      panelSelect.add( myComboType[item], new GridBagConstraints(
//          xp++, item, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
//          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
    return panelSelect;
  }

  private void cancel() {
    dispose();
  }
  //----------------------------------------------------------------------
  //
  private boolean apply() {
    return true;
  }
  //-------------------------------------------------------------------
  //
  public void setTraceHeaders( csHeaderDef[] traceHeaders ) {
    if( traceHeaders == null ) {
      return;
    }
    if( myTraceHeaderDef == traceHeaders ) return;

    myIsUpdating = true;
    myTraceHeaderDef = traceHeaders;
    csHeaderDef[] newHeaders = new csHeaderDef[traceHeaders.length];
    for( int ihdr = 0; ihdr < traceHeaders.length; ihdr++ ) {
      newHeaders[ihdr] = new csHeaderDef( traceHeaders[ihdr] );
    }
    java.util.Arrays.sort( newHeaders );
    removeAllCurves();
    for( int item = 0; item < myNumItems; item++ ) {
      if( myComboHdrNames.get(item) != null && myComboHdrNames.get(item).getSelectedIndex() >= 0 ) {
        String hdrName = myComboHdrNames.get(item).getSelectedItem().toString();
        myComboHdrNames.get(item).setModel( new DefaultComboBoxModel<csHeaderDef>(newHeaders) );
        myComboHdrNames.get(item).setSelectedIndex(-1);
        boolean found = false;
        for( int ihdr = 0; ihdr < myTraceHeaderDef.length; ihdr++ ) {
          if( myTraceHeaderDef[ihdr].name.compareTo( hdrName ) == 0 ) {
            for( int ihdr2 = 0; ihdr2 < newHeaders.length; ihdr2++ ) {
              if( newHeaders[ihdr2].name.compareTo( hdrName ) == 0 ) {
                myCheckItem[item].setSelected(true);
                myComboHdrNames.get(item).setSelectedIndex( ihdr2 );
                addCurve(item);
                found = true;
                break;
              } // END: if compareTo hdrName
            } // END: for ihdr2
            break;
          } // END: if compareTo hdrName
        } // END: for ihdr
        if( !found ) { // Trace header no longer exists --> Remove curve from graph
          myCheckItem[item].setSelected(false);
        }
      }
      else {
        myComboHdrNames.get(item).setModel( new DefaultComboBoxModel<csHeaderDef>(newHeaders) );
        myComboHdrNames.get(item).setSelectedIndex(-1);
      }
    } // END: for item
    myIsUpdating = false;
    updateYLabel();
  }
}



/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import cseis.swing.csColorButton;
import cseis.swing.csColorChangeListener;

/**
 * Attributes defining style of 2D curve.
 * @author Bjorn Olofsson
 */
public class csCurveAttributes {
  public static final int POINT_TYPE_NONE   = 0;
  public static final int POINT_TYPE_CIRCLE = 1;
  public static final int POINT_TYPE_CROSS  = 2;
  public static final int POINT_TYPE_SQUARE = 3;
  public static final int LINE_TYPE_NONE   =  0;
  public static final int LINE_TYPE_SOLID   = 1;
  public static final int LINE_TYPE_DASH    = 2;
  public static final int LINE_TYPE_DOT     = 3;
  public static final int FILLED_TYPE_NONE  = 0;
  public static final int FILLED_TYPE_ZERO  = 1;
  public static final int FILLED_TYPE_MIN   = 2;

  public static final int[] POINT_TYPES = {
    POINT_TYPE_NONE,
    POINT_TYPE_CIRCLE,
    POINT_TYPE_CROSS,
    POINT_TYPE_SQUARE
  };
  public static final String[] TEXT_POINT_TYPES = {
    "None",
    "Circle",
    "Cross",
    "Square"
  };
  public static final int[] LINE_TYPES = {
    LINE_TYPE_NONE,
    LINE_TYPE_SOLID,
    LINE_TYPE_DASH,
    LINE_TYPE_DOT
  };
  public static final String[] TEXT_LINE_TYPES = {
    "None",
    "Solid",
    "Dashed",
    "Dotted"
  };
  public static final int[] FILL_TYPES = {
    FILLED_TYPE_NONE,
    FILLED_TYPE_ZERO,
    FILLED_TYPE_MIN
  };
  public static final String[] TEXT_FILL_TYPES = {
    "None",
    "Zero line",
    "Minimum"
  };
  
  public static final int Y_AXIS_LEFT   = 11;
  public static final int Y_AXIS_RIGHT  = 12;
  public static final int X_AXIS_TOP    = 13;
  public static final int X_AXIS_BOTTOM = 14;

  private static int staticCurveCounter = 1;

  public String name;
  public int pointType;
  public int lineType;
  public int filledType;
  public int pointSize;
  public int lineSize;
  public Color pointColor;
  public Color lineColor;
  public Color fillColor;
  public int axisXType;
  public int axisYType;

  private float myMinX;
  private float myMaxX;
  private float myMinY;
  private float myMaxY;
  private float myMaxTotal;
  private boolean myIsMinMaxSet;
  private JPanel myPanel;

  private ArrayList<csCurveAttrChangeListener> myListeners;
  
  public csCurveAttributes() {
    this( "Curve #" + csCurveAttributes.staticCurveCounter );
  }
  public csCurveAttributes( String name_in ) {
    csCurveAttributes.staticCurveCounter += 1;
    name    = name_in;
    myPanel = null;
    myListeners = new ArrayList<csCurveAttrChangeListener>();

    //    pointType = POINT_TYPE_NONE;
//    filledType = FILLED_TYPE_NONE;
    pointType  = POINT_TYPE_SQUARE;
    lineType   = LINE_TYPE_SOLID;
    filledType = FILLED_TYPE_MIN;
    pointSize  = 10;
    lineSize   = 1;
    pointColor = Color.red;
    lineColor  = Color.black;
    fillColor  = Color.green;
    axisXType  = X_AXIS_BOTTOM;
    axisYType  = Y_AXIS_LEFT;

    // Set private members:
    myMinX = 0;
    myMaxX = 0;
    myMinY = 0;
    myMaxY = 0;
    myMaxTotal = 0;
    myIsMinMaxSet = false;
  }
  public boolean isMinMaxSet() {
    return myIsMinMaxSet;
  }
  public float maxTotal() { return myMaxTotal; }
  public float minX() { return myMinX; }
  public float maxX() { return myMaxX; }
  public float minY() { return myMinY; }
  public float maxY() { return myMaxY; }
  public void setMinMax( float minX, float maxX, float minY, float maxY ) {
    myIsMinMaxSet = true;
    myMinX = minX;
    myMaxX = maxX;
    myMinY = minY;
    myMaxY = maxY;
    myMaxTotal = Math.max(
        Math.max( Math.abs(myMaxX), Math.abs(myMaxY)),
        Math.max( Math.abs(myMinX), Math.abs(myMinY)) );
  }
  public void addListener( csCurveAttrChangeListener listener ) {
    myListeners.add( listener );
  }
  public void removeListener( csCurveAttrChangeListener listener ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      if( myListeners.get(i).equals(listener) ) {
        myListeners.remove( listener );
        return;
      }
    }
  }
  public void fireCurveAttrChangeEvent() {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).curveAttrChanged();
    }
  }
  public JPanel getPanel() {
    if( myPanel == null ) {
      myPanel = new CurveAttrPanel();
    }
    return myPanel;
  }
  /**
   * JPanel providing functionality to set curve attributes
   */
  @SuppressWarnings("serial")
protected class CurveAttrPanel extends JPanel {
    private JComboBox<String> myComboPointType;
    private JComboBox<String> myComboLineType;
    private JComboBox<String> myComboFillType;
//    private JComboBox<String> myComboAxisXType;
//    private JComboBox<String> myComboAxisYType;
    private JTextField myTextPointSize;
    private JTextField myTextLineSize;
    private csColorButton myButtonPointColor;
    private csColorButton myButtonLineColor;
    private csColorButton myButtonFillColor;

    CurveAttrPanel() {
      super( new GridBagLayout() );
      myComboPointType = new JComboBox<String>(TEXT_POINT_TYPES);
      myComboLineType  = new JComboBox<String>(TEXT_LINE_TYPES);
      myComboFillType  = new JComboBox<String>(TEXT_FILL_TYPES);
      myComboPointType.setSelectedIndex( pointType );
      myComboLineType.setSelectedIndex( lineType );
      myComboFillType.setSelectedIndex( filledType );
//      JComboBox myComboAxisXType;
//      JComboBox myComboAxisYType;
      myTextPointSize = new JTextField("" + pointSize);
      myTextLineSize = new JTextField("" + lineSize);
      myButtonPointColor = new csColorButton(this,pointColor);
      myButtonLineColor = new csColorButton(this,lineColor);
      myButtonFillColor = new csColorButton(this,fillColor);
//    JComboBox myComboAxisXType;
//    JComboBox myComboAxisYType;
      myTextPointSize = new JTextField("" + pointSize);
      myTextLineSize = new JTextField("" + lineSize);
      myButtonPointColor = new csColorButton(this,pointColor);
      myButtonLineColor = new csColorButton(this,lineColor);
      myButtonFillColor = new csColorButton(this,fillColor);      

      myTextPointSize.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          try {
            pointSize = Integer.parseInt( myTextPointSize.getText() );
            if( pointSize <= 0 ) throw( new NumberFormatException() );
            csCurveAttributes.this.fireCurveAttrChangeEvent();
          }
          catch( NumberFormatException exc ) {
            myTextPointSize.setText("" + pointSize );
          }
        }
      });
      myTextLineSize.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          try {
            lineSize = Integer.parseInt( myTextLineSize.getText() );
            if( lineSize <= 0 ) throw( new NumberFormatException() );
            csCurveAttributes.this.fireCurveAttrChangeEvent();
          }
          catch( NumberFormatException exc ) {
            myTextLineSize.setText("" + lineSize );
          }
        }
      });
      myTextPointSize.addFocusListener( new FocusListener() {
        @Override
        public void focusGained(FocusEvent e) {
        }
        @Override
        public void focusLost(FocusEvent e) {
          try {
            pointSize = Integer.parseInt( myTextPointSize.getText() );
            if( pointSize <= 0 ) throw( new NumberFormatException() );
            csCurveAttributes.this.fireCurveAttrChangeEvent();
          }
          catch( NumberFormatException exc ) {
            myTextPointSize.setText("" + pointSize );
          }
        }
      });
      myTextLineSize.addFocusListener( new FocusListener() {
        @Override
        public void focusGained(FocusEvent e) {
        }
        @Override
        public void focusLost(FocusEvent e) {
          try {
            lineSize = Integer.parseInt( myTextLineSize.getText() );
            if( lineSize <= 0 ) throw( new NumberFormatException() );
            csCurveAttributes.this.fireCurveAttrChangeEvent();
          }
          catch( NumberFormatException exc ) {
            myTextLineSize.setText("" + lineSize );
          }
        }
      });
      
      myButtonPointColor.addColorChangeListener( new csColorChangeListener() {
        @Override
        public void colorChanged( Object obj, Color color ) {
          pointColor = color;
          csCurveAttributes.this.fireCurveAttrChangeEvent();
        }
      });
      myButtonLineColor.addColorChangeListener( new csColorChangeListener() {
        @Override
        public void colorChanged( Object obj, Color color ) {
          lineColor = color;
          csCurveAttributes.this.fireCurveAttrChangeEvent();
        }
      });
      myButtonFillColor.addColorChangeListener( new csColorChangeListener() {
        @Override
        public void colorChanged( Object obj, Color color ) {
          fillColor = color;
          csCurveAttributes.this.fireCurveAttrChangeEvent();
        }
      });

      myComboPointType.addItemListener( new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent e) {
          if( e.getStateChange() == ItemEvent.SELECTED ) {
            pointType = myComboPointType.getSelectedIndex();
            csCurveAttributes.this.fireCurveAttrChangeEvent();
          }
        }
      });
      myComboLineType.addItemListener( new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent e) {
          if( e.getStateChange() == ItemEvent.SELECTED ) {
            lineType = myComboLineType.getSelectedIndex();
            csCurveAttributes.this.fireCurveAttrChangeEvent();
          }
        }
      });
      myComboFillType.addItemListener( new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent e) {
          if( e.getStateChange() == ItemEvent.SELECTED ) {
            filledType = myComboFillType.getSelectedIndex();
            csCurveAttributes.this.fireCurveAttrChangeEvent();
          }
        }
      });
      
      createPanel();
    }
    private void createPanel() {
      int yp = 0;
      this.add( Box.createHorizontalGlue(), new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( new JLabel("Size"), new GridBagConstraints(
        1, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( Box.createHorizontalGlue(), new GridBagConstraints(
        2, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( new JLabel("Type"), new GridBagConstraints(
        3, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

      this.add( new JLabel("Point:"), new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( myTextPointSize, new GridBagConstraints(
        1, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( myButtonPointColor, new GridBagConstraints(
        2, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( myComboPointType, new GridBagConstraints(
        3, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

      this.add( new JLabel("Line:"), new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( myTextLineSize, new GridBagConstraints(
        1, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( myButtonLineColor, new GridBagConstraints(
        2, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( myComboLineType, new GridBagConstraints(
        3, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

      this.add( new JLabel("Fill:"), new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( Box.createHorizontalGlue(), new GridBagConstraints(
        1, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( myButtonFillColor, new GridBagConstraints(
        2, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
      this.add( myComboFillType, new GridBagConstraints(
        3, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    }
  }
}



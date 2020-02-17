/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.graph;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.*;

import cseis.general.csStandard;
import java.util.Collection;
import java.util.Iterator;

/**
 * Dialog where graph settings can be changed.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csGraphDispDialog extends JDialog implements csIGraph2DListener {
  private static final int TEXT_MIN_WIDTH = 80;
  
  protected csGraphAttributes myAttr;
  
  // Text fields
  private final JTextField myTextInset;
  private final JTextField myTextMinAxisX;
  private final JTextField myTextMaxAxisX;
  private final JTextField myTextMinAxisY;
  private final JTextField myTextMaxAxisY;

  private final JTextField myBorderSize;
  private final JTextField myBorderPadding;
  private final JTextField myInnerBorderSize;
  private final JTextField myInnerBorderPadding;

  private final JTextField myTextGridIncX;
  private final JTextField myTextGridIncY;
  private final JTextField myTextMinorRatioX;
  private final JTextField myTextMinorRatioY;

  private final JTextField myGraphTitle;
  private final JTextField myXAxisLabel;
  private final JTextField myYAxisLabel;
  private final JTextField myTextRefValueY;

  private final JCheckBox myAutoScaleAxes;
  private final JCheckBox myCentreAxisX;
  private final JCheckBox myCentreAxisY;

  private final JCheckBox myShowBorder;
  private final JCheckBox myShowInnerBorder;
  private final JCheckBox myShowGrid;
  private final JCheckBox myShowAxisX;
  private final JCheckBox myShowAxisY;
  private final JCheckBox myShowAxisZero;
  private final JCheckBox myShowXTics;
  private final JCheckBox myShowYTics;
  private final JCheckBox myShowAnnotationX;
  private final JCheckBox myShowAnnotationY;

  private final JRadioButton myButtonYLinearScale;
  private final JRadioButton myButtonYLogScale;
  private final JRadioButton myButtonXLinearScale;
  private final JRadioButton myButtonXLogScale;
  private final JLabel myLabelXAxisScale;
  private final JLabel myLabelYAxisScale;

  private final JCheckBox myIsTrueScale;

//  private final JLabel myLabelMinAxisX;
//  private final JLabel myLabelMaxAxisX;
//  private final JLabel myLabelMinAxisY;
//  private final JLabel myLabelMaxAxisY;
//  private final JLabel myLabelGraphTitle;
  
  private final csGraph2D myGraph;
  private final JTabbedPane myPanelCurves;

  
  private JCheckBox myUseRefValueY;
  private boolean myIsUpdating;

  //----------------------
  private final JButton myButtonOK;
  private final JButton myButtonApply;
  private final JButton myButtonCancel;
  
  public csGraphDispDialog( csGraph2D graph ) {
    super.setTitle( "Graph settings" );
    super.setModal( false );
    
    myGraph = graph;
    myGraph.addGraph2DListener( this);

    myAttr = myGraph.getGraphAttributes();
    myIsUpdating = false;
    
    // Text fields
    myTextInset = new JTextField(" ");
    myTextGridIncX = new JTextField(" ");
    myTextMinAxisX = new JTextField(" ");
    myTextMaxAxisX = new JTextField(" ");
    myTextMinAxisY = new JTextField(" ");
    myTextMaxAxisY = new JTextField(" ");

    myTextGridIncY = new JTextField(" ");
    myTextMinorRatioX = new JTextField(" ");
    myTextMinorRatioY = new JTextField(" ");

    myInnerBorderSize = new JTextField(" ");
    myInnerBorderPadding = new JTextField(" ");
    myBorderSize = new JTextField(" ");
    myBorderPadding = new JTextField(" ");

    myGraphTitle = new JTextField(" ");
    myXAxisLabel = new JTextField(" ");
    myYAxisLabel = new JTextField(" ");
    myTextRefValueY = new JTextField(" ");

    int height = myTextInset.getPreferredSize().height;
    myTextInset.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextGridIncX.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );

    myTextMinAxisX.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMinAxisY.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMaxAxisX.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMaxAxisY.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextGridIncY.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMinorRatioY.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMinorRatioX.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myBorderSize.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myBorderPadding.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myInnerBorderSize.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myInnerBorderPadding.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myGraphTitle.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myXAxisLabel.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myYAxisLabel.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextRefValueY.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );

    myAutoScaleAxes   = new JCheckBox("Auto scale",myAttr.autoScaleAxes);
    myCentreAxisX     = new JCheckBox("Centre X axis",myAttr.centreAxisX);
    myCentreAxisY     = new JCheckBox("Centre Y axis",myAttr.centreAxisY);
    // Axes setup
    myShowAxisX       = new JCheckBox("",myAttr.showAxisX);
    myShowAxisY       = new JCheckBox("",myAttr.showAxisY);
    myShowXTics       = new JCheckBox("",myAttr.showXTics);
    myShowYTics       = new JCheckBox("",myAttr.showYTics);
    myShowAnnotationX = new JCheckBox("",myAttr.showAxisXAnnotation);
    myShowAnnotationY = new JCheckBox("",myAttr.showAxisYAnnotation);
    // General settings
    myShowAxisZero    = new JCheckBox("Zero axis",myAttr.showZeroAxis);
    myShowBorder      = new JCheckBox("Outer border",myAttr.showBorder);
    myShowInnerBorder = new JCheckBox("Inner border",myAttr.showInnerBorder);
    myShowGrid        = new JCheckBox("Grid",myAttr.showGrid);

    myButtonYLinearScale  = new JRadioButton("Linear",myAttr.axisScaleY == csGraphAttributes.AXIS_SCALE_LINEAR);
    myButtonYLogScale  = new JRadioButton("Logarithmic",myAttr.axisScaleY == csGraphAttributes.AXIS_SCALE_LOG);
    myButtonXLinearScale  = new JRadioButton("Linear",myAttr.axisScaleX == csGraphAttributes.AXIS_SCALE_LINEAR);
    myButtonXLogScale  = new JRadioButton("Logarithmic",myAttr.axisScaleY == csGraphAttributes.AXIS_SCALE_LOG);
    myLabelXAxisScale = new JLabel("X axis");
    myLabelYAxisScale = new JLabel("Y axis");
    ButtonGroup groupXAxis = new ButtonGroup();
    groupXAxis.add(myButtonXLinearScale);
    groupXAxis.add(myButtonXLogScale);
    ButtonGroup groupYAxis = new ButtonGroup();
    groupYAxis.add(myButtonYLinearScale);
    groupYAxis.add(myButtonYLogScale);

    myIsTrueScale = new JCheckBox("True scale",myAttr.isTrueScale);
    myUseRefValueY = new JCheckBox("Use ref value",myAttr.useRefValueY);

    // Labels
//    myLabelMinAxisX = new JLabel("Min X value");
//    myLabelMaxAxisX = new JLabel("Max X value");
//    myLabelMinAxisY = new JLabel("Min Y value");
//    myLabelMaxAxisY = new JLabel("Max Y value");
//    myLabelGraphTitle = new JLabel("Title");

    //=============================================================================
    // Creating panels
    //
    JPanel panelCurves = new JPanel(new GridBagLayout());
    panelCurves.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Curve attributes"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelOther = new JPanel(new GridBagLayout());
    panelOther.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Settings 3"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    //--------------------------------------------------------------
    int yp = 0;

/*
    panelText.add( myLabelMinAxisX, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelText.add( myTextMinAxisX, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelText.add( myLabelMaxAxisX, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelText.add( myTextMaxAxisX, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelText.add( myLabelMinAxisY, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelText.add( myTextMinAxisY, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelText.add( myLabelMaxAxisY, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelText.add( myTextMaxAxisY, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
*/

//----------------------------------------------------------------------
    JPanel panelAxes = new JPanel(new GridBagLayout());
    panelAxes.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Axes setup"),
        csStandard.INNER_EMPTY_BORDER ) );

    yp = 0;
    int xp = 0;
    panelAxes.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( new JLabel("Label"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.8, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( new JLabel("Axis"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( new JLabel("Tics"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelAxes.add( new JLabel("Increment"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelAxes.add( new JLabel("Min/maj ratio"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( new JLabel("Values"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    xp = 0;
    yp += 1;
    panelAxes.add( new JLabel("X:"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( myXAxisLabel, new GridBagConstraints(
        xp++, yp, 1, 1, 0.8, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( myShowAxisX, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( myShowXTics, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelAxes.add( myTextGridIncX, new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelAxes.add( myTextMinorRatioX, new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( myShowAnnotationX, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    xp = 0;
    yp += 1;
    panelAxes.add( new JLabel("Y:"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( myYAxisLabel, new GridBagConstraints(
        xp++, yp, 1, 1, 0.8, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( myShowAxisY, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( myShowYTics, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelAxes.add( myTextGridIncY, new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelAxes.add( myTextMinorRatioY, new GridBagConstraints(
        xp++, yp, 1, 1, 0.1, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 10 ), 0, 0 ) );
    panelAxes.add( myShowAnnotationY, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    //------------------------------------------------------------------------
    //
    JPanel panelGeneral = new JPanel(new GridBagLayout());
    panelGeneral.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("General settings"),
        csStandard.INNER_EMPTY_BORDER ) );
    xp = 0;
    yp = 0;
    panelGeneral.add( myShowGrid, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelGeneral.add( myShowAxisZero, new GridBagConstraints(
        xp++, yp++, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelGeneral.add( new JLabel("Title:"), new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelGeneral.add( myGraphTitle, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelGeneral.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.VERTICAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    //------------------------------------------------------------------------
    //
    JPanel panelBorder = new JPanel(new GridBagLayout());
    panelBorder.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Border setup"),
        csStandard.INNER_EMPTY_BORDER ) );
    xp = 0;
    yp = 0;
    panelBorder.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelBorder.add( new JLabel("Size"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelBorder.add( new JLabel("Padding"), new GridBagConstraints(
        xp++, yp++, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    xp = 0;
    panelBorder.add( myShowBorder, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelBorder.add( myBorderSize, new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelBorder.add( myBorderPadding, new GridBagConstraints(
        xp++, yp++, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    xp = 0;
    panelBorder.add( myShowInnerBorder, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelBorder.add( myInnerBorderSize, new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelBorder.add( myInnerBorderPadding, new GridBagConstraints(
        xp++, yp++, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    xp = 0;
    panelBorder.add( new JLabel("Inset:"), new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelBorder.add( myTextInset, new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelBorder.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, yp, 1, 1, 0.5, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    //------------------------------------------------------------------------
    //
    yp = 0;

    panelOther.add( myLabelXAxisScale, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelOther.add( myButtonXLinearScale, new GridBagConstraints(
        1, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelOther.add( myButtonXLogScale, new GridBagConstraints(
        2, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelOther.add( myLabelYAxisScale, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelOther.add( myButtonYLinearScale, new GridBagConstraints(
        1, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelOther.add( myButtonYLogScale, new GridBagConstraints(
        2, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

//    panelOther.add( myUseRefValueY, new GridBagConstraints(
//        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
//        GridBagConstraints.NONE, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
//    panelOther.add( myTextRefValueY, new GridBagConstraints(
//        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
//        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    //  ----------------------------
    //  Button panel
     myButtonOK = new JButton("OK");
     myButtonApply = new JButton("Apply");
     myButtonCancel = new JButton("Cancel");
     this.getRootPane().setDefaultButton(myButtonApply);

     //  ----------------------------
     myPanelCurves = new JTabbedPane( JTabbedPane.TOP, JTabbedPane.SCROLL_TAB_LAYOUT );
     updateCurvePanel();

     JPanel panelMany = new JPanel(new GridBagLayout());
     panelMany.setBorder( BorderFactory.createCompoundBorder(
         BorderFactory.createTitledBorder("General settings"),
         csStandard.INNER_EMPTY_BORDER ) );
     panelMany.add( panelAxes, new GridBagConstraints(
         0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.NORTH,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelMany.add( panelGeneral, new GridBagConstraints(
         0, 1, 1, 1, 1.0, 1.0, GridBagConstraints.NORTH,
         GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelMany.add( panelBorder, new GridBagConstraints(
         0, 2, 1, 1, 1.0, 1.0, GridBagConstraints.NORTH,
         GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     JPanel panelButtons = new JPanel( new GridBagLayout() );
     xp = 0;
     panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
         xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelButtons.add( myButtonOK, new GridBagConstraints(
         xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
         GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
     panelButtons.add( myButtonApply, new GridBagConstraints(
         xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
         GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
     panelButtons.add( myButtonCancel, new GridBagConstraints(
         xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
         GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
     panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
         xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     JPanel panelAll = new JPanel(new BorderLayout());
     panelAll.setBorder( csStandard.DIALOG_BORDER );

     yp = 0;

//     panelCurves.add( myPanelCurves, BorderLayout.WEST );
     panelCurves.add( myPanelCurves, new GridBagConstraints(
         0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 5, 5, 5, 5 ), 0, 0 ) );
//     panelCurves.add( Box.createHorizontalGlue(), new GridBagConstraints(
//         1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
//         GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     panelAll.add(panelCurves,BorderLayout.NORTH);
     panelAll.add(panelMany,BorderLayout.CENTER);
     panelAll.add(panelButtons,BorderLayout.SOUTH);
     this.getContentPane().add(panelAll);

     this.changedSettings( myAttr );

     // Other settings
     myTextRefValueY.setEnabled( myUseRefValueY.isSelected() );
     
     this.pack();

     //------------------------------------------------------------
     //
     myButtonOK.addActionListener(new ActionListener() {
       @Override
       public void actionPerformed(ActionEvent e) {
         apply();
         dispose();
       }
     });
     myButtonApply.addActionListener(new ActionListener() {
       @Override
       public void actionPerformed(ActionEvent e) {
         apply();
       }
     });
     myButtonCancel.addActionListener(new ActionListener() {
       @Override
       public void actionPerformed(ActionEvent e) {
         cancel();
       }
     });

     myUseRefValueY.addActionListener(new ActionListener() {
       @Override
       public void actionPerformed(ActionEvent e) {
         myTextRefValueY.setEnabled( myUseRefValueY.isSelected() );
         apply();
       }
     });

      myAutoScaleAxes.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myCentreAxisX.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myCentreAxisY.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowAxisX.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowAxisY.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowAxisZero.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowXTics.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowYTics.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowAnnotationX.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowAnnotationY.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowBorder.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowInnerBorder.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myShowGrid.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myIsTrueScale.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myButtonYLinearScale.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myButtonYLogScale.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myButtonXLinearScale.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });
      myButtonXLogScale.addActionListener(new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          apply();
        }
      });

     this.addWindowListener(new WindowAdapter() {
       @Override
       public void windowClosing(WindowEvent e) {
         cancel();
       }
     });
  }

  //*****************************************************************************************
  //
  //
  public void updateCurvePanel() {
    myPanelCurves.removeAll();
    Collection<csCurve> curveList = myGraph.getCurves();
    Iterator<csCurve> iter = curveList.iterator();
    int counter = 0;
    while( iter.hasNext() ) {
      csCurve curve = iter.next();
      myPanelCurves.add( curve.attr.getPanel() );
      myPanelCurves.setTitleAt( counter++, curve.attr.name );
    }
  }
  private void cancel() {
    dispose();
  }
  public void apply() {
    if( myGraph == null ) return;
    csGraphAttributes ga = new csGraphAttributes( myAttr );  // Copy current settings

    // Text fields
    try {
      String text;
      text = myTextInset.getText();
      ga.insetLeft   = Integer.parseInt( text );
      ga.insetBottom = ga.insetLeft;
      ga.insetTop    = ga.insetLeft;
      ga.insetRight  = ga.insetLeft;
      text = myTextMinAxisX.getText();
      ga.minAxisX = Float.parseFloat( text );
      text = myTextMaxAxisX.getText();
      ga.maxAxisX = Float.parseFloat( text );
      text = myTextMinAxisY.getText();
      ga.minAxisY = Float.parseFloat( text );
      text = myTextMaxAxisY.getText();
      ga.maxAxisY = Float.parseFloat( text );
      text = myTextGridIncX.getText();
      ga.gridXPixelInc = Integer.parseInt( text );
      text = myTextGridIncY.getText();
      ga.gridYPixelInc = Integer.parseInt( text );
      text = myTextMinorRatioX.getText();
      ga.minorXTicRatio = Integer.parseInt( text );
      if( ga.minorXTicRatio < 1 ) {
        ga.minorXTicRatio = 1;
        myTextMinorRatioX.setText("1");
      }
      text = myTextMinorRatioY.getText();
      ga.minorYTicRatio = Integer.parseInt( text );
      if( ga.minorYTicRatio < 1 ) {
        ga.minorYTicRatio = 1;
        myTextMinorRatioY.setText("1");
      }

      text = myBorderSize.getText();
      ga.borderSize = Integer.parseInt( text );
      text = myBorderPadding.getText();
      ga.borderPadding = Integer.parseInt( text );
      text = myInnerBorderSize.getText();
      ga.innerBorderSize = Integer.parseInt( text );
      text = myInnerBorderPadding.getText();
      ga.innerBorderPadding = Integer.parseInt( text );
      
      ga.title  = myGraphTitle.getText();
      ga.xLabel = myXAxisLabel.getText();
      ga.yLabel = myYAxisLabel.getText();

      text = myTextRefValueY.getText();
      ga.maxRefValueY = Float.parseFloat( text );
    }
    catch( NumberFormatException e ) {
//      System.out.println("Error in text fields...");
      return;
    }
    ga.autoScaleAxes = myAutoScaleAxes.isSelected();
    ga.showAxisX = myShowAxisX.isSelected();
    ga.showAxisY = myShowAxisY.isSelected();
    ga.showXTics = myShowXTics.isSelected();
    ga.showYTics = myShowYTics.isSelected();
    ga.showAxisXAnnotation = myShowAnnotationX.isSelected();
    ga.showAxisYAnnotation = myShowAnnotationY.isSelected();
    ga.showZeroAxis = myShowAxisZero.isSelected();
    ga.centreAxisX = myCentreAxisX.isSelected();
    ga.centreAxisY = myCentreAxisY.isSelected();
    ga.showBorder = myShowBorder.isSelected();
    ga.showInnerBorder = myShowInnerBorder.isSelected();
    ga.showGrid = myShowGrid.isSelected();
    ga.isTrueScale = myIsTrueScale.isSelected();
    ga.useRefValueY = myUseRefValueY.isSelected();
    ga.axisScaleX = myButtonXLinearScale.isSelected() ? csGraphAttributes.AXIS_SCALE_LINEAR : csGraphAttributes.AXIS_SCALE_LOG;
    ga.axisScaleY = myButtonYLinearScale.isSelected() ? csGraphAttributes.AXIS_SCALE_LINEAR : csGraphAttributes.AXIS_SCALE_LOG;

    myAttr = ga;
    if( myGraph != null ) {
      myIsUpdating = true;
      myGraph.setGraphAttributes( myAttr );
      myIsUpdating = false;
    }
  }

  //-------------------------------------------------------------------------------
  //
  //
  public void changedSettings( csGraphAttributes ga ) {
    if( myIsUpdating ) return;

    myAttr = ga;

    myTextInset.setText(""+ga.insetTop);
    myTextMinAxisX.setText(""+ga.minAxisX);
    myTextMaxAxisX.setText(""+ga.maxAxisX);
    myTextMinAxisY.setText(""+ga.minAxisY);
    myTextMaxAxisY.setText(""+ga.maxAxisY);
    myTextGridIncX.setText(""+ga.gridXPixelInc);
    myTextGridIncY.setText(""+ga.gridYPixelInc);
    myTextMinorRatioX.setText(""+ga.minorXTicRatio);
    myTextMinorRatioY.setText(""+ga.minorYTicRatio);

    myBorderSize.setText(""+ga.borderSize);
    myBorderPadding.setText(""+ga.borderPadding);
    myInnerBorderSize.setText(""+ga.innerBorderSize);
    myInnerBorderPadding.setText(""+ga.innerBorderPadding);

    myGraphTitle.setText(""+ga.title);
    myXAxisLabel.setText(""+ga.xLabel);
    myYAxisLabel.setText(""+ga.yLabel);
    myTextRefValueY.setText(""+ga.maxRefValueY);

    myAutoScaleAxes.setSelected(ga.autoScaleAxes);
    myShowAxisX.setSelected(ga.showAxisX);
    myShowAxisY.setSelected(ga.showAxisY);
    myShowXTics.setSelected(ga.showXTics);
    myShowYTics.setSelected(ga.showYTics);
    myShowAnnotationX.setSelected(ga.showAxisXAnnotation);
    myShowAnnotationY.setSelected(ga.showAxisYAnnotation);
    myShowAxisZero.setSelected(ga.showZeroAxis);
    myShowBorder.setSelected(ga.showBorder);
    myShowGrid.setSelected(ga.showGrid);
    myCentreAxisX.setSelected(ga.centreAxisX);
    myCentreAxisY.setSelected(ga.centreAxisY);
    myIsTrueScale.setSelected(ga.isTrueScale);
    myUseRefValueY.setSelected(ga.useRefValueY);
    myButtonXLinearScale.setSelected(ga.axisScaleX == csGraphAttributes.AXIS_SCALE_LINEAR);
    myButtonXLogScale.setSelected(ga.axisScaleX == csGraphAttributes.AXIS_SCALE_LOG);
    myButtonYLinearScale.setSelected(ga.axisScaleY == csGraphAttributes.AXIS_SCALE_LINEAR);
    myButtonYLogScale.setSelected(ga.axisScaleY == csGraphAttributes.AXIS_SCALE_LOG);
  }
  //--------------------------
  //
  //
  
  public static void main( String[] args ) {
    try {
//      UIManager.setLookAndFeel( "com.sun.java.swing.plaf.windows.WindowsLookAndFeel" );
//      UIManager.setLookAndFeel( new MetalLookAndFeel() );
    }
    catch( Exception e ) {
      e.printStackTrace();
    }
    
    csGraph2D graph = new csGraph2D();
    csGraphDispDialog dialog = new csGraphDispDialog( graph );
    dialog.setVisible(true);
    
    dialog.addWindowListener(new WindowAdapter() {
      public void windowClosing( WindowEvent e ) {
        System.exit( 0 );
      }
    });    
  }

  @Override
  public void graph2DValues(float xModel, float yModel) {
  }
  @Override
  public void graphChanged() {
    updateCurvePanel();
  }
}



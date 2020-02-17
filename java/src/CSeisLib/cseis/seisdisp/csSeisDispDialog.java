/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import cseis.cmapgen.csColorMapComboBox;
import cseis.cmapgen.csCustomColorMapModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.text.DecimalFormat;

import javax.swing.*;

import cseis.seis.csISeismicTraceBuffer;
import cseis.general.csColorMap;
import cseis.general.csStandard;
import cseis.general.csSquareIcon;

/**
 * Seismic display settings dialog.<br>
 * Dialog where settings for seismic display can be changed.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csSeisDispDialog extends JDialog implements csISeisViewListener {
  private static final int TEXT_MIN_WIDTH = 90;
  
  private csSeisDispSettings mySettings;

  // Text fields
  private JTextField myTextScalar;
  private JTextField myTextTraceScalar;
  private JTextField myTextMinValue;
  private JTextField myTextMaxValue;
  private JTextField myTextZoomVert;
  private JTextField myTextZoomHorz;
  private JTextField myTextZoomVertInch;
  private JTextField myTextZoomHorzInch;
  private JTextField myTextMinorInc;
  private JTextField myTextMajorInc;
  private JTextField myTextTraceClip;

  private JCheckBox myBoxShowHorzLines;
  private JCheckBox myBoxAutoHorzLines;
  private JCheckBox myBoxShowZeroLines;
  private JCheckBox myBoxTraceClip;
  private JCheckBox myBoxWiggle;
  private JCheckBox myBoxPosFill;
  private JCheckBox myBoxNegFill;
  private JCheckBox myBoxVI;
  private JButton   myButtonAutoRange;
  
  private JCheckBox myBoxUseVarColor;
  
  // Group
  private JRadioButton myButtonScalar;
  private JRadioButton myButtonRange;
  private JRadioButton myButtonTrace;

  // Group
  private JRadioButton myButtonWiggleLinear;
  private JRadioButton myButtonWiggleCubic;
  
  // Group
  private JRadioButton myButtonPolarityNormal;
  private JRadioButton myButtonPolarityReversed;

  // Group
  private JRadioButton myButtonVIDiscrete;
  private JRadioButton myButtonVIVertical;
  private JRadioButton myButtonVI2DSpline;

  // group
  private JRadioButton myButtonMaximum;
  private JRadioButton myButtonAverage;
  
  // Color/color map selectors
  private JComboBox<Object> myComboColorMapWiggle;
  private JComboBox<Object> myComboColorMapVI;
  private JButton myButtonColorWigglePos;
  private JButton myButtonColorWiggleNeg;
  private JButton myButtonHighlightColor;
  
  private csSeisView mySeisView;

  //----------------------
  private JButton myButtonOK;
  private JButton myButtonApply;
  private JButton myButtonCancel;

  private JComboBox<String> myComboPlotDir;
  private JLabel myLabelPlotDir;
  //----------------------
  private JLabel myLabelMinValue;
  private JLabel myLabelMaxValue;
  private JLabel myLabelMinorInc;
  private JLabel myLabelMajorInc;

  private JLabel myLabelZoomVert;
  private JLabel myLabelZoomHorz;

  private JLabel myLabelColorMap;
//  private JLabel myLabelHighlightColor;
  private JPanel myPanelAllDefault;
  
  private boolean myIsUpdating = false;
  private DecimalFormat myFormat0000 = new DecimalFormat("0.0000");
  
  public csSeisDispDialog( csSeisView seisview, int colorMapType ) {
    this( seisview, colorMapType, seisview.getDispSettings() );
  }
  public csSeisDispDialog( csSeisView seisview, int colorMapType, csSeisDispSettings settings ) {
    this( seisview, settings, new csCustomColorMapModel(colorMapType) );
  }
  public csSeisDispDialog( csSeisView seisview, csSeisDispSettings settings,
          csCustomColorMapModel cmapModel )
  {
    super( seisview.getParentFrame(), "Seismic Display Settings" );
    super.setModal( false );
    
    mySeisView = seisview;
    mySettings = settings;

    // Text fields
    myTextScalar      = new JTextField("1");
    myTextTraceScalar  = new JTextField("1");
    myTextMinValue    = new JTextField("0");
    myTextMaxValue    = new JTextField("1");
    myTextZoomVert    = new JTextField("1");
    myTextZoomHorz    = new JTextField("1");
    myTextZoomVertInch  = new JTextField("1");
    myTextZoomHorzInch  = new JTextField("1");
    myTextTraceClip   = new JTextField("1");
    myTextMinorInc    = new JTextField("1");
    myTextMajorInc    = new JTextField("2");
    //    myTextAGC         = new JTextField("500");
    int height = myTextScalar.getPreferredSize().height;
    myTextScalar.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextTraceScalar.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMinValue.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMaxValue.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextZoomVert.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextZoomHorz.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextZoomVertInch.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextZoomHorzInch.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextTraceClip.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMinorInc.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    myTextMajorInc.setPreferredSize( new Dimension( TEXT_MIN_WIDTH, height ) );
    
    myBoxWiggle        = new JCheckBox( "Show wiggle" );
    myBoxTraceClip     = new JCheckBox();
    myBoxShowZeroLines = new JCheckBox( "Show zero lines" );
    myBoxPosFill       = new JCheckBox( "Positive fill" );
    myBoxNegFill       = new JCheckBox( "Negative fill" );
    myBoxVI            = new JCheckBox( "Show VI plot" );
    myBoxShowHorzLines = new JCheckBox( "Show time lines" );
    myBoxAutoHorzLines = new JCheckBox( "Automatic settings" );
//    myBoxAGC           = new JCheckBox( "AGC" );

    myBoxUseVarColor   = new JCheckBox("Use variable fill color");

    // Group
    myButtonWiggleLinear = new JRadioButton("Linear");
    myButtonWiggleCubic  = new JRadioButton("Cubic");
    ButtonGroup group = new ButtonGroup();
    group.add(myButtonWiggleLinear);
    group.add(myButtonWiggleCubic);

    // group
    myButtonPolarityNormal   = new JRadioButton( "Normal" );
    myButtonPolarityReversed = new JRadioButton( "Reversed" );
    group = new ButtonGroup();
    group.add( myButtonPolarityNormal );
    group.add( myButtonPolarityReversed );

    // group
    myButtonScalar = new JRadioButton("Scalar");
    myButtonRange  = new JRadioButton("Range");
    myButtonTrace  = new JRadioButton("Full trace");
    group = new ButtonGroup();
    group.add( myButtonScalar );
    group.add( myButtonRange );
    group.add( myButtonTrace );

    // group
    myButtonVIDiscrete = new JRadioButton("Discrete");
    myButtonVIVertical = new JRadioButton("Vertical interpolation");
    myButtonVI2DSpline = new JRadioButton("2D spline interpolation");
    group = new ButtonGroup();
    group.add( myButtonVIDiscrete );
    group.add( myButtonVIVertical );
    group.add( myButtonVI2DSpline );
    
    myButtonMaximum = new JRadioButton("Maximum");
    myButtonAverage = new JRadioButton("Average");
    group = new ButtonGroup();
    group.add( myButtonAverage );
    group.add( myButtonMaximum );
    
    myButtonAutoRange = new JButton("Auto");
    
    // Color/color map selectors
    myButtonColorWigglePos = new JButton("Select");
    myButtonColorWiggleNeg = new JButton("Select");
    myButtonHighlightColor = new JButton("Select");

    myComboPlotDir = new JComboBox<String>();
    myComboPlotDir.addItem( "VERTICAL" );
    myComboPlotDir.addItem( "HORIZONTAL" );
//    myComboPlotDir.setEnabled(false);
    
    myComboColorMapVI = new csColorMapComboBox( cmapModel );
    myComboColorMapWiggle = new csColorMapComboBox( cmapModel );
    // Labels
    myLabelMinValue = new JLabel("Min value");
    myLabelMaxValue = new JLabel("Max value");
    myLabelMinorInc = new JLabel("Minor inc");
    myLabelMajorInc = new JLabel("Major inc");

    JLabel labelPixels     = new JLabel("pixel units");
    JLabel labelInch       = new JLabel("inch units");
    Font font = labelPixels.getFont();
    labelPixels.setFont( new Font( font.getName(), Font.PLAIN, font.getSize() ) );
    labelInch.setFont( new Font( font.getName(), Font.PLAIN, font.getSize() ) );
    myLabelZoomVert = new JLabel("Sample axis");
    myLabelZoomHorz = new JLabel("Trace axis");

    myLabelPlotDir  = new JLabel("Plot direction");    
    myLabelColorMap = new JLabel("Map");
    //=============================================================================
    // Putting all together
    //

    JPanel panelPreview = new JPanel(new GridBagLayout());
    panelPreview.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Preview"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelScaling = new JPanel(new GridBagLayout());
    panelScaling.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Scaling"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelWiggle = new JPanel(new GridBagLayout());
    panelWiggle.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Wiggle plot"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelVADisplay = new JPanel(new GridBagLayout());
    panelVADisplay.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Variable intensity plot"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelInterpolation = new JPanel(new GridBagLayout());
    panelInterpolation.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Vertical interpolation"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelPolarity = new JPanel(new GridBagLayout());
    panelPolarity.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Polarity"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelDisplay = new JPanel(new GridBagLayout());
    panelDisplay.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Display"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelZoom = new JPanel(new GridLayout(3,3));
    panelZoom.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Zoom settings"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JPanel panelTimeAxis = new JPanel(new GridBagLayout());
    panelTimeAxis.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Time axis annotation"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelClipping = new JPanel(new GridBagLayout());
    panelClipping.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Trace clipping"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelPlotDir = new JPanel(new GridBagLayout());
    panelPlotDir.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Plot direction"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    //--------------------------------------------------------------
    int yp = 0;
    panelScaling.add( myButtonScalar, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelScaling.add( myTextScalar, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelScaling.add( myButtonRange, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 5, 0, 0, 0 ), 0, 0 ) );
    panelScaling.add( myButtonAutoRange, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 5, 0, 0, 0 ), 0, 0 ) );
    panelScaling.add( myLabelMinValue, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 15, 0, 5 ), 0, 0 ) );
    panelScaling.add( myTextMinValue, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelScaling.add( myLabelMaxValue, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 15, 0, 5 ), 0, 0 ) );
    panelScaling.add( myTextMaxValue, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelScaling.add( myButtonTrace, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 5, 0, 0, 5 ), 0, 0 ) );
    panelScaling.add( myTextTraceScalar, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 0, 0, 0 ), 0, 0 ) );
    panelScaling.add( myButtonAverage, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 15, 0, 0 ), 0, 0 ) );
    panelScaling.add( myButtonMaximum, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 15, 0, 0 ), 0, 0 ) );

    panelScaling.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    //-------------------------------
    
    yp = 0;
    panelWiggle.add( myBoxWiggle, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelWiggle.add( myBoxPosFill, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelWiggle.add( myButtonColorWigglePos, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelWiggle.add( myBoxNegFill, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelWiggle.add( myButtonColorWiggleNeg, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    panelWiggle.add( myBoxUseVarColor, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelWiggle.add( myComboColorMapWiggle, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 11, 0, 0 ), 0, 0 ) );
    panelWiggle.add( myBoxShowZeroLines, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 0, 0, 0 ), 0, 0 ) );

    panelWiggle.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    
    yp = 0;
    panelInterpolation.add( myButtonWiggleLinear, new GridBagConstraints(
        0, yp, 1, 1, 0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelInterpolation.add( myButtonWiggleCubic, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelInterpolation.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    yp = 0;
    panelPolarity.add( myButtonPolarityNormal, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelPolarity.add( myButtonPolarityReversed, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );    
    panelPolarity.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    yp = 0;
    panelZoom.add( Box.createHorizontalGlue() );
    panelZoom.add( labelPixels );
    panelZoom.add( labelInch );
    panelZoom.add( myLabelZoomVert );
    panelZoom.add( myTextZoomVert );
    panelZoom.add( myTextZoomVertInch );
    panelZoom.add( myLabelZoomHorz );
    panelZoom.add( myTextZoomHorz );
    panelZoom.add( myTextZoomHorzInch );
    
    yp = 0;
    panelVADisplay.add( myBoxVI, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelVADisplay.add( myButtonVIDiscrete, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 11, 0, 0 ), 0, 0 ) );
    panelVADisplay.add( myButtonVIVertical, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 11, 0, 0 ), 0, 0 ) );
    panelVADisplay.add( myButtonVI2DSpline, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 11, 0, 0 ), 0, 0 ) );
    panelVADisplay.add( myLabelColorMap, new GridBagConstraints(
        0, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 0, 0, 0 ), 0, 0 ) );
    panelVADisplay.add( myComboColorMapVI, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 5, 0, 0 ), 0, 0 ) );
    panelVADisplay.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    yp = 0;
    panelTimeAxis.add( myBoxAutoHorzLines, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTimeAxis.add( myLabelMinorInc, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 15, 0, 5 ), 0, 0 ) );
    panelTimeAxis.add( myTextMinorInc, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTimeAxis.add( myLabelMajorInc, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 15, 0, 5 ), 0, 0 ) );
    panelTimeAxis.add( myTextMajorInc, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelTimeAxis.add( myBoxShowHorzLines, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    
  /*  
    panelOther.add( myBoxShowHorzLines, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelOther.add( myBoxAutoHorzLines, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 15, 0, 0 ), 0, 0 ) );
    panelOther.add( myLabelMinorInc, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 15, 0, 5 ), 0, 0 ) );
    panelOther.add( myTextMinorInc, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelOther.add( myLabelMajorInc, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 15, 0, 5 ), 0, 0 ) );
    panelOther.add( myTextMajorInc, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
*/
    yp = 0;
    panelClipping.add( myBoxTraceClip, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelClipping.add( myTextTraceClip, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelClipping.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    yp = 0;
//    panelPlotDir.add( myLabelPlotDir, new GridBagConstraints(
//        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
 //       GridBagConstraints.NONE, new Insets( 5, 15, 0, 5 ), 0, 0 ) );
    panelPlotDir.add( myComboPlotDir, new GridBagConstraints(
        0, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelPlotDir.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 1, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    yp = 0;
/*    panelColors.add( myLabelColorWigglePos, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelColors.add( myButtonColorWigglePos, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelColors.add( myLabelColorWiggleNeg, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelColors.add( myButtonColorWiggleNeg, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelColors.add( myLabelHighlightColor, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelColors.add( myButtonHighlightColor, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
*/
/*    panelColors.add( myButtonColorWigglePos, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelColors.add( myButtonColorWigglePos, new GridBagConstraints(
        0, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelColors.add( myButtonHighlightColor, new GridBagConstraints(
        1, yp++, 2, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) ); */
/*    panelColors.add( myLabelColorMap, new GridBagConstraints(
        0, yp, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.NONE, new Insets( 0, 0, 0, 5 ), 0, 0 ) );
    panelColors.add( myComboColorMap, new GridBagConstraints(
        1, yp++, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelColors.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, yp, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
*/
    //  ----------------------------
    //  Button panel
     myButtonOK = new JButton("OK");
     myButtonApply = new JButton("Apply");
     myButtonCancel = new JButton("Cancel");

     this.getRootPane().setDefaultButton(myButtonApply);

    // remove the binding for pressed
    getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
            .put(KeyStroke.getKeyStroke("ENTER"), "none");
    // retarget the binding for released
    getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
            .put(KeyStroke.getKeyStroke("released ENTER"), "press");

//    JRootPane buttonRootPane = SwingUtilities.getRootPane( myButtonApply ); 
//    buttonRootPane.setDefaultButton( myButtonApply );
     
     //  ----------------------------
     int xp = 0;
     JPanel panelUpper = new JPanel(new GridBagLayout());
/*
     panelUpper.add( panelWiggle, new GridBagConstraints(
         0, 0, 1, 1, 0.5, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelUpper.add( panelScaling, new GridBagConstraints(
         0, 1, 1, 2, 0.5, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     panelUpper.add( panelVADisplay, new GridBagConstraints(
         1, 0, 1, 1, 0.5, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelUpper.add( panelZoom, new GridBagConstraints(
         1, 1, 1, 1, 0.5, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelUpper.add( panelTimeAxis, new GridBagConstraints(
         1, 2, 1, 1, 0.5, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     
     panelUpper.add( panelInterpolation, new GridBagConstraints(
         0, 3, 1, 1, 0.5, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelUpper.add( panelPolarity, new GridBagConstraints(
         0, 4, 1, 1, 0.5, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );


     panelUpper.add( panelOther, new GridBagConstraints(
         1, 3, 1, 2, 0.5, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
*/
     
     panelUpper.add( panelWiggle, new GridBagConstraints(
         0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelUpper.add( panelScaling, new GridBagConstraints(
         1, 0, 2, 1, 1.0, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     panelUpper.add( panelVADisplay, new GridBagConstraints(
         0, 1, 1, 2, 0.0, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     panelUpper.add( panelZoom, new GridBagConstraints(
         1, 1, 2, 1, 1.0, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelUpper.add( panelClipping, new GridBagConstraints(
         1, 2, 1, 1, 0.5, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelUpper.add( panelPlotDir, new GridBagConstraints(
         2, 2, 1, 1, 0.0, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.VERTICAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     panelUpper.add( panelInterpolation, new GridBagConstraints(
         0, 3, 1, 1, 0.0, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelUpper.add( panelPolarity, new GridBagConstraints(
         0, 4, 1, 1, 0.0, 0.5, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     
     panelUpper.add( panelTimeAxis, new GridBagConstraints(
         1, 3, 2, 2, 1.0, 1.0, GridBagConstraints.WEST,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

     
     JPanel panelButtons = new JPanel( new GridBagLayout() );
     xp = 0;
     panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
         xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
     panelButtons.add( myButtonOK, new GridBagConstraints(
         xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
         GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
     panelButtons.add( myButtonApply, new GridBagConstraints(
         xp++, 0, 1, 1, 0.1, 0.1, GridBagConstraints.SOUTHWEST,
         GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
     panelButtons.add( myButtonCancel, new GridBagConstraints(
         xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
         GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
     panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
         xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
         GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

//  Nuts and bolts
     myPanelAllDefault = new JPanel(new BorderLayout());
     myPanelAllDefault.setBorder( csStandard.DIALOG_BORDER );

     yp = 0;

     myPanelAllDefault.add(panelUpper,BorderLayout.CENTER);
     myPanelAllDefault.add(panelButtons,BorderLayout.SOUTH);
     this.getContentPane().add(myPanelAllDefault,BorderLayout.CENTER);

     this.changedSettings( mySettings );

     this.pack();
     this.setLocationRelativeTo( seisview );

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

     //-------------------------------------------------
     myBoxUseVarColor.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_Wiggle();
       }
     });
     myBoxNegFill.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_Wiggle();
       }
     });
     myBoxPosFill.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_Wiggle();
       }
     });

     myBoxVI.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_VA();
       }
     });
     myButtonVIDiscrete.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_VA();
       }
     });
     myButtonVIVertical.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_VA();
       }
     });
     myButtonVI2DSpline.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_VA();
       }
     });
     myBoxShowHorzLines.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_HorzLine();
       }
     });
     myBoxAutoHorzLines.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_HorzLine();
       }
     });
     myBoxTraceClip.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_TraceClip();
       }
     });
     myButtonScalar.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_Scaling();
       }
     });
     myButtonRange.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         changedSettings_Scaling();
       }
     });
     myButtonTrace.addActionListener(new ActionListener() {
       public void actionPerformed(ActionEvent e) {
         changedSettings_Scaling();
       }
     });
     myButtonAutoRange.addActionListener(new ActionListener() {
       public void actionPerformed(ActionEvent e) {
         csISeismicTraceBuffer buffer = mySeisView.getTraceBuffer();
         myTextMinValue.setText( "" + buffer.minTotalAmplitude() );
         myTextMaxValue.setText( "" + buffer.maxTotalAmplitude() );
       }
     });
     myButtonColorWigglePos.addActionListener(new ActionListener() {
       public void actionPerformed(ActionEvent e) {
         JButton b = (JButton)e.getSource();
         Color color = b.getBackground();
         Color colorNew = JColorChooser.showDialog( csSeisDispDialog.this, "Select color", color );
         if( colorNew != null ) {
           ((csSquareIcon)b.getIcon()).setColor(colorNew);
           b.repaint();
         }
       }
     });
     myButtonColorWiggleNeg.addActionListener(new ActionListener() {
       public void actionPerformed(ActionEvent e) {
         JButton b = (JButton)e.getSource();
         Color color = b.getBackground();
         Color colorNew = JColorChooser.showDialog( csSeisDispDialog.this, "Select color", color );
         if( colorNew != null ) {
           ((csSquareIcon)b.getIcon()).setColor(colorNew);
           b.repaint();
         }
       }
     });
     myButtonHighlightColor.addActionListener(new ActionListener() {
       public void actionPerformed(ActionEvent e) {
         JButton b = (JButton)e.getSource();
         Color color = b.getBackground();
         Color colorNew = JColorChooser.showDialog( csSeisDispDialog.this, "Select color", color );
         if( colorNew != null ) {
           ((csSquareIcon)b.getIcon()).setColor(colorNew);
           b.repaint();
         }
       }
     });
    myTextZoomVert.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        try {
          float zoomVert = Float.parseFloat( myTextZoomVert.getText() );
          myTextZoomVertInch.setText( myFormat0000.format( convertZoomVertToInch(zoomVert)) );
        } catch( NumberFormatException exc ) {}
      }
    });
     myTextZoomHorz.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         try {
          float zoomHorz = Float.parseFloat( myTextZoomHorz.getText() );
          myTextZoomHorzInch.setText( myFormat0000.format( convertZoomHorzToInch(zoomHorz)) );
         } catch( NumberFormatException exc ) {}
       }
     });
     myTextZoomVertInch.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         try {
          float zoomVertInch = Float.parseFloat( myTextZoomVertInch.getText() );
          myTextZoomVert.setText( myFormat0000.format( convertZoomVertToPixels(zoomVertInch)) );
         } catch( NumberFormatException exc ) {}
       }
     });
     myTextZoomHorzInch.addActionListener(new ActionListener() {
      @Override
       public void actionPerformed(ActionEvent e) {
         try {
          float zoomHorzInch = Float.parseFloat( myTextZoomHorzInch.getText() );
          myTextZoomHorz.setText( myFormat0000.format( convertZoomHorzToPixels(zoomHorzInch)) );
         } catch( NumberFormatException exc ) {}
       }
     });
     this.addWindowListener(new WindowAdapter() {
       public void windowClosing(WindowEvent e) {
         cancel();
       }
       public void windowActivated(WindowEvent evt) {
//         myBoxWiggle.requestFocusInWindow();
       }
     });
//     myBoxWiggle.requestFocusInWindow();
     setToolTips();
  }

  //*****************************************************************************************
  //
  //
  private void cancel() {
    dispose();
  }
  public void setSeisView( csSeisView seisview ) {
    mySeisView = seisview;
    // Add call to update function
  }
  public void apply() {
    if( mySeisView == null ) return;
    csSeisDispSettings ds = new csSeisDispSettings( mySettings );  // Copy current settings

    // Text fields
    try {
      String text;
      text = myTextScalar.getText();
      ds.dispScalar = Float.parseFloat( text );
      text = myTextZoomVert.getText();
      ds.zoomVert = Float.parseFloat( text );
      text = myTextZoomHorz.getText();
      ds.zoomHorz = Float.parseFloat( text );

      text = myTextTraceScalar.getText();
      ds.fullTraceScalar = Float.parseFloat( text );
      
      text = myTextTraceClip.getText();
      ds.traceClip = Float.parseFloat( text );

      text = myTextMinorInc.getText();
      ds.timeLineMinorInc = Float.parseFloat( text );
      text = myTextMajorInc.getText();
      ds.timeLineMajorInc = Float.parseFloat( text );

      text = myTextMinValue.getText();
      ds.minValue = Float.parseFloat( text );
      text = myTextMaxValue.getText();
      ds.maxValue = Float.parseFloat( text );
    }
    catch( NumberFormatException e ) {
      System.out.println("Error in text fields...");
      return;
    }

    // Group
    if( myButtonScalar.isSelected() ) {
      ds.scaleType = csSeisDispSettings.SCALE_TYPE_SCALAR;
    }
    else if( myButtonRange.isSelected() ) {
      ds.scaleType = csSeisDispSettings.SCALE_TYPE_RANGE;
    }
    else { //if( myButtonTrace.isSelected() ) {
      ds.scaleType = csSeisDispSettings.SCALE_TYPE_TRACE;
    }

    if( myButtonPolarityNormal.isSelected() ) {
      ds.polarity = csSeisDispSettings.POLARITY_NORMAL;
    }
    else {//if( myButtonPolarityreversed.isSelected() ) {
      ds.polarity = csSeisDispSettings.POLARITY_REVERSED;
    }

    if( myButtonVIDiscrete.isSelected() ) {
      ds.viType = csSeisDispSettings.VA_TYPE_DISCRETE;
    }
    else if( myButtonVIVertical.isSelected() ) {
      ds.viType = csSeisDispSettings.VA_TYPE_VERTICAL;
    }
    else { //if( myButtonVADiscrete.isSelected() ) {
      ds.viType = csSeisDispSettings.VA_TYPE_2DSPLINE;
    }

    ds.isVariableColor = myBoxUseVarColor.isSelected();

    if( myButtonMaximum.isSelected() ) {
      ds.traceScaling = csSeisDispSettings.TRACE_SCALING_MAXIMUM;
    }
    else { //if( myButtonAverage.isSelected() ) {
      ds.traceScaling = csSeisDispSettings.TRACE_SCALING_AVERAGE;
    }

    ds.showTimeLines = myBoxShowHorzLines.isSelected();
    ds.isTimeLinesAuto     = myBoxAutoHorzLines.isSelected();
    ds.showZeroLines       = myBoxShowZeroLines.isSelected();
    ds.doTraceClipping     = myBoxTraceClip.isSelected();
    ds.showWiggle          = myBoxWiggle.isSelected();
    ds.isVIDisplay         = myBoxVI.isSelected();
    ds.isPosFill           = myBoxPosFill.isSelected();
    ds.isNegFill           = myBoxNegFill.isSelected();

    if( myButtonWiggleCubic.isSelected() ) {
      ds.wiggleType = csSeisDispSettings.WIGGLE_TYPE_CUBIC;
    }
    else if( myButtonWiggleLinear.isSelected() ) {
      ds.wiggleType = csSeisDispSettings.WIGGLE_TYPE_LINEAR;
    }
        
    ds.plotDirection = myComboPlotDir.getSelectedIndex();
    // Color/color map selectors
    ds.wiggleColorPos = ((csSquareIcon)myButtonColorWigglePos.getIcon()).getColor();
    ds.wiggleColorNeg = ((csSquareIcon)myButtonColorWiggleNeg.getIcon()).getColor();
    ds.highlightColor = ((csSquareIcon)myButtonHighlightColor.getIcon()).getColor();
    ds.viColorMap.resetColors( ((csColorMapListItem)myComboColorMapVI.getSelectedItem()).map );
    ds.wiggleColorMap.resetColors( ((csColorMapListItem)myComboColorMapWiggle.getSelectedItem()).map );

    // Check integrity
    if( ds.timeLineMinorInc > ds.timeLineMajorInc ) {
      JOptionPane.showMessageDialog( this,
          "Wrong horizontal minor/major increment specified:\n\n" +
          "Minor increment for horizontal lines must be smaller than\n" +
          "or equal to major increment.",
          null, JOptionPane.ERROR_MESSAGE);
      return;
    }
    if( ds.timeLineMinorInc < mySeisView.getSampleInt() ) {
      JOptionPane.showMessageDialog( this,
          "Too small horizontal minor increment specified:\n\n" +
          "Minor increment for horizontal lines must be larger than\n" +
          "or equal to sample interval (=" + mySeisView.getSampleInt() + ").",
          null, JOptionPane.ERROR_MESSAGE);
      return;
    }

    mySettings = ds;
    if( mySeisView != null ) {
      myIsUpdating = true;
      mySeisView.updateDispSettings( mySettings );
      myIsUpdating = false;
    }
    setZoomInchUnits( mySettings.zoomVert, mySettings.zoomHorz );    
  }

  public void vertScrollChanged( int scrollValue ) {}
  public void horzScrollChanged( int scrollValue ) {}
  public void sizeChanged( Dimension size ) {}
  public void addCustomPanel( JPanel panelCustom ) {
    myPanelAllDefault.add( panelCustom, BorderLayout.NORTH );
    pack();
    invalidate();
    repaint();
  }
  //-------------------------------------------------------------------------------
  //
  //
  private double convertZoomVertToInch( double zoomVert ) {    
    double inchesPerPixel = 1.0 / (double)Toolkit.getDefaultToolkit().getScreenResolution();
    double inchesPerSecond = inchesPerPixel  * zoomVert * 1000.0 / mySeisView.getSampleInt();
    return inchesPerSecond;
  }
  private double convertZoomHorzToInch( double zoomHorz ) {    
    double inchesPerPixel = 1.0 / (double)Toolkit.getDefaultToolkit().getScreenResolution();
    double tracesPerInch   = 1.0 / ( inchesPerPixel  * zoomHorz );
    return tracesPerInch;
  }
  private double convertZoomVertToPixels( double zoomVertInch ) {    
    double inchesPerPixel = 1.0 / (double)Toolkit.getDefaultToolkit().getScreenResolution();
    double zoomVertInPixels = mySeisView.getSampleInt() * zoomVertInch / (1000.0 * inchesPerPixel);
    return zoomVertInPixels;
  }
  private double convertZoomHorzToPixels( double zoomHorzInch ) {    
    double inchesPerPixel = 1.0 / (double)Toolkit.getDefaultToolkit().getScreenResolution();
    double zoomHorzInPixels =  1.0 / ( inchesPerPixel  * zoomHorzInch );
    return zoomHorzInPixels;
  }
  private void setZoomInchUnits( double zoomVert, double zoomHorz ) {    
    myTextZoomVertInch.setText( myFormat0000.format( convertZoomVertToInch(zoomVert) ) );
    myTextZoomHorzInch.setText( myFormat0000.format( convertZoomHorzToInch(zoomHorz) ) );
  }
  @Override
  public void traceBufferChanged( csISeismicTraceBuffer traceBuffer ) {
    // Nothing to be done
  }
  @Override
  public void changedSettings( csSeisDispSettings ds ) {
    if( myIsUpdating ) return;

    mySettings = ds;
    mySettings.dispScalar = ds.dispScalar;
    mySettings.zoomVert = ds.zoomVert;
    
    myTextScalar.setText("" + mySettings.dispScalar );
    myTextTraceScalar.setText("" + mySettings.fullTraceScalar );
    myTextMinValue.setText("" + mySettings.minValue );
    myTextMaxValue.setText("" + mySettings.maxValue );

    myTextZoomVert.setText("" + myFormat0000.format(mySettings.zoomVert) );
    myTextZoomHorz.setText("" + myFormat0000.format(mySettings.zoomHorz) );
    myTextZoomVertInch.setText( myFormat0000.format( convertZoomVertToInch(mySettings.zoomVert)) );
    myTextZoomHorzInch.setText( myFormat0000.format( convertZoomHorzToInch(mySettings.zoomHorz)) );

    myTextTraceClip.setText("" + mySettings.traceClip );
    myTextMinorInc.setText("" + mySettings.timeLineMinorInc );
    myTextMajorInc.setText("" + mySettings.timeLineMajorInc );

    myBoxShowHorzLines.setSelected( mySettings.showTimeLines );
    myBoxAutoHorzLines.setSelected( mySettings.isTimeLinesAuto );
    myBoxShowZeroLines.setSelected( mySettings.showZeroLines );
    myBoxWiggle.setSelected(mySettings.showWiggle);
    myBoxPosFill.setSelected( mySettings.isPosFill );
    myBoxNegFill.setSelected( mySettings.isNegFill );
    myBoxVI.setSelected(mySettings.isVIDisplay);
    myBoxTraceClip.setSelected(mySettings.doTraceClipping);

    myButtonWiggleCubic.setSelected( mySettings.wiggleType == csSeisDispSettings.WIGGLE_TYPE_CUBIC );
    myButtonWiggleLinear.setSelected( mySettings.wiggleType == csSeisDispSettings.WIGGLE_TYPE_LINEAR );
    myBoxUseVarColor.setSelected( mySettings.isVariableColor );

//    myButtonFixedColor.setSelected( !mySettings.isVariableColor );
//    myButtonVariableColor.setSelected( mySettings.isVariableColor );

    myButtonPolarityNormal.setSelected( mySettings.polarity == csSeisDispSettings.POLARITY_NORMAL );
    myButtonPolarityReversed.setSelected( mySettings.polarity == csSeisDispSettings.POLARITY_REVERSED );

    myButtonScalar.setSelected( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_SCALAR );
    myButtonRange.setSelected( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_RANGE );
    myButtonTrace.setSelected( mySettings.scaleType == csSeisDispSettings.SCALE_TYPE_TRACE );

    // group
    myButtonVIDiscrete.setSelected( mySettings.viType == csSeisDispSettings.VA_TYPE_DISCRETE );
    myButtonVIVertical.setSelected( mySettings.viType == csSeisDispSettings.VA_TYPE_VERTICAL );
    myButtonVI2DSpline.setSelected( mySettings.viType == csSeisDispSettings.VA_TYPE_2DSPLINE );
    
    myButtonMaximum.setSelected( mySettings.traceScaling == csSeisDispSettings.TRACE_SCALING_MAXIMUM );
    myButtonAverage.setSelected( mySettings.traceScaling == csSeisDispSettings.TRACE_SCALING_AVERAGE );

    myComboPlotDir.setSelectedIndex(mySettings.plotDirection);

    myButtonColorWigglePos.setIcon( new csSquareIcon(csStandard.ICON_SIZE, mySettings.wiggleColorPos ) );
    myButtonColorWiggleNeg.setIcon( new csSquareIcon(csStandard.ICON_SIZE, mySettings.wiggleColorNeg ) );
    myButtonHighlightColor.setIcon( new csSquareIcon(csStandard.ICON_SIZE, mySettings.highlightColor ) );
    for( int i = 0; i < myComboColorMapVI.getItemCount(); i++ ) {
      if( mySettings.viColorMap.toString().compareTo( ((csColorMapListItem)myComboColorMapVI.getItemAt(i)).map.toString()) == 0 ) {
        myComboColorMapVI.setSelectedIndex(i);
        break;
      }
    }
    if( myComboColorMapVI.getSelectedIndex() < 0 ) myComboColorMapVI.setSelectedIndex(0);
    for( int i = 0; i < myComboColorMapWiggle.getItemCount(); i++ ) {
      if( mySettings.wiggleColorMap.toString().compareTo( ((csColorMapListItem)myComboColorMapWiggle.getItemAt(i)).map.toString()) == 0 ) {
        myComboColorMapWiggle.setSelectedIndex(i);
        break;
      }
    }
    if( myComboColorMapWiggle.getSelectedIndex() < 0 ) myComboColorMapWiggle.setSelectedIndex(0);
//    myComboColorMapVA.setSelectedIndex( mySettings.vaColorMap.getDefaultMapIndex() );
//    myComboColorMapWiggle.setSelectedIndex( mySettings.wiggleColorMap.getDefaultMapIndex() );

    setZoomInchUnits( mySettings.zoomVert, mySettings.zoomHorz );    

    changedSettings_VA();
    changedSettings_Wiggle();
    changedSettings_HorzLine();
    changedSettings_TraceClip();
    changedSettings_Scaling();
  }
  private void changedSettings_HorzLine() {
//    myBoxAutoHorzLines.setEnabled( myBoxShowHorzLines.isSelected() );
    boolean isSelected = !myBoxAutoHorzLines.isSelected();
    myLabelMinorInc.setEnabled( isSelected );
    myLabelMajorInc.setEnabled( isSelected );
    myTextMinorInc.setEnabled( isSelected );
    myTextMajorInc.setEnabled( isSelected );
  }
  private void changedSettings_Wiggle() {
    boolean isSelectedVA  = myBoxUseVarColor.isSelected();
    boolean isSelectedPos = myBoxPosFill.isSelected();
    boolean isSelectedNeg = myBoxNegFill.isSelected();
    myBoxUseVarColor.setEnabled(isSelectedNeg||isSelectedPos);
    myComboColorMapWiggle.setEnabled(isSelectedVA&&(isSelectedNeg||isSelectedPos));
    if( !isSelectedVA ) {
      myButtonColorWigglePos.setEnabled( isSelectedPos );
      myButtonColorWiggleNeg.setEnabled( isSelectedNeg );
    }
    else {
      myButtonColorWigglePos.setEnabled( !isSelectedVA );
      myButtonColorWiggleNeg.setEnabled( !isSelectedVA );
    }
  }
  private void changedSettings_VA() {
    boolean isSelected = myBoxVI.isSelected();
    myButtonVIDiscrete.setEnabled( isSelected );
    myButtonVIVertical.setEnabled( isSelected );
    myButtonVI2DSpline.setEnabled( isSelected && !myTextTraceScalar.isEnabled() );
    myLabelColorMap.setEnabled( isSelected );
    myComboColorMapVI.setEnabled( isSelected );
    myButtonTrace.setEnabled( !myButtonVI2DSpline.isSelected() );
  }
  private void changedSettings_TraceClip() {
    myTextTraceClip.setEnabled( myBoxTraceClip.isSelected() );
  }
  private void changedSettings_Scaling() {
    boolean selectedScalar = myButtonScalar.isSelected();
    boolean selectedRange  = myButtonRange.isSelected();
    boolean selectedTrace  = myButtonTrace.isSelected();
    myTextScalar.setEnabled(selectedScalar);
    myTextTraceScalar.setEnabled(selectedTrace);
    myButtonAutoRange.setEnabled(selectedRange);
    myLabelMinValue.setEnabled(selectedRange);
    myLabelMaxValue.setEnabled(selectedRange);
    myTextMinValue.setEnabled(selectedRange);
    myTextMaxValue.setEnabled(selectedRange);
    myButtonMaximum.setEnabled(selectedTrace);
    myButtonAverage.setEnabled(selectedTrace);

    if( selectedTrace ) {
      myButtonVI2DSpline.setEnabled( false );
    }
    else if( myBoxVI.isSelected() ) {
      myButtonVI2DSpline.setEnabled( true );
    }
  }
  //--------------------------
  //
  //
  private void setToolTips() {
    myTextScalar.setToolTipText("Constant scalar to apply to seismic amplitudes (for display only)");
    myTextTraceScalar.setToolTipText("Constant scalar to apply after full trace equalisation");
    myTextMinValue.setToolTipText("Minimum value for range setting");
    myTextMaxValue.setToolTipText("maximum value for range setting");
    myTextZoomVert.setToolTipText("Screen pixels per sample");
    myTextZoomHorz.setToolTipText("Screen pixels per trace");
    myTextZoomVertInch.setToolTipText( "Inches per second" );
    myTextZoomHorzInch.setToolTipText( "Traces per inch" );
    myTextMinorInc.setToolTipText("Minor increment for time axis annotation");
    myTextMajorInc.setToolTipText("Major increment for time axis annotation");
    myTextTraceClip.setToolTipText("Wiggles are clipped when extending over more than the specified number of traces");

    myBoxShowHorzLines.setToolTipText("Show time lines on top of seismic display");
    myBoxAutoHorzLines.setToolTipText("Automatically set time annotation increments");
    myBoxShowZeroLines.setToolTipText("Display straight line at each seismic trace (wiggle trace with zero amplitude)");
    myBoxTraceClip.setToolTipText("Perform trace clipping");
    myBoxWiggle.setToolTipText("Show wiggle trace");
    myBoxPosFill.setToolTipText("Fill positive wiggle");
    myBoxNegFill.setToolTipText("Full negative wiggle");
    myBoxVI.setToolTipText("Show variable intensity display");

    myButtonAutoRange.setToolTipText("Automatically set min/max range");
    myBoxUseVarColor.setToolTipText("Use variable color scale to fill wiggle traces");
    
    myButtonScalar.setToolTipText("Scale data by constant scalar");
    myButtonRange.setToolTipText("Scale data by specified min/max range");
    myButtonTrace.setToolTipText("Equalise each trace to an amplitude of 1, then scale by the specified trace scalar.");

    myButtonWiggleLinear.setToolTipText("Use linear interpolation between samples");
    myButtonWiggleCubic.setToolTipText("Use cubic interpolation between samples");
    
    myButtonPolarityNormal.setToolTipText("<html>'Normal' polarity: <b>Right</b>-hand side of wiggle is <b>positive</b>");
    myButtonPolarityReversed.setToolTipText("<html>'Reverse' polarity: <b>Left</b>-hand side of wiggle is <b>positive</b>");

    myButtonVIDiscrete.setToolTipText("Plot each sample as a discrete rectangle with constant color");
    myButtonVIVertical.setToolTipText("Interpolate each trace vertically, using the specified vertical interpolation (linear/cubic)");
    myButtonVI2DSpline.setToolTipText("Use 2D spline interpolation");

    myButtonMaximum.setToolTipText("Equalise data by trace maximum value");
    myButtonAverage.setToolTipText("Equalise data by trace average value");
    
    // Color/color map selectors
//    myComboColorMapWiggle.setToolTipText("");
//    myComboColorMapVA.setToolTipText("");
    myButtonColorWigglePos.setToolTipText("Select color for positive wiggle fill");
    myButtonColorWiggleNeg.setToolTipText("Select color for negative wiggle fill");
    myButtonHighlightColor.setToolTipText("Select color for trace highlighting");

//    myComboPlotDir.setToolTipText("Plot direction");
    myLabelPlotDir.setToolTipText("Plot direction");

    myLabelMinValue.setToolTipText("Minimum value for range setting");
    myLabelMaxValue.setToolTipText("Maximum value for range setting");
    myLabelMinorInc.setToolTipText("Minor increment for time axis annotation");
    myLabelMajorInc.setToolTipText("Major increment for time axis annotation");

    myLabelZoomVert.setToolTipText("Sample spacing in pixels per sample");
    myLabelZoomHorz.setToolTipText("Trace spacing in pixels per trace");

    myLabelColorMap.setToolTipText("Color map for variable intensity display");
  }
  public void setVisible( boolean isVisible ) {
     this.getRootPane().setDefaultButton(myButtonApply);

    // remove the binding for pressed
    getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
            .put(KeyStroke.getKeyStroke("ENTER"), "none");
    // retarget the binding for released
    getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW)
            .put(KeyStroke.getKeyStroke("released ENTER"), "press");

  	super.setVisible(isVisible);
    if( isVisible ) {
      SwingUtilities.invokeLater( new Runnable() { 
        public void run() { 
//          requestFocusInWindow();
//          myBoxWiggle.requestFocusInWindow();
        }
      });
    }
  }

  public static void main( String[] args ) {
    try {
      UIManager.setLookAndFeel( "com.sun.java.swing.plaf.windows.WindowsLookAndFeel" );
//      UIManager.setLookAndFeel( new MetalLookAndFeel() );
    }
    catch( Exception e ) {
      e.printStackTrace();
    }
    
    csSeisView seisview = new csSeisView(null,null);
    csSeisDispDialog dialog = new csSeisDispDialog( seisview, csColorMap.COLOR_MAP_TYPE_32BIT );
    dialog.setVisible(true);
    
    dialog.addWindowListener(new WindowAdapter() {
      public void windowClosing( WindowEvent e ) {
        System.exit( 0 );
      }
    });    
  }
  @Override
  public void mouseExited() {
  }
}



/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import javax.swing.*;
import java.awt.*;

/**
 * Color map with static list of default color maps.<br>
 * <ul>
 * <li> Defines a number of fixed colors at key positions
 * <li> Interpolates colors to every position inside a given range
 * </ul>
 * @author Bjorn Olofsson
 */
public class csColorMap {
  public static final int DEFAULT = 0;
  public static final int GRAY_WB = 1;
  public static final int GRAY_BW = 2;
  public static final int BLUE_WHITE_RED = 3;
  public static final int RAINBOW = 4;
  public static final int BLACK_WHITE_ORANGE = 5;
  public static final int GRAY_BWB = 6;
  public static final int GRAY_WBW = 7;
  public static final int RAINBOW_BLACK = 8;
  public static final int RAINBOW_BLACK_REV = 9;
  public static final int RAINBOW_MIRROR = 10;
  public static final int COLD_WARM  = 11;
  public static final int BLUE_WHITE_RED2 = 12;
  public static final int BLACK_WHITE_RED = 13;
  public static final int BROWN = 14;
  public static final int NUM_DEFAULT_MAPS = 15;

  public static final String TEXT_DEFAULT        = "default";
  public static final String TEXT_GRAY_WB        = "gray_w2b";
  public static final String TEXT_GRAY_BW        = "gray_b2w";
  public static final String TEXT_BLUE_WHITE_RED = "blue_white_red";
  public static final String TEXT_RAINBOW        = "rainbow";
  public static final String TEXT_BLACK_WHITE_ORANGE = "black_white_orange";
  public static final String TEXT_GRAY_BWB        = "gray_bwb";
  public static final String TEXT_GRAY_WBW        = "gray_wbw";
  public static final String TEXT_RAINBOW_BLACK   = "rainbow_black";
  public static final String TEXT_RAINBOW_BLACK_REV = "rainbow_2";
  public static final String TEXT_RAINBOW_MIRROR  = "rainbow_mirror";
  public static final String TEXT_COLD_WARM   = "cold_warm";
  public static final String TEXT_BLUE_WHITE_RED2 = "blue_white_red2";
  public static final String TEXT_BLACK_WHITE_RED = "black_white_red";
  public static final String TEXT_BROWN = "brown";
  
  public static final String[] TEXT = {
    TEXT_DEFAULT,
    TEXT_GRAY_WB,
    TEXT_GRAY_BW,
    TEXT_BLUE_WHITE_RED,
    TEXT_RAINBOW,
    TEXT_BLACK_WHITE_ORANGE,
    TEXT_GRAY_BWB,
    TEXT_GRAY_WBW,
    TEXT_RAINBOW_BLACK,
    TEXT_RAINBOW_BLACK_REV,
    TEXT_RAINBOW_MIRROR,
    TEXT_COLD_WARM,
    TEXT_BLUE_WHITE_RED2,
    TEXT_BLACK_WHITE_RED,
    TEXT_BROWN
  };

  public static final int NUM_COLORS_DEFAULT = 101;
  
  public static final int METHOD_SCALAR = 1;  // Global scalar
  public static final int METHOD_RANGE  = 2;  // Global range
  public static final int METHOD_TRACE  = 3;  // Trace-by-trace scalar

  public static final int COLOR_MAP_TYPE_32BIT = 91;
  public static final int COLOR_MAP_TYPE_8BIT  = 92;

  public static final int SMOOTH_MODE_DISCRETE     = 101;
  public static final int SMOOTH_MODE_CONTINUOUS   = 102;

  protected int myNumColors;
  protected double myValueRange;
  protected double myScalar;
  protected double myMinValue;
  protected double myMaxValue;
  /// Map index of color map. If the color map is user defined, the index is set to NUM_DEFAULT_MAPS
  protected int myDefaultMapIndex;
  protected int myScaleMethod;
  protected int mySmoothMode = csColorMap.SMOOTH_MODE_CONTINUOUS;
  protected int myDiscreteStep = 1;
  
  protected Color[] myColorKneePoints;
  protected double[] myWeightKneePoints;

  protected int myColorMapType;
  protected Color[] myColors;
  protected byte[]  myColors8bit;
  protected String myColorMapName;
  
  //----------------------------------------------------------------------
  // "Copy" constructor
  //
  public csColorMap( csColorMap map ) {
    resetAll( map );
  }
  //----------------------------------------------------------------------
  // Constructors with colors defined in input arguments
  //
  public csColorMap( Color color ) {
    Color [] colors = new Color[2];
    colors[0] = color;
    colors[1] = color;
    double[] weights = new double[2];
    weights[0] = 0.0;
    weights[1] = 1.0;
    myColorMapName = "USER DEFINED";
    setColors( colors, weights, colors.length );
  }
  public csColorMap( Color [] colors, int numColorResolution ) {
    myColorMapType = csColorMap.COLOR_MAP_TYPE_32BIT;
    double[] weights = new double[colors.length];
    if( colors.length == 1 ) {
      weights[0] = 1.0;
    }
    else {
      for( int i = 0; i < colors.length; i++ ) {
        weights[i] = (double)i/(double)(colors.length-1);
      }
    }
    myColorMapName = "USER DEFINED";
    setColors( colors, weights, numColorResolution );
  }
  public csColorMap( Color [] colors, double[] weights ) {
    this( colors, weights, weights.length );
  }
  public csColorMap( Color [] colors, double[] weights, int numColorResolution ) {
    myColorMapType = csColorMap.COLOR_MAP_TYPE_32BIT;
    myColorMapName = "USER DEFINED";
    setColors( colors, weights, numColorResolution );
  }
  //----------------------------------------------------------------------
  // Constructors with predefined color map index:
  //
  public csColorMap() {
    this( DEFAULT, COLOR_MAP_TYPE_32BIT );
  }
  public csColorMap( int predefColorMapIndex ) {
    this( predefColorMapIndex, COLOR_MAP_TYPE_32BIT, NUM_COLORS_DEFAULT );
  }
  public csColorMap( int predefColorMapIndex, int colorMapType ) {
    this( predefColorMapIndex, colorMapType, NUM_COLORS_DEFAULT );
  }
  public csColorMap( int predefColorMapIndex, int colorMapType, int numColors_in ) {
    myNumColors   = 0;
    myScalar      = 1.0;
    myMinValue    = -1.0;
    myMaxValue    = 1.0;
    myValueRange  = myMaxValue-myMinValue;
    myScaleMethod = METHOD_RANGE;
    myColorMapType = colorMapType;
    myColors8bit  = null;
    myColors      = null;
    
    Color[] colors;
    double[] weights = null;

    if( predefColorMapIndex == DEFAULT ) {
      colors = new Color[5];
      colors[0] = Color.cyan;
      colors[1] = Color.black;
      colors[2] = Color.white;
      colors[3] = Color.red;
      colors[4] = Color.yellow;
      weights = new double[5];
      weights[0] = 0.0;
      weights[1] = 0.05;
      weights[2] = 0.5;
      weights[3] = 0.95; 
      weights[4] = 1.0;
    }
    else if( predefColorMapIndex == GRAY_WB ) {
      colors = new Color[2];
      colors[0] = Color.white;
      colors[1] = Color.black;
    }
    else if( predefColorMapIndex == GRAY_BW ) {
      colors = new Color[2];
      colors[0] = Color.black;
      colors[1] = Color.white;
    }
    else if( predefColorMapIndex == BLUE_WHITE_RED ) {
      colors  = new Color[3];
      colors[0] = Color.blue;
      colors[1] = Color.white;
      colors[2] = Color.red;
    }
    else if( predefColorMapIndex == BLUE_WHITE_RED2 ) {
      colors  = new Color[15];
      colors[0] = new Color(0,0,144);
      colors[1] = new Color(5,5,154);
      colors[2] = new Color(10,10,164);
      colors[3] = new Color(16,16,174);
      colors[4] = new Color(120,120,184);
      colors[5] = new Color(165,165,216);
      colors[6] = new Color(215,215,245);
      colors[7] = new Color(255,255,255);
      colors[8] = new Color(255,215,215);
      colors[9] = new Color(246,165,165);
      colors[10] = new Color(224,120,120);
      colors[11] = new Color(214,16,16);
      colors[12] = new Color(204,10,10);
      colors[13] = new Color(194,5,5);
      colors[14] = new Color(184,0,0);
      weights = new double[15];
      weights[0] = 0.0;
      weights[1] = 0.1;
      weights[2] = 0.2;
      weights[3] = 0.3;
      weights[4] = 0.35;
      weights[5] = 0.4;
      weights[6] = 0.45;
      weights[7] = 0.5;
      weights[8] = 0.55;
      weights[9] = 0.6;
      weights[10] = 0.65;
      weights[11] = 0.7;
      weights[12] = 0.8;
      weights[13] = 0.9;
      weights[14] = 1.0;
    }
    else if( predefColorMapIndex == RAINBOW ) {
      colors  = new Color[6];
      colors[0] = Color.blue;
      colors[1] = Color.cyan;
      colors[2] = Color.green;
      colors[3] = Color.yellow;
      colors[4] = Color.red;
      colors[5] = Color.magenta;
    }
    else if( predefColorMapIndex == RAINBOW_BLACK ) {
      colors  = new Color[7];
      colors[0] = Color.black;
      colors[1] = Color.blue;
      colors[2] = Color.cyan;
      colors[3] = Color.green;
      colors[4] = Color.yellow;
      colors[5] = Color.red;
      colors[6] = Color.magenta;
    }
    else if( predefColorMapIndex == COLD_WARM ) {
      colors  = new Color[6];
      colors[0] = Color.black;
      colors[1] = Color.blue;
      colors[2] = Color.white;
      colors[3] = Color.yellow;
      colors[4] = Color.red;
      colors[5] = Color.magenta;
      int numColors = colors.length;
      weights = new double[numColors];
      weights[0] = 0.0;
      weights[1] = 0.24;
      weights[2] = 0.48;
      weights[3] = 0.76;
      weights[4] = 0.90;
      weights[5] = 1.0;
    }
    else if( predefColorMapIndex == BLACK_WHITE_ORANGE ) {
      colors  = new Color[3];
      colors[0] = Color.black;
      colors[1] = Color.white;
      colors[2] = new Color( 255, 130, 0 );
    }
    else if( predefColorMapIndex == BLACK_WHITE_RED ) {
      colors  = new Color[6];
      colors[0] = Color.black;
      colors[1] = Color.darkGray;
      colors[2] = Color.white;
      colors[3] = Color.yellow;
      colors[4] = new Color(224,30,30);
      colors[5] = new Color(204,10,10);
      weights = new double[6];
      weights[0] = 0.0;
      weights[1] = 0.4;
      weights[2] = 0.5;
      weights[3] = 0.55;
      weights[4] = 0.6;
      weights[5] = 1.0;
    }
    else if( predefColorMapIndex == GRAY_BWB ) {
      colors = new Color[3];
      colors[0] = Color.black;
      colors[1] = Color.white;
      colors[2] = Color.black;
    }
    else if( predefColorMapIndex == GRAY_WBW ) {
      colors = new Color[3];
      colors[0] = Color.white;
      colors[1] = Color.black;
      colors[2] = Color.white;
    }
    else if( predefColorMapIndex == RAINBOW_MIRROR ) {
      colors  = new Color[11];
      colors[0] = Color.magenta;
      colors[1] = Color.red;
      colors[2] = Color.yellow;
      colors[3] = Color.green;
      colors[4] = Color.cyan;
      colors[5] = Color.blue;
      colors[6] = Color.cyan;
      colors[7] = Color.green;
      colors[8] = Color.yellow;
      colors[9] = Color.red;
      colors[10] = Color.magenta;
    }
    else if( predefColorMapIndex == RAINBOW_BLACK_REV ) {
      colors  = new Color[6];
      colors[0] = Color.blue;
      colors[1] = Color.green;
      colors[2] = Color.yellow;
      colors[3] = Color.red;
      colors[4] = Color.magenta;
      colors[5] = Color.black;
    }
    else if( predefColorMapIndex == BROWN ) {
      colors  = new Color[64];
      colors[0] = new Color(255,255,0);
      colors[1] = new Color(238,224,3);
      colors[2] = new Color(220,193,6);
      colors[3] = new Color(203,162,9);
      colors[4] = new Color(185,131,13);
      colors[5] = new Color(168,100,16);
      colors[6] = new Color(150,69,19);
      colors[7] = new Color(155,76,28);
      colors[8] = new Color(159,84,38);
      colors[9] = new Color(163,91,47);
      colors[10] = new Color(167,98,56);
      colors[11] = new Color(171,106,66);
      colors[12] = new Color(175,113,75);
      colors[13] = new Color(179,121,84);
      colors[14] = new Color(184,128,94);
      colors[15] = new Color(188,135,103);
      colors[16] = new Color(192,143,112);
      colors[17] = new Color(196,150,122);
      colors[18] = new Color(200,157,131);
      colors[19] = new Color(204,165,140);
      colors[20] = new Color(208,172,150);
      colors[21] = new Color(213,179,159);
      colors[22] = new Color(217,187,169);
      colors[23] = new Color(221,194,178);
      colors[24] = new Color(225,202,187);
      colors[25] = new Color(229,209,197);
      colors[26] = new Color(233,216,206);
      colors[27] = new Color(237,224,215);
      colors[28] = new Color(242,231,225);
      colors[29] = new Color(246,238,234);
      colors[30] = new Color(250,246,243);
      colors[31] = new Color(254,253,253);
      colors[32] = new Color(247,247,247);
      colors[33] = new Color(237,237,237);
      colors[34] = new Color(226,226,226);
      colors[35] = new Color(216,216,216);
      colors[36] = new Color(206,206,206);
      colors[37] = new Color(195,195,195);
      colors[38] = new Color(185,185,185);
      colors[39] = new Color(174,174,174);
      colors[40] = new Color(164,164,164);
      colors[41] = new Color(153,153,153);
      colors[42] = new Color(143,143,143);
      colors[43] = new Color(133,133,133);
      colors[44] = new Color(122,122,122);
      colors[45] = new Color(112,112,112);
      colors[46] = new Color(101,101,101);
      colors[47] = new Color(91,91,91);
      colors[48] = new Color(81,81,81);
      colors[49] = new Color(70,70,70);
      colors[50] = new Color(60,60,60);
      colors[51] = new Color(49,49,49);
      colors[52] = new Color(39,39,39);
      colors[53] = new Color(29,29,29);
      colors[54] = new Color(18,18,18);
      colors[55] = new Color(8,8,8);
      colors[56] = new Color(0,4,10);
      colors[57] = new Color(0,21,49);
      colors[58] = new Color(0,38,88);
      colors[59] = new Color(0,55,127);
      colors[60] = new Color(0,72,167);
      colors[61] = new Color(0,89,206);
      colors[62] = new Color(0,107,245);
      colors[63] = new Color(0,111,255);
    }
    else {
      // Must be a custom color map...
      colors = new Color[2];
      colors[0] = Color.white;
      colors[1] = Color.black;
      myColorMapName = "UNKNOWN";
    }
    if( predefColorMapIndex >= 0 && predefColorMapIndex < NUM_DEFAULT_MAPS ) {
      myColorMapName = TEXT[predefColorMapIndex];
    }
    else {
      myColorMapName = "UNKNOWN";
    }
    if( weights == null ) {
      int numColors = colors.length;
      weights = new double[numColors];
      for( int i = 0; i < numColors; i++ ) {
        weights[i] = (double)i/(double)(numColors-1);
      }
    }
    setColors( colors, weights, numColors_in );
    myDefaultMapIndex = predefColorMapIndex;
  }
  public void dump() {
    System.out.println(" *** Color map dump *** ");
    System.out.println(" Color map type:    " + myColorMapType + " numColors: " + myColors.length);
    System.out.println(" Color smooth mode: " + mySmoothMode);
    for( int i = 0; i < myColorKneePoints.length; i++ ) {
      System.out.println("Weight #" + (i+1) + ": " + myWeightKneePoints[i] + " color: " + myColorKneePoints[i]);
    }
  }
  public void setName( String name ) {
    myColorMapName = name;
  }
  public void setSmoothModeDiscrete( int discreteStep ) {
    mySmoothMode   = csColorMap.SMOOTH_MODE_DISCRETE;
    myDiscreteStep = discreteStep;
  }
  public void setSmoothModeContinuous() {
    mySmoothMode = csColorMap.SMOOTH_MODE_CONTINUOUS;
  }
  public int getSmoothMode() {
    return mySmoothMode;
  }
  public int getDiscreteStep() {
    return myDiscreteStep;
  }
  public void setColorMapType( int type ) {
    if( myColorMapType != type ) {
      myColorMapType = type;
      if( myColorMapType == csColorMap.COLOR_MAP_TYPE_32BIT ) {
        myColors8bit = null;
      }
      else {
        computeColors8bit();
      }
    }
  }
  public int getNumColorResolution() {
    if( myColors == null ) return 0;
    return myColors.length;
  }
  public void setNumColorResolution( int numColors ) {
    if( numColors < myColorKneePoints.length ) numColors = myColorKneePoints.length;
    setColors( myColorKneePoints, myWeightKneePoints, numColors );
  }
  /**
   * Reset all parameters to the values of the specified color map.
   * @param map Color map
   */
  protected void resetAll( csColorMap map ) {
    resetColors( map );
    myScalar     = map.myScalar;
    myMinValue   = map.myMinValue;
    myMaxValue   = map.myMaxValue;
    myValueRange = map.myValueRange;
    myScaleMethod= map.myScaleMethod;
    mySmoothMode = map.mySmoothMode;
  }
  /**
   * Reset colors to the colors of the specified color map. Do not reset any other parameters.
   * @param map Color map giving the new colors
   */
  public void resetColors( csColorMap map ) {
    if( map == this ) return;
    myNumColors       = map.myNumColors;
    myColorMapType    = map.myColorMapType;
    myDefaultMapIndex = map.myDefaultMapIndex;
    myColorKneePoints = map.myColorKneePoints;
    myWeightKneePoints = map.myWeightKneePoints;
    myColorMapName = map.myColorMapName;

    Color[] colors = new Color[myNumColors];
    System.arraycopy(map.myColors, 0, colors, 0, myNumColors);
    myColors = colors;
    if( myColorMapType == COLOR_MAP_TYPE_8BIT) {
      computeColors8bit();
    }
  }
  /**
   * Set specified 
   * @param colors
   * @param colorWeights
   * @param numColors
   */
  protected void setColors( Color[] colors, double[] colorWeights, int numColors ) {
    myDefaultMapIndex  = NUM_DEFAULT_MAPS;
    myNumColors        = numColors;
    myColorKneePoints  = colors;
    myWeightKneePoints = colorWeights;
    myColors = csColorMap.computeColors( myColorKneePoints, myWeightKneePoints, myNumColors );
    if( myColorMapType == COLOR_MAP_TYPE_8BIT) {
      computeColors8bit();
    }
  }
  /**
   * Compute the specified number of interpolated colors for the specified set of color knee points
   * and associated weights/distances.
   * 
   * @param colors       N colors
   * @param colorWeights N weights, values increasing from 0-1
   * @param numColors    Number of output colors
   * @return Computed colors
   */
  protected static Color[] computeColors( Color[] colors, double[] colorWeights, int numColors ) {
    Color[] colorsOut  = new Color[numColors];

    int index = 0;
    for( int i = 0; i < numColors; i++ ) {
      float weight = (float)i / (float)(numColors-1);
      while( index < (colorWeights.length-1) && weight >= colorWeights[index] ) index += 1;
      int redMin = colors[index-1].getRed();
      int blueMin = colors[index-1].getBlue();
      int greenMin = colors[index-1].getGreen();
      int redMax = colors[index].getRed();
      int blueMax = colors[index].getBlue();
      int greenMax = colors[index].getGreen();
      float color_weight_step = (float)(colorWeights[index] - colorWeights[index-1]);
      
      if( color_weight_step > 0 ) weight  = ( weight - (float)colorWeights[index-1] ) / color_weight_step;
      else weight = 0.0f;
      int red   = (int) ( redMin + weight * ( redMax - redMin ) );
      int green = (int) ( greenMin + weight * ( greenMax - greenMin ) );
      int blue  = (int) ( blueMin + weight * ( blueMax - blueMin ) );
  
      colorsOut[i] = new Color( red, green, blue );
    }
    return colorsOut;
  }
  /**
   * Set min/max values
   * @param minValue
   * @param maxValue
   */
  public void setMinMax( double minValue, double maxValue ) {
    myScaleMethod = METHOD_RANGE;
    setMinMax_internal( minValue, maxValue );
  }
  private void setMinMax_internal( double minValue, double maxValue ) {
    myMinValue = minValue;
    myMaxValue = maxValue;
    myValueRange = myMaxValue-myMinValue;
    if( myValueRange == 0 ) myValueRange = 1;
  }
  /**
   * Set constant scalar. This scalar is applied to each value for which the corresponding color shall be retrieved.
   * @param scalar
   */
  public void setScalar( double scalar ) {
    myScaleMethod = METHOD_SCALAR;
    myScalar = Math.abs( scalar );
    setMinMax_internal( -myScalar, myScalar );
  }
  public int getDefaultMapIndex() {
    return myDefaultMapIndex;
  }
  public int getMethod() {
    return myScaleMethod;
  }
  public int getColorMapType() {
    return myColorMapType;
  }
  /**
   * 
   * @param value
   * @return RGB color at given value
   */
  public int getColorRGB( float value ) {
    return myColors[ bufferIndex(value) ].getRGB();
  }
  /**
   * 
   * @param value
   * @return Color at given value
   */
  public Color getColor( float value ) {
    return myColors[ bufferIndex(value) ];
  }
  public Color[] getColors() {
    return myColors;
  }
  public int getNumKneePoints() {
    return myColorKneePoints.length;
  }
  public Color getColorKneePoint( int index ) {
    return myColorKneePoints[index];
  }
  public double getWeightKneePoint( int index ) {
    return myWeightKneePoints[index];
  }
  public Color[] getColorKneePoints() {
    return myColorKneePoints;
  }
  public double[] getWeightKneePoints() {
    return myWeightKneePoints;
  }
  /**
   * 
   * @param value
   * @return RGB color at given value
   */
  protected int bufferIndex( float value ) {
    if( mySmoothMode == csColorMap.SMOOTH_MODE_DISCRETE ) {
      value = (float)( myDiscreteStep * Math.round( value / (float)myDiscreteStep ) );
    }
    double weight;
    switch( myScaleMethod ) {
      case METHOD_RANGE:
        weight = ( value - myMinValue ) / myValueRange;
        break;
      case METHOD_SCALAR:
        weight = (1 + value*myScalar)*0.5;
        break;
      default:
        weight = 0;
    }
    int index = (int)( weight * (double)(myNumColors-1) + 0.5 );
    if( index < 0 ) index = 0;
    if( index >= myNumColors ) index = myNumColors-1;
    return index;
  }
  public double getScalar() {
    return myScalar;
  }
  public double getMinValue() {
    return myMinValue;
  }
  public double getMaxValue() {
    return myMaxValue;
  }
  public static String textDefaultMap( int i ) {
    if( i >= 0 && i < NUM_DEFAULT_MAPS ) {
      return TEXT[i];
    }
    return "UNKNOWN";
  }
  //*******************************************************************
  //
  // 8bit methods
  //
  //*******************************************************************
  protected void computeColors8bit() {
    myColors8bit = new byte[myNumColors];
    for( int i = 0; i < myNumColors; i++ ) {
      int red   = myColors[i].getRed();
      int green = myColors[i].getGreen();
      int blue  = myColors[i].getBlue();
      myColors8bit[i] = csColorMap.convertColorTo8bit( red, green, blue );
    }
  }
  /**
   * 
   * @param value
   * @return Indexed byte color at given value
   */
  public byte getColor8bit( float value ) {
    return myColors8bit[bufferIndex( value )];
  }
  /**
   * Convert color to indexed byte value, according to default indexed byte color model (see BufferedImage TYPE_BYTE_INDEXED)
   * 
   * @param color
   * @return
   */
  public static byte convertColorTo8bit( Color color ) {
    return convertColorTo8bit( color.getRed(), color.getGreen(), color.getBlue() );
  }
  public static byte convertColorTo8bit( int red, int green, int blue ) {
    if( red != green || red != blue ) {
      return (byte)( 36*(int)(red/51.0 + 0.5) + 6*(int)(green/51.0 + 0.5) + (int)(blue/51.0 + 0.5) );
    }
    else {
      return (byte)( 216 + (int)(red/6.54 + 0.5) );
    }
  }
  public static int getDefaultMapIndex( String text ) {
    for( int i = 0; i < TEXT.length; i++ ) {
      if( text.compareTo( TEXT[i] ) == 0 ) {
        return i;
      }
    }
    return -1;
  }
  @Override
  public String toString() {
    return myColorMapName;
  }
  //-------------------------------------------------
  //
  public static void main( String[] args ) {
    JDialog dialog = new JDialog();
    @SuppressWarnings("serial")
    JPanel panel = new JPanel( new BorderLayout() ) {
      @Override
      public void paintComponent( Graphics g ) {
        int height = getHeight();
        int width = getWidth();
        csColorMap map = new csColorMap( csColorMap.RAINBOW, csColorMap.COLOR_MAP_TYPE_32BIT );
        map.setMinMax( 0, height-1 );
        map.setSmoothModeDiscrete(height/4);
        for( int i = 0; i < height; i++ ) {
          int colorRGB = map.getColorRGB(i);
          g.setColor( new Color(colorRGB) );
          g.drawLine(0, i, width-1, i);
        }
      }
    };
    dialog.getContentPane().add( panel );
    dialog.setSize( 100, 600 );
    dialog.setVisible(true);
    
  }
}



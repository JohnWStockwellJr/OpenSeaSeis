/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import java.awt.Color;

/**
 * 8-bit version of color map.
 * @author Bjorn Olofsson
 */
public class csColorMap8bit {
  private csColorMap myMap;
  private byte[] myColorsByte;
  
  public csColorMap8bit() {
    myMap = new csColorMap();
    computeColors();
  }
  public csColorMap8bit( csColorMap map ) {
    myMap = map;
    computeColors();
  }
  public void reset( csColorMap map ) {
    myMap = map;
    computeColors();
  }
  protected void computeColors() {
    Color[] colors = myMap.getColors();
    int numColors = colors.length;
    myColorsByte = new byte[numColors];
    for( int i = 0; i < numColors; i++ ) {
      int red   = colors[i].getRed();
      int green = colors[i].getGreen();
      int blue  = colors[i].getBlue();
      myColorsByte[i] = csColorMap8bit.convertColorTo8bit( red, green, blue );
    }
  }
  /**
   * 
   * @param value
   * @return Indexed byte color at given value
   */
  public byte getColor8bit( float value ) {
    int index = myMap.bufferIndex( value );
    return myColorsByte[index];
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

}



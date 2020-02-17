/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;

/**
 * Custom color map.
 * @author 2013 Felipe Punto
 */
public class csCustomColorMap extends csColorMap {
  private static final String TEXT_IDENTIFIER = "#CMAP";
  private static final String DEFAULT_NAME = "custom_map01";

  public csCustomColorMap( BufferedReader reader ) throws IOException {
    super( -1 );
    load( reader );
  }
  public csCustomColorMap() {
    super( -1 );
    myColorMapName = DEFAULT_NAME;
  }
  public csCustomColorMap( String name ) {
    super( -1 );
    myColorMapName = name;
  }
  public csCustomColorMap( Color color ) {
    super( color );
    myColorMapName = DEFAULT_NAME;
  }
  public csCustomColorMap( Color [] colors, int numColorResolution ) {
    super( colors, numColorResolution );
    myColorMapName = DEFAULT_NAME;
  }
  public csCustomColorMap( Color [] colors, double[] weights, int numColorResolution ) {
    super( colors, weights, numColorResolution );
    myColorMapName = DEFAULT_NAME;
  }
  public csCustomColorMap( csColorMap map ) {
    super( map );
  }
  public void write( BufferedWriter writer ) throws IOException {
    writer.write( TEXT_IDENTIFIER + "\n");
    writer.write(myColorMapName + "\n");
    writer.write("" + getNumKneePoints() + " " + getNumColorResolution() + "\n");
    for( int i = 0; i < getNumKneePoints(); i++ ) {
      writer.write("" + myWeightKneePoints[i] + " " + myColorKneePoints[i].getRGB() + "\n");
    }
  }
  public void load( BufferedReader reader ) throws IOException {
    String line = reader.readLine();
    if( line == null ) throw new IOException("Input file is empty.");
    if( !line.startsWith( TEXT_IDENTIFIER ) ) {
      throw new IOException("Input file is not a valid CMAP file. First line must start with " + TEXT_IDENTIFIER);
    }
    myColorMapName = reader.readLine();
    line = reader.readLine();
    String[] tokens = line.trim().split(" ");
    if( tokens.length < 2 ) throw new IOException("Incorrect CMAP file format: Missing number of color");
    int numKneePoints = Integer.parseInt( tokens[0] );
    int numColorResolution = Integer.parseInt( tokens[1] );
    Color[] colorKneePoints = new Color[numKneePoints];
    double[] weightKneePoints = new double[numKneePoints];
    for( int i = 0; i < numKneePoints; i++ ) {
      line = reader.readLine();
      tokens = line.trim().split(" ");
      if( tokens.length < 2 ) throw new IOException("Incorrect CMAP file format: Missing weight/color");
      weightKneePoints[i] = Double.parseDouble( tokens[0] );
      colorKneePoints[i] = new Color( Integer.parseInt( tokens[1] ) );
    }
    setColors( colorKneePoints, weightKneePoints, numColorResolution );
  }
}


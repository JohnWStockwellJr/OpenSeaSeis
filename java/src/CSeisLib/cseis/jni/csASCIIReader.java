/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

import cseis.seis.csHeader;
import cseis.seis.csSeismicTrace;
import cseis.seis.csISeismicReader;
import cseis.seis.csHeaderDef;

import java.io.*;
import java.util.*;

import java.io.IOException;

/**
 * ASCII file reader.
 * @author 2011 Bjorn Olofsson
 */
 // TODO: ..this class should be moved to C++ and called through JNI

public class csASCIIReader implements csISeismicReader {
  private File myFile;
  private BufferedReader myReader;

  private int myNumHeaders;
  private int myNumSamples;
  private int myNumTraces;
  private float mySampleInt;
  private Number[][] myHeaderValues;
  private String[] myHeaderNames;
  private String[] myHeaderDesc;
  private float[][] myTraceSamples;

  private int myCurrentTrace;
  private csITraceHeaderScanNotifier mySelectionNotifier;
  private int mySelectionHdrIndex;

  public csASCIIReader( String filename ) throws IOException, Exception {
    try {
      myFile = new File(filename);
      myReader = new BufferedReader(new FileReader(myFile));
    }
    catch( Exception e ) {
      throw(e);
    }
    mySampleInt = 1000.0f;
    myNumHeaders = 1;
    myHeaderNames = new String[myNumHeaders];
    myHeaderDesc  = new String[myNumHeaders];
    myHeaderNames[0] = "trcno";
    myHeaderDesc[0]  = "Trace number";
    myNumSamples = 0;
    myNumTraces  = 0;
    myHeaderValues = null;
    myTraceSamples = null;
    myCurrentTrace = 0;
    readFile();
  }
  private void readFile() throws Exception {
    myNumTraces = 0;
    String line;
    ArrayList<float[]> sampleList = new ArrayList<float[]>();
    try {
      while ((line = myReader.readLine()) != null) {
        String[] tokens = line.replaceAll("\t"," ").trim().split(" ");
        if( tokens.length == 0 ) {
          continue; // Empty line
        }
        if( tokens[0].charAt(0) == '#' ) continue; // Comment line 
        if( myNumTraces == 0 ) {
          myNumTraces = tokens.length;
        }
        else if( myNumTraces != tokens.length ) {
          throw new Exception("Syntax error in input file: Unequal number of columns/traces (" +
                  tokens.length + ") in line " + sampleList.size() );
        }
        float[] samples = new float[myNumTraces];
        try {
          for( int i = 0; i < tokens.length; i++ ) {
            samples[i] = Float.parseFloat(tokens[i]);
          }
        } catch (NumberFormatException e) {
          e.printStackTrace();
        }
        sampleList.add(samples);
      }
    }
    catch( IOException e ) {
    }
    myNumSamples = sampleList.size();
    myHeaderValues = new Number[myNumHeaders][myNumTraces];
    myTraceSamples = new float[myNumTraces][myNumSamples];
    for( int itrc = 0; itrc < myNumTraces; itrc++ ) {
      myHeaderValues[0][itrc] = itrc;
    }
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      float[] samples = sampleList.get(isamp);
      for( int itrc = 0; itrc < myNumTraces; itrc++ ) {
        myTraceSamples[itrc][isamp] = samples[itrc];
      }
    }
  }
  @Override
  public void closeFile() {
    try {
      myReader.close();
    }
    catch( IOException e ) {
    }
  }
  public int sampleByteSize() {
    return 4;
  }
  @Override
  public int numTraces() {
    return myNumTraces;
  }
  @Override
  public int numSamples() {
    return myNumSamples;
  }
  @Override
  public int numHeaders() {
    return myNumHeaders;
  }
  @Override
  public float sampleInt() {
    return mySampleInt;
  }
  //-------------------------------------------------------------
  @Override
  public boolean getNextTrace( csSeismicTrace trace ) throws Exception {
    if( myNumTraces <= 0 ) return false;
    if( myCurrentTrace == myNumTraces-1 ) return false;

    //    System.arraycopy(this, myNumTraces, this, myNumTraces, myNumTraces);
    float[] samples = trace.samples();
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      samples[isamp] = myTraceSamples[myCurrentTrace][isamp];
    }
    if( myCurrentTrace < myNumTraces-1 ) myCurrentTrace += 1;
    return true;
  }
  @Override
  public boolean hasRandomFileAccess() {
    return true;
  }
  @Override
  public boolean moveToTrace( int traceIndex, int numTracesToRead ) throws Exception {
    if( myNumTraces <= 0 ) return false;
    myCurrentTrace = traceIndex;
    if( myCurrentTrace > myNumTraces-1 ) myCurrentTrace = myNumTraces-1;
    if( myCurrentTrace < 0 ) myCurrentTrace = 0;
    return true;
  }
  //-------------------------------------------------------------
  @Override
  public int hdrIntValue( int hdrIndex ) {
    return myHeaderValues[hdrIndex][myCurrentTrace].intValue();
  }
  @Override
  public float hdrFloatValue( int hdrIndex ) {
    return myHeaderValues[hdrIndex][myCurrentTrace].floatValue();
  }
  @Override
  public double hdrDoubleValue( int hdrIndex ) {
    return myHeaderValues[hdrIndex][myCurrentTrace].doubleValue();
  }
  @Override
  public String headerName( int hdrIndex ) {
    return myHeaderNames[hdrIndex];
  }
  @Override
  public String headerDesc( int hdrIndex ) {
    return myHeaderDesc[hdrIndex];
  }
  @Override
  public int headerType( int hdrIndex ) {
    return csHeaderDef.TYPE_INT;
  }
  @Override
  public boolean isFrequencyDomain() {
    return false;
  }
  @Override
  public int verticalDomain() {
    return cseis.general.csUnits.DOMAIN_TIME;
  }
  @Override
  public void setHeaderToPeek( String headerName ) throws Exception {
  }
  @Override
  public boolean peekHeaderValue( int traceIndex, csHeader value ) {
    int hdrIndex = 0;
    value.setValue( myHeaderValues[hdrIndex][traceIndex] );
    return true;
  }
  @Override
  public boolean setSelection( String hdrValueSelectionText, String headerName, int sortOrder, int sortMethod,
          csITraceHeaderScanNotifier notifier ) {
    // NOTE: Ignore sort method etc
    mySelectionNotifier = notifier;
    mySelectionHdrIndex = -1;
    
    myHeaderValues = new Number[myNumHeaders][myNumTraces];
    mySelectionHdrIndex = 0;  // There is only one header in this class which is trcno
    
    for( int ihdr = 0; ihdr < myNumHeaders; ihdr++ ) {
      if( myHeaderNames[ihdr].compareTo(headerName) == 0 ) {
        mySelectionHdrIndex = ihdr;
        mySelectionNotifier.traceHeaderScanNotify( numTraces()-1 );
        return true;
      }
    }
    mySelectionNotifier.traceHeaderScanNotify( numTraces()-1 );
    return false;
  }
  @Override
  public int getNumSelectedTraces() {
    return numTraces();
  }
  @Override
  public void getSelectedValues( csSelectedHeaderBundle hdrBundle ) {
    if( mySelectionHdrIndex < 0 ) return;
    for( int itrc = 0; itrc < numTraces(); itrc++ ) {
      hdrBundle.hdrValues[itrc].setValue( myHeaderValues[mySelectionHdrIndex][itrc] );
      hdrBundle.traceIndexList[itrc] = itrc;  // No sorting
    }
  }
  //----------------------------------------------
  public void disposeNative() {
  }
  public void finalize() {
    disposeNative();
  }
}



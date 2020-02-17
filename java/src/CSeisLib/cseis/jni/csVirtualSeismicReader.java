/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicReader;
import cseis.seis.csSeismicTrace;
import cseis.seis.csTraceBuffer;

/**
 * Virtual seismic reader.<br>
 * Actually stores all data in memory = trace buffer, instead of linking to a disk file or reader object.
 * @author 2013 Felipe Punto
 */
public class csVirtualSeismicReader implements csISeismicReader {
  protected csTraceBuffer myTraceBuffer;
  protected csHeaderDef[] myTraceHeaderDef;
  protected float mySampleInt;
  protected int myCurrentTraceIndex;
  protected int myPeekHeaderIndex;
  protected int myVerticalDomain;
  
  protected csITraceHeaderScanNotifier mySelectionNotifier;
  protected int mySelectionHdrIndex;
  
  public csVirtualSeismicReader() {   
  }
  /**
   * 
   * @param numSamples  Number of samples per trace
   * @param numHeaders  Number of trace headers per trace
   * @param sampleInt   Sample interval [ms]/[Hz]/[m]
   * @param headerDef   Trace header definitions: Names, descriptions, types of all trace headers
   * @param verticalDomain     Vertical domain of data: Time, depth, frequency
   */
  public csVirtualSeismicReader( int numSamples, int numHeaders, float sampleInt,
          csHeaderDef[] headerDef, int verticalDomain )
  {
    myTraceBuffer = new csTraceBuffer( numSamples, numHeaders );
    mySampleInt = sampleInt;
    myCurrentTraceIndex = 0;
    myPeekHeaderIndex = 0;
    mySelectionNotifier = null;
    myVerticalDomain = verticalDomain;

    myTraceHeaderDef = new csHeaderDef[headerDef.length];
    for( int i = 0; i < myTraceHeaderDef.length; i++ ) {
      myTraceHeaderDef[i] = new csHeaderDef( headerDef[i] );
    }
  }
//  public csVirtualSeismicReader( csISeismicTraceBuffer traceBuffer_in ) {
  public csTraceBuffer retrieveTraceBuffer() {
    return myTraceBuffer;
  }
  @Override
  public int numTraces() {
    return myTraceBuffer.numTraces();
  }
  @Override
  public int numSamples() {
    return myTraceBuffer.numSamples();
  }
  @Override
  public int numHeaders() {
    return myTraceBuffer.numHeaders();
  }
  @Override
  public float sampleInt() {
    return mySampleInt;
  }

  @Override
  public int verticalDomain() {
    return myVerticalDomain;
  }

  @Override
  public boolean isFrequencyDomain() {
    return( myVerticalDomain == cseis.general.csUnits.DOMAIN_FREQ );
  }

  @Override
  public boolean getNextTrace( csSeismicTrace trace ) throws Exception {
    if( myCurrentTraceIndex < numTraces() ) {
      System.arraycopy( myTraceBuffer.samples(myCurrentTraceIndex), 0, trace.samples(), 0, numSamples() );
      myCurrentTraceIndex += 1;
      return true;
    }
    return false;
  }
  @Override
  public boolean hasRandomFileAccess() {
    return true;
  }

  @Override
  public boolean moveToTrace( int traceIndex, int numTracesToRead ) throws Exception {
    if( traceIndex >= 0 && traceIndex < numTraces() ) {
      myCurrentTraceIndex = traceIndex;
      return true;
    }
    return false;
  }

  @Override
  public int hdrIntValue(int hdrIndex) {
    csHeader header = myTraceBuffer.headerValues(myCurrentTraceIndex)[hdrIndex];
    return header.intValue();
  }

  @Override
  public float hdrFloatValue(int hdrIndex) {
    csHeader header = myTraceBuffer.headerValues(myCurrentTraceIndex)[hdrIndex];
    return header.floatValue();
  }

  @Override
  public double hdrDoubleValue(int hdrIndex) {
    csHeader header = myTraceBuffer.headerValues(myCurrentTraceIndex)[hdrIndex];
    return header.doubleValue();
  }

  @Override
  public String headerName(int hdrIndex) {
    return myTraceHeaderDef[hdrIndex].name;
  }

  @Override
  public String headerDesc(int hdrIndex) {
    return myTraceHeaderDef[hdrIndex].desc;
  }

  @Override
  public int headerType(int hdrIndex) {
    return myTraceHeaderDef[hdrIndex].type;
  }

  @Override
  public void setHeaderToPeek( String headerName ) throws Exception {
    myPeekHeaderIndex = 0;
  }

  @Override
  public boolean peekHeaderValue( int traceIndex, csHeader value ) {
    if( traceIndex >= 0 && traceIndex < numTraces() ) {
      value = myTraceBuffer.headerValues(traceIndex)[myPeekHeaderIndex];
      return true;
    }
    return false;
  }
  //
  @Override
  public boolean setSelection( String hdrValueSelectionText, String headerName, int sortOrder, int sortMethod,
          csITraceHeaderScanNotifier notifier ) {
    // NOTE: Ignore sort method etc
    mySelectionNotifier = notifier;
    mySelectionHdrIndex = -1;
    for( int ihdr = 0; ihdr < myTraceHeaderDef.length; ihdr++ ) {
      if( myTraceHeaderDef[ihdr].name.compareTo(headerName) == 0 ) {
        mySelectionHdrIndex = ihdr;
        mySelectionNotifier.traceHeaderScanNotify( myTraceBuffer.numTraces()-1 );
        return true;
      }
    }
    mySelectionNotifier.traceHeaderScanNotify( myTraceBuffer.numTraces()-1 );
    return false;
  }
  @Override
  public int getNumSelectedTraces() {
    return myTraceBuffer.numTraces();
  }
  @Override
  public void getSelectedValues( csSelectedHeaderBundle hdrBundle ) {
    if( mySelectionHdrIndex < 0 ) return;
    for( int itrc = 0; itrc < myTraceBuffer.numTraces(); itrc++ ) {
      hdrBundle.hdrValues[itrc] = myTraceBuffer.headerValues(itrc)[mySelectionHdrIndex];
      hdrBundle.traceIndexList[itrc] = itrc;  // No sorting
    }
  }

  @Override
  public void closeFile() {
    // Nothing really to do
    myTraceBuffer = null;
  }
}


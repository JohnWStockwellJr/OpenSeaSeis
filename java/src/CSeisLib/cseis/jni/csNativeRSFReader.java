/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

import cseis.seis.csHeader;
import cseis.seis.csSeismicTrace;
import cseis.seis.csISeismicReader;

import java.io.IOException;

/**
 * JNI Interface to RSF reader (C++ code)
 * @author Bjorn Olofsson
 */

public class csNativeRSFReader implements csISeismicReader {
  private native long native_createInstance( String filename );
  private native void native_closeFile( long ptr );
  private native int native_numTraces( long ptr );
  private native int native_numSamples( long ptr );
  private native int native_numHeaders( long ptr );
  private native float native_sampleInt( long ptr );
  private native float native_hdrFloatValue( long ptr, int hdrIndex );
  private native double native_hdrDoubleValue( long ptr, int hdrIndex );
  private native int native_hdrIntValue( long ptr, int hdrIndex );

  private native boolean native_getNextTrace( long ptr, float[] samples, csSeismicTrace trace ) throws Exception;
  private native boolean native_moveToTrace( long ptr, int traceIndex, int numTracesToRead ) throws Exception;

  private native void native_setHeaderToPeek( long ptr, String headerName );
  private native void native_peekHeaderValue( long ptr, int traceIndex, csHeader value );

  private native String native_headerName( long ptr, int hdrIndex );
  private native String native_headerDesc( long ptr, int hdrIndex );
  private native int native_headerType( long ptr, int hdrIndex );

  private native int native_verticalDomain( long ptr );

  private native boolean native_setSelection(
          long ptr,
          String hdrValueSelectionText,
          String headerName,
          int sortOrder,
          int sortMethod,
          csITraceHeaderScanNotifier notifier );
  private native int native_getNumSelectedTraces( long ptr );
  private native void native_getSelectedValues( long ptr, csSelectedHeaderBundle hdrBundle );
  
  private native void native_freeInstance( long nativePtr );

  protected long myNativePtr = 0;

  public csNativeRSFReader( String filename ) throws IOException, Exception {
    try {
      myNativePtr = native_createInstance( filename );
    }
    catch( Exception e ) {
      throw(e);
    }
    if( myNativePtr == 0 ) {
      throw new IOException("Unknown error occured when opening RSF file " + filename);
    }
  }

  @Override
  public void closeFile() {
    native_closeFile( myNativePtr );
  }
  @Override
  public int numTraces() {
    return native_numTraces( myNativePtr );
  }
  @Override
  public int numSamples() {
    return native_numSamples( myNativePtr );
  }
  @Override
  public int numHeaders() {
    return native_numHeaders( myNativePtr );
  }
  @Override
  public float sampleInt() {
    return native_sampleInt( myNativePtr );
  }
  //-------------------------------------------------------------
  @Override
  public boolean getNextTrace( csSeismicTrace trace ) throws Exception {
    try {
      return native_getNextTrace( myNativePtr, trace.samples(), trace );
    }
    catch( Exception e ) {
//      System.out.println("Error occurred while getting next trace: " + e.getMessage());
      throw(e);
    }
  }
  @Override
  public boolean hasRandomFileAccess() {
    return true;
  }
  @Override
  public boolean moveToTrace( int traceIndex, int numTracesToRead ) throws Exception {
    return native_moveToTrace( myNativePtr, traceIndex, numTracesToRead );
  }
  //-------------------------------------------------------------
  @Override
  public int hdrIntValue( int hdrIndex ) {
    return native_hdrIntValue( myNativePtr, hdrIndex );
  }
  @Override
  public float hdrFloatValue( int hdrIndex ) {
    return native_hdrFloatValue( myNativePtr, hdrIndex );
  }
  @Override
  public double hdrDoubleValue( int hdrIndex ) {
    return native_hdrDoubleValue( myNativePtr, hdrIndex );
  }
  @Override
  public String headerName( int hdrIndex ) {
    String name = native_headerName( myNativePtr, hdrIndex );
    return name;
  }
  @Override
  public String headerDesc( int hdrIndex ) {
    String desc = native_headerDesc( myNativePtr, hdrIndex );
    return desc;
  }
  @Override
  public int headerType( int hdrIndex ) {
    return native_headerType( myNativePtr, hdrIndex );
  }
  @Override
  public boolean isFrequencyDomain() {
    return false;
  }
  @Override
  public int verticalDomain() {
    int jniDomain = native_verticalDomain( myNativePtr );
    int domain = csJNIDef.convertVerticalDomain_fromJNI( jniDomain );
    return domain;
  }
  @Override
  public void setHeaderToPeek( String headerName ) throws Exception {
//    throw( new Exception("Current implementation of RSF reader does not support header selection") );
    native_setHeaderToPeek( myNativePtr, headerName );
  }
  @Override
  public boolean peekHeaderValue( int traceIndex, csHeader value ) {
    try {
      native_peekHeaderValue( myNativePtr, traceIndex, value );
      return true;
    }
    catch( Exception e ) {
      System.out.println("Error: " + e.getMessage());
      return false;
    }
  }
  //
  @Override
  public boolean setSelection( String hdrValueSelectionText, String headerName, int sortOrder, int sortMethod,
          csITraceHeaderScanNotifier notifier ) {
    if( notifier == null ) return false;
    return native_setSelection( myNativePtr, hdrValueSelectionText, headerName, sortOrder, sortMethod, notifier );
  }
  @Override
  public int getNumSelectedTraces() {
    return native_getNumSelectedTraces( myNativePtr );
  }
  @Override
  public void getSelectedValues( csSelectedHeaderBundle hdrBundle ) {
    native_getSelectedValues( myNativePtr, hdrBundle );
  }
  //----------------------------------------------
  public void disposeNative() {
    if( myNativePtr != 0 ) {
      native_freeInstance( myNativePtr );
      myNativePtr = 0;
    }
  }
  public void finalize() {
    disposeNative();
  }
}



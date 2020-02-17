/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

import cseis.seis.csHeader;
import cseis.seis.csSeismicTrace;
import cseis.seis.csISeismicReader;

import java.io.IOException;

/**
 * JNI Interface to SEGY reader (C++ code)
 * @author Bjorn Olofsson
 */

public class csNativeSegyReader implements csISeismicReader {
  public static final int HDR_MAP_STANDARD = 1110;
  public static final int HDR_MAP_OBC      = 1111;
  public static final int HDR_MAP_SEND     = 1112;
  public static final int HDR_MAP_ARMSS    = 1113;
  public static final int HDR_MAP_PSEGY    = 1114;
  public static final int HDR_MAP_NODE1    = 1115;
  public static final int HDR_MAP_NODE2    = 1116;
  public static final int HDR_MAP_SU       = 1117;
  public static final int HDR_MAP_SU_ONLY  = 1118;
  public static final int HDR_MAP_SU_BOTH  = 1119;
  public static final int HDR_MAP_DAS      = 1120;
  
  public static final int NUM_DEFAULT_MAPS = 11;
  public static final int[] DEFAULT_MAPS = {
    HDR_MAP_STANDARD,
    HDR_MAP_OBC,
    HDR_MAP_SEND,
    HDR_MAP_ARMSS,
    HDR_MAP_PSEGY,
    HDR_MAP_NODE1,
    HDR_MAP_NODE2,
    HDR_MAP_SU,
    HDR_MAP_SU_ONLY,
    HDR_MAP_SU_BOTH,
    HDR_MAP_DAS
  };
  public static final String[] NAME_DEFAULT_MAPS = {
    "Standard",
    "OBC",
    "SEND",
    "ARMSS",
    "PSEGY",
    "NODE1",
    "NODE2",
    "SU (cseis)",
    "SU (orig)",
    "SU (cseis+orig)",
    "DAS"
  };

  private native long native_createInstance(
      String filename,
      int nTracesBuffer,
      int segyHeaderSet,
      boolean reverseByteOrderData,
      boolean reverseByteOrderHdr,
      boolean autoscale_hdrs );
  private native void native_closeFile( long ptr );
  private native int native_sampleByteSize( long ptr ) ;
  private native int native_numTraces( long ptr );
  private native int native_numSamples( long ptr );
  private native int native_numHeaders( long ptr );
  private native float native_sampleInt( long ptr );
  private native int native_hdrIntValue( long ptr, int hdrIndex );
  private native float native_hdrFloatValue( long ptr, int hdrIndex );
  private native double native_hdrDoubleValue( long ptr, int hdrIndex );

  private native boolean native_getNextTrace( long ptr, float[] samples, csSeismicTrace trace ) throws Exception;
  private native boolean native_moveToTrace( long ptr, int traceIndex, int numTracesToRead ) throws Exception;

  private native void native_setHeaderToPeek( long ptr, String headerName );
  private native void native_peekHeaderValue( long ptr, int traceIndex, csHeader value );

  private native void native_binHeader( long ptr, cseis.segy.csSegyBinHeader binHdr );
  private native String native_charHeader( long ptr );
  
  private native String native_headerName( long ptr, int hdrIndex );
  private native String native_headerDesc( long ptr, int hdrIndex );
  private native int native_headerType( long ptr, int hdrIndex );

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

  private static native void native_trcHeaderMap( int hdrMap, csSegyTrcHeaderDefinition hdrDef );

  protected long myNativePtr = 0;

  public csNativeSegyReader( String filename ) throws IOException, Exception {
    this( filename, HDR_MAP_STANDARD );
  }
  public csNativeSegyReader( String filename, int segyHeaderSet ) throws IOException, Exception {
    this( filename, segyHeaderSet, 200, false, false, true );
  }
  public csNativeSegyReader( String filename, int segyHeaderSet, int nTracesBuffer,
    boolean reverseByteOrderData, boolean reverseByteOrderHdr, boolean autoscale_hdrs ) throws IOException, Exception {
    try {
      myNativePtr = native_createInstance( filename, nTracesBuffer, segyHeaderSet, reverseByteOrderData, reverseByteOrderHdr, autoscale_hdrs );
    }
    catch( Exception e ) {
      throw(e);
    }
    if( myNativePtr == 0 ) {
      throw new IOException("Unknown error occured when opening SEGY file " + filename);
    }
  }

  @Override
  public void closeFile() {
    native_closeFile( myNativePtr );
  }
  public int sampleByteSize() {
    return native_sampleByteSize( myNativePtr );
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
  public String charHeader() {
    if( myNativePtr != 0 ) {
      return native_charHeader( myNativePtr );
    }
    else {
      return "";
    }
  }
  public cseis.segy.csSegyBinHeader binHeader() {
    cseis.segy.csSegyBinHeader binHdr = new cseis.segy.csSegyBinHeader();
    native_binHeader( myNativePtr, binHdr );
    return binHdr;
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
    return cseis.general.csUnits.DOMAIN_TIME;
  }
  @Override
  public void setHeaderToPeek( String headerName ) throws Exception {
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
  //-------------------------------------------------------------
  public static void trcHdrMap( int hdrMap, csSegyTrcHeaderDefinition hdrDef ) {
    native_trcHeaderMap( hdrMap, hdrDef );
  }

//  public boolean setSelection( String hdrValueSelectionText, String headerName, int sortOrder, int sortMethod ) {
//    return native_setSelection( myNativePtr, hdrValueSelectionText, headerName, sortOrder, sortMethod );
//  }  

  
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



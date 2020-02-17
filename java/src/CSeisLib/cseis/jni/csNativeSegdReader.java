/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

import cseis.seis.csHeader;
import cseis.seis.csSeismicTrace;
import cseis.seis.csISeismicReader;

import java.io.IOException;

/**
 * JNI Interface to SEGD reader (C++ code)
 * @author Bjorn Olofsson
 */

public class csNativeSegdReader implements csISeismicReader {
  /// Manufacturer/recording system IDs
  public static final int RECORDING_SYSTEM_SEAL    = 121;
  public static final int RECORDING_SYSTEM_GEORES  = 122;
  public static final int RECORDING_SYSTEM_GECO    = 123;

  /// External header/nav system IDs
  public static final int NAV_HEADER_NONE     = 200;
  public static final int NAV_HEADER_PSI      = 201;
  public static final int NAV_HEADER_LABO     = 202;
  public static final int NAV_HEADER_HYDRONAV_VER1 = 203;
  public static final int NAV_HEADER_HYDRONAV_VER6 = 204;

  public static final int  UNKNOWN                  = -1;

  /// Nav interface IDs
  public static final int  CM_DIGI_COMP_A = 101;
  public static final int  CM_DIGI_PSIB   = 102;
  public static final int  CM_DIGI_TS     = 103;
  public static final int  CM_DIGI_BIRD_DEPTHS = 104;

  public static final int  DUMP_NONE     = 0;
  public static final int  DUMP_GENERAL  = 2;
  public static final int  DUMP_CHANSET  = 4;
  public static final int  DUMP_EXTENDED = 8;
  public static final int  DUMP_EXTERNAL = 16;
  public static final int  DUMP_TRACE    = 32;
  public static final int  DUMP_EXTERNAL_BIRDS = 64;
  public static final int  DUMP_ALL      = DUMP_GENERAL + DUMP_CHANSET + DUMP_EXTENDED + DUMP_EXTERNAL + DUMP_TRACE;


  private native long native_createInstance(
      String filename,
      boolean readAuxTraces,
      boolean isDebug,
      boolean thisIsRev0,
      int navInterfaceID,
      int navSystemID );
  private native void native_closeFile( long ptr );
  private native int native_numTraces( long ptr );
  private native int native_numSamples( long ptr );
  private native int native_numHeaders( long ptr );
  private native float native_sampleInt( long ptr );
  private native int native_hdrIntValue( long ptr, int hdrIndex );
  private native float native_hdrFloatValue( long ptr, int hdrIndex );
  private native double native_hdrDoubleValue( long ptr, int hdrIndex );

  private native boolean native_getNextTrace( long ptr, float[] samples, csSeismicTrace trace ) throws Exception;
  private native boolean native_moveToTrace( long ptr, int traceIndex, int numTracesToRead ) throws Exception;

//  private native void native_setHeaderToPeek( long ptr, String headerName );
//  private native long native_peekHeaderValue( long ptr, int traceIndex );

  private native String native_headerName( long ptr, int hdrIndex );
  private native String native_headerDesc( long ptr, int hdrIndex );
  private native int native_headerType( long ptr, int hdrIndex );

  private native void native_freeInstance( long nativePtr );

  protected long myNativePtr = 0;

  public csNativeSegdReader( String filename ) throws IOException, Exception {
    boolean readAuxTraces = true;
    boolean isDebug    = false;
    boolean thisIsRev0 = false;
    int navInterfaceID = 0;
    int navSystemID    = csNativeSegdReader.NAV_HEADER_NONE;
    init( filename, readAuxTraces, isDebug, thisIsRev0, navInterfaceID, navSystemID );
  }
  private void init(
          String filename,
          boolean readAuxTraces,
          boolean isDebug,
          boolean thisIsRev0,
          int navInterfaceID,
          int navSystemID ) throws IOException, Exception {
    try {
      myNativePtr = native_createInstance( filename, readAuxTraces, isDebug, thisIsRev0, navInterfaceID, navSystemID );
    }
    catch( Exception e ) {
      throw(e);
    }
    if( myNativePtr == 0 ) {
      throw new IOException("Unknown error occured when opening SEG-D file " + filename);
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
    return false;
  }
  @Override
  public boolean moveToTrace( int traceIndex, int numTracesToRead ) throws Exception {
    throw( new Exception("Current implementation of SEG-D reader does not support move operation.") );
//    return native_moveToTrace( myNativePtr, traceIndex, numTracesToRead );
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
    throw( new Exception("Current implementation of SEG-D reader does not support header selection") );
    // Not implemented
    //    native_setHeaderToPeek( myNativePtr, headerName );
  }
  @Override
  public boolean peekHeaderValue( int traceIndex, csHeader value ) {
    return false;
//    try {
//      long retValue = native_peekHeaderValue( myNativePtr, traceIndex );
//      val.value = retValue;
//      return true;
//    }
//    catch( Exception e ) {
//      System.out.println("Error: " + e.getMessage());
//      return false;
//    }
  }
  //
  @Override
  public boolean setSelection( String hdrValueSelectionText, String headerName, int sortOrder, int sortMethod,
          csITraceHeaderScanNotifier notifier ) {
    return false;
  }
  @Override
  public int getNumSelectedTraces() {
    return 0;
  }
  @Override
  public void getSelectedValues( csSelectedHeaderBundle hdrBundle ) {
    return;
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



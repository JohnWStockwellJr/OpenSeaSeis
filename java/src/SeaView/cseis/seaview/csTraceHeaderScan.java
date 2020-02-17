/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.jni.csITraceHeaderScanNotifier;
import cseis.jni.csJNIDef;
import cseis.jni.csSeismicReaderUtils;
import cseis.jni.csSelectedHeaderBundle;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicReader;
import cseis.tools.csIStoppable;

/**
 * csTraceHeaderScan manages one trace header scanning operation.<br>
 * An object of this class is alive while a trace header scanning operation is ongoing. The scanning operation
 * is performed in a separate thread.
 * @author 2013 Felipe Punto
 */
public class csTraceHeaderScan implements csIStoppable {
  private csISeismicReader myReader;
  private String myFilename;
  private csITraceHeaderScanListener myListener;
  private Thread myThread;
  private csHeaderDef myHeaderDef;
  private int mySortOrder;
  private boolean myIsComplete;
  private boolean myIsRunning;
  private csITraceHeaderScanNotifier myNotifier;
  
  public csTraceHeaderScan( String filename, csHeaderDef headerDef, csITraceHeaderScanListener listener,
          csITraceHeaderScanNotifier notifier ) {
    myListener = listener;
    myNotifier = notifier;
    mySortOrder = csJNIDef.SORT_NONE;
    myHeaderDef = headerDef;
    myFilename = filename;
    myReader = null;
    myIsComplete = false;
    myIsRunning = false;
    myThread = null;
  }
  public boolean isRunning() {
    return myIsRunning;
  }
  public boolean start( csISeismicReader refReader ) {
    if( myIsRunning ) return false;
    try {
      myReader = csSeismicReaderUtils.createSeismicReader( refReader, myFilename );
      if( myReader == null ) throw new Exception("Trace header scanning is not supported for the current data format");
    }
    catch( Exception ex ) {
      myReader = null;
      myThread = null;
      return false;
    }
    myThread = new Thread( new Runnable() {
      @Override
      public void run() {
        myIsRunning = true;
        if( myReader == null ) return;
        boolean success = myReader.setSelection( "*", myHeaderDef.name, mySortOrder, csJNIDef.SIMPLE_SORT, myNotifier );
        myIsComplete = true;
        myIsRunning = false;
        myThread = null;
        if( success ) myListener.scanCompleted( myHeaderDef );
      }
    });
    myThread.start();
    return true;
  }
  public boolean isComplete() {
    return myIsComplete;
  }
  public csSelectedHeaderBundle retrieveSelectedValues() {
    int numTraces = myReader.getNumSelectedTraces();
    csSelectedHeaderBundle hdrBundle = new csSelectedHeaderBundle( numTraces );
    myReader.getSelectedValues( hdrBundle );
    return hdrBundle;
  }
  public void close() {
    myReader.closeFile();
    myThread = null;
  }
  @Override
  public void stopOperation() {
    if( myThread != null && myThread.isAlive() ) {
      myThread.interrupt();
      myThread = null;
    }
    myIsRunning = false;
  }
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seis;

import cseis.jni.csITraceHeaderScanNotifier;
import cseis.jni.csSelectedHeaderBundle;

/**
 * Interface to be implemented by seismic data readers.
 * @author Bjorn Olofsson
 */
public interface csISeismicReader {
  /**
   * @return Number of traces
   */
  public int numTraces();
  /**
   * @return Number of samples
   */
  public int numSamples();
  /**
   * @return Number of trace headers
   */
  public int numHeaders();
  /**
   * @return Sample interval in [ms], [m] or [Hz] depending on vertical domain
   */
  public float sampleInt();
  /**
   * @return Vertical domain
   */
  public int verticalDomain();
  /**
   * @return true if data is in frequency domain
   */
  public boolean isFrequencyDomain();

  /**
   * Retrieve next trace from input file
   * @param trace New seismic trace
   * @return true if another trace was read in
   */
  public boolean getNextTrace( csSeismicTrace trace ) throws Exception;

  /**
   * Move (file pointer) to specified trace.
   * The next call to getNextTrace() retrieves the trace with the specified trace index.
   * 
   * @param traceIndex (i) Trace index, starting at 0 for first trace
   * @param numTracesToRead (i) Number of traces to be read in at once (optimization to speed up reading)
   * @return true if operation was successful, false if error occurred (e.g. specified trace does not exist)
   */
  public boolean moveToTrace( int traceIndex, int numTracesToRead ) throws Exception;

  /**
   * @param hdrIndex Header index
   * @return Value of integer trace header
   */
  public int hdrIntValue( int hdrIndex );
  /**
   * @param hdrIndex Header index
   * @return Value of float trace header
   */
  public float hdrFloatValue( int hdrIndex );
  /**
   * @param hdrIndex Header index
   * @return Value of double trace header
   */
  public double hdrDoubleValue( int hdrIndex );
  /**
   * @param hdrIndex Header index
   * @return Name of trace header
   */
  public String headerName( int hdrIndex );
  /**
   * @param hdrIndex Header index
   * @return Description of trace header
   */
  public String headerDesc( int hdrIndex );
  /**
   * @param hdrIndex Header index
   * @return Type of trace header
   */
  public int headerType( int hdrIndex );

  public void setHeaderToPeek( String headerName ) throws Exception;
  public boolean peekHeaderValue( int traceIndex, csHeader value );

  /**
   * Prepare trace header selection
   * This should retrieve all trace header values of the given trace header which fulfil
   * the given selection criteria (first argument). The selected trace header values
   * can be retrieved by a subsequent call to getSelectedValues().
   * 
   * @param hdrValueSelectionText  String defining header value selection, following Seaseis' value selection syntax rules
   *                               Examples: "1600" "1600-1700" "1600-1700(2)" ">=1650"  "*" [=ALL]
   * @param headerName   Trace header name to use for selection
   * @param sortOrder    SORT_NONE: No sorting,  otherwise SORT_INCREASING or SORT_DECREASING
   * @param sortMethod   SIMPLE_SORT or TREE_SORT
   * @param notifier     Object which will be notified regularly during selection operation. Pass 'null' if not required
   * @return true if selection preparation was successful
   */
  public boolean setSelection(
          String hdrValueSelectionText,
          String headerName,
          int sortOrder,
          int sortMethod,
          csITraceHeaderScanNotifier notifier );
  /**
   * @return Number of selected traces (requires initial call to setSelection())
   */
  public int getNumSelectedTraces();
  /**
   * Retrieve selected trace values as part of the csSelectedHeaderBundle class
   * (requires initial call to setSelection())
   * @param hdrBundle Retrieved trace header values & more
   */
  public void getSelectedValues( csSelectedHeaderBundle hdrBundle );
  
  /**
   * Close file
   */
  public void closeFile();
  /**
   * 
   * @return true is random file access is supported by this reader
   */
  public boolean hasRandomFileAccess();
}



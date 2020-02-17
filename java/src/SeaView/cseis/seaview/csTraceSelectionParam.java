/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

/**
 * Bundles all parameters defining the trace selection method.<br>
 * <br>
 * Defines how and which traces are going to be selected for display
 * Possible selection modes: By sequential trace number, by ensemble, and by trace header value
 * 
 * @author 2013 Felipe Punto
 */
public class csTraceSelectionParam {
  public static final int NUM_TRACES_DEFAULT = 5000;
  public static int NUM_TRACES = NUM_TRACES_DEFAULT;
  
  /// Select traces by sequential trace number
  public static int SELECT_TRACE  = 0;
  /// Select traces by ensemble number
  public static int SELECT_ENS    = 1;
  /// Select traces by trace header value
  public static int SELECT_HEADER = 2;
  
  /// Index of first trace to display (starting at 0)
  public int firstTraceIndex;
  /// Number of traces to display
  public int numTraces;
  /// Trace increment. Only display traces with trace index = N x traceStep (N=0,1,2...)
  public int traceStep;
  
  /// Trace selection option: SELECT_TRACE, SELECT_ENS or SELECT_HEADER
  public int selectOption;
  
  /// For ensemble selection: Index of first ensemble to display (starting at 0)
  public int firstEnsIndex;
  /// For ensemble selection: Number of ensembles to display
  public int numEns;
  /// For ensemble selection: ???
//  public int ensHdrSelectedIndex;

  /// For header selection: Value of selected header to display
  double selectedHdrValue;
  /// For header selection: Number of traces to display
  public int numTracesHdr;
  
  // For ensemble & header selection: Selected header name
  public String selectedHdrName;
  // For ensemble & header selection: Index of selected header
  public int selectedHdrIndex;
  
  public csTraceSelectionParam() {
    selectOption    = SELECT_TRACE;
    firstTraceIndex = 0;
    numTraces       = NUM_TRACES;
    traceStep       = 1;
    firstEnsIndex   = 0;
    numEns          = 1;
    selectedHdrName = "";

    selectedHdrIndex = -1;
    selectedHdrValue = 0;
    numTracesHdr     = NUM_TRACES;
  }
  public csTraceSelectionParam( csTraceSelectionParam sp ) {
    selectOption    = sp.selectOption;
    firstTraceIndex = sp.firstTraceIndex;
    numTraces       = sp.numTraces;
    traceStep       = sp.traceStep;
    firstEnsIndex   = sp.firstEnsIndex;
    selectedHdrName = sp.selectedHdrName;
    numEns          = sp.numEns;
    selectedHdrIndex = sp.selectedHdrIndex;
    selectedHdrValue   = sp.selectedHdrValue;
    numTracesHdr       = sp.numTracesHdr;
  }
  public void dump() {
    System.out.println("-----------------------------------");
    System.out.println("Dump of csTraceSelectionParam");
    System.out.println("selectOption:     " + selectOption);
    System.out.println("firstTraceIndex:  " + firstTraceIndex);
    System.out.println("numTraces:        " + numTraces);
    System.out.println("traceStep:        " + traceStep);

    System.out.println("firstEnsIndex:    " + firstEnsIndex);
    System.out.println("numEns:           " + numEns);
    System.out.println("selectedHdrName:  " + selectedHdrName);
    System.out.println("selectedHdrValue: " + selectedHdrValue);
    System.out.println("selectedHdrIndex: " + selectedHdrIndex);
    System.out.println("numTracesHdr:     " + numTracesHdr);
  }
}



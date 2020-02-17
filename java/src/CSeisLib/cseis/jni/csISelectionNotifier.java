/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

/**
 * Listener interface for 'selection' notifications across JNI interface.<br>
 * Notifications allow C++ code to communicate which trace is currently being scanned
 * to Java code. 'Selection' functionality means trace selection based on certain trace header
 * values. To enable this functionality, all traces need to be scanned.
 * @author 2013 Felipe Punto
 */
public interface csISelectionNotifier {
  /**
   * Notify listener about the trace which is currently being worked on (=scanned).
   * @param traceIndex Current trace index
   */
  public void notify( int traceIndex );
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.jni.csSelectedHeaderBundle;
import cseis.seis.csHeaderDef;
import java.util.ArrayList;

/**
 * Class defining parameters for scanned trace headers.
 * @author 2013 Felipe Punto
 */
public class csScanHeaderInfo {
  private csHeaderDef myHeaderDef;
  private ArrayList<Number> ensValueList;
  private ArrayList<Integer> ensFirstTraceIndexList;
  private Object[] mySortedHdrValues;
  private int[] mySortedHdrFirstTraceIndex;
  private int[] mySortedHdrLastTraceIndex;
  private csSelectedHeaderBundle mySelectedHeaderBundle;
  private int myTotalNumTraces;

  /// Search for traces
  private Number myCurrentHdrValue;
  private int myCurrentTraceIndex;
  private int myCurrentSortedHdrIndex;
  private boolean myCurrentMoveForwards;
  
  public csScanHeaderInfo() {
    myHeaderDef = null;
    mySortedHdrValues = null;
    ensValueList           = new ArrayList<Number>();
    ensFirstTraceIndexList = new ArrayList<Integer>();
    mySelectedHeaderBundle = null;
    myTotalNumTraces = 0;
    myCurrentHdrValue = null;
    myCurrentTraceIndex = -1;
    mySortedHdrFirstTraceIndex = null;
    mySortedHdrLastTraceIndex = null;
    myCurrentSortedHdrIndex = -1;
  }
  public csScanHeaderInfo( csHeaderDef headerDef, csSelectedHeaderBundle bundle ) {
    this();
    myHeaderDef = headerDef;
    mySelectedHeaderBundle = bundle;
    myTotalNumTraces = mySelectedHeaderBundle.hdrValues.length;
    ArrayList<Number> uniqueValueList = new ArrayList<Number>();
    int numTraces = bundle.hdrValues.length;
    if( numTraces == 0 ) {
      return;
    }
    Number valuePrev = (Number)bundle.hdrValues[0].value();
    ensValueList.add( valuePrev );
    ensFirstTraceIndexList.add( bundle.traceIndexList[0] );
    uniqueValueList.add( valuePrev );
    for( int i = 1; i < numTraces; i++ ) {
      Number value = (Number)bundle.hdrValues[i].value();
      if( !uniqueValueList.contains( value ) ) {
        uniqueValueList.add( value );
      }
      if( !value.equals(valuePrev) ) {
        ensValueList.add( value );
        ensFirstTraceIndexList.add( bundle.traceIndexList[i] );
        valuePrev = value;
      }
    }
    mySortedHdrValues = uniqueValueList.toArray();
    java.util.Arrays.sort( mySortedHdrValues );
    mySortedHdrFirstTraceIndex = new int[mySortedHdrValues.length];
    mySortedHdrLastTraceIndex = new int[mySortedHdrValues.length];
    for( int ihdr = 0; ihdr < mySortedHdrValues.length; ihdr++ ) {
      Number value = (Number)mySortedHdrValues[ihdr];
      for( int i = 0; i < numTraces; i++ ) {
        if( value.equals(mySelectedHeaderBundle.hdrValues[i].value()) ) {
          mySortedHdrFirstTraceIndex[ihdr] = i;
          break;
        }
      }
      for( int i = numTraces-1; i >= 0; i-- ) {
        if( value.equals(mySelectedHeaderBundle.hdrValues[i].value()) ) {
          mySortedHdrLastTraceIndex[ihdr] = i;
          break;
        }
      }
    }
  }
  public String getHeaderName() {
    return myHeaderDef.name;
  }
  public int getHeaderType() {
    return myHeaderDef.type;
  }
  public int getNumEnsembles() {
    return ensValueList.size();
  }
  public Number getEnsValue( int ensIndex ) {
    return ensValueList.get(ensIndex);
  }
  public int getEnsFirstTraceIndex( int ensIndex ) {
    return ensFirstTraceIndexList.get(ensIndex).intValue();
  }

  public int getNumSortedHeaderValues() {
    return mySortedHdrValues.length;
  }
  public Number getSortedHeaderValue( int sortedHdrIndex ) {
    return (Number)mySortedHdrValues[sortedHdrIndex];
  }
  
  public void reset_getNextHeaderTraceIndex( double hdrValue, int traceIndex1, boolean moveForwards ) throws Exception {
    Number numberValue = null;
    myCurrentMoveForwards = moveForwards;
    switch( myHeaderDef.type ) {
      case csHeaderDef.TYPE_INT:
        numberValue = (int)hdrValue;
        break;
      case csHeaderDef.TYPE_DOUBLE:
        numberValue = hdrValue;
        break;
      case csHeaderDef.TYPE_FLOAT:
        numberValue = hdrValue;
        break;
      case csHeaderDef.TYPE_LONG:
        numberValue = (long)hdrValue;
        break;
      default:
        numberValue = (int)hdrValue;
        break;
    }
    int sortedHdrIndex = java.util.Arrays.binarySearch(mySortedHdrValues, numberValue );
    if( sortedHdrIndex >= 0 ) {
      myCurrentHdrValue = (Number)mySortedHdrValues[sortedHdrIndex];
      myCurrentTraceIndex = traceIndex1;
      myCurrentSortedHdrIndex = sortedHdrIndex;
    }
    else {
      throw new Exception("Trace header value " + numberValue.toString() + " not found.");
    }
//    for( int i = 0; i < mySortedHdrValues.length; i++ ) {
//      if( mySortedHdrValues[i].equals(hdrValue) ) {
//        return;
//      }
//    }
//    throw( new Exception("non-existent header value provided") );
  }
  public int getCurrentHdrFirstTraceIndex() {
    return mySortedHdrFirstTraceIndex[myCurrentSortedHdrIndex];
  }
  public int getCurrentHdrLastTraceIndex() {
    return mySortedHdrLastTraceIndex[myCurrentSortedHdrIndex];
  }
  public int getNextHeaderTraceIndex() {
    if( myCurrentMoveForwards ) {
      while( myCurrentTraceIndex < myTotalNumTraces &&
             !mySelectedHeaderBundle.hdrValues[myCurrentTraceIndex].value().equals(myCurrentHdrValue) ) {
        myCurrentTraceIndex += 1;
      }
      if( myCurrentTraceIndex < myTotalNumTraces ) {
        return myCurrentTraceIndex++;
      }
    }
    else {
      while( myCurrentTraceIndex >= 0 &&
             !mySelectedHeaderBundle.hdrValues[myCurrentTraceIndex].value().equals(myCurrentHdrValue) ) {
        myCurrentTraceIndex -= 1;
      }
      if( myCurrentTraceIndex >= 0 ) {
        return myCurrentTraceIndex--;
      }
    }
    return -1;
  }
}


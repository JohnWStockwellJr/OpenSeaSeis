/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import java.util.ArrayList;

/**
 * Class containing time picks for one seismic horizon.
 * Horizon key is the consecutive trace number.
 * Maximum one pick per trace.
 */
public class csHorizon {
  /// Pick 'arrays' - Each pick 'array' contains an arbitrary number of picks made on a contiguous set of traces
  /// Each pick (on a particular trace) is unique and is only contained in one array
  private final ArrayList<csHorizonPickArray> myPickArrayList;
  /// Left-most trace which has been picked in this horizon
  private int myMinTrace;
  /// Right-most trace which has been picked in this horizon
  private int myMaxTrace;
  private int myIteratorIndex;

  csHorizon() {
    myPickArrayList = new ArrayList<csHorizonPickArray>();
    myIteratorIndex = 0;
    myMinTrace = Integer.MAX_VALUE;
    myMaxTrace = 0;
  }
  public void startIteration() {
    myIteratorIndex = 0;
  }
  public boolean hasNext() {
    return( myIteratorIndex < myPickArrayList.size() );
  }
  public Object next() {
    if( myIteratorIndex >= myPickArrayList.size() ) return null;
    return myPickArrayList.get(myIteratorIndex++);
  }
  public csHorizonPickArray getPickArray( int trace1Out, int trace2Out ) {
    csHorizonPickArray pickArrayOut = new csHorizonPickArray( trace1Out, trace2Out );
    if( pickArrayOut.numTraces == 0 ) return null;
    for( csHorizonPickArray pickArrayIn : myPickArrayList ) {
      int trace1In = pickArrayIn.trace1;
      int trace2In = pickArrayIn.trace2;
      if( trace2In < trace1Out ) continue;
      if( trace1In > trace2Out ) break;
      int traceFirst = Math.max( trace1In, trace1Out );
      int traceLast = Math.min( trace2In, trace2Out );
      for( int itrc = traceFirst; itrc <= traceLast; itrc++ ) {
        pickArrayOut.picks[itrc-trace1Out] = pickArrayIn.picks[itrc-trace1In];
      }
    }
    return pickArrayOut;
  }
  public void updatePickArray( int trace1New, int trace2New, float[] pickTimesNew ) {
    updatePickArray( new csHorizonPickArray( trace1New, trace2New, pickTimesNew ) );
  }
  public void updatePickArray( csHorizonPickArray pickArrayNew ) {
    if( pickArrayNew.trace1 < myMinTrace ) myMinTrace = pickArrayNew.trace1;
    if( pickArrayNew.trace2 > myMaxTrace ) myMaxTrace = pickArrayNew.trace2;
    for( int iarray = 0; iarray < myPickArrayList.size(); iarray++ ) {
      csHorizonPickArray pickArrayCurrent = myPickArrayList.get(iarray);
      // a) New array is fully on the right-hand side from the current array
      if( pickArrayCurrent.trace2 < pickArrayNew.trace1 ) continue;
      // b) New array is fully on the left-hand side from the current array --> Add to array list without further changes
      if( pickArrayCurrent.trace1 > pickArrayNew.trace2 ) {
        myPickArrayList.add( iarray, pickArrayNew ); // Add new array at current location into array list
        break;
      }
      // --> New array starts at or above start of current array
      csHorizonPickArray leftArray = pickArrayCurrent;
      int leftArrayIndex  = iarray;
      int rightArrayIndex = iarray;
      while( pickArrayCurrent.trace2 < pickArrayNew.trace2 ) {
        iarray += 1;
        if( iarray == myPickArrayList.size() ) break;
        rightArrayIndex = iarray;
        pickArrayCurrent = myPickArrayList.get(iarray); // Replace current array with the next in list
      }
      csHorizonPickArray rightArray = pickArrayCurrent;
      // c) Remove all arrays from the array list which partly (or entirely) overlap the new array
      iarray = leftArrayIndex;
      for( int i = leftArrayIndex; i <= rightArrayIndex; i++ ) {
        myPickArrayList.remove( iarray );
      }
      // d) If there are remaining picks in left-most array, copy them to new subarray
      if( leftArray.trace1 < pickArrayNew.trace1 ) {
        csHorizonPickArray subArrayLeft = leftArray.createSubArray(leftArray.trace1, pickArrayNew.trace1-1 );
        myPickArrayList.add( iarray++, subArrayLeft );
      }
      myPickArrayList.add( iarray++, pickArrayNew );
      // e) If there are remaining picks in right-most array, copy them to new subarray
      if( rightArray.trace2 > pickArrayNew.trace2 ) {
        csHorizonPickArray subArrayRight = rightArray.createSubArray(pickArrayNew.trace2+1, rightArray.trace2 );
        myPickArrayList.add( iarray++, subArrayRight );
      }
      return;
    }
    myPickArrayList.add( pickArrayNew );
  }
  public void resetAllPicks( ArrayList<Integer> traceNumList, ArrayList<Float> pickList ) {
    myPickArrayList.clear();
    int numTraces = traceNumList.size();
    int minTrace  = traceNumList.get(0);
    int maxTrace = traceNumList.get(numTraces-1);
    csHorizonPickArray array = new csHorizonPickArray( minTrace, maxTrace );
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      int traceNum = traceNumList.get(itrc);
      float pick = pickList.get(itrc);
      array.picks[traceNum-minTrace] = pick;
    }
    myPickArrayList.add( array );
  }
  public void dump() {
    for( csHorizonPickArray a : myPickArrayList ) {
      System.out.println("" + a.trace1 + " " + a.picks[0] + "\n" + a.trace2 + " " + a.picks[a.numTraces-1] + " \n" );
    }
  }

}

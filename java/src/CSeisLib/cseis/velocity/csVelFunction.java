/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.velocity;

import java.util.ArrayList;

/**
 * Velocity picks for one 1D/2D/3D ensemble/gather.
 */
public class csVelFunction {
  public static final int MIN_TIME_DISTANCE = 50; // [ms]
  private final csEnsemble myEns;
  private final ArrayList<csVelFunction.Item> myPicks;
  /**
   * 
   * @param idVal1 Unique identifier of picked gather, e.g. CMP
   */
  public csVelFunction( int idVal1 ) {
    myEns = new csEnsemble( idVal1 );
    myPicks = new ArrayList<csVelFunction.Item>();
  }
  /**
   * 
   * @param idVal1 Unique identifier of picked gather (together with idVal2), e.g. INLINE
   * @param idVal2 Unique identifier of picked gather (together with idVal1), e.g. XLINE 
   */
  public csVelFunction( int idVal1, int idVal2 ) {
    myEns = new csEnsemble( idVal1, idVal2 );
    myPicks = new ArrayList<csVelFunction.Item>();
  }
  /**
   * 
   * @param ens Unique ensemble = identifier of picked gather
   */
  public csVelFunction( csEnsemble ens ) {
    myEns = ens;
    myPicks = new ArrayList<csVelFunction.Item>();
  }
  public int getIDVal1() {
    return myEns.id1();
  }
  public int getIDVal2() {
    return myEns.id2();
  }
  public csEnsemble getEns() {
    return myEns;
  }
  public int numPicks() {
    return myPicks.size();
  }
  public csVelFunction.Item getPick( int index ) {
    return myPicks.get( index );
  }
  public java.util.List<csVelFunction.Item> getPicks() {
    return myPicks;
  }
  public void addPick( float time, float vel ) {
    Item velPickNew = new Item(time,vel);
    if( myPicks.isEmpty() ) {
      myPicks.add( velPickNew );
      return;
    }
    
    int index = getNearestPickIndex(time);
    float dist;
    if( index < myPicks.size() ) {
      dist = myPicks.get(index).time - time;
      if( Math.abs(dist) < MIN_TIME_DISTANCE ) {
        myPicks.set( index, velPickNew );
      }
      else if( dist > 0 ) {
        myPicks.add( index, velPickNew );
      }
      else {
        myPicks.add( index+1, velPickNew );
      }
    }
    else {
      dist = time - myPicks.get(myPicks.size()-1).time;
      if( dist < MIN_TIME_DISTANCE ) {
        myPicks.set( index-1, velPickNew );
      }
      else {
        myPicks.add( velPickNew );
      }
    }
  }
  public boolean removePick( float time ) {
    if( myPicks.isEmpty() ) return false;
    int index = getNearestPickIndex( time );
    if( index < myPicks.size() ) {
      myPicks.remove(index);
    }
    else {
      myPicks.remove(myPicks.size()-1);
    }
    return true;
  }
  public void clearPicks() {
    myPicks.clear();
  }
  private int getNearestPickIndex( float time ) {
    for( int i = 0; i < myPicks.size(); i++ ) {
      Item p = myPicks.get(i);
      if( p.time >= time ) {
        if( i == 0 ) return i;
        float dist1 = p.time - time;
        float dist2 = time - myPicks.get(i-1).time;
        if( dist1 < dist2 ) return i;
        else return i-1;
      }
    }
    return myPicks.size();
  }
  public String toString() {
    String text = myEns.toString() + "\n";
    for( Item item : myPicks ) {
      text += " " + item.time + " " + item.vel + "\n";
    }
    return text;
  }
  //-----------------------------------------------------------
  public class Item {
    public Item( float t, float v ) {
      time = t;
      vel = v;
    }
    /// Time in milliseconds
    public float time;
    /// Velocity in meters per second (or any other unit)
    public float vel;
  }
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.velocity;

import java.util.Objects;

public class csEnsemble {
  public static final int DIM_1D = 11;
  public static final int DIM_2D = 22;
  public static final int DIM_3D = 33;

  private Integer myId1;
  private Integer myId2;
  private int myDim;

  public csEnsemble() {
    this( -1, -1 );
    myDim = DIM_1D;
  }
  public csEnsemble( int id1_in ) {
    this( id1_in, -1 );
    myDim = DIM_2D;
  }
  public csEnsemble( int id1_in, int id2_in ) {
    myId1 = id1_in;
    myId2 = id2_in;
    myDim = DIM_3D;
  }
  public boolean is1D() {
    return( myDim == csEnsemble.DIM_1D );
  }
  public boolean is2D() {
    return( myDim == csEnsemble.DIM_2D );
  }
  public boolean is3D() {
    return( myDim == csEnsemble.DIM_3D );
  }
  public int getDim() {
    return( myDim );
  }
  public int id1() {
    return myId1;
  }
  public int id2() {
    return myId2;
  }
  @Override
  public boolean equals( Object o ) {
    csEnsemble ens2 = (csEnsemble)o;
    if( is1D() || ens2.is1D() ) return true;
    if( is2D() || ens2.is2D() ) return ( Objects.equals(myId1, ens2.myId1) );
    return( Objects.equals(myId1, ens2.myId1) && Objects.equals(myId2, ens2.myId2) );
  }
  public int compareTo( csEnsemble ens2 ) {
    if( is1D() || ens2.is1D() ) return 0;
    if( is2D() || ens2.is2D() || !Objects.equals(myId1, ens2.myId1) ) return( myId1 - ens2.myId1 );
    return( myId2 - ens2.myId2 );
  }
  @Override
  public String toString() {
    if( myDim == DIM_2D ) return "" + myId1;
    if( myDim == DIM_3D ) return "" + myId1 + " "  + myId2;
    return "";
  }
  public String dimToString() {
    if( myDim == DIM_1D ) return "1D";
    else if( myDim == DIM_2D ) return "2D";
    else return "3D";
  }
}

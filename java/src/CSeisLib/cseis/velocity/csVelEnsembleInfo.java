/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.velocity;

public class csVelEnsembleInfo {
  public int hdrIndexVel;
  public int hdrIndexEns1;
  public int hdrIndexEns2;
  public String hdrNameVel;
  public String hdrNameEns1;
  public String hdrNameEns2;
  private final int myDim;
  public csVelEnsembleInfo( String hdrNameVel_in ) {
    hdrNameVel  = hdrNameVel_in;
    myDim = csEnsemble.DIM_1D;
  }
  public csVelEnsembleInfo( String hdrNameVel_in, String hdrNameEns1_in ) {
    hdrNameVel   = hdrNameVel_in;
    hdrNameEns1  = hdrNameEns1_in;
    myDim = csEnsemble.DIM_2D;
  }
  public csVelEnsembleInfo( String hdrNameVel_in, String hdrNameEns1_in, String hdrNameEns2_in ) {
    hdrNameVel   = hdrNameVel_in;
    hdrNameEns1  = hdrNameEns1_in;
    hdrNameEns2  = hdrNameEns2_in;
    myDim = csEnsemble.DIM_3D;
  }
  public int getDim() {
    return myDim;
  }
  public void dump() {
    System.out.println("Ensemble info:");
    System.out.println("Velocity: " + hdrNameVel + " " + hdrIndexVel);
    if( myDim != csEnsemble.DIM_1D ) System.out.println("Ensemble1: " + hdrNameEns1 + " " + hdrIndexEns1);
    if( myDim == csEnsemble.DIM_3D ) System.out.println("Ensemble2: " + hdrNameEns2 + " " + hdrIndexEns2);
    System.out.println("Dimension: " + myDim );
  }
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.velocity;

public class csVelPickBundle {
  csEnsemble ens;
  int indexOfTrace1;
  float velOfTrace1;
  float velInc;
  csVelFunction velFunc;

  public csVelPickBundle( csEnsemble e ) {
    ens = e;
  }
  public csVelPickBundle( csEnsemble e, int trace, float vel, float velInc_in ) {
    ens = e;
    indexOfTrace1 = trace;
    velOfTrace1 = vel;
    velInc = velInc_in;
    velFunc = new csVelFunction( ens );
  }
  public float traceToVel( float traceIndex ) {
    return( velOfTrace1 + velInc * (traceIndex-indexOfTrace1) );
  }
  @Override
  public String toString() {
    return "Ensemble keys: " + ens.toString() + ", vel/index/velinc: " + velOfTrace1 + "/" + indexOfTrace1 + "/" + velInc;
  }
  @Override
  public boolean equals( Object o ) {
    if( o instanceof csEnsemble ) {
      return( ens.equals( o ) );
    }
    return( ens.equals( ((csVelPickBundle)o).ens ) );
  }
  
}

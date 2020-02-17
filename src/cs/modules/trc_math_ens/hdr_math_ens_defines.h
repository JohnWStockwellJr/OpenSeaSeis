/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef HDR_MATH_ENS_DEFINES_H
#define HDR_MATH_ENS_DEFINES_H

#include <cstdio>

using namespace std;

namespace mod_trc_math_ens {
  class Ens {
  public:
    Ens( int numEns ) {
      myNumEnsMax = numEns;
      myNumTraces = new int[myNumEnsMax];
      myEnsValues = new double[myNumEnsMax];
      myNumEnsCurrent = 0;
    }
    ~Ens() { delete [] myNumTraces; delete [] myEnsValues; }
    int numEns() const { return myNumEnsCurrent; }
    int numTraces( int ens ) const { return myNumTraces[ens]; }
    int numAllTraces() const {
      int ntraces = 0;
      for( int i = 0; i < myNumEnsCurrent; i++ ) {
        ntraces += myNumTraces[i];
      }
      return ntraces;
    }
    double ensValue( int ens ) const { return myEnsValues[ens]; }
    int addEns( int numTraces, double ensValue ) {
      //      fprintf(stdout,"Ensemble: %d (max: %d)  %f\n", myNumEnsCurrent, myNumEnsMax, ensValue);
      if( myNumEnsCurrent == myNumEnsMax ) {
        fprintf(stdout,"Incorrect ensemble: %d (max: %d)  %f\n", myNumEnsCurrent, myNumEnsMax, ensValue);
        exit(-1);
      }
      myNumTraces[myNumEnsCurrent] = numTraces;
      myEnsValues[myNumEnsCurrent] = ensValue;
      myNumEnsCurrent += 1;
      return myNumEnsCurrent;
    }
    int releaseFirstEns() {
      if( myNumEnsCurrent == 0 ) return 0;
      int ntraces = myNumTraces[0];
      for( int i = 1; i < myNumEnsCurrent; i++ ) {
        myNumTraces[i-1] = myNumTraces[i];
        myEnsValues[i-1] = myEnsValues[i];
      }
      myNumEnsCurrent -= 1;
      return ntraces;
    }
  private:
    int* myNumTraces;
    double* myEnsValues;
    int myNumEnsMax;
    int myNumEnsCurrent;
  };

  class Dim {
  public:
    Dim() {}
    ~Dim() {}
    int num;
    int hdrId;
    double valCurrent;
    int numCurrent;
  };

} // END namespace
#endif

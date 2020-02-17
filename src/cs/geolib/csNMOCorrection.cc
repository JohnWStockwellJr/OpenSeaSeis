/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <cstring>
#include <iostream>
#include "csNMOCorrection.h"
#include "csInterpolation.h"
#include "csTimeFunction.h"
#include "csException.h"

using namespace cseis_geolib;

namespace cseis_geolib {
  float getQuadAmplitudeAtSample( float const* traceData, double sample, int numSamples );
}

csNMOCorrection::csNMOCorrection( double sampleInt_ms, int numSamples, int method_nmo ) {
  myNMOMethod          = method_nmo;
  mySampleInt_sec      = sampleInt_ms/1000.0;
  myNumSamples         = numSamples;
  myTimeOfLastSample   = (myNumSamples-1)*mySampleInt_sec;

  myBuffer = new float[myNumSamples];
  myModeOfApplication = csNMOCorrection::NMO_APPLY;
  myOffsetApex = 4000;
  myIsHorizonBasedNMO = false;
  myVelocityTrace = NULL;
  myETATrace      = NULL;
  myTimeTrace     = NULL;
  myTimeTraceInverse = NULL;
  myTimeTraceDiff    = NULL;
  myInterpol = NULL;
  myTimeSample1_s = 0.0;

  allocateVelocity();
}
//--------------------------------------------------------------------------------
//
csNMOCorrection::~csNMOCorrection() {
  if( myBuffer != NULL ) {
    delete [] myBuffer;
    myBuffer = NULL;
  }
  if( myVelocityTrace != NULL ) {
    delete [] myVelocityTrace;
    myVelocityTrace = NULL;
  }
  if( myETATrace != NULL ) {
    delete [] myETATrace;
    myETATrace = NULL;
  }
  if( myTimeTrace != NULL ) {
    delete [] myTimeTrace;
    myTimeTrace = NULL;
  }
  if( myTimeTraceInverse != NULL ) {
    delete [] myTimeTraceInverse;
    myTimeTraceInverse = NULL;
  }
  if( myTimeTraceDiff != NULL ) {
    delete [] myTimeTraceDiff;
    myTimeTraceDiff = NULL;
  }
  if( myInterpol != NULL ) {
    delete myInterpol;
    myInterpol = NULL;
  }
}
//--------------------------------------------------------------------------------
//
void csNMOCorrection::setTimeSample1( float timeSample1_ms ) {
  if( timeSample1_ms > 0 ) throw( csException("Inconsistent time of first sample provided: %f. Must be <= 0\n", timeSample1_ms) );
  myTimeSample1_s = timeSample1_ms / 1000.0f;
}
//--------------------------------------------------------------------------------
//
void csNMOCorrection::setHorizonBasedNMO( bool setHorizonBased, int method_interpolation ) {
  if( myIsHorizonBasedNMO == setHorizonBased ) return;
  myIsHorizonBasedNMO = setHorizonBased;
  myHorInterpolationMethod = method_interpolation;
  allocateVelocity();
}
//--------------------------------------------------------------------------------
//
void csNMOCorrection::allocateVelocity() {
  if( myVelocityTrace != NULL ) {
    delete [] myVelocityTrace;
    myVelocityTrace = NULL;
  }
  if( myETATrace != NULL ) {
    delete [] myETATrace;
    myETATrace = NULL;
  }
  if( myTimeTrace != NULL ) {
    delete [] myTimeTrace;
    myTimeTrace = NULL;
  }
  if( myIsHorizonBasedNMO ) {
    myVelocityTrace = NULL;
  }
  else {
    myVelocityTrace = new float[myNumSamples];
    if( myNMOMethod == csNMOCorrection::PP_NMO_VTI ) {
      myETATrace = new float[myNumSamples];
    }
    myTimeTrace     = new float[myNumSamples];
    // Probably unnecessary to set time trace here...will be overwritten anyway
    for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
      myTimeTrace[isamp] = (float)( mySampleInt_sec * isamp ) + myTimeSample1_s;
    }
    myInterpol = new csInterpolation( myNumSamples, mySampleInt_sec );
  }
}
//--------------------------------------------------------------------------------
//
void csNMOCorrection::setModeOfApplication( int mode ) {
  myModeOfApplication = mode;
}
//--------------------------------------------------------------------------------
//
void csNMOCorrection::setEmpiricalNMO( double offsetApex_m, double zeroOffsetDamping ) {
  myNMOMethod  = csNMOCorrection::EMPIRICAL_NMO;
  myOffsetApex = offsetApex_m;
  myZeroOffsetDamping = zeroOffsetDamping;
}
//--------------------------------------------------------------------------------
//
void csNMOCorrection::perform_nmo( float const* samplesIn, int numVelocities_in, float const* time_in, float const* vel_rms_in, double offset, float* samplesOut ) {
  memcpy(samplesOut, samplesIn, myNumSamples*sizeof(float) );
  perform_nmo( numVelocities_in, time_in, vel_rms_in, offset, samplesOut );
}

void csNMOCorrection::perform_nmo( int numVelocities_in, float const* time_in, float const* vel_rms_in, double offset, float* samplesOut ) {
  if( myIsHorizonBasedNMO ) {
    return( perform_nmo_horizonBased_prepare( numVelocities_in, time_in, vel_rms_in, offset, samplesOut ) );
  }
  else {
    return( perform_nmo_internal( numVelocities_in, time_in, vel_rms_in, offset, samplesOut, NULL ) );
  }
}
//--------------------------------------------------------------------
//
//
void csNMOCorrection::perform_nmo( csTimeFunction<double> const* velTimeFunc, double offset, float* samplesOut ) {
  int numVelocities_in = velTimeFunc->numValues();
  // Horizon based: Add one vel location at top and bottom of trace (temp fix to make sure that interpolation works)
  int numVelocities    = myIsHorizonBasedNMO ? numVelocities_in+2 : numVelocities_in;

  float* t0      = new float[numVelocities];
  float* vel_rms = new float[numVelocities];

  if( myIsHorizonBasedNMO ) {
    t0[0]      = 0;
    vel_rms[0] = (float)velTimeFunc->valueAtIndex(0);
    for( int i = 0; i < numVelocities_in; i++ ) {
      t0[i+1]      = (float)(velTimeFunc->timeAtIndex(i)/1000.0);  // Convert to seconds
      vel_rms[i+1] = (float)velTimeFunc->valueAtIndex(i);
    }
    t0[numVelocities-1]      = myTimeOfLastSample;
    vel_rms[numVelocities-1] = vel_rms[numVelocities-2];
    
    perform_nmo_horizonBased( numVelocities, t0, vel_rms, offset, samplesOut );
  }
  else {
    for( int i = 0; i < numVelocities_in; i++ ) {
      t0[i]      = (float)velTimeFunc->timeAtIndex(i)/1000.0;  // Convert to seconds
      vel_rms[i] = (float)velTimeFunc->valueAtIndex(i);
    }
    if( myNMOMethod != csNMOCorrection::PP_NMO_VTI ) {
      perform_nmo_internal( numVelocities_in, t0, vel_rms, offset, samplesOut, NULL );
    }
    else {
      if( velTimeFunc->numSpatialValues() < 2 ) throw( csException("csNMOCorrection::perform_nmo: VTI velocity function is lacking an eta function") );
      float* eta = new float[numVelocities];
      for( int i = 0; i < numVelocities_in; i++ ) {
        eta[i] = velTimeFunc->valueAtIndex(i,1);
      }
      perform_nmo_internal( numVelocities_in, t0, vel_rms, offset, samplesOut, eta );
      delete [] eta;
    }
  }
  delete [] t0;
  delete [] vel_rms;
}
//--------------------------------------------------------------------------------
//
void csNMOCorrection::perform_nmo_internal( int numVels, float const* times, float const* vel_rms, double offset, float* samplesOut, float const* eta ) {

  // Linearly interpolate input velocities
  csInterpolation::linearInterpolation( numVels, times, vel_rms, myNumSamples, mySampleInt_sec, myVelocityTrace );
  if( myNMOMethod == csNMOCorrection::PP_NMO_VTI ) {
    csInterpolation::linearInterpolation( numVels, times, eta, myNumSamples, mySampleInt_sec, myETATrace );
  }

  memcpy( myBuffer, samplesOut, myNumSamples*sizeof(float) );
  double offset_sq = offset*offset;
  double tmp;
  int isampTimeZero = (int)round( -myTimeSample1_s / mySampleInt_sec );
  for( int isamp = isampTimeZero; isamp < myNumSamples; isamp++ ) {
    double timeOut = (double)isamp*mySampleInt_sec + myTimeSample1_s;
    double timeOut_sq = timeOut * timeOut;
    int isampVel = isamp - isampTimeZero;
    if( isampVel < 0 ) isampVel = 0;
    double vel = myVelocityTrace[isampVel];
    if( vel <= 0.0 ) {
      throw( csException("csNMOCorrection::perform_nmo_internal(): Velocity <= 0 encountered. Wrong input velocity function?") );
    }
    switch( myNMOMethod ) {
    case csNMOCorrection::PP_NMO:
      myTimeTrace[isamp] = sqrt( timeOut_sq + offset_sq / (vel*vel) );
      break;
    case csNMOCorrection::PS_NMO:
      myTimeTrace[isamp] = 0.5 * sqrt( timeOut_sq ) + 0.5 * sqrt( timeOut_sq + 2.0 * offset_sq / (vel*vel) );
      break;
    case csNMOCorrection::PP_NMO_VTI:
      tmp = vel*vel * (timeOut_sq*vel*vel + (1+2*myETATrace[isampVel])*offset_sq);
      if( tmp != 0 ) myTimeTrace[isamp] = sqrt( timeOut_sq + offset_sq/(vel*vel) - (2*myETATrace[isampVel]*offset_sq*offset_sq) / tmp );
      else myTimeTrace[isamp] = 0;
      break;
    case csNMOCorrection::EMPIRICAL_NMO:
      myTimeTrace[isamp] = timeOut + ( ( vel * ( pow( (offset-myOffsetApex)/1000,2) - pow(myOffsetApex/1000,2) ) -  (vel/(vel+0.1)) * myZeroOffsetDamping * exp(-0.5*pow(offset/1000,2)) ) / 1000);
      break;
    case csNMOCorrection::OUTPUT_VEL:
      samplesOut[isamp] = vel;
      break;
    }
  }
  for( int isamp = 0; isamp < isampTimeZero; isamp++ ) {
    myTimeTrace[isamp] = myTimeTrace[isampTimeZero] + myTimeTrace[isampTimeZero] - myTimeTrace[2*isampTimeZero-isamp+1];
  }

  if( myNMOMethod == csNMOCorrection::OUTPUT_VEL ) {
    // Do nothing
  }
  else if( myModeOfApplication == csNMOCorrection::NMO_APPLY ) {
    myInterpol->process( mySampleInt_sec, myTimeSample1_s, myBuffer, myTimeTrace, samplesOut );
  }
  else { // if( myModeOfApplication == csNMOCorrection::NMO_REMOVE ) {
    if( myTimeTraceInverse == NULL ) {
      myTimeTraceInverse  = new float[myNumSamples];
    }
    if( myTimeSample1_s != 0 ) {
      for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
        myTimeTrace[isamp] -= myTimeSample1_s;
      }
    }
    csInterpolation::xy2yxInterpolation( myTimeTrace, myTimeTraceInverse, myNumSamples, mySampleInt_sec );
    myInterpol->process( mySampleInt_sec, 0, myBuffer, myTimeTraceInverse, samplesOut );
  }
}

//--------------------------------------------------------------------------------
//
void csNMOCorrection::perform_nmo_horizonBased_prepare( int numVelocities_in, float const* t0_in, float const* vel_rms_in, double offset, float* samplesOut ) {
  // Set fields. Add t0 & vel_rms at top and bottom if not there yet
  float* t0      = new float[numVelocities_in+2];
  float* vel_rms = new float[numVelocities_in+2];
  int numVelocities = numVelocities_in;
  if( t0_in[0] > 0 ) {
    t0[0] = 0;
    vel_rms[0] = vel_rms_in[0];
    memcpy( &(t0[1]), t0_in, numVelocities_in*sizeof(float) );
    memcpy( &(vel_rms[1]), vel_rms_in, numVelocities_in*sizeof(float) );
    numVelocities += 1;
  }
  else {
    memcpy( t0, t0_in, numVelocities_in*sizeof(float) );
    memcpy( vel_rms, vel_rms_in, numVelocities_in*sizeof(float) );
  }
  if( t0_in[numVelocities_in-1] < myTimeOfLastSample ) {
    numVelocities += 1;
    t0[numVelocities-1]      = myTimeOfLastSample;
    vel_rms[numVelocities-1] = vel_rms_in[numVelocities_in-1];
  }

  //  fprintf(stderr,"Velocities: %d  time: %f %f, vel: %f %f %f   num: %d\n", numVelocities, t0[0], t0[1], vel_rms[0], vel_rms[1], offset, numVelocities_in);

  perform_nmo_horizonBased( numVelocities, t0, vel_rms, offset, samplesOut );
}
//--------------------------------------------------------------------------------
//
void csNMOCorrection::perform_nmo_horizonBased( int numVelocities, float const* times, float const* vel_rms, double offset, float* samplesInOut ) {
  double offset_sq = offset*offset;

  int currentVelIndex = 0;
  double t0Top    = 0;
  double t0Bot    = times[currentVelIndex];
  double velTop_sq = vel_rms[currentVelIndex]*vel_rms[currentVelIndex];
  double velBot_sq = velTop_sq;
  double time_sqDifference = t0Bot*t0Bot - t0Top*t0Top;
  //  double time_difference   = t0Bot - t0Top;

  double nmo_sign = ( myModeOfApplication == csNMOCorrection::NMO_APPLY ) ? 1.0 : -1.0;

  // Go through all output samples, retrieve sample value at NMO time corresponding to zero-offset t0
  for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    double timeOut = (double)isamp*mySampleInt_sec;
    while( currentVelIndex < numVelocities && t0Bot <= timeOut ) {
      t0Top      = t0Bot;
      velTop_sq  = velBot_sq;

      currentVelIndex += 1;
      t0Bot      = times[currentVelIndex];
      velBot_sq  = vel_rms[currentVelIndex] * vel_rms[currentVelIndex];
      //      time_difference   = t0Bot - t0Top;
      time_sqDifference = t0Bot*t0Bot - t0Top*t0Top;
    }

    double timeOut_sq = timeOut * timeOut;
    double weight;
    double v_sq;
    if( myHorInterpolationMethod == csNMOCorrection::HORIZON_METHOD_LINEAR ) {
      weight = (t0Bot - timeOut) / (t0Bot - t0Top);
      v_sq   = weight * (velTop_sq - velBot_sq) + velBot_sq;
    }
    else {
      weight = (t0Bot*t0Bot - timeOut_sq)/time_sqDifference;
      v_sq   = weight * (velTop_sq - velBot_sq) + velBot_sq;
    }
    double vel = sqrt( v_sq );

    double timeIn = 0;
    switch( myNMOMethod ) {
    case csNMOCorrection::PP_NMO:
      timeIn = sqrt( timeOut_sq + nmo_sign * offset_sq / v_sq );
      break;
    case csNMOCorrection::PS_NMO:
      timeIn = 0.5 * sqrt( timeOut_sq ) + 0.5 * sqrt( timeOut_sq + nmo_sign *  2.0 * offset_sq / v_sq );
      break;
    case csNMOCorrection::EMPIRICAL_NMO:
      timeIn = timeOut + ( vel * ( pow( (offset-myOffsetApex)/1000,2) - pow(myOffsetApex/1000,2) ) -  (vel/(vel+0.1)) * myZeroOffsetDamping * exp(-0.5*pow(offset/1000,2)) ) / 1000;
      break;
    case csNMOCorrection::OUTPUT_VEL:
      myBuffer[isamp] = vel;
      break;
    }

    if( myNMOMethod != csNMOCorrection::OUTPUT_VEL ) {
      double sampleIndex = timeIn / mySampleInt_sec;
      myBuffer[isamp] = getQuadAmplitudeAtSample( samplesInOut, sampleIndex, myNumSamples );
    }
  }

  memcpy( samplesInOut, myBuffer, myNumSamples*sizeof(float) );
}

//--------------------------------------------------------------------------------------------------
// Differential NMO:
//
void csNMOCorrection::perform_differential_nmo( float const* samplesIn, int numVelocities_in, float const* time_in, float const* vel_rms_in, double offsetIn, double offsetOut, float* samplesOut ) {
  memcpy(samplesOut, samplesIn, myNumSamples*sizeof(float) );
  perform_differential_nmo( numVelocities_in, time_in, vel_rms_in, offsetIn, offsetOut, samplesOut );
}

void csNMOCorrection::perform_differential_nmo( int numVelocities_in, float const* time_in, float const* vel_rms_in, double offsetIn, double offsetOut, float* samplesOut ) {
  if( myIsHorizonBasedNMO ) {
    throw( csException("Horizon based NMO is not supported for differential NMO") );
  }
  else {
    return( perform_differential_nmo_internal( numVelocities_in, time_in, vel_rms_in, offsetIn, offsetOut, samplesOut, NULL ) );
  }
}

void csNMOCorrection::perform_differential_nmo( csTimeFunction<double> const* velTimeFunc, double offsetIn, double offsetOut, float* samplesOut ) {
  if( myTimeSample1_s != 0 ) throw( csException("csNMOCorrection::perform_differential_nmo(): Time of first sample != 0 not supported yet for differential NMO") );
  int numVelocities_in = velTimeFunc->numValues();
  // Horizon based: Add one vel location at top and bottom of trace (temp fix to make sure that interpolation works)
  int numVelocities    = myIsHorizonBasedNMO ? numVelocities_in+2 : numVelocities_in;

  float* t0      = new float[numVelocities];
  float* vel_rms = new float[numVelocities];

  if( myIsHorizonBasedNMO ) {
    throw( csException("Horizon based NMO is not supported for differential NMO") );
  }
  else {
    for( int i = 0; i < numVelocities_in; i++ ) {
      t0[i]      = velTimeFunc->timeAtIndex(i)/1000.0;  // Convert to seconds
      vel_rms[i] = velTimeFunc->valueAtIndex(i);
    }
    if( myNMOMethod != csNMOCorrection::PP_NMO_VTI ) {
      perform_differential_nmo_internal( numVelocities_in, t0, vel_rms, offsetIn, offsetOut, samplesOut, NULL );
      //      perform_nmo_internal( numVelocities_in, t0, vel_rms, offset, samplesOut, NULL );
    }
    else {
      if( velTimeFunc->numSpatialValues() < 2 ) throw( csException("csNMOCorrection::perform_differential_nmo: VTI velocity function is lacking an eta function") );
      float* eta = new float[numVelocities];
      for( int i = 0; i < numVelocities_in; i++ ) {
        eta[i] = velTimeFunc->valueAtIndex(i,1);
      }
      perform_differential_nmo_internal( numVelocities_in, t0, vel_rms, offsetIn, offsetOut, samplesOut, eta );
      delete [] eta;
    }


  }
  delete [] t0;
  delete [] vel_rms;
}

void csNMOCorrection::perform_differential_nmo_ORIG( csTimeFunction<double> const* velTimeFunc, double offsetIn, double offsetOut, float* samplesOut ) {
  int numVelocities_in = velTimeFunc->numValues();
  // Horizon based: Add one vel location at top and bottom of trace (temp fix to make sure that interpolation works)
  int numVelocities    = myIsHorizonBasedNMO ? numVelocities_in+2 : numVelocities_in;

  float* t0      = new float[numVelocities];
  float* vel_rms = new float[numVelocities];

  if( myIsHorizonBasedNMO ) {
    throw( csException("Horizon based NMO is not supported for differential NMO") );
  }
  else {
    for( int i = 0; i < numVelocities_in; i++ ) {
      t0[i]      = velTimeFunc->timeAtIndex(i)/1000.0;  // Convert to seconds
      vel_rms[i] = velTimeFunc->valueAtIndex(i);
    }

    perform_differential_nmo_internal_ORIG( numVelocities_in, t0, vel_rms, offsetIn, offsetOut, samplesOut );    
  }
  delete [] t0;
  delete [] vel_rms;
}


void csNMOCorrection::perform_differential_nmo_internal( int numVels, float const* times, float const* vel_rms, double offsetIn, double offsetOut, float* samplesOut, float const* eta ) {

  // Linearly interpolate input velocities
  csInterpolation::linearInterpolation( numVels, times, vel_rms, myNumSamples, mySampleInt_sec, myVelocityTrace );
  if( myNMOMethod == csNMOCorrection::PP_NMO_VTI ) {
    csInterpolation::linearInterpolation( numVels, times, eta, myNumSamples, mySampleInt_sec, myETATrace );
  }

  double offsetIn_sq = offsetIn*offsetIn;
  double offsetOut_sq = offsetOut*offsetOut;

  if( myTimeTraceInverse == NULL ) {
    myTimeTraceInverse  = new float[myNumSamples];
  }
  if( myTimeTraceDiff == NULL ) {
    myTimeTraceDiff  = new float[myNumSamples];
  }         

  int isampTimeZero = (int)round( -myTimeSample1_s / mySampleInt_sec );
  for( int isamp = isampTimeZero; isamp < myNumSamples; isamp++ ) {
    double timeOut    = (double)isamp*mySampleInt_sec + myTimeSample1_s;
    double timeOut_sq = timeOut * timeOut;
    int isampVel = isamp - isampTimeZero;
    if( isampVel < 0 ) isampVel = 0;
    double vel = (isampVel >= 0) ? myVelocityTrace[isampVel] : myVelocityTrace[0];
    if( vel <= 0.0 ) {
      throw( csException("csNMOCorrection::perform_nmo_internal(): Velocity <= 0 encountered. Wrong input velocity function?") );
    }
    switch( myNMOMethod ) {
    case csNMOCorrection::PP_NMO:
      myTimeTrace[isamp]     = sqrt( timeOut_sq + offsetIn_sq / (vel*vel) );
      myTimeTraceDiff[isamp] = sqrt( timeOut_sq + offsetOut_sq / (vel*vel) );
      break;
    case csNMOCorrection::PP_NMO_VTI:
      myTimeTrace[isamp]     = sqrt( timeOut_sq + offsetIn_sq/(vel*vel) - (2*myETATrace[isampVel]*offsetIn_sq*offsetIn_sq) / (vel*vel * (timeOut_sq*vel*vel + (1+2*myETATrace[isampVel])*offsetIn_sq) ) );
      myTimeTraceDiff[isamp] = sqrt( timeOut_sq + offsetOut_sq/(vel*vel) - (2*myETATrace[isampVel]*offsetOut_sq*offsetOut_sq) / (vel*vel * (timeOut_sq*vel*vel + (1+2*myETATrace[isampVel])*offsetOut_sq) ) );
      break;
    case csNMOCorrection::PS_NMO:
      myTimeTrace[isamp] = 0.5 * sqrt( timeOut_sq ) + 0.5 * sqrt( timeOut_sq + 2.0 * offsetIn_sq / (vel*vel) );
      myTimeTraceDiff[isamp] = 0.5 * sqrt( timeOut_sq ) + 0.5 * sqrt( timeOut_sq + 2.0 * offsetOut_sq / (vel*vel) );
      break;
    case csNMOCorrection::EMPIRICAL_NMO:
      myTimeTrace[isamp] = timeOut + ( ( vel * ( pow( (offsetIn-myOffsetApex)/1000,2) - pow(myOffsetApex/1000,2) ) -  (vel/(vel+0.1)) * myZeroOffsetDamping * exp(-0.5*pow(offsetIn/1000,2)) ) / 1000);
      myTimeTraceDiff[isamp] = timeOut + ( ( vel * ( pow( (offsetOut-myOffsetApex)/1000,2) - pow(myOffsetApex/1000,2) ) -  (vel/(vel+0.1)) * myZeroOffsetDamping * exp(-0.5*pow(offsetOut/1000,2)) ) / 1000);
      break;
    }
  }
  /*
...review when implementing timesample1_s != 0
  for( int isamp = 0; isamp < isampTimeZero; isamp++ ) {
    myTimeTrace[isamp] = myTimeTrace[isampTimeZero] + myTimeTrace[isampTimeZero] - myTimeTrace[2*isampTimeZero-isamp+1];
    myTimeTraceDiff[isamp] = myTimeTraceDiff[isampTimeZero] + myTimeTraceDiff[isampTimeZero] - myTimeTraceDiff[2*isampTimeZero-isamp+1];
  }
  */

  // Ver2: Brute force, apply NMO, then inverse NMO
  memcpy( myBuffer, samplesOut, myNumSamples*sizeof(float) );
  if( myModeOfApplication == csNMOCorrection::NMO_APPLY ) {
    csInterpolation::xy2yxInterpolation( myTimeTraceDiff, myTimeTraceInverse, myNumSamples, mySampleInt_sec );
    myInterpol->process( mySampleInt_sec, 0.0, myBuffer, myTimeTrace, samplesOut );
    memcpy( myBuffer, samplesOut, myNumSamples*sizeof(float) );
    myInterpol->process( mySampleInt_sec, 0.0, myBuffer, myTimeTraceInverse, samplesOut );
  }
  else {
    csInterpolation::xy2yxInterpolation( myTimeTrace, myTimeTraceInverse, myNumSamples, mySampleInt_sec );
    myInterpol->process( mySampleInt_sec, 0.0, myBuffer, myTimeTraceDiff, samplesOut );
    memcpy( myBuffer, samplesOut, myNumSamples*sizeof(float) );
    myInterpol->process( mySampleInt_sec, 0.0, myBuffer, myTimeTraceInverse, samplesOut );
  }
  
  //  fprintf(stderr,"Diff NMO %.2f --> %.2f    %f\n", offsetIn, offsetOut, myTimeTrace[1000]);

  // Ver1: Compute time squeeze function, then apply NMO once. Only coded for APPLY yet. Does not keep data above direct arrival
  /*
  for( int isamp = 0; isamp < isampTimeZero; isamp++ ) {
    myTimeTrace[isamp] = myTimeTrace[isampTimeZero] + myTimeTrace[isampTimeZero] - myTimeTrace[2*isampTimeZero-isamp+1];
    myTimeTraceDiff[isamp] = myTimeTraceDiff[isampTimeZero] + myTimeTraceDiff[isampTimeZero] - myTimeTraceDiff[2*isampTimeZero-isamp+1];
  }
  for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    myTimeTraceInverse[isamp] = 0;
  }
  int index1 = (int)myTimeTraceDiff[0];
  int index2 = (int)myTimeTraceDiff[myNumSamples-1];

  for( int isamp = 0; isamp < index1; isamp++ ) {
    myTimeTraceInverse[isamp] = 0;
  }
  for( int isamp = index2; isamp < myNumSamples; isamp++ ) {
    myTimeTraceInverse[isamp] = 0;
  }
  for( int isamp = 1; isamp < myNumSamples; isamp++ ) {
    float timeIn1 = myTimeTrace[isamp-1];
    float timeIn2 = myTimeTrace[isamp];
    float timeOut1 = myTimeTraceDiff[isamp-1];
    float timeOut2 = myTimeTraceDiff[isamp];
    float indexOut1 = timeOut1/mySampleInt_sec;
    float indexOut2 = timeOut2/mySampleInt_sec;
    index1 = (int)indexOut1+1;
    index2 = (int)indexOut2;
    if( (float)index2 < indexOut1 ) continue; // Skip too finely sampled time steps
    if( index1 >= 0 && index2 <= myNumSamples-1 ) {
      for( int index = index1; index <= index2; index++ ) {
        float ratio = ( (float)index - indexOut1 ) / ( indexOut2 - indexOut1 );
        myTimeTraceInverse[index] = timeIn1 + ratio * ( timeIn2 - timeIn1 );
        //        fprintf(stdout,"%d %f\n", index, index*mySampleInt_sec );
      }
    }
  }
  index1 = (int)myTimeTraceDiff[0];
  for( int isamp = 0; isamp < index1; isamp++ ) {
    int indexTmp = 2*index1-isamp+1;
    if( indexTmp < myNumSamples ) {
      myTimeTraceInverse[isamp] = 2*myTimeTraceInverse[index1] - myTimeTraceInverse[indexTmp];
    }
    else {
      myTimeTraceInverse[isamp] = 0;
    }
  }

  memcpy( myBuffer, samplesOut, myNumSamples*sizeof(float) );
  myInterpol->process( mySampleInt_sec, 0.0, myBuffer, myTimeTraceInverse, samplesOut );
  */
  

  // Assume APPLY: ...another option
  /*
  csInterpolation::xy2yxInterpolation( myTimeTraceDiff, myTimeTraceInverse, myNumSamples, mySampleInt_sec );
  for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    fprintf(stdout,"%d  %.4f %.4f %.4f %.4f\n", isamp, isamp*mySampleInt_sec, myTimeTrace[isamp], myTimeTraceDiff[isamp], myTimeTraceInverse[isamp]);
  }
  for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    float time = myTimeTraceInverse[isamp];
    int index1 = (int)(time/mySampleInt_sec);
    float timeIndex1 = index1*mySampleInt_sec;
    if( index1 < 0 || index1 >= myNumSamples-2 ) {
      myTimeTraceDiff[isamp] = 0;
    }
    else {
      float time1 = myTimeTrace[index1];
      float time2 = myTimeTrace[index1+1];
      myTimeTraceDiff[isamp] = time1 + ( (time-timeIndex1)/mySampleInt_sec ) * ( time2 - time1 );
    }
  }
  myInterpol->process( mySampleInt_sec, 0.0, myBuffer, myTimeTraceDiff, samplesOut );
  */


  /*
  if( myModeOfApplication == csNMOCorrection::NMO_APPLY ) {
    myInterpol->process( mySampleInt_sec, myTimeSample1_s, myBuffer, myTimeTrace, samplesOut );
  }
  else { // if( myModeOfApplication == csNMOCorrection::NMO_REMOVE ) {
    if( myTimeTraceInverse == NULL ) {
      myTimeTraceInverse  = new float[myNumSamples];
    }
    if( myTimeSample1_s != 0 ) {
      for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
        myTimeTrace[isamp] -= myTimeSample1_s;
      }
    }
    csInterpolation::xy2yxInterpolation( myTimeTrace, myTimeTraceInverse, myNumSamples, mySampleInt_sec );
    myInterpol->process( mySampleInt_sec, 0, myBuffer, myTimeTraceInverse, samplesOut );
  }
*/
}

void csNMOCorrection::perform_differential_nmo_internal_ORIG( int numVels, float const* times, float const* vel_rms, double offsetIn, double offsetOut, float* samplesOut ) {

  // Linearly interpolate input velocities
  csInterpolation::linearInterpolation( numVels, times, vel_rms, myNumSamples, mySampleInt_sec, myVelocityTrace );

  memcpy( myBuffer, samplesOut, myNumSamples*sizeof(float) );
  double offsetIn_sq = offsetIn*offsetIn;
  double offsetOut_sq = offsetOut*offsetOut;

  for( int isamp = 0; isamp < myNumSamples; isamp++ ) {
    double timeOut = (double)isamp*mySampleInt_sec;
    double timeOut_sq = timeOut * timeOut;
    double vel        = myVelocityTrace[isamp];
    if( vel <= 0.0 ) {
      throw( csException("csNMOCorrection::perform_nmo_internal(): Velocity <= 0 encountered. Wrong input velocity function?") );
    }
    switch( myNMOMethod ) {
    case csNMOCorrection::PP_NMO:
      myTimeTrace[isamp] = timeOut + sqrt( timeOut_sq + offsetIn_sq / (vel*vel) ) - sqrt( timeOut_sq + offsetOut_sq / (vel*vel) );
      break;
    case csNMOCorrection::PS_NMO:
      myTimeTrace[isamp] = timeOut + 0.5 * sqrt( timeOut_sq ) + 0.5 * sqrt( timeOut_sq + 2.0 * offsetIn_sq / (vel*vel) ) - ( 0.5 * sqrt( timeOut_sq ) + 0.5 * sqrt( timeOut_sq + 2.0 * offsetOut_sq / (vel*vel) ) );
      break;
    case csNMOCorrection::EMPIRICAL_NMO:
      myTimeTrace[isamp] = timeOut + ( ( vel * ( pow( (offsetIn-myOffsetApex)/1000,2) - pow(myOffsetApex/1000,2) ) -  (vel/(vel+0.1)) * myZeroOffsetDamping * exp(-0.5*pow(offsetIn/1000,2)) ) / 1000) -
        (( ( vel * ( pow( (offsetOut-myOffsetApex)/1000,2) - pow(myOffsetApex/1000,2) ) -  (vel/(vel+0.1)) * myZeroOffsetDamping * exp(-0.5*pow(offsetOut/1000,2)) ) / 1000));
      break;
    }
    //    fprintf(stdout,"%f %f %f\n", timeOut, myTimeTrace[isamp], myVelocityTrace[isamp]);
  }

  if( myModeOfApplication == csNMOCorrection::NMO_APPLY ) {
    myInterpol->process( mySampleInt_sec, 0.0, myBuffer, myTimeTrace, samplesOut );
  }
  else { // if( myModeOfApplication == csNMOCorrection::NMO_REMOVE ) {
    if( myTimeTraceInverse == NULL ) {
      myTimeTraceInverse  = new float[myNumSamples];
    }
    csInterpolation::xy2yxInterpolation( myTimeTrace, myTimeTraceInverse, myNumSamples, mySampleInt_sec );
    myInterpol->process( mySampleInt_sec, 0.0, myBuffer, myTimeTraceInverse, samplesOut );
  }
}



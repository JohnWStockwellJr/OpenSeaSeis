/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_DESPIKE_H
#define CS_DESPIKE_H

#include <cstdio>
#include <cstdlib>
#include <cmath>

namespace cseis_geolib {

template <typename T> class csVector;
class csDespike;
struct DespikeConfig;


/**
 * Despike class.
 *
 * Removes spikes and/or short noise bursts from input data.
 *
 * @author Bjorn Olofsson
 * @date 2008
 */
class csDespike {
private:
  struct DespikePoint {
    DespikePoint( int theIndex, float theValue ) {
      index = theIndex;
      value = theValue;
    }
    int   index;
    float value;
  };

public:
  /// COSINE_TAPER: Apply cosine taper around each identified spike event
  static int const COSINE_TAPER = 212;  /// SET_TO_ZERO: Set al samples around each identified spike event to zero
  static int const SET_TO_ZERO  = 213;
  /// LINEAR_INTERPOLATION: Interpolate samples around each identified spike event, between first and last sample
  static int const LINEAR_INTERPOLATION = 214;

  static int const METHOD_WIN_MEAN   = 301;
  static int const METHOD_WIN_MEDIAN = 302;

public:
  csDespike( int numSamples, double sampleInt, DespikeConfig& config );
  ~csDespike();
  static void getDefaultFrequencySpikeConfig( DespikeConfig& config );
  static void getDefaultTimeNoiseBurstConfig( DespikeConfig& config );
  void set( DespikeConfig const& config );
  /**
  * Apply despiking operation to data
  */
  void apply( float* samples, int numSamples, int& numSpikesFound, int& numSamplesAffected );
  int method() const { return myMethod; }

private:
  int    myNumSamples;
  double mySampleInt;
  bool   myPerformDebias;
  int    myMethod;
  int    myMethodRefWin;
  int    myMethodSpikeWin;

  // Widths, increments etc are all in number of samples
  int myIncWin;
  int myWidthRefWin;   // Must be odd number of samples
  int myWidthMeanWin;
  int myWidthSpikeWin;
  int myStartSample;
  int myStopSample;
  int myNumWindows;

  float myMaxRatio;
  float* myRatios;
  float* myValueRefWin;
  
  // Special advanced settings
  int myRatioAmplifier;
  int myWindowAmplifier;

  void init( int numSamples, double sampleInt );
  bool computeRatios( float* samples, int numSamples );
  void insertNewValue( DespikePoint* p, cseis_geolib::csVector<DespikePoint*>* sortedList );
};

/**
 * Configuration values for Despike class
 * All values are in units of the current trace data,
 * for example [ms] or [Hz]
 */
struct DespikeConfig {
  DespikeConfig() {
    widthRefWin = 0;
    incWin      = 0;
    start       = 0;
    stop        = 0;
    minWidthSpikeWin = 0;
    advanced.widthMeanWin_inSamples = 11;
    advanced.ratioAmplifier  = 2;
    advanced.windowAmplifier = 4;
    performDebias = false;
    maxRatio   = 0.0;
    method     = csDespike::COSINE_TAPER;
    methodRefWin = csDespike::METHOD_WIN_MEDIAN;
    methodSpikeWin = csDespike::METHOD_WIN_MEAN;
  }
  /// Width of reference window (for example in [ms] or [Hz], depending on data domain)
  /// The reference window is used to determine the energy level of the background
  float widthRefWin;
  /// Increment of window. Set to zero to use one sample step = sliding window
  float incWin;
  /// Start time/frequency for despiking. All samples before 'start' will not be despiked
  double start;
  /// Stop time/frequency for despiking. All samples after 'stop' will not be despiked
  double stop;
  /// Minimum width of spike removal window.
  /// All samples in the spike window will be despiked if a spike has been found.
  /// If anomalous amplitudes are found in a wider window, a wider window will be used
  float minWidthSpikeWin;
  /// Remove DC bias before determining spikes
  bool performDebias;
  float maxRatio;
  int method;
  int methodRefWin;
  int methodSpikeWin;

  struct Advanced {
    /// Width of window over which mean (spike) amplitude is determined
    /// Keeping it too short and each single sample value will be considered.
    /// Making it too long and short spikes will be overlooked
    int widthMeanWin_inSamples;
    int ratioAmplifier;
    int windowAmplifier;
  } advanced;
};

}  // end namespace

#endif


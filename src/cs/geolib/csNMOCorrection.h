/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_NMO_CORRECTION_H
#define CS_NMO_CORRECTION_H

namespace cseis_geolib {

 template<typename T> class csTimeFunction;
 class csInterpolation;

/**
  * Performs NMO correction for a single trace.
  *
  * Uses one of three NMO interpolation methods.
  * NOTE: This is an adaptation for model based NMO correction, not entirely suitable for NMO & stacking
  *
  * @author Bjorn Olofsson
  * @date 2005
  */
class csNMOCorrection
{
public:
  static const int NMO_APPLY  = 11;
  static const int NMO_REMOVE = 12;
  static const int PP_NMO = 1;
  static const int PS_NMO = 2;
  static const int EMPIRICAL_NMO = 3;
  static const int OUTPUT_VEL = 4;
  static const int PP_NMO_VTI = 5;

  static const int HORIZON_METHOD_LINEAR = 21;
  static const int HORIZON_METHOD_QUAD   = 22;
  //  static const int METHOD_DIFF  = 23;

public:
  /**
  * sampleInt [ms]
  * numSamples Number of samples
  * method_interpolation NMO interpolation method
  */
  csNMOCorrection( double sampleInt_ms, int numSamples, int method_nmo = csNMOCorrection::PP_NMO );
  ~csNMOCorrection();

  void setModeOfApplication( int mode );
  void setEmpiricalNMO( double offsetApex_m, double zeroOffsetDamping );
  void setHorizonBasedNMO( bool setHorizonBased, int method_interpolation = csNMOCorrection::HORIZON_METHOD_LINEAR );
  /**
   * Set time of first sample in input data
   */
  void setTimeSample1( float timeSample1_ms );

  /**
  * samples    Trace samples
  * vel_rms    RMS velocities, one function valid for this trace [m/s]
  * t0         Zero-offset times at which RMS velocities are defined [s]
  * numVelocities  Number of velocity/t0 values
  * offset     Trace offset [m]
  */
  //void perform_nmo( float* samplesOut, float const* vel_rms, float const* t0_s, int numVelocities, double offset );
  //void perform_nmo( float* samplesOut, float const* vel_rms, float const* t0_s, int numVelocities, double offset, float const* samplesIn );
  //void perform_nmo( float* samplesOut, csTimeFunction<double> const*, double offset );

  void perform_nmo( int numVelocities_in, float const* time_in, float const* vel_rms_in, double offset, float* samplesOut );
  void perform_nmo( float const* samplesIn, int numVelocities_in, float const* time_in, float const* vel_rms_in, double offset, float* samplesOut );
  void perform_nmo( csTimeFunction<double> const* velTimeFunc, double offset, float* samplesOut );

  void perform_differential_nmo( int numVelocities_in, float const* time_in, float const* vel_rms_in, double offsetIn, double offsetOut, float* samplesOut );
  void perform_differential_nmo( float const* samplesIn, int numVelocities_in, float const* time_in, float const* vel_rms_in, double offsetIn, double offsetOut, float* samplesOut );
  void perform_differential_nmo( csTimeFunction<double> const* velTimeFunc, double offsetIn, double offsetOut, float* samplesOut );
  void perform_differential_nmo_ORIG( csTimeFunction<double> const* velTimeFunc, double offsetIn, double offsetOut, float* samplesOut );

private:
  void perform_nmo_internal( int numVels, float const* times, float const* vel_rms, double offset, float* samplesOut, float const* eta );
  void perform_nmo_horizonBased_prepare( int numVelocities_in, float const* times, float const* vel_rms_in, double offset, float* samplesOut );
  void perform_nmo_horizonBased( int numVelocities, float const* times, float const* vel_rms, double offset, float* samplesInOut );

  void perform_differential_nmo_internal( int numVels, float const* times, float const* vel_rms, double offsetIn, double offsetOut, float* samplesOut );
  void perform_differential_nmo_internal( int numVels, float const* times, float const* vel_rms, double offsetIn, double offsetOut, float* samplesOut, float const* eta );
  void perform_differential_nmo_internal_ORIG( int numVels, float const* times, float const* vel_rms, double offsetIn, double offsetOut, float* samplesOut );


  //  void perform_nmo_horizonBased( float* samplesOut, float const* vel_rms, float const* t0_s, int numVelocities, double offset, float const* samplesIn );
  // void perform_nmo_internal( double offset, float* samplesOut );
  // void perform_nmo_internal( int numVels, float const* times, float* vels, double offset, float* samplesOut );
  void allocateVelocity();

  /// NMO 'method': PP, PS ...
  int myNMOMethod;
  /// NMO interpolation method for horizon based interpolation
  int myHorInterpolationMethod;
  /// Time of first sample, typically specified in seconds [s] (otehrwise: unit of input data and velocity field).
  float myTimeSample1_s;

  float* myBuffer;
  /// Number of samples expected in the traces to be processed
  int    myNumSamples;
  /// Sample interval in seconds
  double mySampleInt_sec;    // [s]
  /// Time [s] of last sample
  double myTimeOfLastSample;  // [s]
  /// 'Mode' of application: Apply or remove NMO
  int myModeOfApplication;

  float* myVelocityTrace;
  float* myETATrace;
  float* myTimeTrace;
  float* myTimeTraceInverse;
  float* myTimeTraceDiff;
  csInterpolation* myInterpol;

  /// 'Horizon-based' NMO (stretch happens mostly between horizons, which can lead to sharp 'jumps' at the horizons)
  bool myIsHorizonBasedNMO;

  /// Parameters for 'empirical' NMO function:
  /// Offset at which residual NMO apex appears
  double myOffsetApex;
  /// Amount of 'damping' to apply at zero offset. Set to high value to reduce residual correction at zero offset
  double myZeroOffsetDamping;

};

} // namespace
#endif


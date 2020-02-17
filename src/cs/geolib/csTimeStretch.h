/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TIME_STRETCH
#define CS_TIME_STRETCH

namespace cseis_geolib {

template<typename T> class csTimeFunction;
class csInterpolation;

/**
  * Performs NMO correction for a single trace.
  * Uses one of three NMO interpolation methods.
  *
  * @author Bjorn Olofsson
  * @date 2005
  */
class csTimeStretch
{
public:
  static const int LAYER_INTERPOL_LIN   = 21;
  static const int LAYER_INTERPOL_QUAD  = 22;
  static const int LAYER_INTERPOL_SLOPE = 23;

  static const int SAMPLE_INTERPOL_LIN  = 30;
  static const int SAMPLE_INTERPOL_QUAD = 31;
  static const int SAMPLE_INTERPOL_SINC = 32;

public:
  csTimeStretch( double sampleInt_ms, int numSamples );
  csTimeStretch( double sampleInt_ms, int numSamples, int methodSampleInterpolation );
  ~csTimeStretch();

  void setLayerInterpolationMethod( int method );

  /**
   * Set stretch behaviour for bottom layer (below last specified layer)
   * @param stretchBottom  true: Apply same stretch as in last layer, false: Only apply static shift to bottom of trace
   */
  void setBottomStretch( bool stretchBottom ) { myIsBottomStretch = stretchBottom; }

  /**
   * Apply stretch function to input data
   * @param samplesIn   Input data samples
   * @param timesIn_ms  List of times defining each layer. First time is top of layer 1, second time is bottom of layer 1 (=top of layer 2) etc 
   * @param stretch_ms  List of stretch factors to apply to each layer. This list contains one value less than the time list, one for each layer.
   * @param numTimes    Number of time values.
   * @param samplesOut  Output data samples
   */
  void applyStretchFunction( float const* samplesIn, float const* timesIn_ms, float const* stretch_ms, int numTimes, float* samplesOut );
  /**
   * Apply stretch to input data.
   * Same as applyStretchFunction() but providing pairs of input/output horizon times instead of layer stretch function
   */
  void applyTimeInterval( float const* samplesIn, float const* timesIn_ms, float const* timesOut_ms, int numTimes, float* samplesOut );
  /**
   */
  void initialize_timeDepthConversion( int numSamplesOut, float sampleIntOut );
  void apply_timeDepthConversion( float const* samplesIn, float const* velFunc, float sampleIntVel_m, int numSamplesVel, float* samplesOut, bool isTime2Depth );
  //  void convert_depth2time( float const* samplesIn, float const* velFunc, float sampleIntVel_m, int numSamplesVel, float* samplesOut );
  //  void convert_time2depth( float const* samplesIn, float const* velFunc, float sampleIntVel_m, int numSamplesVel, float* samplesOut );
private:
  void init( double sampleInt_ms, int numSamples, int methodSampleInterpolation );
  csInterpolation* mySampleInterpolation;
  int    myLayerInterpolationMethod;
  int    mySampleInterpolationMethod;
  int    myNumSamples;
  double mySampleInt;        // [ms]
  double myTimeOfLastSample;  // [s]
  bool   myIsBottomStretch;

  float* myIndexBufferOut;
  float* myTimeDepthFunc;
  int myNumSamplesOut;
  float mySampleIntOut;
  cseis_geolib::csInterpolation* myInterpol_timeDepthConversion;
};

} // namespace
#endif


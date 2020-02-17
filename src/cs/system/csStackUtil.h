/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_STACK_UTIL_H
#define CS_StACK_UTIL_H

#include <map>

namespace cseis_geolib {
  template <typename T> class csVector;
}

namespace cseis_system {
class csTraceGather;
class csTrace;
class csTraceHeader;

/**
 * Utility class for seismic stacking operation
 * Helps with the logistics of stacking seismic traces
 *
 * @author Bjorn Olofsson
 * @date 2008
 */
class csStackUtil {
public:
  static int const OUTPUT_FIRST   = 1;
  static int const OUTPUT_LAST    = 2;
  static int const OUTPUT_AVERAGE = 3;

public:
  /**
  * Constructor
  * @param nSamples    Number of samples in traces to be stacked
  * @param normFactor  Normalisation factor: Stacked traces will be normalised by (1 / ntraces^normFactor)
  * @param outputOption  Option defining what output is required: OUTPUT_FIRST, OUTPUT_LAST, OUTPUT_AVERAGE
  */
  csStackUtil( int nSamples, float normFactor, int outputOption );
  csStackUtil( int nSamples, float normFactor, int outputOption, bool normTimeVariant );
  ~csStackUtil();
  void stackAndNormalizeTraces( csTraceGather* traceGather );

  void stackTraceOld( csTrace* stackedTrace, csTrace const* traceIn );
  void stackTrace( csTrace* stackedTrace );
  void stackTrace( csTrace* stackedTrace, csTrace const* traceIn );
  void normStackedTraceOld( csTrace* trace, int nTraces );
  void normStackedTrace( csTrace* trace, int nTraces );

  void setOutputNormTrace( bool doOutputNormTrace );
  void setTimeVariantNorm( bool doTimeVariantNorm, int hdrId_keyValue, float zeroThreshold );
private:
  void init( int nSamples, float normFactor, int outputOption );
  void stackTrace( csTrace* stackedTrace, csTrace const* traceIn, bool newTrace );
  int* getNormTrace_internal( int keyValue, bool createNew );
  void clearNormTrace_internal();
  int myOutputOption;
  int myNumSamples;
  float myNormFactor;
  bool myNormTimeVariant;
  //  int* myNormTrace;
  bool myOutputNormTrace;
  void stackHeaders( csTraceHeader* trcHdrOut, csTraceHeader const* trcHdrIn );
  void normHeaders( csTraceHeader* trcHdrOut, int nTraces );

  cseis_geolib::csVector<int*>* myNormTraceList;
  std::map<int,int>* myNormTraceIndexMap;
  int myHdrId_keyValue;
  float myZeroThreshold;
};

} // end namespace

#endif


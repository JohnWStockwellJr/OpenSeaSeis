/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEGY_BIN_HEADER_H
#define CS_SEGY_BIN_HEADER_H

#include "geolib_defines.h"
#include <cstdio>

namespace cseis_geolib {

/**
* SEGY binary header
*
* @author Bjorn Olofsson
* @date 2006
*/
class csSegyBinHeader {
public:
  csSegyBinHeader( bool doSwapEndian );
  ~csSegyBinHeader();
  void extractHeaders( byte_t const* buffer );
  void encodeHeaders( byte_t* buffer ) const;
  void dump( std::FILE* fout ) const;
public:
  static const int HDR_SIZE = 400;
  int jobID;
  int lineNum;
  int reelNum;
  int numTraces;
  int numAuxTraces;
  int sampleIntUS;
  int sampleIntOrigUS;
  int numSamples;
  int numSamplesOrig;

  /// SEGY data format
  /// 1 = 32-bit IBM floating point
  /// 2 = 32-bit fixed-point (integer)
  /// 3 = 16-bit fixed-point (integer)
  /// 4 = 32-bit fixed-point with gain code (obsolete)
  /// 5 = 32-bit IEEE floating point
  unsigned short dataSampleFormat;
  unsigned short fold;
  /**
  * Sorting code
  * -1 = Other
  * 0 = Unknown
  * 1 = As recorded (no sorting)
  * 2 = CDP ensemble
  * 3 = Single fold continuous profile
  * 4 = horizontally stacked
  * 5 = Common source point
  * 6 = Common receiver point
  * 7 = Common offset point
  * 8 = Common mid-point
  * 9 = Common conversion point
  */
  unsigned short sortCode;
  /// Vertical sum code. 1 = no sum, 2 = two sum ... N=M-1 (M=2-32767)
  unsigned short vertSumCode;
  unsigned short sweepFreqStart; // Hz
  unsigned short sweepFreqEnd;   // Hz
  unsigned short sweepCode;
  unsigned short taperType;
  unsigned short correlatedTraces;
  unsigned short gainRecovered;
  /// Amplitude recovery method. 1 = none, 2 = spherical divergence, 3 = AGC, 4 = other
  unsigned short ampRecoveryMethod;
  /// Measurement system. 1 = Meters, 2 = Feet
  unsigned short unitSystem;
  /// Impule signal polarity. 1 = Increase in pressure is negative number, 2 = otherwise
  unsigned short polarity;
  unsigned short vibPolarityCode;
  /// SEGY revision number. 0100 means version 1.00.
  unsigned short revisionNum;
  unsigned short fixedTraceLengthFlag;
  unsigned short numExtendedBlocks;

private:
  csSegyBinHeader();
  bool myDoSwapEndian;
};

} // end namespace

#endif


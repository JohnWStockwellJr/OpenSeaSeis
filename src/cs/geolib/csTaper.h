/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TAPER_H
#define CS_TAPER_H

namespace cseis_geolib {

class csTaper {
 public:
  csTaper();
  csTaper( int numSamplesIn, int numSamplesTaper, int taperType, int normFlag );
  ~csTaper();
  void set( int numSamplesIn, int numSamplesTaper, int taperType, int normFlag );
  /**
   * Apply taper to data samples
   *
   * @param dataSamples (i/o) : 
   * @param numSamplesIn (i)  : Number of samples
   * @param numSamplesTaper (i) : Taper length in samples
   * @param taperType  (i)    : Taper type
   * @param normFlag   (i)    : Normalization flag
   * @param normScalar (o)    : Normalization scalar
   */
  void applyTaper( float* dataSamples, int numSamplesIn, int numSamplesTaper, int taperType, int normFlag, float& normScalar );
  void apply( float* dataSamples, float& normScalar );

 public:
  static const int NONE    = -1;
  static const int COSINE  = 1;
  static const int HANNING = 2;
  static const int BLACKMAN = 3;

  static const int NORM_NO   = 1;
  static const int NORM_YES  = 2;
  static const int NORM_ZERO = 3;
  static const int NORM_NSAMP = 10;
  static const int NORM_SQRT  = 11;

 private:
  int myTaperType;
  int myNormFlag;
  int myNumSamplesIn;
  int myNumSamplesTaper;
};

} // END namespace

#endif

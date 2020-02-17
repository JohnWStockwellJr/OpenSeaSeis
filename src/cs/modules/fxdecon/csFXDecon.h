/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_FX_DECON_H
#define CS_FX_DECON_H

#include <complex>

namespace cseis_geolib {
  class csFFT;
}

namespace mod_fxdecon {

struct Attr {
  float fmin;
  float fmax;
  int winLen_samp;
  int taperLen_samp;
  int ntraces_design;
  int ntraces_filter;
  int numWin;
  float taperLen_s;
};

class csFXDecon {
 public:
  csFXDecon();
  ~csFXDecon();
  void initialize( float sampleInt_ms, int numSamplesIn, mod_fxdecon::Attr const& attr );
  void apply( float** samplesIn, float** samplesOut, int numTraces, int numSamples );
  void dump() const;

  static void cxcor( int num1, int index1, std::complex<float> *in1,
                     int num2, int index2, std::complex<float> *in2, 
                     int numCorr, int indexCorr, std::complex<float> *corr );
  static void cconv( int num1, int index1, std::complex<float>* in1,
                     int num2, int index2, std::complex<float>* in2, 
                     int numCorr, int indexCorr, std::complex<float>* corr );
  static void LUDecomposition( float** AA, int nIn, int* indx, float* dd );
  static void LUBackSub( float** AA, int nIn, int* indx, float* bb );

 private:
  void initialize_internal( int numTraces, int numSamples );
  void freeMem( int numTraces );

 private:
  cseis_geolib::csFFT* myFFT;
  int myNumSamplesFFT;
  int myNumFreq;
  float myFreqStep_hz;
  float mySampleInt_s;
  int myNumSamplesWinF;
  int myNumSamplesWinI;
  int myNumTraces;

  float fmin;
  float fmax;
  int winLen_samp;
  int taperLen_samp;
  int ntraces_design;
  int ntraces_filter;
  int numWin;
  float taperLen_s;

  int ntraces_fdataw;

 private:
  std::complex<float>** fdata;
  std::complex<float>** fdataw;
  std::complex<float>** ffreq;
  float** rmatrix;

  std::complex<float>* fvector;
  std::complex<float>*  sfreq;
  std::complex<float>* autocorr;
  float* rautoc;
  float* iautoc;
  float* gvector;

  float* info;
  int* ipvt;
  float* bufferRealImag;

  std::complex<float>* sfreqout;
  float* tidataw;
  float* ttodataw;

};

} // end namespace

#endif

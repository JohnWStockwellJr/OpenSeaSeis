/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csGeolibUtils.h"
#include "csFFT.h"
#include <cmath>
#include <cstring>
#include <limits>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: FFT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_fft {
  struct VariableStruct {
    csFFT* fft;

    int direction;       // Direction of transform: 'forward' or 'reverse'

    int normMethod;      // Normalization method
    bool apply2x_doubleSided; // Apply factor 2x for double-sided amplitude/psd spectra
    float normValue;     // Normalisation value to be applied after forward FFT (or before inverse FFT)
    bool nonZeroSamplesOnly; // Use non-zero samples only for time averaging normalization?
    
    int fftDataType;        // Datatype of the transformed data after forward transform 

    int taperType;            // Type of taper tp apply to input 
    int taperLengthInSamples; // Taper length in number of samples (from 0 to 1)

    int numSamplesIn;    // Number of samples input (4-byte float)
    float sampleIntIn;   // Sample rate of input data (ms or hz).
  };
  static const int TAPER_NONE    = -1;
  static const int TAPER_COSINE  = 1;
  static const int TAPER_HANNING = 2;
  static const int TAPER_BLACKMAN = 3;

  static const int FORWARD = 11;
  static const int INVERSE = 12;

  static const int NORM_NONE         = 1;
  static const int NORM_STANDARD     = 2;
  static const int NORM_TIME_AVERAGE = 3;
  static const int NORM_TIME_AVERAGE_SQRT = 4;
}
using mod_fft::VariableStruct;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_fft_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader*  shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  edef->setMPISupport( true ); // MPI is supported by this module

  vars->direction      = mod_fft::FORWARD;
  vars->normValue      = 1.0f;
  vars->normMethod     = mod_fft::NORM_STANDARD;
  vars->apply2x_doubleSided = true;
  vars->nonZeroSamplesOnly  = false;
  vars->fftDataType    = FX_NONE;
  vars->taperType      = mod_fft::TAPER_NONE;
  vars->taperLengthInSamples = 20;

  vars->numSamplesIn   = shdr->numSamples;
  vars->sampleIntIn    = shdr->sampleInt;
  
  vars->fft = NULL;

  std::string text;

  if( param->exists("direction") ) {
    param->getString( "direction", &text );
    text = toLowerCase( text );

    if( text.compare("forward") == 0 ) {
      vars->direction = mod_fft::FORWARD;
    }
    else if( text.compare("inverse") == 0 ) {
      vars->direction = mod_fft::INVERSE;
    }
    else {
      writer->line("ERROR: Unknown 'direction' option: %s", text.c_str());
      env->addError();
    }
  }

  // Ignore domain stored in the superheader if it conflicts with the user parameters.
  if( param->exists("override") ) {
    bool overrideFX = false;
    param->getString("override", &text);
    text = toLowerCase( text );
    if( text.compare("xt") == 0 ) {
      shdr->domain = DOMAIN_XT;
    }
    else if( text.compare("xd") == 0 ) {
      shdr->domain = DOMAIN_XD;
    }
    else if( text.compare("fx")  == 0 ) {
      overrideFX = true;
      shdr->domain = DOMAIN_FX;
    }
    else if( text.compare("no")  == 0 ) {
      // Nothing
    }
    else {
      writer->line("ERROR:Unknown 'override' option: %s", text.c_str());
      env->addError();
    }
    if( overrideFX ) {
      param->getString("override", &text, 1);
      if( text.compare("amp_phase") == 0 ) {
        shdr->fftDataType = FX_AMP_PHASE;
      }
      else if( text.compare("real_imag") == 0 ) {
        shdr->fftDataType = FX_REAL_IMAG;
      }
      else if( text.compare("amp") == 0 ) {
        shdr->fftDataType = FX_AMP;
      }
      else if( text.compare("psd") == 0 ) {
        shdr->fftDataType = FX_PSD;
      }
      else {
        writer->line("ERROR: Unknown option '%s'", text.c_str() );
        env->addError();
      }
    }
  }

  // Normalize the output.
  if( param->exists("norm") ) {
    param->getString("norm", &text);
    if( text.compare("none") == 0 ) {
      vars->normMethod = mod_fft::NORM_NONE;
    }
    else if( text.compare("standard") == 0 ) {
      vars->normMethod = mod_fft::NORM_STANDARD;
    }
    else if( text.compare("time_average") == 0 ) {
      vars->normMethod = mod_fft::NORM_TIME_AVERAGE;
    }
    else if( text.compare("time_average_sqrt") == 0 ) {
      vars->normMethod = mod_fft::NORM_TIME_AVERAGE_SQRT;
    }
    else {
      writer->line("Error: Unknown 'norm' option: %s", text.c_str());
      env->addError();
    }
    if( param->getNumValues("norm") > 1 ) {
      param->getString("norm", &text, 1);
      if( text.compare("all") == 0 ) {
        vars->nonZeroSamplesOnly = false;
      }
      else if( text.compare("nonzero") == 0 ) {
        vars->nonZeroSamplesOnly = true;
        if( vars->normMethod != mod_fft::NORM_TIME_AVERAGE || vars->normMethod != mod_fft::NORM_TIME_AVERAGE_SQRT ) {
          writer->warning("Ignored user-specified 'non-zero sample usage': Only applies to time averaging normalization.");
          vars->nonZeroSamplesOnly = false;
        }
      }
      else {
        writer->line("Error: Unknown 'norm' option: %s", text.c_str());
        env->addError();
      }
    }
  }
  if( param->exists("sidedness") ) {
    param->getString("sidedness", &text, 1);
    if( text.compare("double") == 0 ) {
      vars->apply2x_doubleSided = true;
    }
    else if( text.compare("single") == 0 ) {
      vars->apply2x_doubleSided = false;
    }
    else {
      writer->line("Error: Unknown 'sidedness' option: %s", text.c_str());
      env->addError();
    }
  }
  
  //--------------------------------------------------------------------------------
  // Taper the input
  if( param->exists("taper_type") ) {

    if( vars->direction == mod_fft::INVERSE ) {
      writer->warning("User parameter 'taper_type' will be ignored for INVERSE FFT transform");
    }
    else {
      param->getString("taper_type", &text);
      text = toLowerCase( text );

      if( text.compare("none") == 0 ) {
        vars->taperType = mod_fft::TAPER_NONE;
      }
      else if( text.compare("cos") == 0 ) {
        vars->taperType = mod_fft::TAPER_COSINE;
      }
      else if( text.compare("hanning") == 0 ) {
        vars->taperType = mod_fft::TAPER_HANNING;
        vars->taperLengthInSamples = shdr->numSamples/2;
      }
      else if( text.compare("blackman") == 0 ) {
        vars->taperType = mod_fft::TAPER_BLACKMAN;
        vars->taperLengthInSamples = shdr->numSamples/2;
      }
      else {
        writer->line("ERROR:Unknown 'taper_type' option: %s", text.c_str());
        env->addError();
      }
    }
  }
  if( param->exists("taper_len") ) {
    if( vars->direction == mod_fft::INVERSE ) {
      writer->warning("User parameter 'taper_len' will be ignored for INVERSE FFT transform");
    }
    else {
      if( vars->taperType == mod_fft::TAPER_HANNING ) {
        writer->line("ERROR:Cannot specify 'taper_len' for Hanning taper, taper length is fixed to half the trace length");
        env->addError();
      }
      else {
        float taperLength;
        param->getFloat("taper_len", &taperLength);
        vars->taperLengthInSamples = (int)round( taperLength / shdr->sampleInt );
      }
    }
  }
  if( vars->taperType != mod_fft::TAPER_NONE && vars->taperLengthInSamples > vars->numSamplesIn ) {
    writer->error("Specified taper is too long. Number of samples in input trace: %d, taper length: %d", vars->numSamplesIn, vars->taperLengthInSamples);
  }

  // Output type of the FORWARD transform data
  bool outputEvenPSD = false;
  if( param->exists("output") ) {
    if( vars->direction == mod_fft::INVERSE ) {
      writer->warning("User parameter 'output' will be ignored for INVERSE FFT transform");
    }
    else {
      param->getString("output", &text);
      if( text.compare("amp_phase") == 0 ) {
        vars->fftDataType = FX_AMP_PHASE;
      }
      else if( text.compare("amp") == 0 ) {
        vars->fftDataType = FX_AMP;
      }
      else if( text.compare("phase") == 0 ) {
        vars->fftDataType = FX_PHASE;
      }
      else if( text.compare("psd") == 0 ) {
        vars->fftDataType = FX_PSD;
      }
      else if( text.compare("psd_even") == 0 ) {
        vars->fftDataType = FX_PSD;
        outputEvenPSD = true;
      }
      else if( text.compare("real_imag") == 0 ) {
        vars->fftDataType = FX_REAL_IMAG;
      }
      else {
        writer->line("ERROR:Unknown 'output' option: %s", text.c_str());
        env->addError();
      }
    }

  }
  else if( vars->direction == mod_fft::FORWARD ) {
    writer->line("ERROR:User parameter 'output' is required for forward FFT transform");
    env->addError();
  }

  int numSamples_time = ( vars->direction == mod_fft::FORWARD ) ? vars->numSamplesIn : shdr->numSamplesXT;
  if( param->exists("nsamples_fft") ) {
    int numSamples_fft;
    param->getInt("nsamples_fft", &numSamples_fft);
    if( numSamples_fft < shdr->numSamples ) writer->error("Number of samples in FFT (%d) must exceed number of samples in input data (%d)", numSamples_fft, shdr->numSamples);
    vars->fft = new cseis_geolib::csFFT( numSamples_time, numSamples_fft );
  }
  else {
    vars->fft = new cseis_geolib::csFFT( numSamples_time );
  }
  //--------------------------------------------------------------------------------
  // Set normalization value
  //
  if( vars->normMethod == mod_fft::NORM_STANDARD ) {
    if( vars->direction == mod_fft::FORWARD ) {
      vars->normValue = shdr->sampleInt / 1000.0f; // 1/1000 to convert from [ms] to [s]
    }
    else { // INVERSE
      vars->normValue = 1000.0f / shdr->sampleIntXT; // x1000 to convert from [ms] to [s]
    }
  }
  else if( vars->normMethod == mod_fft::NORM_TIME_AVERAGE ) {
    if( vars->direction == mod_fft::FORWARD ) {
      vars->normValue = 1.0f / (float)shdr->numSamples; // Number of samples in X-T domain
    }
    else { // INVERSE
      vars->normValue = (float)shdr->numSamplesXT; // Number of samples in X-T domain
    }
  }
  else if( vars->normMethod == mod_fft::NORM_TIME_AVERAGE_SQRT ) {
    if( vars->direction == mod_fft::FORWARD ) {
      vars->normValue = 1.0f / (float)sqrt( shdr->numSamples ); // Number of samples in X-T domain
    }
    else { // INVERSE
      vars->normValue = (float)sqrt( shdr->numSamplesXT ); // Number of samples in X-T domain
    }
  }
  
  // FORWARD FFT
  if( vars->direction == mod_fft::FORWARD ) {
    writer->line("FFT FORWARD transform." );

    if( shdr->domain != DOMAIN_XT && shdr->domain != DOMAIN_XD ) {
      writer->line("ERROR: Current trace is not in XT (or XD) domain. FFT forward transform not possible. Actual domain: %s (%d)", csGeolibUtils::domain2Text(shdr->domain) );
      env->addError();
    }

    int numSamplesOut = 0;
    if( vars->fftDataType == FX_AMP_PHASE ) {
      numSamplesOut = 2 * vars->fft->numFreqValues(); 
    }
    else if( vars->fftDataType == FX_AMP || vars->fftDataType == FX_PSD || vars->fftDataType == FX_PHASE ) {
      numSamplesOut = vars->fft->numFreqValues(); 
      if( outputEvenPSD ) numSamplesOut = vars->fft->numFFTValues()/2;
    }
    else if( vars->fftDataType == FX_REAL_IMAG ) {
      numSamplesOut = 2 * vars->fft->numFreqValues(); 
    }
    else {
      writer->error("Unsupported FFT data type: ", vars->fftDataType);
    }

    // Evaluate frequency info
    float freqNyquist = 1.0f / ( 2.0f * ( vars->sampleIntIn / 1000.0f ) );
    shdr->sampleInt = freqNyquist / (float)(vars->fft->numFreqValues()-1);
    writer->line("Input data sample rate = %.6f   Nyquist = %.6f    Delta f = %.8f", vars->sampleIntIn, freqNyquist, shdr->sampleInt );

    shdr->numSamplesXT = vars->numSamplesIn;
    shdr->sampleIntXT  = vars->sampleIntIn;
    shdr->numSamples   = numSamplesOut;

    shdr->domain      = DOMAIN_FX;
    shdr->fftDataType = vars->fftDataType;
  }
  else if( vars->direction == mod_fft::INVERSE ) {
    writer->line("FFT INVERSE transform." );

    if( shdr->domain != DOMAIN_FX ) {
      writer->line("ERROR:Current trace is not in FX domain. FFT inverse transform not possible. Actual domain: %s", csGeolibUtils::domain2Text(shdr->domain) );
      env->addError();
    }
    else if( shdr->fftDataType != FX_REAL_IMAG && shdr->fftDataType != FX_AMP_PHASE ) {
      writer->line("ERROR: Cannot perform inverse transform due to lack of information: Data is neither of type AMP_PHASE or REAL_IMAG. Actual type defined in super header: %s",   csGeolibUtils::fftDataType2Text( shdr->fftDataType ) );
      env->addError();              
    }
    vars->fftDataType  = shdr->fftDataType;

    shdr->numSamples   = shdr->numSamplesXT;
    shdr->sampleInt    = shdr->sampleIntXT;
    shdr->numSamplesXT = 0;
    shdr->sampleIntXT  = 0;

    shdr->domain      = DOMAIN_XT;  // CHANGE: May need to set to XD domain...
    shdr->fftDataType = FX_NONE; 
  }

  writer->line("Number of input samples:  %d", vars->numSamplesIn);
  writer->line("Number of output samples: %d", shdr->numSamples);
  writer->line("Number of time samples:   %d", numSamples_time);
  writer->line("Number of FFT values:     %d", vars->fft->numFFTValues());
  writer->line("Number of frequencies:    %d", vars->fft->numFreqValues());
  writer->line("Normalization value:      %f", vars->normValue);
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_fft_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);
  float* samples = trace->getTraceSamples();
  float normValueCurrent = vars->normValue;
  int numFreq = vars->fft->numFreqValues();

  // Forward transform
  if( vars->direction == mod_fft::FORWARD ) {

    if( vars->nonZeroSamplesOnly ) {
      int numSamplesNonZero = vars->numSamplesIn;
      for( int i = 0; i < vars->numSamplesIn; i++ ) {
        if( samples[i] == 0.0 ) numSamplesNonZero -= 1;
      }
      if( numSamplesNonZero == 0 ) numSamplesNonZero = 1;
      if( vars->taperType == mod_fft::TAPER_NONE ) {
        normValueCurrent = 1.0f / (float)numSamplesNonZero;  // Number of non-zero samples in X-T domain
      }
      else if( vars->taperType == mod_fft::TAPER_COSINE || vars->taperType == mod_fft::TAPER_HANNING ) {
        double weightSum = 0.0;
        int weightCounter = 0;
        for( int i = 0; i < vars->taperLengthInSamples; i++ ) {
          if( samples[i] == 0 ) continue;
          double scalar = cos( M_PI_2 * (float)(vars->taperLengthInSamples-i)/(float)vars->taperLengthInSamples );
          weightSum += scalar;
          weightCounter += 1;
        }
        for( int i = vars->numSamplesIn-vars->taperLengthInSamples; i < vars->numSamplesIn; i++ ) {
          if( samples[i] == 0 ) continue;
          double scalar = cos( M_PI_2 * (float)(vars->taperLengthInSamples-vars->numSamplesIn+i+1)/(float)vars->taperLengthInSamples );
          weightSum += scalar;
          weightCounter += 1;
        }

        normValueCurrent = 1.0f / (float)( numSamplesNonZero + 2.0 * ( weightSum - (double)weightCounter ) );
      }
      else if( vars->taperType == mod_fft::TAPER_BLACKMAN ) {
        double weightSum  = 0.0;
        int weightCounter = 0;
        double alpha = 0.16;
        double a0 = 0.5 * (1.0 - alpha);
        double a1 = 0.5;
        double a2 = 0.5 * alpha;
        for( int i = 0; i < vars->numSamplesIn; i++ ) {
          double piFactor = (2.0 * M_PI) * (double)i / (double)(vars->numSamplesIn - 1);
          float weight = (float)( a0 - a1*cos( piFactor ) + a2*cos( 2.0 * piFactor ) );
          weightSum += weight;
          weightCounter += 1;
        }
        if( weightCounter < 1 ) weightSum = 1.0;
        normValueCurrent = 1.0f / (float)weightSum;
      }
    } // END: NORM_...NONZERO
    
    // Apply taper to input trace
    if( vars->taperType == mod_fft::TAPER_COSINE || vars->taperType == mod_fft::TAPER_HANNING ) {
      for( int i = 0; i < vars->taperLengthInSamples; i++ ) {
        if( samples[i] == 0 ) continue;
        double scalar = cos( M_PI_2 * (float)(vars->taperLengthInSamples-i)/(float)vars->taperLengthInSamples );
        samples[i] *= scalar;
        samples[vars->numSamplesIn-1-i] *= scalar;
      }
    }
    else if( vars->taperType == mod_fft::TAPER_BLACKMAN ) {
      double alpha = 0.16;
      double a0 = 0.5 * (1.0 - alpha);
      double a1 = 0.5;
      double a2 = 0.5 * alpha;
      for( int i = 0; i < vars->numSamplesIn; i++ ) {
        double piFactor = (2.0 * M_PI) * (double)i / (double)(vars->numSamplesIn - 1);
        float weight = (float)( a0 - a1*cos( piFactor ) + a2*cos( 2.0 * piFactor ) );
        samples[i] *= weight;
      }
    }
    
    vars->fft->forwardTransform( samples, samples, vars->fftDataType );

    if( vars->normMethod != mod_fft::NORM_NONE ) {
      if( vars->fftDataType == FX_PSD ) {
        normValueCurrent = pow( normValueCurrent, 2 );
      }
      if( vars->apply2x_doubleSided ) {
        normValueCurrent *= 2.0f;
      }
    }
    

    //----------------------------------------------------------------------
    if( vars->normMethod != mod_fft::NORM_NONE ) {
      // FORWARD - PSD
      if( vars->fftDataType == FX_PSD ) {
        for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
          samples[isamp] *= normValueCurrent;
        }
      }
      // FORWARD - AMP_PHASE, AMP
      else if( vars->fftDataType == FX_AMP_PHASE || vars->fftDataType == FX_AMP ) {
        for( int isamp = 0; isamp < numFreq; isamp++ ) {
          samples[isamp] *= normValueCurrent;
        }
      }
      // FORWARD - REAL_IMAG
      // FORWARD - FX_PHASE
      else {
        // Nothing to be done
      }
    } // END normMethod != NORM_NONE
  }
  //--------------------------------------------------------------------------------
  else { // Inverse transform.  if( vars->direction == mod_fft::INVERSE ) {
    if( vars->fftDataType == FX_AMP_PHASE && vars->normMethod != mod_fft::NORM_NONE ) {
      if( vars->apply2x_doubleSided ) {
        normValueCurrent /= 2.0f;
      }
      int numSamplesNormalise = vars->numSamplesIn/2;
      for( int is = 0; is < numSamplesNormalise; is++ ) {
        samples[is] *= normValueCurrent;
      }
    }

    vars->fft->inverseTransform( samples, samples, vars->fftDataType );
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_fft_( csParamDef* pdef ) {
  pdef->setModule( "FFT", "FFT transform", "Transform input data from x-t to f-x (x-w) domain (forward) or vice versa (inverse)");
  pdef->setVersion(1,0);
  pdef->addDoc("Depending on whether the FFTW library was available at build time, the output spectrum will have a slightly different sample step in the frequency domain.");
  pdef->addDoc("  a) Normalization method 'standard' (= default method)");
  pdef->addDoc("  - For active source signal, this normalization is invariant of trace length and sample interval (assumes that data contains same, noise-free, signal)");
  pdef->addDoc("  - Exploration geophysicist's amplitude unit:  'dB re (uPa/Hz)' (amplitude)  or  'dB re (uPa/Hz)^2' (energy)");
  pdef->addDoc("  - ...to adher to the SEG standard for a source signal (dB re (uPa/Hz) @ 1m), must compensate for spherical spreading to move receiver to 1m distance from source");
  pdef->addDoc("  b) Normalization method 'time_average'");
  pdef->addDoc("  - Ambient noise: Use power = energy / time (instead of energy) --> divide above normalization by time");
  pdef->addDoc("  - For (stationary) signal or noise, this normalization is invariant of trace length and sample interval (assumes that data consists of time-invariant continuous noise)");
  pdef->addDoc("  - Acousticians standard amplitude unit:  'dB re (uPa/sqrt(Hz))' (amplitude)  or  'dB re (uPa/sqrt(Hz))^2' (energy)");
  pdef->addDoc("  c) Normalization method 'none'");
  pdef->addDoc("  - Keep raw fft values. No normalization.");
  
  pdef->addParam( "direction", "Direction of transform", NUM_VALUES_FIXED );
  pdef->addValue( "forward", VALTYPE_OPTION );
  pdef->addOption( "forward", "Forward transform from x-t to f-x(x-w)" );
  pdef->addOption( "inverse", "Inverse transform from f-x(x-w) to x-t", "Inverse transform only works if a forward transform was applied earlier" );

  pdef->addParam( "taper_type", "Type of taper to apply to input trace", NUM_VALUES_FIXED );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "none", "Do not apply any taper to input trace" );
  pdef->addOption( "cos", "Apply cosine taper to input trace" );
  pdef->addOption( "hanning", "Apply 'Hanning' cosine taper to input trace. Taper length is 1/2 trace length" );
  pdef->addOption( "blackman", "Apply 'Blackman' taper (alpha=0.16) to input trace" );

  pdef->addParam( "taper_len", "Taper length [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "20", VALTYPE_NUMBER, "Taper length in [ms]" );

  pdef->addParam( "output", "Output options for forward transform.", NUM_VALUES_FIXED );
  pdef->addValue( "amp_phase", VALTYPE_OPTION );
  pdef->addOption( "amp_phase", "Output amplitude and phase spectrum concatenated into one trace, amplitudes at the beginning of the trace, phase at the end" );
  pdef->addOption( "amp", "Output amplitude spectrum", "For input data unit [U], output unit is [U/Hz]" );
  pdef->addOption( "real_imag", "Output real and imaginary values concatenated into one trace, real values at the beginning of the trace, imaginary values at the end" );
  pdef->addOption( "psd", "Output power spectral density", "For input data unit [U], output unit is [U/Hz]^2" );
  pdef->addOption( "psd_even", "Output PSD spectrum. Omit value at Nyquist frequency, i.e. output 2^N samples" );
  pdef->addOption( "phase", "Output phase spectrum" );

  pdef->addParam( "nsamples_fft", "For FORWARD transform, override the number of values in FFT transform.", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER );

  pdef->addParam( "override", "Override domain specified in superheader.", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Acknowledge domain found in super header" );
  pdef->addOption( "fx", "Override domain found in super header. Set to FX." );
  pdef->addOption( "xt", "Override domain found in super header. Set to XT (time)." );
  pdef->addOption( "xd", "Override domain found in super header. Set to XD (depth)." );
  pdef->addValue( "amp_phase", VALTYPE_OPTION, "Override data type" );
  pdef->addOption( "amp_phase", "Set FFT data type to amplitude/phase spectrum. Input data consists of amplitude and phase spectra" );
  pdef->addOption( "real_imag", "Set FFT data type to real/imaginary spectrum" );

  pdef->addParam( "norm", "Normalization method for amplitude and PSD spectra", NUM_VALUES_VARIABLE );
  pdef->addValue( "standard", VALTYPE_OPTION );
  pdef->addOption( "standard", "Standard normalization (=multiply raw fft result by <sample interval>. Appropriate for active source signal. In case of strong short-termed signal, keeps amplitude spectrum invariant of trace length and sample interval" );
  pdef->addOption( "none", "No normalization (=do not apply any factors to raw fft result)" );
  pdef->addOption( "time_average", "Normalize by 1/tim (=divide raw fft result by <number of samples>). Appropriate for (statinary) signal or noise. In case of strong ambient noise, keeps amplitude spectrum invariant of trace length and sample interval" );
  pdef->addOption( "time_average_sqrt", "Same as 'time_average' but taking square root of number of samples: This can work better when comparing random noise in windows of different length" );
  pdef->addValue( "all", VALTYPE_OPTION, "How to treat non-zero values" );
  pdef->addOption( "all", "Applies to time averaging only: Count all samples");
  pdef->addOption( "nonzero", "Applies to time averaging only: Count non-zero samples only");

  pdef->addParam( "sidedness", "Single- or double-sided spectrum?", NUM_VALUES_FIXED );
  pdef->addValue( "double", VALTYPE_OPTION );
  pdef->addOption( "double", "Scale output amplitude/PSD spectra by 2x: Compensate spectra to account for positive and negative frequencies");
  pdef->addOption( "single", "Do not apply additional scaling to output amplitude/PSD spectra: Output single-sided spectra, positive frequencies only");
}

//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_fft_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_fft::VariableStruct* vars = reinterpret_cast<mod_fft::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_fft_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_fft::VariableStruct* vars = reinterpret_cast<mod_fft::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->fft != NULL ) {
    delete vars->fft;
    vars->fft = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_fft_( csParamDef* pdef ) {
  params_mod_fft_( pdef );
}
extern "C" void _init_mod_fft_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_fft_( param, env, writer );
}
extern "C" bool _start_exec_mod_fft_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_fft_( env, writer );
}
extern "C" void _exec_mod_fft_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_fft_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_fft_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_fft_( env, writer );
}

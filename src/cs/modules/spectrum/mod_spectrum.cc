/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csGeolibUtils.h"
#include "csFFT.h"
#include "csTaper.h"
#include <cmath>
#include <cstring>
#include <limits>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPECTRUM
 *
 * @author Felipe Punto
 * @date   2015
 */
namespace mod_spectrum {
  struct VariableStruct {
    float* bufferAuto;
    int numSamplesWin;
    int numSamplesLag;
    int numSamplesAuto;
    bool isTimeDomain;
    int scaleTest;

    int normFlag;       // Normalize data: Yes or no?
    int normFunction;
    double normScalar;   // Normalization scalar to be applied before inverse FFT.
                         // This scalar accounts for zeros that may have been padded to input trace
    int numSamplesToOutput;   // Number of samples output as equivalent # of floats (ex: 1 complex = 2 floats)
    int fftDataType;        // Datatype of the transformed data after FORWARD. 

    int taperType;            // Type of taper tp apply to input 
    int numSamplesTaper; // Taper length in number of samples (from 0 to 1)
    int autoTaperNumSamples;
    int autoTaperStartSample;
    bool isAutoTaper;

    int numSamplesToFFT;      // Number of samples after padding *before* FOWARD
    int nFreqToFFT;         // Number of samples in frequency domain (complex == 2*float)

    int option; // Spectrum computation option
    int lag_samples;  // Auto-correlation lag

    int startSamp; // Start sample index of input window
    int endSamp;   // End sample index of input window

    cseis_geolib::csFFT* fftObj;
   };
  static const int OPTION_ACOR = 34;
  static const int OPTION_XCOR = 35;
}
using mod_spectrum::VariableStruct;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spectrum_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader*  shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->fftDataType = FX_AMP;
  vars->taperType   = cseis_geolib::csTaper::NONE;
  vars->normFlag    = cseis_geolib::csTaper::NORM_NO;
  vars->normFunction = cseis_geolib::csTaper::NORM_NSAMP;
  vars->option       = mod_spectrum::OPTION_ACOR;
  vars->bufferAuto  = NULL;
  vars->numSamplesWin = shdr->numSamples;
  vars->numSamplesLag = vars->numSamplesWin - 1;
  vars->numSamplesAuto = 2*vars->numSamplesLag + 1;
  vars->numSamplesTaper = (int)( 20.0f / shdr->sampleInt );
  vars->startSamp       = 0;
  vars->endSamp         = shdr->numSamples-1;
  vars->isTimeDomain = true;
  vars->scaleTest = 1;
  vars->isAutoTaper = false;
  vars->autoTaperStartSample = 0;
  vars->autoTaperNumSamples = 0;

  std::string text;

  vars->normScalar     = 1.0;

  if( param->exists("norm") ) {
    param->getString("norm", &text);
    text = toLowerCase( text );
    if( text.compare("yes") == 0 ) {
      vars->normFlag = cseis_geolib::csTaper::NORM_YES;
    }
    else if( text.compare("no") == 0 ) {
      vars->normFlag = cseis_geolib::csTaper::NORM_NO;
    }
    else {
      writer->line("ERROR:Unknown 'norm' option: %s", text.c_str());
      env->addError();
    }
    if( param->getNumValues("norm") > 1 ) {
      param->getString("norm", &text, 1);
      if( !text.compare("nsamp") ) {
        vars->normFunction = cseis_geolib::csTaper::NORM_NSAMP;
      }
      else if( !text.compare("sqrt") ) {
        vars->normFunction = cseis_geolib::csTaper::NORM_SQRT;
      }
      else {
        writer->error("Unknown option: %s", text.c_str());
      }
    }
  }

  if( param->exists("option") ) {
    param->getString("option", &text);
    if( text.compare("auto") == 0 ) {
      vars->option = mod_spectrum::OPTION_ACOR;
    }
    else if( text.compare("cross") == 0 ) {
      vars->option = mod_spectrum::OPTION_ACOR;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }


  if( param->exists("output") ) {
    param->getString("output", &text);
    if( text.compare("amp") == 0 ) {
      vars->fftDataType = FX_AMP;
    }
    else if( text.compare("psd") == 0 ) {
      vars->fftDataType = FX_PSD;
    }
    else {
      writer->line("ERROR:Unknown 'output' option: %s", text.c_str());
      env->addError();
    }
  }
  //--------------------------------------------------------------------------------
  // Set up two-sided auto-correlation

  if( param->exists( "domain" ) ) {
    param->getString( "domain", &text );
    if( !text.compare( "sample" ) ) {
      vars->isTimeDomain = false;
    }
    else if( !text.compare( "time" ) ) {
      vars->isTimeDomain = true;
    }
    else {
      writer->line("Domain not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }
  float startTime = 0.0;
  float endTime   = shdr->sampleInt * (float)vars->endSamp;
  
  if( param->exists("start") ) {
    string textStart;
    param->getString( "start", &textStart );
    if( isDigit( textStart.at(0) ) ) {
      float start = (float)atof( textStart.c_str() );
      if( vars->isTimeDomain ) {
        startTime = start;
        vars->startSamp = (int)round( startTime / shdr->sampleInt );
      }
      else {
        vars->startSamp = (int)round( start ) - 1;  // User provides sample index starting at 1. Convert to start index = 0
        startTime = (float)vars->startSamp * shdr->sampleInt;
      }
    }
    else {
      writer->error("Start time/sample not a number: %s", text.c_str());
      //      vars->hdrId_start   = hdef->headerIndex(textStart);
      //   vars->hdrType_start = hdef->headerType(textStart);
    }
  }
  if( param->exists("end") ) {
    string textEnd;
    param->getString( "end", &textEnd );
    if( isDigit( textEnd.at(0) ) ) {
      float end = (float)atof( textEnd.c_str() );
      if( vars->isTimeDomain ) {
        endTime = end;
        vars->endSamp = (int)round( endTime / shdr->sampleInt );
      }
      else {
        vars->endSamp = (int)end - 1;  // User provides sample index starting at 1
        endTime = (float)vars->endSamp * shdr->sampleInt;
      }
    }
    else {
      writer->error("End time/sample not a number: %s", text.c_str());
      //      vars->hdrId_end   = hdef->headerIndex(textEnd);
      //  vars->hdrType_end = hdef->headerType(textEnd);
    }
  }

  if( vars->isTimeDomain ) {
    if( startTime < 0.0 ) writer->error("Start time (%f) needs to be greater or equal to 0.0.", startTime);
    if( startTime > endTime ) writer->error("Start time (%f) needs to be smaller than end time (%f).", startTime, endTime);
    if( endTime >= (float)shdr->numSamples*shdr->sampleInt ) writer->error("Specified end time (%fms) is greater than trace length (%fms).", endTime, (float)shdr->numSamples*shdr->sampleInt);
  }
  else {
    // NOTE: User input is '1' for first sample. Internally, '0' is used
    if( vars->startSamp < 0 ) writer->error("Start sample (%d) needs to be greater or equal to 1.", vars->startSamp+1);
    if( vars->startSamp > vars->endSamp ) writer->error("Start sample (%d) needs to be smaller than end sample (%d).", vars->startSamp+1, vars->endSamp+1);
    if( vars->endSamp > shdr->numSamples-1 ) writer->error("End sample/time greater than number samples in trace (%d).", vars->endSamp+1, shdr->numSamples);
  }

  vars->numSamplesWin  = vars->endSamp - vars->startSamp + 1;
  vars->numSamplesLag  = vars->numSamplesWin - 1;
  vars->numSamplesAuto = 2*vars->numSamplesLag + 1;
  vars->bufferAuto = new float[vars->numSamplesAuto];
  if( vars->normFlag != cseis_geolib::csTaper::NORM_NO ) {
    vars->normScalar = 1.0 / vars->numSamplesAuto;  // Scale by inverse of number of samples input to FFT
  }
  vars->fftObj = new cseis_geolib::csFFT( vars->numSamplesAuto );

  // Taper the input
  if( param->exists("taper_type") ) {
    param->getString("taper_type", &text);
    text = toLowerCase( text );

    if( text.compare("none") == 0 ) {
      vars->taperType = cseis_geolib::csTaper::NONE;
    }
    else if( text.compare("cos") == 0 ) {
      vars->taperType = cseis_geolib::csTaper::COSINE;
    }
    else if( text.compare("hanning") == 0 ) {
      vars->taperType = cseis_geolib::csTaper::HANNING;
      vars->numSamplesTaper = vars->numSamplesAuto/2;
    }
    else if( text.compare("blackman") == 0 ) {
      vars->taperType = cseis_geolib::csTaper::BLACKMAN;
      vars->numSamplesTaper = vars->numSamplesAuto/2;
    }
    else {
      writer->line("ERROR:Unknown 'taper_type' option: %s", text.c_str());
      env->addError();
    }
  }

  // Taper applied to input data
  if( param->exists("taper_len") ) {
    if( vars->taperType == cseis_geolib::csTaper::HANNING ) {
      writer->line("ERROR: Do not specify 'taper_len' for Hanning taper: Taper length is fixed to half the trace length");
      env->addError();
    }
    else {
      float taperLength;
      param->getFloat("taper_len", &taperLength);
      vars->numSamplesTaper = (int)round( taperLength / shdr->sampleInt );
    }
  }
  if( vars->numSamplesTaper >= vars->numSamplesWin ) {
    writer->error("Specified input taper is too long. Number of samples in windowed input trace: %d, and in taper: %d", vars->numSamplesWin, vars->numSamplesTaper);
  }

  // Taper applied to auto-correlation trace before FFT
  if( param->exists("taper_auto") ) {
    float autoTaperStart_ms;
    float autoTaperWidth_ms;
    param->getFloat("taper_auto", &autoTaperStart_ms, 0);
    param->getFloat("taper_auto", &autoTaperWidth_ms, 1);
    vars->autoTaperStartSample = (int)round( autoTaperStart_ms / shdr->sampleInt );
    vars->autoTaperNumSamples  = (int)round( autoTaperWidth_ms / shdr->sampleInt );
    if( (vars->autoTaperStartSample + vars->autoTaperNumSamples) > (int)((vars->numSamplesAuto)/2) ) {
      writer->error("Auto-correlation taper is too long.\nAuto-correlation half length in samples: %d.\nTaper start sample and length in samples: %d %d (sum=%d).",
                 (int)((vars->numSamplesAuto)/2), vars->autoTaperStartSample, vars->autoTaperNumSamples, vars->autoTaperStartSample+vars->autoTaperNumSamples );
    }
    vars->isAutoTaper = true;
  }
  //--------------------------------------------------------------------------------
  // Set up FFT
  vars->numSamplesToFFT = vars->fftObj->numFFTValues();
  vars->nFreqToFFT      = vars->fftObj->numFreqValues();

  vars->numSamplesToOutput = vars->nFreqToFFT; 

  // Save time-domain information
  shdr->numSamplesXT = vars->numSamplesAuto;
  shdr->sampleIntXT  = shdr->sampleInt;

  float nyquist = (float)( 1.0/(2.0*(shdr->sampleIntXT/1000.0)) );
  
  // Reset samples & sample rate for output
  shdr->numSamples  = vars->numSamplesToOutput;
  shdr->sampleInt   = nyquist/(float)(vars->nFreqToFFT-1);  // New sample interval in frequency domain
  shdr->domain      = DOMAIN_FX;
  shdr->fftDataType = vars->fftDataType;


  writer->line("  Input data length:          %d  samples", shdr->numSamples );
  writer->line("  Input data window length:   %d  samples", vars->numSamplesWin );
  writer->line("  Auto-correlation lag:       %d  samples", vars->numSamplesLag);
  writer->line("  Auto-correlation length:    %d  samples", vars->numSamplesAuto);
  writer->line("  FFT:");
  writer->line("  Number of frequencies:      %d", vars->nFreqToFFT );
  writer->line("  Number of output samples:   %d", vars->numSamplesToOutput );
  writer->line("  Input sample rate = %f -> Nyquist is %f @ delta f = %f.", shdr->sampleIntXT, nyquist, shdr->sampleInt );
  writer->line("");
  shdr->dump( writer->getFile() );
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spectrum_(
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

  float normScalarCurrent = (float)vars->normScalar;

  //--------------------------------------------------------------------------------
  // (0) Apply taper to input trace
  if( vars->taperType != cseis_geolib::csTaper::NONE ) {
  }

  //--------------------------------------------------------------------------------
  // (1) Compute auto-corrlelation
  //
  float* tracePtr = trace->getTraceSamples();
  float* bufferPtr;
  bufferPtr = &vars->bufferAuto[vars->numSamplesLag];  // Set pointer to zero-lag sample
  compute_onesided_auto_correlation(
                                    &tracePtr[vars->startSamp],
                                    vars->numSamplesWin,
                                    bufferPtr,
                                    vars->numSamplesLag,
                                    false );

  // Copy one-sided auto-correlation to other side:
  for( int isamp = vars->numSamplesLag-1; isamp >= 0; isamp-- ) {
    vars->bufferAuto[isamp] = vars->bufferAuto[ vars->numSamplesAuto - 1 - isamp ];
  }

  //--------------------------------------------------------------------------------
  // (2) Compute FFT spectrum from auto-correlation trace
  //
  float* samples = trace->getTraceSamples();

  // (2.a) Apply taper to auto-correlation trace
  if( vars->isAutoTaper ) {
    int startSampTaper = vars->numSamplesLag - (vars->autoTaperStartSample + vars->autoTaperNumSamples);
    int numSamplesTemp = vars->numSamplesAuto - 2*startSampTaper;
    for( int isamp = 0; isamp < startSampTaper; isamp++ ) {
      vars->bufferAuto[isamp] = 0;
      vars->bufferAuto[vars->numSamplesAuto-isamp-1] = 0;
    }
    cseis_geolib::csTaper taper;
    taper.applyTaper( &vars->bufferAuto[startSampTaper], numSamplesTemp, vars->autoTaperNumSamples, csTaper::COSINE, vars->normFlag, normScalarCurrent );
  }

  vars->fftObj->forwardTransform( vars->bufferAuto, samples, vars->fftDataType );

  for ( int isamp=0; isamp < shdr->numSamples; isamp++) {
    samples[isamp] = shdr->sampleIntXT * sqrt(samples[isamp]);   // Double-sqrt application to revert back to amplitude spectrum. Normalize by sample interval.
  }

  if( vars->normFlag != cseis_geolib::csTaper::NORM_NO ) {
    if( vars->normFunction == cseis_geolib::csTaper::NORM_SQRT ) {
      normScalarCurrent = sqrt(normScalarCurrent);
    }
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] *= normScalarCurrent;
    }
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_spectrum_( csParamDef* pdef ) {
  pdef->setModule( "SPECTRUM", "Compute frequency spectrum", "For straight-forward FFT transform, use $FFT");
  pdef->setVersion(1,0);
 
  pdef->addParam( "option", "Option for computing frequency spectrum.", NUM_VALUES_FIXED );
  pdef->addValue( "auto", VALTYPE_OPTION );
  pdef->addOption( "auto", "Compute frequency spectrum from auto-correlation" );
  //  pdef->addOption( "cross", "Compute frequency spectrum from cross-correlation of adjacent traces" );

  pdef->addParam( "domain", "Is correlation window given in time or in samples?", NUM_VALUES_FIXED );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Window is specified in time [ms]" );
  pdef->addOption( "sample", "Window is specified in samples (1 for first sample)" );

  pdef->addParam( "start", "Start time/sample", NUM_VALUES_FIXED, "Start time or sample, this depends on the 'domain' setting" );
  pdef->addValue( "0", VALTYPE_STRING, "Start time/sample of cross-correlation window", "Alternatively, name of trace header containing start time/sample" );

  pdef->addParam( "end", "End time/sample", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "End time/sample of cross-correlation window", "Alternatively, name of trace header containing end time/sample" );

  pdef->addParam( "taper_type", "Type of taper to apply to auto-correlation trace.", NUM_VALUES_FIXED );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "none", "Do not apply any taper to input trace" );
  pdef->addOption( "cos", "Apply cosine taper to input trace" );
  pdef->addOption( "hanning", "Apply 'Hanning' cosine taper to input trace. Taper length is 1/2 trace" );
  pdef->addOption( "blackman", "Apply 'Blackman' taper (alpha=0.16) to input trace" );

  pdef->addParam( "taper_len", "For FORWARD transform, taper length [ms]. THIS FEATURE IS UNTESTED!", NUM_VALUES_FIXED );
  pdef->addValue( "20", VALTYPE_NUMBER, "Taper length [ms]" );

  pdef->addParam( "taper_auto", "Apply cosine taper to auto-correlation trace", NUM_VALUES_FIXED, "Taper is applied before FFT" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Taper start time [ms] (relative to lag zero)" );
  pdef->addValue( "500", VALTYPE_NUMBER, "Taper length [ms]" );

  pdef->addParam( "norm", "Normalize output", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Normalize output values" );
  pdef->addOption( "no", "Do not normalize output values" );
  pdef->addValue( "nsamp", VALTYPE_OPTION );
  pdef->addOption( "nsamp", "Normalize by number of samples: Appropriate for stationary signal" );
  pdef->addOption( "sqrt", "Normalize by square root of number of samples: Appropriate for non-stationary noise" );

  pdef->addParam( "output", "Output options for frequency spectrum.", NUM_VALUES_FIXED );
  pdef->addValue( "amp", VALTYPE_OPTION );
  pdef->addOption( "amp", "Output amplitude spectrum" );
  pdef->addOption( "psd", "Output power spectrum (PSD - power spectral density)" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_spectrum_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_spectrum::VariableStruct* vars = reinterpret_cast<mod_spectrum::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_spectrum_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_spectrum::VariableStruct* vars = reinterpret_cast<mod_spectrum::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->bufferAuto != NULL ) {
    delete [] vars->bufferAuto;
    vars->bufferAuto = NULL;
  }
  if( vars->fftObj != NULL ) {
    delete vars->fftObj;
    vars->fftObj = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spectrum_( csParamDef* pdef ) {
  params_mod_spectrum_( pdef );
}
extern "C" void _init_mod_spectrum_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spectrum_( param, env, writer );
}
extern "C" bool _start_exec_mod_spectrum_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spectrum_( env, writer );
}
extern "C" void _exec_mod_spectrum_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spectrum_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spectrum_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spectrum_( env, writer );
}

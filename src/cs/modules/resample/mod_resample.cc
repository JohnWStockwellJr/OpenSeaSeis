/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_methods.h"
#include "csFFTTools.h"
#include "csInterpolation.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: RESAMPLE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_resample {
  struct VariableStruct {
    int   numSamplesOld;
    float sampleIntOld;  // [ms]
    float* buffer;
    bool debias;
    int filterFlag;
    int hdrID_scalar;
    cseis_geolib::csInterpolation* interpol;
    int normOption;
    float normScalar;
    int upSampleMethod;

    cseis_geolib::csFFTTools* fftTool;
  };
  static int const FILTER_NONE   = 0;
  static int const FILTER_FIR    = 11;
  static int const FILTER_BUTTER = 22;
  static int const NORM_YES   = 101;
  static int const NORM_NO    = 102;
  static int const NORM_RMS   = 103;
  static int const UPSAMPLE_SINC   = 301;
  static int const UPSAMPLE_QUAD   = 302;
  static int const UPSAMPLE_LINEAR   = 303;
}
using namespace mod_resample;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_resample_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->buffer = NULL;
  vars->debias = false;
  vars->filterFlag = mod_resample::FILTER_NONE;
  vars->hdrID_scalar = -1;
  vars->interpol = NULL;
  vars->normOption = mod_resample::NORM_YES;
  vars->normScalar = 1.0f;
  vars->upSampleMethod = mod_resample::UPSAMPLE_SINC;
 
  //---------------------------------------------
  vars->numSamplesOld = shdr->numSamples;
  vars->sampleIntOld  = shdr->sampleInt;

  float sampleIntNew;
  param->getFloat( "sample_int", &sampleIntNew );  // [ms]

  if( sampleIntNew <= 0 ) {
    writer->error("Non-physical sample interval specified: %fms", sampleIntNew);
  }

  if( param->exists( "debias" ) ) {
    std::string text;
    param->getString( "debias", &text );
    text = toLowerCase( text );
    if( !text.compare( "yes" ) ) {
      vars->debias = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->debias = false;
    }
    else {
      writer->line("Unknown argument for user parameter 'debias': '%s'.", text.c_str());
      env->addError();
    }
  }

  float cutOffFreq = 120.0f;
  float slope      = 200.0f;
  if( param->exists( "anti_alias" ) ) {
    vars->filterFlag = mod_resample::FILTER_BUTTER;
    param->getFloat("anti_alias", &cutOffFreq, 0 );
    param->getFloat("anti_alias", &slope, 1 );
  }

  bool isSinc = true;
  if( param->exists( "upsampling" ) ) {
    std::string text;
    param->getString( "upsampling", &text );
    if( !text.compare( "sinc" ) ) {
      isSinc = true;
      vars->upSampleMethod = mod_resample::UPSAMPLE_SINC;
    }
    else if( !text.compare( "quad" ) ) {
      isSinc = false;
      vars->upSampleMethod = mod_resample::UPSAMPLE_QUAD;
    }
    else if( !text.compare( "linear" ) ) {
      vars->upSampleMethod = mod_resample::UPSAMPLE_LINEAR;
    }
    else {
      writer->line("Unknown option: '%s'.", text.c_str());
      env->addError();
    }
  }
  if( param->exists( "norm" ) ) {
    std::string text;
    param->getString( "norm", &text );
    if( !text.compare( "rms" ) ) {
      vars->normOption = mod_resample::NORM_RMS;
    }
    else if( !text.compare( "yes" ) ) {
      vars->normOption = mod_resample::NORM_YES;
    }
    else if( !text.compare( "no" ) ) {
      vars->normOption = mod_resample::NORM_NO;
    }
    else {
      writer->line("Unknown option: '%s'.", text.c_str());
      env->addError();
    }
  }

  int numSamplesNew = shdr->numSamples;
  if( sampleIntNew < shdr->sampleInt ) {
    numSamplesNew = (int)( ((float)vars->numSamplesOld)*(shdr->sampleInt/sampleIntNew) );
    vars->buffer = new float[numSamplesNew];
    if( isSinc ) vars->interpol = new csInterpolation( vars->numSamplesOld, vars->sampleIntOld, 8 );
  }
  else if( sampleIntNew > shdr->sampleInt ) {
    float freqNy  = 500.0/shdr->sampleInt;
    double ratio  = (int)((double)sampleIntNew / (double)shdr->sampleInt );
    if( fabs( sampleIntNew - shdr->sampleInt*ratio) > 0.001 ) {
      writer->warning("Output sample interval was set to %f. Current implementation requires new_sample_int to be N*old_sample_int, where N=2^M", shdr->sampleInt*ratio);
    }
    sampleIntNew = shdr->sampleInt*ratio;
    numSamplesNew     = (int)ceil((double)vars->numSamplesOld / ratio);  // Workaround BUGFIX 100504: Avoids clash for num samples = 2^N+1
    vars->fftTool     = new csFFTTools( vars->numSamplesOld, numSamplesNew, vars->sampleIntOld, sampleIntNew );

    if( vars->filterFlag != mod_resample::FILTER_NONE ) {
      vars->fftTool->setFilter( cutOffFreq, slope );
    }
    writer->line("Sample int old/new: %f/%f ms\nNyquist: %f Hz\nCut-off frequency: %f Hz\nRatio: %f\nSlope[db/oct]: %f, #samples old/new/fft/fftout: %d/%d/%d/%d\n",
                 vars->sampleIntOld, sampleIntNew, freqNy, cutOffFreq, ratio, slope, vars->numSamplesOld, numSamplesNew, vars->fftTool->numFFTSamples(), vars->fftTool->numFFTSamplesOut() );

    if( edef->isDebug() ) {
      fprintf(stderr,"SampleInt: %f, Nyquist: %f, cutOff: %f, ratio: %f\n", shdr->sampleInt, freqNy, cutOffFreq, ratio );
    }
  }
  else {
    writer->warning("Specified sample interval of %fms is the same as in the input data.\nData will stay unchanged.", sampleIntNew );
  }

  shdr->numSamples = numSamplesNew;
  shdr->sampleInt  = sampleIntNew;

  if( edef->isDebug() ) {
    float ratio = shdr->sampleInt/vars->sampleIntOld;
    int ratio_int = (int)round( ratio );
    writer->line("nSamples (old,new): %d %d, sampleInt (old,new): %8.4f %8.4f, ratio: %f (%d)",
                 vars->numSamplesOld, shdr->numSamples, vars->sampleIntOld, shdr->sampleInt, ratio, ratio_int);
    fprintf(stderr,"nSamples (old,new): %d %d, sampleInt (old,new): %8.4f %8.4f, ratio: %f (%d)\n",
            vars->numSamplesOld, shdr->numSamples, vars->sampleIntOld, shdr->sampleInt, ratio, ratio_int);
  }

  if( vars->normOption == mod_resample::NORM_RMS ) {
    if( !hdef->headerExists("resample_scalar") ) {
      hdef->addHeader( cseis_geolib::TYPE_FLOAT, "resample_scalar", "Inverse scalar applied during RESAMPLE" );
    }
    vars->hdrID_scalar = hdef->headerIndex("resample_scalar");
    writer->line("Apply normalization?  Yes, using trace RMS. Store applied scalar value to trace header 'resample_scalar'\n");
  }
  else if( vars->normOption == mod_resample::NORM_YES ) {
    float ratio = shdr->sampleInt / vars->sampleIntOld;
    vars->normScalar = 1.0/ratio;
    writer->line("Apply normalization? Yes, using ratio between sample intervals Scalar = %.4f / %.4f = %f\n", shdr->sampleInt, vars->sampleIntOld, vars->normScalar);
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_resample_(
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
  int numSamplesNew = shdr->numSamples;

  float resampleScalar = 1.0f;
  // Case #1: Interpolate to smaller sample interval
  if( shdr->sampleInt < vars->sampleIntOld ) {
    if( vars->upSampleMethod == mod_resample::UPSAMPLE_SINC ) {
      for( int isamp = 0; isamp < numSamplesNew; isamp++ ) {
        float time_ms = (float)isamp * shdr->sampleInt;
        vars->buffer[isamp] = vars->interpol->valueAt( time_ms, samples );
      }
    }
    else if( vars->upSampleMethod == mod_resample::UPSAMPLE_QUAD ) {
      for( int isamp = 0; isamp < numSamplesNew; isamp++ ) {
        double sampleIndexOld = (double)(isamp*shdr->sampleInt) / vars->sampleIntOld;
        vars->buffer[isamp] = getQuadAmplitudeAtSample( samples, sampleIndexOld, vars->numSamplesOld );
      }
    }
    else { // if( vars->upSampleMethod == mod_resample::UPSAMPLE_LINEAR ) {
      for( int isamp = 0; isamp < numSamplesNew; isamp++ ) {
        double sampleIndexOld = (double)(isamp*shdr->sampleInt) / vars->sampleIntOld;
        vars->buffer[isamp] = getLinAmplitudeAtSample( samples, sampleIndexOld, vars->numSamplesOld );
      }
    }

    memcpy( samples, vars->buffer, numSamplesNew*sizeof(float) );
  }
  // Case #2: Resample to larger sample interval. Apply anti-alias filter if specified
  else if( shdr->sampleInt > vars->sampleIntOld ) {
    float mean = 0.0f;
    if( vars->debias ) { // Remove DC bias first
      double sum = 0.0;
      int sampleCounter = 0;
      for( int isamp = 0; isamp < vars->numSamplesOld; isamp++ ) {
        sampleCounter += 1;
        sum += samples[isamp];
      }
      mean = sum / (double)sampleCounter;
      for( int isamp = 0; isamp < vars->numSamplesOld; isamp++ ) {
        samples[isamp] -= mean;
      }
    }
    resampleScalar = vars->fftTool->resample( samples, vars->filterFlag != mod_resample::FILTER_NONE, vars->normOption == mod_resample::NORM_RMS );
    if( vars->normOption == mod_resample::NORM_YES ) {
      for( int isamp = 0; isamp < vars->numSamplesOld; isamp++ ) {
        samples[isamp] *= vars->normScalar;
      }      
    }
    if( vars->debias ) { // Reapply DC bias afterwards
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] += mean;
      }
    }
  }
  if( vars->normOption == mod_resample::NORM_RMS ) {
    if( resampleScalar > 0 ) resampleScalar = 1/resampleScalar;
    trace->getTraceHeader()->setFloatValue( vars->hdrID_scalar, resampleScalar );
  }
  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_resample_( csParamDef* pdef ) {
  pdef->setModule( "RESAMPLE", "Resample trace to different sample interval" );

  pdef->addParam( "sample_int", "New sample interval [ms]", NUM_VALUES_FIXED, "Must be even multiplier of current sample interval" );
  pdef->addValue( "", VALTYPE_NUMBER, "Sample interval [ms]" );

  pdef->addParam( "debias", "Remove DC bias before resampling?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Remove DC bias before resampling, reapply afterwards.", "This should be done to avoid FFT related artefacts due to DC bias" );
  pdef->addOption( "no", "Do not remove DC bias before resampling." );

  pdef->addParam( "anti_alias", "Anti-alias filter", NUM_VALUES_FIXED, "Do not specify if no anti-alias filter is required" );
  pdef->addValue( "", VALTYPE_NUMBER, "Cut-off frequency [Hz] (-3dB point)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Filter slope [dB/oct]" );

  pdef->addParam( "upsampling", "Interpolation method used for upsampling (output sample interval < input sample interval)", NUM_VALUES_FIXED );
  pdef->addValue( "sinc", VALTYPE_OPTION );
  pdef->addOption( "sinc", "Use sinc interpolation" );
  pdef->addOption( "quad", "Use quadratic interpolation" );
  pdef->addOption( "linear", "Linear interpolation" );
 
  pdef->addParam( "norm", "Apply normalization?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Apply normalization: Ratio of input/output sample intervals" );
  pdef->addOption( "no", "Do not apply normalization" );
  pdef->addOption( "rms", "Apply trace-by-trace RMS normalization. Inverse scalar is stored in output trace header 'resample_scalar'. Applying this scalar reverses the normalization" );
}



//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_resample_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_resample::VariableStruct* vars = reinterpret_cast<mod_resample::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_resample_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_resample::VariableStruct* vars = reinterpret_cast<mod_resample::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->buffer != NULL ) {
    delete [] vars->buffer;
    vars->buffer = NULL;
  }
  if( vars->fftTool != NULL ) {
    delete vars->fftTool;
    vars->fftTool = NULL;
  }
  if( vars->interpol != NULL ) {
    delete vars->interpol;
    vars->interpol = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_resample_( csParamDef* pdef ) {
  params_mod_resample_( pdef );
}
extern "C" void _init_mod_resample_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_resample_( param, env, writer );
}
extern "C" bool _start_exec_mod_resample_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_resample_( env, writer );
}
extern "C" void _exec_mod_resample_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_resample_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_resample_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_resample_( env, writer );
}

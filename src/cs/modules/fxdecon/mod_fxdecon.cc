/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstring>
#include "cseis_includes.h"
#include "csFXDecon.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: FXDECON
 *
 */
namespace mod_fxdecon {
  struct VariableStruct {
    csFXDecon* fxdecon;
    int output;
    Attr attr;
  };
  static int const MODE_ENSEMBLE = 11;
  static int const MODE_TRACE    = 12;
  static int const MODE_ALL      = 13;

  static int const OUTPUT_FILT = 41;
  static int const OUTPUT_DIFF = 42;
}
using namespace mod_fxdecon;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_fxdecon_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  //  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->fxdecon = NULL;
  vars->output = OUTPUT_FILT;

  env->execPhaseDef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->attr.fmin           = 0;
  vars->attr.fmax           = 0;
  vars->attr.winLen_samp    = 0;
  vars->attr.taperLen_samp  = 0;
  vars->attr.taperLen_s     = 0;
  vars->attr.numWin         = 0;
  vars->attr.ntraces_design = 0;
  vars->attr.ntraces_filter = 0;

  if( param->exists( "freq_range" ) ) {
    param->getFloat( "freq_range", &vars->attr.fmin, 0 );
    param->getFloat( "freq_range", &vars->attr.fmax, 1 );
  }

  float windowLength_ms = shdr->numSamples * shdr->sampleInt;
  if( param->exists( "win_len" ) ) {
    param->getFloat( "win_len", &windowLength_ms, 0 );
    vars->attr.winLen_samp = (int)round( windowLength_ms / shdr->sampleInt ) + 1;
  }
  vars->attr.numWin = (int)round( (shdr->numSamples-1) * shdr->sampleInt / windowLength_ms );

  vars->attr.taperLen_samp  = (int)round( 0.1 * vars->attr.winLen_samp );
  if( param->exists( "taper_len" ) ) {
    float taperLength_ms;
    param->getFloat( "taper_len", &taperLength_ms, 0 );
    vars->attr.taperLen_samp = (int)round( taperLength_ms / shdr->sampleInt );
    vars->attr.taperLen_samp = ( vars->attr.taperLen_samp%2 ? vars->attr.taperLen_samp+1 : vars->attr.taperLen_samp ); 
  }

  if( param->exists( "win_traces" ) ) {
    param->getInt( "win_traces", &vars->attr.ntraces_design, 0 );
    param->getInt( "win_traces", &vars->attr.ntraces_filter, 1 );
  }

  if( param->exists("output") ) {
    std::string text;
    param->getString("output", &text );
    if( !text.compare("filt") ) {
      vars->output = mod_fxdecon::OUTPUT_FILT;
    }
    else if( !text.compare("diff") ) {
      vars->output = mod_fxdecon::OUTPUT_DIFF;
    }
  }

  vars->fxdecon = new mod_fxdecon::csFXDecon();
  vars->attr.taperLen_s = vars->attr.taperLen_samp * shdr->sampleInt / 1000.0;
  if( vars->attr.numWin == 0 ) vars->attr.taperLen_s = 0;

  vars->fxdecon->initialize( shdr->sampleInt, shdr->numSamples, vars->attr );

  if( edef->isDebug() ) {
    vars->fxdecon->dump();
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_fxdecon_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;


  int numSamples = traceGather->trace(0)->numSamples();
  int numTraces  = traceGather->numTraces();

  if( numTraces < vars->attr.ntraces_design ) return; // Do not apply FX decon if input ensemble has too few traces

  float** samplesIn  = new float*[numTraces];
  float** samplesOut = new float*[numTraces];
  for( int itrc = 0; itrc < numTraces; itrc++ ) {
    samplesIn[itrc]  = traceGather->trace(itrc)->getTraceSamples();
    samplesOut[itrc] = new float[numSamples];
    memset( samplesOut[itrc], 0, numSamples*sizeof(float) );
  }

  vars->fxdecon->apply( samplesIn, samplesOut, numTraces, numSamples );
  if( vars->output == mod_fxdecon::OUTPUT_FILT ) {
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      memcpy( samplesIn[itrc], samplesOut[itrc], numSamples*sizeof(float) );
    }
  }
  else if( vars->output == mod_fxdecon::OUTPUT_DIFF ) {
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samplesIn[itrc][isamp] -= samplesOut[itrc][isamp];
      }
    }
  }

  for( int itrc = 0; itrc < numTraces; itrc++ ) {
    delete [] samplesOut[itrc];
  }
  delete [] samplesIn;
  delete [] samplesOut;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_fxdecon_( csParamDef* pdef ) {
  pdef->setModule( "FXDECON", "FX deconvolution");

  pdef->addParam( "freq_range", "Frequency range to filter", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Minimum frequency [Hz]" );
  pdef->addValue( "12", VALTYPE_NUMBER, "Maximum frequency [Hz]" );

  pdef->addParam( "win_len", "Window length [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "2000", VALTYPE_NUMBER );

  pdef->addParam( "win_traces", "Spatial/filter window", NUM_VALUES_FIXED );
  pdef->addValue( "10", VALTYPE_NUMBER, "Number of traces in design window" );
  pdef->addValue( "4", VALTYPE_NUMBER, "Number of traces for filter (smaller than window size)" );

  pdef->addParam( "taper_len", "Taper length [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "100", VALTYPE_NUMBER );

  pdef->addParam( "output", "Output data", NUM_VALUES_FIXED );
  pdef->addValue( "filt", VALTYPE_OPTION );
  pdef->addOption( "filt", "Output filtered data" );
  pdef->addOption( "diff", "Output difference between input and filtered data" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_fxdecon_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_fxdecon::VariableStruct* vars = reinterpret_cast<mod_fxdecon::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_fxdecon_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_fxdecon::VariableStruct* vars = reinterpret_cast<mod_fxdecon::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->fxdecon != NULL ) {
    delete vars->fxdecon;
    vars->fxdecon = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_fxdecon_( csParamDef* pdef ) {
  params_mod_fxdecon_( pdef );
}
extern "C" void _init_mod_fxdecon_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_fxdecon_( param, env, writer );
}
extern "C" bool _start_exec_mod_fxdecon_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_fxdecon_( env, writer );
}
extern "C" void _exec_mod_fxdecon_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_fxdecon_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_fxdecon_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_fxdecon_( env, writer );
}

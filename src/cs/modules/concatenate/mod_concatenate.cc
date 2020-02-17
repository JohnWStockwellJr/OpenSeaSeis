/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: CONCATENATE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_concatenate {
  struct VariableStruct {
    int mode;
    int numTraces;
    int numSamplesOrig;
  };
  static const int MODE_TRACES   = 10;
  static const int MODE_ENSEMBLE = 11;
}
using namespace mod_concatenate;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_concatenate_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );


  vars->mode           = MODE_TRACES;
  vars->numTraces      = 0;
  vars->numSamplesOrig = shdr->numSamples;

//---------------------------------------------
// Retrieve sampIndex and scalar
//
  csVector<std::string> valueList;

  if( param->exists("mode") ) {
    string text;
    param->getString( "mode", &text );
    if( !text.compare("traces") ) {
      vars->mode = MODE_TRACES;
    }
    else if( !text.compare("ensemble") ) {
      vars->mode = MODE_ENSEMBLE;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  param->getInt( "ntraces", &vars->numTraces );

  if( vars->mode == MODE_TRACES ) {
    edef->setTraceSelectionMode( TRCMODE_FIXED, vars->numTraces );
    //edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }
  else {
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }

  shdr->numSamples = vars->numTraces*shdr->numSamples;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_concatenate_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;


  int nTraces = traceGather->numTraces();
//  fprintf(stderr,"Concatenate: Input traces: %d, concatenate %d traces\n", nTraces, vars->numTraces );
  if( nTraces == 0 ) return;

  int numTracesToPass = 1;

  if( vars->mode == MODE_TRACES ) {
    float* samplesOutPtr = traceGather->trace(0)->getTraceSamples();
    // Concatenate all traces into trace 1
    for( int itrc = 1; itrc < vars->numTraces; itrc++ ) {
      float* tracePtr = traceGather->trace(itrc)->getTraceSamples();
      memcpy( &samplesOutPtr[itrc*vars->numSamplesOrig], tracePtr, vars->numSamplesOrig*sizeof(float) );
    }
    for( int itrc = nTraces; itrc < vars->numTraces; itrc++ ) {
      float* samples = &samplesOutPtr[itrc*vars->numSamplesOrig];
      for( int isamp = 0; isamp < vars->numSamplesOrig; isamp++ ) {
        samples[isamp] = 0.0;
      }
    }
  }
  else {  // MODE_ENSEMBLE
    if( nTraces <= vars->numTraces ) {
      // Actual number of ensemble traces is smaller than specified number of traces
      // --> Concatenate all traces, fill rest with zeros
      float* samplesOutPtr = traceGather->trace(0)->getTraceSamples();
      for( int itrc = 1; itrc < nTraces; itrc++ ) {
        float* samplesInPtr = traceGather->trace(itrc)->getTraceSamples();
        memcpy( &samplesOutPtr[itrc*vars->numSamplesOrig], samplesInPtr, vars->numSamplesOrig*sizeof(float) );
      }
      for( int itrc = nTraces; itrc < vars->numTraces; itrc++ ) {
        float* samples = &samplesOutPtr[itrc*vars->numSamplesOrig];
        for( int isamp = 0; isamp < vars->numSamplesOrig; isamp++ ) {
          samples[isamp] = 0.0;
        }
      }
    }
    else {
      int nLoops = (int)( (nTraces-1) / vars->numTraces ) + 1;
      for( int iloop = 0; iloop < nLoops; iloop++ ) {
        float* samplesOutPtr = traceGather->trace(iloop)->getTraceSamples();

        int reduce = iloop*vars->numTraces;
        int traceIndexStart = iloop*vars->numTraces;
        int traceIndexEnd = traceIndexStart + vars->numTraces;
        if( traceIndexStart == 0 ) traceIndexStart = 1;
        if( traceIndexEnd > nTraces ) traceIndexEnd = nTraces;

        for( int itrc = traceIndexStart; itrc < traceIndexEnd; itrc++ ) {
          float* samplesInPtr = traceGather->trace(itrc)->getTraceSamples();
          memcpy( &samplesOutPtr[(itrc-reduce)*vars->numSamplesOrig], samplesInPtr, vars->numSamplesOrig*sizeof(float) );
        }

        if( iloop == (nLoops-1) ) {
          for( int itrc = nTraces; itrc < nLoops*vars->numTraces; itrc++ ) {
            float* samples = &samplesOutPtr[(itrc-reduce)*vars->numSamplesOrig];
            for( int isamp = 0; isamp < vars->numSamplesOrig; isamp++ ) {
              samples[isamp] = 0.0;
            }
          }
        }
      }

    }
    
  }
  
  // Free all remaining traces
  traceGather->freeTraces( numTracesToPass, nTraces-numTracesToPass );
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_concatenate_( csParamDef* pdef ) {
  pdef->setModule( "concatenate", "Concatenate adjacent traces" );

  pdef->addDoc("This module concatenates traces to form a longer output trace.");
  pdef->addDoc("For example, if user parameter 'ntraces' is set to 3, three sequential traces are concatenated from top to bottom, forming a new output trace with 3x the trace length");

  pdef->addParam( "mode", "Mode of concatenation", NUM_VALUES_FIXED );
  pdef->addValue( "traces", VALTYPE_OPTION );
  pdef->addOption( "traces", "Concatenate fixed number of sequential traces, specified in user parameter 'ntraces'" );
  pdef->addOption( "ensemble", "Concatenate all traces in input ensemble. Specify the maximum number of traces expected in input ensembles in user parameter 'ntraces'." );

  pdef->addParam( "ntraces", "Number of adjacent traces to concatenate", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Number of traces to concatenate" );

  //  pdef->addParam( "absolute_time", "Take into account absolute time?", NUM_VALUES_VARIABLE );
  // pdef->addValue( "no", VALTYPE_OPTION );
  // pdef->addOption( "no", "Concatenate traces last sample trace N / first sample trace N+1" );
  // pdef->addOption( "yes", "Concatenate traces according to absolute time found in standard trace headers time_samp1 & time_samp1_us" );
  // pdef->addValue( "first", VALTYPE_OPTION );
  // pdef->addOption( "first", "In case of overlapping times, use sample values from first trace" );
  // pdef->addOption( "stack", "Stack overlapping parts" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_concatenate_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_concatenate::VariableStruct* vars = reinterpret_cast<mod_concatenate::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_concatenate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_concatenate::VariableStruct* vars = reinterpret_cast<mod_concatenate::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_concatenate_( csParamDef* pdef ) {
  params_mod_concatenate_( pdef );
}
extern "C" void _init_mod_concatenate_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_concatenate_( param, env, writer );
}
extern "C" bool _start_exec_mod_concatenate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_concatenate_( env, writer );
}
extern "C" void _exec_mod_concatenate_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_concatenate_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_concatenate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_concatenate_( env, writer );
}

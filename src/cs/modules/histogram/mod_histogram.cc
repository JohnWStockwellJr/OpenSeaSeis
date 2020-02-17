/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csStandardHeaders.h"
#include "csVector.h"
#include "csStackUtil.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: HISTOGRAM
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_histogram {
  struct VariableStruct {
    float ampMin;
    float ampMax;
    float ampStep;
    int numSamplesIn;
    float sampleIntIn;
    int numSamplesOut;
    float sampleIntOut;
    int* histBuffer;
  }; 
}
using namespace mod_histogram;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_histogram_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  //  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );


  vars->numSamplesIn  = shdr->numSamples;
  vars->numSamplesOut = shdr->numSamples;
  vars->sampleIntIn   = shdr->sampleInt;
  vars->sampleIntOut  = shdr->sampleInt;
  vars->ampMin = 0.0;
  vars->ampMax = 0.0;
  vars->ampStep = 1.0;
  vars->histBuffer = NULL;

  //  vars->rms = 0.0;
  //  vars->counter = 0;
  //  vars->traceCounter = 0;
  //  vars->rms_prev = 0.0;

  if( param->exists("hist") ) {
    param->getFloat("hist", &vars->ampMin, 0 );
    param->getFloat("hist", &vars->ampMax, 1 );
    param->getFloat("hist", &vars->ampStep, 2 );
    vars->numSamplesOut = (int)((vars->ampMax - vars->ampMin)/vars->ampStep + 0.5);
    vars->sampleIntOut  = (vars->ampMax - vars->ampMin) / (float)vars->numSamplesOut;
    vars->histBuffer = new int[vars->numSamplesOut];
    writer->line("Amplitudes: %f - %f (step %f)", vars->ampMin, vars->ampMax, vars->sampleIntOut);
  }
  else {
    //    vars->hdrID_rms = hdef->addHeader( TYPE_DOUBLE, "rms", "rms value" );
    //    vars->hdrID_counter = hdef->addHeader( TYPE_INT, "counter", "Counter" );
  }
  if( param->exists("ntraces") ) {
    int ntraces = 0;
    param->getInt("ntraces", &ntraces );
    edef->setTraceSelectionMode( TRCMODE_FIXED, ntraces );
  }
  else {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }
  shdr->sampleInt  = vars->sampleIntOut;
  shdr->numSamples = vars->numSamplesOut;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_histogram_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  //  csSuperHeader const*    shdr = env->superHeader;


  //--------------------------------------------------------------------------------
  for( int isamp = 0; isamp < vars->numSamplesOut; isamp++ ) {
    vars->histBuffer[isamp] = 0;
  }
  int ntraces = traceGather->numTraces();
  for( int itrc = 0; itrc < ntraces; itrc++ ) {
    float* samples = traceGather->trace(itrc)->getTraceSamples();
    for( int isamp = 0; isamp < vars->numSamplesIn; isamp++ ) {
      int index = (int)(( samples[isamp] - vars->ampMin ) / vars->sampleIntOut + 0.5 );
      if( index < 0 ) index = 0;
      else if( index >= vars->numSamplesOut ) index = vars->numSamplesOut-1;
      vars->histBuffer[index] += 1;
    }
  }
  if( ntraces > 1 ) traceGather->freeTraces( 1, ntraces-1 );
  float* samples = traceGather->trace(0)->getTraceSamples();
  for( int isamp = 0; isamp < vars->numSamplesOut; isamp++ ) {
    samples[isamp] = (float)vars->histBuffer[isamp];
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_histogram_( csParamDef* pdef ) {
  pdef->setModule( "HISTOGRAM", "Histogram" );

  pdef->addParam( "ntraces", "Number of input traces", NUM_VALUES_FIXED);
  pdef->addValue( "0", VALTYPE_NUMBER );

  pdef->addParam( "hist", "Amplitude range for histogram", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_NUMBER, "Minimum amplitude" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum amplitude" );
  pdef->addValue( "", VALTYPE_NUMBER, "Step" );
}



//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_histogram_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_histogram::VariableStruct* vars = reinterpret_cast<mod_histogram::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_histogram_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_histogram::VariableStruct* vars = reinterpret_cast<mod_histogram::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->histBuffer != NULL ) {
    delete [] vars->histBuffer;
    vars->histBuffer = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_histogram_( csParamDef* pdef ) {
  params_mod_histogram_( pdef );
}
extern "C" void _init_mod_histogram_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_histogram_( param, env, writer );
}
extern "C" bool _start_exec_mod_histogram_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_histogram_( env, writer );
}
extern "C" void _exec_mod_histogram_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_histogram_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_histogram_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_histogram_( env, writer );
}

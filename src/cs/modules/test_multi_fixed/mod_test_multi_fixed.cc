/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cstdio>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_test_multi_fixed {
  struct VariableStruct {
    int callCounter;
    int numTracesToRoll;
    int numTracesToAdd;
    int numTraces;
    int repeat;
    int skip;
    int hdrId_bias;
  };
}
using mod_test_multi_fixed::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_test_multi_fixed_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  //  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->callCounter     = 0;
  vars->numTracesToRoll = 0;
  vars->numTracesToAdd  = 0;
  vars->numTraces       = 1;
  vars->repeat          = 0;
  vars->skip            = 0;
  vars->hdrId_bias      = -1;


  if( param->exists("repeat") ) {
    param->getInt("repeat", &vars->repeat );
    if( vars->repeat < 0 ) {
      writer->error("Incorrect entry for number of times to repeat proessing: %d", vars->repeat);
    }
  }
  if( param->exists("ntraces_roll") ) {
    param->getInt("ntraces_roll", &vars->numTracesToRoll );
    if( vars->numTracesToRoll < 0 ) {
      writer->error("Incorrect entry for number of traces to roll: %d", vars->numTracesToRoll);
    }
  }
  if( param->exists("skip") ) {
    param->getInt("skip", &vars->skip );
    if( vars->skip < 0 ) {
      writer->error("Incorrect entry for number of times to skip: %d", vars->skip);
    }
  }
  if( param->exists("ntraces_in") ) {
    param->getInt("ntraces_in", &vars->numTraces );
    if( vars->numTraces <= 0 ) {
      writer->error("Incorrect entry for number of input traces: %d", vars->numTraces);
    }
  }
  if( param->exists("ntraces_add") ) {
    param->getInt("ntraces_add", &vars->numTracesToAdd );
    if( vars->numTracesToAdd <= 0 ) {
      writer->error("Incorrect entry for number of traces to add: %d", vars->numTracesToAdd);
    }
  }

  if( vars->numTracesToAdd != 0 && vars->numTracesToRoll != 0 ) {
    writer->error("Cannot add and roll traces at the same time... Specify one at a time");
  }

  if( !hdef->headerExists("dc") ) {
    hdef->addStandardHeader( "dc" );
  }
  vars->hdrId_bias = hdef->headerIndex("dc");

  //  int numSamples = shdr->numSamples;
  //  float sampleInt = shdr->sampleInt;

  edef->setTraceSelectionMode( TRCMODE_FIXED, vars->numTraces );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_test_multi_fixed_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;


  int numTracesIn = traceGather->numTraces();
  int numTracesOut;
  if( numTracesIn == vars->numTraces && !edef->isLastCall() ) {
    numTracesOut = numTracesIn - vars->numTracesToRoll;
    *numTrcToKeep = vars->numTracesToRoll;
  }
  else {  // Number of input traces is less than specified and/or last call: This is the last call to this module
    numTracesOut = numTracesIn;
    *numTrcToKeep = 0;
  }

  if( edef->isDebug() ) {
    writer->line("Number of input traces: %d, last call: %d", numTracesIn, edef->isLastCall());
    fprintf(stdout,"Number of input traces: %d, last call: %d\n", numTracesIn, edef->isLastCall());
  }

  vars->callCounter += 1;
  if( vars->callCounter > vars->skip ) {
    double mean = 0.0;
    for( int itrc = 0; itrc < numTracesIn; itrc++ ) {
      float* samples = traceGather->trace(itrc)->getTraceSamples();
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        mean += samples[isamp];
      }
    }
    mean /= numTracesIn*shdr->numSamples;
    if( edef->isDebug() ) writer->line("Mean amplitude: %f", mean);
    for( int itrc = 0; itrc < numTracesOut; itrc++ ) {
      traceGather->trace(itrc)->getTraceHeader()->setFloatValue( vars->hdrId_bias, mean );
      float* samples = traceGather->trace(itrc)->getTraceSamples();
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] -= (float)mean;
      }
    }
    if( vars->numTracesToAdd > 0 ) {
      traceGather->createTraces( numTracesIn, vars->numTracesToAdd, hdef, shdr->numSamples );
      csTraceData* trcDataOrig  = traceGather->trace(0)->getTraceDataObject();
      csTraceHeader* trcHdrOrig = traceGather->trace(0)->getTraceHeader();
      for( int itrc = numTracesIn; itrc < traceGather->numTraces(); itrc++ ) {
        traceGather->trace(itrc)->getTraceDataObject()->setData( trcDataOrig );
        traceGather->trace(itrc)->getTraceHeader()->copyFrom( trcHdrOrig );
        traceGather->trace(itrc)->getTraceHeader()->setFloatValue( vars->hdrId_bias, mean );
      }
    }
    vars->callCounter = 0;
  }
  else {
    traceGather->freeAllTraces();
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_test_multi_fixed_( csParamDef* pdef ) {
  pdef->setModule( "TEST_MULTI_FIXED", "Test module - multi-trace, fixed number of input traces" );

  pdef->addParam( "ntraces_in", "Number of input traces", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of input traces" );

  pdef->addParam( "skip", "Number of times to skip processing and drop all input traces", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of times to skip processing and drop all input traces" );

  pdef->addParam( "ntraces_roll", "Number of traces to roll", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to roll" );

  pdef->addParam( "ntraces_add", "Number of traces to add", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to add each time module is called" );

  pdef->addParam( "repeat", "Number of times to repeat", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of times to repeat each processing" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_test_multi_fixed_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_test_multi_fixed::VariableStruct* vars = reinterpret_cast<mod_test_multi_fixed::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_test_multi_fixed_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_test_multi_fixed::VariableStruct* vars = reinterpret_cast<mod_test_multi_fixed::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_test_multi_fixed_( csParamDef* pdef ) {
  params_mod_test_multi_fixed_( pdef );
}
extern "C" void _init_mod_test_multi_fixed_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_test_multi_fixed_( param, env, writer );
}
extern "C" bool _start_exec_mod_test_multi_fixed_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_test_multi_fixed_( env, writer );
}
extern "C" void _exec_mod_test_multi_fixed_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_test_multi_fixed_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_test_multi_fixed_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_test_multi_fixed_( env, writer );
}

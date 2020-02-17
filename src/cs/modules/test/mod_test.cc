/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_test {
  struct VariableStruct {
    int traceCounter;
    int numTracesToSkip;
    int hdrID;
  };
}
using mod_test::VariableStruct;

//*************************************************************************************************
// Init phase
//
void init_mod_test_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  //  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  vars->traceCounter    = 0;
  vars->numTracesToSkip = 0;
  vars->hdrID = -1;

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  if( param->exists("skip") ) {
    param->getInt("skip", &vars->numTracesToSkip );
    if( vars->numTracesToSkip < 0 ) {
      writer->error("Incorrect entry for number of traces to skip: %d", vars->numTracesToSkip);
    }
  }
  if( param->exists("header") ) {
    std::string text;
    param->getString("header", &text );
    vars->hdrID = hdef->headerIndex( text );
  }
  else {
    vars->hdrID = hdef->addHeader( TYPE_FLOAT, "mean", "Mean value" );
  }
  //  int numSamples  = shdr->numSamples;
  //  float sampleInt = shdr->sampleInt;
}

//*************************************************************************************************
// Exec phase
//
void exec_mod_test_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>(env->execPhaseDef->variables());
  //  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  vars->traceCounter += 1;
  if( vars->traceCounter > vars->numTracesToSkip ) {
    float* samples = trace->getTraceSamples();
    double mean = 0.0;
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      mean += samples[isamp];
    }
    mean /= shdr->numSamples;
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] -= (float)mean;
    }
    csTraceHeader* trcHdr = trace->getTraceHeader();
    trcHdr->setFloatValue( vars->hdrID, mean );
    vars->traceCounter = 0;
    return;
  }
  else {
    traceGather->freeAllTraces();
    return;
  }
}

//*************************************************************************************************
// Parameter definition
//
void params_mod_test_( csParamDef* pdef ) {
  pdef->setModule( "TEST", "Demonstration single trace module", "Extra description..." );

  pdef->addParam( "skip", "Number of traces to skip", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to skip" );

  pdef->addParam( "header", "Trace header to store mean trace value", NUM_VALUES_FIXED );
  pdef->addValue( "mean", VALTYPE_STRING, "Trace header name" );
}

//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_test_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_test::VariableStruct* vars = reinterpret_cast<mod_test::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_test_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_test::VariableStruct* vars = reinterpret_cast<mod_test::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_test_( csParamDef* pdef ) {
  params_mod_test_( pdef );
}
extern "C" void _init_mod_test_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_test_( param, env, writer );
}
extern "C" bool _start_exec_mod_test_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_test_( env, writer );
}
extern "C" void _exec_mod_test_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_test_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_test_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_test_( env, writer );
}

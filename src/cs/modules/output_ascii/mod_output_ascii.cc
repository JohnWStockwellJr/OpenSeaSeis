/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: OUTPUT_ASCII
 *
 */
namespace mod_output_ascii {
  struct VariableStruct {
    std::string filename;
    FILE* fout;
    int nTracesOut;
    bool isFirstCall;
  };
}
using namespace mod_output_ascii;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_output_ascii_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  //  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->filename    = "";
  vars->isFirstCall = true;
  vars->nTracesOut  = 0;
  vars->fout = NULL;

  param->getString( "filename", &vars->filename );

  writer->line("");
  writer->line("  File name:             %s", vars->filename.c_str());
  writer->line("  Sample interval [ms]:  %f", shdr->sampleInt);
  writer->line("  Number of samples:     %d", shdr->numSamples);
  writer->line("");
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_output_ascii_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep, 
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;


  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    vars->fout = fopen( vars->filename.c_str(), "w" );
    if( !vars->fout ) writer->error("Error occurred when opening ASCII output file: %s\nDoes folder exist?\n", vars->filename.c_str() );
  }

  int numTraces = traceGather->numTraces();
  for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      float* samples = traceGather->trace(itrc)->getTraceSamples();
      fprintf(vars->fout,"%.7e ", samples[isamp]);
    }
    fprintf(vars->fout,"\n");
  }
  vars->nTracesOut += numTraces;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_output_ascii_( csParamDef* pdef ) {
  pdef->setModule( "OUTPUT_ASCII", "Output data to ASCII file", "Writes input ensemble sample values to ASCII table in matrix form" );

  pdef->addParam( "filename", "Output file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_output_ascii_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_output_ascii::VariableStruct* vars = reinterpret_cast<mod_output_ascii::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_output_ascii_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_output_ascii::VariableStruct* vars = reinterpret_cast<mod_output_ascii::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->fout != NULL ) {
    fclose(vars->fout);
    vars->fout = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_output_ascii_( csParamDef* pdef ) {
  params_mod_output_ascii_( pdef );
}
extern "C" void _init_mod_output_ascii_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_output_ascii_( param, env, writer );
}
extern "C" bool _start_exec_mod_output_ascii_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_output_ascii_( env, writer );
}
extern "C" void _exec_mod_output_ascii_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_output_ascii_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_output_ascii_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_output_ascii_( env, writer );
}

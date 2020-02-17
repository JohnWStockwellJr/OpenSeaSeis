/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cstddef>
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: FFT_2D
 *
 * @date   2011
 */

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_fft_2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  writer->error("This is a dummy module which does not function");
}

//*******************************************************************************************
// Exec phase
//*******************************************************************************************

void exec_mod_fft_2d_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{

}

//*************************************************************************************************
// Parameter definition
//*************************************************************************************************
void params_mod_fft_2d_( csParamDef* pdef ) {
  pdef->setModule( "FFT_2D", "2D FFT - DUMMY MODULE");
  pdef->setVersion(0,5);
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_fft_2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_fft_2d::VariableStruct* vars = reinterpret_cast<mod_fft_2d::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_fft_2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_fft_2d::VariableStruct* vars = reinterpret_cast<mod_fft_2d::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  delete vars; vars = NULL;
}

extern "C" void _params_mod_fft_2d_( csParamDef* pdef ) {
  params_mod_fft_2d_( pdef );
}
extern "C" void _init_mod_fft_2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_fft_2d_( param, env, writer );
}
extern "C" bool _start_exec_mod_fft_2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_fft_2d_( env, writer );
}
extern "C" void _exec_mod_fft_2d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_fft_2d_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_fft_2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_fft_2d_( env, writer );
}

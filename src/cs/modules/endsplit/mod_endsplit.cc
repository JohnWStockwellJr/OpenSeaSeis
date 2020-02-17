/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ENDSPLIT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */

using namespace cseis_system;

void init_mod_endsplit_( csParamManager* userParams, csInitPhaseEnv* env, csLogWriter* writer ) {
  env->execPhaseDef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
}
void exec_mod_endsplit_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  traceGather->freeAllTraces();
}
void params_mod_endsplit_( csParamDef* pdef ) {
  pdef->setModule( "ENDSPLIT", "Endsplit statement", "Marks the end of a Split-endsplit block" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_endsplit_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_endsplit_( csExecPhaseEnv* env, csLogWriter* writer ) {
  // Nothing to do
}

extern "C" void _params_mod_endsplit_( csParamDef* pdef ) {
  params_mod_endsplit_( pdef );
}
extern "C" void _init_mod_endsplit_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_endsplit_( param, env, writer );
}
extern "C" bool _start_exec_mod_endsplit_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_endsplit_( env, writer );
}
extern "C" void _exec_mod_endsplit_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_endsplit_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_endsplit_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_endsplit_( env, writer );
}

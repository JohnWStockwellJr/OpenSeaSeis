/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"
#include "csSelectionManager.h"

using namespace cseis_geolib;
using namespace cseis_system;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ELSEIF
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_elseif {
  struct VariableStruct {
    csSelectionManager* selectionManager;
  };
}
using mod_elseif::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//*************************************************************************************************
void init_mod_elseif_( csParamManager* userParams, csInitPhaseEnv* env, csLogWriter* writer ) {
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  //  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
    
  vars->selectionManager = NULL;

  csVector<std::string> valueList;

  userParams->getAll( "header", &valueList );
  if( valueList.size() == 0 ) {
    writer->warning("%s: Wrong number of parameters for option 'HEADER'. Expected: > 0, found: %d.", edef->moduleName().c_str(), valueList.size());
    env->addError();
  }

  std::string text;
  userParams->getString( "select", &text );

  try {
    vars->selectionManager = new csSelectionManager();
    vars->selectionManager->set( &valueList, &text, hdef );
  }
  catch( csException& e ) {
    vars->selectionManager = NULL;
    writer->error( "%s: %s", edef->moduleName().c_str(), e.getMessage() );
  }
  if( edef->isDebug() ) vars->selectionManager->dump();
}

//*************************************************************************************************
// Exec phase
//
//
//*************************************************************************************************
void exec_mod_elseif_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;

  csTrace* trace = traceGather->trace(0);


  if( vars->selectionManager->contains( trace->getTraceHeader() ) ) {
    *port = 0;
  }
  else {
    *port = 1;
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_elseif_( csParamDef* pdef ) {
  pdef->setModule( "ELSEIF", "Elseif statement", "Branch traces that match specified header selection" );

  pdef->addParam( "header", "Names of trace headers used for trace selection", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of trace header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING,
      "List of selection strings, one for each specified header. See documentation for more detailed description of selection syntax" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_elseif_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_elseif::VariableStruct* vars = reinterpret_cast<mod_elseif::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_elseif_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_elseif::VariableStruct* vars = reinterpret_cast<mod_elseif::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->selectionManager ) {
    delete vars->selectionManager;
    vars->selectionManager = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_elseif_( csParamDef* pdef ) {
  params_mod_elseif_( pdef );
}
extern "C" void _init_mod_elseif_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_elseif_( param, env, writer );
}
extern "C" bool _start_exec_mod_elseif_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_elseif_( env, writer );
}
extern "C" void _exec_mod_elseif_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_elseif_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_elseif_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_elseif_( env, writer );
}

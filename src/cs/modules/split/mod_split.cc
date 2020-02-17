/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSelectionManager.h"
#include "csSelection.h"

using namespace cseis_geolib;
using namespace cseis_system;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPLIT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_split {
  struct VariableStruct {
    csSelectionManager* selectionManager;
    int  outputPort;
    bool acceptAll;
    bool isTraceStored;
  };
}
using namespace mod_split;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_split_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  //  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->outputPort    = 0;
  vars->acceptAll     = false;
  vars->isTraceStored = false;

  if( param->numParameters() == 0 ) {
    vars->acceptAll = true;
  }
  else {
    csVector<std::string> valueList;
    
    param->getAll( "header", &valueList );
    if( valueList.size() == 0 ) {
      writer->warning("%s: Wrong number of parameters for option 'header'. Expected: > 0, found: %d.", edef->moduleName().c_str(), valueList.size());
      env->addError();
    }

    std::string text;
    param->getString( "select", &text );

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

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_split_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{

  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csTraceHeaderDef const* hdef = env->headerDef;


  int numTraces = traceGather->numTraces();
  if( numTraces == 0 ) return;

  if( edef->isDebug() ) {
    int trcno_id = hdef->headerIndex("trcno");
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      fprintf(stderr,"Trace %d, trcno: %d\n", itrc, traceGather->trace(itrc)->getTraceHeader()->intValue(trcno_id) );
    }
  }

  csTrace* traceIn = traceGather->trace(0);
  if( vars->isTraceStored ) {
    if( edef->isDebug() ) {
      int trcno = traceGather->trace(0)->getTraceHeader()->intValue(hdef->headerIndex("trcno"));
      fprintf(stderr,"SECOND CALL %d, trcno: %d\n", numTraces, trcno);
    }
    *port = vars->outputPort;
    vars->isTraceStored = false;
    return;
  }
  else {
    bool isContained = true;
    if( !vars->acceptAll ) {
      isContained = vars->selectionManager->contains( traceIn->getTraceHeader() );
    }
    if( isContained ) {
      if( edef->isDebug() ) {
        int trcno = traceGather->trace(0)->getTraceHeader()->intValue(hdef->headerIndex("trcno"));
        fprintf(stderr,"FIRST CALL %d, trcno: %d\n", numTraces, trcno);
      }
      traceGather->copyTraceTo( 0, traceGather );
      *port = 0;
      vars->outputPort = 1;
      vars->isTraceStored = true;
      *numTrcToKeep = 1;
    }
    else {
      if( edef->isDebug() ) fprintf(stderr,"NOT CONTAINED\n");
      *port = 1;
    }
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_split_( csParamDef* pdef ) {
  pdef->setModule( "SPLIT", "Split/branch trace flow", "Create copy of traces that match specified header selection and process them within SPLIT-ENDSPLIT block" );

  pdef->addParam( "header", "Names of trace headers used for trace selection", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING,
    "List of selection strings, one for each specified header. See documentation for more detailed description of selection syntax" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_split_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_split::VariableStruct* vars = reinterpret_cast<mod_split::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_split_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_split::VariableStruct* vars = reinterpret_cast<mod_split::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->selectionManager ) {
    delete vars->selectionManager;
    vars->selectionManager = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_split_( csParamDef* pdef ) {
  params_mod_split_( pdef );
}
extern "C" void _init_mod_split_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_split_( param, env, writer );
}
extern "C" bool _start_exec_mod_split_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_split_( env, writer );
}
extern "C" void _exec_mod_split_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_split_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_split_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_split_( env, writer );
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSelectionManager.h"
#include "csSelection.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: KILL
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_kill {
  struct VariableStruct {
    std::string* hdrNames;
    csSelectionManager* selectionManager;
    int nHeaders;
    bool isModeInclude;
    bool killZeroTraces;
    bool killOneZeroTrace;
  };
}
using mod_kill::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_kill_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->killZeroTraces  = false;
  vars->killOneZeroTrace = false;
  vars->selectionManager = NULL;
  vars->hdrNames = NULL;
  vars->nHeaders = 0;
  vars->isModeInclude = true;

  //---------------------------------------------------------
  std::string text;

  if( param->exists( "zero_traces" ) ) {
    param->getString( "zero_traces", &text );
    if( !text.compare("yes") ) {
      vars->killZeroTraces = true;
    }
    else if( !text.compare("no") ) {
      vars->killZeroTraces = false;
    }
    else if( !text.compare("one") ) {
      vars->killZeroTraces = true;
      vars->killOneZeroTrace = true;
    }
    else {
      writer->error( "Option not recognised: '%s'", text.c_str() );
    }
  }

  if( !vars->killZeroTraces ) {
    //---------------------------------------------------------
    // Create new headers
    int nLines = param->getNumLines( "header" );
    if( nLines > 1 ) {
      writer->line( "More than one line encountered for user parameter '%s'. Only one line is supported.", "header" );
      env->addError();
    }

    csVector<std::string> valueList;
    param->getAll( "header", &valueList );

    if( valueList.size() == 0 ) {
      writer->line("Wrong number of arguments for user parameter '%s'. Expected: >0, found: %d.", "header", valueList.size());
      env->addError();
    }
    else {
      vars->nHeaders = valueList.size();
      param->getString( "select", &text );

      vars->selectionManager = NULL;
      try {
        vars->selectionManager = new csSelectionManager();
        vars->selectionManager->set( &valueList, &text, hdef );
      }
      catch( csException& e ) {
        vars->selectionManager = NULL;
        writer->error( "%s", e.getMessage() );
      }
      if( edef->isDebug() ) vars->selectionManager->dump();
    }
  }

  if( param->exists( "mode" ) ) {
    param->getString( "mode", &text );
    if( !text.compare("include") ) {
      vars->isModeInclude = true;
    }
    else if( !text.compare("exclude") ) {
      vars->isModeInclude = false;
    }
    else {
      writer->error( "Mode option not recognised: '%s'", text.c_str() );
    }
  }

  if( edef->isDebug() ) writer->line("Mode include: %d", vars->isModeInclude );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_kill_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);

  bool returnValue = true;

  if( !vars->killZeroTraces ) {
    if( vars->selectionManager->contains( trace->getTraceHeader() ) ) {
      returnValue = !vars->isModeInclude;
    }
    else {
      returnValue = vars->isModeInclude;
    }
  }
  else {
    float* samples = trace->getTraceSamples();
    if( !vars->killOneZeroTrace ) {
      returnValue = !vars->isModeInclude;
      for( int i = 0; i < shdr->numSamples; i++ ) {
        if( samples[i] != 0.0 ) {
          returnValue = vars->isModeInclude;
          break;
        }
      }
    }      
    else {
      returnValue = vars->isModeInclude;
      for( int i = 0; i < shdr->numSamples; i++ ) {
        if( samples[i] == 0.0 ) {
          returnValue = !vars->isModeInclude;// false;
          break;
        }
      }
    }
  }
  if( !returnValue ) {
    traceGather->freeAllTraces();
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//
//*************************************************************************************************
void params_mod_kill_( csParamDef* pdef ) {
  pdef->setModule( "KILL", "Kill traces", "Kill traces that match specified header selection" );

  pdef->addParam( "header", "Names of trace headers used for trace selection", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING,
      "List of selection strings, one for each specified header. See documentation for more detailed description of selection syntax" );

  pdef->addParam( "mode", "Mode of selection", NUM_VALUES_FIXED );
  pdef->addValue( "include", VALTYPE_OPTION );
  pdef->addOption( "include", "Kill traces matching the specified selection criteria" );
  pdef->addOption( "exclude", "Kill traces NOT matching the specified selection criteria" );

  pdef->addParam( "zero_traces", "Kill zero traces", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not kill zero traces" );
  pdef->addOption( "yes", "Kill zero traces" );
  pdef->addOption( "one", "Kill traces that contains at least one zero" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_kill_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_kill::VariableStruct* vars = reinterpret_cast<mod_kill::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_kill_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_kill::VariableStruct* vars = reinterpret_cast<mod_kill::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->selectionManager ) {
    delete vars->selectionManager; vars->selectionManager = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_kill_( csParamDef* pdef ) {
  params_mod_kill_( pdef );
}
extern "C" void _init_mod_kill_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_kill_( param, env, writer );
}
extern "C" bool _start_exec_mod_kill_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_kill_( env, writer );
}
extern "C" void _exec_mod_kill_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_kill_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_kill_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_kill_( env, writer );
}

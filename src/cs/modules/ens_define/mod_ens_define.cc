/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ENS_DEFINE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_ens_define {
  struct VariableStruct {
    int hdrID_all;
  }; 
}
using namespace mod_ens_define;

//*******************************************************************
//
//  Init phase
//
//*******************************************************************

void init_mod_ens_define_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->hdrID_all = -1;

  shdr->clearEnsembleKeys();

  if( param->exists("all") ) {
    string text;
    param->getString("all", &text);
    if( !text.compare("no") ) {
      // Nothing
    }
    else if( !text.compare("yes") ) {
      if( hdef->headerExists("all") ) {
        type_t type = hdef->headerType("all");
        if( type != TYPE_INT ) {
          writer->error("Trace header 'all' already exists but has wrong type. Should be integer");
        }
        else {
          writer->warning("Trace header 'all' already exists and will be overwritten");
        }
      }
      else {
        hdef->addHeader(TYPE_INT,"all","Created by ENS_DEFINE: Defines all traces as one ensemble");
      }
      vars->hdrID_all = hdef->headerIndex("all");
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
    shdr->setEnsembleKey( "all" );
  }

  if( param->exists("header") ) {
    csVector<std::string> valueList;
    param->getAll( "header", &valueList );
    if( vars->hdrID_all >= 0 ) {
      writer->line("Error: User parameter 'header' cannot be used in conjunction with user parameter 'all'. Specify one of these two parameters only.");
      env->addError();
    }
    else {
      for( int i = 0; i < valueList.size(); i++ ) {
        string name = valueList.at(i);
        if( !hdef->headerExists( name ) ) {
          writer->line("Error: Unknown trace header '%s'", name.c_str() );
          env->addError();
        }
        shdr->setEnsembleKey( name, i );
      }
    }
  }
  else if( !param->exists("all") ) {
    writer->error("No user parameter specified: 'header' or 'all'");
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_ens_define_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csTrace* trace = traceGather->trace(0);


  if( vars->hdrID_all >= 0 ) {
    trace->getTraceHeader()->setIntValue( vars->hdrID_all, 1 );
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_ens_define_( csParamDef* pdef ) {
  pdef->setModule( "ENS_DEFINE", "Define ensemble trace headers",
    "An 'ensemble' consists of all consecutive traces for which the ensemble trace headers do not change" );

  pdef->addParam( "header", "List of trace header names", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "all", "Define all input traces as one ensemble", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Define all input traces as one ensemble. This is done by adding a new integer trace header called 'all' which is set to '1'", "If trace header 'all' already exists, a warning is displayed in the log file" );
  pdef->addOption( "no", "Do not define all input traces as one ensemble" );

}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_ens_define_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_ens_define::VariableStruct* vars = reinterpret_cast<mod_ens_define::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_ens_define_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_ens_define::VariableStruct* vars = reinterpret_cast<mod_ens_define::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_ens_define_( csParamDef* pdef ) {
  params_mod_ens_define_( pdef );
}
extern "C" void _init_mod_ens_define_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_ens_define_( param, env, writer );
}
extern "C" bool _start_exec_mod_ens_define_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_ens_define_( env, writer );
}
extern "C" void _exec_mod_ens_define_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_ens_define_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_ens_define_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_ens_define_( env, writer );
}

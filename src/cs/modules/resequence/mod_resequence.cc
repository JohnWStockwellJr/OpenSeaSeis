/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: RESEQUENCE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_resequence {
  struct VariableStruct {
    long traceCounter;
    int  hdrId;
    char hdrType;
    int mode;
    int nTraces;
    double valueStart;
    double valueInc;
    double groupInc;
    double valueCurrent;
    int valueStartInt;
    int valueIncInt;
    int groupIncInt;
    int valueCurrentInt;

    int groupCurrent;
  };
  static const int MODE_ENSEMBLE = 11;
  static const int MODE_ALL      = 12;
  static const int MODE_FIXED    = 13;
}
using namespace mod_resequence;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_resequence_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );


  vars->traceCounter = 0;
  vars->hdrId        = -1;
  vars->hdrType      = 0;
  vars->mode         = 0;
  vars->nTraces      = 0;
  vars->valueStart   = 0;
  vars->valueInc     = 0;
  vars->groupInc     = 0;
  vars->valueCurrent  = 0;
  vars->valueStartInt = 0;
  vars->valueIncInt   = 0;
  vars->groupIncInt   = 0;
  vars->valueCurrentInt = 0;
  vars->groupCurrent  = 0;

  //---------------------------------------------------------
  std::string text;
  param->getString( "header", &text );

  vars->hdrType = 0;
  if( hdef->headerExists( text ) ) {
    vars->hdrId   = hdef->headerIndex( text );
    vars->hdrType = hdef->headerType( text );
  }
  else {
    writer->warning("Unknown trace header: %s", text.c_str());
    env->addError();
  }
  //---------------------------------------------------------
  
  int nValues;
  nValues = param->getNumValues( "set" );
  if( nValues >= 2 ) {
    if( vars->hdrType == TYPE_INT ) {
      param->getInt( "set", &vars->valueStartInt, 0 );
      param->getInt( "set", &vars->valueIncInt, 1 );
      vars->valueCurrentInt = vars->valueStartInt;
      vars->groupIncInt = 0;
      if( nValues == 3 ) {
        param->getInt( "set", &vars->groupIncInt, 2 );
      }
    }
    else if( vars->hdrType == TYPE_FLOAT || vars->hdrType == TYPE_DOUBLE || vars->hdrType == TYPE_INT64 ) {
      param->getDouble( "set", &vars->valueStart, 0 );
      param->getDouble( "set", &vars->valueInc, 1 );
      vars->valueCurrent = vars->valueStart;
      vars->groupInc = 0;
      if( nValues == 3 ) {
        param->getDouble( "set", &vars->groupInc, 2 );
      }
    }
    else {
      writer->line("Trace header used for resequencing must be of number type.");
      env->addError();
    }
  }
  else {
    writer->line("Error: User parameter 'set': Expected 2 values, found %d.", nValues);
    env->addError();
  }

  //---------------------------------------------------------
  param->getString("mode", &text);
  if( !text.compare("ensemble") ) {
    vars->mode = MODE_ENSEMBLE;
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }
  else if( !text.compare("all") ) {
    vars->mode = MODE_ALL;
    edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }
  else if( !text.compare("fixed") ) {
    vars->mode = MODE_FIXED;
    param->getInt( "ntraces", &vars->nTraces );
    if( vars->nTraces <= 0 ) {
      writer->error("Number of traces lower or equal to zero: %d", vars->nTraces);
    }
    edef->setTraceSelectionMode( TRCMODE_FIXED, vars->nTraces );
  }
  else {
    writer->error("Unknown option: '%s'", text.c_str());
  }

  vars->traceCounter = 0;
  vars->groupCurrent = 0;
  
  if( edef->isDebug() ) {
    writer->line("Mode: %d", vars->mode);
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_resequence_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;


  int nTraces = traceGather->numTraces();
  if( edef->isDebug() ) writer->line("Number of input traces: %d, current group: %d", traceGather->numTraces(), vars->groupCurrent );
  if( vars->mode == MODE_ENSEMBLE || vars->mode == MODE_FIXED ) {
    if( vars->hdrType == TYPE_INT ) {
      int valueStart = vars->valueStartInt + vars->groupIncInt * vars->groupCurrent;
      for( int itrc = 0; itrc < nTraces; itrc++ ) {
        csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
        int value = valueStart + itrc * vars->valueIncInt;
        trcHdr->setIntValue( vars->hdrId, value );
      }
    }
    else if( vars->hdrType == TYPE_INT64 ) {
      double valueStart = vars->valueStart + vars->groupInc * (double)vars->groupCurrent;
      for( int itrc = 0; itrc < nTraces; itrc++ ) {
        csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
        double value = valueStart + itrc * vars->valueInc;
        trcHdr->setInt64Value( vars->hdrId, (csInt64_t)value );
      }
    }
    else {
      double valueStart = vars->valueStart + vars->groupInc * (double)vars->groupCurrent;
      for( int itrc = 0; itrc < nTraces; itrc++ ) {
        csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
        double value = valueStart + itrc * vars->valueInc;
        trcHdr->setDoubleValue( vars->hdrId, value );
      }
    }
    vars->groupCurrent += 1;
  }
  else {  // MODE_ALL
    vars->traceCounter += 1;
    csTraceHeader* trcHdr = traceGather->trace(0)->getTraceHeader();
    if( vars->hdrType == TYPE_INT ) {
      trcHdr->setIntValue( vars->hdrId, vars->valueCurrentInt );
      vars->valueCurrentInt += vars->valueIncInt;
    }
    else if( vars->hdrType == TYPE_INT ) {
      trcHdr->setInt64Value( vars->hdrId, (csInt64_t)vars->valueCurrent );
      vars->valueCurrent += vars->valueInc;
    }
    else {  // DOUBLE
      trcHdr->setDoubleValue( vars->hdrId, vars->valueCurrent );
      vars->valueCurrent += vars->valueInc;
    }
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_resequence_( csParamDef* pdef ) {
  pdef->setModule( "RESEQUENCE", "Resequence trace header", "Reset a trace header to increasing (or decreasing) values for each trace group or ensemble" );

  pdef->addParam( "header", "Name of trace header to resequence", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "set", "Value to set trace header to.", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Start value" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Value increment", "Amount by which value is incremented for each consecutive trace" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Group increment", "Amount by which value is incremented for each consecutive trace group or ensemble" );

  pdef->addParam( "mode", "Mode of resequencing", NUM_VALUES_FIXED );
  pdef->addValue( "all", VALTYPE_OPTION );
  pdef->addOption( "ensemble", "Resequence inside each ensemble" );
  pdef->addOption( "all", "Resequence all traces" );
  pdef->addOption( "fixed", "Resequence fixed number of traces (specified in 'ntraces')" );

  pdef->addParam( "ntraces", "Number of consecutive traces to resequence", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Number of traces" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_resequence_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_resequence::VariableStruct* vars = reinterpret_cast<mod_resequence::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_resequence_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_resequence::VariableStruct* vars = reinterpret_cast<mod_resequence::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_resequence_( csParamDef* pdef ) {
  params_mod_resequence_( pdef );
}
extern "C" void _init_mod_resequence_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_resequence_( param, env, writer );
}
extern "C" bool _start_exec_mod_resequence_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_resequence_( env, writer );
}
extern "C" void _exec_mod_resequence_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_resequence_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_resequence_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_resequence_( env, writer );
}

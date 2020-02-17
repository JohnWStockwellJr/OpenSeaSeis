/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: REPEAT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_repeat {
  struct VariableStruct {
    int hdrId_repeat;
    int repeat;
    int mode;
    int traceCounter;

    int numTracesOrig;
    int newTraceCountdown;
  };
  static int const MODE_ENSEMBLE = 11;
  static int const MODE_TRACE    = 12;
  static int const MODE_ALL      = 13;
}
using namespace mod_repeat;

//*******************************************************************
//
//  Init phase
//
//*******************************************************************

void init_mod_repeat_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );


  vars->hdrId_repeat = -1;
  vars->repeat       = 0;
  vars->mode         = 0;
  vars->traceCounter = 0;

  vars->numTracesOrig = 0;
  vars->newTraceCountdown  = 0;

  if( param->exists( "mode" ) ) {
    std::string text;
    param->getString( "mode", &text );
    text = toLowerCase( text );
    if( !text.compare( "ensemble" ) ) {
      vars->mode = MODE_ENSEMBLE;
      env->execPhaseDef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
    }
    else if( !text.compare( "trace" ) ) {
      vars->mode = MODE_TRACE;
      env->execPhaseDef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
//      writer->line("Mode option not supported yet in this version: '%s'.", text.c_str());
//      env->addError();
    }
    else if( !text.compare( "all" ) ) {
      vars->mode = MODE_ALL;
      writer->line("Mode option not supported yet in this version: '%s'.", text.c_str());
      env->addError();
    }
    else {
      writer->line("Unknown argument for user parameter 'mode': '%s'.", text.c_str());
      env->addError();
    }
  }
  else {
    vars->mode = MODE_TRACE;
    env->execPhaseDef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }

  param->getInt( "repeat", &vars->repeat );
  if( vars->repeat <= 1 ) {
    writer->line("Wrong entry. Number of repeats must be larger than 1. Number specified: %d.", vars->repeat);
    env->addError();
  }
  if( !hdef->headerExists( "repeat" ) ) {
    vars->hdrId_repeat = hdef->addHeader( TYPE_INT, "repeat", "Repeat index" );
  }
  else {
    if( hdef->headerType("repeat") != TYPE_INT ) {
      writer->line("Trace header 'repeat' exists but has wrong header type. Must be INT.");
      env->addError();
    }
    else {
      vars->hdrId_repeat = hdef->headerIndex("repeat");
    }
  }
  vars->traceCounter = 0;

  if( edef->isDebug() ) {
    writer->line("Mode: %d", vars->mode );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_repeat_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;


  if( vars->newTraceCountdown > 0 ) {
    // Do nothing
  }
  else if( vars->numTracesOrig == 0 ) {
    vars->numTracesOrig = traceGather->numTraces();
    vars->newTraceCountdown = vars->numTracesOrig * ( vars->repeat - 1);

    // Set header 'repeat' for original traces to maximum repeat index, because
    // these traces will be output last
    for( int itrc = 0; itrc < vars->numTracesOrig; itrc++ ) {
      csTraceHeader* trcHdrOrig = traceGather->trace(itrc)->getTraceHeader();
      trcHdrOrig->setIntValue( vars->hdrId_repeat, vars->repeat );
    }
  }
  else {
    vars->numTracesOrig = 0;
    *numTrcToKeep = 0;
    return;
  }

  int traceIndexNew = 0;  // New trace is always the first one in gather

  traceGather->createTrace( traceIndexNew, hdef, shdr->numSamples );

  int repeat        = vars->repeat - (int)( (vars->newTraceCountdown-1)/ vars->numTracesOrig ) - 1;
  int traceIndex    = vars->numTracesOrig - (int)( (vars->newTraceCountdown-1) % vars->numTracesOrig );
//  fprintf(stdout,"Trace: %d / %d    %d, old trace %d\n", vars->numTracesOrig, vars->newTraceCountdown, repeat, traceIndex);
//  fflush(stdout);
  vars->newTraceCountdown -= 1;

  csTraceData* trcDataOrig  = traceGather->trace(traceIndex)->getTraceDataObject();
  csTraceHeader* trcHdrOrig = traceGather->trace(traceIndex)->getTraceHeader();
  csTraceData* trcDataCopy  = traceGather->trace(traceIndexNew)->getTraceDataObject();
  csTraceHeader* trcHdrCopy = traceGather->trace(traceIndexNew)->getTraceHeader();

  // Copy seismic & header data
  trcDataCopy->setData( trcDataOrig );
  trcHdrCopy->copyFrom( trcHdrOrig );
  // Set value for trace header 'repeat'
  trcHdrCopy->setIntValue( vars->hdrId_repeat, repeat );

  *numTrcToKeep = vars->numTracesOrig;   // Keep all original traces, only pass newly created trace

  //------------------------------------------------------------------------------------------
  // This module could be made more efficient in its memory requirements if
  // only one trace was allocated each time the module is called.
  // This would require some more logistics, by keeping one trace (or gather in case of 'ensemble' mode) (numTrcKeep)
  // until all requested traces have been output. At last, output the input trace(s) and set numTrcToKeep=0.
  // This will decrease the amount of traces to be allocated by 1/numRepeatedTraces
  //------------------------------------------------------------------------------------------
/*
  int numTracesOrig = traceGather->numTraces();
  int numTracesNew  = numTracesOrig * ( vars->repeat - 1);

  if( edef->isDebug() ) {
    writer->line("-----------------------------\n");
    writer->line("numTraces (orig,new): %6d %6d", numTracesOrig, numTracesNew );
    for( int i = 0; i < numTracesOrig; i++ ) {
      csTraceHeader* trcHdrOrig = traceGather->trace(i)->getTraceHeader();
      writer->line("   sensor %5d, source %5d,  trace: %d ",
        trcHdrOrig->intValue( hdef->headerIndex("sensor") ),
        trcHdrOrig->intValue( hdef->headerIndex("source") ), ++vars->traceCounter );
    }
  }

  traceGather->createTraces( numTracesOrig, numTracesNew, hdef, shdr->numSamples );

  for( int itrc = 0; itrc < numTracesOrig; itrc++ ) {
    csTraceData* trcDataOrig  = traceGather->trace(itrc)->getTraceDataObject();
    csTraceHeader* trcHdrOrig = traceGather->trace(itrc)->getTraceHeader();
    // Set value for trace header 'repeat'
    trcHdrOrig->setIntValue( vars->hdrId_repeat, 1 );

    for( int irepeat = 1; irepeat < vars->repeat; irepeat++ ) {
      int trcIndexNew = itrc + irepeat*numTracesOrig;
      csTraceData* trcDataCopy  = traceGather->trace(trcIndexNew)->getTraceDataObject();
      csTraceHeader* trcHdrCopy = traceGather->trace(trcIndexNew)->getTraceHeader();

      // Copy seismic & header data
      trcDataCopy->setData( trcDataOrig );
      trcHdrCopy->copyFrom( trcHdrOrig );
      // Set value for trace header 'repeat'
      trcHdrCopy->setIntValue( vars->hdrId_repeat, irepeat+1 );
    }
  }

  if( edef->isDebug() ) {
    for( int i = 0; i < numTracesNew; i++ ) {
      csTraceHeader* trcHdrNew = traceGather->trace(i+numTracesOrig)->getTraceHeader();
      writer->line("NEW sensor %5d, source %5d,  trace: %d ",
        trcHdrNew->intValue( hdef->headerIndex("sensor") ),
        trcHdrNew->intValue( hdef->headerIndex("source") ), ++vars->traceCounter );
    }
  }
*/
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_repeat_( csParamDef* pdef ) {
  pdef->setModule( "REPEAT", "Repeat/duplicate traces", "Copy traces and all trace headers. Set new trace header called 'repeat'" );

  pdef->addParam( "mode", "Mode of operation", NUM_VALUES_FIXED );
  pdef->addValue( "trace", VALTYPE_OPTION );
  pdef->addOption( "trace", "Repeat individual trace, trace-by-trace" );
  pdef->addOption( "ensemble", "Repeat whole ensemble, ensemble-by-ensemble" );
  pdef->addOption( "all", "Repeat all traces: Collect all traces as one ensemble, then repeat" );

  pdef->addParam( "repeat", "Total number of output traces per input trace", NUM_VALUES_FIXED, "repeat = 2 --> Output input trace and one copy" );
  pdef->addValue( "2", VALTYPE_NUMBER );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_repeat_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_repeat::VariableStruct* vars = reinterpret_cast<mod_repeat::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_repeat_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_repeat::VariableStruct* vars = reinterpret_cast<mod_repeat::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_repeat_( csParamDef* pdef ) {
  params_mod_repeat_( pdef );
}
extern "C" void _init_mod_repeat_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_repeat_( param, env, writer );
}
extern "C" bool _start_exec_mod_repeat_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_repeat_( env, writer );
}
extern "C" void _exec_mod_repeat_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_repeat_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_repeat_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_repeat_( env, writer );
}

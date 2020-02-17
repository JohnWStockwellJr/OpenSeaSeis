/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csStandardHeaders.h"
#include "csVector.h"
#include "csStackUtil.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: STACK
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_stack {
  struct VariableStruct {
    int mode;
    csTraceGather* stackedTraces;
    cseis_geolib::csVector<cseis_system::csTrace*>* stackTraceList;
    cseis_geolib::csVector<double>* hdrValueList;  // List of stack header values for all stacked/buffered traces so far
    cseis_geolib::csVector<int>* numStackedTracesList;
    int hdrId_stack;
    int hdrId_fold;
    float normFactor;
    cseis_system::csStackUtil* stackUtil;
    bool isFirstCall;
    bool doResetFold;
    bool foldHdrExists;
  }; 
  static int const MODE_ENSEMBLE = 11;
  static int const MODE_SORTED   = 12;
  static int const MODE_UNSORTED = 13;
  static int const MODE_ALL      = 14;
}
using namespace mod_stack;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_stack_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );


  vars->mode              = -1;
  vars->stackedTraces     = NULL;
  vars->stackTraceList    = NULL;
  vars->stackUtil         = NULL;
//  vars->stackBufferList = NULL;
  vars->numStackedTracesList = NULL;
  vars->hdrValueList    = NULL;
  vars->hdrId_stack     = -1;
  vars->hdrId_fold      = -1;
  vars->normFactor      = 1;
  vars->isFirstCall     = true;
  vars->doResetFold = true;
  vars->foldHdrExists = false;

  float zeroThreshold = 0;

  if( param->exists("zero_threshold") ) {
    param->getFloat("zero_threshold", &zeroThreshold );
    zeroThreshold = abs(zeroThreshold);
  }

  std::string text;
  int outputOption = csStackUtil::OUTPUT_FIRST;

  if( param->exists("output_hdr") ) {
    param->getString("output_hdr", &text);
    if( !text.compare("first") ) {
      outputOption = csStackUtil::OUTPUT_FIRST;
    }
    else if( !text.compare("last") ) {
      outputOption = csStackUtil::OUTPUT_LAST;
    }
    else if( !text.compare("average") ) {
//      writer->error("Option AVERAGE is not supported yet...");
      outputOption = csStackUtil::OUTPUT_AVERAGE;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
      env->addError();
    }
  }

  param->getString("mode", &text);
  if( !text.compare("ensemble") ) {
    vars->mode = MODE_ENSEMBLE;
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }
  else if( !text.compare("all") ) {
    vars->mode = MODE_ALL;
    edef->setTraceSelectionMode( TRCMODE_FIXED, 2 );
    vars->stackedTraces        = new csTraceGather();
    vars->numStackedTracesList = new csVector<int>();
    vars->numStackedTracesList->insert( 0 );
  }
  else {
    if( !text.compare("sorted") ) {
      vars->mode = MODE_SORTED;
    }
    else if( !text.compare("unsorted") ) {
      vars->mode = MODE_UNSORTED;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
    edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
    vars->stackedTraces        = new csTraceGather();
    vars->hdrValueList         = new csVector<double>();
    vars->numStackedTracesList = new csVector<int>();
  }

  bool normTimeVariant = false;
  bool outputNormTrace = false;
  if( param->exists("norm_time_variant") ) {
    param->getString("norm_time_variant", &text);
    if( !text.compare("no") ) {
      normTimeVariant = false;
    }
    else if( !text.compare("yes") ) {
      normTimeVariant = true;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
    if( param->getNumValues("norm_time_variant" ) > 1 ) {
      param->getString("norm_time_variant", &text, 1);
      if( !text.compare("norm_trace") ) {
        outputNormTrace = true;
      }
      else if( !text.compare("stack_trace") ) {
        outputNormTrace = false;
      }
      else {
        writer->line("Unknown option: '%s'", text.c_str());
      }
    }
  }

  if( param->exists("reset_fold") ) {
    param->getString("resetFold", &text);
    if( !text.compare("no") ) {
      vars->doResetFold = false;
    }
    else if( !text.compare("yes") ) {
      vars->doResetFold = true;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
  }

  if( param->exists("norm") ) {
    param->getFloat("norm", &vars->normFactor);
  }
  if( !hdef->headerExists( HDR_FOLD.name ) ) {
    vars->foldHdrExists = false;
    hdef->addStandardHeader( HDR_FOLD.name );
  }
  else {
    vars->foldHdrExists = true;
  }
  vars->hdrId_fold = hdef->headerIndex( HDR_FOLD.name );

  if( vars->mode != MODE_ENSEMBLE && vars->mode != MODE_ALL ) {
    param->getString("header", &text);
    vars->hdrId_stack = hdef->headerIndex( text );
  }

  vars->stackUtil = new csStackUtil( shdr->numSamples, vars->normFactor, outputOption );
  vars->stackUtil->setTimeVariantNorm( normTimeVariant, vars->hdrId_stack, zeroThreshold );
  vars->stackUtil->setOutputNormTrace( outputNormTrace );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_stack_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  //  csSuperHeader const*    shdr = env->superHeader;


  int numTracesIn = traceGather->numTraces();

  if( numTracesIn == 0 || edef->isLastCall() ) {
    // This is the very last call to this module --> output all buffered traces
    if( vars->mode == MODE_ALL ) {
      csTrace* traceOut = traceGather->trace(0);
      int numStackedTraces = vars->numStackedTracesList->at(0)+1;
      vars->stackUtil->normStackedTrace( traceOut, numStackedTraces );
      traceOut->getTraceHeader()->setIntValue( vars->hdrId_fold, numStackedTraces );
      vars->numStackedTracesList->clear();
      *numTrcToKeep = 0;
      return;
    }
    else if( vars->mode == MODE_SORTED ) {
      if( vars->stackedTraces->numTraces() == 0 || vars->numStackedTracesList->size() == 0 ) {
        // Must be a single input trace: Do nothing except setting fold header
        if( numTracesIn > 0 ) traceGather->trace(0)->getTraceHeader()->setIntValue( vars->hdrId_fold, 1 );
        return;
      }
      csTrace* traceOut = vars->stackedTraces->trace(0);
      int numStackedTraces = vars->numStackedTracesList->at(0);
      if( numTracesIn != 0 ) {
        csTrace* traceIn  = traceGather->trace(0);
        double hdrValueIn = traceIn->getTraceHeader()->doubleValue(vars->hdrId_stack);
        if( vars->hdrValueList->size() == 0 ) {
          writer->error("Program bug in module STACK: Header value list is zero.");
        }
        if( vars->hdrValueList->at(0) == hdrValueIn ) {  // Same header value? --> Stack input trace into buffered trace
          vars->stackUtil->stackTrace( traceOut, traceIn );
          vars->numStackedTracesList->set( numStackedTraces+1, 0 );
          traceGather->freeTrace( 0 );
          numStackedTraces += 1;
        }
        else {
          // Do nothing, just leave last trace as it is...
        }
      }
      vars->stackUtil->normStackedTrace( traceOut, numStackedTraces );
      traceOut->getTraceHeader()->setIntValue( vars->hdrId_fold, numStackedTraces );
      // Swap buffered and input trace
      vars->stackedTraces->moveTraceTo( 0, traceGather );  // Move buffered stacked trace to traceGather (new traceIndex '1')
      return;
    }
    else if( vars->mode == MODE_ENSEMBLE && numTracesIn != 0 ) {
      vars->stackUtil->stackAndNormalizeTraces( traceGather );
      traceGather->trace(0)->getTraceHeader()->setIntValue( vars->hdrId_fold, numTracesIn );
      return;
    }
    else if( vars->mode == MODE_UNSORTED && numTracesIn == 0 ) {
      int numTracesOut = vars->stackedTraces->numTraces();
      vars->stackedTraces->moveTracesTo( 0, numTracesOut, traceGather );
      for( int itrcOut = 0; itrcOut < numTracesOut; itrcOut++ ) {
        csTrace* traceOut    = traceGather->trace(itrcOut);
        int numStackedTraces = vars->numStackedTracesList->at(itrcOut);
        vars->stackUtil->normStackedTrace( traceOut, numStackedTraces );
        traceOut->getTraceHeader()->setIntValue( vars->hdrId_fold, numStackedTraces );
      }
      vars->hdrValueList->clear();
      vars->numStackedTracesList->clear();
      return;
    }
  }

  //---------------------------------------------------------------------
  // Ensemble mode
  // Simply stack all input traces into one stacked output trace
  //
  if( vars->mode == MODE_ENSEMBLE ) {
    vars->stackUtil->stackAndNormalizeTraces( traceGather );
    traceGather->trace(0)->getTraceHeader()->setIntValue( vars->hdrId_fold, numTracesIn );
  }
  //---------------------------------------------------------------------
  // More complicated modes: Stacked trace(s) have to be buffered before outputting
  //
  else if( vars->mode == MODE_ALL ) {
    csTrace* traceOut = traceGather->trace(0);
    csTrace* traceIn  = traceGather->trace(1);
    vars->stackUtil->stackTrace( traceOut, traceIn );      
    traceGather->freeTraces(1,1);
    vars->numStackedTracesList->set( vars->numStackedTracesList->at(0)+1, 0 );
    *numTrcToKeep = 1;
  }
  else {
    if( numTracesIn > 1 ) {
      writer->error("STACK: Wrong number of input traces: %d (expected: 1). This is probably due to a program bug within this module.", numTracesIn);
    }

    // One new input trace to stack...
    csTrace* traceIn  = traceGather->trace(0);
    double hdrValueIn = 0.0;
    if( vars->mode != MODE_ALL ) {
      hdrValueIn = traceIn->getTraceHeader()->doubleValue(vars->hdrId_stack);
    }
    // This is the very first trace. --> Save trace in stacked trace buffer
    if( vars->stackedTraces->numTraces() == 0 ) {
      traceGather->moveTraceTo( 0, vars->stackedTraces );
      vars->stackUtil->stackTrace( vars->stackedTraces->trace(0) );  // Bug fix 150528: Stack in new trace only once (was done twice)
      vars->numStackedTracesList->insertEnd( 1 );
      if( vars->mode != MODE_ALL ) {
        vars->hdrValueList->insertEnd( hdrValueIn );
      }
    }
    // Input data is sorted --> if new stack ensemble is found, output current stack buffer first, otherwise just stack into buffer
    else if( vars->mode == MODE_SORTED ) {
      csTrace* traceOut = vars->stackedTraces->trace(0);
      int numStackedTraces = vars->numStackedTracesList->at(0);
      if( vars->hdrValueList->at(0) == hdrValueIn ) {  // Same header value? --> Stack input trace into buffered trace
        vars->stackUtil->stackTrace( traceOut, traceIn );
        vars->numStackedTracesList->set( numStackedTraces+1, 0 );
        traceGather->freeTrace( 0 );
      }
      else { // New header value. --> Normalise and output current stacked trace, buffer new trace
        vars->stackUtil->normStackedTrace( traceOut, numStackedTraces );
        traceOut->getTraceHeader()->setIntValue( vars->hdrId_fold, numStackedTraces );
        // Swap buffered and input trace
        vars->stackedTraces->moveTraceTo( 0, traceGather );  // Move buffered stacked trace to traceGather (new traceIndex '1')
        traceGather->moveTraceTo( 0, vars->stackedTraces );  // Move input trace (traceGather traceIndex '0') to buffer
        vars->hdrValueList->set( hdrValueIn, 0 );
        vars->numStackedTracesList->set( 1, 0 );
        // Stack in new trace:
        vars->stackUtil->stackTrace( vars->stackedTraces->trace(0) );  // Bug fix 150528: Stack in new trace only once (was done twice)
      }
    }
    //----------------------------------------------------------------------------
    else {  // MODE_UNSORTED && MODE_ALL
      int stackedTraceIndex = 0;
      if( vars->mode != MODE_ALL ) {
        stackedTraceIndex = -1;
        // Search for header value, if available
        for( int i = 0; i < vars->hdrValueList->size(); i++ ) {
          if( vars->hdrValueList->at(i) == hdrValueIn ) {
            stackedTraceIndex = i;
            break;
          }
        }
      }
      if( stackedTraceIndex < 0 ) {  // Header value not found --> add new stack trace to buffer
        traceGather->moveTraceTo( 0, vars->stackedTraces );
        vars->hdrValueList->insertEnd( hdrValueIn );
        vars->numStackedTracesList->insertEnd( 1 ); 
        vars->stackUtil->stackTrace( vars->stackedTraces->trace(vars->stackedTraces->numTraces()-1) );
      }
      else {  // Header value found --> stack input trace into buffer
        csTrace* traceOut = vars->stackedTraces->trace(stackedTraceIndex);
        vars->stackUtil->stackTrace( traceOut, traceIn );
        vars->numStackedTracesList->set( vars->numStackedTracesList->at(stackedTraceIndex)+1, stackedTraceIndex );

        traceGather->freeTrace( 0 );
      }
      if( edef->isLastCall() ) {
        int numTracesOut = vars->stackedTraces->numTraces();
        vars->stackedTraces->moveTracesTo( 0, numTracesOut, traceGather );
        for( int itrcOut = 0; itrcOut < numTracesOut; itrcOut++ ) {
          csTrace* traceOut    = traceGather->trace(itrcOut);
          int numStackedTraces = vars->numStackedTracesList->at(itrcOut);
          vars->stackUtil->normStackedTrace( traceOut, numStackedTraces );
          traceOut->getTraceHeader()->setIntValue( vars->hdrId_fold, numStackedTraces );
        }
        vars->hdrValueList->clear();
        vars->numStackedTracesList->clear();

        return;
      }
    }
    if( vars->numStackedTracesList->size() != 0 ) {
      edef->setTracesAreWaiting();
    }
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_stack_( csParamDef* pdef ) {
  pdef->setModule( "STACK", "Ensemble stack", "Stack all traces in ensemble" );

  pdef->addParam( "output_hdr", "How shall trace headers of output trace be determined?", NUM_VALUES_FIXED );
  pdef->addValue( "first", VALTYPE_OPTION );
  pdef->addOption( "first", "Use trace header values from first input trace" );
  pdef->addOption( "last", "Use trace header values from last input trace" );
  pdef->addOption( "average", "Take average of all input trace headers" );

  pdef->addParam( "mode", "Mode of operation", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "ensemble", "Stack each input ensemble." );
  pdef->addOption( "all", "Stack all incoming traces." );
  pdef->addOption( "sorted",   "Input data have already been pre-sorted by header(s) specified in parameter 'header'." );
  pdef->addOption( "unsorted", "Input data have NOT been pre-sorted.", "This means that the stack module waits until all traces have been input before outputting the first trace. Note that output traces may not be sorted in (increasing,decreasing) order of the stack header.");

  pdef->addParam( "header", "Stack header", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace 'stack' header name. Stack all traces with same stack header value." );

  pdef->addParam( "norm", "Normalisation factor", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Output stack value is normalised by number of stacked traces to the power of the 'norm' factor.", "Specify 0.5 for sqrt(N) normalization, 0 for no normalization" );

  pdef->addParam( "norm_time_variant", "Time variant normalisation?", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not apply time variant normalisation. Apply constant normalisation for each sample value" );
  pdef->addOption( "yes", "Apply time variant normalisation according to stacked non-zero samples.");
  pdef->addValue( "stack_trace", VALTYPE_OPTION );
  pdef->addOption( "stack_trace", "Do not output normalisation trace - output stacked trace" );
  pdef->addOption( "norm_trace", "Output normalisation trace instead of stacked data");

  pdef->addParam( "zero_threshold", "Threshold below which trace amplitude is considered zero", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Absolute zero threshold value. Values below this threshold are considered zero values"  );

  // Not fully implemented yet:
  //  pdef->addParam( "reset_fold", "Reset fold (trace header 'fold') of each input trace to zero?", NUM_VALUES_FIXED );
  //  pdef->addValue( "yes", VALTYPE_OPTION );
  //  pdef->addOption( "yes", "Reset fold of each input trace to zero");
  //  pdef->addOption( "no", "Do not reset fold of input traces. Accumulate fold" );
}



//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_stack_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_stack::VariableStruct* vars = reinterpret_cast<mod_stack::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_stack_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_stack::VariableStruct* vars = reinterpret_cast<mod_stack::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->stackedTraces != NULL) {
    delete vars->stackedTraces;
    vars->stackedTraces = NULL;
  }
  if( vars->stackTraceList != NULL ) {
    delete vars->stackTraceList;
    vars->stackTraceList = NULL;
  }
  if( vars->hdrValueList != NULL ) {
    delete vars->hdrValueList;
    vars->hdrValueList = NULL;
  }
  if( vars->numStackedTracesList != NULL ) {
    delete vars->numStackedTracesList;
    vars->numStackedTracesList = NULL;
  }
  if( vars->stackUtil != NULL ) {
    delete vars->stackUtil;
    vars->stackUtil = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_stack_( csParamDef* pdef ) {
  params_mod_stack_( pdef );
}
extern "C" void _init_mod_stack_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_stack_( param, env, writer );
}
extern "C" bool _start_exec_mod_stack_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_stack_( env, writer );
}
extern "C" void _exec_mod_stack_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_stack_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_stack_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_stack_( env, writer );
}

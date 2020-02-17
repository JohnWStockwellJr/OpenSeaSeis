/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/*
 * This module has a known bug when padding:
 * If header values outside specified range exists, a segmentation error may occur
 * Workaround: delete all traces outside specified pad range before calling this module
 */

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TRC_ADD_ENS
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_trc_add_ens {
  struct VariableStruct {
    int numTraces;
    float value;
    int method;
    int hdrID_pad;
    int hdrType_pad;
    int hdrID_padFlag;
    csFlexNumber start;
    csFlexNumber stop;
    csFlexNumber inc;
    bool deleteInconsistentTraces;
    bool isPadRange;
  };
  static const int METHOD_PAD     = 71;
  static const int METHOD_NTRACES = 72;
}
using mod_trc_add_ens::VariableStruct;

//*******************************************************************
//
//  Init phase
//
//*******************************************************************

void init_mod_trc_add_ens_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  env->execPhaseDef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->numTraces      = 0;
  vars->value          = 0;
  vars->start          = 0;
  vars->stop           = 0;
  vars->inc            = 0;
  vars->hdrID_padFlag  = -1;
  vars->hdrID_pad      = -1;
  vars->method         = mod_trc_add_ens::METHOD_NTRACES;
  vars->deleteInconsistentTraces = false;
  vars->isPadRange = false;
  
  //-----------------------------------------------
  string text;
  if( param->exists("method") ) {
    param->getString( "method", &text );
    if( !text.compare("ntraces") ) {
      vars->method = mod_trc_add_ens::METHOD_NTRACES;
    }
    else if( !text.compare("pad") ) {
      vars->method = mod_trc_add_ens::METHOD_PAD;
      if( param->exists("delete") ) {
        string text;
        param->getString("delete",&text);
        if( !text.compare("yes") ) {
          vars->deleteInconsistentTraces = true;
        }
        else if( !text.compare("no") ) {
          vars->deleteInconsistentTraces = false;
        }
        else {
          writer->error("Unknown option: %s", text.c_str());
        }
      }
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  //-----------------------------------------------
  if( vars->method == mod_trc_add_ens::METHOD_NTRACES ) {
    param->getInt( "ntraces", &vars->numTraces );
  }
  else { // method PAD
    param->getString( "header", &text );
    if( !hdef->headerExists(text) ) {
      writer->error("Trace header does not exist: %s", text.c_str() );
    }
    if( !hdef->headerExists("pad_flag") ) {
      hdef->addHeader(cseis_geolib::TYPE_INT, "pad_flag", "Padded trace? 0:no, 1:yes");
    }
    vars->hdrID_padFlag = hdef->headerIndex("pad_flag");
    vars->hdrID_pad   = hdef->headerIndex(text);
    vars->hdrType_pad = hdef->headerType(text);
    if( vars->hdrType_pad == TYPE_STRING || vars->hdrType_pad == TYPE_CHAR ) {
      writer->error("Trace header does not have a number type. This is currently not supported.");
    }
    else if( vars->hdrType_pad != TYPE_INT ) {
      writer->error("Trace header has floating point type, or 64bit integer. Only 32bit integer type trace headers are currently supported for trace padding.");
    }
    //    if( vars->hdrType_pad == TYPE_INT ) {
    int tmp3;
    param->getInt( "pad_inc", &tmp3 );
    vars->inc.setIntValue(tmp3);
    if( param->exists("pad_range") ) {
      vars->isPadRange = true;
      int tmp1, tmp2;
      param->getInt( "pad_range", &tmp1, 0 );
      vars->start.setIntValue(tmp1);
      param->getInt( "pad_range", &tmp2, 1 );
      vars->stop.setIntValue(tmp2);
      vars->numTraces = (tmp2 - tmp1) / tmp3 + 1;
    }
    if( vars->isPadRange ) writer->line("Number of output traces per ensemble: %d", vars->numTraces );
  }

  if( param->exists("value") ) {
    param->getFloat( "value", &vars->value );
  }

  if( edef->isDebug() ) {
    writer->line("Start/stop/inc: %f %f %f", vars->start.doubleValue(), vars->stop.doubleValue(), vars->inc.doubleValue() );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_trc_add_ens_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;

  int numTracesIn = traceGather->numTraces();
  if( numTracesIn == 0 ) return;

  //------------------------------------------------
  if( vars->method == mod_trc_add_ens::METHOD_NTRACES ) {
    traceGather->createTraces( numTracesIn, vars->numTraces, hdef, shdr->numSamples );

    csTraceHeader* trcHdrOrig = traceGather->trace(numTracesIn-1)->getTraceHeader();
    for( int itrc = 0; itrc < vars->numTraces; itrc++ ) {
      int trcIndexNew = numTracesIn + itrc;

      // Copy header data
      csTraceHeader* trcHdrCopy = traceGather->trace(trcIndexNew)->getTraceHeader();
      trcHdrCopy->copyFrom( trcHdrOrig );

      float* samples = traceGather->trace(trcIndexNew)->getTraceSamples();
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] = vars->value;
      }
    }
  }
  //------------------------------------------------
  else {  // method PAD
    // if( vars->hdrType == TYPE_INT {
    int increment  = vars->inc.intValue();
    int valueStart = traceGather->trace(0)->getTraceHeader()->intValue( vars->hdrID_pad );
    int valueEnd   = traceGather->trace(numTracesIn-1)->getTraceHeader()->intValue( vars->hdrID_pad );
    if( vars->isPadRange ) {
      valueStart = vars->start.intValue();
      valueEnd   = vars->stop.intValue();
    }

    // 1) Check for duplicate trace values --> Remove if specified by user
    int valuePrev = traceGather->trace( numTracesIn-1 )->getTraceHeader()->intValue(vars->hdrID_pad);
    for( int itrc = numTracesIn-2; itrc >= 0; itrc-- ) {
      csTraceHeader* trcHdrCurrent = traceGather->trace(itrc)->getTraceHeader();
      int valueCurrent = trcHdrCurrent->intValue( vars->hdrID_pad );
      if( valueCurrent == valuePrev ) {
        // Duplicate trace header value found --> Remove
        writer->warning("Duplicate trace header value: %d. Trace number: %d (out of %d)", valueCurrent, itrc+1, numTracesIn );
        traceGather->freeTrace( itrc );
      }
      else {
        int valSign = increment * (valuePrev - valueCurrent);
        if( valSign < 0 ) {
          writer->error("Inconsistent trace values: %d %d (traces #%d and #%d). Different sign of increment than specified 'pad_inc' (=%d)", valueCurrent, valuePrev, itrc+1, itrc+2, increment );
        }
        if( vars->isPadRange && (valueCurrent-valueStart) % increment != 0 ) {
          writer->error("Inconsistent trace header value: %d. Expected values should fit into the specified padding range frm %d to %d, increment %d. Input data must be sorted accordingly",
                        valueCurrent, valueStart, valueEnd, increment );
        }
        trcHdrCurrent->setIntValue( vars->hdrID_padFlag, 0 );
      }
      valuePrev = valueCurrent;
    }
    int numTracesInUnique = traceGather->numTraces();
    if( numTracesInUnique != numTracesIn && !vars->deleteInconsistentTraces ) {
      writer->error("Job aborted due to duplicate trace header values. To avoid job termination, specify user parameter 'delete yes'");
    }
    int numTracesOut = vars->isPadRange ? vars->numTraces : ( (valueEnd - valueStart)/increment + 1 );
    int numTracesToPad = numTracesOut - numTracesInUnique;
    if( edef->isDebug() ) {
      writer->line("Number of input traces/after deletion of duplicates/to be padded/total: %d / %d / %d / %d", numTracesIn, numTracesInUnique, numTracesToPad, numTracesOut );
    }
    if( numTracesToPad == 0 ) return;
    
    numTracesIn = numTracesInUnique;
    csTraceHeader* trcHdr0 = traceGather->trace(0)->getTraceHeader();
    for( int itrcOut = 0; itrcOut < numTracesOut; itrcOut++ ) {
      int valueOut = valueStart + itrcOut * increment;
      int valueCurrentTrace = -valueOut; // Make up value
      if( itrcOut < traceGather->numTraces() ) {
        valueCurrentTrace = traceGather->trace(itrcOut)->getTraceHeader()->intValue(vars->hdrID_pad);
      }
      if( valueOut != valueCurrentTrace ) {
        csTrace* traceNew = traceGather->createTrace( itrcOut, hdef, shdr->numSamples );
        csTraceHeader* trcHdrNew = traceNew->getTraceHeader();
        trcHdrNew->copyFrom( trcHdr0 );
        trcHdrNew->setIntValue( vars->hdrID_padFlag, 1 );
        trcHdrNew->setIntValue( vars->hdrID_pad, valueOut );
        float* samples = traceNew->getTraceSamples();
        for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
          samples[isamp] = vars->value;
        }
      }
    }
  } // END pad method
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_trc_add_ens_( csParamDef* pdef ) {
  pdef->setModule( "TRC_ADD_ENS", "Ensemble trace adding" );

  pdef->addParam( "method", "Method for trace adding", NUM_VALUES_FIXED );
  pdef->addValue( "ntraces", VALTYPE_OPTION );
  pdef->addOption( "ntraces", "Add specified number of traces at beginning and end of ensemble" );
  pdef->addOption( "pad", "Pad traces as specified by trace header name and padding values" );

  pdef->addParam( "ntraces", "Number of traces to add at end of each ensemble", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER );

  pdef->addParam( "value", "Initialize trace samples to the given value", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER );

  pdef->addParam( "header", "Trace header name for trace padding", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "pad_range", "Pad traces within the specified (enforced) range of trace header values", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "First value" );
  pdef->addValue( "", VALTYPE_NUMBER, "Last value" );

  pdef->addParam( "pad_inc", "Pad traces with the specified trace header value increment between consecutive traces", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Increment" );

  pdef->addParam( "delete", "Delete inconsistent or out-of-range traces when trace padding?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Delete inconsistent traces" );
  pdef->addOption( "no", "Do not delete inconsistent traces" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_trc_add_ens_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_trc_add_ens::VariableStruct* vars = reinterpret_cast<mod_trc_add_ens::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_trc_add_ens_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_trc_add_ens::VariableStruct* vars = reinterpret_cast<mod_trc_add_ens::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_trc_add_ens_( csParamDef* pdef ) {
  params_mod_trc_add_ens_( pdef );
}
extern "C" void _init_mod_trc_add_ens_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_trc_add_ens_( param, env, writer );
}
extern "C" bool _start_exec_mod_trc_add_ens_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_trc_add_ens_( env, writer );
}
extern "C" void _exec_mod_trc_add_ens_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_trc_add_ens_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_trc_add_ens_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_trc_add_ens_( env, writer );
}

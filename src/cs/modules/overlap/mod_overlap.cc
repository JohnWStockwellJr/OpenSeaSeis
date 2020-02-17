/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: OVERLAP
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_overlap {
  struct VariableStruct {
    int hdrID_time_samp1_s;
    int hdrID_time_samp1_us;
    bool honorAbsoluteTime;
    int timePrev_s;
    int timePrev_us;
    bool isFirstCall;
    int overlap_numSamples;
    int numSamplesIn;
  };
}
using namespace mod_overlap;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_overlap_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 2 );

  vars->hdrID_time_samp1_s  = -1;
  vars->hdrID_time_samp1_us = -1;
  vars->honorAbsoluteTime   = false;
  vars->isFirstCall         = true;
  vars->overlap_numSamples  = 0;
  vars->numSamplesIn        = shdr->numSamples;

  //---------------------------------------------------------
  if( param->exists("absolute_time") ) {
    std::string text;
    param->getString( "absolute_time", &text );
    if( !text.compare( "yes" ) ) {
      vars->honorAbsoluteTime = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->honorAbsoluteTime = false;
    }
    else {
      writer->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }
  //---------------------------------------------------------
  float overlap_ms;
  param->getFloat( "overlap", &overlap_ms );
  float traceLength_ms = (float)shdr->numSamples * shdr->sampleInt;
  if( overlap_ms <= 0 || overlap_ms > traceLength_ms ) {
    writer->line("Error in user parameter 'overlap': Inconsistent overlap specified: %fms. Overlap has to be given in [ms]. Current trace length: %fms", overlap_ms, traceLength_ms);
    env->addError();
  }

  vars->hdrID_time_samp1_s  = hdef->headerIndex( HDR_TIME_SAMP1.name );
  vars->hdrID_time_samp1_us = hdef->headerIndex( HDR_TIME_SAMP1_US.name );

  vars->overlap_numSamples = (int)round( overlap_ms / shdr->sampleInt );

  shdr->numSamples += 2*vars->overlap_numSamples;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_overlap_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;


  int numTracesIn = traceGather->numTraces();
  if( numTracesIn == 0 ) {
    return;
  }
  float* samples1 = traceGather->trace(0)->getTraceSamples();
  if( vars->isFirstCall ) {
    for( int isamp = vars->numSamplesIn-1; isamp >= 0; isamp-- ) {
      samples1[isamp+vars->overlap_numSamples] = samples1[isamp];
    }
    for( int isamp = 0; isamp < vars->overlap_numSamples; isamp++ ) {
      samples1[isamp] = 0.0f;
    }
    vars->isFirstCall = false;
    if( vars->honorAbsoluteTime ) {
      csTraceHeader* trcHdr = traceGather->trace(0)->getTraceHeader();
      vars->timePrev_s  = trcHdr->intValue( vars->hdrID_time_samp1_s );
      vars->timePrev_us = trcHdr->intValue( vars->hdrID_time_samp1_us );
    }
  }
  if( numTracesIn == 2 ) {
    bool doCopy = true;
    if( vars->honorAbsoluteTime ) {
      csTraceHeader* trcHdr = traceGather->trace(1)->getTraceHeader();
      int time_samp1_s  = trcHdr->intValue( vars->hdrID_time_samp1_s );
      int time_samp1_us = trcHdr->intValue( vars->hdrID_time_samp1_us );
      int diff_s  = time_samp1_s  - vars->timePrev_s;
      int diff_us = time_samp1_us - vars->timePrev_us;
      int diff_samples = (int)(round((double)diff_s*1000.0 + (double)diff_us/1000.0 )/((double)shdr->sampleInt) );
      vars->timePrev_s  = time_samp1_s;
      vars->timePrev_us = time_samp1_us;
      //      fprintf(stderr,"Diff samples: %d   %d %d   %d %d\n", diff_samples, time1_samp1_s, time2_samp1_s, time1_samp1_us, time2_samp1_us );
      if( diff_samples > vars->numSamplesIn ) {
        doCopy = false;
      }
    }
    float* samples2 = traceGather->trace(1)->getTraceSamples();
    // 1) Copy first # overlap samples from trace 2 to bottom of trace 1
    if( doCopy ) {
      memcpy( &samples1[shdr->numSamples-vars->overlap_numSamples], samples2, vars->overlap_numSamples*sizeof(float) );
    }
    // 2) Move samples in trace 2 down by # overlap samples
    for( int isamp = vars->numSamplesIn-1; isamp >= 0; isamp-- ) {
      samples2[isamp+vars->overlap_numSamples] = samples2[isamp];
    }
    // 3) Copy last (original!) # overlap samples from trace 1 to top of trace 2
    if( doCopy ) {
      memcpy( samples2, &samples1[shdr->numSamples-2*vars->overlap_numSamples], vars->overlap_numSamples*sizeof(float) );
    }
    else { // Zero out top of trace #2
      for( int isamp = 0; isamp < vars->overlap_numSamples; isamp++ ) {
        samples2[isamp] = 0.0f;
      }
    }
    *numTrcToKeep = 1; // Last trace is kept. First trace is output
  }
  else {
    for( int isamp = shdr->numSamples-vars->overlap_numSamples; isamp < shdr->numSamples; isamp++ ) {
      samples1[isamp] = 0.0f;
    }
    *numTrcToKeep = 0;
  }

  if( edef->isDebug() ) {
    for( int itrc = 0; itrc < traceGather->numTraces(); itrc++ ) {
      csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
      int time_samp1_s  = trcHdr->intValue( vars->hdrID_time_samp1_s );
      int time_samp1_us = trcHdr->intValue( vars->hdrID_time_samp1_us );
      writer->line("Trace #-5d,  Time:  %12ds   %10dus", itrc, time_samp1_s, time_samp1_us);
    }
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_overlap_( csParamDef* pdef ) {
  pdef->setModule( "OVERLAP", "Create data overlap between adjacent traces", "Duplicates data from adjacent traces and pads it at start and end of trace" );

  pdef->addParam( "overlap", "Size of overlap [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Size of overlap [ms]", "Data of length 'overlap'ms from adjacent traces is added to start and end of trace" );

  pdef->addParam( "absolute_time", "Acknowledge absolute time?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Absolute time is acknowledged. No overlap will be created from adjacent traces that do not have adjacent absolute time stamps" );
  pdef->addOption( "no", "Absolute time is not acknowledged. Overlap will be created for all adjacent traces" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_overlap_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_overlap::VariableStruct* vars = reinterpret_cast<mod_overlap::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_overlap_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_overlap::VariableStruct* vars = reinterpret_cast<mod_overlap::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_overlap_( csParamDef* pdef ) {
  params_mod_overlap_( pdef );
}
extern "C" void _init_mod_overlap_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_overlap_( param, env, writer );
}
extern "C" bool _start_exec_mod_overlap_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_overlap_( env, writer );
}
extern "C" void _exec_mod_overlap_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_overlap_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_overlap_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_overlap_( env, writer );
}

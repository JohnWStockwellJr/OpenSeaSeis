/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csGeolibUtils.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: RMS
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_rms {
  struct VariableStruct {
    int startSamp;
    int endSamp;
    int hdrId_rms;
    int hdrType_rms;
    int mode;
    int hdrId_trcno;

    int dataWinStepInSamples;
    int dataWinLengthInSamples;
    int dataWinHalfLengthInSamples;
    int dataLastMidSample;
    float* dataBuffer;
  };
  static int const MODE_HEADER = 1;
  static int const MODE_DATA   = 2;
  static int const MODE_ENSEMBLE = 3;
}

using mod_rms::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_rms_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->startSamp   = 0;
  vars->endSamp     = 0;
  vars->hdrId_rms   = -1;
  vars->hdrType_rms = 0;
  vars->mode = mod_rms::MODE_HEADER;

  vars->dataWinStepInSamples   = 0;
  vars->dataWinLengthInSamples = 0;
  vars->dataWinHalfLengthInSamples = 0;
  vars->dataLastMidSample      = 0;
  vars->dataBuffer = NULL;

  //---------------------------------------------------------

  if( param->exists("mode") ) {
    std::string text;
    param->getString( "mode", &text );
    if( !text.compare( "header" ) ) {
      vars->mode = mod_rms::MODE_HEADER;
    }
    else if( !text.compare( "data" ) ) {
      vars->mode = mod_rms::MODE_DATA;
    }
    else if( !text.compare( "ensemble" ) ) {
      vars->mode = mod_rms::MODE_ENSEMBLE;
    }
    else {
      writer->error("Option not recognized: '%s'.", text.c_str());
    }
  }


  bool isTimeDomain = true;

  if( param->exists("domain") ) {
    std::string text;
    param->getString( "domain", &text );
    if( !text.compare( "sample" ) ) {
      isTimeDomain = false;
    }
    else if( !text.compare( "time" ) ) {
      isTimeDomain = true;
    }
    else {
      writer->line("Domain option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  if( vars->mode == mod_rms::MODE_HEADER || vars->mode == mod_rms::MODE_ENSEMBLE ) {
    float startTime = 0.0f;
    float endTime   = 0.0f;
    double start_in = 0.0f;
    double end_in   = 0.0f;
    if( param->exists("start") ) {
      param->getDouble( "start", &start_in );
    }
    if( param->exists("end") ) {
      param->getDouble( "end", &end_in );
    }
    if( isTimeDomain ) {
      startTime = start_in;
      endTime   = end_in;
      if( startTime < 0 ) {
        writer->warning("Start time (%f) needs to be greater or equal to 0.0.", startTime);
        startTime = 0;
      }
      if( endTime > (shdr->numSamples-1)*shdr->sampleInt ) {
        writer->warning("End time (%f) exceeds length of trace (%f).", endTime, (shdr->numSamples-1)*shdr->sampleInt );
        endTime = (shdr->numSamples-1)*shdr->sampleInt;
      }
      else if( endTime == 0 ) {
        endTime = (shdr->numSamples-1)*shdr->sampleInt;
      }
      if( startTime > endTime ) writer->error("Start time (%f) needs to be smaller than end time (%f).", startTime, endTime);
      
      vars->startSamp = (int)(startTime / shdr->sampleInt);  // All in milliseconds
      vars->endSamp   = (int)(endTime / shdr->sampleInt);
    }
    else {
      //
      // NOTE: User input is '1' for first sample. Internally, '0' is used!!
      //
      vars->startSamp = (int)start_in;
      vars->endSamp   = (int)end_in;
      if( vars->startSamp < 1 ) writer->error("Start sample (%d) needs to be greater or equal to 1.", vars->startSamp);
      if( vars->startSamp > vars->endSamp ) writer->error("Start sample (%d) needs to be smaller than end sample (%d).", vars->startSamp, vars->endSamp);
      if( vars->endSamp > shdr->numSamples ) writer->error("End sample (%d) exceeds number of samples (%d).", vars->endSamp, shdr->numSamples );
      else if( vars->endSamp == 0 ) {
        vars->endSamp = shdr->numSamples;
      }
      vars->startSamp -= 1;   // see note above..
      vars->endSamp   -= 1;
      startTime = (float)vars->startSamp * shdr->sampleInt;
      endTime   = (float)vars->endSamp * shdr->sampleInt;
    }

    //---------------------------------------------
    //
    std::string headerName_rms("rms");
    if( param->exists("hdr_rms") ) {
      param->getString("hdr_rms", &headerName_rms, 0);
    }
    if( !hdef->headerExists( headerName_rms ) ) {
      hdef->addHeader( TYPE_FLOAT, headerName_rms, "RMS value" );
    }

    vars->hdrId_rms = hdef->headerIndex(headerName_rms);
    vars->hdrId_trcno = hdef->headerIndex("trcno");
    vars->hdrType_rms = hdef->headerType(headerName_rms);
    if( vars->hdrType_rms != TYPE_FLOAT && vars->hdrType_rms != TYPE_DOUBLE ) {
      writer->error("Trace header '%s' exists but has wrong type (%s). Type should be 'float' or 'double'.",
                    cseis_geolib::csGeolibUtils::typeText(vars->hdrType_rms), headerName_rms.c_str());
    }
  } // END: if mode == mode_header
  else if( vars->mode == mod_rms::MODE_DATA ) {
    float length = 0;
    float step   = 0;
    param->getFloat( "data_param", &length, 0 );
    param->getFloat( "data_param", &step, 1 );
    if( isTimeDomain ) {
      vars->dataWinLengthInSamples = (int)round( length / shdr->sampleInt );
      vars->dataWinStepInSamples   = (int)round( step / shdr->sampleInt );
    }
    else {
      vars->dataWinLengthInSamples = (int)round( length );
      vars->dataWinStepInSamples   = (int)round( step );
    }
    vars->dataWinStepInSamples   = std::max( vars->dataWinStepInSamples, 1 );
    vars->dataWinLengthInSamples = std::min( std::max( vars->dataWinLengthInSamples, 2 ), shdr->numSamples-1 );
    if( (int)(vars->dataWinLengthInSamples/2)*2 == vars->dataWinLengthInSamples ) vars->dataWinLengthInSamples += 1; // Window must be odd number of samples
    vars->dataWinHalfLengthInSamples = vars->dataWinLengthInSamples/2;
    vars->dataLastMidSample = shdr->numSamples - vars->dataWinHalfLengthInSamples - 1;
    int numSteps = ( vars->dataLastMidSample - vars->dataWinHalfLengthInSamples ) / vars->dataWinStepInSamples;
    vars->dataLastMidSample = vars->dataWinHalfLengthInSamples + numSteps * vars->dataWinStepInSamples;
    vars->dataBuffer = new float[shdr->numSamples];
  }

  if( vars->mode == mod_rms::MODE_ENSEMBLE ) {
    int numSamplesOut = 0;
    param->getInt( "nsamples_out", &numSamplesOut );
    shdr->numSamples = numSamplesOut;
    shdr->sampleInt  = 1.0f;
    env->execPhaseDef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }
  else {
    env->execPhaseDef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_rms_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);

  if( vars->mode == mod_rms::MODE_HEADER ) {
    int nSamples = vars->endSamp - vars->startSamp + 1;
    float rms = compute_rms( &(trace->getTraceSamples()[vars->startSamp]), nSamples );
    // Set double value to account for different header types
    trace->getTraceHeader()->setDoubleValue( vars->hdrId_rms, rms );
  }
  else if( vars->mode == mod_rms::MODE_DATA ) {
    float* samplesIn = trace->getTraceSamples();
    for( int isamp = vars->dataWinHalfLengthInSamples; isamp <= vars->dataLastMidSample; isamp += vars->dataWinStepInSamples ) {
      vars->dataBuffer[isamp] = compute_rms( &(samplesIn[isamp-vars->dataWinHalfLengthInSamples]), vars->dataWinLengthInSamples );
    }
    float rms = vars->dataBuffer[vars->dataWinHalfLengthInSamples];
    for( int isamp = 0; isamp < vars->dataWinHalfLengthInSamples; isamp++ ) {
      vars->dataBuffer[isamp] = rms;
    }
    //--------- MID -------------
    float rms1 = vars->dataBuffer[vars->dataWinHalfLengthInSamples];
    for( int isamp = vars->dataWinHalfLengthInSamples; isamp <= vars->dataLastMidSample; isamp += vars->dataWinStepInSamples ) {
      float rms2 = vars->dataBuffer[isamp+vars->dataWinStepInSamples];
      float rms1Ratio = rms1 / (float)vars->dataWinStepInSamples;
      float rms2Ratio = rms2 / (float)vars->dataWinStepInSamples;
      for( int counter = 1; counter < vars->dataWinStepInSamples; counter++ ) {
        vars->dataBuffer[isamp+counter] = rms1 + (float)counter*( rms2Ratio - rms1Ratio );
      }
      rms1 = rms2;
    }
    //------------- LAST ----------------
    rms = vars->dataBuffer[vars->dataLastMidSample];
    for( int outsamp = vars->dataLastMidSample+1; outsamp < shdr->numSamples; outsamp++ ) {
      vars->dataBuffer[outsamp] = rms;
    }
    memcpy( samplesIn, vars->dataBuffer, shdr->numSamples*sizeof(float) );
  }
  else if( vars->mode == mod_rms::MODE_ENSEMBLE ) {
    float* samplesOut = trace->getTraceSamples();
    int numTracesIn = traceGather->numTraces();
    int numTracesToProcess = std::min( shdr->numSamples, numTracesIn );
    int nSamples = vars->endSamp - vars->startSamp + 1;
    for( int itrc = 0; itrc < numTracesToProcess; itrc++ ) {
      float* samplesCurrent = traceGather->trace( itrc )->getTraceSamples();
      float rms = compute_rms( &samplesCurrent[vars->startSamp], nSamples );
      samplesOut[itrc] = rms;
    }
    for( int isamp = numTracesToProcess; isamp < shdr->numSamples; isamp++ ) {
      samplesOut[isamp] = 0.0f;
    }
    traceGather->freeTraces( 1, numTracesIn-1 );
  }
  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_rms_( csParamDef* pdef ) {
  pdef->setModule( "RMS", "Compute RMS value in given time window", "Computed RMS value is stored in trace header 'rms'" );

  pdef->addParam( "mode", "RMS mode", NUM_VALUES_FIXED );
  pdef->addValue( "header", VALTYPE_OPTION );
  pdef->addOption( "header", "Compute RMS over time window, store value in trace header" );
  pdef->addOption( "data", "Compute RMS over sliding window along trace, store in trace samples", "Specify user parameter 'data_param'" );
  pdef->addOption( "ensemble", "Compute RMS over all traces of input gather, output one RMS trace", "Specify user parameter 'nsamples_out'" );

  pdef->addParam( "domain", "Time or sample domain", NUM_VALUES_FIXED );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Window is specified in time [ms] (or frequency [Hz])" );
  pdef->addOption( "sample", "Window is specified in samples (1 for first sample)" );

  pdef->addParam( "start", "Start time/sample", NUM_VALUES_FIXED, "Start time or sample, this depends on the 'domain' setting" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Start time/sample of cross-correlation window" );

  pdef->addParam( "end", "End time/sample", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "End time/sample of cross-correlation window", "= 0 : Until end of trace" );

  pdef->addParam( "hdr_rms", "Name of trace header where RMS value shall be stored", NUM_VALUES_FIXED, "...may be new or existing trace header" );
  pdef->addValue( "rms", VALTYPE_STRING );

  pdef->addParam( "data_param", "Parameters for RMS 'data' mode: RMS window length and increment", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Window length [ms]" );
  pdef->addValue( "0", VALTYPE_STRING, "Time step [ms]", "0: Single sample step. Interpolate linearly between samples for step > sample interval" );

  pdef->addParam( "nsamples_out", "Number of samples in output RMS trace. Only applies to mode 'ensemble'", NUM_VALUES_FIXED, "Should match number of traces per input ensemble" );
  pdef->addValue( "", VALTYPE_NUMBER );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_rms_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_rms::VariableStruct* vars = reinterpret_cast<mod_rms::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;

  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_rms_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_rms::VariableStruct* vars = reinterpret_cast<mod_rms::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->dataBuffer != NULL ) {
    delete vars->dataBuffer;
    vars->dataBuffer = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_rms_( csParamDef* pdef ) {
  params_mod_rms_( pdef );
}
extern "C" void _init_mod_rms_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_rms_( param, env, writer );
}
extern "C" bool _start_exec_mod_rms_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_rms_( env, writer );
}
extern "C" void _exec_mod_rms_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_rms_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_rms_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_rms_( env, writer );
}

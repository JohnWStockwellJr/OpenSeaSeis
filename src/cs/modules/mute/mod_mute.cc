/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csTableAll.h"
#include "csTableManagerNew.h"
#include "csFlexNumber.h"
#include "csStandardHeaders.h"
#include "csSort.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: MUTE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_mute {
  struct VariableStruct {
    int mode;
    float mute_time;
    int mute_start_samp;
    int mute_end_samp;
    csTableManagerNew* tableManager;
    bool killZeroTraces;
    int taperLengthSamples;
    bool indicate;
    float indicateValue;
    int   indicateWidthSamples;
    int taperType;
    int hdrId_time;
    int hdrId_taper;
    int hdrId_mute;
    int windowStartSample;
    int windowEndSample;
    bool isMedian;
    float* sortValues;
    cseis_geolib::csSort<float>* sortObj;
    bool isBackwardCompatible;
  };
  static int const TAPER_LINEAR = 1;
  static int const TAPER_COSINE = 2;
}
using mod_mute::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_mute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->mode         = MUTE_FRONT;
  vars->tableManager        = NULL;
  vars->mute_start_samp = 0;
  vars->mute_end_samp   = 0;
  vars->killZeroTraces  = false;
  vars->indicate        = false;
  vars->indicateValue   = 0.0;
  vars->indicateWidthSamples = 1;
  vars->taperType       = mod_mute::TAPER_COSINE;
  float taperLength     = 0; // See below
  vars->hdrId_time  = -1;
  vars->hdrId_mute  = -1;
  vars->hdrId_taper = -1;
  vars->windowStartSample = 0; // By default, apply mute to full trace length
  vars->windowEndSample   = shdr->numSamples-1; // By default, apply mute to full trace length
  vars->isMedian = false;
  vars->sortValues = NULL;
  vars->sortObj = NULL;
  vars->isBackwardCompatible = false;

//---------------------------------------------
// Retrieve mute table
//
  std::string text;

  if( param->exists("table") ) {

    param->getString("table", &text );

    try {
      vars->tableManager = new csTableManagerNew( text, csTableAll::TABLE_TYPE_UNIQUE_KEYS, hdef );
      if( vars->tableManager->valueName().compare("time") ) {
        writer->error("Mute table must contain 'value' column labelled 'time', for example: '@%s time'. Value label found: '%s'",
                   vars->tableManager->numKeys() > 0 ? vars->tableManager->keyName(0).c_str() : "offset",
                   vars->tableManager->valueName().c_str() );
      }
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s': %s\n", text.c_str(), exc.getMessage() );
    }
    if( vars->tableManager->numKeys() < 1 ) {
      writer->error("Number of table keys = %d. Specify table key by placing the character '%c' in front of the key name. Example: %csource time vel  (source is a table key)",
                 csTableAll::KEY_CHAR, csTableAll::KEY_CHAR);
    }

    if( edef->isDebug() ) {
      vars->tableManager->dump();
    }

  }

//---------------------------------------------
// Retrieve mute time
//
  if( param->exists("mode") ) {
    param->getString("mode",&text);
    if( !text.compare("front") ) {
      vars->mode = MUTE_FRONT;
    }
    else if( !text.compare("end") ) {
      vars->mode = MUTE_END;
    }
    else {
      writer->error("Unknown mode option: '%s'", text.c_str());
    }
  }
  if( param->exists("window") ) {
    float startTime;
    float endTime;
    param->getFloat( "window", &startTime, 0 );
    param->getFloat( "window", &endTime, 1 );
    vars->windowStartSample = std::max( 0, (int)round( startTime / (float)shdr->sampleInt ) );
    vars->windowEndSample   = std::min( shdr->numSamples-1, (int)round( endTime / (float)shdr->sampleInt ) );
    if( vars->windowStartSample >= shdr->numSamples ) writer->error("Window start time (=%f) exceeds trace length (=%f)", startTime, shdr->numSamples );
    if( vars->windowEndSample <= vars->windowStartSample ) writer->error("Window start time (=%f) exceeds or equals window end time (=%f)", startTime, endTime );
  }

  if( param->exists("kill") ) {
    param->getString("kill",&text);
    if( !text.compare("yes") ) {
      vars->killZeroTraces = true;
    }
    else if( !text.compare("no") ) {
      vars->killZeroTraces = false;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  if( param->exists("indicate") ) {
    vars->indicate = true;
    param->getFloat( "indicate", &vars->indicateValue );
    if( param->getNumValues("indicate") > 1 ) {
      param->getInt( "indicate", &vars->indicateWidthSamples, 1 );
    }
  }

  if( param->exists("taper_len") ) {
    param->getString("taper_len",&text);
    csFlexNumber number;
    if( !number.convertToNumber( text ) ) {
      vars->hdrId_taper = hdef->headerIndex(text);
    }
    else {
      taperLength = number.floatValue();
    }
  }
  vars->taperLengthSamples = (int)( taperLength/shdr->sampleInt + 0.5 );

  if( param->exists("taper_type") ) {
    param->getString("taper_type", &text);
    if( !text.compare("cos") ) {
      vars->taperType = mod_mute::TAPER_COSINE;
    }
    else if( !text.compare("linear") ) {
      vars->taperType = mod_mute::TAPER_LINEAR;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  if( param->exists("back_comp") ) {
    param->getString("back_comp", &text);
    if( !text.compare("yes") ) {
      vars->isBackwardCompatible = true;
    }
    else if( !text.compare("no") ) {
      vars->isBackwardCompatible = false;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  if( !vars->tableManager ) {
    param->getString("time",&text);
    csFlexNumber number;
    if( !number.convertToNumber( text ) ) {
      vars->hdrId_time = hdef->headerIndex(text);
    }
    else {
      vars->mute_time = number.floatValue();
      float recordLength = (float)shdr->numSamples * shdr->sampleInt;
      int sampleIndex = (int)( vars->mute_time / shdr->sampleInt );

      if( sampleIndex < 0 ) {
        writer->warning("Specified mute time (%.1fms) is negative. Set to 0ms.", vars->mute_time );
        vars->mute_time = 0;
        sampleIndex = 0;
      }
      else if( sampleIndex >= shdr->numSamples ) {
        writer->warning("Specified mute time (%.1fms) is outside of valid range (0-%.1fms). Set to %.1fms.",
                     vars->mute_time, recordLength, recordLength );
        vars->mute_time = recordLength;
        sampleIndex = shdr->numSamples-1;
      }
      if( vars->mode == MUTE_FRONT ) {
        vars->mute_start_samp = 0;
        vars->mute_end_samp   = sampleIndex;
        if( vars->isBackwardCompatible ) vars->mute_end_samp -= 1;
      }
      else {
        vars->mute_start_samp = sampleIndex;
        vars->mute_end_samp   = shdr->numSamples-1;
        if( vars->isBackwardCompatible ) vars->mute_start_samp += 1;
      }
    }
  }

  if( param->exists("mute_value") ) {
    param->getString("mute_value", &text );
    if( !text.compare("zero") ) {
      vars->isMedian = false;
    }
    else if( !text.compare("median") ) {
      vars->isMedian = true;
      vars->sortObj = new csSort<float>();
      vars->sortValues = new float[shdr->numSamples];
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }


//----------------------------------------------
  csHeaderInfo info;
  if( vars->mode == MUTE_FRONT ) {
    info = cseis_geolib::HDR_MUTE_START;
  }
  else {
    info = cseis_geolib::HDR_MUTE_END;
  }
  if( !hdef->headerExists(info.name) ) {
    hdef->addStandardHeader(info.name);
  }
  vars->hdrId_mute = hdef->headerIndex(info.name);
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_mute_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;

  csTrace* trace = traceGather->trace(0);


  float* samples = trace->getTraceSamples();
  csTraceHeader* trcHdr = trace->getTraceHeader();

  if( vars->tableManager != NULL ) {
    vars->mute_time = vars->tableManager->getValue( trace->getTraceHeader() );
    int sampleIndex = (int)( vars->mute_time / shdr->sampleInt );
    if( vars->mode == MUTE_FRONT ) {
      vars->mute_start_samp = 0;
      vars->mute_end_samp   = std::min( shdr->numSamples-1,sampleIndex);
      if( vars->isBackwardCompatible ) vars->mute_end_samp = std::min( shdr->numSamples-1,sampleIndex-1);
    }
    else { // MUTE_END
      vars->mute_start_samp = std::max( 0, sampleIndex );
      vars->mute_end_samp   = shdr->numSamples-1;
      if( vars->isBackwardCompatible ) vars->mute_start_samp = std::max( 0, sampleIndex+1 );
    }
    if( edef->isDebug() ) writer->line("Mute time: %f  Start/end sample: %d %d", vars->mute_time, vars->mute_start_samp, vars->mute_end_samp );
  }
  else {
    // Read in values from trace headers if necessary
    if( vars->hdrId_time >= 0 ) {
      vars->mute_time = trcHdr->floatValue( vars->hdrId_time );
      float recordLength = (float)shdr->numSamples * shdr->sampleInt;
      if( vars->mute_time < 0 ) {
        vars->mute_time = 0;
      }
      else if( vars->mute_time > recordLength ) {
        vars->mute_time = recordLength;
      }
      int sampleIndex = (int)( vars->mute_time / shdr->sampleInt );
      if( vars->mode == MUTE_FRONT ) {
        vars->mute_start_samp = 0;
        vars->mute_end_samp   = sampleIndex;
        if( vars->isBackwardCompatible ) vars->mute_end_samp -= 1;
      }
      else {
        vars->mute_start_samp = sampleIndex;
        vars->mute_end_samp   = shdr->numSamples-1;
        if( vars->isBackwardCompatible ) vars->mute_start_samp += 1;
      }
    }
  }
  if( vars->hdrId_taper >= 0 ) {
    float taperLength = trcHdr->floatValue( vars->hdrId_taper );
    vars->taperLengthSamples = (int)( taperLength/shdr->sampleInt + 0.5 );
    if( vars->taperLengthSamples < 0 ) vars->taperLengthSamples = 0;
    if( edef->isDebug() ) fprintf(stdout,"Taper %f\n", taperLength);
  }

  if( edef->isDebug() ) {
    writer->line("Mute time sample index: %d %d (nsamp: %d)", vars->mute_start_samp, vars->mute_end_samp, shdr->numSamples);
  }

  if( vars->indicate ) {
    int sampleIndex = (int)( vars->mute_time / shdr->sampleInt );
    int startSamp = std::max( vars->windowStartSample, std::max( sampleIndex-vars->indicateWidthSamples, 0) );
    int endSamp   = std::min( vars->windowEndSample+1, std::min( sampleIndex+vars->indicateWidthSamples, shdr->numSamples ) );
    for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
      samples[isamp] = vars->indicateValue;
    }
    if( vars->mode == MUTE_FRONT ) {
      trcHdr->setFloatValue( vars->hdrId_mute, startSamp * shdr->sampleInt );
    }
    else {
      trcHdr->setFloatValue( vars->hdrId_mute, endSamp * shdr->sampleInt );
    }
    return;
  }
  else {
    if( vars->killZeroTraces && vars->mute_end_samp >= shdr->numSamples-1 && vars->mute_start_samp <= 1 ) {
      traceGather->freeAllTraces();
      return;
    }
    int startSamp = std::max( vars->windowStartSample, vars->mute_start_samp);
    int endSamp   = std::min( vars->windowEndSample, vars->mute_end_samp );
    float muteValue = 0.0f;
    for( int isamp = startSamp; isamp <= endSamp; isamp++ ) {
      samples[isamp] = 0.0;
    }
    if( vars->mode == MUTE_FRONT ) {
      int startSampTaper = vars->mute_end_samp+1;
      startSamp = std::max( std::max( vars->windowStartSample, startSampTaper ), 0 );
      endSamp   = std::min( vars->windowEndSample+1, std::min( vars->mute_end_samp+1+vars->taperLengthSamples, shdr->numSamples ) );
      if( edef->isDebug() ) writer->line("Start/end sample of taper: %d %d", startSamp, endSamp);

      if( vars->isMedian ) {
        int ns = endSamp-startSamp+1;
        memcpy( vars->sortValues, &samples[startSamp], ns*sizeof(float) );
        vars->sortObj->treeSort( vars->sortValues, ns );
        muteValue = vars->sortValues[ns/2];
        for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
          samples[isamp] -= muteValue;
        }
      }

      trcHdr->setFloatValue( vars->hdrId_mute, startSamp * shdr->sampleInt );
      if( vars->taperType == mod_mute::TAPER_LINEAR ) {
        for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
          float weight = (float)(isamp-startSampTaper+1)/(float)vars->taperLengthSamples;
          samples[isamp] *= weight;
        }
      }
      else { // COSINE taper
        for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
          float phase  = ( ( (float)(isamp-startSampTaper) / (float)vars->taperLengthSamples ) - 1.0 ) * M_PI;
          float weight = 0.5 * (cos(phase) + 1.0);
          samples[isamp] *= weight;
        }
      }
    }
    else {
      startSamp = std::max( vars->windowStartSample, std::max( 0, vars->mute_start_samp-vars->taperLengthSamples ) );
      int endSampTaper = std::min( vars->mute_start_samp, shdr->numSamples );
      endSamp   = std::min( vars->windowEndSample+1, endSampTaper );

      if( vars->isMedian ) {
        int ns = endSamp-startSamp+1;
        memcpy( vars->sortValues, &samples[startSamp], ns*sizeof(float) );
        vars->sortObj->treeSort( vars->sortValues, ns );
        muteValue = vars->sortValues[ns/2];
        for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
          samples[isamp] -= muteValue;
        }
      }

      trcHdr->setFloatValue( vars->hdrId_mute, endSamp * shdr->sampleInt );
      if( edef->isDebug() ) writer->line("Start/end sample of taper: %d %d", startSamp, endSamp);
      if( vars->taperType == mod_mute::TAPER_LINEAR ) {
        for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
          float weight = (float)(endSampTaper-isamp)/(float)vars->taperLengthSamples;
          samples[isamp] *= weight;
        }
      }
      else { // COSINE taper
        for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
          float phase  = ( ( (float)(endSampTaper-isamp) / (float)vars->taperLengthSamples ) - 1.0 ) * M_PI;
          float weight = 0.5 * (cos(phase) + 1.0);
          samples[isamp] *= weight;
        }
      }
    } // if END mute

    if( vars->isMedian ) {
      for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
        samples[isamp] += muteValue;
      }
      if( vars->mode == MUTE_FRONT ) {
        for( int isamp = 0; isamp < startSamp; isamp++ ) {
          samples[isamp] = muteValue;
        }
      }
      else {
        for( int isamp = endSamp; isamp < shdr->numSamples; isamp++ ) {
          samples[isamp] = muteValue;
        }
      }
    } // END if isMedian
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_mute_( csParamDef* pdef ) {
  pdef->setModule( "MUTE", "Mute trace data" );

  pdef->addParam( "mode", "Front mute or end mute", NUM_VALUES_FIXED );
  pdef->addValue( "front", VALTYPE_OPTION );
  pdef->addOption( "front", "Front mute. Specify time of first unmuted sample.", "Samples are muted from 0ms to specified time" );
  pdef->addOption( "end", "End mute. Specify time of last unmuted sample.", "Samples are muted from specified time to the end of the trace" );

  pdef->addParam( "table", "Mute table", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Mute table file, full path.",
                  "The mute table must have at least two columns, one giving a key and the second giving the mute time in [ms]" );
  //  pdef->addValue( "", VALTYPE_STRING, "Mute table. This must be specified within the flow file, using directive &table <tablename> <filename> {methodInterpolation}",
  //    "The mute table must have at least two columns, one giving a key and the second giving the mute time in [ms]" );

  pdef->addParam( "time", "Mute time [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "Mute time [ms] (or header name containing mute time)" );

  pdef->addParam( "taper_len", "Taper length [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_HEADER_NUMBER, "Mute taper length [ms] (or header name containing taper length)" );

  pdef->addParam( "taper_type", "Type of mute taper", NUM_VALUES_FIXED );
  pdef->addValue( "cos", VALTYPE_OPTION );
  pdef->addOption( "linear", "Apply linear taper" );
  pdef->addOption( "cos", "Apply cosine taper" );

  pdef->addParam( "window", "Restrict mute to certain window?", NUM_VALUES_FIXED, "If specified, mute will not be applied outside of the specified window, including no tapering" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Start time" );
  pdef->addValue( "99999", VALTYPE_NUMBER, "End time" );

  pdef->addParam( "kill", "Kill zero traces?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Kill zero traces" );
  pdef->addOption( "yes", "Do not kill zero traces" );

  pdef->addParam( "indicate", "Do not mute, indicate mute time by setting samples to the given value", NUM_VALUES_VARIABLE,
                  "Input data will not be muted. Instead, spikes are placed at the mute times, with the given amplitude." );
  pdef->addValue( "0", VALTYPE_NUMBER, "Value that mute samples are set to." );
  pdef->addValue( "1", VALTYPE_NUMBER, "Width in samples to indicate with given value." );

  pdef->addParam( "mute_value", "Mute value", NUM_VALUES_FIXED );
  pdef->addValue( "zero", VALTYPE_OPTION );
  pdef->addOption( "zero", "Set fully muted values to zero" );
  pdef->addOption( "median", "Use median value determined over taper zone" );

  pdef->addParam( "back_comp", "Backward compatibility?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Use latest implementation" );
  pdef->addOption( "yes", "Use 'previous' implementation" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_mute::VariableStruct* vars = reinterpret_cast<mod_mute::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_mute::VariableStruct* vars = reinterpret_cast<mod_mute::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->tableManager != NULL ) {
    delete vars->tableManager;
    vars->tableManager = NULL;
  }
  if( vars->sortObj != NULL ) {
    delete vars->sortObj;
    vars->sortObj = NULL;
  }
  if( vars->sortValues != NULL ) {
    delete [] vars->sortValues;
    vars->sortValues = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_mute_( csParamDef* pdef ) {
  params_mod_mute_( pdef );
}
extern "C" void _init_mod_mute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_mute_( param, env, writer );
}
extern "C" bool _start_exec_mod_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_mute_( env, writer );
}
extern "C" void _exec_mod_mute_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_mute_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_mute_( env, writer );
}

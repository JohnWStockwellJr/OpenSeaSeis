/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csInterpolation.h"
#include <string>
#include <cstring>
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: STATICS
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_statics {
  struct VariableStruct {
    bool   doHeaderStatic;
    int    hdrId;
    int    hdrType;
    float  bulkShift_ms;
    float* buffer;
    bool   isApplyMode;
    int hdrID_time_samp1_s;
    int hdrID_time_samp1_us;
    int hdrID_stat_res;
    cseis_geolib::csInterpolation* interpol;
    bool  applyFraction;
    bool isSampleDomain;
    bool updateTimeHdr;
  };
}
using mod_statics::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_statics_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->doHeaderStatic = false;
  vars->hdrId     = -1;
  vars->hdrType   = 0;
  vars->bulkShift_ms = 0;
  vars->buffer    = NULL;
  vars->isApplyMode = true;
  vars->interpol  = NULL;
  vars->hdrID_time_samp1_s = -1;
  vars->hdrID_time_samp1_us = -1;
  vars->hdrID_stat_res     = -1;
  vars->applyFraction = true;
  vars->isSampleDomain = false;
  vars->updateTimeHdr = false;

//---------------------------------------------

  string text;
  if( param->exists("domain") ) {
    param->getString( "domain", &text );
    if( !text.compare( "sample" ) ) {
      vars->isSampleDomain = true;
    }
    else if( !text.compare( "time" ) ) {
      vars->isSampleDomain = false;
    }
    else {
      writer->line("Domain option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  if( param->exists("update_time_hdr") ) {
    param->getString( "update_time_hdr", &text );
    if( !text.compare( "yes" ) ) {
      vars->updateTimeHdr = true;
    }
    else if( !text.compare( "time" ) ) {
      vars->updateTimeHdr = false;
    }
    else {
      writer->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  // For now, stick to the default interpolation re number of coefficients etc
  int numCoefficients = 8;
  if( param->exists("ncoef") ) {
    param->getInt("ncoef", &numCoefficients );
    if( numCoefficients < 2 ) {
      writer->error("Number of interpolation coefficients too small: %d. Minimum is 2", numCoefficients);
    }
    if( numCoefficients > 32 ) {
      writer->warning("Number of interpolation coefficients is very large: %d.", numCoefficients);
    }
    if( numCoefficients > 128 ) {
      writer->warning("Number of interpolation coefficients too large: %d.", numCoefficients);
    }
  }

  if( param->exists("header") ) {
    std::string headerName;
    vars->doHeaderStatic = true;
    param->getString( "header", &headerName );
    if( ( headerName.c_str() ) ) {
      vars->hdrType = hdef->headerType( headerName.c_str() );
      vars->hdrId   = hdef->headerIndex( headerName.c_str() );
    }
    else {
      writer->warning("Unknown trace header name: '%s'", headerName.c_str());
      env->addError();
    }
    writer->line("User specified header: '%s', type#: %d", headerName.c_str(), vars->hdrType );
  }

  if( param->exists("mode") ) {
    param->getString("mode",&text);
    if( !text.compare("apply") ) {
      vars->isApplyMode = true;
    }
    else if( !text.compare("remove") ) {
      vars->isApplyMode = false;
    }
    else {
      writer->error("Unknown mode option: '%s'", text.c_str());
    }
  }

  if( param->exists("apply_fraction") ) {
    string text;
    param->getString("apply_fraction",&text);
    if( !text.compare("yes") ) {
      vars->applyFraction = true;
    }
    else if( !text.compare("no") ) {
      vars->applyFraction = false;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }

  if( param->exists("bulk_shift") ) {
    param->getFloat( "bulk_shift", &vars->bulkShift_ms );  // [ms]
    if( vars->isSampleDomain ) vars->bulkShift_ms *= shdr->sampleInt;
  }
  else {
    if( !vars->doHeaderStatic ) {
      writer->line("Error: No statics option selected. Empty module.");
      env->addError();
    }
    vars->bulkShift_ms = 0.0;
  }

  if( !vars->isApplyMode ) {
    vars->bulkShift_ms *= -1.0;
  }

  if( vars->applyFraction ) {
    vars->interpol = new csInterpolation( shdr->numSamples, shdr->sampleInt, numCoefficients );
  }
  if( !hdef->headerExists( HDR_STAT_RES.name ) ) {
    hdef->addStandardHeader( HDR_STAT_RES.name );
  }
  vars->hdrID_stat_res = hdef->headerIndex( HDR_STAT_RES.name );
  
  vars->buffer = new float[shdr->numSamples];
  if( vars->updateTimeHdr ) {
    vars->hdrID_time_samp1_s  = hdef->headerIndex( HDR_TIME_SAMP1.name );
    vars->hdrID_time_samp1_us = hdef->headerIndex( HDR_TIME_SAMP1_US.name );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_statics_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  float shift_ms = vars->bulkShift_ms;  // Shift in [ms]
  if( vars->doHeaderStatic ) {
    float stat_hdr = trace->getTraceHeader()->floatValue(vars->hdrId);
    if( vars->isSampleDomain ) stat_hdr *= shdr->sampleInt;
    if( vars->isApplyMode ) {
      shift_ms += stat_hdr;
    }
    else {
      shift_ms -= stat_hdr;
    }
    if( edef->isDebug() ) writer->line("Static total shift = %f ms (...incl header shift: %f ms)", shift_ms, stat_hdr);
  }
  float* samples = trace->getTraceSamples();

  if( edef->isDebug() ) { writer->line("Apply static %f ms", shift_ms); }

  csTraceHeader* trcHdr = trace->getTraceHeader();
  bool isFullSample =  ( (float)( (int)round( fabs(shift_ms) / shdr->sampleInt ) ) * shdr->sampleInt == fabs(shift_ms) );

  if( vars->applyFraction && !isFullSample ) {
    memcpy( vars->buffer, trace->getTraceSamples(), shdr->numSamples*sizeof(float) );
    vars->interpol->static_shift( shift_ms, vars->buffer, samples );
    trcHdr->setFloatValue( vars->hdrID_stat_res, 0 );
  }
  else {
    float* samples = trace->getTraceSamples();
    int shiftSamples = (int)cseis_geolib::SIGN(shift_ms)*(int)( fabs(shift_ms) / shdr->sampleInt );
    shift_ms = shdr->sampleInt * shiftSamples;
    float stat_res = shift_ms - shift_ms;
    trcHdr->setFloatValue( vars->hdrID_stat_res, stat_res );
    if( shiftSamples < 0 ) {
      for( int isamp = -shiftSamples; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp+shiftSamples] = samples[isamp];
      }
      for( int isamp = std::max(0,shdr->numSamples+shiftSamples); isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] = 0;
      }
    }
    else if( shiftSamples > 0 ) {
      shiftSamples = std::min( shiftSamples, shdr->numSamples );
      for( int isamp = shdr->numSamples-1; isamp >= shiftSamples; isamp-- ) {
        samples[isamp] = samples[isamp-shiftSamples];
      }
      for( int isamp = 0; isamp < shiftSamples; isamp++ ) {
        samples[isamp] = 0;
      }
    }
  }

  if( vars->updateTimeHdr ) {
    int time_samp1_s  = trcHdr->intValue( vars->hdrID_time_samp1_s );
    int time_samp1_us = trcHdr->intValue( vars->hdrID_time_samp1_us );
    time_samp1_us -= (int)round(shift_ms*1000.0);
    if( time_samp1_us < 0 ) {
      time_samp1_us += 1000000;
      time_samp1_s  -= 1;
    }
    else if( time_samp1_us >= 1000000 ) {
      time_samp1_us -= 1000000;
      time_samp1_s  += 1;
    }
    trcHdr->setIntValue( vars->hdrID_time_samp1_s, time_samp1_s );
    trcHdr->setIntValue( vars->hdrID_time_samp1_us, time_samp1_us );
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_statics_( csParamDef* pdef ) {
  pdef->setModule( "STATICS", "Apply trace statics" );

  pdef->addParam( "bulk_shift", "Apply static bulk shift to all traces", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Static shift [ms]. Positive value shifts samples downwards" );

  pdef->addParam( "header", "Apply static shift from trace header", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header containing static shift [ms]. Positive value shifts samples downwards" );

  pdef->addParam( "mode", "Mode of application", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply all specified statics" );
  pdef->addOption( "remove", "Remove all specified statics", "Apply inverse/negative statics" );

  pdef->addParam( "ncoef", "Number of interpolation coefficients", NUM_VALUES_FIXED );
  pdef->addValue( "8", VALTYPE_STRING, "Number of interpolation coefficients" );

  pdef->addParam( "apply_fraction", "Apply fractional static beyond full sample interval?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Apply full static including fraction, use sinc interpolation" );
  pdef->addOption( "no", "Shift samples by N x sample interval without interpolation. Store ,fractional residual static [ms] in trace header 'stat_res'" );

  pdef->addParam( "domain", "Time or sample domain", NUM_VALUES_FIXED, "Values in trace header and/or user parameter 'bulk_shift' are provided in these units" );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Statics are specified in unit of trace (e.g. [ms] or [hz])" );
  pdef->addOption( "sample", "Statics are specified in number of samples (may be fractions)" );

  pdef->addParam( "update_time_hdr", "Update time headers according to static time shift?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Change internal time headers (time_samp1 and time_samp1_us) according to static shift in order to reflect changed time of first sample" );
  pdef->addOption( "no", "Do not update time headers" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_statics_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_statics::VariableStruct* vars = reinterpret_cast<mod_statics::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_statics_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_statics::VariableStruct* vars = reinterpret_cast<mod_statics::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->buffer ) {
    delete [] vars->buffer;
    vars->buffer = NULL;
  }
  if( vars->interpol != NULL ) {
    delete vars->interpol;
    vars->interpol = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_statics_( csParamDef* pdef ) {
  params_mod_statics_( pdef );
}
extern "C" void _init_mod_statics_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_statics_( param, env, writer );
}
extern "C" bool _start_exec_mod_statics_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_statics_( env, writer );
}
extern "C" void _exec_mod_statics_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_statics_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_statics_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_statics_( env, writer );
}

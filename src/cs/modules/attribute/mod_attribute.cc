/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_methods.h"
#include "csTableNew.h"
#include "csTableValueList.h"
#include "csGeolibUtils.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: attribute
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_attribute {
  struct VariableStruct {
    int method;
    int startSample;
    int endSample;
    int rmsWinWidthSample;
    int rmsWinStepSample;
    int hdrID_attr1;
    int hdrID_attr2;
    cseis_geolib::type_t type_attr1;
    cseis_geolib::type_t type_attr2;
    float* buffer;
    bool interpolate;

    csTableNew* table;
    int* hdrId_keys;
    int tableValueIndex;
  };
  static int const METHOD_NONE = -1;
  static int const METHOD_MINIMUM = 1;
  static int const METHOD_MAXIMUM = 2;
  static int const METHOD_RMS_MAXIMUM = 3;
  static int const METHOD_ZRUNS = 4;
  static int const METHOD_HORIZON = 5;
  static int const METHOD_AMPLITUDE = 6;
}
using mod_attribute::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_attribute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->method      = mod_attribute::METHOD_NONE;
  vars->hdrID_attr1 = -1;
  vars->hdrID_attr2 = -1;
  vars->type_attr1 = cseis_geolib::TYPE_DOUBLE;
  vars->type_attr2 = cseis_geolib::TYPE_DOUBLE;
  vars->startSample = 0;
  vars->endSample   = shdr->numSamples-1;
  vars->buffer = NULL;
  vars->rmsWinWidthSample = 0;
  vars->rmsWinStepSample  = 0;
  vars->interpolate = true;
  vars->table = NULL;
  vars->hdrId_keys = NULL;

  std::string text;

  if( param->exists( "method" ) ) {
    param->getString( "method", &text );
    text = toLowerCase( text );
    if( !text.compare( "minimum" ) ) {
      vars->method = mod_attribute::METHOD_MINIMUM;
    }
    else if( !text.compare( "maximum" ) ) {
      vars->method = mod_attribute::METHOD_MAXIMUM;
    }
    else if( !text.compare( "max_rms" ) ) {
      vars->method = mod_attribute::METHOD_RMS_MAXIMUM;
    }
    else if( !text.compare( "zruns" ) ) {
      vars->method = mod_attribute::METHOD_ZRUNS;
    }
    else if( !text.compare( "horizon" ) ) {
      vars->method = mod_attribute::METHOD_HORIZON;
    }
    else if( !text.compare( "amplitude" ) ) {
      vars->method = mod_attribute::METHOD_AMPLITUDE;
    }
    else {
      writer->line("Unknown argument for user parameter 'method': '%s'.", text.c_str());
      env->addError();
    }
  }

 //--------------------------------------------------
  if( vars->method == mod_attribute::METHOD_HORIZON ) {
    vars->table = new csTableNew( csTableNew::TABLE_TYPE_UNIQUE_KEYS );

    int numKeys = param->getNumLines("table_key");
    if( numKeys == 0 ) {
      writer->error("No table key(s) specified.");
    }
    vars->hdrId_keys  = new int[numKeys];
    for( int ikey = 0; ikey < numKeys; ikey++ ) {
      std::string headerName;
      int col;
      bool interpolate = false;
      param->getStringAtLine( "table_key", &headerName, ikey, 0 );
      param->getIntAtLine( "table_key", &col, ikey, 1 );
      if( param->getNumValues( "table_key", ikey ) > 2 ) {
        param->getStringAtLine( "table_key", &text, ikey, 2 );
        if( !text.compare("yes") ) {
          interpolate = true;
        }
        else if( !text.compare("no") ) {
          interpolate = false;
        }
        else {
          writer->error("Unknown option: %s", text.c_str() );
        }
      }
      vars->table->addKey( col-1, interpolate );  // -1 to convert from 'user' column to 'C' column
      if( !hdef->headerExists( headerName ) ) {
        writer->error("No matching trace header found for table key '%s'", headerName.c_str() );
      }
      vars->hdrId_keys[ikey] = hdef->headerIndex( headerName );
    } // END for ikey
    
    int col;
    param->getInt("table_time",&col);
    vars->table->addValue( col-1 );  // -1 to convert from 'user' column to 'C' column
    vars->tableValueIndex = 0;

    param->getString("table", &text );
    bool sortTable = false;
    try {
      vars->table->initialize( text, sortTable );
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s': %s\n", text.c_str(), exc.getMessage() );
    }
  }

  //--------------------------------------------------

  if( param->exists( "interpolate" ) ) {
    param->getString( "interpolate", &text );
    text = toLowerCase( text );
    if( !text.compare( "yes" ) ) {
      vars->interpolate = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->interpolate = false;
    }
    else {
      writer->line("Unknown option: '%s'.", text.c_str());
      env->addError();
    }
  }

  //---------------------------------------------------------

  if( param->exists("window") ) {
    double start;
    double end;
    param->getDouble( "window", &start, 0 );
    param->getDouble( "window", &end, 1 );
    vars->startSample = (int)( start / shdr->sampleInt + 0.5 );
    vars->endSample   = (int)( end / shdr->sampleInt + 0.5 );

    if( vars->startSample < 0 ) vars->startSample = 0;
    if( vars->endSample >= shdr->numSamples ) vars->endSample = shdr->numSamples-1;
  }
  if( vars->method == mod_attribute::METHOD_RMS_MAXIMUM ) {
    double width;
    double inc;
    int widthFull = vars->endSample - vars->startSample + 1;
    param->getDouble("rms",&width,0);
    param->getDouble("rms",&inc,1);
    vars->rmsWinWidthSample = (int)( width / shdr->sampleInt + 0.5 );
    vars->rmsWinStepSample  = (int)( inc / shdr->sampleInt + 0.5 );
    if( vars->rmsWinWidthSample >= widthFull ) {
      writer->warning("Specified window is shorter than specified time window");
      vars->rmsWinWidthSample = widthFull-1;
    }
    if( vars->rmsWinStepSample > vars->rmsWinWidthSample ) vars->rmsWinStepSample = vars->rmsWinWidthSample;
    if( vars->rmsWinStepSample <= 0 ) vars->rmsWinStepSample = 1;
    vars->endSample = (int)((widthFull - vars->rmsWinWidthSample) / vars->rmsWinStepSample ) * vars->rmsWinStepSample + vars->startSample;
    if( edef->isDebug() ) {
      writer->line(" Start/end sample of analysis: %d / %d, RMS window width/increment: %d / %d, Total numSamples: %d",
        vars->startSample, vars->endSample, vars->rmsWinWidthSample, vars->rmsWinStepSample, shdr->numSamples );
    }
  }

  std::string headerName_attr1("attr1");
  vars->type_attr1 = cseis_geolib::TYPE_DOUBLE;
  if( param->exists("hdr_attr1") ) {
    param->getString("hdr_attr1", &headerName_attr1, 0);
    if( param->getNumValues("hdr_attr1") > 1 ) {
      param->getString("hdr_attr1", &text, 1);
      vars->type_attr1 = csGeolibUtils::text2Type( text );
    }
  }
  std::string headerName_attr2("attr2");
  vars->type_attr2 = cseis_geolib::TYPE_DOUBLE;
  if( param->exists("hdr_attr2") ) {
    param->getString("hdr_attr2", &headerName_attr2, 0);
    if( param->getNumValues("hdr_attr2") > 1 ) {
      param->getString("hdr_attr2", &text, 1);
      vars->type_attr2 = csGeolibUtils::text2Type( text );
    }
  }

  if( !hdef->headerExists( headerName_attr1 ) ) {
    hdef->addHeader( vars->type_attr1, headerName_attr1, "Attribute 1" );
  }
  else {
    vars->type_attr1 = hdef->headerType(headerName_attr1);
  }
  if( vars->type_attr1 != cseis_geolib::TYPE_INT && vars->type_attr1 != cseis_geolib::TYPE_INT64 &&
      vars->type_attr1 != cseis_geolib::TYPE_FLOAT && vars->type_attr1 != cseis_geolib::TYPE_DOUBLE ) {
    writer->error("Attribute 1 (%s): Supported types are 'float', 'double', 'int' and 'int64'", headerName_attr1.c_str() );
  }
  if( !hdef->headerExists( headerName_attr2 ) ) {
    hdef->addHeader( vars->type_attr2, headerName_attr2, "Attribute 2" );
  }
  else {
    vars->type_attr2 = hdef->headerType(headerName_attr2);
  }
  if( vars->type_attr2 != cseis_geolib::TYPE_INT && vars->type_attr2 != cseis_geolib::TYPE_INT64 &&
      vars->type_attr2 != cseis_geolib::TYPE_FLOAT && vars->type_attr2 != cseis_geolib::TYPE_DOUBLE ) {
    writer->error("Attribute 2 (%s): Supported types are 'float', 'double', 'int' and 'int64'", headerName_attr2.c_str() );
  }

  vars->hdrID_attr1 = hdef->headerIndex( headerName_attr1 );
  vars->hdrID_attr2 = hdef->headerIndex( headerName_attr2 );

  vars->buffer = new float[5];
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_attribute_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  float* samples = trace->getTraceSamples();
  csTraceHeader* trcHdr = trace->getTraceHeader();
  double attr1 = 0.0;
  double attr2 = 0.0;

  if( vars->method == mod_attribute::METHOD_MINIMUM ) {
    int sampMin = vars->startSample;
    float min = samples[vars->startSample];
    for( int i = vars->startSample+1; i <= vars->endSample; i++ ) {
      if( samples[i] < min ) {
        min = samples[i];
        sampMin = i;
      }
    }
    attr1 = min;
    attr2 = sampMin*shdr->sampleInt;

    if( vars->interpolate ) {
      int start = std::max(sampMin-2,0);
      int sampMinNew = sampMin-start;
      int end   = std::min(start+5,shdr->numSamples-1);
      for( int i = start; i < end; i++ ) {
        vars->buffer[i-start] = -samples[i];
      }
      float sampleIndex;
      float maxAmplitude;
        sampleIndex = getQuadMaxSample( vars->buffer, sampMinNew, 5, &maxAmplitude );
      attr1 = -maxAmplitude;
      attr2 = (sampleIndex+start)*shdr->sampleInt;
    }
  }
  else if( vars->method == mod_attribute::METHOD_MAXIMUM ) {
    int sampMax = vars->startSample;
    float max = samples[vars->startSample];
    for( int i = vars->startSample+1; i <= vars->endSample; i++ ) {
      if( samples[i] > max ) {
        max = samples[i];
        sampMax = i;
      }
    }
    attr1 = max;
    attr2 = sampMax*shdr->sampleInt;

    if( vars->interpolate ) {
      float sampleIndex = sampMax;
      float maxAmplitude = samples[sampMax];
    
      sampleIndex = getQuadMaxSample( samples, sampMax, shdr->numSamples, &maxAmplitude );

      attr1 = maxAmplitude;
      attr2 = sampleIndex*shdr->sampleInt;
    }
  }
  else if( vars->method == mod_attribute::METHOD_AMPLITUDE ) {
    attr1 = samples[vars->startSample];
    attr2 = samples[vars->endSample];
  }
  //--------------------------------------------------------------
  else if( vars->method == mod_attribute::METHOD_RMS_MAXIMUM ) {
    int sampMax = vars->startSample;
    float max = 0;
    for( int isamp = vars->startSample; isamp <= vars->endSample; isamp += vars->rmsWinStepSample ) {
      double rms = 0;
      for( int i = isamp; i < isamp+vars->rmsWinWidthSample; i++ ) {
        rms += samples[i]*samples[i];
      }
      rms = sqrt(rms/(double)vars->rmsWinWidthSample);
      if( rms > max ) {
        max = rms;
        sampMax = isamp + vars->rmsWinWidthSample/2;
      }
    }
    attr1 = max;
    attr2 = sampMax*shdr->sampleInt;
  }
  else if( vars->method == mod_attribute::METHOD_ZRUNS ) {
    int count1 = 0;
    int count2 = 0;
    for( int isamp = vars->startSample+1; isamp <= vars->endSample; isamp += 1 ) {
      if( fabs(samples[isamp]) == fabs(samples[isamp-1]) ) {
        count2 += 1;
        if( samples[isamp] == samples[isamp-1] ) count1 += 1;
      }
    }
    attr1 = (double)count1;
    attr2 = (double)count2;
  }
  else if( vars->method == mod_attribute::METHOD_HORIZON ) {
    double* keyValueBuffer = new double[vars->table->numKeys()];
    for( int ikey = 0; ikey < vars->table->numKeys(); ikey++ ) {
      keyValueBuffer[ikey] = trcHdr->doubleValue( vars->hdrId_keys[ikey] );
    }
    attr1 = 0;
    attr2 = 0;
    try {
      attr2 = vars->table->getValue( keyValueBuffer );
    }
    catch( csException& e ) {
      delete [] keyValueBuffer;
      writer->error("Error occurred in INSERT_DATA: %s", e.getMessage());
      throw(e);
    }
    delete [] keyValueBuffer;
    float sampleIndex = attr2 / shdr->sampleInt;
    attr1 = getQuadAmplitudeAtSample( samples, sampleIndex, shdr->numSamples );
  }

  trcHdr->setDoubleValue( vars->hdrID_attr1, attr1 );
  trcHdr->setDoubleValue( vars->hdrID_attr2, attr2 );
/*
  int start = MAX(sampMin-2,0);
  int sampMinNew = sampMin-start;
  int end   = MIN(start+5,shdr->numSamples-1);
  for( int i = start; i <= end; i++ ) {
    vars->buffer[i] = -samples[i];
  }

  float maxAmplitude;
  float sampleIndex = getQuadMaxSample( vars->buffer, sampMinNew, 5, &maxAmplitude );

  trace->getTraceHeader()->setDoubleValue( vars->hdrID_attr1, -maxAmplitude );
  trace->getTraceHeader()->setDoubleValue( vars->hdrID_attr2, sampleIndex*shdr->sampleInt );
*/

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_attribute_( csParamDef* pdef ) {
  pdef->setModule( "ATTRIBUTE", "Extract attribute from data" );

  pdef->addParam( "method", "Method", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "minimum", "Extract minimum from data. Store minimum amplitude (attr1) and time/frequency (attr2)." );
  pdef->addOption( "maximum", "Extract maximum from data. Store maximum amplitude (attr1) and time/frequency (attr2)." );
  pdef->addOption( "max_rms", "Extract maximum RMS from data. Store maximum RMS amplitude (attr1) and time/frequency (attr2).", "RMS is computed over sliding window, given by user parameter 'rms'" );
  pdef->addOption( "zruns", "Count how many adjacent samples have the same sample value." );
  pdef->addOption( "horizon", "Extract amplitude at horizon time.", "Specify input horizon in user parameter 'table'" );
  pdef->addOption( "amplitude", "Extract amplitude at two constant times (at window start and end time)" );

  pdef->addParam( "window", "Computation window", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Start time/frequency [ms] or [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "End time/frequency [ms] or [Hz]" );

  pdef->addParam( "rms", "Definition of RMS window", NUM_VALUES_FIXED, "...in units of trace, e.g. [ms] or [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Width of RMS window." );
  pdef->addValue( "0", VALTYPE_NUMBER, "Increment of RMS window. Set to 0 for sliding window (one sample)." );

  pdef->addParam( "interpolate", "Interpolate max/min between full samples...?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Find min/max amplitude using quadratic interpolation" );
  pdef->addOption( "no", "Pick min/max amplitude only on full samples without interpolation" );

  pdef->addParam( "table", "Table containing 'horizon'", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Name of table containing horizon definition");

  pdef->addParam( "table_key", "Key trace header used to match values found in specified table columns", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name of key header" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table" );
  pdef->addValue( "no", VALTYPE_OPTION, "Interpolate based to this key?" );
  pdef->addOption( "yes", "Use this key for interpolation of value" );
  pdef->addOption( "no", "Do not use this key for interpolation", "The input table is expected to contain the exact key values for this trace header" );

  pdef->addParam( "table_time", "Table column containing seismic time (or depth, frequency..)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing vertical dimension value" );

  pdef->addParam( "hdr_attr1", "Name of trace header where first attribute shall be stored", NUM_VALUES_VARIABLE, "...may be new or existing trace header" );
  pdef->addValue( "attr1", VALTYPE_STRING );
  pdef->addValue( "double", VALTYPE_OPTION, "Trace header type" );
  pdef->addOption( "int", "Integer, 4 byte integer type" );
  pdef->addOption( "int64", "Integer, 8 byte integer type" );
  pdef->addOption( "float", "Float, 4 byte floating point type" );
  pdef->addOption( "double", "Double, 8 byte floating point type" );

  pdef->addParam( "hdr_attr2", "Name of trace header where second attribute shall be stored (if applicable)", NUM_VALUES_VARIABLE, "...may be new or existing trace header" );
  pdef->addValue( "attr2", VALTYPE_STRING );
  pdef->addValue( "double", VALTYPE_OPTION, "Trace header type" );
  pdef->addOption( "int", "Integer, 4 byte integer type" );
  pdef->addOption( "int64", "Integer, 8 byte integer type" );
  pdef->addOption( "float", "Float, 4 byte floating point type" );
  pdef->addOption( "double", "Double, 8 byte floating point type" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_attribute_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_attribute::VariableStruct* vars = reinterpret_cast<mod_attribute::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_attribute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_attribute::VariableStruct* vars = reinterpret_cast<mod_attribute::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->buffer != NULL ) {
    delete [] vars->buffer;
    vars->buffer = NULL;
  }
  if( vars->table != NULL ) {
    delete vars->table;
    vars->table = NULL;
  }
  if( vars->hdrId_keys != NULL ) {
    delete [] vars->hdrId_keys;
    vars->hdrId_keys = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_attribute_( csParamDef* pdef ) {
  params_mod_attribute_( pdef );
}
extern "C" void _init_mod_attribute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_attribute_( param, env, writer );
}
extern "C" bool _start_exec_mod_attribute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_attribute_( env, writer );
}
extern "C" void _exec_mod_attribute_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_attribute_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_attribute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_attribute_( env, writer );
}

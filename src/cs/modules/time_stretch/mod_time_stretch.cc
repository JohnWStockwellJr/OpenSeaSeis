/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csTimeStretch.h"
#include "csTable.h"
#include "csRSFHeader.h"
#include "csRSFReader.h"
#include <cmath>
#include <cstring>
#include <cstdio>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TIME_STRETCH
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_time_stretch {
  struct VariableStruct {
    float* times;    // [s]
    float* values;  // [ms]
    csTimeStretch* timeStretchObj;
    int numTimes;
    float* traceBuffer;
    int mode;
    int method;

    int hdrId_stretch;
    int* hdrId_time_in;
    int* hdrId_time_out;
    int hdrId_vel_dim2;
    int hdrId_vel_dim3;

    csTable const* table; // Do not allocate or free!
    int numKeys;
    int* hdrId_keys;
    double* keyValueBuffer;

    cseis_io::csRSFReader* rsfReader;
    cseis_io::csRSFHeader* rsfHdr;
    std::string filename_vel;
    float sampleInt_in;
    int numSamples_in;
    bool isOutVel;
  };
  static int const METHOD_STRETCH_FACTOR = 1;
  static int const METHOD_HORIZON_TIME   = 2;
  static int const METHOD_TIME2DEPTH     = 3;
  static int const METHOD_DEPTH2TIME     = 4;

  static const int MODE_APPLY  = 11;
  static const int MODE_REMOVE = 12;
}
using mod_time_stretch::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_time_stretch_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

//---------------------------------------------
//
  vars->numTimes       = 0;
  vars->times          = NULL;
  vars->values          = NULL;
  vars->timeStretchObj = NULL;
  vars->traceBuffer    = NULL;
  vars->hdrId_stretch  = -1;
  vars->method      = mod_time_stretch::METHOD_STRETCH_FACTOR;

  vars->table          = NULL;
  vars->hdrId_keys     = NULL;
  vars->keyValueBuffer = NULL;
  vars->numKeys        = 0;
  vars->mode           = mod_time_stretch::MODE_APPLY;

  vars->hdrId_time_in  = NULL;
  vars->hdrId_time_out = NULL;

  vars->rsfReader = NULL;
  vars->rsfHdr    = NULL;
  vars->hdrId_vel_dim2  = -1;
  vars->hdrId_vel_dim3  = -1;
  vars->numSamples_in = shdr->numSamples;
  vars->sampleInt_in  = shdr->sampleInt;
  vars->isOutVel = false;

  csVector<std::string> valueList;

  //-------------------------------------
  std::string text;

  if( param->exists("mode") ) {
    param->getString( "mode", &text );
    if( !text.compare("apply") ) {
      vars->mode = mod_time_stretch::MODE_APPLY;
    }
    else if( !text.compare("remove") ) {
      vars->mode = mod_time_stretch::MODE_REMOVE;
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }
  }

  param->getString( "method", &text );
  if( !text.compare("stretch") ) {
    vars->method = mod_time_stretch::METHOD_STRETCH_FACTOR;
  }
  else if( !text.compare("time") ) {
    vars->method = mod_time_stretch::METHOD_HORIZON_TIME;
  }
  else { // depth2time or time2depth conversion
    param->getString( "filename_vel", &vars->filename_vel );
    try {
      vars->rsfReader = new cseis_io::csRSFReader( vars->filename_vel, 1, false );
    }
    catch( csException& e ) {
      vars->rsfReader = NULL;
      writer->error("Error when opening RSF file '%s'.\nSystem message: %s", vars->filename_vel.c_str(), e.getMessage() );
    }
    vars->rsfHdr    = new cseis_io::csRSFHeader();
    vars->rsfReader->initialize( vars->rsfHdr );
    vars->rsfReader->dump( writer->getFile() );

    if( !text.compare("time2depth") ) {
      vars->method = mod_time_stretch::METHOD_TIME2DEPTH;
    }
    else if( !text.compare("depth2time") ) {
      vars->method = mod_time_stretch::METHOD_DEPTH2TIME;
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }

    param->getString( "hdr_vel", &text, 0 );
    vars->hdrId_vel_dim2  = hdef->headerIndex( text );
    param->getString( "hdr_vel", &text, 1 );
    vars->hdrId_vel_dim3  = hdef->headerIndex( text );

    param->getFloat( "sample_int_out", &shdr->sampleInt );
    param->getInt( "nsamples_out", &shdr->numSamples );
  }

//---------------------------------------------
// Retrieve velocity table
//
/*
  if( param->exists("table") ) {
    param->getString("table", &text );
    vars->table = env->getTable( text );
    if( vars->table == NULL ) {
      writer->error("No stretch table found with name '%s'. The table must be specified within the flow file with directive &table <tablename> <filename> {methodInterpolation}", text.c_str() );
    }
    vars->numKeys = vars->table->numKeys();
    if( vars->table->dimension() != csTable::TABLE_DIM_1D ) {
      writer->error("Stretch table '%s' has a time column which makes it a 2D table.", text.c_str() );
    }
    vars->hdrId_keys     = new int[vars->numKeys];
    vars->keyValueBuffer = new double[vars->numKeys];
    for( int ikey = 0; ikey < vars->numKeys; ikey++ ) {
      if( !hdef->headerExists( vars->table->keyName(ikey) ) ) {
        writer->error("No matching trace header found for table key '%s'", vars->table->keyName(ikey).c_str() );
      }
      vars->hdrId_keys[ikey] = hdef->headerIndex( vars->table->keyName(ikey) );
    }
  }
*/
//---------------------------------------------
// Retrieve time-stretch pairs
//
  if( param->exists( "time" ) ) {
    param->getAll( "time", &valueList );
    if( valueList.size() < 1 ){
      writer->error("No times specified in user parameter 'time'");
    }
  }
  vars->numTimes = valueList.size();

  if( vars->numTimes > 0 ) {
    vars->times = new float[vars->numTimes];
    vars->values = new float[vars->numTimes];
    for( int i = 0; i < vars->numTimes; i++ ) {
      vars->times[i] = (float)atof( valueList.at(i).c_str() );
    }
  }

  if( param->exists( "hdr_stretch" ) ) {
    if( vars->method != mod_time_stretch::METHOD_STRETCH_FACTOR ) {
      writer->error("User parameter 'hdr_stretch' can only be specified in conjunction with value type 'stretch'");
    }
    param->getString( "hdr_stretch", &text );
    if( !hdef->headerExists( text ) ) {
      writer->error("Unknown trace header: '%s'", text.c_str());
    }
    vars->hdrId_stretch = hdef->headerIndex( text );
  }
  else if( param->exists( "hdr_time_in" ) ) {
    if( vars->numTimes > 0 ) writer->error("Specify either a list of time/value pairs, or parameter 'hdr_time', not both.");
    vars->numTimes = param->getNumValues("hdr_time_in");
    int numtmp = param->getNumValues("hdr_time_out");
    if( numtmp != vars->numTimes ) writer->error("Specify same number of input and output time headers: 'hdr_time_in' and 'hdr_time_out'");
    vars->hdrId_time_in = new int[vars->numTimes];
    vars->hdrId_time_out = new int[vars->numTimes];
    for( int i = 0; i < vars->numTimes; i++ ) {
      param->getString("hdr_time_in", &text, i);
      if( !hdef->headerExists( text ) ) {
        writer->error("Unknown trace header: '%s'", text.c_str());
      }
      vars->hdrId_time_in[i] = hdef->headerIndex( text );
      param->getString("hdr_time_out", &text, i);
      if( !hdef->headerExists( text ) ) {
        writer->error("Unknown trace header: '%s'", text.c_str());
      }
      vars->hdrId_time_out[i] = hdef->headerIndex( text );
    }
    vars->times  = new float[vars->numTimes];
    vars->values = new float[vars->numTimes];
  }
  else if( vars->method == mod_time_stretch::METHOD_STRETCH_FACTOR || vars->method == mod_time_stretch::METHOD_HORIZON_TIME ) {  // No stretch/time headers provided
    valueList.clear();
    if( vars->method == mod_time_stretch::METHOD_STRETCH_FACTOR ) {
      param->getAll( "stretch", &valueList );
    }
    else if( vars->method == mod_time_stretch::METHOD_HORIZON_TIME ) {
      param->getAll( "time_out", &valueList );
    }
    int numValues = valueList.size();
    if( valueList.size() < 1 ) {
      writer->error("No values supplied for user parameter 'value'");
    }
    else if( vars->method == mod_time_stretch::METHOD_STRETCH_FACTOR && numValues != vars->numTimes-1 ) {
      writer->error("Non-matching number of stretch factors (=%d) and horizon times (=%d). Each stretch factor applies to one layer in between two horizon times, so the numberr of stretch factors should be one less",
        numValues, vars->numTimes );
    }
    else if( vars->method == mod_time_stretch::METHOD_HORIZON_TIME && numValues != vars->numTimes ) {
      writer->error("Non-matching number of output (=%d) and input horizon times (=%d).", numValues, vars->numTimes );
    }
    float sign_scalar = ( vars->mode == mod_time_stretch::MODE_APPLY ) ? 1.0 : -1.0;
    csFlexNumber number;
    for( int i = 0; i < numValues; i++ ) {
      if( !number.convertToNumber( valueList.at(i) ) ) {
        writer->error("Specified value is not a valid number: '%s'", valueList.at(i).c_str() );
      }
      vars->values[i] = sign_scalar * number.floatValue();
      if( edef->isDebug() ) writer->line("Value #%d: '%s' --> %f", i, valueList.at(i).c_str(), vars->values[i] );
    }
  }
  else { // Velocity model for time2depth or depth2time conversion
    //
  }

  //---------------------------------------------
  // Set TimeStretch object
  //
  vars->timeStretchObj = new csTimeStretch( vars->sampleInt_in, vars->numSamples_in );
  vars->traceBuffer    = new float[vars->numSamples_in];

  if( vars->hdrId_time_in == NULL ) {
    if( vars->method == mod_time_stretch::METHOD_STRETCH_FACTOR ) {
      for( int i = 0; i < vars->numTimes-1; i++ ) {
        writer->line("Horizon: %d  Top: %f  Bottom: %f   Stretch: %f", i, vars->times[i], vars->times[i+1], vars->values[i] );
      }
    }
    else if( vars->method == mod_time_stretch::METHOD_HORIZON_TIME ) {
      for( int i = 0; i < vars->numTimes; i++ ) {
        writer->line("Horizon: %d  Input: %f   Output: %f", i, vars->times[i], vars->values[i] );
      }
    }
  }
  if( vars->method == mod_time_stretch::METHOD_TIME2DEPTH || vars->method == mod_time_stretch::METHOD_DEPTH2TIME ) {
    vars->timeStretchObj->initialize_timeDepthConversion( shdr->numSamples, shdr->sampleInt );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_time_stretch_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);
  csTraceHeader* trcHdr = trace->getTraceHeader();

  if( vars->hdrId_time_in != NULL ) {
    for( int i = 0; i < vars->numTimes; i++ ) {
      vars->times[i]  = trcHdr->floatValue( vars->hdrId_time_in[i] );
      vars->values[i] = trcHdr->floatValue( vars->hdrId_time_out[i] );
    }
  }
  float* samples = trace->getTraceSamples();
  memcpy( vars->traceBuffer, samples, vars->numSamples_in*sizeof(float) );

  if( vars->hdrId_stretch > 0 ) {
    float s = trace->getTraceHeader()->floatValue( vars->hdrId_stretch );
    float sign_scalar = ( vars->mode == mod_time_stretch::MODE_APPLY ) ? 1.0 : -1.0;
    for( int i = 0; i < vars->numTimes-1; i++ ) {
      vars->values[i] = sign_scalar * s;
      if( edef->isDebug() ) writer->line("TIME_STRETCH DEBUG: Stretch #%d: %f", i, vars->values[i]);
    }
  }
  if( vars->method == mod_time_stretch::METHOD_STRETCH_FACTOR ) {
    vars->timeStretchObj->applyStretchFunction( vars->traceBuffer, vars->times, vars->values, vars->numTimes, samples );
  }
  else if( vars->method == mod_time_stretch::METHOD_HORIZON_TIME ) {
    vars->timeStretchObj->applyTimeInterval( vars->traceBuffer, vars->times, vars->values, vars->numTimes, samples );
  }
  else {
    double val_dim2 = trcHdr->doubleValue( vars->hdrId_vel_dim2 );
    double val_dim3 = trcHdr->doubleValue( vars->hdrId_vel_dim3 );
    if( val_dim2 < vars->rsfHdr->o2 ) val_dim2 = vars->rsfHdr->o2;
    if( val_dim2 > vars->rsfHdr->e2 ) val_dim2 = vars->rsfHdr->e2;
    if( val_dim3 < vars->rsfHdr->o3 ) val_dim3 = vars->rsfHdr->o3;
    if( val_dim3 > vars->rsfHdr->e3 ) val_dim3 = vars->rsfHdr->e3;
    //    int traceIndex = vars->rsfReader->computeTrace( val_dim2, val_dim3 );
    //    fprintf(stderr,"Trace index %d (%d):  %.2f %.2f\n", traceIndex, vars->rsfReader->numTraces(), val_dim2, val_dim3);
    vars->rsfReader->moveToComputedTrace( val_dim2, val_dim3 );
    float const* samples_vel = vars->rsfReader->getNextTracePointer();
    float sampleInt_velFunc = vars->rsfReader->sampleInt();
    int numSamples_velFunc  = vars->rsfReader->numSamples();

    vars->timeStretchObj->apply_timeDepthConversion( vars->traceBuffer, samples_vel, sampleInt_velFunc, numSamples_velFunc, samples, (vars->method == mod_time_stretch::METHOD_TIME2DEPTH) );
  }
  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_time_stretch_( csParamDef* pdef ) {
  pdef->setModule( "TIME_STRETCH", "Stretch/squeeze data trace, time/depth conversion" );

  pdef->addDoc("Specify either pairs of horizon input/output times, or otherwise the input horizon times and a stretch factor for each layer.");
  pdef->addDoc("For time/depth conversion, specify a velocity field, two trace headers defining the trace location (e.g. row and col), and the output sample interval & number of samples.");

  pdef->addParam( "method", "Time/depth stretch method to apply", NUM_VALUES_FIXED );
  pdef->addValue( "stretch_factor", VALTYPE_OPTION );
  pdef->addOption( "stretch_factor", "Apply time stretch using a stretch factor. Provide time/stretch factor pairs in user parameter 'time' and 'stretch' or through trace header specified in 'hdr_stretch'" );
  pdef->addOption( "horizon_time", "Apply time stretch using horizon times. Provide input/output horizon time(s) in user parameters 'time' and 'time_out' or through trace header specified in 'hdr_time_in' and 'hdr_time_out'" );
  pdef->addOption( "time2depth", "Apply time-to-depth stretch using velocty function specified in user parammeter 'filename_vel'" );
  pdef->addOption( "depth2time", "Apply depth-to-time stretch using velocty function specified in user parammeter 'filename_vel'" );

  pdef->addParam( "time", "List of input horizon times [ms]", NUM_VALUES_VARIABLE,
                  "For a model with N layers, specify N+1 time values. First value gives top of first layer, each subsequent value defines the bottom of one layer and the top of the next layer" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of input horizon times [ms]." );

  pdef->addParam( "time_out", "List of output horizon times [ms]", NUM_VALUES_VARIABLE, "Number of output times must be the same as the number of input horizon times" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of output horizon times [ms]." );

  pdef->addParam( "stretch", "List of stretch factors [ms]", NUM_VALUES_VARIABLE, "Number of stretch factors must be one less than number of times, one for each model layer.");
  pdef->addValue( "", VALTYPE_NUMBER, "Stretch factor [ms]. Polarity for stretch factor: Stretch(+) or squeeze(-)" );

  pdef->addParam( "hdr_stretch", "Trace header containing stretch value", NUM_VALUES_FIXED, "This is an alternative to providing a list of values. This only applies to a single layer" );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header containing stretch value in [ms] (or [Hz])" );

  pdef->addParam( "hdr_time_in", "(List of) Trace header containing time value for one horizon", NUM_VALUES_FIXED, "This is an alternative to providing a list of times and values." );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header containing time value in [ms] (or [Hz])" );
  pdef->addParam( "hdr_time_out", "(List of) Trace header containing time value for one horizon", NUM_VALUES_FIXED, "This is an alternative to providing a list of times and values." );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header containing time value in [ms] (or [Hz])" );

  pdef->addParam( "hdr_vel", "Trace headers defining location in 3D velocity field", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header / rsf dimension 2" );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header / rsf dimension 3" );

  pdef->addParam( "sample_int_out", "Output sample interval (depth[m] or time[ms]) for depth/time conversion", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER );

  pdef->addParam( "nsamples_out", "Output number of samples for depth/time conversion", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER );

  pdef->addParam( "filename_vel", "Velocity field for time<-->depth conversion", NUM_VALUES_FIXED, "Currently, only 'rsf' format is supported for velocity field input" );
  pdef->addValue( "", VALTYPE_STRING );

  pdef->addParam( "mode", "Mode of application", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply" );
  pdef->addOption( "remove", "Remove" );

  pdef->addParam( "output", "Output data (only applicable for time/depth conversion)", NUM_VALUES_FIXED );
  pdef->addValue( "default", VALTYPE_OPTION );
  pdef->addOption( "default", "Output default data set" );
  pdef->addOption( "velocity", "Output velocity function", "This overrides the specified output sample interval & number of samples" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_time_stretch_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_time_stretch::VariableStruct* vars = reinterpret_cast<mod_time_stretch::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_time_stretch_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_time_stretch::VariableStruct* vars = reinterpret_cast<mod_time_stretch::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->times != NULL ) {
    delete [] vars->times;
    vars->times = NULL;
  }
  if( vars->values != NULL ) {
    delete [] vars->values;
    vars->values = NULL;
  }
  if( vars->timeStretchObj == NULL ) {
    delete vars->timeStretchObj;
    vars->timeStretchObj = NULL;
  }
  if( vars->hdrId_time_in != NULL ) {
    delete [] vars->hdrId_time_in;
    vars->hdrId_time_in = NULL;
  }
  if( vars->hdrId_time_out != NULL ) {
    delete [] vars->hdrId_time_out;
    vars->hdrId_time_out = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_time_stretch_( csParamDef* pdef ) {
  params_mod_time_stretch_( pdef );
}
extern "C" void _init_mod_time_stretch_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_time_stretch_( param, env, writer );
}
extern "C" bool _start_exec_mod_time_stretch_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_time_stretch_( env, writer );
}
extern "C" void _exec_mod_time_stretch_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_time_stretch_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_time_stretch_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_time_stretch_( env, writer );
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include <cmath>

#include "csTimeFunction.h"
#include "csTableNew.h"
#include "csTableValueList.h"


using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SCALING
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_scaling {
  struct VariableStruct {
    int option;
    int numTimes;
    float* sampIndex;
    float* scalar;
    int hdrId;
    bool isFile;
    bool isScalar;
    float* bufferScalar;
    int windowStartSample;
    int windowEndSample;

    csTableNew* table;
    int* hdrId_keys;
  };
}
using mod_scaling::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_scaling_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef  = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->option    = 0;
  vars->numTimes  = 0;
  vars->sampIndex = NULL;
  vars->scalar    = NULL;
  vars->hdrId     = -1;
  vars->bufferScalar = NULL;
  vars->isFile = false;
  vars->isScalar = false;
  vars->windowStartSample = 0; // By default, apply mute to full trace length
  vars->windowEndSample   = shdr->numSamples-1; // By default, apply mute to full trace length

  vars->table          = NULL;
  vars->hdrId_keys     = NULL;

//---------------------------------------------
// Retrieve sampIndex and scalar
//
  csVector<std::string> valueList;

  std::string text;
  if( param->exists( "filename" ) ) {
    param->getString("filename",&text);
    csVector<float> timeList;
    csVector<float> scalarList;
    FILE* fin = fopen(text.c_str(),"r");
    if( fin == NULL ) writer->error("Error occurred when opening input file %s.", text.c_str());
    char buffer[132];
    while( fgets(buffer,132,fin) != NULL ) {
      float time;
      float scalar;
      sscanf(buffer,"%f %f", &time, &scalar);
      timeList.insertEnd(time);
      scalarList.insertEnd(scalar);
      int num = timeList.size();
      if( num > 1 && timeList.at(num-1) <= timeList.at(num-2) ) writer->error("Time #%d not monotonically increasing: %f(#%d) --> %f(#%d)",
                                                                           num, timeList.at(num-2), num-1, timeList.at(num-1), num);
    }
    fclose(fin);
    if( scalarList.size() == 0 ) writer->error("Input file %s does not contain any data", text.c_str());
    int numSamples = scalarList.size();
    vars->bufferScalar = new float[shdr->numSamples];
    int indexCurrent = 0;
    float timeCurrent   = timeList.at(0);
    float timePrev      = timeCurrent;
    float scalarCurrent = scalarList.at(0);
    float scalarPrev    = scalarCurrent;
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      float time = (float)isamp * shdr->sampleInt;
      if( (indexCurrent == 0 && time <= timeCurrent) || (indexCurrent == numSamples-1) ) {
        vars->bufferScalar[isamp] = scalarCurrent;
      }
      else {
        if( time > timeCurrent ) {
          do {
            indexCurrent += 1;
            timeCurrent  = timeList.at(indexCurrent);  
          } while( time > timeCurrent && indexCurrent < numSamples-1 );
          timePrev   = timeList.at(indexCurrent-1);
          scalarPrev = scalarList.at(indexCurrent-1);
          scalarCurrent = scalarList.at(indexCurrent);
        }
        if( timePrev != timeCurrent ) {
          vars->bufferScalar[isamp] = scalarPrev + (time-timePrev)/(timeCurrent-timePrev) * (scalarCurrent - scalarPrev);
        }
        else {
          vars->bufferScalar[isamp] = scalarCurrent;
        }
      }
      if( edef->isDebug() ) fprintf(stdout,"%d %f %e\n", isamp, time, vars->bufferScalar[isamp]);
    }
    vars->isFile = true;

  }

//---------------------------------------------
// Retrieve scalars from ASCII table
//
  if( param->exists("table") ) {
    if( vars->isFile ) writer->error("Specify either user parameter 'table' or 'filename', not both");
    std::string tableFilename;
    int colTime = 0;
    int colScalar  = 1;
    param->getString("table", &tableFilename );
    param->getInt("table_col",&colTime,0);
    param->getInt("table_col",&colScalar,1);
    if( colTime < 1 || colScalar < 1 ) writer->error("Column numbers in table (user parameter 'table_col') must be larger than 0");
    vars->table = new csTableNew( csTableNew::TABLE_TYPE_TIME_FUNCTION, colTime-1 );
    int numKeys = param->getNumLines("table_key");
    if( numKeys > 0 ) {
      int col;
      vars->hdrId_keys  = new int[numKeys];
      for( int ikey = 0; ikey < numKeys; ikey++ ) {
        std::string headerName;
        bool interpolate = true;
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
        vars->table->addKey( col-1, interpolate );  // -1 to convert from 'user' column to 'C++' column
        if( !hdef->headerExists( headerName ) ) {
          writer->error("No matching trace header found for table key '%s'", headerName.c_str() );
        }
        vars->hdrId_keys[ikey] = hdef->headerIndex( headerName );
      } // END for ikey
    }

    // Add scalar 'value' column to input table
    vars->table->addValue( colScalar-1 );  // -1 to convert from 'user' column to 'C++' column
    bool sortTable = false;
    try {
      vars->table->initialize( tableFilename, sortTable );
      //      vars->table->dump();
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s': %s\n", tableFilename.c_str(), exc.getMessage() );
    }
  }
  //--------------------------------------------------------------------------------

  if( param->exists( "time" ) ) {
    param->getAll( "time", &valueList );
    if( valueList.size() < 1 ){
      writer->error("No times specified in user parameter 'TIME'!");
    }
  }
  vars->numTimes = valueList.size();

  if( vars->numTimes > 0 ) {
    vars->sampIndex = new float[vars->numTimes];
    for( int i = 0; i < vars->numTimes; i++ ) {
      vars->sampIndex[i] = atof( valueList.at(i).c_str() ) / shdr->sampleInt;  // Convert to sample index
      //      fprintf(stderr,"Sampindex:  %d %f\n", i, vars->sampIndex[i] );
    }
  }
  else {
    vars->numTimes = 1;
    vars->sampIndex = new float[vars->numTimes];
    vars->sampIndex[0] = 0;
  }

  vars->scalar = new float[vars->numTimes];

  valueList.clear();

  if( param->exists("header") ) {
    string name;
    param->getString("header",&name);
    vars->hdrId = hdef->headerIndex( name );
  }
  else if( param->exists("scalar") ) {
    vars->isScalar = true; 
    param->getAll( "scalar", &valueList );

    if( valueList.size() < 1 ){
      writer->error("Missing user parameter 'scalar'!");
    }
    else if( valueList.size() != vars->numTimes ) {
      writer->error("Unequal number of scalar(%d) and times(%d)", valueList.size(), vars->numTimes );
    }
    csFlexNumber number;
    for( int i = 0; i < vars->numTimes; i++ ) {
//    vars->scalar[i] = atof( valueList.at(i).c_str() );
      if( !number.convertToNumber( valueList.at(i) ) ) {
        writer->error("Specified scalar is not a valid number: '%s'", valueList.at(i).c_str() );
      }
      vars->scalar[i] = number.floatValue();
      if( edef->isDebug() ) writer->line("Scalar #%d: '%s' --> %f", i, valueList.at(i).c_str(), vars->scalar[i] );
    }
  }
  else if( !vars->isFile && vars->table == NULL) {
    writer->error("Need to specify at least one scaling parameter: 'header', 'scalar', 'filename' or 'table'");
  }
//----------------------------------------------
  
  std::string tableName;
  if( param->exists("table") ) {
    param->getString( "table", &tableName );
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
    writer->line("Window start/end samples: %d - %d  (Times: %.3f - %.3f)",
                 vars->windowStartSample, vars->windowEndSample, vars->windowStartSample*shdr->sampleInt, vars->windowEndSample*shdr->sampleInt );
  }

//----------------------------------------------
  if( edef->isDebug() ) {
    for( int i = 0; i < vars->numTimes; i++ ) {
      writer->line("sampleIndex: %8.2f , time: %8.2f, scalar: %8.4f", vars->sampIndex[i], vars->sampIndex[i]*shdr->sampleInt,
        vars->scalar[i] );
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_scaling_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  float* scalar = vars->scalar;
  float* sampIndex = vars->sampIndex;
  float* samples = trace->getTraceSamples();

  int firstSampleToProcess = vars->windowStartSample;
  int lastSampleToProcess  = vars->windowEndSample;
  //  fprintf(stderr,"First/last sample to process: %d %d   %f %f\n", firstSampleToProcess, lastSampleToProcess, sampIndex[0], sampIndex[1]);

  if( vars->hdrId >= 0 ) {
    double scalarValue = trace->getTraceHeader()->doubleValue(vars->hdrId);
    for( int isamp = firstSampleToProcess; isamp <= lastSampleToProcess; isamp++ ) {
      samples[isamp] *= scalarValue;
    }
  }
  else if( vars->isScalar ) {
    float sc;
    if( vars->numTimes == 1 ) {
      sc = scalar[0];
      for( int isamp = firstSampleToProcess; isamp <= lastSampleToProcess; isamp++ ) {
        samples[isamp] *= sc;
      }
    }
    else {
      int sampIndexFirst = (int)( sampIndex[0] ) + 1;
      int sampIndexLast  = (int)( sampIndex[vars->numTimes-1] + 0.01 );
      
      sampIndexFirst = std::min( std::max( sampIndexFirst, firstSampleToProcess ), lastSampleToProcess );
      //      sampIndexLast  = std::max( std::min( sampIndexLast, lastSampleToProcess ), lastSampleToProcess );
      sampIndexLast  = std::min( sampIndexLast, lastSampleToProcess );
      //      fprintf(stderr,"First/last index: %d %d\n", sampIndexFirst, sampIndexLast );

      sc = scalar[0];
      for( int isamp = firstSampleToProcess; isamp < sampIndexFirst; isamp++ ) {
        samples[isamp] *= sc;
      }
      sc = scalar[vars->numTimes-1];
      for( int isamp = sampIndexLast; isamp <= lastSampleToProcess; isamp++ ) {
        samples[isamp] *= sc;
      }
      
      int index = 1;
      float weight;
      float norm = 1.0;
      bool hasChanged = true;
      //      fprintf(stderr,"First/last index: %d %d\n", sampIndexFirst, sampIndexLast );
      for( int isamp = sampIndexFirst; isamp < sampIndexLast; isamp++ ) {
        while( isamp > sampIndex[index] ) {
          hasChanged = true;
          index++;
        }
        if( hasChanged ) {
          norm = sampIndex[index] - sampIndex[index-1];
          hasChanged = false;
        }
        //        fprintf(stdout,"%d %d   %d   - %f %f\n", index, vars->numTimes, hasChanged, sampIndex[index-1], sampIndex[index] );
        weight = ((float)isamp-sampIndex[index-1])/norm;
        sc = ( (1.0-weight) * scalar[index-1] + weight * scalar[index] );
        samples[isamp] *= sc;
      }
    }
  }
  if( vars->isFile ) {
    for( int isamp = firstSampleToProcess; isamp <= lastSampleToProcess; isamp++ ) {
      samples[isamp] *= vars->bufferScalar[isamp];
    }
  }
  else if( vars->table != NULL ) {
    bool dump = false;
    csTimeFunction<double> const* timeFunc;
    if( vars->table->numKeys() > 0 ) {
      double* keyValueBuffer = new double[vars->table->numKeys()];
      for( int ikey = 0; ikey < vars->table->numKeys(); ikey++ ) {
        keyValueBuffer[ikey] = trace->getTraceHeader()->doubleValue( vars->hdrId_keys[ikey] );
      }
      timeFunc = vars->table->getFunction( keyValueBuffer, dump );
      delete [] keyValueBuffer;
    }
    else {
      timeFunc = vars->table->getFunction( NULL, dump );
    }
    for( int isamp = firstSampleToProcess; isamp < lastSampleToProcess; isamp++ ) {
      samples[isamp] *= timeFunc->valueAt( (double)isamp * (double)shdr->sampleInt);
    }
  }


  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_scaling_( csParamDef* pdef ) {
  pdef->setModule( "SCALING", "Scale trace data with linear function", "Apply time variant linear scaling function to seismic trace data" );

  pdef->setVersion( 1, 0 );

  pdef->addParam( "filename", "Input ASCII file name, containing pairs of time/scalar values", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "table", "ASCII able", NUM_VALUES_FIXED,
                  "Table format: The ASCII file should contain only numbers, no text. Required are two columns containing time and scalar pairs. Optional: Up to 2 key values. Lines starting with '#' are considered comment lines." );
  pdef->addValue( "", VALTYPE_STRING, "Filename containing table.");

  pdef->addParam( "table_col", "Table columns containing time in [ms] and scalar value", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing time [ms]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing scalar value" );

  pdef->addParam( "table_key", "Key trace header used to match values found in specified table columns", NUM_VALUES_VARIABLE,
                  "Specify the 'table_key' parameter for each key in the table file" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name of key header" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table" );
  pdef->addValue( "yes", VALTYPE_OPTION, "Interpolate based to this key?" );
  pdef->addOption( "yes", "Use this key for interpolation of value" );
  pdef->addOption( "no", "Do not use this key for interpolation", "The input table is expected to contain the exact key values for this trace header" );

  //  pdef->addParam( "table", "ASCII table, containing up to two keys together with pairs of time/scalar values", NUM_VALUES_FIXED );
  //  pdef->addValue( "", VALTYPE_STRING, "Table file name" );

  pdef->addParam( "time", "List of time value [ms]", NUM_VALUES_VARIABLE, "Time knee points at which specified scalars apply" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of time values [ms]..." );

  pdef->addParam( "scalar", "List of scalars", NUM_VALUES_VARIABLE,
      "Scalar values that apply at specified time knee points. In between time knee points, scalar is linearly interpolated." );
  pdef->addValue( "", VALTYPE_NUMBER, "List of scalars..." );

  pdef->addParam( "header", "Trace header name containing scalar value", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "window", "Restrict application of scalar to time window", NUM_VALUES_FIXED, "If specified, scaling will not be applied outside of the specified window" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Start time" );
  pdef->addValue( "99999", VALTYPE_NUMBER, "End time" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_scaling_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_scaling::VariableStruct* vars = reinterpret_cast<mod_scaling::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_scaling_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_scaling::VariableStruct* vars = reinterpret_cast<mod_scaling::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->sampIndex != NULL ) {
    delete [] vars->sampIndex; vars->sampIndex = NULL;
  }
  if( vars->scalar != NULL ) {
    delete [] vars->scalar; vars->scalar = NULL;
  }
  if( vars->bufferScalar != NULL ) {
    delete [] vars->bufferScalar;
    vars->bufferScalar = NULL;
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

extern "C" void _params_mod_scaling_( csParamDef* pdef ) {
  params_mod_scaling_( pdef );
}
extern "C" void _init_mod_scaling_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_scaling_( param, env, writer );
}
extern "C" bool _start_exec_mod_scaling_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_scaling_( env, writer );
}
extern "C" void _exec_mod_scaling_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_scaling_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_scaling_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_scaling_( env, writer );
}

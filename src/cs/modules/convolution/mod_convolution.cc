/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csASCIIFileReader.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: CONVOLUTION
 *
 * @author Bjorn Olofsson
 * @date   2012
 */
namespace mod_convolution {
  struct VariableStruct {
    int asciiFormat;
    cseis_io::ASCIIParam* asciiParam;
    float* bufferTrace;
    float* wavelet;
    int sampleAtZeroTime;
  };
  static int const UNIT_MS = 3;
  static int const UNIT_S  = 4;
}
using mod_convolution::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_convolution_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->asciiFormat     = cseis_io::csASCIIFileReader::FORMAT_COLUMNS;
  vars->bufferTrace     = NULL;
  vars->wavelet   = NULL;
  vars->asciiParam = new cseis_io::ASCIIParam();
  vars->sampleAtZeroTime = 0;

//---------------------------------------------
//

  std::string filename;
  std::string text;
  int unitTime_ascii = mod_convolution::UNIT_MS;

  bool overrideSampleInt = false;
  if( param->exists("override_sample_int") ) {
    param->getString("override_sample_int", &text);
    if( !text.compare("no") ) {
      overrideSampleInt = false;
    }
    else if( !text.compare("yes") ) {
      overrideSampleInt = true;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  param->getString("input_wavelet", &filename);
  if( param->getNumValues("input_wavelet") > 1 ) {
    param->getString("input_wavelet", &text, 1);
    if( !text.compare("ms") ) {
      unitTime_ascii = mod_convolution::UNIT_MS;
    }
    else if( !text.compare("s") ) {
      unitTime_ascii = mod_convolution::UNIT_S;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  if( param->exists("format") ) {
    param->getString("format",&text,0);
    if( !text.compare("columns") ) {
      vars->asciiFormat = cseis_io::csASCIIFileReader::FORMAT_COLUMNS;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  //----------------------------------------------------------------------
  // Read in signature from external ASCII file
  //

  try {
    cseis_io::csASCIIFileReader asciiFileReader( filename, vars->asciiFormat );
    bool success = asciiFileReader.initialize( vars->asciiParam );
    if( !success ) writer->error("Unknown error occurred during initialization of signature input file. Incorrect or unsupported format?");
    success = asciiFileReader.readNextTrace( vars->asciiParam );
    if( !success ) writer->error("Unknown error occurred when reading in samples from signature input file. Incorrect or unsupported format?");
  }
  catch( csException& e ) {
    writer->error("Error occurred when initializing input ASCII file: %s", e.getMessage() );
  }

  //
  //----------------------------------------------------------------------

  float sampleInt_ms = vars->asciiParam->sampleInt;
  float timeZero_ms  = -(float)vars->asciiParam->timeFirstSamp;

  if( unitTime_ascii == mod_convolution::UNIT_S ) {
    sampleInt_ms *= 1000.0;
    timeZero_ms  *= 1000.0;
  }
  if( sampleInt_ms != shdr->sampleInt ) {
    if( !overrideSampleInt ) {
      writer->error("Wavelet input file has different sample interval (=%f ms) than input data (=%f ms). Unsupported case.", sampleInt_ms, shdr->sampleInt);
    }
    else {
      writer->warning("Wavelet input file has different sample interval (=%f ms) than input data (=%f ms). Ignored.", sampleInt_ms, shdr->sampleInt);
    }
  }
  if( vars->asciiParam->numSamples() <= 0 ) {
    writer->error("Could not read in any sample values from input file. Unsupported or incorrect file format?");
  }

  vars->bufferTrace   = new float[shdr->numSamples];
  vars->wavelet = new float[vars->asciiParam->numSamples()];
  for( int isamp = 0; isamp < vars->asciiParam->numSamples(); isamp++ ) {
    vars->wavelet[isamp] = vars->asciiParam->sample(isamp);
  }

  //--------------------------------------------------

  if( param->exists("zero_time") ) {
    param->getFloat("zero_time",&timeZero_ms);
    if( param->getNumValues("zero_time") > 1 ) {
      string unit;
      param->getString( "zero_time", &unit, 1 );
      if( !unit.compare("samples") ) {
        timeZero_ms *= vars->asciiParam->sampleInt;
      }
      else if( !unit.compare("s") ) {
        timeZero_ms *= 1000.0;
      }
      else if( !unit.compare("ms") ) {
        // Nothing
      }
    }
  }
  vars->sampleAtZeroTime = (int)round( timeZero_ms / shdr->sampleInt );

  if( edef->isDebug() ) {
    for( int isamp = 0; isamp < vars->asciiParam->numSamples(); isamp++ ) {
      fprintf( stdout, "%d %f\n", isamp, vars->asciiParam->sample(isamp) );
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_convolution_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  float* samples = trace->getTraceSamples();

  for( int isampOut = 0; isampOut < shdr->numSamples; isampOut++ ) {
    float sum = 0.0;
    int startSamp = std::min( std::max( 0, isampOut + vars->sampleAtZeroTime - (vars->asciiParam->numSamples()-1) ), shdr->numSamples-1 );
    int endSamp   = std::min( isampOut + vars->sampleAtZeroTime, shdr->numSamples-1 );
    //    fprintf(stderr,"%4d  -  %4d  %4d    %d %d\n", isampOut, startSamp, endSamp, vars->asciiParam->numSamples(), shdr->numSamples);
    for( int isampIn = startSamp; isampIn <= endSamp; isampIn++ ) {
      int index = isampOut + vars->sampleAtZeroTime - isampIn;
      sum += samples[isampIn] * vars->wavelet[index];
    }
    vars->bufferTrace[isampOut] = sum;
  }
  memcpy( samples, vars->bufferTrace, sizeof(float)*shdr->numSamples );

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_convolution_( csParamDef* pdef ) {
  pdef->setModule( "CONVOLUTION", "Convolve trace with input (filter) wavelet" );

  pdef->addParam( "input_wavelet", "Name of file containing input (filter) wavelet", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit of time values found in input file, if any" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );

  pdef->addParam( "format", "Input wavelet ASCII file format", NUM_VALUES_FIXED);
  pdef->addValue( "columns", VALTYPE_OPTION );
  pdef->addOption( "columns", "Simple file format with 2 or 3 columns:  Time[ms]  Amplitude" );

  pdef->addParam( "zero_time", "Zero time in input wavelet", NUM_VALUES_VARIABLE,
                  "Specifying this parameter overrides the zero time that is found in the ASCII signature file" );
  pdef->addValue( "0", VALTYPE_NUMBER );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );
  pdef->addOption( "samples", "Samples" );

  pdef->addParam( "override_sample_int", "Override sample interval?", NUM_VALUES_FIXED);
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "" );
  pdef->addOption( "yes", "Ignore sample interval of input wavelet. Assume it is the same as the input data" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_convolution_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_convolution::VariableStruct* vars = reinterpret_cast<mod_convolution::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_convolution_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_convolution::VariableStruct* vars = reinterpret_cast<mod_convolution::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->asciiParam != NULL ) {
    delete vars->asciiParam;
    vars->asciiParam = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_convolution_( csParamDef* pdef ) {
  params_mod_convolution_( pdef );
}
extern "C" void _init_mod_convolution_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_convolution_( param, env, writer );
}
extern "C" bool _start_exec_mod_convolution_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_convolution_( env, writer );
}
extern "C" void _exec_mod_convolution_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_convolution_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_convolution_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_convolution_( env, writer );
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csStandardHeaders.h"
#include "csASCIIFileReader.h"
#include <cstring>
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: INPUT_ASCII
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_input_ascii {
  struct VariableStruct {
    int inputFormat;

    bool atEOF;
    bool isZeroPhase;
    int traceCounter;

    int hdrId_time_samp1;
    int hdrId_time_samp1_us;
    int hdrId_sou_z;
    int hdrId_binx;
    int hdrId_biny;
    int hdrId_trcno;
    int hdrId_chan;

    double zmap_noValue;
    int zmap_numValues;
    int zmap_numBlocks;
    double zmap_x1;
    double zmap_y1;
    double zmap_x2;
    double zmap_y2;

    int unit;

    cseis_io::csASCIIFileReader* asciiFileReader;
    cseis_io::ASCIIParam* asciiParam;
    float phaseScalar;
  };
  static int const UNIT_MS = 3;
  static int const UNIT_S  = 4;
  static int const UNIT_HZ = 5;
}
using mod_input_ascii::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_ascii_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->inputFormat       = cseis_io::csASCIIFileReader::FORMAT_COLUMNS;
  vars->atEOF        = false;
  vars->asciiParam   = NULL;
  vars->asciiFileReader = NULL;
  vars->traceCounter = 0;
  vars->unit = mod_input_ascii::UNIT_MS;
  vars->zmap_noValue = 0;
  vars->zmap_numValues = 0;
  vars->zmap_numBlocks = 0;
  vars->zmap_x1 = 0;
  vars->zmap_y1 = 0;
  vars->zmap_x2 = 0;
  vars->zmap_y2 = 0;

  vars->hdrId_binx  = -1;
  vars->hdrId_biny  = -1;
  vars->hdrId_trcno = -1;
  vars->hdrId_chan  = -1;
  vars->hdrId_time_samp1 = -1;
  vars->hdrId_time_samp1_us = -1;
  vars->hdrId_sou_z = -1;
  vars->isZeroPhase = true;
  vars->phaseScalar = 1.0f;

  string filename;
  param->getString( "filename", &filename );

  //----------------------------------------------------

  if( param->exists("format") ) {
    string text;
    param->getString( "format", &text );
    if( !text.compare("signature") ) {
      vars->inputFormat = cseis_io::csASCIIFileReader::FORMAT_NUCLEUS_SIGNATURE;
    }
    else if( !text.compare("columns") ) {
      vars->inputFormat = cseis_io::csASCIIFileReader::FORMAT_COLUMNS;
    }
    else if( !text.compare("zmap") ) {
      vars->inputFormat = cseis_io::csASCIIFileReader::FORMAT_ZMAP;
    }
    else {
      writer->error("Option not recognised: %s", text.c_str());
    }
  }
  int colIndexTime  = 0;
  int colIndexValue = 1;
  int colIndexTrace = -1;
  if( vars->inputFormat == cseis_io::csASCIIFileReader::FORMAT_COLUMNS ) {
    if( param->exists("format") && param->getNumValues("format") > 1 ) {
      param->getInt( "format", &colIndexTime, 1 );
      colIndexTime  -= 1;
      if( param->getNumValues("format") > 2 ) {
        param->getInt( "format", &colIndexValue, 2 );
        colIndexValue -= 1;
        if( param->getNumValues("format") > 3 ) {
          param->getInt( "format", &colIndexTrace, 3 );
          colIndexTrace -= 1;
        }
      }
    }
  }

  if( param->exists("phase") ) {
    string text;
    param->getString( "phase", &text );
    if( !text.compare("zero") ) {
      vars->isZeroPhase = true;
    }
    else if( !text.compare("linear") ) {
      vars->isZeroPhase = false;
      if( param->getNumValues("phase") > 1 ) {
        param->getFloat( "phase", &vars->phaseScalar, 1 );
      }
    }
    else {
      writer->error("Option not recognised: %s", text.c_str());
    }
  }

  if( param->exists("unit") ) {
    string text;
    param->getString( "unit", &text );
    if( !text.compare("ms") ) {
      vars->unit = cseis_io::csASCIIFileReader::FORMAT_NUCLEUS_SIGNATURE;
    }
    else if( !text.compare("s") ) {
      vars->unit = mod_input_ascii::UNIT_S;
    }
    else if( !text.compare("hz") ) {
      vars->unit = mod_input_ascii::UNIT_HZ;
    }
    else {
      writer->error("Option not recognised: %s", text.c_str());
    }
  }

  bool isOverrideSampleInt = false;
  if( param->exists("sample_int") ) {
    param->getFloat( "sample_int", &shdr->sampleInt, 0 );
    isOverrideSampleInt = true;
  }

  bool isOverrideNumSamples = false;
  if( param->exists("nsamples") ) {
    param->getInt( "nsamples", &shdr->numSamples, 0 );
    isOverrideNumSamples = true;
  }

  //----------------------------------------------------------------------
  //

  try {
    vars->asciiParam = new cseis_io::ASCIIParam();
    vars->asciiFileReader = new cseis_io::csASCIIFileReader( filename, vars->inputFormat );
    bool success = false;
    if( vars->inputFormat != cseis_io::csASCIIFileReader::FORMAT_ZMAP ) {
      success = vars->asciiFileReader->initialize( vars->asciiParam, colIndexTime, colIndexValue, colIndexTrace );
    }
    else {
      success = vars->asciiFileReader->initializeZMap( vars->asciiParam,
                                                       vars->zmap_numBlocks,
                                                       vars->zmap_x1,
                                                       vars->zmap_y1,
                                                       vars->zmap_x2,
                                                       vars->zmap_y2 );
    }
    if( !success ) writer->error("Unknown error occurred during initialization of signature input file. Incorrect or unsupported format?");
  }
  catch( csException& e ) {
    writer->error("Error occurred when initializing input ASCII file: %s", e.getMessage() );
  }
  //
  //----------------------------------------------------------------------


  if( vars->asciiParam->sampleInt <= 0 ) {
    writer->error("Inconsistent sample interval (=%f).", vars->asciiParam->sampleInt);
  }
  if( vars->asciiParam->numSamples() <= 0 ) {
    writer->error("Inconsistent number of samples (=%d).", vars->asciiParam->numSamples());
  }

  if( vars->unit == mod_input_ascii::UNIT_HZ ) {
    if( !isOverrideSampleInt ) {
      shdr->sampleInt  = (float)( vars->asciiParam->timeLastSamp - vars->asciiParam->timeFirstSamp ) / (float)(vars->asciiParam->numSamples()-1);
    }
    shdr->domain       = DOMAIN_FX;
    shdr->fftDataType  = FX_AMP_PHASE;
    // Split trace into amplitude & phase. Read in amplitude, set phase to zero.
    // Input data must have for example 1024+1 samples. 
    if( isOverrideNumSamples ) {
      writer->error("Unsupported option for frequency input data: Number of samples overriden by user. Do not specify input parameter 'nsamples'.");
    }
    shdr->numSamples   = 2*vars->asciiParam->numSamples();
    shdr->numSamplesXT = shdr->numSamples-2; // Remove 2 for Nyquist frequency
    int numFFTSamples  = shdr->numSamplesXT;
    shdr->sampleIntXT  = round( 1000.0 / (float)( numFFTSamples * shdr->sampleInt ) );
  }
  else {
    shdr->domain = DOMAIN_XT;
    if( !isOverrideNumSamples ) {
      shdr->numSamples = vars->asciiParam->numSamples();
    }
    if( !isOverrideSampleInt ) {
      if( vars->unit == mod_input_ascii::UNIT_MS ) {
        shdr->sampleInt = vars->asciiParam->sampleInt;
      }
      else {
        shdr->sampleInt = 1000.0 * vars->asciiParam->sampleInt;
      }
    }
  }

  if( vars->inputFormat == cseis_io::csASCIIFileReader::FORMAT_NUCLEUS_SIGNATURE ) {
    if( !hdef->headerExists(cseis_geolib::HDR_SOU_Z.name) ) {
      hdef->addStandardHeader(cseis_geolib::HDR_SOU_Z.name);
    }
    vars->hdrId_sou_z = hdef->headerIndex(cseis_geolib::HDR_SOU_Z.name);
  }
  //-----------------------------------------------------------------------------
  else if( vars->inputFormat == cseis_io::csASCIIFileReader::FORMAT_COLUMNS ) {
  }
  //--------------------------------------------------------------------------------
  // ZMAP data
  //
  else if( vars->inputFormat == cseis_io::csASCIIFileReader::FORMAT_ZMAP ) {
    vars->hdrId_biny = hdef->addStandardHeader(cseis_geolib::HDR_BIN_Y.name);
  }

  vars->hdrId_chan          = hdef->addStandardHeader(cseis_geolib::HDR_CHAN.name);
  vars->hdrId_trcno         = hdef->addStandardHeader(cseis_geolib::HDR_TRCNO.name);
  vars->hdrId_time_samp1    = hdef->headerIndex(cseis_geolib::HDR_TIME_SAMP1.name);
  vars->hdrId_time_samp1_us = hdef->headerIndex(cseis_geolib::HDR_TIME_SAMP1_US.name);

  if( edef->isDebug() ) {
    writer->line("Input file parameters:");
    writer->line("  Time/Freq of first sample: %f", vars->asciiParam->timeFirstSamp);
    writer->line("  Time/Freq of last sample:  %f", vars->asciiParam->timeLastSamp);
    writer->line("  Source depth:              %f", vars->asciiParam->srcDepth);
    writer->line("  Sample interval:           %f", vars->asciiParam->sampleInt);
    if( isOverrideSampleInt ) {
      writer->line("    ..overriden by user specified sample interval:       %f", shdr->sampleInt);
    }
    writer->line("  Number of samples:         %d", vars->asciiParam->numSamples());
  }

  writer->line("Input file            :   %s", filename.c_str() );
  writer->line("Output sample interval:   %f %s", shdr->sampleInt, (shdr->domain == DOMAIN_FX) ? "Hz" : "ms");
  writer->line("Number of output samples: %d", shdr->numSamples );

  vars->traceCounter = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_input_ascii_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const*  shdr = env->superHeader;
//  csTraceHeaderDef const*  hdef = env->headerDef;

  csTrace* trace = traceGather->trace(0);

  if( vars->atEOF ) {
    traceGather->freeAllTraces();
    return;
  }

  vars->traceCounter += 1;

  // Read in next trace
  try {
    bool success = vars->asciiFileReader->readNextTrace( vars->asciiParam );
    if( !success ) { // No further trace could be read in
      vars->atEOF = true;
      traceGather->freeAllTraces();
      return;
    }
  }
  catch( csException& e ) {
    writer->error("Error occurred when reading from input ASCII file: %s", e.getMessage() );
  }

  //--------------------------------------------------------------------------------
  //
  if( vars->inputFormat == cseis_io::csASCIIFileReader::FORMAT_NUCLEUS_SIGNATURE ) {
    // Only one trace to read in
    vars->atEOF = true;
  }
  //--------------------------------------------------------------------------------
  //
  else if( vars->inputFormat == cseis_io::csASCIIFileReader::FORMAT_COLUMNS ) {
    int numSamples = vars->asciiParam->numSamples();
    if( vars->unit != mod_input_ascii::UNIT_HZ ) {
      if( numSamples != shdr->numSamples ) {
        writer->error("Incorrect number of samples for trace %d: %d != %d", vars->traceCounter, numSamples, shdr->numSamples);
      }
      if( edef->isDebug() ) writer->line("Number of samples for trace %d: %d", vars->traceCounter, shdr->numSamples);
    }
    else { // FX DOMAIN
      if( 2*numSamples != shdr->numSamples ) {
        writer->error("Incorrect number of samples for trace %d: %d != %d", vars->traceCounter, numSamples, (int)(0.5*shdr->numSamples));
      }
    }
  }
  //--------------------------------------------------------------------------------
  //
  else if( vars->inputFormat == cseis_io::csASCIIFileReader::FORMAT_ZMAP ) {
    trace->getTraceHeader()->setDoubleValue(vars->hdrId_biny, vars->traceCounter * (vars->zmap_y1-vars->zmap_y2)/(vars->zmap_numBlocks-1) + vars->zmap_y2 );
  }

  // Copy sample values to output trace
  float* samples = trace->getTraceSamples();
  if( vars->unit != mod_input_ascii::UNIT_HZ ) {
    memcpy( samples, vars->asciiParam->getSamples(), vars->asciiParam->numSamples()*sizeof(float) );
    for( int isamp = vars->asciiParam->numSamples(); isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] = 0.0;
    }
  }
  else {
    memcpy( samples, vars->asciiParam->getSamples(), (shdr->numSamples/2)*sizeof(float) ); // Copy amplitudes up to Nyquist
    if( vars->isZeroPhase ) {
      for( int i = shdr->numSamples/2; i < shdr->numSamples; i++ ) {
        samples[i] = 0.0;
      }
    }
    else { // Set linear phase shift
      float scalar = vars->phaseScalar;
      for( int i = shdr->numSamples/2; i < shdr->numSamples; i++ ) {
        float freqHz = (i - shdr->numSamples/2) * shdr->sampleInt;
        float phase = freqHz * M_PI * scalar;
        samples[i] = phase;
      }
    }
  }

  // Set trace headers
  csTraceHeader* trcHdr = trace->getTraceHeader();
  trcHdr->setIntValue( vars->hdrId_trcno, vars->traceCounter );
  trcHdr->setIntValue( vars->hdrId_chan, vars->asciiParam->traceNumber );

  double scalar = ( vars->unit == mod_input_ascii::UNIT_MS ) ? 0.001 : 1.0;
  int timeFirstSamp_s  = (int)( vars->asciiParam->timeFirstSamp * scalar );
  int timeFirstSamp_us = (int)( (vars->asciiParam->timeFirstSamp * scalar - (double)timeFirstSamp_s)*1e6 );
  trcHdr->setIntValue( vars->hdrId_time_samp1, timeFirstSamp_s );
  trcHdr->setIntValue( vars->hdrId_time_samp1_us, timeFirstSamp_us );
  if( vars->hdrId_sou_z >= 0 ) {
    trcHdr->setFloatValue(vars->hdrId_sou_z, vars->asciiParam->srcDepth);
  }

  //------------------------------------------------------------------------------------

  return;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_ascii_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_ASCII", "Input ASCII file", "Reads seismic data from ASCII file" );
  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "format", "Input ASCII file format", NUM_VALUES_VARIABLE);
  pdef->addValue( "columns", VALTYPE_OPTION );
  pdef->addOption( "signature", "Read in source signature from Nucleus ASCII file", "Supported ASCII formats are: Far-field signature, Nucleus+ wavelet file" );
  pdef->addOption( "columns", "Simple file format with 2 or 3 columns:  Time[ms]  Amplitude  (Trace number)",
                   "Trace number column is optional" );
  pdef->addOption( "zmap", "Input grid data in ZMAP format");
  //  pdef->addOption( "nucleus_wav", "Read in source signature from Nucleus+ wavelet (*.wav) ASCII file" );
  //  pdef->addOption( "spikogram", "Read in pairs of time/amplitude values");
  pdef->addValue( "1", VALTYPE_NUMBER, "For 'column' format: Column number specifying time (or other sample unit)" );
  pdef->addValue( "2", VALTYPE_NUMBER, "For 'column' format: Column number specifying sample value" );
  pdef->addValue( "-1", VALTYPE_NUMBER, "Optional: Column number specifying trace number. -1: Trace number not applicable" );

  pdef->addParam( "sample_int", "Sample interval [ms]", NUM_VALUES_FIXED,
                  "Set/override sample interval found in input file with this value");
  pdef->addValue( "", VALTYPE_NUMBER, "Sample interval [ms]" );

  pdef->addParam( "nsamples", "Number of samples in output trace", NUM_VALUES_FIXED, "Override number of samples found in input file");
  pdef->addValue( "", VALTYPE_NUMBER, "Number of samples" );

  pdef->addParam( "unit", "Unit of vertical dimension found in input file", NUM_VALUES_FIXED,
                  "Only required when input data is amplitude spectrum (specify 'hz')");
  pdef->addValue( "ms", VALTYPE_OPTION );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "seconds" );
  pdef->addOption( "Hz", "Hertz", "Input data is amplitude spectrum, phase = 0" );

  pdef->addParam( "phase", "Phase behaviour in case input is amplitude spectrum", NUM_VALUES_VARIABLE);
  pdef->addValue( "zero", VALTYPE_OPTION );
  pdef->addOption( "linear", "Set to linear phase" );
  pdef->addOption( "zero", "Set to zero phase" );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Linear phase shift (scalar)" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_input_ascii_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_input_ascii::VariableStruct* vars = reinterpret_cast<mod_input_ascii::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_input_ascii_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_input_ascii::VariableStruct* vars = reinterpret_cast<mod_input_ascii::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->asciiFileReader != NULL ) {
    delete vars->asciiFileReader;
    vars->asciiFileReader = NULL;
  }
  if( vars->asciiParam != NULL ) {
    delete vars->asciiParam;
    vars->asciiParam = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_input_ascii_( csParamDef* pdef ) {
  params_mod_input_ascii_( pdef );
}
extern "C" void _init_mod_input_ascii_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_ascii_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_ascii_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_ascii_( env, writer );
}
extern "C" void _exec_mod_input_ascii_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_ascii_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_ascii_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_ascii_( env, writer );
}

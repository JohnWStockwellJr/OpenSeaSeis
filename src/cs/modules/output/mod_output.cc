/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csIODefines.h"
#include "csSeismicWriter.h"
#include "csTimer.h"
#include "csFileUtils.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: OUTOUT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_output {
  struct VariableStruct {
    std::string filename;
    cseis_system::csSeismicWriter* writer;
    int nTracesOut;
    bool isFirstCall;
    int numTracesBuffer;
    int sampleByteSize; //, doOverwrite
  };
}
using namespace mod_output;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_output_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  //  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->filename   = "";
  vars->writer     = NULL;
  vars->nTracesOut = 0;
  vars->numTracesBuffer = 20;
  vars->sampleByteSize = 4;
  vars->isFirstCall = true;

  bool doOverwrite = true;
  if( param->exists("overwrite") ) {
    std::string text;
    param->getString("overwrite", &text);
    if( !text.compare("yes") ) {
      doOverwrite = true;
    }
    else if( !text.compare("no") ) { 
      doOverwrite = false;
    }
    else {
      writer->line("Unknown option for user parameter overwrite: '%s'", text.c_str());
      env->addError();
    }
  }

  if( param->exists("compress") ) {
    std::string text;
    param->getString("compress", &text);
    if( !text.compare("no") || !text.compare("32bit") ) {
      vars->sampleByteSize = 4;
    }
    else if( !text.compare("16bit") ) { 
      vars->sampleByteSize = 2;
    }
    else if( !text.compare("8bit") ) { 
      vars->sampleByteSize = 1;
    }
    else {
      writer->line("Unknown option for user parameter compress: '%s'", text.c_str());
      env->addError();
    }
  }

  param->getString( "filename", &vars->filename );
  int length = (int)vars->filename.length();
  if( length < 6 || (vars->filename.substr( length-6, 6 ).compare(".cseis")) ) {
    writer->warning("Output file name does not have the standard SeaSeis extension '.cseis': '%s'\n", vars->filename.c_str());
    if( vars->filename.substr( length-6, 6 ).compare(".oseis") ) {
      //    if( filename.find('.') != string::npos ) {
      writer->error("SeaSeis standard file name extension '.cseis' was omitted or spelled wrong. Please add extension '.cseis' to output file name and re-run.");
    }
  }

  if( param->exists( "ntraces_buffer" ) ) {
    param->getInt( "ntraces_buffer", &vars->numTracesBuffer );
    if( vars->numTracesBuffer <= 0 || vars->numTracesBuffer > 999999 ) {
      writer->warning("Number of buffered traces out of range (=%d). Changed to 1.", vars->numTracesBuffer);
      vars->numTracesBuffer = 1;
    }
  }

  if( !doOverwrite && csFileUtils::fileExists( vars->filename ) ) {
    writer->error("File %s already exists but user parameter set to 'overwrite no'.", vars->filename.c_str() );
  }
  if( !csFileUtils::createDoNotOverwrite( vars->filename ) ) {
    writer->error("Cannot open SeaSeis output file '%s'", vars->filename.c_str() );
  }

  writer->line("");
  writer->line("  File name:             %s", vars->filename.c_str());
  writer->line("  Sample interval [ms]:  %f", shdr->sampleInt);
  writer->line("  Number of samples:     %d", shdr->numSamples);
  writer->line("");

  vars->nTracesOut = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_output_(
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

  csTrace* trace = traceGather->trace(0);


  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    try {
      vars->writer = new csSeismicWriter( vars->filename, vars->numTracesBuffer, vars->sampleByteSize, true );
    }
    catch( csException& exc ) {
      writer->error("Error occurred when opening SeaSeis file. System message:\n%s", exc.getMessage() );
    }
  }


  if( edef->isDebug() ) writer->line("Output SeaSeis trace #%d", vars->nTracesOut+1);

  float* samples = trace->getTraceSamples();
  char const* hdrValueBlock   = trace->getTraceHeader()->getTraceHeaderValueBlock();

  try {
    bool success = true;
    if( vars->nTracesOut == 0 ) {
      success = vars->writer->writeFileHeader( shdr, hdef );
    }
    if( success ) success = vars->writer->writeTrace( samples, hdrValueBlock );
    if( !success ) {
      writer->error("Unknown error occurred when writing to SeaSeis file.\n" );
    }
  }
  catch( csException& exc ) {
    writer->error("Error occurred when writing to SeaSeis file. System message:\n%s", exc.getMessage() );
  }

  vars->nTracesOut += 1;

  return;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_output_( csParamDef* pdef ) {
  pdef->setModule( "OUTPUT", "Output SeaSeis data", "Writes SeaSeis formatted data to disk file, file extension '.cseis'" );

  pdef->addParam( "filename", "Output file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );
  
  pdef->addParam( "ntraces_buffer", "Number of traces to buffer before write operation", NUM_VALUES_FIXED,
                  "Writing a large number of traces at once enhances performance, but requires more memory" );
  pdef->addValue( "20", VALTYPE_NUMBER, "Number of traces to buffer before writing" );

  pdef->addParam( "overwrite", "Overwrite exisiting file?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Overwrite file if it already exists" );
  pdef->addOption( "no", "Do not overwrite file if it already exists");


  pdef->addParam( "compress", "Compress data before output?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "No compression. Save data samples as 32bit floating point" );
  pdef->addOption( "32bit", "No compression. Same as option 'no'");
  pdef->addOption( "16bit", "Compress data samples to 16bit");
  pdef->addOption( "8bit", "Compress data samples to 8bit");
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_output_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_output::VariableStruct* vars = reinterpret_cast<mod_output::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_output_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_output::VariableStruct* vars = reinterpret_cast<mod_output::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->writer != NULL ) {
    delete vars->writer;
    vars->writer = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_output_( csParamDef* pdef ) {
  params_mod_output_( pdef );
}
extern "C" void _init_mod_output_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_output_( param, env, writer );
}
extern "C" bool _start_exec_mod_output_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_output_( env, writer );
}
extern "C" void _exec_mod_output_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_output_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_output_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_output_( env, writer );
}

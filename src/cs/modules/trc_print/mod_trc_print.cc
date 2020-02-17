/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFileUtils.h"
#include <cmath>
#include <string>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TRC_PRINT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_trc_print {
  struct VariableStruct {
    int option;
    cseis_geolib::type_t hdrType_user;
    int hdrId_user;
    bool print_hdr;
    std::string format_time;
    std::string format_value;
    std::string printFormat;
    FILE* fout;
    std::string filename;
    int traceCounter;
    bool addBlankLine;
    float timeSamp1;
    int hdrId_timeSamp1;
    bool isFirstCall;
  };
}
using mod_trc_print::VariableStruct;

const int OPTION_LIST  = 10;
const int OPTION_TABLE = 11;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_trc_print_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->option       = 0;
  vars->hdrId_user   = -1;
  vars->hdrType_user = TYPE_UNKNOWN;
  vars->print_hdr    = false;
  vars->fout         = NULL;
  vars->format_time  = "%10.2f";
  vars->format_value = "%16.10e";
  vars->printFormat = "";
  vars->traceCounter = 0;
  vars->addBlankLine = true;
  vars->hdrId_timeSamp1 = -1;
  vars->timeSamp1    = 0;
  vars->isFirstCall = true;

  if( param->exists("filename") ) {
    param->getString("filename",&vars->filename);
    if( !csFileUtils::createDoNotOverwrite( vars->filename ) ) {
      writer->error("Unable to open output file %s. Wrong path name..?", vars->filename.c_str() );
    }
  }
  else {
    vars->fout = writer->getFile();
  }

  if( param->exists("format") ) {
    param->getString("format", &vars->format_time,0);
    param->getString("format", &vars->format_value,1);
  }
  if( param->exists("header") ) {
    vars->print_hdr  = true;
    string name;
    param->getString("header", &name);
    vars->hdrId_user   = env->headerDef->headerIndex(name);
    vars->hdrType_user = env->headerDef->headerType(name);
  }
  if( param->exists("blank_lines") ) {
    std::string text;
    param->getString("blank_lines", &text);
    if( !text.compare( "yes" ) ) {
      vars->addBlankLine = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->addBlankLine = false;
    }
    else {
      writer->line("Unknown option:: '%s'.", text.c_str());
      env->addError();
    }
  }
  if( param->exists("time_samp1") ) {
    csFlexNumber number;
    std::string text;
    param->getString( "time_samp1", &text );
    if( !number.convertToNumber( text ) ) {
      if( !hdef->headerExists(text) ) {
        writer->error("Specified trace header name containing time of first sample does not exist: '%s'", text.c_str());
      }
      vars->hdrId_timeSamp1 = hdef->headerIndex(text);
    }
    else { // User specified a constant value
      vars->timeSamp1 = number.floatValue();
      writer->line("Constant time of first sample = %.2f", vars->timeSamp1);
    }
  }
  vars->printFormat = "%d  " + vars->format_time + "  " + vars->format_value + "\n";
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_trc_print_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    if( vars->fout == NULL ) {
      vars->fout = fopen(vars->filename.c_str(),"w");
      if( vars->fout == NULL ) writer->error("Unable to open output file %s. Wrong path name..?", vars->filename.c_str() );
    }
  }

  //  writer->line("TRC_PRINT: numSamples: %d, superheader: %d",trace->numSamples(), shdr->numSamples);
  float* samples = trace->getTraceSamples();
  vars->traceCounter += 1;

  if( vars->addBlankLine && vars->traceCounter != 1 ) fprintf(vars->fout,"\n");

  if( vars->hdrId_timeSamp1 >= 0 ) {
    vars->timeSamp1 = trace->getTraceHeader()->floatValue( vars->hdrId_timeSamp1 );
  }

  if( vars->print_hdr ) {
    if( vars->hdrType_user == TYPE_INT ) {
      string printFormat = "%d  %8d " + vars->format_time + "  " + vars->format_value + "\n";
      int value = trace->getTraceHeader()->intValue( vars->hdrId_user );
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        float time = isamp*shdr->sampleInt + vars->timeSamp1;
        fprintf(vars->fout,printFormat.c_str(), isamp, value, time, samples[isamp] );
        //        fprintf(vars->fout,"%d  %8d %10.2f  %16.10e\n", isamp, value, time, samples[isamp] );
      }
    }
    else if( vars->hdrType_user == TYPE_FLOAT || vars->hdrType_user == TYPE_DOUBLE ) {
      string printFormat = "%d  %14.6f " + vars->format_time + "  " + vars->format_value + "\n";
      double value = trace->getTraceHeader()->doubleValue( vars->hdrId_user );
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        float time = isamp*shdr->sampleInt + vars->timeSamp1;
        fprintf(vars->fout,printFormat.c_str(), isamp, value, time, samples[isamp] );
        //        fprintf(vars->fout,"%d  %14.6f %10.2f  %16.10e\n", isamp, value, time, samples[isamp] );
      }
    }
    else if( vars->hdrType_user == TYPE_STRING ) {
      string printFormat = "%d  %s " + vars->format_time + "  " + vars->format_value + "\n";
      string value = trace->getTraceHeader()->stringValue( vars->hdrId_user );
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        float time = isamp*shdr->sampleInt + vars->timeSamp1;
        fprintf(vars->fout,printFormat.c_str(), isamp, value.c_str(), time, samples[isamp] );
        //        fprintf(vars->fout,"%d  %s %10.2f  %16.10e\n", isamp, value.c_str(), time, samples[isamp] );
      }
    }
  }
  else {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      float time = isamp*shdr->sampleInt + vars->timeSamp1;
      fprintf(vars->fout,vars->printFormat.c_str(), isamp, time, samples[isamp] );
      //      fprintf(vars->fout,"%d  %10.2f  %16.10e\n", isamp, time, samples[isamp] );
    }
  }
  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_trc_print_( csParamDef* pdef ) {
  pdef->setModule( "TRC_PRINT", "Print trace samples", "Print trace samples to log file" );

  pdef->addParam( "filename", "Output file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );

  pdef->addParam( "header", "Trace header to print", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "format", "Specify floating point format (C style)", NUM_VALUES_FIXED );
  pdef->addValue( "%10.2f", VALTYPE_STRING, "Format for sample time/frequency" );
  pdef->addValue( "%16.10e", VALTYPE_STRING, "Format for sample value" );

  pdef->addParam( "blank_lines", "Add blank lines between traces?", NUM_VALUES_FIXED, "Adds blank lines in output file between all traces" );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Add blank lines between traces" );
  pdef->addOption( "no", "Do NOT add blank lines between traces" );

  pdef->addParam( "time_samp1", "Time of first sample (in units of data)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_HEADER_NUMBER, "Time of first sample: Value or name of trace header containing the value");
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_trc_print_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_trc_print::VariableStruct* vars = reinterpret_cast<mod_trc_print::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_trc_print_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_trc_print::VariableStruct* vars = reinterpret_cast<mod_trc_print::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->fout != NULL ) {
    fclose( vars->fout );
    vars->fout = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_trc_print_( csParamDef* pdef ) {
  params_mod_trc_print_( pdef );
}
extern "C" void _init_mod_trc_print_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_trc_print_( param, env, writer );
}
extern "C" bool _start_exec_mod_trc_print_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_trc_print_( env, writer );
}
extern "C" void _exec_mod_trc_print_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_trc_print_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_trc_print_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_trc_print_( env, writer );
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFFT.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TRC_INTERPOL
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_trc_interpol {
  struct VariableStruct {
    int method;
    int mode;
    int hdrID_in;
    int hdrID_out;
    double value;
    cseis_geolib::csFFT* fftTool1;
    cseis_geolib::csFFT* fftTool2;
    float* buffer;
  };
  static int const METHOD_SIMPLE_AVERAGE = 1;
  static int const METHOD_FX_AVERAGE     = 2;
  
  static int const MODE_ENSEMBLE = 11;
  static int const MODE_FIXED    = 12;
  static int const MODE_HEADER   = 13;
}

using namespace mod_trc_interpol;

//*******************************************************************
//
//  Init phase
//
//*******************************************************************

void init_mod_trc_interpol_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->method = METHOD_SIMPLE_AVERAGE;
  vars->mode = MODE_FIXED;
  vars->hdrID_out = -1;
  vars->hdrID_in = -1;
  vars->value = 0.0;
  vars->fftTool1 = NULL;
  vars->fftTool2 = NULL;
  vars->buffer   = NULL;
  
  std::string text;
  if( param->exists("mode") ) {
    param->getString("mode",&text);
    if( !text.compare("fixed") ) {
      vars->mode = MODE_FIXED;
    }
    else if( !text.compare("ensemble") ) {
      vars->mode = MODE_ENSEMBLE;
    }
    else if( !text.compare("header") ) {
      vars->mode = MODE_HEADER;
      param->getString("header",&text);
      param->getDouble("header",&vars->value,1);
      vars->hdrID_in = hdef->headerIndex( text.c_str() );
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  if( param->exists("method") ) {
    param->getString("method",&text);
    if( !text.compare("simple_average") ) {
      vars->method = METHOD_SIMPLE_AVERAGE;
    }
    else if( !text.compare("fx_average") ) {
      vars->method = METHOD_FX_AVERAGE;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  if( vars->mode == MODE_ENSEMBLE ) {
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }
  else if( vars->mode == MODE_FIXED ) {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 2 );
  }
  else {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 3 );
  }

  if( !hdef->headerExists("interpolated") ) {
    hdef->addHeader( TYPE_INT, "interpolated", "1: Trace is interpolated" );
  }
  vars->hdrID_out = hdef->headerIndex("interpolated");

  if( vars->method == METHOD_FX_AVERAGE ) {
    vars->fftTool1 = new cseis_geolib::csFFT( shdr->numSamples );
    vars->fftTool2 = new cseis_geolib::csFFT( shdr->numSamples );
    vars->buffer = new float[shdr->numSamples];
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************

void exec_mod_trc_interpol_(
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


  if( edef->isDebug() ) writer->line("Number of input traces: %d, num trc to keep: %d", traceGather->numTraces(), *numTrcToKeep);

  int numTraces = traceGather->numTraces();
  if( numTraces < 2 ) {  // This must be the last trace of the data set. Simply return, do nothing
    return;
  }

  int nSamples = shdr->numSamples;
  if( vars->mode == MODE_FIXED ) {
    // 1) Create new trace in between the two traces
    traceGather->createTraces( 1, 1, hdef, nSamples );

    // 2) Interpolate trace samples
    float* samples1 = traceGather->trace(0)->getTraceSamples();
    float* samples2 = traceGather->trace(2)->getTraceSamples();
    float* samplesNew = traceGather->trace(1)->getTraceSamples();
  
    if( vars->method == METHOD_SIMPLE_AVERAGE ) {
      for( int isamp = 0; isamp < nSamples; isamp++ ) {
        samplesNew[isamp] = (samples1[isamp] + samples2[isamp])*0.5;
      }
    }

    // 2) Interpolate trace header values
    csTraceHeader* trcHdr1 = traceGather->trace(0)->getTraceHeader();
    csTraceHeader* trcHdr2 = traceGather->trace(2)->getTraceHeader();
    csTraceHeader* trcHdrNew = traceGather->trace(1)->getTraceHeader();

    for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
      char type = hdef->headerType(ihdr);
      switch( type ) {
      case TYPE_INT:
        trcHdrNew->setIntValue( ihdr, (trcHdr1->intValue(ihdr)+trcHdr2->intValue(ihdr))/2 );
        break;
      case TYPE_FLOAT:
        trcHdrNew->setFloatValue( ihdr, (trcHdr1->floatValue(ihdr)+trcHdr2->floatValue(ihdr))/2 );
        break;
      case TYPE_DOUBLE:
        trcHdrNew->setDoubleValue( ihdr, (trcHdr1->doubleValue(ihdr)+trcHdr2->doubleValue(ihdr))/2 );
        break;
      case TYPE_INT64:
        trcHdrNew->setInt64Value( ihdr, (trcHdr1->int64Value(ihdr)+trcHdr2->int64Value(ihdr))/2 );
        break;
      case TYPE_STRING:
        trcHdrNew->setStringValue( ihdr, trcHdr1->stringValue(ihdr) );
        break;
      default:
        writer->error("mod_trc_interpol: Unknown trace header type, code: %d", type);
      }
    }
    trcHdrNew->setIntValue( vars->hdrID_out, 1 );

    *numTrcToKeep = 1;  // Keep last trace for next pair of traces to interpolate in between
  } // END FIXED
  else if( vars->mode == MODE_ENSEMBLE ) {
    for( int itrc = numTraces-1; itrc > 0; itrc-- ) {
      int trcIndex1 = itrc-1;
      int trcIndex2 = itrc+1;
      // 1) Create new trace in between the two traces
      traceGather->createTrace( itrc, hdef, nSamples );

      // 2) Interpolate trace samples
      float* samples1 = traceGather->trace(trcIndex1)->getTraceSamples();
      float* samples2 = traceGather->trace(trcIndex2)->getTraceSamples();
      float* samplesNew = traceGather->trace(itrc)->getTraceSamples();
  
      if( vars->method == METHOD_SIMPLE_AVERAGE ) {
        for( int isamp = 0; isamp < nSamples; isamp++ ) {
          samplesNew[isamp] = (samples1[isamp] + samples2[isamp])*0.5;
        }
      }
      else {
        int fftDataType = cseis_geolib::FX_REAL_IMAG;
        vars->fftTool1->forwardTransform( samples1, samples1, fftDataType );
        vars->fftTool2->forwardTransform( samples2, samples2, fftDataType );
        float* real1Ptr = &samples1[0];
        float* imag1Ptr = &samples1[vars->fftTool1->numFreqValues()];
        float const* real2Ptr = &samples2[0];
        float const* imag2Ptr = &samples2[vars->fftTool1->numFreqValues()];
        for( int ifreq = 0; ifreq < vars->fftTool1->numFreqValues(); ifreq++ ) {
          real1Ptr[ifreq] = (real1Ptr[ifreq] + real2Ptr[ifreq])*0.5f;
          imag1Ptr[ifreq] = (imag1Ptr[ifreq] + imag2Ptr[ifreq])*0.5f;
        }
        vars->fftTool1->inverseTransform( samples1, samplesNew, fftDataType );
      }

      // 2) Interpolate trace header values
      csTraceHeader* trcHdr1 = traceGather->trace(trcIndex1)->getTraceHeader();
      csTraceHeader* trcHdr2 = traceGather->trace(trcIndex2)->getTraceHeader();
      csTraceHeader* trcHdrNew = traceGather->trace(itrc)->getTraceHeader();

      for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
        char type = hdef->headerType(ihdr);
        switch( type ) {
        case TYPE_INT:
          trcHdrNew->setIntValue( ihdr, (trcHdr1->intValue(ihdr)+trcHdr2->intValue(ihdr))/2 );
          break;
        case TYPE_FLOAT:
          trcHdrNew->setFloatValue( ihdr, (trcHdr1->floatValue(ihdr)+trcHdr2->floatValue(ihdr))/2 );
          break;
        case TYPE_DOUBLE:
          trcHdrNew->setDoubleValue( ihdr, (trcHdr1->doubleValue(ihdr)+trcHdr2->doubleValue(ihdr))/2 );
          break;
        case TYPE_INT64:
          trcHdrNew->setInt64Value( ihdr, (trcHdr1->int64Value(ihdr)+trcHdr2->int64Value(ihdr))/2 );
          break;
        case TYPE_STRING:
          trcHdrNew->setStringValue( ihdr, trcHdr1->stringValue(ihdr) );
          break;
        default:
          writer->error("mod_trc_interpol: Unknown trace header type, code: %d", type);
        }
      }
      trcHdrNew->setIntValue( vars->hdrID_out, 1 );

    } // END for itrc
  } // END if MODE_ENSEMBLE
  else {
    if( numTraces < 3 ) {  // This must be the last trace of the data set. Simply return, do nothing
      return;
    }
    double value1 = traceGather->trace(1)->getTraceHeader()->doubleValue( vars->hdrID_in );
    double value2 = traceGather->trace(2)->getTraceHeader()->doubleValue( vars->hdrID_in );
    if( value1 == vars->value ) {
      // 2) Interpolate trace samples
      float* samples1 = traceGather->trace(0)->getTraceSamples();
      float* samples2 = traceGather->trace(2)->getTraceSamples();
      float* samplesMid = traceGather->trace(1)->getTraceSamples();
  
      if( vars->method == METHOD_SIMPLE_AVERAGE ) {
        for( int isamp = 0; isamp < nSamples; isamp++ ) {
          samplesMid[isamp] = (samples1[isamp] + samples2[isamp])*0.5;
        }
      }
      traceGather->trace(1)->getTraceHeader()->setIntValue( vars->hdrID_out, 1 );
    }
    if( value2 == vars->value ) {
      *numTrcToKeep = 2;
    }
    else {
      *numTrcToKeep = 1;
    }
  } // END if MODE_HEADER
}


//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_trc_interpol_( csParamDef* pdef ) {
  pdef->setModule( "TRC_INTERPOL", "Interpolate traces (EXPERIMENTAL MODULE)", "Interpolate 1 trace between every two adjacent traces, e.g. if 10 traces are input, 19 traces will be output" );

  pdef->addParam( "method", "", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "simple_average", "Interpolate new trace by simply averaging adjacent traces. Header values are averaged." );
  pdef->addOption( "fx_average", "..." );

  pdef->addParam( "mode", "", NUM_VALUES_FIXED );
  pdef->addValue( "fixed", VALTYPE_OPTION );
  pdef->addOption( "fixed", "Input 2 sequential traces at a time" );
  pdef->addOption( "ensemble", "Input one ensemble at a time" );
  pdef->addOption( "header", "Trace header indicates traces which shall be replaced by an interpolated trace" );

  pdef->addParam( "header", "Trace header indicating trace which shall be interpolated", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Value indicating this trace shall be replaced/interpolated" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_trc_interpol_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_trc_interpol::VariableStruct* vars = reinterpret_cast<mod_trc_interpol::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_trc_interpol_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_trc_interpol::VariableStruct* vars = reinterpret_cast<mod_trc_interpol::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->fftTool1 != NULL ) {
    delete vars->fftTool1;
    vars->fftTool1 = NULL;
  }
  if( vars->fftTool2 != NULL ) {
    delete vars->fftTool2;
    vars->fftTool2 = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_trc_interpol_( csParamDef* pdef ) {
  params_mod_trc_interpol_( pdef );
}
extern "C" void _init_mod_trc_interpol_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_trc_interpol_( param, env, writer );
}
extern "C" bool _start_exec_mod_trc_interpol_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_trc_interpol_( env, writer );
}
extern "C" void _exec_mod_trc_interpol_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_trc_interpol_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_trc_interpol_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_trc_interpol_( env, writer );
}

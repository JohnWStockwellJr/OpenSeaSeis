/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csGeolibUtils.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TIME_SLICE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_time_slice {
  class Dimension {
  public:
    Dimension() : val1(0), val2(0), inc(0) {}
    int val1;
    int val2;
    int inc;
    int nVal;
  };
  struct VariableStruct {
    int nSlices;
    int* sampleIndexSlice;
    int* hdrId_slice;
    cseis_system::csTraceGather* gather;

    Dimension dim1;
    Dimension dim2;
    int hdrId_dim1;
    int hdrId_dim2;

    float sampleIntIn;
    int numSamplesIn;
    int mode;
  };
  static int const MODE_HEADER = 1;
  static int const MODE_DATA   = 2;
}

using mod_time_slice::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_time_slice_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->nSlices     = 0;
  vars->hdrId_slice = NULL;
  vars->sampleIndexSlice = NULL;
  vars->gather = NULL;
  vars->mode   = mod_time_slice::MODE_HEADER;
  vars->sampleIntIn  = 0;
  vars->numSamplesIn = 0;
  vars->hdrId_dim1  = -1;
  vars->hdrId_dim2  = -1;

  //---------------------------------------------------------
  bool isTimeDomain = true;

  if( param->exists("domain") ) {
    std::string text;
    param->getString( "domain", &text );
    if( !text.compare( "sample" ) ) {
      isTimeDomain = false;
    }
    else if( !text.compare( "time" ) ) {
      isTimeDomain = true;
    }
    else {
      writer->line("Domain option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  //---------------------------------------------------------
  std::string text;

  if( param->exists("mode") ) {
    param->getString( "mode", &text );
    if( !text.compare( "header" ) ) {
      vars->mode = mod_time_slice::MODE_HEADER;
    }
    else if( !text.compare( "data" ) ) {
      vars->mode = mod_time_slice::MODE_DATA;
    }
    else {
      writer->error("Option not recognized: '%s'.", text.c_str());
    }
  }


  vars->nSlices = param->getNumLines("slice");
  if( vars->nSlices == 0 && vars->mode == mod_time_slice::MODE_HEADER ) {
    writer->error("No time slice specified, user parameter 'slice'");
  }
  if( vars->nSlices > 0 ) {
    int nHdrSlices = vars->nSlices;
    if( vars->mode == mod_time_slice::MODE_DATA ) nHdrSlices = 1;
    //    if( vars->mode == mod_time_slice::MODE_HEADER ) 
    vars->hdrId_slice = new int[nHdrSlices];
    for( int islice = 0; islice < nHdrSlices; islice++ ) {
      string hdrName;
      param->getStringAtLine("slice",&hdrName,islice,1);
      if( hdef->headerExists(hdrName) ) {
        type_t type = hdef->headerType(hdrName);
        if( type != TYPE_FLOAT && type != TYPE_DOUBLE) {
          writer->error("Time slice trace header '%s' should be of floating point type. Found: %s", hdrName.c_str(), csGeolibUtils::typeText(type) );
        }
      }
      else {
        hdef->addHeader( TYPE_FLOAT, hdrName, "Time slice");
      }
      vars->hdrId_slice[islice] = hdef->headerIndex(hdrName);
    }

    vars->sampleIndexSlice = new int[vars->nSlices];

    for( int islice = 0; islice < vars->nSlices; islice++ ) {
      float slice_in;
      param->getFloatAtLine("slice",&slice_in,islice,0);
      
      if( isTimeDomain ) {
        float sliceTime = slice_in;
        if( sliceTime < 0 ) {
          writer->error("Slice time (%f) needs to be greater or equal to 0.0.", sliceTime);
        }
        if( sliceTime > (float)(shdr->numSamples-1)*shdr->sampleInt ) {
          writer->error("Slice time (%f) exceeds length of trace (%f).", sliceTime, (float)(shdr->numSamples-1)*shdr->sampleInt );
        }
        vars->sampleIndexSlice[islice] = (int)(sliceTime / shdr->sampleInt);  // All in milliseconds / Hz
      }
      else {
        //
        // NOTE: User input is '1' for first sample. Internally, '0' is used!!
        //
        int sampleIndex = (int)slice_in;
        if( sampleIndex < 1 ) writer->error("Slice sample (%d) needs to be greater or equal to 1.", sampleIndex);
        if( sampleIndex > shdr->numSamples ) writer->error("Slice sample (%d) exceeds number of samples (%d).", sampleIndex, shdr->numSamples );
        vars->sampleIndexSlice[islice] = sampleIndex-1;   // see note above..
      }
    } // END for
  } // END: if nSlices > 0

  if( vars->mode == mod_time_slice::MODE_DATA ) {
    param->getString("data_dim1", &text, 0);
    if( hdef->headerType(text) != TYPE_INT ) writer->error("Trace header %s is not of type INT. This is currently not supported", text.c_str());
    vars->hdrId_dim1 = hdef->headerIndex(text);
    param->getInt("data_dim1", &vars->dim1.val1, 1);
    param->getInt("data_dim1", &vars->dim1.val2, 2);
    param->getInt("data_dim1", &vars->dim1.inc, 3);

    param->getString("data_dim2", &text, 0);
    if( hdef->headerType(text) != TYPE_INT ) writer->error("Trace header %s is not of type INT. This is currently not supported", text.c_str());
    vars->hdrId_dim2 = hdef->headerIndex(text);
    param->getInt("data_dim2", &vars->dim2.val1, 1);
    param->getInt("data_dim2", &vars->dim2.val2, 2);
    param->getInt("data_dim2", &vars->dim2.inc, 3);

    if( param->exists("data_slice") ) {
      double slice1 = 0;
      double slice2 = 0;
      double sliceInc = 0;

      if( vars->nSlices > 0 ) writer->error("Found both user parameters 'slice' and 'data_slice'. Only one of these can be specified at once.");
      param->getString("data_slice", &text, 0);
      if( !hdef->headerExists(text) ) {
        hdef->addHeader(TYPE_FLOAT,text,"Time slice");
      }
      vars->hdrId_slice = new int[1];
      vars->hdrId_slice[0] = hdef->headerIndex(text);
      param->getDouble("data_slice", &slice1, 1);
      slice2 = slice1;
      sliceInc = 1;
      if( param->getNumValues("data_slice") > 2 ) {
        param->getDouble("data_slice", &slice2, 2);
        sliceInc = slice2 - slice1;
        if( param->getNumValues("data_slice") > 3 ) {
          param->getDouble("data_slice", &sliceInc, 3);
        }
      }
      vars->nSlices = (int)round( (slice2 - slice1) / sliceInc ) + 1;
      vars->sampleIndexSlice = new int[vars->nSlices];
      for( int islice = 0; islice < vars->nSlices; islice++ ) {
        float sliceTime = (float)( (double)islice * sliceInc + slice1 );
        vars->sampleIndexSlice[islice] = (int)round(sliceTime / shdr->sampleInt);  // All in milliseconds / Hz
      }
    }

    vars->dim1.nVal = (int)round( (vars->dim1.val2 - vars->dim1.val1) / vars->dim1.inc ) + 1;
    vars->dim2.nVal = (int)round( (vars->dim2.val2 - vars->dim2.val1) / vars->dim2.inc ) + 1;

    // Clean up trace headers:
    for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
      int index = hdef->headerIndex(hdef->headerName(ihdr));
      if( index != vars->hdrId_slice[0] && index != vars->hdrId_dim2 ) {
        hdef->deleteHeader( hdef->headerName(ihdr) );
      }
    }
    hdef->resetByteLocation();

    vars->sampleIntIn  = shdr->sampleInt;
    vars->numSamplesIn = shdr->numSamples;
    shdr->numSamples = vars->dim1.nVal;
    shdr->sampleInt  = (float)vars->dim1.inc * 1000.0f;

    vars->gather = new csTraceGather( hdef );
    vars->gather->createTraces( 0, vars->dim2.nVal*vars->nSlices, hdef, shdr->numSamples );
    for( int islice = 0; islice < vars->nSlices; islice++ ) {
      float time = (float)vars->sampleIndexSlice[islice] * vars->sampleIntIn;
      for( int itrc = 0; itrc < vars->dim2.nVal; itrc++ ) {
        int indexTrace = islice * vars->dim2.nVal + itrc;
        csTraceHeader* trcHdr = vars->gather->trace(indexTrace)->getTraceHeader();
        trcHdr->setFloatValue( vars->hdrId_slice[0], time );
        trcHdr->setIntValue( vars->hdrId_dim2, itrc * vars->dim2.inc + vars->dim2.val1 );
        float* samples = vars->gather->trace(indexTrace)->getTraceSamples();
        for( int isamp = 0; isamp < vars->dim1.nVal; isamp++ ) {
          samples[isamp] = 0.0;
        }
      }
    }

  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_time_slice_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;


  float* samplesIn = traceGather->trace(0)->getTraceSamples();
  csTraceHeader* trcHdr = traceGather->trace(0)->getTraceHeader();

  if( vars->mode == mod_time_slice::MODE_HEADER ) {
    for( int islice = 0; islice < vars->nSlices; islice++ ) {
      trcHdr->setFloatValue( vars->hdrId_slice[islice], samplesIn[vars->sampleIndexSlice[islice]] );
    }
  }
  else {
    if( edef->isLastCall() ) {
      vars->gather->moveTracesTo( 0, vars->gather->numTraces(), traceGather );
    }
    else {
      int val_dim1 = trcHdr->intValue( vars->hdrId_dim1 );
      int val_dim2 = trcHdr->intValue( vars->hdrId_dim2 );
      int indexDim1 = (int)( (val_dim1-vars->dim1.val1)/vars->dim1.inc );
      int indexDim2 = (int)( (val_dim2-vars->dim2.val1)/vars->dim2.inc );
      bool isOK = ( (val_dim1 >= vars->dim1.val1) && (val_dim1 <= vars->dim1.val2) && (indexDim1*vars->dim1.inc+vars->dim1.val1 == val_dim1) );
      isOK = isOK && ( (val_dim2 >= vars->dim2.val1) && (val_dim2 <= vars->dim2.val2) && (indexDim2*vars->dim2.inc+vars->dim2.val1 == val_dim2) );
      if( isOK ) {
        if( indexDim1 >= shdr->numSamples ) {
          throw csException("Program bug: Incorrect samle index: %d  (numSamples = %d)\n", indexDim1, shdr->numSamples );
        }
        for( int islice = 0; islice < vars->nSlices; islice++ ) {
          int indexTrace  = islice * vars->dim2.nVal + indexDim2;
          if( indexTrace >= vars->gather->numTraces() ) {
            throw csException("Program bug: Incorrect trace index: %d  (numTraces = %d)\n", indexTrace, vars->gather->numTraces() );
          }
          float* samples = vars->gather->trace(indexTrace)->getTraceSamples();
          samples[indexDim1] = samplesIn[vars->sampleIndexSlice[islice]];
        }
      } // END isOK
      // else {
      //        writer->warning("Throwing out trace with trace header values (dim1/dim2):  %d / %d", val_dim1, val_dim2);
      // }
      traceGather->freeAllTraces();
      edef->setTracesAreWaiting();
    }
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_time_slice_( csParamDef* pdef ) {
  pdef->setModule( "TIME_SLICE", "Extract time slice(s) from input data, write to trace header field or output data" );

  pdef->addParam( "mode", "Time slice mode", NUM_VALUES_FIXED );
  pdef->addValue( "header", VALTYPE_OPTION );
  pdef->addOption( "header", "Place amplitudes of time slice in trace header only", "Specify user parameter 'slice'" );
  pdef->addOption( "data", "Write time slices to output data traces" );

  pdef->addParam( "slice", "Time slice", NUM_VALUES_VARIABLE, "This parameter can be repeated" );
  pdef->addValue( "", VALTYPE_NUMBER, "Time [ms] / frequency [Hz] or sample at which time slice shall be extracted" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name where time slice amplitude will be stored", "In case of 'data' mode, the same trace header should be specified for all slices" );

  pdef->addParam( "data_slice", "Time slice (output data trace & store in trace header)", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name where slice time will be stored" );
  pdef->addValue( "", VALTYPE_NUMBER, "Time/depth of first slice", "..in given units of trace (typically milliseconds, meters or Hertz)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Time/depth of last slice (optional)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Increment (optional)" );

  pdef->addParam( "data_dim1", "Definition of first dimension (vertical axis) of output time slice", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name for first dimension (vertical axis on output trace)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Start value" );
  pdef->addValue( "", VALTYPE_NUMBER, "End value" );
  pdef->addValue( "", VALTYPE_NUMBER, "Increment" );

  pdef->addParam( "data_dim2", "Definition of second dimension (horizontal axis) of output time slice", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name for first dimension (horizontal axis on output trace)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Start value" );
  pdef->addValue( "", VALTYPE_NUMBER, "End value" );
  pdef->addValue( "", VALTYPE_NUMBER, "Increment" );

  pdef->addParam( "domain", "Time or sample domain", NUM_VALUES_FIXED );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Window is specified in time [ms] (or frequency [Hz])" );
  pdef->addOption( "sample", "Window is specified in samples (1 for first sample)" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_time_slice_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_time_slice::VariableStruct* vars = reinterpret_cast<mod_time_slice::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_time_slice_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_time_slice::VariableStruct* vars = reinterpret_cast<mod_time_slice::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->gather != NULL ) {
    delete vars->gather;
    vars->gather = NULL;
  }
  if( vars->hdrId_slice != NULL ) {
    delete [] vars->hdrId_slice;
    vars->hdrId_slice = NULL;
  }
  if( vars->sampleIndexSlice != NULL ) {
    delete [] vars->sampleIndexSlice;
    vars->sampleIndexSlice = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_time_slice_( csParamDef* pdef ) {
  params_mod_time_slice_( pdef );
}
extern "C" void _init_mod_time_slice_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_time_slice_( param, env, writer );
}
extern "C" bool _start_exec_mod_time_slice_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_time_slice_( env, writer );
}
extern "C" void _exec_mod_time_slice_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_time_slice_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_time_slice_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_time_slice_( env, writer );
}

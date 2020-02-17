/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csFFTTools.h"
#include <cmath>
#include <cstring>
#include "csTimeFunction.h"
#include "csTableNew.h"
#include "csTableValueList.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: FILTER
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_filter {
  struct VariableStruct {
    csFFTTools* fftTool;
    csFFTTools** fftToolWin;
    int filterType;
    bool isPad;
    int numPaddedSamples;
    int numSamplesInclPad;
    float* bufferPaddedTrace;
    int output;
    float* bufferInput;
    int hdrID_both;

    // Windowed filter application
    int numWin; // Number of windows
    //    float* winTimesSamp; // Window end time in samples
    int* winStartSample; // Windows start samples (inclusive)
    int* winEndSample;   // Window end samples (inclusive)
    int* winNumSamples;  // Number of samples in windows
    int  winOverlapSamples; // Window overlap in number of samples 
    float* freqLowPass;
    float* freqHighPass;
    float* slopeLowPass;
    float* slopeHighPass;
    float* winBufferIn;
    float* winBufferOut;
    int numLowPassFilters;
    int numHighPassFilters;

    bool isNotchFilter;
    bool isNotchCosineTaper;
    float notchFreqHz;
    float notchWidthHz;
    float notchFilterSlope;

    csTableNew* table;
    int* hdrId_keys;
    float filterSampleInt;
    int filterNumSamples;
    int spatialMaxTime_nsamp;
  };
  static int const UNIT_HZ   = 1;
  static int const UNIT_PERCENT = 2;

  static int const TYPE_BUTTER  = 10;
  static int const TYPE_TABLE   = 11;
  static int const TYPE_SPATIAL = 12;

  static int const OUTPUT_FILT = 41;
  static int const OUTPUT_DIFF = 42;
  static int const OUTPUT_BOTH = 43;
  static int const OUTPUT_IMPULSE = 44;
  static int const OUTPUT_NOTCH = 45;

  static int const PAD_SELF       = 51;
  static int const PAD_CONTINUOUS = 52;

  void padTrace( mod_filter::VariableStruct* vars, float* samplesInOut, int numSamples );
  void padTrace2x( mod_filter::VariableStruct* vars, float* samplesInOut, int numSamples );

}
using namespace mod_filter;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_filter_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  csTraceHeaderDef* hdef = env->headerDef;
  edef->setVariables( vars );

  vars->fftTool       = NULL;
  vars->filterType = mod_filter::TYPE_BUTTER;

  vars->isPad             = false;
  vars->numPaddedSamples     = 0;
  vars->numSamplesInclPad = 0;
  vars->bufferPaddedTrace = NULL;
  vars->output = OUTPUT_FILT;
  vars->bufferInput = NULL;
  vars->hdrID_both = -1;

  vars->numWin = 0;
  vars->fftToolWin    = NULL;
  vars->winOverlapSamples = 0;
  vars->winStartSample = NULL;
  vars->winEndSample = NULL;
  vars->winNumSamples = NULL;
  vars->freqLowPass  = NULL;
  vars->freqHighPass = NULL;
  vars->slopeLowPass  = NULL;
  vars->slopeHighPass = NULL;
  vars->winBufferIn = NULL;
  vars->winBufferOut = NULL;

  vars->notchFreqHz = 0;
  vars->notchWidthHz = 0;
  vars->isNotchFilter = false;
  vars->isNotchCosineTaper = false;
  vars->notchFilterSlope = 72;

  vars->table          = NULL;
  vars->hdrId_keys     = NULL;
  vars->filterSampleInt = shdr->sampleInt;
  vars->filterNumSamples = shdr->numSamples;

  edef->setMPISupport( true ); // MPI is supported by this module
//---------------------------------------------
//
  std::string text;
  if( param->exists("type") ) {
    param->getString("type", &text );
    if( !text.compare("butterworth") ) {
      vars->filterType = mod_filter::TYPE_BUTTER;
    }
    else if( !text.compare("table") ) {
      vars->filterType = mod_filter::TYPE_TABLE;
    }
    else if( !text.compare("spatial") ) {
      vars->filterType = mod_filter::TYPE_SPATIAL;
      vars->filterSampleInt = 1.0f;
      param->getInt("spatial_param", &vars->filterNumSamples, 0 );
      float maxTime;
      param->getFloat("spatial_param", &maxTime, 1 );
      vars->spatialMaxTime_nsamp = (int)round( maxTime / shdr->sampleInt );
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  if( vars->filterType != mod_filter::TYPE_SPATIAL ) {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }
  else {
    edef->setTraceSelectionMode( TRCMODE_FIXED, vars->filterNumSamples );
  }

  bool doRestoreFilter = false;
  if( param->exists("mode") ) {
    param->getString("mode", &text );
    if( !text.compare("restore") ) {
      doRestoreFilter = true;
    }
    else if( !text.compare("apply") ) {
      doRestoreFilter = false;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

//---------------------------------------------
//
  if( param->exists("pad") ) {
    float padLength_ms;
    param->getFloat("pad", &padLength_ms, 0 );
    vars->numPaddedSamples = (int)round( padLength_ms / vars->filterSampleInt );
    vars->numSamplesInclPad = 2 * vars->numPaddedSamples + vars->filterNumSamples;
    vars->bufferPaddedTrace = new float[vars->numSamplesInclPad];
    vars->isPad = true;
    if( vars->filterNumSamples/2 < vars->numPaddedSamples ) {
      writer->error("Too much padding: Pad length cannot exceed 1/2 trace length. Maximum pad length = %fms", (float)vars->filterNumSamples*vars->filterSampleInt/2.0f);
    }
  }

//---------------------------------------------
//
  int freqUnit = mod_filter::UNIT_HZ;
  if( param->exists("unit") ) {
    std::string text;
    param->getString("unit", &text);
    if( !text.compare("hz") ) {
      freqUnit = mod_filter::UNIT_HZ;
    }
    else if( !text.compare("percent") ) {
      freqUnit = mod_filter::UNIT_PERCENT;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str() );
    }
  }

//---------------------------------------------
// Retrieve filter scalars from ASCII table
//
  bool isWindowedApp = param->exists("win_times");
  if( param->exists("table") ) {
    if( vars->filterType != mod_filter::TYPE_TABLE ) writer->error("When reading in filter scalars from a table, you need to specify user paremeter 'type table'");
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
        if( !isWindowedApp || ikey > 0 ) { // First key in windowed application must be "window" which is not a trace header
          if( !hdef->headerExists( headerName ) ) {
            writer->error("No matching trace header found for table key '%s'", headerName.c_str() );
          }
          vars->hdrId_keys[ikey] = hdef->headerIndex( headerName );
        }
      } // END for ikey
    }

    // Add scalar 'value' column to input table
    vars->table->addValue( colScalar-1 );  // -1 to convert from 'user' column to 'C++' column
    bool sortTable = false;
    try {
      vars->table->initialize( tableFilename, sortTable );
      if( edef->isDebug() ) vars->table->dump();
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s': %s\n", tableFilename.c_str(), exc.getMessage() );
    }
  }

  vars->numLowPassFilters  = param->getNumLines("lowpass");
  vars->numHighPassFilters = param->getNumLines("highpass");
  float freqNy = 500.0f/vars->filterSampleInt;
  float slope;
  if( vars->numLowPassFilters  != 0 ) {
    vars->freqLowPass   = new float[vars->numLowPassFilters];
    vars->slopeLowPass  = new float[vars->numLowPassFilters];
    for( int i = 0; i < vars->numLowPassFilters; i++ ) {
      param->getFloatAtLine("lowpass", &vars->freqLowPass[i], i, 0 );
      param->getFloatAtLine("lowpass", &slope, i, 1 );
      vars->slopeLowPass[i] = fabs(slope);
      if( vars->slopeLowPass[i] < 0.6 ) {
        writer->error("Specified lowpass slope %f outside of valid range", slope);
      }
      if( doRestoreFilter ) {
        vars->slopeLowPass[i]  = -vars->slopeLowPass[i];
      }
      if( freqUnit == mod_filter::UNIT_PERCENT ) {
        vars->freqLowPass[i] *= freqNy / 100.0f;
      }
      if( vars->freqLowPass[i] < 0 || vars->freqLowPass[i] > freqNy ) {
        writer->error("Lowpass filter frequency exceeds valid range (0-%fHz): %fHz", freqNy, vars->freqLowPass[i] );
      }
    }
  }
  if( vars->numHighPassFilters != 0 ) {
    vars->freqHighPass   = new float[vars->numHighPassFilters];
    vars->slopeHighPass  = new float[vars->numHighPassFilters];
    for( int i = 0; i < vars->numHighPassFilters; i++ ) {
      param->getFloatAtLine("highpass", &vars->freqHighPass[i], i, 0 );
      param->getFloatAtLine("highpass", &slope, i, 1 );
      vars->slopeHighPass[i] = fabs(slope);
      if( vars->slopeHighPass[i] < 0.6 ) {
        writer->error("Specified highpass slope %f outside of valid range", slope);
      }
      if( doRestoreFilter ) {
        vars->slopeHighPass[i]  = -vars->slopeHighPass[i];
      }
      if( freqUnit == mod_filter::UNIT_PERCENT ) {
        vars->freqHighPass[i] *= freqNy / 100.0f;
      }
      if( vars->freqHighPass[i] < 0 || vars->freqHighPass[i] > freqNy ) {
        writer->error("Highpass filter frequency exceeds valid range (0-%fHz): %fHz", freqNy, vars->freqHighPass[i] );
      }
    }
  }

  if( param->exists("output") ) {
    param->getString("output", &text );
    if( !text.compare("filt") ) {
      vars->output = OUTPUT_FILT;
    }
    else if( !text.compare("diff") ) {
      vars->output = OUTPUT_DIFF;
    }
    else if( !text.compare("both") ) {
      vars->output = OUTPUT_BOTH;
    }
    else if( !text.compare("impulse") ) {
      vars->output = OUTPUT_IMPULSE;
    }
    else if( !text.compare("notch_filter") ) {
      vars->output = OUTPUT_NOTCH;
      shdr->numSamples = vars->fftTool->numFFTSamples();
      shdr->sampleInt  = vars->fftTool->sampleIntFreq();
      vars->filterSampleInt = shdr->sampleInt;
      vars->filterNumSamples = shdr->numSamples;
      shdr->domain = DOMAIN_FX;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
    if( vars->output == mod_filter::OUTPUT_DIFF || vars->output == mod_filter::OUTPUT_BOTH ) {
      vars->bufferInput = new float[vars->filterNumSamples];
    }
    if( vars->output == mod_filter::OUTPUT_BOTH ) {
      text = "trcno";
      if( param->getNumValues("output") > 1 ) {
        param->getString("output",&text,1);
      }
      vars->hdrID_both = hdef->headerIndex( text.c_str() );
    }
  }


  if( param->exists("win_times") ) {
    if( vars->isPad ) writer->error("Windowed filter application does not work in conjunction with padding. Remove user parameter 'pad' and retry.");
    if( doRestoreFilter ) writer->error("Windowed filter application does not suppport  the 'mode restore' option. Remove user parameter 'mode' and retry.");
    if( vars->output == mod_filter::OUTPUT_IMPULSE ) writer->error("Windowed filter application does not work in conjunction with impulse output option.");
    if( freqUnit == mod_filter::UNIT_PERCENT ) writer->error("Windowed filter application does not support percent option. Remove user parameter 'unit' and retry.");
    vars->numWin = param->getNumValues("win_times");

    if( vars->filterType != mod_filter::TYPE_TABLE ) {
      if( vars->numWin == 0 || (vars->numLowPassFilters!=vars->numWin && vars->numHighPassFilters!=vars->numWin) ) {
        writer->error("Inconsistent window definition. Need to specify same number of low (=%d) and/or high (=%d) pass points as window start times (=%d)", vars->numLowPassFilters, vars->numHighPassFilters, vars->numWin);
      }
    }
    else { // Table input
      if( vars->table->numKeys() == 0 ) {
        writer->error("Windowed filter application using an input tablerequires a key column in the input table containing the number of the application window, starting with 1 for the first window.");
      }
    }
    vars->winOverlapSamples = (int)(500 / vars->filterSampleInt);
    if( param->exists("win_overlap") ) {
      float overlap_ms = 0;
      param->getFloat("win_overlap",&overlap_ms);
      vars->winOverlapSamples = (int)round( 0.5*overlap_ms / vars->filterSampleInt ) * 2;
    }
    if( ((int)(vars->winOverlapSamples/2))*2 != vars->winOverlapSamples ) vars->winOverlapSamples += 1; // Must be even

    vars->winStartSample = new int[vars->numWin];
    vars->winEndSample   = new int[vars->numWin];
    for( int i = 0; i < vars->numWin; i++ ) {
      float time;
      param->getFloat("win_times", &time, i );
      vars->winStartSample[i] = (int)(round(time / vars->filterSampleInt));
    }

    vars->winStartSample[0] = 0;
    vars->winEndSample[vars->numWin-1] = vars->filterNumSamples-1;
    for( int i = 1; i < vars->numWin; i++ ) {
      vars->winEndSample[i-1] = vars->winStartSample[i] + vars->winOverlapSamples;
      vars->winStartSample[i] -= vars->winOverlapSamples;
      if( vars->winStartSample[i] < 0 ) writer->error("Window #%d, including overlap (=%.2fms), is too large.", i+1, vars->winOverlapSamples*vars->filterSampleInt);
      if( vars->winEndSample[i-1] >= vars->filterNumSamples ) writer->error("Window #%d, including overlap (=%.2fms), is too large.", i+1, vars->winOverlapSamples*vars->filterSampleInt);
    }
    vars->fftToolWin = new csFFTTools*[vars->numWin];
    for( int i = 0; i < vars->numWin; i++ ) {
      int numSamplesWin = vars->winEndSample[i] - vars->winStartSample[i] + 1;
      vars->fftToolWin[i] = new csFFTTools( numSamplesWin, vars->filterSampleInt );
    }
    vars->winBufferIn  = new float[vars->filterNumSamples];
    vars->winBufferOut = new float[vars->filterNumSamples];

    writer->line("Window overlap time/samples:  %8.2fms/%-5d", vars->winOverlapSamples*vars->filterSampleInt, vars->winOverlapSamples);
    for( int iwin = 0; iwin < vars->numWin; iwin++ ) {
      writer->line("Window #%d: Start time/sample %8.2fms/%-5d, end time/sample: %8.2fms/%-5d, number of samples: %d", iwin+1,
              vars->filterSampleInt*vars->winStartSample[iwin], vars->winStartSample[iwin], vars->filterSampleInt*vars->winEndSample[iwin], vars->winEndSample[iwin], vars->winEndSample[iwin]-vars->winStartSample[iwin]+1 );
    }
  } // END win_times
  //--------------------------------------------------------------------------------
  //
  if( vars->numWin == 0 && (vars->filterType == mod_filter::TYPE_BUTTER || vars->filterType == mod_filter::TYPE_SPATIAL) ) {

    if( param->exists("notch_filter") ) {
      param->getFloat("notch_filter", &vars->notchFreqHz, 0);
      param->getFloat("notch_filter", &vars->notchWidthHz, 1);
      param->getFloat("notch_filter", &vars->notchFilterSlope, 2);
      //      vars->notchFilterOrder = fabs(vars->notchFilterOrder) / 6.0;
      std::string text;
      param->getString("notch_filter", &text, 3);
      if( !text.compare("cosine") ) {
        vars->isNotchCosineTaper = true;
      }
      else if( !text.compare("butter") ) {
        vars->isNotchCosineTaper = false;
      }
      else {
        writer->error("Unknown option: '%s'", text.c_str() );
      }
      vars->isNotchFilter = true;
    }

    if( vars->numLowPassFilters == 0 && vars->numHighPassFilters == 0 && !vars->isNotchFilter ) {
      writer->error("No filter option specified. Specify for user parameter 'lowpass' and/or 'highpass' and/or 'notch_filter'");
    }
  } // END: Setup bandpass filter

  if( vars->numWin == 0 ) {
    if( !vars->isPad ) {
      vars->fftTool = new csFFTTools( vars->filterNumSamples, vars->filterSampleInt );
    }
    else {
      vars->fftTool = new csFFTTools( vars->numSamplesInclPad, vars->filterSampleInt );
    }
  }


  if( vars->numWin == 0 ) {
    if( vars->isNotchFilter ) {
      vars->fftTool->prepareNotchFilter( vars->notchFreqHz, vars->notchWidthHz, vars->notchFilterSlope, vars->isNotchCosineTaper ); 
    }
    if( vars->numLowPassFilters > 0 || vars->numHighPassFilters > 0 ) {
      vars->fftTool->prepareBandpassFilter( vars->numLowPassFilters, vars->freqLowPass, vars->slopeLowPass, vars->numHighPassFilters, vars->freqHighPass, vars->slopeHighPass );
    }
  }
  else {
    int numLP = ( vars->numLowPassFilters > 0 ) ? 1 : 0;
    int numHP = ( vars->numHighPassFilters > 0 ) ? 1 : 0;
    float* freqLP  = NULL;
    float* slopeLP = NULL;
    float* freqHP  = NULL;
    float* slopeHP = NULL;
    for( int iwin = 0; iwin < vars->numWin; iwin++ ) {
      if( numLP > 0 ) {
        freqLP  = &vars->freqLowPass[iwin];
        slopeLP = &vars->slopeLowPass[iwin];
      }
      if( numHP > 0 ) {
        freqHP  = &vars->freqHighPass[iwin];
        slopeHP = &vars->slopeHighPass[iwin];
      }
      vars->fftToolWin[iwin]->prepareBandpassFilter( numLP, freqLP, slopeLP, numHP, freqHP, slopeHP );
    }
  }
  writer->line("Number of user specified low/high pass filters: %d + %d\n", vars->numLowPassFilters, vars->numHighPassFilters );

  if( edef->isDebug() ) {
    for( int i = 0; i < vars->numLowPassFilters; i++ ) {
      writer->line("LOWPASS  #%d:  %.2f Hz  %.2f db/oct", i+1, vars->freqLowPass[i], vars->slopeLowPass[i] );
    }
    for( int i = 0; i < vars->numHighPassFilters; i++ ) {
      writer->line("HIGHPASS #%d:  %.2f Hz  %.2f db/oct", i+1, vars->freqHighPass[i], vars->slopeHighPass[i] );
    }
  }
}
//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_filter_(
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
  float* samplesInOut = trace->getTraceSamples();
  int numSamples = shdr->numSamples;
  int numTraces = traceGather->numTraces();
  
  if( vars->filterType == mod_filter::TYPE_SPATIAL ) {
    if( numTraces != vars->filterNumSamples ) writer->error("Inconsistent number of traces in gather: %d != %d(expected)", numTraces, vars->filterNumSamples);

    float* samplesSpatial = new float[numSamples];
    for( int isampActual = 0; isampActual < vars->spatialMaxTime_nsamp; isampActual++ ) {
      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        samplesInOut = traceGather->trace(itrc)->getTraceSamples();
        samplesSpatial[itrc] = samplesInOut[isampActual];
      }
      if( vars->isPad ) {
        padTrace( vars, samplesSpatial, numSamples );
        samplesSpatial = vars->bufferPaddedTrace;
      }
      
      vars->fftTool->applyBandpassFilter( samplesSpatial, false );

      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        samplesInOut = traceGather->trace(itrc)->getTraceSamples();
        samplesInOut[isampActual] = samplesSpatial[itrc];
      }
    }
    delete [] samplesSpatial;
    return;
  }

  if( vars->output == mod_filter::OUTPUT_DIFF || vars->output == mod_filter::OUTPUT_BOTH ) {
    memcpy( vars->bufferInput, samplesInOut, numSamples*sizeof(float) );
  }

  // Padding: Pad samples and extrapolate data samples. Apply cosine taper to padded data.
  if( vars->isPad ) {
    padTrace( vars, samplesInOut, numSamples );
    samplesInOut = vars->bufferPaddedTrace; // Direct samples pointer to padded buffer
    numSamples   = vars->numSamplesInclPad;
    if( edef->isDebug() ) {
      for( int isamp = 0; isamp < vars->numSamplesInclPad; isamp++ ) {
        fprintf(stdout,"%f %e\n", (isamp - vars->numPaddedSamples)*vars->filterSampleInt, samplesInOut[isamp]);
      }
    }
  }

  //----------------------------------------------------------------------
  // Windowed filter application
  if( vars->numWin != 0 ) {
    // Zero output buffer
    for( int isamp = 0; isamp < numSamples; isamp++ ) {
      vars->winBufferOut[isamp] = 0.0;
    }
    for( int iwin = 0; iwin < vars->numWin; iwin++ ) {
      int numSamplesWin = vars->winEndSample[iwin] - vars->winStartSample[iwin] + 1;
      // Copy windowed data to input buffer
      memcpy( vars->winBufferIn, &samplesInOut[vars->winStartSample[iwin]], numSamplesWin*sizeof(float) );
      // Apply cosine taper to 1/2 overlap samples at each end
      int numOverlapSamplesHalf = vars->winOverlapSamples/2;
      for( int isamp = 0; isamp < numOverlapSamplesHalf; isamp++ ) {
        float taper = 0.5 * ( 1.0 - cos( M_PI * ( (float)isamp/(float)numOverlapSamplesHalf ) ) );
        vars->winBufferIn[isamp] *= taper;
        vars->winBufferIn[numSamplesWin-1-isamp] *= taper;
      }
      // Apply low-pass and high-pass filters:
      vars->fftToolWin[iwin]->applyBandpassFilter( vars->winBufferIn, (vars->output == OUTPUT_IMPULSE) );
      if( vars->table != NULL ) {
        csTimeFunction<double> const* freqFunc;
        if( vars->table->numKeys() > 0 ) {
          double* keyValueBuffer = new double[vars->table->numKeys()];
          keyValueBuffer[0] = iwin + 1; // +1 to convert to 'user' number
          for( int ikey = 1; ikey < vars->table->numKeys(); ikey++ ) {
            keyValueBuffer[ikey] = trace->getTraceHeader()->doubleValue( vars->hdrId_keys[ikey] );
          }
          freqFunc = vars->table->getFunction( keyValueBuffer, false );
          delete [] keyValueBuffer;
        }
        else {
          freqFunc = vars->table->getFunction( NULL, false );
        }
        vars->fftToolWin[iwin]->applyFilter( vars->winBufferIn, freqFunc );
      }
      int samp1Copy = vars->winStartSample[iwin] + vars->winOverlapSamples + numOverlapSamplesHalf;
      int samp2Copy = vars->winEndSample[iwin]   - vars->winOverlapSamples - numOverlapSamplesHalf;
      int sampFrom  = vars->winOverlapSamples + numOverlapSamplesHalf;
      if( iwin == 0 ) {
        samp1Copy = 0;
        sampFrom = 0;
      }
      else if( iwin == vars->numWin-1 ) {
        samp2Copy = vars->winEndSample[iwin];
      }
      if( samp2Copy >= samp1Copy ) memcpy( &vars->winBufferOut[samp1Copy], &vars->winBufferIn[sampFrom], (samp2Copy-samp1Copy+1)*sizeof(float) );
      // b) Add overlap data to output buffer using linear taper
      if( iwin > 0 ) {
        for( int isamp = 0; isamp < vars->winOverlapSamples; isamp++ ) {
          float scalar = (float)(isamp) / (float)(vars->winOverlapSamples-1);
          vars->winBufferOut[isamp+numOverlapSamplesHalf+vars->winStartSample[iwin]] += scalar*vars->winBufferIn[isamp+numOverlapSamplesHalf];
        }
      }
      if( iwin < vars->numWin-1 ) {
        for( int isamp = 1; isamp <= vars->winOverlapSamples; isamp++ ) {
          float scalar = (float)(isamp-1) / (float)(vars->winOverlapSamples-1);
          vars->winBufferOut[vars->winEndSample[iwin]-isamp-numOverlapSamplesHalf+1] += scalar*vars->winBufferIn[numSamplesWin-isamp-numOverlapSamplesHalf];
        }
      }
      
    } // END: Windowed filter application
    memcpy( samplesInOut, vars->winBufferOut, numSamples*sizeof(float) );
  }
  //----------------------------------------------------------------------
  // Non-windowed filter application
  else if( vars->numLowPassFilters > 0 || vars->numHighPassFilters > 0 ) {
    if( vars->filterType == mod_filter::TYPE_BUTTER ) {
      vars->fftTool->applyBandpassFilter( samplesInOut, (vars->output == OUTPUT_IMPULSE) );
    }
    else if( vars->filterType == mod_filter::TYPE_TABLE ) {
      csTimeFunction<double> const* freqFunc;

      if( vars->table->numKeys() > 0 ) {
        double* keyValueBuffer = new double[vars->table->numKeys()];
        for( int ikey = 0; ikey < vars->table->numKeys(); ikey++ ) {
          keyValueBuffer[ikey] = trace->getTraceHeader()->doubleValue( vars->hdrId_keys[ikey] );
        }
        freqFunc = vars->table->getFunction( keyValueBuffer, false );
        delete [] keyValueBuffer;
      }
      else {
        freqFunc = vars->table->getFunction( NULL, false );
      }
      vars->fftTool->applyFilter( samplesInOut, freqFunc );
    }
  }
  if( vars->isNotchFilter ) {
    if( vars->output == mod_filter::OUTPUT_NOTCH ) {
      float const* ptr = vars->fftTool->prepareNotchFilter( vars->notchFreqHz, vars->notchWidthHz, vars->notchFilterSlope, vars->isNotchCosineTaper );
      for( int isamp = 0; isamp < numSamples; isamp++ ) {
        samplesInOut[isamp] = (float)ptr[isamp];
      }
      if( vars->isPad ) {
        memcpy( trace->getTraceSamples(), &samplesInOut[vars->numPaddedSamples], shdr->numSamples * sizeof(float) );    
      }
      return;
    }
    vars->fftTool->applyNotchFilter( samplesInOut, (vars->output == OUTPUT_IMPULSE) );
  }

  // Padding: Copy non-padded & filtered samples back to input trace
  if( vars->isPad ) {
    memcpy( trace->getTraceSamples(), &samplesInOut[vars->numPaddedSamples], shdr->numSamples * sizeof(float) );    
    numSamples = shdr->numSamples; // Reset numSamples to input/output trace length
  }

  //--------------------------------------------------------------------------------
  //
  if( vars->output == mod_filter::OUTPUT_DIFF || vars->output == mod_filter::OUTPUT_BOTH ) {
    float* samplesOut = samplesInOut;
    if( vars->output == mod_filter::OUTPUT_BOTH ) {
      traceGather->createTrace( hdef, numSamples);
      traceGather->trace(1)->getTraceHeader()->copyFrom( trace->getTraceHeader() );
      samplesOut = traceGather->trace(1)->getTraceSamples();
      traceGather->trace(0)->getTraceHeader()->setIntValue( vars->hdrID_both,1 );
      traceGather->trace(1)->getTraceHeader()->setIntValue( vars->hdrID_both,2 );
    }
    // Compute difference:
    for( int isamp = 0; isamp < numSamples; isamp++ ) {
      samplesOut[isamp] = vars->bufferInput[isamp] - samplesInOut[isamp];
    }
  }
  
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_filter_( csParamDef* pdef ) {
  pdef->setModule( "FILTER", "Frequency filter" );
  //  pdef->setVersion(1,0);

  pdef->addParam( "type", "Filter type", NUM_VALUES_FIXED );
  pdef->addValue( "butterworth", VALTYPE_OPTION );
  pdef->addOption( "butterworth", "Butterworth filter" );
  pdef->addOption( "table", "Zero-phase filter. Read in coefficients/scalars from ASCII table. Specify user parameters 'table' and 'table_col'" );
  pdef->addOption( "spatial", "Apply filter spatially across ensemble, sample by sample" );

  pdef->addParam( "lowpass", "Lowpass filter", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Cutoff frequency for low-pass filter [Hz]", "The cutoff frequency will be damped by -3db" );
  pdef->addValue( "12", VALTYPE_NUMBER, "Filter slope [dB/oct]" );

  pdef->addParam( "highpass", "Highpass filter", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Cutoff frequency for high-pass filter [Hz]", "The cutoff frequency will be damped by -3db" );
  pdef->addValue( "12", VALTYPE_NUMBER, "Filter slope [dB/oct]" );

  pdef->addParam( "table", "ASCII able", NUM_VALUES_FIXED, "Table format: The ASCII file should contain only numbers, no text. Required are two columns containing frequency and scalar pairs. Optional: Up to 2 key values. Lines starting with '#' are considered comment lines. For windowed application, the first specified key must be the window number" );
  pdef->addValue( "", VALTYPE_STRING, "Filename containing table");

  pdef->addParam( "table_col", "Table columns containing frequency in [Hz] and scalar value", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing frequency [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing scalar value" );

  pdef->addParam( "table_key", "Key trace header used to match values found in specified table columns", NUM_VALUES_VARIABLE,
                  "Specify the 'table_key' parameter for each key in the table file" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name of key header" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table" );
  pdef->addValue( "yes", VALTYPE_OPTION, "Interpolate based to this key?" );
  pdef->addOption( "yes", "Use this key for interpolation of value" );
  pdef->addOption( "no", "Do not use this key for interpolation", "The input table is expected to contain the exact key values for this trace header" );

  pdef->addParam( "unit", "Unit of frequency values supplied in parameters", NUM_VALUES_FIXED );
  pdef->addValue( "hz", VALTYPE_OPTION );
  pdef->addOption( "hz", "Frequencies are specified in [Hz]" );
  pdef->addOption( "percent", "Frequencies are specified as percent of Nyquist" );

  pdef->addParam( "pad", "Pad trace to avoid filter edge effects", NUM_VALUES_FIXED, "Pad and extrapolate trace at top and bottom before filter application. Remove padded samples after filter." );
  pdef->addValue( "0", VALTYPE_NUMBER, "Pad length at top and bottom, in units of trace (for example [ms])" );

  pdef->addParam( "mode", "Mode of operation: Apply specified filter or inverse filter", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply filter" );
  pdef->addOption( "restore", "'Restore' filter: Apply inverse filter" );

  pdef->addParam( "output", "Output data", NUM_VALUES_VARIABLE );
  pdef->addValue( "filt", VALTYPE_OPTION );
  pdef->addOption( "filt", "Output filtered data" );
  pdef->addOption( "impulse", "Output filter impulse response, zero phase, placed at trace centre." );
  pdef->addOption( "diff", "Output difference between input and filtered data" );
  pdef->addOption( "both", "Output both filtered data and difference between input and filtered data" );
  pdef->addOption( "notch_filter", "Output notch filter" );
  pdef->addValue( "trcno", VALTYPE_STRING, "Optional, used for option 'both': Trace header to distinguish between filtered data (1) and difference (2)" );

  pdef->addParam( "win_times", "List of N window start times for windowed filter application", NUM_VALUES_VARIABLE, "Specify N lines with lowpass and/or highpass filter points, user parameters 'lowpass' and 'highpass', or an input table with filter scalars" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of times [ms]" );

  pdef->addParam( "win_overlap", "Overlap between windows", NUM_VALUES_FIXED, "Overlap is centred around window start times" );
  pdef->addValue( "500", VALTYPE_NUMBER, "Window time overlap [ms]" );

  pdef->addParam( "notch_filter", "Apply notch filter", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Notch frequency [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Width of notch filter [Hz]" );
  pdef->addValue( "72", VALTYPE_NUMBER, "Notch filter slope [dB/oct]" );
  pdef->addValue( "butter", VALTYPE_OPTION, "Filter type" );
  pdef->addOption( "butter", "Use Butterworth filter" );
  pdef->addOption( "cosine", "Use cosine taper" );

  pdef->addParam( "spatial_param", "Specify for spatial filter application", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Number of traces in ensemble" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum time to apply filter" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_filter_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_filter::VariableStruct* vars = reinterpret_cast<mod_filter::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_filter_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_filter::VariableStruct* vars = reinterpret_cast<mod_filter::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->table != NULL ) {
    delete vars->table;
    vars->table = NULL;
  }
  if( vars->hdrId_keys != NULL ) {
    delete [] vars->hdrId_keys;
    vars->hdrId_keys = NULL;
  }
  if( vars->fftTool != NULL ) {
    delete vars->fftTool;
    vars->fftTool = NULL;
  }
  if( vars->bufferPaddedTrace != NULL ) {
    delete [] vars->bufferPaddedTrace;
    vars->bufferPaddedTrace = NULL;
  }
  if( vars->bufferInput != NULL ) {
    delete [] vars->bufferInput;
    vars->bufferInput = NULL;
  }
  if( vars->winBufferIn != NULL ) {
    delete [] vars->winBufferIn;
    vars->winBufferIn = NULL;
  }
  if( vars->winBufferOut != NULL ) {
    delete [] vars->winBufferOut;
    vars->winBufferOut = NULL;
  }
  if( vars->fftToolWin != NULL ) {
    for( int iwin = 0; iwin < vars->numWin; iwin++ ) {
      delete vars->fftToolWin[iwin];
    }
    delete [] vars->fftToolWin;
    vars->fftToolWin = NULL;
  }
  if( vars->winStartSample != NULL ) {
    delete [] vars->winStartSample;
    vars->winStartSample= NULL;
  }
  if( vars->winEndSample != NULL ) {
    delete [] vars->winEndSample;
    vars->winEndSample = NULL;
  }
  if( vars->winNumSamples != NULL ) {
    delete [] vars->winNumSamples;
    vars->winNumSamples = NULL;
  }
  if( vars->freqLowPass  != NULL ) {
    delete [] vars->freqLowPass;
    vars->freqLowPass  = NULL;
  }
  if( vars->freqHighPass != NULL ) {
    delete [] vars->freqHighPass;
    vars->freqHighPass = NULL;
  }
  if( vars->slopeLowPass != NULL ) {
    delete [] vars->slopeLowPass;
    vars->slopeLowPass  = NULL;
  }
  if( vars->slopeHighPass != NULL ) {
    delete [] vars->slopeHighPass;
    vars->slopeHighPass = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_filter_( csParamDef* pdef ) {
  params_mod_filter_( pdef );
}
extern "C" void _init_mod_filter_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_filter_( param, env, writer );
}
extern "C" bool _start_exec_mod_filter_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_filter_( env, writer );
}
extern "C" void _exec_mod_filter_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_filter_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_filter_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_filter_( env, writer );
}


namespace mod_filter {

void padTrace( mod_filter::VariableStruct* vars, float* samplesInOut, int numSamples ) {
  for( int isamp = 0; isamp < vars->numPaddedSamples; isamp++ ) {
    vars->bufferPaddedTrace[isamp] = 0;
  }
  for( int isamp = vars->numSamplesInclPad-vars->numPaddedSamples; isamp < vars->numSamplesInclPad; isamp++ ) {
    vars->bufferPaddedTrace[isamp] = 0;
  }
  // Copy input data
  memcpy( &vars->bufferPaddedTrace[vars->numPaddedSamples], samplesInOut, numSamples * sizeof(float) );
  // Extrapolate start of trace
  for( int isamp = vars->numPaddedSamples-1; isamp >= 0; isamp-- ) {
    float ratio = (float)isamp/(float)(vars->numPaddedSamples-1);
    float taper = 0.5 * ( 1 + cos( M_PI*(ratio-1.0) ) );
    vars->bufferPaddedTrace[isamp] = taper * ( samplesInOut[vars->numPaddedSamples-isamp] );
  }
  // Extrapolate end of trace
  int num = 2*numSamples+vars->numPaddedSamples-2;
  for( int isamp = numSamples+vars->numPaddedSamples; isamp < vars->numSamplesInclPad; isamp++ ) {
    float ratio = (float)(vars->numSamplesInclPad-isamp-1)/(float)(vars->numPaddedSamples-1);
    float taper = 0.5 * ( 1 + cos( M_PI*(ratio-1.0) ) );
    vars->bufferPaddedTrace[isamp] = taper * ( samplesInOut[num-isamp] );
  }
}

void padTrace2x( mod_filter::VariableStruct* vars, float* samplesInOut, int numSamples ) {
  for( int isamp = 0; isamp < vars->numPaddedSamples; isamp++ ) {
    vars->bufferPaddedTrace[isamp] = 0;
  }
  for( int isamp = vars->numSamplesInclPad-vars->numPaddedSamples; isamp < vars->numSamplesInclPad; isamp++ ) {
    vars->bufferPaddedTrace[isamp] = 0;
  }
  // Copy input data
  memcpy( &vars->bufferPaddedTrace[vars->numPaddedSamples], samplesInOut, numSamples * sizeof(float) );
  // Extrapolate start of trace
  float value2x = 2 * samplesInOut[0];
  for( int isamp = vars->numPaddedSamples-1; isamp >= 0; isamp-- ) {
    float ratio = (float)isamp/(float)(vars->numPaddedSamples-1);
    float taper = 0.5 * ( 1 + cos( M_PI*(ratio-1.0) ) );
    vars->bufferPaddedTrace[isamp] = taper * ( value2x - samplesInOut[vars->numPaddedSamples-isamp] );
  }
  // Extrapolate end of trace
  value2x = 2 * samplesInOut[numSamples - 1];
  int num = 2*numSamples+vars->numPaddedSamples-2;
  for( int isamp = numSamples+vars->numPaddedSamples; isamp < vars->numSamplesInclPad; isamp++ ) {
    float ratio = (float)(vars->numSamplesInclPad-isamp-1)/(float)(vars->numPaddedSamples-1);
    float taper = 0.5 * ( 1 + cos( M_PI*(ratio-1.0) ) );
    vars->bufferPaddedTrace[isamp] = taper * ( value2x - samplesInOut[num-isamp] );
  }
}


} // END namespace


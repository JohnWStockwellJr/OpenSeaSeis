/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"
#include "csSort.h"
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: trc_split
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_trc_split {
  struct VariableStruct {
    int inputMode;
    int numTraces;
    int oldTraceLengthMS;  // New trace length in [ms]
    int newTraceLengthS;   // New trace length in [s]
    int newTraceLengthMS;  // New trace length in [ms]
    int hdrID_subTrace;
    int numSamplesIn;
    int numSamplesOut;

    // Absolute reference time: Start times are offset from reference time
    bool isRefTime;
    int hdrID_ens;
    int hdrType_ens;
    int hdrID_time_samp1_s;
    int hdrID_time_samp1_us;
    int hdrID_endtime_us;
    //    cseis_geolib::csFlexHeader* ensValue;
    int numStartTimes;
    int* startTimes_s;

    int currentSplitTime;
    int currentSplitTrace;

    int reference_time_s;
    int residual_time_ms;
    int tred_samp1_trc1_ms_out;
    int tred_samp1_trcN_ms_out;
    int numTracesOut;
    int numLeadingZeros;
  };
  static const int MODE_NONE    = -1;
  static const int MODE_LENGTH  = 1;
  static const int MODE_NTRACES = 2;
  
  static const int INPUT_TRACE    = 11;
  static const int INPUT_ENSEMBLE = 12;
}
using namespace mod_trc_split;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_trc_split_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct*   vars = new VariableStruct();
  edef->setVariables( vars );


  vars->numTraces        = 0;
  vars->numSamplesIn     = shdr->numSamples;
  vars->numSamplesOut    = 0;
  vars->oldTraceLengthMS = (int)(vars->numSamplesIn * shdr->sampleInt + 0.5);
  vars->newTraceLengthMS = 0;
  vars->hdrID_subTrace   = -1;
  vars->hdrID_ens        = -1;
  vars->isRefTime        = false;
  //  vars->ensValue         = NULL;
  vars->hdrID_time_samp1_s  = -1;
  vars->hdrID_time_samp1_us = -1;
  vars->hdrID_endtime_us    = -1;
  vars->numStartTimes       = 0;
  vars->startTimes_s        = NULL;

  vars->inputMode           = INPUT_TRACE;

  int mode = MODE_NONE;
  float traceLengthMS_in = shdr->numSamples*shdr->sampleInt;

  if( (int)(shdr->sampleInt * 1000 + 0.5) % 1000 != 0 ) {
    writer->error("Sample interval is not an integer multiple of 1ms. This can not be handled by the current version of this module.");
  }

  //---------------------------------------------------------
  if( param->exists("ntraces") ) {
    mode = MODE_NTRACES;
    param->getInt("ntraces", &vars->numTraces);
    if( vars->numTraces <= 1 ) {
      writer->error("Specified number of traces ('ntraces') too small: %d", vars->numTraces );
    }
    vars->numSamplesOut = (int)( shdr->numSamples / vars->numTraces );
    vars->newTraceLengthMS = (int)(shdr->sampleInt * (float)vars->numSamplesOut);
  }
  //---------------------------------------------------------
  if( param->exists("mode") ) {
    string text;
    param->getString("mode", &text);
    if( !text.compare("trace") ) {
      vars->inputMode = INPUT_TRACE;
    }
    else if( !text.compare("ensemble") ) {
      vars->inputMode = INPUT_ENSEMBLE;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }
  if( vars->inputMode == INPUT_TRACE ) {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }
  else {
//    edef->setTraceSelectionMode( TRCMODE_FIXED, 4 );
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }
  //---------------------------------------------------------
  if( param->exists("length") ) {
    if( mode == MODE_NTRACES ) {
      writer->error("Only one of the two user parameters 'length' or 'ntraces' can be specified");
    }
    mode = MODE_LENGTH;
    param->getInt("length", &vars->newTraceLengthMS);
    vars->numSamplesOut = (int)(vars->newTraceLengthMS/shdr->sampleInt);
    if( vars->newTraceLengthMS < 0 || vars->numSamplesOut >= vars->numSamplesIn ) {
      writer->error("Inconsistent new trace length specified: %f. Current trace length: %f", vars->newTraceLengthMS, traceLengthMS_in );
    }
    vars->numTraces = (int)((vars->numSamplesIn-1) / vars->numSamplesOut) + 1;
  }
  //---------------------------------------------------------
  if( param->exists("ref_time") ) {
    if( mode != MODE_LENGTH ) {
      //      writer->error("Reference time only works in combination with selection mode 'length'. Mode 'ntraces' was specified");
    }
    string text;
    param->getString("ref_time", &text);
    if( !text.compare("no") ) {
      vars->isRefTime = false;
    }
    else if( !text.compare("yes") ) {
      vars->isRefTime = true;
      vars->numSamplesOut = (int)(vars->newTraceLengthMS/shdr->sampleInt + 0.5);
      if( vars->newTraceLengthMS < 0 || vars->numSamplesOut >= vars->numSamplesIn ) {
        writer->error("Inconsistent output trace length specified: %f. Current trace length: %f", vars->newTraceLengthMS, traceLengthMS_in );
      }
      vars->numTraces = (int)((vars->numSamplesIn-1) / vars->numSamplesOut) + 1;
    }
  }

  //---------------------------------------------------------
  if( param->exists("start_times") ) {
    if( mode != MODE_LENGTH ) {
      writer->line("Error: List of start times only works in combination with selection mode 'length'. Mode 'ntraces' was specified");
      env->addError();
    }
    if( vars->isRefTime ) {
      writer->line("Error: List of start times cannot be used in conjunction with UNIX reference time (user parameter 'ref_time').");
      env->addError();
    }
    int numLines = param->getNumLines("start_times");
    vars->numStartTimes = 0;
    for( int iline = 0; iline < numLines; iline++ ) {
      vars->numStartTimes += param->getNumValues("start_times",iline);
    }

    if( vars->numStartTimes <= 0 ) {
      writer->line("no start times specified for user parameter 'start_times'. List is empty.");
      env->addError();
    }
    else {
      vars->startTimes_s = new int[vars->numStartTimes];
      int counter = 0;
      for( int iline = 0; iline < numLines; iline++ ) {
        int numValues = param->getNumValues("start_times",iline);
        for( int ival = 0; ival < numValues; ival++ ) {
          param->getIntAtLine( "start_times", &vars->startTimes_s[counter++], iline, ival );
        }
      }
      vars->numTraces = vars->numStartTimes;
      csSort<int>().simpleSort( vars->startTimes_s, vars->numStartTimes );
      if( edef->isDebug() ) {
        writer->line("Start times: ");
        for( int i = 0; i < vars->numStartTimes; i++ ) {
          writer->line("Start time #%3d: %d", i+1, vars->startTimes_s[i] );
        }
      }
    }
  }

  vars->hdrID_time_samp1_s  = hdef->headerIndex("time_samp1");
  vars->hdrID_time_samp1_us = hdef->headerIndex("time_samp1_us");

  //---------------------------------------------------------
  /*  if( param->exists("ens_hdr") ) {
    if( mode != MODE_LENGTH ) {
      writer->warning("Ensemble header will be ignored for mode 'ntraces'");
    }
    else {
      string headerName;
      param->getString("ens_hdr", &headerName);
      if( !hdef->headerExists( headerName ) ) {
        writer->error("Specified ensemble header does not exist: '%s'.", headerName.c_str() );
      }
      vars->hdrType_ens = hdef->headerType( headerName );
    }
    //    if( vars->hdrType_ens != TYPE_INT && vars->hdrType_ens != TYPE_FLOAT && vars->hdrType_ens != TYPE_DOUBLE ) {
    //  writer->error("Ensemble header can only be of number type, not a string or array type.");
    // }
    } */

  if( mode == MODE_NONE ) {
    writer->error("None of the two mandatory user paramters 'length' or 'ntraces' have been specified");
  }
  if( !vars->isRefTime && vars->inputMode != INPUT_TRACE ) {
    writer->error("Ensemble input mode is currently only available in conjunction with reference time. This is a limitation of the current implementation");
  }
  vars->currentSplitTime  = 0;
  vars->currentSplitTrace = -1;

  //---------------------------------------------------------
  if( !hdef->headerExists( "sub_trcno" ) ) {
    hdef->addHeader( TYPE_INT, "sub_trcno", "Intra-trace number" );
  }
  vars->hdrID_subTrace = hdef->headerIndex( "sub_trcno" );

  shdr->numSamples = vars->numSamplesOut;
  
  writer->line("Old/new number of samples:  %d --> %d", vars->numSamplesIn, vars->numSamplesOut );
  writer->line("Number of output traces:    %d", vars->numTraces);
  writer->line("Number of residual samples: %d  (samples at end will be set to zero)", vars->numTraces*vars->numSamplesOut - vars->numSamplesIn );
  writer->line("Sample int:    %f", shdr->sampleInt );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_trc_split_(
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


  if( traceGather->numTraces() == 0 ) return;
  int numTracesIn = traceGather->numTraces();
  csTraceHeader* trcHdrOrig = traceGather->trace(0)->getTraceHeader();
  int time_samp1_s_in  = trcHdrOrig->intValue( vars->hdrID_time_samp1_s );
  int time_samp1_us_in = trcHdrOrig->intValue( vars->hdrID_time_samp1_us );

  if( time_samp1_us_in % 1000 != 0 ) {
    writer->error("TRC_SPLIT: Time of first sample in input file must be an integer number of milliseconds. This is a limitation of the current implementation");
  }
  else if( (time_samp1_us_in / 1000 ) % (int)(shdr->sampleInt + 0.5) != 0 ) {
    writer->warning("TRC_SPLIT: Time of first sample is not an integer number of sample intervals (%dms). This is a limitation of the current implementation.\nOutput data may be shifted by up to 1/2 sample interval due to that.", (int)(shdr->sampleInt+0.5));
  }
  //---------------------------------------------------------------------
  // No reference time used:
  //
  if( !vars->isRefTime ) {
    int reference_time_s     = time_samp1_s_in;
    int tref_firstSamp_ms_in = (time_samp1_s_in-reference_time_s)*1000 + time_samp1_us_in/1000;
    float* samplesIn = traceGather->trace( 0 )->getTraceSamples();
    trcHdrOrig->setIntValue( vars->hdrID_subTrace, 1 );

    if( vars->numStartTimes == 0 ) {
      traceGather->createTraces( 1, vars->numTraces-1, hdef, shdr->numSamples );
      for( int itrc = 1; itrc < vars->numTraces; itrc++ ) {
        int tref_firstSamp_currentTrace_ms = tref_firstSamp_ms_in + itrc*vars->newTraceLengthMS;
        float* samplesOut = traceGather->trace( itrc )->getTraceSamples();
        memcpy( samplesOut, &samplesIn[itrc*shdr->numSamples], 4*shdr->numSamples );
        csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
        trcHdr->copyFrom( trcHdrOrig );
        trcHdr->setIntValue( vars->hdrID_subTrace, itrc+1 );
        trcHdr->setIntValue( vars->hdrID_time_samp1_s,  (int)( tref_firstSamp_currentTrace_ms / 1000 ) + reference_time_s );
        trcHdr->setIntValue( vars->hdrID_time_samp1_us, ( tref_firstSamp_currentTrace_ms % 1000 ) * 1000 );
      }
    }
    // Start times were given. Cut input data into traces only for matching start times
    else {
      if( edef->isDebug() ) writer->line( "Start times: %d, time_samp1: %d", vars->numStartTimes, time_samp1_s_in );
      int counterFirst = 0;
      int timeFirstSamplePlus_s = time_samp1_s_in;
      if( time_samp1_us_in > 0 ) timeFirstSamplePlus_s += 1;
      while( counterFirst < vars->numStartTimes && vars->startTimes_s[counterFirst] < timeFirstSamplePlus_s ) {
        counterFirst += 1;
      }
      if( counterFirst == vars->numStartTimes ) {
        traceGather->freeAllTraces();
        return;
      }
      int timeLastSampleMinus_s = time_samp1_s_in + (int)( ((double)(vars->numSamplesIn-shdr->numSamples)*shdr->sampleInt + (double)time_samp1_us_in/1000.0 )/1000.0 + 0.5 );
      int counterLast = vars->numStartTimes - 1;
      while( counterLast >= 0 && vars->startTimes_s[counterLast] > timeLastSampleMinus_s ) {
        counterLast -= 1;
      }
      if( counterLast < 0 ) {
        traceGather->freeAllTraces();
        return;
      }
      int numActualTraces     = counterLast - counterFirst + 1;
      int indexFirstStartTime = counterFirst;

      if( numActualTraces <= 1 ) {
        traceGather->freeAllTraces();
        return;
      }
      if( numActualTraces > 1 ) {
        traceGather->createTraces( 1, numActualTraces-1, hdef, shdr->numSamples );
      }

      if( edef->isDebug() ) {
        writer->line( "Start time indices: %d - %d", counterFirst, counterLast );
        writer->line( "Number of actual output traces: %d, input traces: %d", numActualTraces, numTracesIn );
        writer->line( "First/last start time: %d %d", vars->startTimes_s[indexFirstStartTime], vars->startTimes_s[indexFirstStartTime+numActualTraces-1] );
        writer->flush();
      }

      for( int itrc = 0; itrc < numActualTraces; itrc++ ) {
        int tref_firstSamp_currentTrace_ms = tref_firstSamp_ms_in + itrc*vars->newTraceLengthMS;
        float* samplesOut = traceGather->trace( itrc )->getTraceSamples();
        tref_firstSamp_currentTrace_ms = ( vars->startTimes_s[itrc+indexFirstStartTime] - reference_time_s ) * 1000;
        int startSampleIndex = (int)( (double)tref_firstSamp_currentTrace_ms / (double)shdr->sampleInt + 0.5 );
        startSampleIndex -= (int)( (double)( time_samp1_us_in/1000 ) / (double)shdr->sampleInt + 0.5 );
        if( itrc == 0 && startSampleIndex < shdr->numSamples ) {
          for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
            samplesOut[isamp] = samplesIn[startSampleIndex+isamp];
          }
        }
        else if( startSampleIndex+shdr->numSamples < vars->numSamplesIn ) {
          memcpy( samplesOut, &samplesIn[startSampleIndex], 4*shdr->numSamples );
        }
        else {
          int endSamp = vars->numSamplesIn-startSampleIndex;
          fprintf(stdout,"End sample: %d %d %d\n\n", endSamp, vars->numSamplesIn, startSampleIndex );
          for( int isamp = 0; isamp < endSamp; isamp++ ) {
            samplesOut[isamp] = samplesIn[startSampleIndex+isamp];
          }
          for( int isamp = endSamp; isamp < shdr->numSamples; isamp++ ) {
            samplesOut[isamp] = 0.0;
          }
        }

        csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
        if( itrc != 0 ) trcHdr->copyFrom( trcHdrOrig );
        trcHdr->setIntValue( vars->hdrID_subTrace, itrc+1 );
        trcHdr->setIntValue( vars->hdrID_time_samp1_s,  (int)( tref_firstSamp_currentTrace_ms / 1000 ) + reference_time_s );
        trcHdr->setIntValue( vars->hdrID_time_samp1_us, ( tref_firstSamp_currentTrace_ms % 1000 ) * 1000 );
      }
    } // END: if( StartTimes )
  } // END: NO reference time used
  //---------------------------------------------------------------------
  //
  else if ( vars->currentSplitTrace <= 0 ) { // New input trace(s). Compute all variables and output first split traces
    for( int itrc = 1; itrc < traceGather->numTraces(); itrc++ ) {
      csTraceHeader* trcHdrTmp = traceGather->trace(itrc)->getTraceHeader();
      int time_samp1_s_tmp  = trcHdrTmp->intValue( vars->hdrID_time_samp1_s );
      int time_samp1_us_tmp = trcHdrTmp->intValue( vars->hdrID_time_samp1_us );
      if( time_samp1_s_tmp != time_samp1_s_in || time_samp1_us_tmp != time_samp1_us_in ) {
        writer->error("Absolute time of first sample not equal for all traces in input ensemble. This is currently not supported");
      }
    }

    // Absolute time of first sample, input trace [ms]
    csInt64_t t64_samp1_ms_in        = (csInt64_t)time_samp1_s_in*1000LL + (csInt64_t)time_samp1_us_in/1000LL;
    // Absolute time of last sample, input trace [ms]
    csInt64_t t64_sampN_ms_in        = t64_samp1_ms_in + (csInt64_t)vars->oldTraceLengthMS - (csInt64_t)shdr->sampleInt;
    // Absolute time of first sample, first output trace [ms]
    csInt64_t t64_samp1_trc1_ms_out  = (csInt64_t)( t64_samp1_ms_in / (csInt64_t)vars->newTraceLengthMS ) * (csInt64_t)vars->newTraceLengthMS;
    // Absolute time of first sample, last output trace [ms]
    csInt64_t t64_samp1_trcN_ms_out  = (csInt64_t)( t64_sampN_ms_in / (csInt64_t)vars->newTraceLengthMS ) * (csInt64_t)vars->newTraceLengthMS;
    // Reference time [s]. Time values are temporarily reduced by reference time for easier computation using 32bit integers:
    // Reference time is nearest absolute time in full seconds before time of first samples of first output trace
    vars->reference_time_s = (int)( t64_samp1_trc1_ms_out / 1000LL );
    vars->residual_time_ms = (int)( t64_samp1_ms_in - t64_samp1_trc1_ms_out );

    // Reduced time of first sample, first output trace [ms]
    vars->tred_samp1_trc1_ms_out = (int)( t64_samp1_trc1_ms_out - (csInt64_t)vars->reference_time_s*1000LL );
    //   = int tred_samp1_trc1_ms_out = (int)( t64_samp1_trc1_ms_out % 1000LL );
    // Reduced time of first sample, last output trace [ms]
    vars->tred_samp1_trcN_ms_out = (int)( t64_samp1_trcN_ms_out - (csInt64_t)vars->reference_time_s*1000LL );

    vars->numTracesOut = (vars->tred_samp1_trcN_ms_out - vars->tred_samp1_trc1_ms_out) / vars->newTraceLengthMS + 1;

    vars->numLeadingZeros = (int) ( (double)vars->residual_time_ms / shdr->sampleInt + 0.5 );
    if( vars->numLeadingZeros < 0 ) writer->error("Module TRC_SPLIT: Unknown program bug... numLeadingZeros: %d", vars->numLeadingZeros);

    vars->currentSplitTrace = 1;
    
    traceGather->createTraces( 0, numTracesIn, hdef, shdr->numSamples );  // Add numTracesIn traces at beginning
    for( int itrc = 0; itrc < numTracesIn; itrc++ ) {
      csTraceHeader* trcHdrIn  = traceGather->trace(numTracesIn+itrc)->getTraceHeader();
      csTraceHeader* trcHdrOut = traceGather->trace(itrc)->getTraceHeader();
      float* samplesIn  = traceGather->trace( numTracesIn+itrc )->getTraceSamples();
      float* samplesOut = traceGather->trace( itrc )->getTraceSamples();
      for( int isamp = shdr->numSamples-1; isamp >= vars->numLeadingZeros; isamp-- ) {
        samplesOut[isamp] = samplesIn[isamp-vars->numLeadingZeros];
      }
      for( int isamp = 0; isamp < vars->numLeadingZeros; isamp++ ) {
        samplesOut[isamp] = 0.0;
      }
      trcHdrOut->copyFrom( trcHdrIn );
      trcHdrOut->setIntValue( vars->hdrID_time_samp1_s,  (int)( vars->tred_samp1_trc1_ms_out / 1000 ) + vars->reference_time_s );
      trcHdrOut->setIntValue( vars->hdrID_time_samp1_us, ( vars->tred_samp1_trc1_ms_out % 1000 ) * 1000 );
      trcHdrOut->setIntValue( vars->hdrID_subTrace, vars->currentSplitTrace );
    }
    vars->currentSplitTrace += 1;


    if( edef->isDebug() ) {
      fprintf(stderr,"t64 time of first & last sample:     %lldms / %lldms   (%dms)\n", t64_samp1_ms_in, t64_sampN_ms_in, vars->newTraceLengthMS );
      fprintf(stderr,"Reference times:      %ds / %ds  --  %dms / %dms\n", vars->reference_time_s, time_samp1_s_in, vars->tred_samp1_trc1_ms_out, vars->tred_samp1_trcN_ms_out );
      fprintf(stderr,"Old/New trace length: %dms / %dms\n", vars->oldTraceLengthMS, vars->newTraceLengthMS );
      fprintf(stderr,"numTracesOut: %d\n", vars->numTracesOut );
      writer->line("Start in: %d, start1: %d, startN: %d, num traces: %d, leading zeros: %d\n",
                time_samp1_s_in, vars->tred_samp1_trc1_ms_out, vars->tred_samp1_trcN_ms_out, vars->numTracesOut, vars->numLeadingZeros );
    }
    *numTrcToKeep = numTracesIn;
  }
  else if( vars->currentSplitTrace < vars->numTracesOut ) {
    traceGather->createTraces( 0, numTracesIn, hdef, shdr->numSamples );  // Add numTracesIn traces at beginning
    int traceIndex = vars->currentSplitTrace - 1;
    for( int itrc = 0; itrc < numTracesIn; itrc++ ) {
      csTraceHeader* trcHdrIn  = traceGather->trace(numTracesIn+itrc)->getTraceHeader();
      csTraceHeader* trcHdrOut = traceGather->trace(itrc)->getTraceHeader();
      float* samplesIn  = traceGather->trace( numTracesIn+itrc )->getTraceSamples();
      float* samplesOut = traceGather->trace( itrc )->getTraceSamples();
      memcpy( samplesOut, &samplesIn[traceIndex*shdr->numSamples-vars->numLeadingZeros], 4*shdr->numSamples );
      trcHdrOut->copyFrom( trcHdrIn );
      trcHdrOut->setIntValue( vars->hdrID_subTrace, vars->currentSplitTrace );
      int tred_samp1_currentTrace_ms =  vars->tred_samp1_trc1_ms_out + traceIndex*vars->newTraceLengthMS;
      trcHdrOut->setIntValue( vars->hdrID_time_samp1_s,  (int)( tred_samp1_currentTrace_ms / 1000 ) + vars->reference_time_s );
      trcHdrOut->setIntValue( vars->hdrID_time_samp1_us, ( tred_samp1_currentTrace_ms % 1000 ) * 1000 );
    }
    vars->currentSplitTrace += 1;
    *numTrcToKeep = numTracesIn;
  }
  else { // Last split trace
    // Special treatment for last trace because potentially need to set last samples to zero
    traceGather->createTraces( 0, numTracesIn, hdef, shdr->numSamples );  // Add numTracesIn traces at beginning
    int traceIndex = vars->currentSplitTrace - 1;
    for( int itrc = 0; itrc < numTracesIn; itrc++ ) {
      csTraceHeader* trcHdrIn  = traceGather->trace(numTracesIn+itrc)->getTraceHeader();
      csTraceHeader* trcHdrOut = traceGather->trace(itrc)->getTraceHeader();
      float* samplesIn  = traceGather->trace( numTracesIn+itrc )->getTraceSamples();
      float* samplesOut = traceGather->trace( itrc )->getTraceSamples();
      int from = traceIndex*shdr->numSamples-vars->numLeadingZeros;
      int to   = MIN( vars->numSamplesIn, from+shdr->numSamples );
      memcpy( samplesOut, &samplesIn[from], 4*(to-from) );
      for( int isamp = to-from; isamp < shdr->numSamples; isamp++ ) {
        samplesOut[isamp] = 0.0;
      }
      trcHdrOut->copyFrom( trcHdrIn );
      trcHdrOut->setIntValue( vars->hdrID_time_samp1_s,  (int)( vars->tred_samp1_trcN_ms_out / 1000 ) + vars->reference_time_s );
      trcHdrOut->setIntValue( vars->hdrID_time_samp1_us, ( vars->tred_samp1_trcN_ms_out % 1000 ) * 1000 );
      trcHdrOut->setIntValue( vars->hdrID_subTrace, vars->currentSplitTrace );
    }

    vars->currentSplitTrace = -1;
    traceGather->freeTraces( numTracesIn, numTracesIn );
    *numTrcToKeep = 0;
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_trc_split_( csParamDef* pdef ) {
  pdef->setModule( "TRC_SPLIT", "Trace split", "Split trace into N traces of equal length" );

  pdef->addParam( "mode", "Mode", NUM_VALUES_FIXED);
  pdef->addValue( "trace", VALTYPE_OPTION );
  pdef->addOption( "trace", "Split one trace after the other..." );
  pdef->addOption( "ensemble", "Input ensemble, N traces. Output split trace 1 from input trace 1-N, then split trace 2 etc." );

  pdef->addParam( "ref_time", "Use reference time (1-Jan-1970 00:00:00)", NUM_VALUES_FIXED,
                  "Start times of output traces are offset by integer multiple of trace lengths relative to reference time");
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not use reference time." );
  pdef->addOption( "yes", "Make sure that trace start time is offset from reference start time by an integer multiple of trace lengths");

  pdef->addParam( "start_times", "Give start times of new traces in list", NUM_VALUES_VARIABLE,
                  "When this parameter is selected, the output trace length must be specified in parameter 'length'");
  pdef->addValue( "", VALTYPE_NUMBER, "List of absolute start times in full UNIX seconds [s]" );

  pdef->addParam( "ntraces", "Number of traces to split into", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Number of new traces" );

  pdef->addParam( "length", "Length of new traces", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Length of new traces [ms]", "...must be smaller than the current trace length" );
}



//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_trc_split_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_trc_split::VariableStruct* vars = reinterpret_cast<mod_trc_split::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_trc_split_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_trc_split::VariableStruct* vars = reinterpret_cast<mod_trc_split::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  //    if( vars->ensValue == NULL ) {
  //   delete vars->ensValue;
  //  vars->ensValue = NULL;
  // }
  if( vars->startTimes_s != NULL ) {
    delete [] vars->startTimes_s;
    vars->startTimes_s = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_trc_split_( csParamDef* pdef ) {
  params_mod_trc_split_( pdef );
}
extern "C" void _init_mod_trc_split_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_trc_split_( param, env, writer );
}
extern "C" bool _start_exec_mod_trc_split_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_trc_split_( env, writer );
}
extern "C" void _exec_mod_trc_split_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_trc_split_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_trc_split_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_trc_split_( env, writer );
}

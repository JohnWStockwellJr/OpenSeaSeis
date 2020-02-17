/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csTableAll.h"
#include "csSort.h"
#include "csTimer.h"
#include "csVector.h"
#include "csGeolibUtils.h"
#include "csInterpolation.h"
#include <cmath>
#include <cstring>
#include <iostream>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SELECT_TIME
 *
 * @author Bjorn Olofsson
 * @date   2007
 *
 * Bug fixes
 *  2009-03-25 - Correctly zero out new samples
 */
namespace mod_select_time {
  struct ShotAttr {
    csInt64_t time_us;
    int       id;
  };
  struct VariableStruct {
    int numSamplesOrig;
    int startSamp;
    int endSamp;
    int mode;
    int hdrID_time_samp1_s;
    int hdrID_time_samp1_us;
    int hdrID_key;
    cseis_geolib::type_t hdrType_key;
    int numLocations;
    double* keyValues;
    int* numLines;
    double** startTimes;
    double** endTimes;
    bool isTable;
    bool isDelTrace;
    float percentDelTrace_ratio;
    bool doTrim;
    bool isShiftTrace;
    cseis_geolib::csVector<ShotAttr>* shotAttrList;
    float* shotBuffer;
    float* shotBufferTemp;
    csInt64_t startTimeCurrentShot_us;
    csInt64_t endTimeCurrentShot_us;
    int idCurrentShot;
    cseis_geolib::csInterpolation* interpol;
    bool isNonZeroShotInMemory;
    cseis_system::csTraceGather* traceGatherKeep;
    int indexCurrentInternalTrace;
    //    int overlapMethod;
    int traceCounter;
    int numSampEdgeEffectReduction;
  };
  static int const MODE_ABSOLUTE = 1;
  static int const MODE_RELATIVE = 2;
  static int const MODE_SHOTS    = 3;
  static int const OVERLAP_FIRST_TRACE = 20;
  static int const OVERLAP_STACK       = 21;
}
using namespace mod_select_time;

void addShot( VariableStruct* vars, cseis_system::csTraceGather* traceGather, csInt64_t& startTimeCurrentTrace_us, csInt64_t& endTimeCurrentTrace_us, cseis_system::csSuperHeader const* shdr );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_select_time_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->isShiftTrace = false;
  vars->startSamp = 0;
  vars->endSamp   = 0;
  vars->mode      = MODE_RELATIVE;
  vars->hdrID_time_samp1_s = -1;
  vars->hdrID_time_samp1_us = -1;
  vars->hdrID_key       = -1;
  vars->hdrType_key     = TYPE_UNKNOWN;
  vars->startTimes      = NULL;
  vars->endTimes        = NULL;
  vars->numLines        = NULL;
  vars->numLocations    = 0;
  vars->isTable         = false;
  vars->isDelTrace      = false;
  vars->percentDelTrace_ratio = 1.0;
  vars->numSamplesOrig = shdr->numSamples;
  vars->doTrim = false;
  vars->shotAttrList = NULL;
  vars->shotBuffer = NULL;
  vars->shotBufferTemp = NULL;
  vars->startTimeCurrentShot_us = 0LL;
  vars->endTimeCurrentShot_us   = 0LL;
  vars->idCurrentShot = 0;
  vars->interpol  = NULL;
  vars->isNonZeroShotInMemory = false;
  vars->traceGatherKeep = NULL;
  //  vars->overlapMethod = mod_concatenate::OVERLAP_FIRST_TRACE;
  vars->traceCounter = 0;
  vars->numSampEdgeEffectReduction = 2;

  std::string text;
  //---------------------------------------------------------
  if( param->exists("mode") ) {
    param->getString( "mode", &text );
    if( !text.compare( "relative" ) ) {
      vars->mode = MODE_RELATIVE;
    }
    else if( !text.compare( "absolute" ) ) {
      vars->mode = MODE_ABSOLUTE;
    }
    else if( !text.compare( "extract_shots" ) ) {
      vars->mode = MODE_SHOTS;
    }
    else {
      writer->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }
  //---------------------------------------------------------
  if( param->exists("free_mem") ) {
    param->getString( "free_mem", &text );
    if( !text.compare( "yes" ) ) {
      vars->doTrim = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->doTrim = false;
    }
    else {
      writer->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }
  if( param->exists("shift") ) {
    param->getString( "shift", &text );
    if( !text.compare( "yes" ) ) {
      vars->isShiftTrace = true;
      if( vars->mode == MODE_ABSOLUTE ) {
        writer->error("Shift trace cannot be applied for absolute time selection");
      }
    }
    else if( !text.compare( "no" ) ) {
      vars->isShiftTrace = false;
    }
    else {
      writer->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }
  //---------------------------------------------------------
  bool isTimeDomain = true;

  if( param->exists("domain") ) {
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
  if( param->exists("del_traces") ) {
    param->getString( "del_traces", &text );
    if( !text.compare( "yes" ) ) {
      vars->isDelTrace = true;
      if( param->getNumValues("del_traces") > 1 ) {
        float percentDelTrace;
        param->getFloat( "del_traces", &percentDelTrace, 1 );
        if( percentDelTrace <= 0 || percentDelTrace > 100 ) {
          writer->line("Error in user parameter 'del_traces': Percent threshold of trace must be in the range of ]0,100]. Given value: %f", percentDelTrace);
          env->addError();
        }
        vars->percentDelTrace_ratio = percentDelTrace / 100.0;
      }
    }
    else if( !text.compare( "no" ) ) {
      vars->isDelTrace = false;
    }
    else {
      writer->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  //---------------------------------------------------------
  if( vars->mode == MODE_SHOTS ) {
    param->getInt("nsamples",&shdr->numSamples);
    // The following restriction makes sure that one input trace will at maximum spread over two output traces, not more.
    if( shdr->numSamples < vars->numSamplesOrig ) writer->error("Unsupported case: Number of output samples (%d) must be greater or equal the number of samples in input data (%d).", shdr->numSamples, vars->numSamplesOrig);
    param->getString( "extract_shots", &text, 0 );
    int numCoefficients = 8;
    
    int numValues = param->getNumValues( "extract_shots" );
    vars->numSampEdgeEffectReduction = 2;
    if( numValues > 1 ) {
      param->getInt( "extract_shots", &vars->numSampEdgeEffectReduction, 1 );
    }

    FILE* fin = fopen(text.c_str(),"r");
    if( fin == NULL ) writer->error("Error opening ASCII file which contains the shot times: '%s'", text.c_str());
    char buffer[132];

    vars->shotAttrList = new cseis_geolib::csVector<ShotAttr>();
    csInt64_t shot_time_prev = 0;
    while( fgets(buffer,132,fin) != NULL ) {
      ShotAttr shotAttr;
      int shot_time_s;
      int shot_time_us;
      if( (int)strlen(buffer) < 4 ) continue;
      sscanf(buffer,"%d %d %d", &shotAttr.id, &shot_time_s, &shot_time_us);
      shotAttr.time_us = (csInt64_t)shot_time_s*1000000LL + (csInt64_t)shot_time_us;
      if( shotAttr.time_us < shot_time_prev ) writer->error("Shot times not in ascending order");
      vars->shotAttrList->insertEnd( shotAttr );
      shot_time_prev = shotAttr.time_us;
      if( edef->isDebug() ) cout << "Shot time " << shot_time_prev << " " << shotAttr.id << "\n";
    }
    int numTimes = vars->shotAttrList->size();
    if( numTimes == 0 ) writer->error("No shot times read in from file %s: Is shot time file empty?\n", text.c_str() );
    int shotTime1_s = (int)vars->shotAttrList->at(0).time_us/1000000LL;
    int shotTime2_s = (int)vars->shotAttrList->at(numTimes-1).time_us/1000000LL;
    writer->line("Read in %d shot times from %s to %s\n", numTimes, csGeolibUtils::UNIXsec2dateString(shotTime1_s).c_str(), csGeolibUtils::UNIXsec2dateString(shotTime2_s).c_str() );
    fclose(fin);
    vars->shotBuffer     = new float[shdr->numSamples];
    vars->shotBufferTemp = new float[vars->numSamplesOrig];
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      vars->shotBuffer[isamp] = 0;
    }
    vars->interpol = new csInterpolation( vars->numSamplesOrig, shdr->sampleInt, numCoefficients );
    vars->idCurrentShot = vars->shotAttrList->at(0).id;
    vars->startTimeCurrentShot_us = vars->shotAttrList->at(0).time_us;

    vars->endTimeCurrentShot_us   = vars->startTimeCurrentShot_us + (csInt64_t)(shdr->numSamples-1) * (csInt64_t)shdr->sampleInt*1000LL;
    vars->isNonZeroShotInMemory = false;
    vars->traceGatherKeep = new cseis_system::csTraceGather( hdef );

    vars->hdrID_time_samp1_s  = hdef->headerIndex( HDR_TIME_SAMP1.name );
    vars->hdrID_time_samp1_us = hdef->headerIndex( HDR_TIME_SAMP1_US.name );
    if( !hdef->headerExists( HDR_FFID.name ) ) hdef->addStandardHeader( HDR_FFID.name );
    vars->hdrID_key = hdef->headerIndex( HDR_FFID.name );

    vars->indexCurrentInternalTrace = 0;
  }
  else if( param->exists("table") ) {
    if( vars->isShiftTrace ) {
      writer->error("Trace shift is not possible to perform when specifying start/end times in a table");
    }

    csSort<double> sortObj;

    vars->isTable = true;
    param->getString( "table", &text );
    csTimer timer;
    timer.start();
    csTableAll table( csTableAll::TABLE_TYPE_DUPLICATE_KEYS );
    try {
      table.initialize( text );
      int numKeys   = table.numKeys();
      int numValues = table.numValues();
      if( numKeys != 1 ) writer->error("Number of key columns in input table file: %d. Only 1 key column is currently supported", numKeys);
      if( numValues != 2 ) writer->error("Number of value columns in input table file: %d. Expected 2 value columns, i.e. 'start' and 'end'", numValues);
      if( table.valueName(0).compare("start") ) {
        writer->warning("Name of first table value is '%s'. Correct/expected value name should be 'start'.", table.valueName(0).c_str() );
      }
      if( table.valueName(1).compare("end") ) {
        writer->warning("Name of second table value is '%s'. Correct/expected value name should be 'end'.", table.valueName(1).c_str() );
      }
      vars->hdrID_key   = hdef->headerIndex( table.keyName(0) );
      vars->hdrType_key = hdef->headerType( table.keyName(0) );

      table.readTableContents();

      if( edef->isDebug() ) table.dump();
      int id_start = 0;
      int id_end   = 1;
      vars->numLocations = table.numLocations();
      vars->keyValues  = new double[vars->numLocations];
      vars->startTimes = new double*[vars->numLocations];
      vars->endTimes   = new double*[vars->numLocations];
      vars->numLines   = new int[vars->numLocations];
      int* indexList   = NULL;
      int maxLines = 0;
      double sampleIntInSeconds = (double)shdr->sampleInt / 1000.0;
      for( int iloc = 0; iloc < vars->numLocations; iloc++ ) {
        vars->keyValues[iloc] = table.getKeyValue( iloc, 0 );
        TableValueList const* tvl = table.getValues( &vars->keyValues[iloc] );
        int numLines = tvl->numLines();
        if( numLines > maxLines ) {
          if( indexList ) delete [] indexList;
          indexList = new int[numLines];
          maxLines  = numLines;
        }
        vars->numLines[iloc] = numLines+1;
        double* startTimes = new double[numLines+1];
        double* endTimes   = new double[numLines+1];
        for( int iline = 0; iline < numLines; iline++ ) {
          startTimes[iline] = tvl->get(id_start,iline);
          indexList[iline]  = iline;
        }
        sortObj.simpleSortIndex( startTimes, numLines, indexList );

        for( int iline = 0; iline < numLines; iline++ ) {
          endTimes[iline] = tvl->get(id_end,indexList[iline]);
        }

        if( edef->isDebug() ) {
          writer->line("Good ranges:");
          for( int iline = 0; iline < numLines; iline++ ) {
            writer->line("  #%-3d:  %f  - %f", iline, startTimes[iline], endTimes[iline]);
          }
        }

        // NOW, convert GOOD selection times into BAD selection times, for easier selection later on
        for( int iline = numLines; iline > 0; iline-- ) {
          endTimes[iline]   = startTimes[iline] - sampleIntInSeconds;          
          startTimes[iline] = endTimes[iline-1] + sampleIntInSeconds;
        }
        endTimes[0]        = startTimes[0] - sampleIntInSeconds;
        startTimes[0]      = 0;
        endTimes[numLines] = 1e40;
        vars->startTimes[iloc] = startTimes;
        vars->endTimes[iloc]   = endTimes;

        if( edef->isDebug() ) {
          writer->line("Bad ranges:");
          for( int iline = 0; iline < numLines+1; iline++ ) {
            writer->line("  #%-3d:  %f  - %f", iline, startTimes[iline], endTimes[iline]);
          }
        }
      }
      if( indexList ) delete [] indexList;

      writer->line("Input table '%s'", table.tableName().c_str() ); 
      writer->line("  Number of key locations:       %d", vars->numLocations);
      writer->line("  CPU time for reading in table: %12.6f seconds\n\n", timer.getElapsedTime() );
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s':\n", text.c_str(), exc.getMessage() );
    }
  }
  else { // No table given
    // No table given
    //---------------------------------------------------------
    float startTime = 0.0;
    float endTime = 0.0;
    
    if( isTimeDomain ) {
      if( vars->mode == MODE_ABSOLUTE ) {
        vars->numLocations = 1;
        vars->numLines = new int[vars->numLocations];
        vars->numLines[0] = param->getNumValues( "start" );
        if( vars->numLines[0] == 0 ) writer->error("No start values given");
        else if( param->getNumValues( "end" ) != vars->numLines[0] ) writer->error("Unequal start times (%d) and end times (%d) given", vars->numLines[0],param->getNumValues( "end" ));
        vars->startTimes = new double*[vars->numLocations];
        vars->endTimes   = new double*[vars->numLocations];
        vars->startTimes[0] = new double[vars->numLines[0]];
        vars->endTimes[0]   = new double[vars->numLines[0]];
        for( int iline = 0; iline < vars->numLines[0]; iline++ ) {
          double startTime;
          double endTime;
          param->getDouble( "start", &startTime, iline );
          param->getDouble( "end", &endTime, iline );
          vars->startTimes[0][iline] = startTime * 1000.0; // Convert to milliseconds
          vars->endTimes[0][iline]   = endTime * 1000.0; // Convert to milliseconds
          if( edef->isDebug() ) {
            writer->line("Start/end times[%d]: %15.3fs  %15.3fs", iline, startTime, endTime);
          }
        }
      }
      else {
        param->getFloat( "start", &startTime );
        endTime = (float)(shdr->numSamples-1) * shdr->sampleInt;
        if( param->exists("end") ) {
          param->getFloat( "end", &endTime );
        }
        if( startTime < 0.0 ) writer->error("Start time (%f) needs to be greater or equal to 0.0.", startTime);
        if( startTime > endTime ) writer->error("Start time (%f) needs to be smaller than end time (%f).", startTime, endTime);
        if( startTime >= shdr->numSamples*shdr->sampleInt ) writer->error("Start time (%f) exceeds input trace end time (%f).", startTime, shdr->numSamples*shdr->sampleInt);
        vars->startSamp = (int)(startTime / shdr->sampleInt);  // All in milliseconds
        vars->endSamp   = (int)(endTime / shdr->sampleInt);
      }
    }
    else { // Sample domain
      if( vars->mode == MODE_ABSOLUTE ) {
        writer->error("Absolute start/end times must be provided in seconds since 1-1-1970, not in samples");
      }
      //
      // NOTE: User input is '1' for first sample. Internally, '0' is used!!
      //
      param->getInt( "start", &vars->startSamp );
      vars->endSamp = shdr->numSamples;
      if( param->exists("end") ) {
        param->getInt( "end", &vars->endSamp );
      }
      if( vars->startSamp < 1 ) writer->error("Start sample (%d) needs to be greater or equal to 1.", vars->startSamp);
      if( vars->startSamp > vars->endSamp ) writer->error("Start sample (%d) needs to be smaller than end sample (%d).", vars->startSamp, vars->endSamp);
      vars->startSamp -= 1;   // see note above..
      vars->endSamp   -= 1;
      startTime = (float)vars->startSamp * shdr->sampleInt;
      endTime   = (float)vars->endSamp * shdr->sampleInt;
    }
    //    if( vars->endSamp > shdr->numSamples ) writer->error("Specified end sample/time (%d/%.2f) is larger than number of input samples (%d)",
    //                                                vars->endSamp, vars->endSamp*shdr->sampleInt, shdr->numSamples);

    // Set new number of samples
    if( vars->mode == MODE_RELATIVE ) {
      if( !vars->isShiftTrace ) {
        shdr->numSamples = vars->endSamp + 1;
      }
      else {
        shdr->numSamples = vars->endSamp - vars->startSamp + 1;
      }
    }
    else {
      vars->hdrID_time_samp1_s  = hdef->headerIndex( HDR_TIME_SAMP1.name );
      vars->hdrID_time_samp1_us = hdef->headerIndex( HDR_TIME_SAMP1_US.name );
    }

    if( edef->isDebug() ) {
      writer->line("time1: %f, time2: %f, sample1: %d, sample2: %d, sampInt: %f\n", startTime, endTime, vars->startSamp, vars->endSamp, shdr->sampleInt );
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_select_time_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;


  vars->traceCounter += traceGather->numTraces();
  //  cout << "************ READ IN TRACE " << vars->traceCounter << "\n";
  if( edef->isLastCall() && vars->mode == mod_select_time::MODE_SHOTS && vars->isNonZeroShotInMemory ) {
    csInt64_t dummy1;
    csInt64_t dummy2;
    addShot( vars, traceGather, dummy1, dummy2, shdr );
    vars->traceGatherKeep->freeAllTraces();
    vars->isNonZeroShotInMemory = false; // Redundant
    vars->indexCurrentInternalTrace = 0;
    return;
  }

  if( traceGather->numTraces() == 0 ) {
    return;
  }
  csTrace* trace = traceGather->trace(0);
  float* samples = trace->getTraceSamples();

  if( edef->isDebug() ) {
    writer->line("SELECT_TIME: Number of samples in trace: %d, Super header: %d", trace->numSamples(), shdr->numSamples );
  }

  if( vars->isTable ) {
    double sampleIntInSeconds = (double)shdr->sampleInt / 1000.0;

    csTraceHeader* trcHdr = trace->getTraceHeader();
    int time_samp1_s  = trcHdr->intValue( vars->hdrID_time_samp1_s );
    int time_samp1_us = trcHdr->intValue( vars->hdrID_time_samp1_us );

    double keyValue = trace->getTraceHeader()->doubleValue( vars->hdrID_key );
    int location = 0;    
    while( location < vars->numLocations && vars->keyValues[location] != keyValue ) {
      location += 1;
    }
    if( location >= vars->numLocations ) {
      writer->warning("Key value %d not found in input table. No time selection applied.", keyValue);
    }
    else {
      double time_firstSamp = (double)time_samp1_s + ((double)time_samp1_us)/1000000.0;   // Absolute time of first sample, in seconds since 01-01-1970
      double time_lastSamp  = time_firstSamp + (double)((shdr->numSamples-1)*shdr->sampleInt)/1000.0;    // Absolute time of last sample

      int iline = 0;
      int numLines = vars->numLines[location];
      double* startTimes = vars->startTimes[location];
      double* endTimes   = vars->endTimes[location];

      int numDeletedSamples = 0;

      if( edef->isDebug() ) writer->line("--- KEY %f ---", keyValue );
      while( iline < numLines && endTimes[iline] < time_firstSamp ) iline += 1;
      while( iline < numLines && startTimes[iline] <= time_lastSamp ) {
        if( startTimes[iline] <= time_firstSamp && endTimes[iline] >= time_lastSamp ) {
          if( edef->isDebug() ) writer->line("   ZERO'ing WHOLE TRACE!");
          memset( samples, 0, shdr->numSamples*4 );
          numDeletedSamples = shdr->numSamples;
          break;
        }
        int startSamp = (int)( ( MAX( startTimes[iline], time_firstSamp ) - time_firstSamp ) / sampleIntInSeconds + 0.5 );
        int endSamp   = MIN( (int)( ( MIN( endTimes[iline], time_lastSamp ) - time_firstSamp ) / sampleIntInSeconds + 0.5 ), shdr->numSamples-1 );
        if( edef->isDebug() ) writer->line("   ZERO'ing samples: %10d to %10d,   %f  to %f  (%f  %f)", startSamp, endSamp, startTimes[iline], endTimes[iline], time_firstSamp, time_lastSamp );
        memset( &samples[startSamp], 0, (endSamp-startSamp+1)*4 );
        numDeletedSamples += (endSamp-startSamp+1);
        iline += 1;
      }
      if( vars->isDelTrace && (float)numDeletedSamples/(float)shdr->numSamples >= vars->percentDelTrace_ratio ) {
        traceGather->freeTrace(0);
      }
      return;
    }
  }
  //--------------------------------------------------------------------------------
  // Splice data into new shot records, based on shot times read in from ASCII file
  else if( vars->mode == MODE_SHOTS ) {
    if( vars->shotAttrList->size() == 0 ) { // No more shots to splice --> Remove traces and return
      traceGather->freeAllTraces();
      return;
    }

    traceGather->moveTraceTo( 0, vars->traceGatherKeep );  // Move all traces to internal buffer for further processing
    if( vars->indexCurrentInternalTrace >= vars->traceGatherKeep->numTraces() ) writer->error("Inconsistent trace index: %d %d", vars->indexCurrentInternalTrace, vars->traceGatherKeep->numTraces() );
    trace = vars->traceGatherKeep->trace( vars->indexCurrentInternalTrace );

    csTraceHeader* trcHdr = trace->getTraceHeader();
    csInt64_t startTimeCurrentTrace_us = (csInt64_t)trcHdr->intValue( vars->hdrID_time_samp1_s )*1000000LL + (csInt64_t)trcHdr->intValue( vars->hdrID_time_samp1_us );
    csInt64_t endTimeCurrentTrace_us   = startTimeCurrentTrace_us + (csInt64_t)( (vars->numSamplesOrig-1) * shdr->sampleInt )*1000LL;
    csInt64_t sampleInt_us = (csInt64_t)round( (double)shdr->sampleInt * 1000.0 );
    
    //    cout << "  Num internal traces: " << vars->traceGatherKeep->numTraces() << " / shot in memory? " << vars->isNonZeroShotInMemory << "\n";
    //    cout << "A Trace time: " <<  startTimeCurrentTrace_us << "  " << endTimeCurrentTrace_us << " / " Shot time:  " <<  vars->startTimeCurrentShot_us << "  " << vars->endTimeCurrentShot_us << "\n";

    //--------------------------------------------------------------------------------
    // (2) Current trace time range is earlier than currently requested shot time --> throw trace away, move on to next trace if available. Return to trace flow if no more trace available
    while( endTimeCurrentTrace_us < vars->startTimeCurrentShot_us && vars->traceGatherKeep->numTraces() > 0 ) {
      //      cout << " ---- (2) trace end time earlier than requested shot ------\n";
      vars->traceGatherKeep->freeTrace(0);
      if( vars->traceGatherKeep->numTraces() == 0 ) {
        return;
      }
      vars->indexCurrentInternalTrace = 0;
      trace  = vars->traceGatherKeep->trace(vars->indexCurrentInternalTrace);
      trcHdr = trace->getTraceHeader();
      startTimeCurrentTrace_us = (csInt64_t)trcHdr->intValue( vars->hdrID_time_samp1_s )*1000000LL + (csInt64_t)trcHdr->intValue( vars->hdrID_time_samp1_us );
      endTimeCurrentTrace_us   = startTimeCurrentTrace_us + (csInt64_t)( (vars->numSamplesOrig-1) * shdr->sampleInt )*1000LL;
    }
    samples = trace->getTraceSamples();

    //    cout << "B Trace time: " <<  startTimeCurrentTrace_us << "  " << endTimeCurrentTrace_us << " / " Shot time:  " <<  vars->startTimeCurrentShot_us << "  " << vars->endTimeCurrentShot_us << "\n";

    // (3) Current trace start time is later than currently requested shot end time --> Output currently buffered shot to trace flow if there is one, and prepare next shot
    while( startTimeCurrentTrace_us > vars->endTimeCurrentShot_us && vars->shotAttrList->size() > 0 ) {
      //      cout << " ---- (3) Trace start time later than requested shot ------\n";
      // Output currently buffered shot to internally stored trace and move this trace out to trace flow
      if( vars->isNonZeroShotInMemory ) {
        addShot( vars, traceGather, startTimeCurrentTrace_us, endTimeCurrentTrace_us, shdr );
      }
      vars->shotAttrList->remove(0); // Remove currently requested shot since this has now been processed. The current and all following traces have later times.
      if( vars->shotAttrList->size() == 0 ) {
        vars->traceGatherKeep->freeAllTraces();
        return; // Nothing more to do. All requested shots have been processed
      }
      vars->idCurrentShot = vars->shotAttrList->at(0).id;
      vars->startTimeCurrentShot_us = vars->shotAttrList->at(0).time_us;

      vars->endTimeCurrentShot_us   = vars->startTimeCurrentShot_us + (csInt64_t)(shdr->numSamples-1) * (csInt64_t)shdr->sampleInt*1000LL;
    } // END: Check if current trace start time is later than requested shot time

    //    cout << "C Trace time: " <<  startTimeCurrentTrace_us << "  " << endTimeCurrentTrace_us << " / " Shot time:  " <<  vars->startTimeCurrentShot_us << "  " << vars->endTimeCurrentShot_us << "\n";

    // (2) Current trace time range overlaps with currently requested shot --> Store overlapping samples in shotBuffer, iterate through all internally stored traces, then return to trace flow.
    // Add shot to output trace flow if needed.
    if( startTimeCurrentTrace_us <= vars->endTimeCurrentShot_us && endTimeCurrentTrace_us >= vars->startTimeCurrentShot_us ) {
      bool loopAgain;
      int shift_us = 0;
      // Loop over internally stored traces (traceGatherKeep):
      do {
        loopAgain = false;
        // Apply sub-sample static shift to align current trace samples to full sample index relative to requested shot time
        shift_us = -( vars->startTimeCurrentShot_us - startTimeCurrentTrace_us ) % sampleInt_us;
        if( shift_us > 0 ) shift_us = shift_us - sampleInt_us;
        float shift_ms = (float)( (double)shift_us/1000.0 );  // Subsample bulk shift in [ms]

        // Shift current trace by shift_us --> shotBufferTemp
        vars->interpol->static_shift( shift_ms, samples, vars->shotBufferTemp ); // Shift trace samples so that they match up on a full sample with requested shot time
        startTimeCurrentTrace_us -= shift_us;
        endTimeCurrentTrace_us   -= shift_us; 
        //        cout << " ---- OVERLAP FOUND " <<  vars->shotAttrList->size() << " " << vars->traceGatherKeep->numTraces()  << " / Shift: " << shift_us << "\n";

        // Reduce overlap range by one sample if overlap zone touches edge of input trace. The reason is that the static shift will have distorted the first or last samples depending on the shift.
        csInt64_t reduceStartSample_us = 0LL; // No reduction needed if time shift was exactly 0
        if( shift_us < 0 ) { // Only reduce trace length when needed.
          reduceStartSample_us    = vars->numSampEdgeEffectReduction * sampleInt_us;
          endTimeCurrentTrace_us -= vars->numSampEdgeEffectReduction * sampleInt_us;
        }

        // Determine overlapping time range
        csInt64_t startTime_us = std::max( vars->startTimeCurrentShot_us, startTimeCurrentTrace_us + reduceStartSample_us );
        csInt64_t endTime_us   = std::min( vars->endTimeCurrentShot_us,   endTimeCurrentTrace_us   );
        int numSampToCopy = (int)( (endTime_us - startTime_us)/sampleInt_us ) + 1;
        int samp1In  = (int)( (startTime_us - startTimeCurrentTrace_us)/sampleInt_us );
        int samp1Out = (int)( (startTime_us - vars->startTimeCurrentShot_us)/sampleInt_us );
        //    cout << "D Trace time: " <<  startTimeCurrentTrace_us << "  " << endTimeCurrentTrace_us << " / " Shot time:  " <<  vars->startTimeCurrentShot_us << "  " << vars->endTimeCurrentShot_us << "\n";
        //    cout << " ..set times: " << startTime_us << " " << endTime_us << " / Copy samples: " << samp1In << " " <<  samp1Out << " " <<  numSampToCopy << "\n";

        if( samp1In < 0 || samp1Out < 0 ) writer->error("SELECT_TIME: Program bug in exec phase. Wrong computed sample index: %d %d\n", samp1In, samp1Out );
        if( samp1In >= vars->numSamplesOrig || samp1Out >= shdr->numSamples ) writer->error("Incorrect computed sample index: samp_in:%d, samp_out:%d (nsamp_in:%d, nsamp_out: %d)\nCheck time stamps (headers time_samp1,time_samp1_us) in input data\n", samp1In, samp1Out, vars->numSamplesOrig, shdr->numSamples );
        
        if( numSampToCopy > 0 ) memcpy( &vars->shotBuffer[samp1Out], &vars->shotBufferTemp[samp1In], numSampToCopy*sizeof(float) );
        vars->isNonZeroShotInMemory = true;
        vars->indexCurrentInternalTrace += 1; // Look at next trace to find missing samples for current on next iteration

        // There are more internally stored traces AND the current shot time range exceeds the current internally stored trace --> iterate through next internally stored trace
        if( vars->indexCurrentInternalTrace < vars->traceGatherKeep->numTraces() && vars->endTimeCurrentShot_us > endTimeCurrentTrace_us ) {
          trace   = vars->traceGatherKeep->trace( vars->indexCurrentInternalTrace );
          trcHdr  = trace->getTraceHeader();
          samples = trace->getTraceSamples();
          startTimeCurrentTrace_us = (csInt64_t)trcHdr->intValue( vars->hdrID_time_samp1_s )*1000000LL + (csInt64_t)trcHdr->intValue( vars->hdrID_time_samp1_us );
          endTimeCurrentTrace_us   = startTimeCurrentTrace_us + (csInt64_t)( (vars->numSamplesOrig-1) * shdr->sampleInt )*1000LL;
          loopAgain = ( vars->endTimeCurrentShot_us > startTimeCurrentTrace_us );
        }
      } while( loopAgain ); // Loop over internally stored traces to seek for more samples for currently requested shot

      if( vars->endTimeCurrentShot_us <= endTimeCurrentTrace_us ) { // Shot has been completely copied --> release this (shot) trace to trace flow, move to next shot
        addShot( vars, traceGather, startTimeCurrentTrace_us, endTimeCurrentTrace_us, shdr );
        vars->shotAttrList->remove(0); // Remove currently requested shot since this has now been processed. This and all following traces have later times.
        if( vars->shotAttrList->size() == 0 ) {
          vars->traceGatherKeep->freeAllTraces();
          return; // Nothing more to do. All requested shots have been processed
        }
        vars->idCurrentShot = vars->shotAttrList->at(0).id;
        vars->startTimeCurrentShot_us = vars->shotAttrList->at(0).time_us;
        vars->endTimeCurrentShot_us   = vars->startTimeCurrentShot_us + (csInt64_t)(shdr->numSamples-1) * (csInt64_t)shdr->sampleInt*1000L;
        vars->indexCurrentInternalTrace = 0; // Start looking at earliest trace for this shot on next iteration
        //        cout << " ADDED SHOT, numtraces in keep gather: " << vars->traceGatherKeep->numTraces() << "\n";
      }
    }
    edef->setTracesAreWaiting( vars->traceGatherKeep->numTraces() > 0 );
    return;
  }
  else if( vars->mode == MODE_ABSOLUTE ) {
    csTraceHeader* trcHdr = trace->getTraceHeader();
    double startTimeCurrent_ms = 1000.0 * (double)trcHdr->intValue( vars->hdrID_time_samp1_s ) + (double)trcHdr->intValue( vars->hdrID_time_samp1_us )/1000.0;
    double endTimeCurrent_ms   = (double)startTimeCurrent_ms + (double)( (shdr->numSamples-1) * shdr->sampleInt );
    if( edef->isDebug() ) {
      writer->line("Start/end times: %15.3fs  %15.3fs", startTimeCurrent_ms/1000, endTimeCurrent_ms/1000);
    }
    for( int iline = 0; iline < vars->numLines[0]; iline++ ) {
      if( endTimeCurrent_ms <= vars->startTimes[0][iline] || startTimeCurrent_ms >= vars->endTimes[0][iline] ) {
        continue;
      }
      else {
        return;  // At least some portion of trace is within specified time window
        break;
      }
    }
    traceGather->freeTrace(0);
    return;
  }
  else if( !vars->isShiftTrace ) {
    //
    // Zero out trace samples outside the selected time gate
    //
    if( vars->startSamp > 0 ) {
      memset( samples, 0, vars->startSamp*4 );
    }
    if( vars->numSamplesOrig < shdr->numSamples ) {
      memset( &samples[vars->numSamplesOrig], 0, (shdr->numSamples - vars->numSamplesOrig)*4 );
    }
  }
  else {
    for( int isamp = vars->startSamp; isamp <= vars->endSamp; isamp++ ) {
      int indexTo = isamp - vars->startSamp;
      samples[indexTo] = samples[isamp];
    }
  }
  // Free any left-over memory
  if( vars->doTrim ) {
    trace->trim();
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_select_time_( csParamDef* pdef ) {
  pdef->setModule( "SELECT_TIME", "Select time of interest", "Samples outside of the chosen time interval are either set to zero, or removed" );

  pdef->addDoc("For active data processing, only the relative time is usually of interest. For continuous data processing, absolute time selection is usually required.");

  pdef->addParam( "domain", "Time or sample domain", NUM_VALUES_FIXED );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Window is specified in time [ms]" );
  pdef->addOption( "sample", "Window is specified in samples (1 for first sample)" );

  pdef->addParam( "mode", "Time selection mode", NUM_VALUES_FIXED );
  pdef->addValue( "relative", VALTYPE_OPTION );
  pdef->addOption( "relative", "Time window(s) are specified in relative time or sample index" );
  pdef->addOption( "absolute", "Time window(s) are specified as absolute times" );
  pdef->addOption( "extract_shots", "Extract shot records by absolute shot time. Read in shot times from table (single column). Specify 'nsamples' for number of output samples" );

  pdef->addParam( "start", "List of start times/samples", NUM_VALUES_FIXED, "Depends on 'domain' parameter" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Start times/samples" );

  pdef->addParam( "end", "List of end times/samples", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "End times/samples" );

  pdef->addParam( "table", "Table with window start/end times or samples", NUM_VALUES_FIXED, "For absolute time mode, columns 'start_time' and 'end_time', given in [ms] since 01-Jan-1970, must be specified" );
  pdef->addValue( "", VALTYPE_STRING, "Full path name of table containing window start/end times");

  pdef->addParam( "extract_shots", "Parameters for shot extraction", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Full path name of ASCII file containing shot ID and times in the form of two integer numbers: SHOTID TIME[s] TIME_FRACTION[us]");
  pdef->addValue( "2", VALTYPE_NUMBER, "Number of samples at start and end of each input trace to skip in order to reduce edge effects from internal sinc interpolation");
  //  pdef->addValue( "8", VALTYPE_NUMBER, "Number of coefficients for internal sinc interpolation (needed for potential relative static shifts between contributions from more than one input trace to one output trace)");


  pdef->addParam( "del_traces", "Delete trace if more then specified amount of trace has been de-selected", NUM_VALUES_VARIABLE,
                  "This option really only makes sense when used in conjunction with an absolute time selection" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not delete any traces. Set de-selected samples to zero." );
  pdef->addOption( "yes", "Delete traces if more than X percent of trace has been de-selected. Otherwise set de-selected samples to zero." );
  pdef->addValue( "100", VALTYPE_NUMBER, "Threshold of removed data that triggers trace deletion, given in percent [%].");

  pdef->addParam( "free_mem", "Free unused data in case output trace is shorter than input trace", NUM_VALUES_VARIABLE,
                  "This can be used to boost memory performance" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "By default, the system decides whether to free any unused memory immediately or reuse it for later." );
  pdef->addOption( "yes", "Free extra memory." );

  pdef->addParam( "shift", "Shift first life sample to start of trace?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not shift trace. This means that output trace length will be (endSample+1) samples long" );
  pdef->addOption( "yes", "Shift first life sample to start of trace. This means that output trace length will be (endSample-startSample+1) samples long" );

  pdef->addParam( "nsamples", "Number of samples per output trace", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Number of output samples", "Specify only in case of 'mode shot_records'" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_select_time_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_select_time::VariableStruct* vars = reinterpret_cast<mod_select_time::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_select_time_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_select_time::VariableStruct* vars = reinterpret_cast<mod_select_time::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->startTimes != NULL ) {
    for( int iloc = 0; iloc < vars->numLocations; iloc++ ) {
      delete [] vars->startTimes[iloc];
      delete [] vars->endTimes[iloc];
    }
    delete [] vars->startTimes;
    delete [] vars->endTimes;
    vars->startTimes = NULL;
    vars->endTimes = NULL;
  }
  if( vars->numLines != NULL ) {
    delete [] vars->numLines;
    vars->numLines = NULL;
  }
  if( vars->shotAttrList != NULL ) {
    delete vars->shotAttrList;
    vars->shotAttrList = NULL;
  }
  if( vars->shotBuffer != NULL ) {
    delete [] vars->shotBuffer;
    vars->shotBuffer = NULL;
  }
  if( vars->shotBufferTemp != NULL ) {
    delete [] vars->shotBufferTemp;
    vars->shotBufferTemp = NULL;
  }
  if( vars->keyValues != NULL ) {
    delete [] vars->keyValues;
    vars->keyValues  = NULL;
  }
  if( vars->interpol != NULL ) {
    delete vars->interpol;
    vars->interpol = NULL;
  }
  if( vars->traceGatherKeep != NULL ) {
    delete vars->traceGatherKeep;
    vars->traceGatherKeep = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_select_time_( csParamDef* pdef ) {
  params_mod_select_time_( pdef );
}
extern "C" void _init_mod_select_time_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_select_time_( param, env, writer );
}
extern "C" bool _start_exec_mod_select_time_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_select_time_( env, writer );
}
extern "C" void _exec_mod_select_time_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_select_time_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_select_time_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_select_time_( env, writer );
}

void addShot( VariableStruct* vars, cseis_system::csTraceGather* traceGather, csInt64_t& startTimeCurrentTrace_us, csInt64_t& endTimeCurrentTrace_us, cseis_system::csSuperHeader const* shdr ) {
  vars->traceGatherKeep->copyTraceTo( 0, traceGather );  // Use earliest trace as template for output trace
  cseis_system::csTrace* trace  = traceGather->trace( traceGather->numTraces()-1 );
  cseis_system::csTraceHeader* trcHdr = trace->getTraceHeader();
  float* samples = trace->getTraceSamples();

  int time_samp1_s  = (int)(vars->startTimeCurrentShot_us/1000000LL);
  int time_samp1_us = (int)(vars->startTimeCurrentShot_us%1000000LL);
  trcHdr->setIntValue( vars->hdrID_time_samp1_s, time_samp1_s );
  trcHdr->setIntValue( vars->hdrID_time_samp1_us, time_samp1_us );
  trcHdr->setIntValue( vars->hdrID_key, vars->idCurrentShot );
  startTimeCurrentTrace_us = (csInt64_t)trcHdr->intValue( vars->hdrID_time_samp1_s )*1000000LL + (csInt64_t)trcHdr->intValue( vars->hdrID_time_samp1_us );
  endTimeCurrentTrace_us   = startTimeCurrentTrace_us + (csInt64_t)( (vars->numSamplesOrig-1) * shdr->sampleInt )*1000LL;
  memcpy( samples, vars->shotBuffer, shdr->numSamples*sizeof(float) );
  for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
    vars->shotBuffer[isamp] = 0;
  }
  vars->isNonZeroShotInMemory = false;
}

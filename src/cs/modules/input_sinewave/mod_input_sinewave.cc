/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: INPUT_SINEWAVE
 * The base code of this module was converted from a demo module of the ProMAX system
 *
 * @date   2007
 */
namespace mod_input_sinewave {
  struct VariableStruct {
    int   nens;            // number of trace ensembles to generate
    int   nTraces;         // Number of traces in ensemble to generate
    int   nfreqs;          // number of frequencies in output traces
    int   nspikes;         // Number of spikes to generate
    int*  spikeSamples;    // Sample index where spikes shall be put onto trace
    float *freqs;          // list of frequencies in output traces 
    float *phase;          // list of phase in output traces 
    int   currEns;         // the current ensemble number 
    int   currTrcInEns;    // current trace number in ensemble 
    bool  firstVisit;      // flag for first visit to exec
    bool  isShot; // true if shot records are to be generated, otherwise cdps
  };
}
using mod_input_sinewave::VariableStruct;

static void sineWaveAdd( float freq, float phase, float sampInt, int numsmp, float *trace );

#define FALSE false;
#define TRUE true;

///************************************************************************************************
// Init phase
//
//
//
///************************************************************************************************
void init_mod_input_sinewave_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  csVector<std::string> valueList;

  vars->freqs = NULL;
  vars->phase = NULL;
  vars->spikeSamples = NULL;

  
// get the primary key name 
  std::string pkeyName;
  param->getString( "pkeynam", &pkeyName );
  if( pkeyName.compare( "source" ) != 0 && pkeyName.compare( "cdp" ) != 0  ){
    writer->error("Only SOURCE and CDP can be primary keys. Specified key is: %s", pkeyName.c_str());
  }

  valueList.clear();
  // get and decode the list of frequencies
  if( param->exists("freq") ) {
    param->getAll( "freq", &valueList );
  }
  //  if( valueList.size() < 1 ){
  //    writer->error("No frequencies were input!");
  // }

// allocate enough space for the frequency list at max length 
  vars->nfreqs = valueList.size();
  int nfreqs = vars->nfreqs;
  if( vars->nfreqs > 0 ) {
    vars->freqs = new float[nfreqs];
    vars->phase = new float[nfreqs];
    for( int i = 0; i < nfreqs; i++ ) {
      vars->freqs[i] = atof( valueList.at(i).c_str() );
      vars->phase[i] = 0.0;
    }
  }

  valueList.clear();
  if( param->exists("phase" ) ) {
    param->getAll( "phase", &valueList );
  }
  int nPhase = valueList.size();
  if( nPhase != 0 && nPhase != nfreqs ) {
    writer->error("Unequal number of frequency/phase pairs");
  }
  if( nPhase > 0 ) {
    for( int i = 0; i < nfreqs; i++ ) {
      vars->phase[i] = 0.0;
    }
    for( int i = 0; i < nPhase; i++ ) {
      vars->phase[i] = DEG2RAD( atof( valueList.at(i).c_str() ) );
    }
  }

  // get the number of traces per ensemble, the sample interval the
  // length of the traces, and the number of ensembles
  int ntrPerEns = 0;
  int nens      = 0;
  float traceLen  = 0.0;
  float sampInt  = 0.0;

  param->getInt( "ntraces", &ntrPerEns );
  param->getInt( "nens", &nens );
  param->getFloat( "tracelen", &traceLen );
  param->getFloat( "sample_int", &sampInt );
  
  // check and set the number of traces per ensemble
  if( ntrPerEns <= 0 ){
    writer->error("Number of traces per ensemble is < 1.");
  }
  vars->nTraces = ntrPerEns;
  // assign the sample interval
  if( sampInt <= 0.0 ){
    writer->error("The sample interval cannot be 0.0 or less.");
  }

  // assign the trace length, use the agfc.h NINT macro
  if( traceLen <= 0.0 ){
    writer->error("The trace length cannot be 0.0 or less.");
  }

  // check the number of ensembles
  if( nens < 1 ){
    writer->error("Number of ensembles is less than 1.");
  }

  valueList.clear();
  if( param->exists( "spikes" ) ) {
    param->getAll( "spikes", &valueList );
  }
  vars->nspikes = valueList.size();
  vars->spikeSamples = NULL;
  if( vars->nspikes > 0 ) {
    vars->spikeSamples = new int[vars->nspikes];
    for( int i = 0; i < vars->nspikes; i++ ) {
      float time  = atof( valueList.at(i).c_str() );
      if( time < 0.0 || time > traceLen ) {
        writer->error("Spike time outside of specified trace length: %f", time);
      }
      vars->spikeSamples[i] = (int)( time/sampInt + 0.5 );
    }
  }


// define the minimum standard header entries 
//  defStdHdr();

// set the primary sort flag depending on primary key name 
  if( strncmp( pkeyName.c_str(), "cdp", 3 ) == 0 ){
    hdef->addStandardHeader( "cmp" );
    hdef->addStandardHeader( "offset" );
    hdef->addStandardHeader( "aoffset" );
    vars->isShot = false;
  }
  else{
    hdef->addStandardHeader( "source" );
    hdef->addStandardHeader( "ffid" );
    hdef->addStandardHeader( "chan" );
    vars->isShot = true;
  }

  shdr->sampleInt = sampInt;
  shdr->numSamples = (int)( traceLen/shdr->sampleInt + 0.5 ) + 1;
  shdr->domain = DOMAIN_XT;

// set the external parameters 
  vars->nens         = nens;
  vars->currEns      = 1;
  vars->currTrcInEns = 1;
  vars->firstVisit   = true;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_input_sinewave_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef*         edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;

  csTrace* trace = traceGather->trace(0);


  csTraceHeader* trcHdr = trace->getTraceHeader();
  float* samples = trace->getTraceSamples();

  float   *freqs       = vars->freqs        ;
  float   *phase       = vars->phase        ;

  if( edef->isDebug() ) writer->line("Traces: %d %d  %d %d", vars->currEns, vars->nens, vars->currTrcInEns, vars->nTraces );

// fill the current trace and header buffers 
//  noDataFound = sineWaveNextTr( samples, ithdr, rthdr, freqs, phase );
  if( vars->currEns > vars->nens ){
  // .. we have no more data to output, this is equivalent to
  // .. failing to read a trace from tape or disk
    if( edef->isDebug() ) writer->line("No more trace to generate...");
    traceGather->freeAllTraces();
    return;
  }
  for( int i = 0; i < shdr->numSamples; i++ ) {
    samples[i] = 0.0;
  }

// add the sine_wave data to the trace
  for( int i = 0; i < vars->nfreqs; i++ ){
    sineWaveAdd( freqs[i], phase[i], shdr->sampleInt, shdr->numSamples, samples );
  }
  
  for( int i = 0; i < vars->nspikes; i++ ){
    samples[vars->spikeSamples[i]] = 1.0;
  }

  if( vars->isShot ){
// .. inputting shots
    trcHdr->setIntValue( hdef->headerIndex("source"), vars->currEns );
    trcHdr->setIntValue( hdef->headerIndex("ffid"), vars->currEns );
    trcHdr->setIntValue( hdef->headerIndex("chan"), vars->currTrcInEns );
  }
  else{
// .. inputting CDPs
    trcHdr->setIntValue( hdef->headerIndex("cdp"), vars->currEns );
    trcHdr->setFloatValue( hdef->headerIndex("offset"), (float)vars->currTrcInEns );
  }
  vars->currTrcInEns += 1;
  if( vars->currTrcInEns > vars->nTraces ){
    vars->currEns += 1;
  }
  return;
}

//-------------------------------------------------------------------
// sineWaveAdd
//
// subroutine to add a sine wave to a trace
// input:
//   freq - frequency of sine wave to add
//   sampInt - float sample interval of trace in ms
//   numsmp  - number of samples in trace
// output:
//   trace - the trace
//
//-------------------------------------------------------------------

static void sineWaveAdd( float freq, float phase, float sampInt, int numsmp, float *trace )
{
  double angle, twoPi, time;

  twoPi = asin(1.0)*4.0;

  for( int i = 0; i < numsmp; i++ ){
    time = (float)i * sampInt/1000.;
    angle = freq*time*twoPi;
    trace[i] += (float)sin( angle + phase );
  }
}

//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_sinewave_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_SINEWAVE", "Create traces with sine waves", "This module was converted from a publically available ProMAX demo module" );

  pdef->addParam( "pkeynam", "Primary key name", NUM_VALUES_FIXED );
  pdef->addValue( "source", VALTYPE_OPTION );
  pdef->addOption( "source", "Create shot gathers" );
  pdef->addOption( "cdp", "Create CDP gathers" );

  pdef->addParam( "nens", "Number of trace ensembles to generate", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of trace ensembles to generate" );

  pdef->addParam( "ntraces", "Number of traces per ensemble", NUM_VALUES_FIXED );
  pdef->addValue( "100", VALTYPE_NUMBER, "Number of traces per ensemble" );

  pdef->addParam( "tracelen", "Trace length for the sinewaveetic traces [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "1000", VALTYPE_NUMBER, "Trace length for the sinewaveetic traces [ms]" );

  pdef->addParam( "sample_int", "Sample interval [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "2", VALTYPE_NUMBER, "Sample interval [ms]" );

  pdef->addParam( "spikes", "List of spike times", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "List of spike times (spikes of amplitude 1.0)" );

  pdef->addParam( "freq", "List of frequencies [Hz]", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "List of frequencies [Hz]" );

  pdef->addParam( "phase", "List of phases [deg]", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "List of phases [deg]" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_input_sinewave_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_input_sinewave::VariableStruct* vars = reinterpret_cast<mod_input_sinewave::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_input_sinewave_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_input_sinewave::VariableStruct* vars = reinterpret_cast<mod_input_sinewave::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->freqs ) {
    delete [] vars->freqs;
    vars->freqs = NULL;
  }
  if( vars->phase ) {
    delete [] vars->phase;
    vars->phase = NULL;
  }
  if( vars->spikeSamples ) {
    delete [] vars->spikeSamples;
    vars->spikeSamples = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_input_sinewave_( csParamDef* pdef ) {
  params_mod_input_sinewave_( pdef );
}
extern "C" void _init_mod_input_sinewave_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_sinewave_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_sinewave_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_sinewave_( env, writer );
}
extern "C" void _exec_mod_input_sinewave_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_sinewave_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_sinewave_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_sinewave_( env, writer );
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csStandardHeaders.h"
#include "csTaper.h"
#include <cmath>
#include <cstring>
#include <ctime>    // For time
#include <cstdlib>  // For srand() and rand()

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: input_create
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_input_create {
  struct Receiver {
    int p;
    double x;
    double y;
  };
  struct Sweep {
    int len_samp;
    int end_samp;
    int start_samp;
    int taper_len_samp;
    double f1;
    double f2;
  };
  struct VariableStruct {
    float value;
    bool isValue;
    int numTraces;
    bool atEOF;

    int numSpikes;
    int* samplesSpikes;
    float* values;

    float maxNoise;

    int hdrId_trcno;
    int hdrId_rcv;
    int hdrId_rec_x;
    int hdrId_rec_y;

    int traceCounter;

    // Plane wave:
    float pw_azimuth;  // in [rad]
    float pw_slowness; // s/m
    float pw_reftime_ms;
    float pw_x0;
    float pw_y0;
    int numPlaneWaves;

    bool isRicker;
    float rickerFreq;
    float rickerAmp;
    float rickerPhase_rad;
    float rickerDamping;

    mod_input_create::Sweep* sweep;

    // Point source:
    float* ps_x0;
    float* ps_y0;
    float ps_slowness; // s/m
    float ps_reftime_ms;
    float ps_spreadingFactor;
    int numPointSources;

    Receiver* receiver;
  };

  double computeSweepSample( double freq, double dt_s, double& phase );
  bool addRickerWavelet( float freq_hz, float dampingFactor, float phaseShift_rad, float sampleInt_ms,
                         float time_ms, float amplitude, float* samples, int nSamples );
  void createRickerWavelet( float freq, float phaseShift_rad, float dampingFactor, float sampleInt_ms, int nSamples, float* samples );
}
using mod_input_create::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_create_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new mod_input_create::VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->hdrId_trcno  = -1;
  vars->hdrId_rcv    = -1;
  vars->hdrId_rec_x  = -1;
  vars->hdrId_rec_y  = -1;
  vars->atEOF        = false;
  vars->traceCounter = 0;
  vars->numTraces    = 1;
  vars->isValue      = false;
  vars->value        = 1.0f;
  vars->numSpikes    = 0;
  vars->samplesSpikes = NULL;
  vars->values       = NULL;
  vars->maxNoise     = 0.0;

  vars->pw_azimuth  = 0.0;
  vars->pw_slowness = 0.0;
  vars->pw_reftime_ms = 0.0;
  vars->pw_x0       = 0.0;
  vars->pw_y0       = 0.0;
  vars->numPlaneWaves = 0;

  vars->ps_slowness = 0.0;
  vars->ps_reftime_ms = 0.0;
  vars->ps_x0       = NULL;
  vars->ps_y0       = NULL;
  vars->ps_spreadingFactor = 0;
  vars->numPointSources = 0;

  vars->receiver = NULL;
  vars->sweep    = NULL;

  vars->isRicker    = false;
  vars->rickerFreq  = 15;
  vars->rickerAmp   = 1;
  vars->rickerPhase_rad = 0;
  vars->rickerDamping   = 1.0;

  //----------------------------------------------------
  double length;
  float sampleInt;

  param->getFloat( "sample_int", &sampleInt, 0 );
  shdr->sampleInt = sampleInt;

  string unit("ms");
  param->getDouble( "length", &length, 0 );
  if( param->getNumValues("length") > 1 ) {
    param->getString( "length", &unit, 1 );
  }
  if( param->exists("value") ) {
    param->getFloat( "value", &vars->value );
    vars->isValue = true;
  }

  if( param->exists("ricker") ) {
    vars->isRicker = true;
    int numValues = param->getNumValues("ricker");
    param->getFloat( "ricker", &vars->rickerFreq, 0 );
    int counter = 1;
    if( numValues > counter ) {
      param->getFloat( "ricker", &vars->rickerAmp, counter++ );
      if( numValues > counter ) {
        param->getFloat( "ricker", &vars->rickerPhase_rad, counter++ );
        vars->rickerPhase_rad *= M_PI / 180.0;
        if( numValues > counter ) {
          param->getFloat( "ricker", &vars->rickerDamping, counter++ );
          if( vars->rickerDamping < 1e-04 ) vars->rickerDamping = 1e-04;
        }
      }
    }
  }

  if( !unit.compare("samples") ) {
    shdr->domain     = DOMAIN_XT;
    shdr->numSamples = (int)(length + 0.5);
  }
  else if( !unit.compare("s") || !unit.compare("seconds") ) {
    shdr->domain     = DOMAIN_XT;
    shdr->numSamples = (int)( 1000*(length/shdr->sampleInt) + 0.5 );
  }
  else if( !unit.compare("ms") || !unit.compare("milliseconds") ) {
    shdr->domain     = DOMAIN_XT;
    shdr->numSamples = (int)( (length/shdr->sampleInt) + 0.5 );
  }
  else if( !unit.compare("m") || !unit.compare("meters") ) {
    shdr->domain     = DOMAIN_XD;
    shdr->numSamples = (int)( (length/shdr->sampleInt) + 0.5 );
  }
  else if( !unit.compare("km") || !unit.compare("kilometers") ) {
    shdr->domain     = DOMAIN_XD;
    shdr->numSamples = (int)( 1000*(length/shdr->sampleInt) + 0.5 );
  }
  else if( !unit.compare("Hz") ) {
    shdr->domain     = DOMAIN_FX;
    shdr->numSamples = (int)( (length/shdr->sampleInt) + 0.5 );
  }
  else {
    writer->error("Unknown option: %s", unit.c_str() );
  }

  //-------------------------------------------------------------------------
  if( param->exists("noise") ) {
    param->getFloat("noise", &vars->maxNoise);
  }
  //-------------------------------------------------------------------------
  if( param->exists( "plane_wave" ) ) {
    int iv = 0;
    vars->numPlaneWaves = 1;
    param->getFloat( "plane_wave", &vars->pw_azimuth, iv++ );
    param->getFloat( "plane_wave", &vars->pw_slowness, iv++ );
    vars->pw_slowness /= 1000.0f;   // Convert from [s/km] to [s/m]
    param->getFloat( "plane_wave", &vars->pw_reftime_ms, iv++ );
    param->getFloat( "plane_wave", &vars->pw_x0, iv++ );
    param->getFloat( "plane_wave", &vars->pw_y0, iv++ );
    vars->pw_azimuth *= (float)(M_PI / 180.0);
    if( vars->value == 0.0 ) {
      writer->error("When computing plane waves, a value (parameter 'value') other than 0 must be specified. Plane wave arrivals will be scaled by this value.");
    }
  }
  //-------------------------------------------------------------------------
  if( param->exists( "point_source" ) ) {
    vars->numPointSources = param->getNumLines("point_source_xy");
    if( vars->numPointSources == 0 ) {
      writer->error("Specify at least one point source location with parameter 'point_source_xy'");
    }
    vars->ps_x0 = new float[vars->numPointSources];
    vars->ps_y0 = new float[vars->numPointSources];
    for( int i = 0; i < vars->numPointSources; i++ ) {
      param->getFloatAtLine( "point_source_xy", &vars->ps_x0[i], i, 0 );
      param->getFloatAtLine( "point_source_xy", &vars->ps_y0[i], i, 1 );
    }
    param->getFloat( "point_source", &vars->ps_slowness, 0 );
    vars->ps_slowness /= 1000.0f;   // Convert from [s/km] to [s/m]
    param->getFloat( "point_source", &vars->ps_reftime_ms, 1 );
    param->getFloat( "point_source", &vars->ps_spreadingFactor, 2 );
    if( vars->value == 0.0 ) {
      writer->error("When computing point sources, a value (parameter 'value') other than 0 must be specified. Arrivals will be scaled by this value.");
    }
  }
  //-------------------------------------------------------------------------
  if( param->exists("spikes") ) {
    vars->numSpikes = param->getNumValues("spikes");
    vars->samplesSpikes = new int[vars->numSpikes];
    for( int i = 0; i < vars->numSpikes; i++ ) {
      float time;
      param->getFloat("spikes", &time, i);
      vars->samplesSpikes[i] = (int)( time / shdr->sampleInt + 0.5);
    }

    vars->values = new float[vars->numSpikes];
    if( param->exists("values") ) {
      if( vars->isValue ) writer->error("Parameter 'values' cannot be specified in conjunction with parameter 'value'. Specify only one of these two");
      int numValues = param->getNumValues("values");
      if( numValues != vars->numSpikes ) writer->error("Parameter 'values': Inconsistent number of values (=%d). Should match number of spikes (=%d)", numValues, vars->numSpikes);
      for( int i = 0; i < vars->numSpikes; i++ ) {
        param->getFloat("values", &vars->values[i], i);
      }
    }
    else {
      if( vars->value == 0.0f ) vars->value = 1.0f;
      for( int i = 0; i < vars->numSpikes; i++ ) {
        vars->values[i] = vars->value;
      }
    }

  }

  //-------------------------------------------------------------------------
  if( param->exists("sweep_lin") ) {
    vars->sweep = new mod_input_create::Sweep();
    float start_time;
    float len_time;
    float taper_len_time;
    param->getFloat("sweep_lin", &start_time, 0);
    param->getFloat("sweep_lin", &len_time, 1);
    param->getFloat("sweep_lin", &taper_len_time, 2);
    param->getDouble("sweep_lin", &vars->sweep->f1, 3);
    param->getDouble("sweep_lin", &vars->sweep->f2, 4);
    vars->sweep->start_samp     = (int)round( start_time / shdr->sampleInt ) + 1;
    vars->sweep->len_samp       = (int)round( len_time / shdr->sampleInt );
    vars->sweep->end_samp       = vars->sweep->start_samp + vars->sweep->len_samp - 1;
    vars->sweep->taper_len_samp = (int)round( taper_len_time / shdr->sampleInt );
    if( vars->sweep->end_samp > shdr->numSamples ) writer->error("Sweep length (%.4f + %.4f) exceeds trace length (%.4f)", len_time, start_time, shdr->numSamples*shdr->sampleInt);
    if( 2*vars->sweep->taper_len_samp > vars->sweep->len_samp ) writer->error("Sweep taper length x2 (%.4f) exceeds sweep length (%.4f)", taper_len_time*2, len_time);
    writer->line("Sweep paramters:");
    writer->line("  Start sample:  %d", vars->sweep->start_samp);
    writer->line("  End sample:    %d", vars->sweep->end_samp);
    writer->line("  Length:        %d", vars->sweep->len_samp);
    writer->line("  Taper length:  %d", vars->sweep->taper_len_samp);
  }

  //-------------------------------------------------------------------------
  if( param->exists("ntraces") ) {
    param->getInt( "ntraces", &vars->numTraces );
  }

  //-------------------------------------------------------------------------
  //
  if( param->exists("rec_geom") || vars->numPlaneWaves > 0 ) {
    std::string text;
    param->getString( "rec_geom", &text );
    FILE* fin = fopen(text.c_str(),"r");
    if( fin == NULL ) {
      writer->error("Problems opening receiver geometry file: %s", text.c_str() );
    }
    csVector<mod_input_create::Receiver> receiverList;
    char line[132];
    while( fgets( line, 132, fin ) != NULL ) {
      mod_input_create::Receiver rcv;
      sscanf(line,"%d %lf %lf", &rcv.p, &rcv.x, &rcv.y);
      receiverList.insert(rcv);
    }
    vars->receiver = new mod_input_create::Receiver[receiverList.size()];
    for( int i = 0; i < receiverList.size(); i++ ) {
      vars->receiver[i].p = receiverList.at(i).p;
      vars->receiver[i].x = receiverList.at(i).x;
      vars->receiver[i].y = receiverList.at(i).y;
    }
    if( param->exists("ntraces") && vars->numTraces != receiverList.size() ) {
      writer->warning("Receiver geometry overrules specified number of traces (=%d): Set number of traces to %d.", vars->numTraces, receiverList.size());
    }
    vars->numTraces = receiverList.size();
    vars->hdrId_rcv   = hdef->addHeader( &cseis_geolib::HDR_RCV );
    vars->hdrId_rec_x = hdef->addHeader( &cseis_geolib::HDR_REC_X );
    vars->hdrId_rec_y = hdef->addHeader( &cseis_geolib::HDR_REC_Y );

    if( edef->isDebug() ) {
      for( int i = 0; i < receiverList.size(); i++ ) {
        writer->line("#%2d   %d %lf %lf", i+1, vars->receiver[i].p, vars->receiver[i].x, vars->receiver[i].y);
      }
    }
  }
  //-------------------------------------------------------------------------
  //
  vars->hdrId_trcno = hdef->addHeader( csStandardHeaders::get("trcno") );
  vars->traceCounter = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_input_create_(
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

  if( vars->traceCounter == vars->numTraces-1 ) {
    vars->atEOF = true;
  }

  float* samples = trace->getTraceSamples();
  for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
    samples[isamp] = 0.0;
  }

  //-------------------------------------------------------------------------------------------
  // Add spikes
  //
  if( vars->samplesSpikes != NULL ) {
    for( int ispike = 0; ispike < vars->numSpikes; ispike++ ) {
      if( vars->samplesSpikes[ispike] < shdr->numSamples && vars->samplesSpikes[ispike] >= 0 ) {
        samples[ vars->samplesSpikes[ispike] ] = vars->values[ispike];
      }
    }
  }
  else if( vars->numPlaneWaves > 0 || vars->numPointSources > 0 ) {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] = 0.0;
    }
    //-------------------------------------------------------------------------------------------
    // Add plane wave
    //
    if( vars->numPlaneWaves > 0 ) {
      double dx = vars->receiver[vars->traceCounter].x - vars->pw_x0;
      double dy = vars->receiver[vars->traceCounter].y - vars->pw_y0;
      double projection = dx * sin(vars->pw_azimuth) + dy * cos(vars->pw_azimuth);
      float delayTime_ms = (float)( projection * vars->pw_slowness ) * 1000.0; // [ms]
      int nearestSampleIndex = (int)round( (delayTime_ms+vars->pw_reftime_ms) / shdr->sampleInt );
      if( nearestSampleIndex >= 0 && nearestSampleIndex < shdr->numSamples ) {
        samples[nearestSampleIndex] += vars->value;
      }
      if( edef->isDebug() ) {
        writer->line("Receiver %3d  %12f %12f  delay time  %10f ms (+%fms reference time)  sample_index %d",
                     vars->receiver[vars->traceCounter].p, vars->receiver[vars->traceCounter].x, vars->receiver[vars->traceCounter].y,
                     delayTime_ms, vars->pw_reftime_ms, nearestSampleIndex);
      }
    }
    //-------------------------------------------------------------------------------------------
    // Add point source
    //
    for( int ip = 0; ip < vars->numPointSources; ip++ ) {
      double dx = vars->receiver[vars->traceCounter].x - vars->ps_x0[ip];
      double dy = vars->receiver[vars->traceCounter].y - vars->ps_y0[ip];
      double offset = sqrt( dx * dx + dy * dy );
      float delayTime_ms = (float)( offset * vars->ps_slowness * 1000.0 ); // [ms]
      int nearestSampleIndex = (int)round( (delayTime_ms+vars->ps_reftime_ms) / shdr->sampleInt );
      if( nearestSampleIndex >= 0 && nearestSampleIndex < shdr->numSamples ) {
        samples[nearestSampleIndex] += vars->value / (float)max(1.0,pow(offset,(double)vars->ps_spreadingFactor));
      }
      if( edef->isDebug() ) {
        writer->line("Receiver %3d  %12f %12f  delay time  %10f ms (+%fms reference time)  sample_index %d, offset: %f",
                     vars->receiver[vars->traceCounter].p, vars->receiver[vars->traceCounter].x, vars->receiver[vars->traceCounter].y,
                     delayTime_ms, vars->ps_reftime_ms, nearestSampleIndex, offset);
      }
    }
  }
  //-------------------------------------------------------------------------------------------
  // Create sweep
  else if( vars->sweep != NULL ) {
    for( int isamp = 0; isamp < vars->sweep->start_samp; isamp++ ) {
      samples[isamp] = 0;
    }
    for( int isamp = vars->sweep->end_samp+1; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] = 0;
    }
    double df = (vars->sweep->f2 - vars->sweep->f1) / (double)( vars->sweep->len_samp-1 );
    int startSamp = vars->sweep->start_samp;
    int endSamp   = vars->sweep->end_samp;
    double phasePrev = 0;
    for( int isamp = startSamp; isamp <= endSamp; isamp++ ) {
      int sweepIndex = isamp - startSamp;
      double freq = vars->sweep->f1 + df * (double)sweepIndex;
      samples[isamp] = vars->value * mod_input_create::computeSweepSample( freq, (double)shdr->sampleInt/1000.0, phasePrev );
      //      fprintf(stdout,"%d %f %f  %f\n", sweepIndex, (double)isamp*(double)shdr->sampleInt, freq, phasePrev*180/M_PI);
    }
    float normScalarDummy;
    cseis_geolib::csTaper taper( vars->sweep->len_samp, vars->sweep->taper_len_samp, cseis_geolib::csTaper::COSINE, cseis_geolib::csTaper::NORM_NO );
    taper.apply( &samples[vars->sweep->start_samp], normScalarDummy );
  }
  //-------------------------------------------------------------------------------------------
  // Set to constant value
  //
  else if( vars->isValue ) {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] = vars->value;
    }
  }
  //-------------------------------------------------------------------------------------------
  // Add noise
  if( vars->maxNoise != 0.0 ) {
    srand( (vars->traceCounter+1)*shdr->numSamples );
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      double random  = 2.0 * ( (double)rand() - 0.5*(double)RAND_MAX ) / (double)RAND_MAX;
      samples[isamp] += (float)random * vars->maxNoise;
    }
  }
  csTraceHeader* trcHdr = trace->getTraceHeader();
  trcHdr->setIntValue( vars->hdrId_trcno, vars->traceCounter+1 );

  if( vars->hdrId_rcv >= 0 ) {
    trcHdr->setIntValue( vars->hdrId_rcv, vars->receiver[vars->traceCounter].p );
    trcHdr->setDoubleValue( vars->hdrId_rec_x, vars->receiver[vars->traceCounter].x );
    trcHdr->setDoubleValue( vars->hdrId_rec_y, vars->receiver[vars->traceCounter].y );
  }

  if( vars->isRicker ) {
    mod_input_create::createRickerWavelet( vars->rickerFreq, vars->rickerPhase_rad, vars->rickerDamping, shdr->sampleInt, shdr->numSamples, samples );
  }

  vars->traceCounter += 1;

  return;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_create_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_CREATE", "Create traces", "Simple synthetic trace generator" );

  pdef->addParam( "length", "Trace length", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Length of traces" );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );
  pdef->addOption( "samples", "Samples" );
  pdef->addOption( "m", "Meters" );
  pdef->addOption( "km", "Kilometers" );

  pdef->addParam( "ntraces", "Number of traces", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of traces" );

  pdef->addParam( "sample_int", "Sample interval", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Sample interval, unit specified under user parameter 'length'" );

  pdef->addParam( "value", "Constant value to set trace samples", NUM_VALUES_FIXED );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Set trace samples to this value" );

  pdef->addParam( "noise", "Add Gaussian random noise", NUM_VALUES_FIXED );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Maximum amplitude of random noise" );

  pdef->addParam( "spikes", "Add spikes", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "List of times/depths/frequencies (unit see user parameter 'length') where a spike shall be added (rest of trace is set to 0)" );

  pdef->addParam( "values", "List of values corresponding to spikes", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Set each spike to the corresponding value (use parameter 'value' to set each spike to the same value)" );

  pdef->addParam( "ricker", "Generate Ricker wavelet centered on output trace", NUM_VALUES_VARIABLE );
  pdef->addValue( "15.0", VALTYPE_NUMBER, "Frequency [Hz]" );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Amplitude" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Phase shift [deg]" );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Damping factor (1: No damping, >1: Damping, <1: Boosting)" );

  pdef->addParam( "sweep_lin", "Add linear sweep. Specify user parameter 'value'", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Start time [ms]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Length [ms] - ...includes taper at both ends" );
  pdef->addValue( "", VALTYPE_NUMBER, "Taper length [ms] - Sweep is tapered at both ends" );
  pdef->addValue( "", VALTYPE_NUMBER, "Start frequency [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "End frequency [Hz]" );

  pdef->addParam( "plane_wave", "Generate plane wave (spikes only)", NUM_VALUES_FIXED, "Requires the definiton of a receiver geometry" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Azimuth [deg]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Slowness [s/km]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Reference time [ms]: Time to add to computed delay time" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Origin X coordinate [m] for delay time computation" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Origin Y coordinate [m] for delay time computation" );

  pdef->addParam( "point_source", "Generate point source (spikes only)", NUM_VALUES_FIXED, "Requires the definiton of a receiver geometry" );
  //  pdef->addValue( "0.0", VALTYPE_NUMBER, "Point source Z coordinate [m]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Slowness [s/km]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Reference time [ms]: Time to add to computed travel time" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Geometric spreading 1/r^N" );

  pdef->addParam( "point_source_xy", "Point source location", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Point source X coordinate [m]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Point source Y coordinate [m]" );

  pdef->addParam( "rec_geom", "Receiver geometry definition", NUM_VALUES_FIXED, "Specify ASCII file name containing listing of receiver XY positions" );
  pdef->addValue( "", VALTYPE_STRING, "ASCII file name. Format:  rcv  rec_x  rec_y" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_input_create_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_input_create::VariableStruct* vars = reinterpret_cast<mod_input_create::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_input_create_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_input_create::VariableStruct* vars = reinterpret_cast<mod_input_create::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->samplesSpikes != NULL ) {
    delete [] vars->samplesSpikes; vars->samplesSpikes = NULL;
  }
  if( vars->receiver != NULL ) {
    delete [] vars->receiver;
    vars->receiver = NULL;
  }
  if( vars->sweep != NULL ) {
    delete [] vars->sweep;
    vars->sweep = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_input_create_( csParamDef* pdef ) {
  params_mod_input_create_( pdef );
}
extern "C" void _init_mod_input_create_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_create_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_create_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_create_( env, writer );
}
extern "C" void _exec_mod_input_create_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_create_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_create_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_create_( env, writer );
}

namespace mod_input_create {
  double computeSweepSample( double freq, double dt_s, double& phase ) {
    phase += 2.0 * M_PI * freq * dt_s;
    return( cos( phase ) );
  }
  bool addRickerWavelet( float freq_hz, float dampingFactor, float phaseShift_rad, float sampleInt_ms,
                         float timeIn_ms, float amplitudeIn, float* samples, int nSamples )
  {
    if( dampingFactor < 0.0001 ) {
      dampingFactor = 2.0f;
    }

    float sampleInt_s = sampleInt_ms / 1000.0f;
    float timeIn_s     = timeIn_ms / 1000.0f;
    float omegaHalf    = (float)( M_PI * freq_hz );
    float omegaHalf_sq = omegaHalf * omegaHalf;
    float exponentDamped_sq = omegaHalf_sq / (dampingFactor * dampingFactor);

    fprintf(stderr,"Damping factor: %f  %f\n", dampingFactor, exponentDamped_sq );

    float waveletWidth_s = 4.0f * (float)( dampingFactor * 3.0 * sqrt(6) / M_PI ) / freq_hz;
    float waveletWidthHalf_s = 0.5f * waveletWidth_s;
    float time1 = timeIn_s - waveletWidthHalf_s;
    float time2 = timeIn_s + waveletWidthHalf_s;
    float maxTime = (float)nSamples * sampleInt_ms * 1000.0f;
    if( time1 >= maxTime || time2 <= 0.0 ) {
      fprintf(stderr,"Outside of valid window. Max time: %.2fms, time1: %.2fms, time2: %.2f\n", maxTime * 1000.0, time1*1000, time2*1000);
      return false;
    }

    int index1 = std::max( (int)round( time1 / sampleInt_s ), 0 );
    int index2 = std::min( (int)round( time2 / sampleInt_s ), nSamples-1 );
    fprintf(stderr,"Wavelet width: %.5fs  %.5fs %.5fs\n", waveletWidth_s, time1, time2);

    for( int index = index1; index <= index2; index++ ) {
      float timeOut_s = sampleInt_s * (float)(index) - timeIn_s;
      float phase     = phaseShift_rad + 2.0f * omegaHalf * timeOut_s;
      float ampOut    = amplitudeIn * cos(phase) * exp( -exponentDamped_sq * timeOut_s * timeOut_s );
      samples[index] += ampOut;
    }

    return true;
  }

  void createRickerWavelet( float freq, float phaseShift_rad, float dampingFactor, float sampleInt_ms, int nSamples, float* samples ) {
    float sampleInt_s = sampleInt_ms / 1000.0f;
    float timeMid_s   = sampleInt_s * (float)(int)(nSamples/2);
    float omegaHalf_sqr = M_PI * M_PI * freq * freq;
    float damping_sqr = 1.0f / (dampingFactor * dampingFactor);
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      float timeRed_s = (float)isamp * sampleInt_s - timeMid_s;
      float timeRed_sqr = timeRed_s * timeRed_s;
      float temp = omegaHalf_sqr * timeRed_sqr;
      float amp = ( 1.0f - 2.0f * temp ) * exp( -temp * damping_sqr ) * cos( phaseShift_rad );
      samples[isamp] = amp;
    }

  }
 
} // END namespace


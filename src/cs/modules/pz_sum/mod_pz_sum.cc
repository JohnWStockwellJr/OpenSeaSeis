/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"
#include "csFlexNumber.h"
#include "geolib_math.h"
#include "geolib_methods.h"
#include "csFFTTools_ORIG.h"
#include "csGeolibUtils.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: PZSUM
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_pz_sum {
  struct VariableStruct {
    int input;
    int hdrId_sensor;
    int hdrId_rec_z;
    int method;
    int windowLengthSamples;
    float percent;
    float water_vel;
    int first_sample;
    int last_sample;
    bool doEqualize;
    float zScalar;

    cseis_geolib::csFFTTools_ORIG* fftTool;
    int numFFTSamples;
    float* buffer;
    float* weightsHyd;
    float* weightsGeo;

    bool isSensorUserDefined;
    int sensorNumP;
    int sensorNumZ;
    int outputOption;
  };

  static int const OUTPUT_PZ = 51;
  static int const OUTPUT_P = 52;
  static int const OUTPUT_Z = 53;
  
  static int const METHOD_SUMMATION = 1;
  static int const METHOD_ZERO = 2;
  static int const METHOD_ANALYSIS = 3;
  static int const METHOD_WEIGHT = 4;

  static int const SENSOR_INDEX_P = 1;
  static int const SENSOR_INDEX_Z = 5;
  static int const SENSOR_INDEX_Z2 = 2;
  static int const NO_VALUE = -999;
}
using namespace mod_pz_sum;

void pz_summation( float* samples_p, float const* samples_z, float percent, int nSamples, int windowLengthSamples, float zScalar, bool doEqualize );
void pz_zero( float* samples_p, float const* samples_z, int nSamples );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_pz_sum_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 2 );

  std::string text;

  vars->water_vel = 1500.0;

  vars->input  = 0;
  vars->hdrId_sensor  = -1;
  vars->hdrId_rec_z   = -1;
  vars->method = 0;
  vars->windowLengthSamples = 0;
  vars->percent      = 0;
  vars->first_sample = 0;
  vars->last_sample  = 0;
  vars->doEqualize   = false;
  vars->zScalar      = 1.0;

  vars->sensorNumP = 1;
  vars->sensorNumZ = 5;
  vars->fftTool = NULL;
  vars->numFFTSamples = 0;
  vars->buffer = NULL;
  vars->weightsHyd = NULL;
  vars->weightsGeo = NULL;
  vars->isSensorUserDefined = false;
  vars->outputOption = mod_pz_sum::OUTPUT_PZ;

//-----------------------------------------------------------------

  param->getString( "method", &text );
  if( !text.compare("sum") ) {
    vars->method = METHOD_SUMMATION;
  }
  else if( !text.compare("subtract") ) {
    vars->zScalar = -1;
    vars->method = METHOD_SUMMATION;
  }
  else if( !text.compare("zero") ) {
    vars->method = METHOD_ZERO;
  }
  else if( !text.compare("analysis") ) {
    vars->method = METHOD_ANALYSIS;
  }
  else if( !text.compare("weight_sum") ) {
    vars->method = METHOD_WEIGHT;
  }
  else {
    writer->error("Unknown option '%s'", text.c_str());
  }

  //-------------------------------------------------------------------------

  std::string sensorName("sensor");
  if( param->exists( "sensor_hdr" ) ) {
    param->getString( "sensor_hdr", &sensorName, 0 );
    param->getInt( "sensor_hdr", &vars->sensorNumP, 1 );
    param->getInt( "sensor_hdr", &vars->sensorNumZ, 2 );
    vars->isSensorUserDefined = true;
  }

  if( vars->method == METHOD_ANALYSIS ) {
    float startTime = -1.0;
    float endTime   = -2.0;
    param->getFloat( "start_time", &startTime );
    param->getFloat( "end_time", &endTime );
    if( startTime == 0.0 && endTime == 0.0 ) {
      vars->first_sample = 0;
      vars->last_sample  = shdr->numSamples-1;
    }
    else if( startTime > endTime ) {
      writer->error("Specified start/end times are invalid. Please check input parameters.");
    }
    else {
      vars->first_sample = (int)(startTime / shdr->sampleInt + 0.5);
      vars->last_sample  = (int)(endTime / shdr->sampleInt + 0.5);
    }

    if( !hdef->headerExists( "rec_z" ) ) {
      writer->error("Required trace header 'rec_z' does not exist.");
    }
    vars->hdrId_rec_z  = hdef->headerIndex( "rec_z" );
  }
  else if( vars->method == METHOD_SUMMATION ) {
    if( param->exists( "equalize" ) ) {
      param->getString( "equalize", &text );
      if( !text.compare("yes") ) {
        vars->doEqualize = true;
      }
      else if( !text.compare("no") ) {
        vars->doEqualize = false;
      }
      else {
        writer->error("Unknown option '%s'", text.c_str());
      }
    }

    //-----------------------------------------------------------------
    if( param->exists("percent") ) {
      param->getFloat( "percent", &vars->percent );
      if( vars->percent < 0.0 || vars->percent > 100.0 ) {
        writer->error("Percentage cannot be in the range [0,100]: %f", vars->percent);
      }
    }
    float windowLength = 0.0;
    if( vars->percent > 0.0 ) {
      param->getFloat( "window_length", &windowLength );
    }
    if( windowLength < 0.0 ) {
      writer->error("Window length cannot be smaller than 0: %f", windowLength);
    }
    else if( windowLength == 0.0 ) {
      vars->windowLengthSamples = shdr->numSamples;
    }
    else {
      vars->windowLengthSamples = (int)(windowLength/shdr->sampleInt + 0.5);
      if( vars->windowLengthSamples > shdr->numSamples ) {
        writer->warning("Specified window length is larger than total number of samples: %f", windowLength);
        vars->windowLengthSamples = shdr->numSamples;
      }
    }
  }
  else if( vars->method == METHOD_WEIGHT ) {
    if( param->exists("weight_output" ) ) {
      param->getString( "weight_output", &text );
      if( !text.compare("pz") ) {
        vars->outputOption = OUTPUT_PZ;
      }
      else if( !text.compare("p") ) {
        vars->outputOption = OUTPUT_P;
      }
      else if( !text.compare("z") ) {
        vars->outputOption = OUTPUT_Z;
      }
      else {
        writer->error("Unknown option '%s'", text.c_str());
      }
    }
    float rec_z;
    int numNotchesHyd;
    double  squeezeFactor;
    param->getFloat("weight_param", &rec_z, 0);
    param->getFloat("weight_param", &vars->water_vel, 1);
    param->getInt("weight_param", &numNotchesHyd, 2);
    param->getDouble("weight_param", &squeezeFactor, 3);
    float weightDamping = 1.0;
    if( param->exists("weight_damping") ) {
      param->getFloat("weight_damping", &weightDamping );
    }
    if( !hdef->headerExists( "rec_z" ) ) {
      writer->error("Required trace header 'rec_z' does not exist.");
    }
    vars->hdrId_rec_z  = hdef->headerIndex( "rec_z" );

    vars->fftTool = new cseis_geolib::csFFTTools_ORIG( shdr->numSamples );
    vars->numFFTSamples = vars->fftTool->numFFTSamples();
    vars->buffer = new float[vars->numFFTSamples + 2];
    vars->weightsHyd = new float[vars->numFFTSamples/2];
    vars->weightsGeo = new float[vars->numFFTSamples/2];
    float T_interval  = (float)vars->numFFTSamples * shdr->sampleInt/1000.0;
    float sampleInt    = 1.0/T_interval;
    float notchFreq = vars->water_vel / (2.0*rec_z);
    float freqMinConst = ( (float)numNotchesHyd - 0.5 ) * notchFreq;

    for( int ifreq = 0; ifreq <= vars->numFFTSamples/2; ifreq++ ) {
      double freq = ifreq * sampleInt;
      //      double omega = 2.0 * M_PI * freq;
      //   if( ifreq == 0 ) omega = 1.0;
      double an_rad = (freq/notchFreq) * 2.0 * M_PI;
      double number = fmod( an_rad, 2*M_PI ) - M_PI;
      if( squeezeFactor < 0.001 ) {
        an_rad = M_PI;
      }
      else if( number >= 0 ) {
        an_rad = pow( number/M_PI, 1.0/squeezeFactor )*M_PI + M_PI;
      }
      else {
        an_rad = -pow( -number/M_PI, 1.0/squeezeFactor )*M_PI + M_PI;
      }
      vars->weightsHyd[ifreq] = (float)(1.0 + -cos( an_rad ) );
      if( freq < freqMinConst ) vars->weightsHyd[ifreq] = 2.0f;
      if( weightDamping < 1.0 ) {
        vars->weightsHyd[ifreq] = weightDamping * (float)(1.0 + -cos( an_rad ) ) + (1.0-weightDamping);
        if( freq < (float)numNotchesHyd * notchFreq ) {
          float w = 0.5*( 1 + weightDamping );
          vars->weightsHyd[ifreq] = w * (float)(1.0 + -cos( an_rad ) ) + (1.0-weightDamping);
        }
        if( freq < freqMinConst ) vars->weightsHyd[ifreq] = 2.0f; // * weightDamping + (1.0-weightDamping);
      }
      if( edef->isDebug() ) fprintf(stdout,"Freq: %f  hydrophone weight: %f  geophone weight: %f\n", freq, vars->weightsHyd[ifreq], 2.0-vars->weightsHyd[ifreq]);
      vars->weightsGeo[ifreq] = 2.0f - vars->weightsHyd[ifreq];
    }
  }
  
//-------------------------------------------------------------------------
// Check headers
  if( !hdef->headerExists( sensorName ) ) {
    writer->error("Required trace header '%s' does not exist.", sensorName.c_str() );
  }
  if( hdef->headerType( sensorName ) != TYPE_INT ) {
    writer->error("Trace header '%s%' has the wrong type: %s. Should be integer type", sensorName.c_str(), csGeolibUtils::typeText( hdef->headerType( sensorName ) ) );
  }  

  vars->hdrId_sensor = hdef->headerIndex( sensorName.c_str() );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_pz_sum_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;


  int nTraces = traceGather->numTraces();
  if( nTraces < 2 ) {
    writer->error("Module %s: Incorrect number of input traces. Expected 2 (or more), found %d", nTraces, edef->moduleName().c_str() );
  }

  int sensor;
  int trace_index_p = NO_VALUE;
  int trace_index_z = NO_VALUE;

  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    sensor = traceGather->trace(itrc)->getTraceHeader()->intValue(vars->hdrId_sensor);
    if( !vars->isSensorUserDefined && sensor == SENSOR_INDEX_P ) {
      if( trace_index_p != NO_VALUE ) {
        writer->error("Found more than one sensor P (%d) trace in input trace pair.", SENSOR_INDEX_P );
      }
      trace_index_p = itrc;
    }
    else if( !vars->isSensorUserDefined && ( sensor == SENSOR_INDEX_Z || sensor == SENSOR_INDEX_Z2 ) ) {
      if( trace_index_z != NO_VALUE ) {
        writer->error("Found more than one sensor Z (%d/%d) trace in input trace pair.", SENSOR_INDEX_Z, SENSOR_INDEX_Z2 );
      }
      trace_index_z = itrc;
    }
    else if( sensor == vars->sensorNumP ) {
      if( trace_index_p != NO_VALUE ) {
        writer->error("Found more than one sensor P (%d) trace in input trace pair.", vars->sensorNumP );
      }
      trace_index_p = itrc;
    }
    else if( sensor == vars->sensorNumZ ) {
      if( trace_index_z != NO_VALUE ) {
        writer->error("Found more than one sensor Z (%d) trace in input trace pair.", vars->sensorNumZ );
      }
      trace_index_z = itrc;
    }

  }
  if( trace_index_p == NO_VALUE ) {
    if( !vars->isSensorUserDefined ) {
      writer->error("Ensemble is missing a sensor P (%d) trace.", SENSOR_INDEX_P );
    }
    else {
      writer->error("Ensemble is missing a sensor P (%d) trace.", vars->sensorNumP );
    }
  }
  if( trace_index_z == NO_VALUE ) {
    if( !vars->isSensorUserDefined ) {
      writer->error("Ensemble is missing a sensor Z (%d/%d) trace.", SENSOR_INDEX_Z, SENSOR_INDEX_Z2 );
    }
    else {
      writer->error("Ensemble is missing a sensor Z (%d) trace.", vars->sensorNumZ );
    }
  }

//---------------------------------------------------
//  csTraceHeader* trcHdr_p = traceGather->trace(trace_index_p)->getTraceHeader();
//  csTraceHeader* trcHdr_z = traceGather->trace(trace_index_z)->getTraceHeader();

  if( vars->method == METHOD_SUMMATION ) {
    pz_summation(
                 traceGather->trace(trace_index_p)->getTraceSamples(),
                 traceGather->trace(trace_index_z)->getTraceSamples(),
                 vars->percent,
                 shdr->numSamples,
                 vars->windowLengthSamples,
                 vars->zScalar,
                 vars->doEqualize );
  }
  else if( vars->method == METHOD_ZERO ) {
    pz_zero(
            traceGather->trace(trace_index_p)->getTraceSamples(),
            traceGather->trace(trace_index_z)->getTraceSamples(),
            shdr->numSamples );
  }
  else if( vars->method == METHOD_ANALYSIS ) {
    int nTraces = 1;
    float** ptr_p_data = new float*[nTraces];
    float** ptr_z_data = new float*[nTraces];
    for( int itrc = 0; itrc < nTraces; itrc++ ) {
      ptr_p_data[itrc] = &traceGather->trace(trace_index_p)->getTraceSamples()[vars->first_sample];
      ptr_z_data[itrc] = &traceGather->trace(trace_index_z)->getTraceSamples()[vars->first_sample];
    }
    float refl_coef;
    float pz_scalar;
    float wdep = traceGather->trace(0)->getTraceHeader()->floatValue(vars->hdrId_rec_z);

    float z_phase;
    /*    pz_analysis(
                ptr_p_data,
                ptr_z_data,
                nTraces,
                (vars->last_sample-vars->first_sample+1),
                &refl_coef,
                &pz_scalar,
                &wdep,
                &z_phase );  */
    pz_analysis2(
                 ptr_p_data,
                 ptr_z_data,
                 nTraces,
                 (vars->last_sample-vars->first_sample+1),
                 vars->water_vel,
                 &refl_coef,
                 &pz_scalar,
                 &wdep,
                 &z_phase );

    delete [] ptr_p_data;
    delete [] ptr_z_data;
  }
  else if( vars->method == METHOD_WEIGHT ) {
    int fftDataType = cseis_geolib::FX_AMP_PHASE;
    float* samplesHyd = traceGather->trace(trace_index_p)->getTraceSamples();
    float* samplesGeo = traceGather->trace(trace_index_z)->getTraceSamples();

    // Apply frequency weights for P data:
    bool success;
    double const* real;
    if( vars->outputOption != mod_pz_sum::OUTPUT_Z ) {
      success = vars->fftTool->fft_forward( samplesHyd, &vars->buffer[0], &vars->buffer[vars->numFFTSamples/2+1] );
      if( !success ) writer->error("FFT transform failed for unknown reasons...");
      for( int ifreq = 0; ifreq <= vars->numFFTSamples/2; ifreq++ ) {
        vars->buffer[ifreq] *= vars->weightsHyd[ifreq];
      }
      success = vars->fftTool->fft_inverse( vars->buffer, fftDataType );
      if( !success ) writer->error("FFT transform failed for unknown reasons...");
      real = vars->fftTool->realData();
      for( int i = 0; i < shdr->numSamples; i++ ) {
        samplesHyd[i] = real[i];
      }
    }
    if( vars->outputOption != mod_pz_sum::OUTPUT_P ) {
      // Apply frequency weights for Z data:
      success = vars->fftTool->fft_forward( samplesGeo, &vars->buffer[0], &vars->buffer[vars->numFFTSamples/2+1] );
      if( !success ) writer->error("FFT transform failed for unknown reasons...");
      for( int ifreq = 0; ifreq <= vars->numFFTSamples/2; ifreq++ ) {
        //      vars->buffer[ifreq] *= (2.0f - vars->weightsHyd[ifreq]);
        vars->buffer[ifreq] *= vars->weightsGeo[ifreq];
      }
      success = vars->fftTool->fft_inverse( vars->buffer, fftDataType );
      if( !success ) writer->error("FFT transform failed for unknown reasons...");
      real = vars->fftTool->realData();
      if( vars->outputOption == mod_pz_sum::OUTPUT_PZ ) {
        for( int i = 0; i < shdr->numSamples; i++ ) {
          samplesHyd[i] = 0.5 * ( samplesHyd[i] + real[i] );
        }
      }
      else if( vars->outputOption == mod_pz_sum::OUTPUT_Z ) {
        for( int i = 0; i < shdr->numSamples; i++ ) {
          samplesHyd[i] = real[i];
        }
      }
    }

  }
  traceGather->freeTrace(trace_index_z);
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_pz_sum_( csParamDef* pdef ) {
  pdef->setModule( "PZ_SUM", "PZ combination (STILL EXPERIMENTAL)", "Combine P and Z trace pair. The combined output trace will keep all trace headers from the P trace" );

  pdef->addDoc("Required input: Pairs of P and Z traces.");
  pdef->addDoc("Required trace headers in input data: sensor, rec_z");
  pdef->addDoc("Use weighted sum method only if input P & Z have been fully deghosted");

  pdef->addParam( "method", "PZ combination method", NUM_VALUES_FIXED );
  pdef->addValue( "sum", VALTYPE_OPTION );
  pdef->addOption( "sum", "Sum P and Z traces in x-t domain, sample by sample." );
  pdef->addOption( "subtract", "Subtract P and Z traces in x-t domain, sample by sample." );
  pdef->addOption( "zero", "Zero samples of opposite polarity." );
  pdef->addOption( "weight_sum", "Frequency-weighted sum. Apply weights in frequency domain, then sum in x-t domain. Use only in input data has been deghosted", "Requires user parameter 'weight_param'." );
  //  pdef->addOption( "analysis", "Analysis. Determine reflection coefficient, PZ scalar, and water depth." );

  pdef->addParam( "window_length", "Length of sliding time window [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "500.0", VALTYPE_NUMBER, "Length of sliding time window [ms] in which normalization factor is computed" );

  pdef->addParam( "start_time", "Start time of analysis window [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Start time of analysis window [ms]" );

  pdef->addParam( "end_time", "End time of analysis window [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "End time of analysis window [ms]" );

  pdef->addParam( "percent", "Percentage of normalization", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "How much normalization shall be done [%]? 0%: No normalization, 100%: Full normalization, scale Z trace by normalization factor" );

  pdef->addParam( "equalize", "Parameter for method 'sum':  Equalize RMS amplitudes between P & Z trace before summation", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Perform full trace equalization on Z trace before summation." );
  pdef->addOption( "no", "Do not equalize Z trace before summation." );

  pdef->addParam( "weight_param", "Parameters for frequency-weighted sum", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Receiver depth [m]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Water velocity [m/s]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Use hydrophone data only for N first notches", "Specify number of last hydrophone notch for which hydrohopne data ony shall be used" );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Squeeze factor, valid range [0,1]", "1.0: Equal weights to P & Z, 0.0: Zero weight to Z" );

  pdef->addParam( "weight_output", "Output data", NUM_VALUES_FIXED );
  pdef->addValue( "pz", VALTYPE_OPTION );
  pdef->addOption( "pz", "Output PZ data." );
  pdef->addOption( "p", "Output P data." );
  pdef->addOption( "z", "Output Z data." );

  pdef->addParam( "weight_damping", "Optional 'damping' of weight curves functions...", NUM_VALUES_FIXED );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Damping (1.0: No damping, 0:0: 100% damping, i.e. P and Z are both weighted 1)" );

  pdef->addParam( "sensor_hdr", "Trace header and number values defining P & Z sensors", NUM_VALUES_FIXED );
  pdef->addValue( "sensor", VALTYPE_STRING, "Trace header name");
  pdef->addValue( "1", VALTYPE_NUMBER, "Value identifying P trace");
  pdef->addValue( "5", VALTYPE_NUMBER, "Value identifying Z trace");
}

//--------------------------------------------------------------------------------
//
//
void pz_summation( float* samples_p, float const* samples_z, float percent, int nSamples, int windowLengthSamples,
                   float zScalar, bool doEqualize ) { // float const* rms_z_in, float const* rms_p_in ) {
  float pz_weight = percent / 100.0;

  float* rms_z;
  float* rms_p;
  rms_z = new float[nSamples];
  rms_p = new float[nSamples];

  if( rms_z == NULL || rms_p == NULL ) throw( csException("pz_summation: Out of memory") );

  for( int isamp = 0; isamp < nSamples; isamp++ ) {
    rms_z[isamp] = samples_z[isamp] * samples_z[isamp];
    rms_p[isamp] = samples_p[isamp] * samples_p[isamp];
  }

  float ratio_all = 1.0;

  if( doEqualize ) {
    float rms_all_p = 0.0;
    float rms_all_z = 0.0;
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      rms_all_z += rms_z[isamp];
      rms_all_p += rms_p[isamp];
    }
    rms_all_p = sqrt( rms_all_p / (float)nSamples );
    rms_all_z = sqrt( rms_all_z / (float)nSamples );
    if( rms_all_z != 0.0 ) {
      ratio_all = rms_all_p / rms_all_z;
    }
  }

  zScalar *= ratio_all;

  int windHalf = (int)( windowLengthSamples / 2 );
  int sampWindStart = 0;
  int sampWindEnd   = windowLengthSamples-1;

  // Compute sum of all RMS values over window, for first sample
  double sum_p = 0.0;
  double sum_z = 0.0;
  double ratio = 1.0;
  for( int isampWind = sampWindStart; isampWind <= sampWindEnd; isampWind++ ) {
    sum_p += rms_p[isampWind];
    sum_z += rms_z[isampWind];
  }
  if( sum_z > 0.0 ) {
    ratio = 1.0 + pz_weight * ( sqrt(sum_p/sum_z)/ratio_all - 1.0 );
  }

  //  fprintf( stdout,"Window 0, PZ weight: %10.5f, ratio: %7.5f   %7.5f\n", pz_weight, ratio, ratio_all );
  for( int isamp = 0; isamp < windHalf; isamp++ ) {
    samples_p[isamp] += zScalar * ratio * samples_z[isamp];
  }
  int sampEnd = nSamples-windHalf;
  for( int isamp = windHalf; isamp < sampEnd; isamp++ ) {
    sampWindStart = isamp-windHalf;
    sum_p -= rms_p[sampWindStart];
    sum_z -= rms_z[sampWindStart];
    if( sum_z < 0 ) sum_z = 0.0;
    if( sum_p < 0 ) sum_p = 0.0;
    sampWindStart += 1;
    sampWindEnd   += 1;
    sum_p = sum_p + rms_p[sampWindEnd];
    sum_z = sum_z + rms_z[sampWindEnd];
    if( sum_z > 0.0 ) {
      ratio = 1.0 + pz_weight * ( sqrt(sum_p/sum_z)/ratio_all - 1.0 );
      samples_p[isamp] += zScalar * ratio * samples_z[isamp];
    }
    else {
      ratio = 1.0;
      samples_p[isamp] = 0.0;
    }
    // fprintf( stdout, "%d %d  %d  %f   %f %f  %f %f\n", isamp, sampWindStart, sampWindEnd, ratio, sum_p, sum_z, rms_p[isamp], rms_z[isamp] );
  }
  sampEnd = nSamples-1;
  for( int isamp = nSamples-windHalf; isamp < sampEnd; isamp++ ) {
    samples_p[isamp] += zScalar * ratio * samples_z[isamp];
  }

  delete [] rms_z;
  delete [] rms_p;
}

void pz_zero( float* samples_p, float const* samples_z, int nSamples ) {
  float rms_p = 0.0;
  for( int isamp = 0; isamp < nSamples; isamp++ ) {
    rms_p += samples_p[isamp] * samples_p[isamp];
  }
  rms_p = sqrt( rms_p );

  for( int isamp = 0; isamp < nSamples; isamp++ ) {
    float pz = samples_p[isamp]*samples_z[isamp];
    if( pz < 0.0 ) {  // Only try to recompute pz_weight if last sample has not been reached yet
      samples_p[isamp] = 0.0;
    }
    else {
      samples_p[isamp] = pz;
      if( samples_z[isamp] < 0.0 ) samples_p[isamp] *= -1.0;
    }
  }
  float rms_pz = 0.0;
  for( int isamp = 0; isamp < nSamples; isamp++ ) {
    rms_pz += samples_p[isamp] * samples_p[isamp];
  }
  rms_pz = sqrt( rms_pz );

  if( rms_pz > 0.0 ) {
    float ratio = rms_p / rms_pz;
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      samples_p[isamp] *= ratio;
    }
  }
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_pz_sum_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_pz_sum::VariableStruct* vars = reinterpret_cast<mod_pz_sum::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_pz_sum_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_pz_sum::VariableStruct* vars = reinterpret_cast<mod_pz_sum::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->fftTool != NULL ) {
    delete vars->fftTool;
    vars->fftTool = NULL;
  }
  if( vars->buffer != NULL ) {
    delete [] vars->buffer;
    vars->buffer = NULL;
  }
  if( vars->weightsHyd != NULL ) {
    delete [] vars->weightsHyd;
    vars->weightsHyd = NULL;
  }
  if( vars->weightsGeo != NULL ) {
    delete [] vars->weightsGeo;
    vars->weightsGeo = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_pz_sum_( csParamDef* pdef ) {
  params_mod_pz_sum_( pdef );
}
extern "C" void _init_mod_pz_sum_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_pz_sum_( param, env, writer );
}
extern "C" bool _start_exec_mod_pz_sum_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_pz_sum_( env, writer );
}
extern "C" void _exec_mod_pz_sum_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_pz_sum_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_pz_sum_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_pz_sum_( env, writer );
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csTimeFunction.h"
#include "csTableNew.h"
#include "csInterpolation.h"
#include "csGeolibUtils.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: GAIN
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_gain {
  struct VariableStruct {
    bool applyTGain;
    bool applyAGC;
    bool applyTraceEqualization;
    bool applySphDiv;

    int sphDivOption;
    csTableNew* table;
    int* hdrId_keys;
    int hdrId_offset;
    float* velocityTrace;
    double* keyValueBuffer;
    bool isVelSet;

    float tgain;
    int agcWindowLengthSamples;
    float* buffer;
    float* scalarTGain;
    float traceAmp;
    int outputOption;
  };
  static int const OUTPUT_DATA = 0;
  static int const OUTPUT_GAIN = 1;
  static int const SPHDIV_NONE = 0;
  static int const SPHDIV_OFFSET = 11;
  static int const SPHDIV_SIMPLE = 22;
}
using mod_gain::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_gain_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->applyTGain = false;
  vars->applyAGC   = false;
  vars->applyTraceEqualization = false;
  vars->applySphDiv = false;
  vars->table = NULL;
  vars->hdrId_keys = NULL;
  vars->hdrId_offset = -1;
  vars->velocityTrace = NULL;
  vars->keyValueBuffer = NULL;
  vars->isVelSet = false;
  vars->sphDivOption = mod_gain::SPHDIV_NONE;

  vars->tgain    = 0;
  vars->agcWindowLengthSamples = 0;
  vars->buffer   = NULL;
  vars->scalarTGain   = NULL;
  vars->traceAmp = 0;
  vars->outputOption = mod_gain::OUTPUT_DATA;
//---------------------------------------------
//
  
  std::string tableName;
  if( param->exists("tgain") ) {
    float tgain;
    param->getFloat("tgain", &tgain );
    vars->applyTGain = true;
    vars->scalarTGain = new float[shdr->numSamples];
    float sampleInt_s = shdr->sampleInt/1000.0f;
    vars->scalarTGain[0] = (tgain >= 0) ? 0.0 : 1.0;
    for( int isamp = 1; isamp < shdr->numSamples; isamp++ ) {
      float time = (float)isamp*sampleInt_s;
      vars->scalarTGain[isamp] = pow( time, tgain );
      if( edef->isDebug() ) writer->line("Gain (sample,time,scalar):  %d  %.2f %.6f", isamp, time, vars->scalarTGain[isamp]);
    }
  }
  else {
    vars->applyTGain = false;
  }
  if( param->exists("output") ) {
    std::string text;
    param->getString("output", &text);
    if( !text.compare("data") ) {
      vars->outputOption = mod_gain::OUTPUT_DATA;
    }
    else if( !text.compare("gain") ) {
      vars->outputOption = mod_gain::OUTPUT_GAIN;
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }
  }

  if( param->exists("sph_div") ) {
    std::string text;
    param->getString("sph_div", &text);
    if( !text.compare("simple") ) {
      vars->applySphDiv = true;
      vars->sphDivOption = mod_gain::SPHDIV_SIMPLE;
    }
    else if( !text.compare("offset") ) {
      vars->applySphDiv = true;
      vars->sphDivOption = mod_gain::SPHDIV_OFFSET;
    }
    else if( !text.compare("no") ) {
      vars->applySphDiv = false;
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }
    if( vars->applySphDiv ) {
      std::string tableFilename;
      int colTime;
      int colVel;
      param->getString("table", &tableFilename );
      param->getInt("table_col",&colTime,0);
      param->getInt("table_col",&colVel,1);
      if( colTime < 1 || colVel < 1 ) writer->error("Column numbers in table ('table_col') must be larger than 0");
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
          vars->table->addKey( col-1, interpolate );  // -1 to convert from 'user' column to 'C' column
          if( !hdef->headerExists( headerName ) ) {
            writer->error("No matching trace header found for table key '%s'", headerName.c_str() );
          }
          vars->hdrId_keys[ikey] = hdef->headerIndex( headerName );
        } // END for ikey
      }
      vars->table->addValue( colVel-1 );  // -1 to convert from 'user' column to 'C' column

      bool sortTable = false;
      try {
        vars->table->initialize( tableFilename, sortTable );
      }
      catch( csException& exc ) {
        writer->error("Error when initializing input table '%s': %s\n", text.c_str(), exc.getMessage() );
      }
      vars->keyValueBuffer = new double[vars->table->numKeys()];
      vars->isVelSet = false;
      if( vars->sphDivOption == mod_gain::SPHDIV_OFFSET ) vars->hdrId_offset = hdef->headerIndex("offset");
      vars->velocityTrace = new float[shdr->numSamples];
    } // if applySphDiv
  }

  if( param->exists("trace_equal") ) {
    param->getFloat("trace_equal", &vars->traceAmp );
    vars->applyTraceEqualization = true;
  }
  else {
    vars->applyTraceEqualization = false;
  }

  float window = 0;
  if( param->exists("agc") ) {
    param->getFloat("agc", &window );
    vars->agcWindowLengthSamples = (int)( window/shdr->sampleInt + 0.5 );
    vars->applyAGC = true;
    vars->buffer = new float[shdr->numSamples];
  }
  else {
    vars->applyAGC = false;
  }
  
  int sum = (int)vars->applyTGain + (int)vars->applyAGC + (int)vars->applyTraceEqualization + (int)vars->applySphDiv;
  //  if( (vars->applyTGain && vars->applyAGC) || (vars->applyTGain && vars->applyTraceEqualization) || (vars->applyAGC && vars->applyTraceEqualization) ) {
  if( sum > 1 ) {
    writer->error("More than one gain option specified. Can only specify one gain option at one time.");
  }
  else if( sum == 0 ) { // !vars->applyTGain && !vars->applyAGC && !vars->applyTraceEqualization ) {
    writer->error("No gain option specified.");
  }

  // TEMP:
  //  for( int i = 0; i < hdef->numHeaders(); i++ ) {
  //   writer->line( "GAIN header #%2d: %s %s\n", i, hdef->headerName(i).c_str(), csGeolibUtils::typeText(hdef->headerType(i)) );
  //  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_gain_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  float* samples = trace->getTraceSamples();
  //  float sampleInt_s = shdr->sampleInt/1000.0;

  if( vars->applyTGain ) {
    if( vars->outputOption == mod_gain::OUTPUT_DATA ) {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] *= vars->scalarTGain[isamp];
      }
    }
    else {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] = vars->scalarTGain[isamp];
      }
    }
  }
  else if( vars->applyAGC ) {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      int firstSamp = MAX(0,isamp-vars->agcWindowLengthSamples);
      int lastSamp  = MIN(shdr->numSamples-1,isamp+vars->agcWindowLengthSamples);
      float rms = 0.0;
      for( int isampSum = firstSamp; isampSum <= lastSamp; isampSum++ ) {
        rms += samples[isampSum]*samples[isampSum];
      }
      rms = sqrt( rms/(float)(lastSamp-firstSamp+1) );
      if( rms != 0 ) vars->buffer[isamp] = 1.0f / rms;
      else vars->buffer[isamp] = 0.0f;
    }
    if( vars->outputOption == mod_gain::OUTPUT_DATA ) {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] *= vars->buffer[isamp];
      }
    }
    else {
      memcpy( samples, vars->buffer, shdr->numSamples*sizeof(float) );
    }
  }
  else if( vars->applyTraceEqualization ) {
    float sum = 0.0;
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      sum += fabs( samples[isamp] );
    }
    float mean = sum / ( (float)shdr->numSamples * vars->traceAmp );
    if( vars->outputOption == mod_gain::OUTPUT_DATA ) {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] /= mean;
      }
    }
    else {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] = 1/mean;
      }
    }
  }
  else if( vars->applySphDiv ) {
    csTimeFunction<double> const* velTimeFunc = NULL;
    if( vars->table->numKeys() > 0 ) {
      double* keyValueBuffer = new double[vars->table->numKeys()];
      for( int ikey = 0; ikey < vars->table->numKeys(); ikey++ ) {
        keyValueBuffer[ikey] = trace->getTraceHeader()->doubleValue( vars->hdrId_keys[ikey] );
      }
      if( vars->isVelSet ) { // If velocity trace has been set previously, check if this trace has different key values
        for( int ikey = 0; ikey < vars->table->numKeys(); ikey++ ) {
          if( keyValueBuffer[ikey] != vars->keyValueBuffer[ikey] ) {
            vars->isVelSet = false;
            break;
          }
        }
      }
      memcpy( vars->keyValueBuffer, keyValueBuffer, vars->table->numKeys() );
      if( !vars->isVelSet ) velTimeFunc = vars->table->getFunction( keyValueBuffer, false );
      delete [] keyValueBuffer;
    }
    else if( !vars->isVelSet ) {
      velTimeFunc = vars->table->getFunction( NULL, false );
    }
    if( !vars->isVelSet ) {
      int numVelocities = velTimeFunc->numValues();
      float* t0      = new float[numVelocities];
      float* vel_rms = new float[numVelocities];
      for( int i = 0; i < numVelocities; i++ ) {
        t0[i]      = velTimeFunc->timeAtIndex(i)/1000.0;  // Convert to seconds
        vel_rms[i] = velTimeFunc->valueAtIndex(i);
      }
      // Linearly interpolate input velocities
      csInterpolation::linearInterpolation( numVelocities, t0, vel_rms, shdr->numSamples, shdr->sampleInt/1000.0, vars->velocityTrace );
      vars->isVelSet = true;
      delete [] t0;
      delete [] vel_rms;
    }
    float vel0     = vars->velocityTrace[0];
    float vel0_sqr = vel0 * vel0;
    if( vars->sphDivOption == mod_gain::SPHDIV_SIMPLE ) {
      // t * V(t)^2 / V(0)
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        float time_s  = (float)isamp * shdr->sampleInt / 1000.0f;
        float vel = vars->velocityTrace[isamp];
        float vel_sqr = vel * vel;
        float sphDivGain = time_s * (vel_sqr / vel0_sqr);
        if( vars->outputOption == mod_gain::OUTPUT_DATA ) {
          samples[isamp] *= sphDivGain;
        }
        else {
          samples[isamp] = sphDivGain;
        }
        //        fprintf(stdout,"%d %f  vel: %f, gain: %f\n", isamp, time_s, vel, sphDivGain );
      }
    }
    else { // Offset dependeint spherical divergence correction:
      float offset = trace->getTraceHeader()->floatValue( vars->hdrId_offset );
      float offset_sqr = offset * offset;
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        float time_s  = (float)isamp * shdr->sampleInt / 1000.0f;
        float time0_s = time_s; //NMO zero-offset time;   // !CHANGE! Compute zero-offset time using inverse NMO
        if( time0_s <= 0 ) time0_s = shdr->sampleInt/1000.0;
        float indexVel = 1000*time0_s / shdr->sampleInt;
        int index1 = (int)indexVel;
        float vel = 0;
        if( index1 >= 0 && index1 < shdr->numSamples-1 ) {
          vel = vars->velocityTrace[index1] + (indexVel-(float)index1) * ( vars->velocityTrace[index1+1] - vars->velocityTrace[index1] );
        }
        else if( index1 < 0 ) {
          vel = vars->velocityTrace[0];
        }
        else if( index1 >= shdr->numSamples-1 ) {
          vel = vars->velocityTrace[shdr->numSamples-1];
        }
        float vel_sqr = vel * vel;
        //    a = ( v(T)^2  - v(0)^2   ) * offset^2   / ( T^2             * v(T)^4 )
        float a = ( vel_sqr - vel0_sqr ) * offset_sqr / ( time0_s*time0_s * vel_sqr * vel_sqr  );
        // gain          = t      * ( v(T)^2 / v(0)^2 ) * sqrt( 1 + a )
        float sphDivGain = time_s * (vel_sqr / vel0_sqr)  * sqrt( 1.0 + a );
        if( vars->outputOption == mod_gain::OUTPUT_DATA ) {
          samples[isamp] *= sphDivGain;
        }
        else {
          samples[isamp] = sphDivGain;
        }
        fprintf(stdout,"%d %f  vel: %f, offset: %f, a: %f, gain: %f\n", isamp, time_s, vel, offset, a, sphDivGain );
      }
    }

  } // END: if applySphDiv
  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_gain_( csParamDef* pdef ) {
  pdef->setModule( "GAIN", "Apply gain function to trace samples" );

  pdef->addParam( "tgain", "Apply time gain function", NUM_VALUES_FIXED );
  pdef->addValue( "2", VALTYPE_NUMBER, "Time gain value: gain = t^value" );

  pdef->addParam( "agc", "Apply automatic gain control (AGC)", NUM_VALUES_FIXED );
  pdef->addValue( "500", VALTYPE_NUMBER, "AGC sliding window length [ms]" );

  pdef->addParam( "trace_equal", "Apply full trace equalisation", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Average output amplitude" );

  pdef->addParam( "sph_div", "Apply offset-dependent spherical divergence correction", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not apply spherical divergence correction" );
  pdef->addOption( "offset", "Apply offset-dependent spherical divergence correction" );
  pdef->addOption( "simple", "Apply simple spherical divergence correction.", "gain = t * vel(t)^2 / vel(0)^2" );

  // gain = t * v(T)^2 / v(0) * sqrt( 1 + a )
  // a    = ( v(T)^2 - v(0)^2 ) * offset^2 / ( T^2 * v(T)^4 )
  // T    = t0 = zero-offset time (NMO)

  pdef->addParam( "table", "Velocity(RMS) table for spherical divergence correction", NUM_VALUES_FIXED,
                  "Velocity table format: The ASCII file should contain only numbers, no text. Required are two columns containing time and velocity pairs. Optional: Up to 2 key values. Lines starting with '#' are considered comment lines." );
  pdef->addValue( "", VALTYPE_STRING, "Filename containing velocity table.");

  pdef->addParam( "table_col", "Table columns containing time in [ms] and velocity in [m/s]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing time [ms]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing velocity [m/s]" );

  pdef->addParam( "table_key", "Key trace header used to match values found in specified table columns", NUM_VALUES_VARIABLE,
                  "Specify the 'key' parameter for each key in the velocity file" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name of key header" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table" );
  pdef->addValue( "yes", VALTYPE_OPTION, "Interpolate based to this key?" );
  pdef->addOption( "yes", "Use this key for interpolation of value" );
  pdef->addOption( "no", "Do not use this key for interpolation", "The input table is expected to contain the exact key values for this trace header" );

  pdef->addParam( "output", "Output data", NUM_VALUES_FIXED );
  pdef->addValue( "data", VALTYPE_OPTION );
  pdef->addOption( "data", "Output data with specified gain applied" );
  pdef->addOption( "gain", "Output gain values" );

  //  pdef->addParam( "mode", "Mode of application", NUM_VALUES_FIXED );
  //  pdef->addValue( "apply", VALTYPE_OPTION );
  //  pdef->addOption( "apply", "Apply gain." );
  //  pdef->addOption( "remove", "Remove gain = Inverse application (if possible)." );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_gain_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_gain::VariableStruct* vars = reinterpret_cast<mod_gain::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_gain_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_gain::VariableStruct* vars = reinterpret_cast<mod_gain::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->scalarTGain ) {
    delete [] vars->scalarTGain;
    vars->scalarTGain = NULL;
  }
  if( vars->buffer ) {
    delete [] vars->buffer;
    vars->buffer = NULL;
  }
  if( vars->keyValueBuffer ) {
    delete [] vars->keyValueBuffer;
    vars->keyValueBuffer = NULL;
  }
  if( vars->velocityTrace ) {
    delete [] vars->velocityTrace;
    vars->velocityTrace = NULL;
  }
  if( vars->hdrId_keys ) {
    delete [] vars->hdrId_keys;
    vars->hdrId_keys = NULL;
  }
  if( vars->table ) {
    delete vars->table;
    vars->table = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_gain_( csParamDef* pdef ) {
  params_mod_gain_( pdef );
}
extern "C" void _init_mod_gain_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_gain_( param, env, writer );
}
extern "C" bool _start_exec_mod_gain_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_gain_( env, writer );
}
extern "C" void _exec_mod_gain_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_gain_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_gain_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_gain_( env, writer );
}

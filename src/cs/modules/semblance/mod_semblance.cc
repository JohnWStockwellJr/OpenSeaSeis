/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csNMOCorrection.h"
#include "csTableManagerNew.h"
#include "csTableAll.h"
#include "csVector.h"
#include "csInterpolation.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SEMBLANCE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_semblance {
  struct VariableStruct {
    csNMOCorrection* nmo;
    int hdrId_offset;
    int hdrId_rec_z;
    int hdrId_sou_z;
    int hdrId_vel_rms;
    int outputOption;

    int windowLengthSamples;
    float* velStart;  // Start velocity at each sample
    float* velEnd;    // End velocity at each sample
    float velInc;
    float velMin;
    float velMax;
    int numVels;

    float* bufferSemblance;

    bool isNMO;
    cseis_geolib::csInterpolation* lmoInterpol;
    float lmoRefVel;  // LMO reference velocity
    float lmoRefVelInverse;

    csTableManagerNew* tableManager;
  };
  static int const OUTPUT_FIRST   = 1;
  static int const OUTPUT_LAST    = 2;
  static int const OUTPUT_AVERAGE = 3;
}
using namespace mod_semblance;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_semblance_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->nmo            = NULL;
  vars->hdrId_offset   = -1;
  vars->hdrId_vel_rms  = -1;
  vars->hdrId_rec_z    = -1;
  vars->hdrId_sou_z   = -1;
  vars->outputOption   = -1;
  vars->windowLengthSamples = 0;
  vars->velStart  = NULL;
  vars->velEnd    = NULL;
  vars->velInc    = 0;
  vars->velMin    = 0;
  vars->velMax    = 0;
  vars->numVels   = 0;
  vars->bufferSemblance  = NULL;
  vars->tableManager = NULL;
  vars->isNMO     = true;
  vars->lmoInterpol  = NULL;
  vars->lmoRefVel = 0;
  vars->lmoRefVelInverse = 0;

  //-------------------------------------
  std::string text;
  int mode_nmo = csNMOCorrection::PP_NMO;
  if( param->exists("wave_mode") ) {
    param->getString( "wave_mode", &text );
    if( !text.compare("pp_iso") ) {
      mode_nmo = csNMOCorrection::PP_NMO;
    }
    else if( !text.compare("ps_iso") ) {
      mode_nmo = csNMOCorrection::PS_NMO;
    }
    else if( !text.compare("p_direct") ) {
      vars->isNMO = false;
      writer->line("Linear moveout correction is performed. Source/receiver depths are read in from trace headers 'rec_z' and 'sou_z'");
      if( param->exists("lmo_refvel") ) {
        param->getFloat("lmo_refvel", &vars->lmoRefVel );
      }
      if( vars->lmoRefVel != 0 ) {
        vars->lmoRefVelInverse = 1.0/vars->lmoRefVel;
      }
    }
    else {
      writer->error("Output option not recognized: %s.", text.c_str());
    }
  }

  //-------------------------------------
  param->getFloat( "vel_inc", &vars->velInc, 0 );

  vars->velStart = new float[shdr->numSamples];
  vars->velEnd   = new float[shdr->numSamples];

  int numVelTimes = param->getNumLines( "vel_range" );
  if( numVelTimes < 1 ) {
    writer->error("User parameter 'vel_range' missing." );
  }
  float* velTime  = new float[numVelTimes];
  float* velStart = new float[numVelTimes];
  float* velEnd   = new float[numVelTimes];
  float velMin = 10e20;
  float velMax = -10e20;
  for( int i = 0; i < numVelTimes; i++ ) {
    param->getFloat( "vel_range", &velTime[i], 0 );
    param->getFloat( "vel_range", &velStart[i], 1 );
    param->getFloat( "vel_range", &velEnd[i], 2 );
    if( velEnd[i] <= velStart[i] || vars->velInc > (velEnd[i]-velStart[i]) ) {
      writer->error("Inconsistent velocity range specified at time %f: Start/End/Inc   %f/%f/%f",
        velTime[i], velStart[i], velEnd[i], vars->velInc );
    }
    if( velStart[i] < velMin ) velMin = velStart[i];
    if( velEnd[i] > velMax )   velMax = velEnd[i];
  }
  vars->velMin = velMin;
  vars->velMax = velMax;

  int isampStart = (int)(velTime[0]/shdr->sampleInt + 0.5);
  int isampEnd   = (int)(velTime[numVelTimes-1]/shdr->sampleInt + 0.5);
  for( int isamp = 0; isamp < isampStart; isamp++ ) {
    vars->velStart[isamp] = velStart[0];
    vars->velEnd[isamp]   = velEnd[0];
  }
  int ivel = 0;
  for( int isamp = isampStart; isamp < isampEnd; isamp++ ) {
    float time = isamp*shdr->sampleInt;
    while( ivel < numVelTimes-1 && velTime[ivel+1] < time ) {
      ivel += 1;
    }
    float weight = (time-velTime[ivel])/(velTime[ivel+1]-velTime[ivel]);
    vars->velStart[isamp] = velStart[ivel] + weight*( velStart[ivel+1] - velStart[ivel] );
    vars->velEnd[isamp]   = velEnd[ivel] + weight*( velEnd[ivel+1] - velEnd[ivel] );
  }
  for( int isamp = isampEnd; isamp < shdr->numSamples; isamp++ ) {
    vars->velStart[isamp] = velStart[numVelTimes-1];
    vars->velEnd[isamp]   = velEnd[numVelTimes-1];
  }

  vars->numVels = (int)( (vars->velMax - vars->velMin) / vars->velInc ) + 1;

  //-----------------------------------------------
  //
  float windowTime = 0;
  param->getFloat( "window", &windowTime, 0 );
  vars->windowLengthSamples = (int)(windowTime/shdr->sampleInt) + 1;

  //-----------------------------------------------
  if( param->exists("mute_table") ) {
    param->getString("mute_table", &text );

    try {
      vars->tableManager = new csTableManagerNew( text, csTableAll::TABLE_TYPE_UNIQUE_KEYS, hdef );
      if( vars->tableManager->valueName().compare("time") ) {
        writer->error("Mute table must contain 'value' column labelled 'time', for example: '@%s time'. Value label found: '%s'",
                   vars->tableManager->numKeys() > 0 ? vars->tableManager->keyName(0).c_str() : "offset",
                   vars->tableManager->valueName().c_str() );
      }
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s': %s\n", text.c_str(), exc.getMessage() );
    }
    if( vars->tableManager->numKeys() < 1 ) {
      writer->error("Number of table keys = %d. Specify table key by placing the character '%c' in front of the key name. Example: %csource time vel  (source is a table key)",
                 csTableAll::KEY_CHAR, csTableAll::KEY_CHAR);
    }

    if( edef->isDebug() ) {
      vars->tableManager->dump();
    }
  }
  //-----------------------------------------------
  if( param->exists("output_hdr") ) {
    std::string text;
    param->getString("output_hdr", &text);
    if( !text.compare("first") ) {
      vars->outputOption = OUTPUT_FIRST;
    }
    else if( !text.compare("last") ) {
      vars->outputOption = OUTPUT_LAST;
    }
    else if( !text.compare("average") ) {
      writer->error("Option AVERAGE is not supported yet...");
      vars->outputOption = OUTPUT_AVERAGE;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
  }
  else {
    vars->outputOption = OUTPUT_AVERAGE;
  }

  //-----------------------------------------------
  csVector<std::string> valueList;

  if( vars->isNMO ) {
    vars->nmo = new csNMOCorrection( shdr->sampleInt, shdr->numSamples, mode_nmo );
  }
  else {
    vars->hdrId_rec_z = hdef->headerIndex( "rec_z" );
    vars->hdrId_sou_z = hdef->headerIndex( "sou_z" );
    vars->lmoInterpol = new csInterpolation( shdr->numSamples, shdr->sampleInt, 8 );
  }

  if( !hdef->headerExists( "offset" ) ) {
    writer->error("Trace header 'offset' does not exist. Cannot perform NMO correction.");
  }
  else if( hdef->headerType( "offset" ) != TYPE_FLOAT && hdef->headerType( "offset" ) != TYPE_DOUBLE ) {
    writer->error("Trace header 'offset' exists but has the wrong number type. Should be FLOAT.");
  }
  vars->hdrId_offset = hdef->headerIndex( "offset" );
  if( !hdef->headerExists( "vel_rms" ) ) {
    hdef->addHeader( TYPE_FLOAT, "vel_rms", "RMS velocity [m/s]" );
  }
  vars->hdrId_vel_rms = hdef->headerIndex( "vel_rms" );

  vars->bufferSemblance = new float[vars->numVels*shdr->numSamples];

  if( edef->isDebug() ) {
    writer->line("Number of velocities: %d   %f %f %f", vars->numVels, vars->velMin, vars->velMax, vars->velInc );
    writer->line("Number of samples in window: %d", vars->windowLengthSamples );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_semblance_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;


  int nTracesIn = traceGather->numTraces();

  // Step 1: Determine mute time. Decimate input data if mute kills full trace
  csVector<int> muteList(nTracesIn);
  if( vars->tableManager != NULL ) {
    for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
      double muteTime = vars->tableManager->getValue( traceGather->trace(itrc)->getTraceHeader() );
      int sampleIndex = (int)( muteTime / shdr->sampleInt + 0.5 );
      muteList.insertEnd( std::max( 0, sampleIndex ) );
    }
    for( int itrc = nTracesIn-1; itrc >= 0; itrc-- ) {
      int muteSample = muteList.at(itrc);
      if( edef->isDebug() ) writer->line(" Trace #%d: mute sample index : %d  (Num samples: %d)", itrc, muteSample, shdr->numSamples );
      if( muteSample >= shdr->numSamples-2 ) {
        if( edef->isDebug() ) writer->line("    ...remove trace");
        muteList.remove(itrc);
        traceGather->freeTrace(itrc);
      }
    }
    if( edef->isDebug() ) writer->line("\nReduced input data from %d traces to %d traces\n", nTracesIn, traceGather->numTraces());
    nTracesIn = traceGather->numTraces();
  }
  int* sampleMuteEnd = new int[nTracesIn];
  if( vars->tableManager != NULL ) {
    for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
      sampleMuteEnd[itrc] = muteList.at(itrc);
    }
    muteList.clear();
  }
  else {
    for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
      sampleMuteEnd[itrc] = 0;
    }
  }

  // !CHANGE! Selection of output trace not coded yet
  /*  int outTrace   = 0;
  if( vars->outputOption == OUTPUT_LAST ) {
    outTrace = nTracesIn-1;
  }
  else if( vars->outputOption == OUTPUT_FIRST ) {
    outTrace = 0;
  }
  */
  traceGather->createTraces( nTracesIn, nTracesIn, env->headerDef, shdr->numSamples );

  float** samplesTmpNMO = new float*[nTracesIn];
  for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
    samplesTmpNMO[itrc] = traceGather->trace(itrc+nTracesIn)->getTraceSamples();
  }

  float* offset = new float[nTracesIn];
  for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
    csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
    offset[itrc]          = trcHdr->doubleValue( vars->hdrId_offset );
  }
  if( !vars->isNMO ) {
    for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
      csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
      double rec_z          = trcHdr->doubleValue( vars->hdrId_rec_z );
      double sou_z          = trcHdr->doubleValue( vars->hdrId_sou_z );
      offset[itrc]          = sqrt( pow(offset[itrc],2) + pow(rec_z-sou_z,2) );
    }
  }
  //--------------------------------------------------------------------
  //
  //
/*
    // TEMP: Compute simple stack response
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      int firstSamp = std::max( 0, isamp-vars->windowLengthSamples );
      int lastSamp  = std::min( shdr->numSamples-1, isamp+vars->windowLengthSamples );
      double sum_stack = 0;
      for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
        for( int isampSum = lastSamp; isampSum <= lastSamp; isampSum++ ) {
          sum_stack += samplesTmpNMO[itrc][isamp];
        }
      }
      semblancePtr[isamp] = sum_stack*sum_stack / nTracesIn;
    }*/

  float time = 0.0;
  for( int ivel = 0; ivel < vars->numVels; ivel++ ) {
    float velocity  = vars->velMin + (float)ivel * vars->velInc;
    if( edef->isDebug() ) writer->line("Compute semblance for velocity #%d: %f", ivel, velocity);
    for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
      float const* samplesIn = traceGather->trace(itrc)->getTraceSamples();
      //      double offset          = traceGather->trace(itrc)->getTraceHeader()->doubleValue( vars->hdrId_offset );
      if( vars->isNMO ) {
        vars->nmo->perform_nmo( samplesIn, 1, &time, &velocity, offset[itrc], samplesTmpNMO[itrc] );
      }
      else {
        float shift_ms = -offset[itrc]*1000.0 * ( (1.0/(vars->lmoRefVel+velocity)) - vars->lmoRefVelInverse ); 
        vars->lmoInterpol->static_shift( shift_ms, samplesIn, samplesTmpNMO[itrc] );
      }
    }
    float* semblancePtr = &vars->bufferSemblance[ivel*shdr->numSamples];

    // Compute semblance

    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      if( vars->velStart[isamp] > velocity || vars->velEnd[isamp] < velocity ) {
        semblancePtr[isamp] = 0.0;
        continue;
      }
      int firstSamp = std::max( 0, isamp-vars->windowLengthSamples );
      int lastSamp  = std::min( shdr->numSamples-1, isamp+vars->windowLengthSamples );
      double sumUpper = 0.0;
      double sumLower = 0.0;
      for( int isampSum = firstSamp; isampSum <= lastSamp; isampSum++ ) {
        double sum    = 0.0;
        double sum_sq = 0.0;
        for( int itrc = 0; itrc < nTracesIn; itrc++ ) {
          if( isampSum <= sampleMuteEnd[itrc] ) continue;
          sum    += samplesTmpNMO[itrc][isampSum];
          sum_sq += samplesTmpNMO[itrc][isampSum]*samplesTmpNMO[itrc][isampSum];
        }
        sumUpper += sum*sum;
        sumLower += sum_sq;
      }
      if( sumLower == 0 ) semblancePtr[isamp] = 0.0;
      else semblancePtr[isamp] = sumUpper / sumLower;
    }
  }

  if( vars->numVels > 2*nTracesIn ) {
    traceGather->createTraces( 2*nTracesIn, vars->numVels-2*nTracesIn, env->headerDef, shdr->numSamples );
  }
  else {
    traceGather->freeTraces( vars->numVels, 2*nTracesIn-vars->numVels );
  }
  for( int ivel = 0; ivel < vars->numVels; ivel++ ) {
    float velocity  = vars->velMin + (float)ivel * vars->velInc;
    csTraceHeader* trcHdr = traceGather->trace(ivel)->getTraceHeader();
    trcHdr->setFloatValue( vars->hdrId_vel_rms, velocity+vars->lmoRefVel );
    float* samplesOut  = traceGather->trace(ivel)->getTraceSamples();
    memcpy( samplesOut, &vars->bufferSemblance[ivel*shdr->numSamples], shdr->numSamples*sizeof(float) );
  }

  delete [] samplesTmpNMO;
  delete [] sampleMuteEnd;
  delete [] offset;
//  *numTrcToKeep = vars->num_vels;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_semblance_( csParamDef* pdef ) {
  pdef->setModule( "SEMBLANCE", "Semblance panel generation" );

  pdef->addParam( "vel_range", "Velocity range to test", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Time at which specified velocity applies [ms]", "Velocity range is interpolated linearly over time" );
  pdef->addValue( "700", VALTYPE_NUMBER, "Start stacking velocity [m/s]" );
  pdef->addValue( "2000", VALTYPE_NUMBER, "End stacking velocity [m/s]" );

  pdef->addParam( "vel_inc", "Velocity increment", NUM_VALUES_FIXED );
  pdef->addValue( "20", VALTYPE_NUMBER, "Velocity increment between tests [m/s]" );

  pdef->addParam( "vel_test", "Maximum deviation from velocity guide function", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Maximum deviation to velocity guide function to test [%]", "'0' means all velocities will be tested" );

  pdef->addParam( "window", "Window length to use for semblance calculation", NUM_VALUES_FIXED );
  pdef->addValue( "100", VALTYPE_NUMBER, "Window length [ms]" );

  pdef->addParam( "wave_mode", "Mode of wave to model", NUM_VALUES_FIXED );
  pdef->addValue( "pp_iso", VALTYPE_OPTION );
  pdef->addOption( "pp_iso", "Isotropic PP mode (downgoing P, upgoing P)" );
  pdef->addOption( "ps_iso", "Isotropic PS mode (downgoing P, upgoing S)" );
  pdef->addOption( "p_direct", "Direct P arrival, straight ray assumption. Essentially to Linear Moveout (LMO)", "For this option, use parameter 'lmo_refvel' to set up a reference velocity in case the input data is already LMO corrected using a constant velocity" );

  pdef->addParam( "lmo_refvel", "Linear moveout reference velocity", NUM_VALUES_FIXED, "Use this parameter to specify the LMO velocity which is already applied to the input data. If a reference velocity other than 0 is specified, the tested velocities are interpreted as differences to the reference velocity, for example +/-10m/s" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Velocity [m/s]", "= 0: Input data is not LMO corrected. If set != 0, the specified semblance velocity range should be a relative velocity difference." );

  pdef->addParam( "mute_table", "Mute table", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Mute table file.",
    "The mute table must have at least two columns, one giving a key and the second giving the mute time in [ms]" );

  pdef->addParam( "output_hdr", "How shall trace header of output trace be determined?", NUM_VALUES_FIXED );
  pdef->addValue( "first", VALTYPE_OPTION );
  pdef->addOption( "first", "Output trace header values of first input trace to output trace" );
  pdef->addOption( "last", "Output trace header values of last input trace to output trace" );
//  pdef->addOption( "average", "Average all input trace headers for output trace" );
}



//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_semblance_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_semblance::VariableStruct* vars = reinterpret_cast<mod_semblance::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_semblance_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_semblance::VariableStruct* vars = reinterpret_cast<mod_semblance::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->bufferSemblance ) {
    delete [] vars->bufferSemblance;
    vars->bufferSemblance = NULL;
  }
  if( vars->tableManager != NULL ) {
    delete vars->tableManager;
    vars->tableManager = NULL;
  }
  if( vars->velStart != NULL ) {
    delete [] vars->velStart;
    vars->velStart = NULL;
  }
  if( vars->velEnd != NULL ) {
    delete [] vars->velEnd;
    vars->velEnd = NULL;
  }
  if( vars->nmo != NULL ) {
    delete vars->nmo;
    vars->nmo = NULL;
  }
  if( vars->lmoInterpol != NULL ) {
    delete vars->lmoInterpol;
    vars->lmoInterpol = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_semblance_( csParamDef* pdef ) {
  params_mod_semblance_( pdef );
}
extern "C" void _init_mod_semblance_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_semblance_( param, env, writer );
}
extern "C" bool _start_exec_mod_semblance_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_semblance_( env, writer );
}
extern "C" void _exec_mod_semblance_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_semblance_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_semblance_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_semblance_( env, writer );
}

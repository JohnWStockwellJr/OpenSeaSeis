/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: BEAM_FORMING
 *
 * @author Bjorn Olofsson
 * @date   2010
 */
namespace mod_beam_forming {
  struct VariableStruct {
    int hdrId_rcv;
    int hdrId_rec_x;
    int hdrId_rec_y;

    int hdrId_azim;
    int hdrId_slowness;
    int hdrId_slowness_x;
    int hdrId_slowness_y;
    int hdrId_beam_sou_x;
    int hdrId_beam_sou_y;

    int numAzimuths;
    float minSlowness; // s/km
    float maxSlowness; // s/km
    int numSlownesses; // s/km

    float ps_slowness; // s/km
    float ps_minOffset; // m
    float ps_maxOffset; // m
    bool isOffsetSet;
    float ps_gain;
    double minX;
    double minY;
    double maxX;
    double maxY;
    double incXY;
    int numStepsX;
    int numStepsY;

    int startSamp;
    int endSamp;
    int numSamples_in;
    int method;
    int outputOption;
  };
  static int const METHOD_PLANE_WAVE_FAN  = 1;
  static int const METHOD_PLANE_WAVE_GRID = 2;
  static int const METHOD_POINT_SOURCE    = 3;
  static int const WINDOW_CUT  = 11;
  static int const WINDOW_ZERO = 12;
}

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_beam_forming_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  mod_beam_forming::VariableStruct* vars = new mod_beam_forming::VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->hdrId_rcv   = -1;
  vars->hdrId_rec_x = -1;
  vars->hdrId_rec_y = -1;
  vars->hdrId_beam_sou_x = -1;
  vars->hdrId_beam_sou_y = -1;

  vars->hdrId_azim     = -1;
  vars->hdrId_slowness = -1;
  vars->numSlownesses = 0;
  vars->numAzimuths   = 0;
  vars->numStepsX   = 0;
  vars->numStepsY   = 0;

  vars->hdrId_slowness_y = -1;
  vars->hdrId_slowness_x = -1;

  vars->ps_slowness    = 0.0;
  vars->ps_minOffset   = 0.0;
  vars->ps_maxOffset   = 0.0;
  vars->isOffsetSet    = false;
  vars->ps_gain        = 0.0;
  vars->minX = 0.0;
  vars->minY = 0.0;
  vars->maxX = 0.0;
  vars->maxY = 0.0;
  vars->incXY = 0.0;
  vars->startSamp = 0;
  vars->endSamp   = shdr->numSamples - 1;
  vars->outputOption = mod_beam_forming::WINDOW_ZERO;
  vars->numSamples_in = shdr->numSamples;

  //------------------------------------------------------------------------------
  if( param->exists("plane_wave") ) {
    vars->method = mod_beam_forming::METHOD_PLANE_WAVE_GRID;
    param->getFloat("plane_wave", &vars->minSlowness, 0);
    param->getFloat("plane_wave", &vars->maxSlowness, 1);
    param->getInt("plane_wave", &vars->numSlownesses, 2);

    vars->hdrId_slowness = hdef->addHeader( TYPE_FLOAT, "beam_slowness", "Beam slowness [s/km]" );
    vars->hdrId_slowness_x = hdef->addHeader( TYPE_FLOAT, "beam_slowness_x", "Beam slowness X component [s/km]" );
    vars->hdrId_slowness_y = hdef->addHeader( TYPE_FLOAT, "beam_slowness_y", "Beam slowness Y component [s/km]" );
  }
  else if( param->exists("plane_wave_fan") ) {
    vars->method = mod_beam_forming::METHOD_PLANE_WAVE_FAN;
    param->getFloat("plane_wave_fan", &vars->minSlowness, 0);
    param->getFloat("plane_wave_fan", &vars->maxSlowness, 1);
    param->getInt("plane_wave_fan", &vars->numSlownesses, 2);
    param->getInt("plane_wave_fan", &vars->numAzimuths, 3);

    vars->hdrId_azim = hdef->addHeader( TYPE_FLOAT, "beam_azim", "Beam azimuth [deg]" );
    vars->hdrId_slowness = hdef->addHeader( TYPE_FLOAT, "beam_slowness", "Beam slowness [s/km]" );
  }
  else if( param->exists("point_source") ) {
    vars->method = mod_beam_forming::METHOD_POINT_SOURCE;
    param->getDouble("point_source", &vars->minX, 0);
    param->getDouble("point_source", &vars->minY, 1);
    param->getDouble("point_source", &vars->maxX, 2);
    param->getDouble("point_source", &vars->maxY, 3);
    param->getDouble("point_source", &vars->incXY, 4);
    vars->hdrId_beam_sou_x = hdef->addHeader( TYPE_FLOAT, "beam_sou_x", "Beam point source X [s/km]" );
    vars->hdrId_beam_sou_y = hdef->addHeader( TYPE_FLOAT, "beam_sou_y", "Beam point source Y [s/km]" );
    vars->numStepsX = (int)round((vars->maxX - vars->minX) / vars->incXY + 1.0 );
    vars->numStepsY = (int)round((vars->maxY - vars->minY) / vars->incXY + 1.0 );
    param->getFloat("point_source_slowness", &vars->ps_slowness);
    if( param->exists("point_source_offset") ) {
      param->getFloat("point_source_offset", &vars->ps_minOffset, 0);
      param->getFloat("point_source_offset", &vars->ps_maxOffset, 1);
      if( vars->ps_maxOffset <= vars->ps_minOffset ) {
        writer->error("Inconsistent min/max point source offset given: min=%f > max=%f", vars->ps_minOffset, vars->ps_maxOffset);
      }
      vars->isOffsetSet = true;
    }
    if( param->exists("point_source_gain") ) {
      param->getFloat("point_source_gain", &vars->ps_gain);
    }
    writer->line("Point source search grid (xmin,xmax,dx,n): %f %f %f %d", vars->minX, vars->maxX, vars->incXY, vars->numStepsX);
    writer->line("Point source search grid (ymin,ymax,dy,n): %f %f %f %d", vars->minY, vars->maxY, vars->incXY, vars->numStepsY);
  }
  else {
    writer->error("Specify parameter 'plane_wave' or 'point_source'.");
  }

  if( param->exists("window") ) {
    float timeStart = 0.0;
    float timeEnd = 0.0;
    param->getFloat("window",&timeStart,0);
    param->getFloat("window",&timeEnd,0);
    vars->startSamp = (int)round(timeStart/shdr->sampleInt);
    vars->endSamp   = (int)round(timeEnd/shdr->sampleInt);
    if( vars->startSamp < 0 ) vars->startSamp = 0;
    if( vars->endSamp > shdr->numSamples-1 ) vars->endSamp = shdr->numSamples-1;
    if( vars->startSamp >= vars->endSamp ) {
      writer->error("Inconsistent window: %f to %f", timeStart, timeEnd);
    }

    if( param->getNumValues("window") > 2 ) {
      std::string text;
      param->getString("window",&text,2);
      if( !text.compare("zero") ) {
        vars->outputOption = mod_beam_forming::WINDOW_ZERO;
      }
      else if( !text.compare("cut") ) {
        vars->outputOption = mod_beam_forming::WINDOW_CUT;
        shdr->numSamples = vars->endSamp - vars->startSamp + 1;
      }
      else {
        writer->error("Unknown option: %s", text.c_str());
      }
    }
  }

//  vars->hdrId_rcv   = hdef->addHeader( csStandardHeaders::get("rcv") );
//  vars->hdrId_rec_x = hdef->addHeader( csStandardHeaders::get("rec_x") );
//  vars->hdrId_rec_y = hdef->addHeader( csStandardHeaders::get("rec_y") );
//  vars->hdrId_rcv   = hdef->headerIndex( "rcv" );
  vars->hdrId_rec_x = hdef->headerIndex( "rec_x" );
  vars->hdrId_rec_y = hdef->headerIndex( "rec_y" );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_beam_forming_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  mod_beam_forming::VariableStruct* vars = reinterpret_cast<mod_beam_forming::VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;

  int nTraces = traceGather->numTraces();
  int nBeams = 0;
  if( vars->method == mod_beam_forming::METHOD_PLANE_WAVE_FAN ) {
    nBeams = vars->numSlownesses * vars->numAzimuths;
  }
  else if( vars->method == mod_beam_forming::METHOD_PLANE_WAVE_GRID ) {
    nBeams = vars->numSlownesses * vars->numSlownesses;
  }
  else {
    nBeams = vars->numStepsX * vars->numStepsY;
  }

  float* rec_dx = new float[nTraces];
  float* rec_dy = new float[nTraces];
  rec_dx[0] = 0.0;
  rec_dy[0] = 0.0;
  for( int itrc = 1; itrc < nTraces; itrc++ ) {
    csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
    rec_dx[itrc] = (float)(trcHdr->doubleValue(vars->hdrId_rec_x) - rec_dx[0]);
    rec_dy[itrc] = (float)(trcHdr->doubleValue(vars->hdrId_rec_y) - rec_dy[0]);
  }
  traceGather->createTraces( nTraces, nBeams, hdef, shdr->numSamples );
  for( int ibeam = 0; ibeam < nBeams; ibeam++ ) {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      float* outSamples = traceGather->trace(nTraces+ibeam)->getTraceSamples();
      outSamples[isamp] = 0.0;
    }
  }
  int beamCounter = 0;
  int sampShift = 0;
  if( vars->outputOption == mod_beam_forming::WINDOW_CUT ) {
    sampShift = vars->startSamp;
  }
//  fprintf(stderr,"%d %d %d, num samples: %d %d\n", vars->startSamp, sampShift, vars->outputOption, shdr->numSamples, vars->numSamples_in);
  if( vars->method == mod_beam_forming::METHOD_PLANE_WAVE_FAN ) {
    float azimStep = (float)( ( 360.0 / (double)vars->numAzimuths ) * M_PI / 180.0 );
    float slownessStep = (vars->maxSlowness - vars->minSlowness) / (float)(vars->numSlownesses-1);
    // Plane wave, search for azimuth and slowness
    for( int iazim = 0; iazim < vars->numAzimuths; iazim++ ) {
      float azim_rad = (float)iazim * azimStep;
      //      fprintf(stderr,"Test azimuth: %f\n", azim_rad * 180.0 / M_PI);
      float sin_azim = sin(azim_rad);
      float cos_azim = cos(azim_rad);
      for( int islow = 0; islow < vars->numSlownesses; islow++ ) {
        float slowness_spkm = (float)islow * slownessStep;  // [s/km]
        float* outSamples = traceGather->trace(nTraces+beamCounter)->getTraceSamples();
        csTraceHeader* trcHdr = traceGather->trace(nTraces+beamCounter)->getTraceHeader();
        trcHdr->setFloatValue(vars->hdrId_azim, (float)(azim_rad*180.0/M_PI));
        trcHdr->setFloatValue(vars->hdrId_slowness, slowness_spkm);
        //  memcpy(outSamples,traceGather->trace(0)->getTraceSamples(),nSamples*sizeof(float));
        for( int itrc = 0; itrc < nTraces; itrc++ ) {
          float* inSamples = traceGather->trace(itrc)->getTraceSamples();
          double projection = rec_dx[itrc] * sin_azim + rec_dy[itrc] * cos_azim;
          float delayTime_ms = 0.0f;
          delayTime_ms = (float)(projection * slowness_spkm);  // [ms]
          int nearestSamp = (int)round( delayTime_ms / shdr->sampleInt );
          int minSamp = min( max(vars->startSamp,-nearestSamp), vars->endSamp );
          int maxSamp = max( min(vars->endSamp,vars->endSamp-nearestSamp), vars->startSamp );
          for( int isamp = minSamp; isamp <= maxSamp; isamp++ ) {
            outSamples[isamp-sampShift] += inSamples[isamp+nearestSamp];
          }
        }
        for( int isamp = vars->startSamp; isamp <= vars->endSamp; isamp++ ) {
          outSamples[isamp-sampShift] /= (float)nTraces;
        }
        beamCounter += 1;
      }
    }
  }  // END plane wave fan
  //----------------------------------------------------------------------------------
  else if( vars->method == mod_beam_forming::METHOD_PLANE_WAVE_GRID ) {
    float slownessStep = (vars->maxSlowness - vars->minSlowness) / (float)(vars->numSlownesses-1);
    // Plane wave, search for slowness XY
    for( int islowx = 0; islowx < vars->numSlownesses; islowx++ ) {
      float slowness_x = (float)islowx * slownessStep;
//      fprintf(stderr,"Test slowness X: %f\n", slowness_x);
      for( int islowy = 0; islowy < vars->numSlownesses; islowy++ ) {
        float slowness_y = (float)islowy * slownessStep;
        float slowness_spkm = sqrt( slowness_x*slowness_x + slowness_y*slowness_y ); // [s/km]
        float azim_rad = atan2( slowness_x, slowness_y );
        float sin_azim = sin(azim_rad);
        float cos_azim = cos(azim_rad);
        float* outSamples = traceGather->trace(nTraces+beamCounter)->getTraceSamples();
        csTraceHeader* trcHdr = traceGather->trace(nTraces+beamCounter)->getTraceHeader();
        trcHdr->setFloatValue(vars->hdrId_azim, (float)(azim_rad*180.0/M_PI));
        trcHdr->setFloatValue(vars->hdrId_slowness, slowness_spkm);
        trcHdr->setFloatValue(vars->hdrId_slowness_x, slowness_x);
        trcHdr->setFloatValue(vars->hdrId_slowness_y, slowness_y);
        //  memcpy(outSamples,traceGather->trace(0)->getTraceSamples(),nSamples*sizeof(float));
        for( int itrc = 0; itrc < nTraces; itrc++ ) {
          float* inSamples = traceGather->trace(itrc)->getTraceSamples();
          double projection = rec_dx[itrc] * sin_azim + rec_dy[itrc] * cos_azim;
          float delayTime_ms = 0.0;
          delayTime_ms = (float)(projection * slowness_spkm); // [ms]
          int nearestSamp = (int)round( delayTime_ms / shdr->sampleInt );
          int minSamp = min( max(vars->startSamp,-nearestSamp), vars->endSamp );
          int maxSamp = max( min(vars->endSamp,vars->endSamp-nearestSamp), vars->startSamp );
          for( int isamp = minSamp; isamp <= maxSamp; isamp++ ) {
            outSamples[isamp-sampShift] += inSamples[isamp+nearestSamp];
          }
        }
        for( int isamp = vars->startSamp; isamp <= vars->endSamp; isamp++ ) {
          outSamples[isamp-sampShift] /= (float)nTraces;
        }
        beamCounter += 1;
      }
    }
  } // END plane wave grid
  //----------------------------------------------------------------------------------
  else if( vars->method == mod_beam_forming::METHOD_POINT_SOURCE ) {
    // Plane wave, search for slowness XY
    for( int ix = 0; ix < vars->numStepsX; ix++ ) {
      float xp = (float)((double)ix * vars->incXY + vars->minX);
//      fprintf(stderr,"Test X coordinate: %f\n", xp);
      for( int iy = 0; iy < vars->numStepsY; iy++ ) {
        float yp = (float)((double)iy * vars->incXY + vars->minY);
        float* outSamples = traceGather->trace(nTraces+beamCounter)->getTraceSamples();
        csTraceHeader* trcHdr = traceGather->trace(nTraces+beamCounter)->getTraceHeader();
        trcHdr->setFloatValue(vars->hdrId_slowness, vars->ps_slowness);
        trcHdr->setFloatValue(vars->hdrId_beam_sou_x, xp);
        trcHdr->setFloatValue(vars->hdrId_beam_sou_y, yp);
        //  memcpy(outSamples,traceGather->trace(0)->getTraceSamples(),nSamples*sizeof(float));
        int numStackedTraces = 0;
        float numStackedGain = 0.0;
        for( int itrc = 0; itrc < nTraces; itrc++ ) {
          float* inSamples = traceGather->trace(itrc)->getTraceSamples();
          double dx = rec_dx[itrc]+rec_dx[0] - xp;
          double dy = rec_dy[itrc]+rec_dy[0] - yp;
          double offset = sqrt( dx*dx + dy*dy );
          if( vars->isOffsetSet ) {
            if( offset < vars->ps_minOffset || offset > vars->ps_maxOffset ) {
              continue;
            }
          }
          numStackedTraces += 1;
          float gain = 1.0f;
          if( vars->ps_gain > 0.0 ) {
            gain = (float)pow( max( 1.0, offset), (double)vars->ps_gain );
          }
          numStackedGain += gain;
          float delayTime_ms = (float)( offset * vars->ps_slowness ); // [ms]
          int nearestSamp = (int)round( delayTime_ms / shdr->sampleInt );
          int minSamp = min( max(vars->startSamp,-nearestSamp), vars->endSamp );
          int maxSamp = max( min(vars->endSamp,vars->endSamp-nearestSamp), vars->startSamp );
          for( int isamp = minSamp; isamp <= maxSamp; isamp++ ) {
            outSamples[isamp-sampShift] += gain*inSamples[isamp+nearestSamp];
          }
        }
        if( numStackedTraces > 0 ) {
          if( vars->ps_gain == 0.0 ) {
            for( int isamp = vars->startSamp; isamp <= vars->endSamp; isamp++ ) {
              outSamples[isamp-sampShift] /= (float)numStackedTraces;
            }
          }
          else {
            for( int isamp = vars->startSamp; isamp <= vars->endSamp; isamp++ ) {
              outSamples[isamp-sampShift] /= numStackedGain;
            }
          }
        }
        beamCounter += 1;
      }
    }
  } // END point source

  traceGather->freeTraces( 0, nTraces );
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_beam_forming_( csParamDef* pdef ) {
  pdef->setModule( "BEAM_FORMING", "Generate beams from receiver array" );

//  pdef->addParam( "method", "Method of beam forming", NUM_VALUES_FIXED);
//  pdef->addValue( "plane_wave", VALTYPE_OPTION );
//  pdef->addOption( "plane_wave", "Assume plane waves.", "Specify slowness range in parameter 'slowness'" );
//  pdef->addOption( "point_source", "Assume point source.", "Specify area in parameter 'area'" );

  pdef->addParam( "plane_wave", "Beam-form plane waves", NUM_VALUES_FIXED,
    "The output is a grid of test beams generated with parameters   min/max slowness(X)  x   min/max slowness(Y)" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Minimum slowness X/Y component [s/km]" );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Maximum slowness X/Y component [s/km]" );
  pdef->addValue( "50", VALTYPE_NUMBER, "Number of slowness steps in each dimension",
    "Number of output traces is  (N*2+1) * (N*2+1)" );

  pdef->addParam( "plane_wave_fan", "Beam-form plane waves in a fan", NUM_VALUES_FIXED,
    "The output is a grid of test beams generated with parameters   min/max slowness  x  min/max azimuth" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Minimum slowness [s/km]" );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Maximum slowness [s/km]" );
  pdef->addValue( "50", VALTYPE_NUMBER, "Number of slowness steps" );
  pdef->addValue( "50", VALTYPE_NUMBER, "Number of azimuth steps" );

  pdef->addParam( "point_source", "Beam-form point sources", NUM_VALUES_FIXED,
    "The output is a areal grid of test beams for the given area with parameter slowness" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Minimum X coordinate [m]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Minimum Y coordinate [m]" );
  pdef->addValue( "10.0", VALTYPE_NUMBER, "Maximum X coordinate [m]" );
  pdef->addValue( "10.0", VALTYPE_NUMBER, "Maximum Y coordinate [m]" );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Increment in each dimension [m]" );

  pdef->addParam( "point_source_slowness", "Constant slowness to use for point source beam-forming", NUM_VALUES_FIXED );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Slowness [s/km]" );

  pdef->addParam( "point_source_offset", "Minimum/maximum offset to use in beam forming", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Minimum offset [m]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Maximum offset [m]" );

  pdef->addParam( "point_source_gain", "Apply gain correction", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Apply gain r^g, where r is the source-receiver offset, and g the specified gain" );

  pdef->addParam( "window", "Analysis window", NUM_VALUES_VARIABLE );
  pdef->addValue( "0", VALTYPE_NUMBER, "Minimum time to analyse [ms]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Maximum time to analyse [ms]" );
  pdef->addValue( "zero", VALTYPE_OPTION, "Output trace option" );
  pdef->addOption( "zero", "Sample values outside of the analysis window are set to 0" );
  pdef->addOption( "cut", "Only samples in analysis window are output to output trace" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_beam_forming_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_beam_forming::VariableStruct* vars = reinterpret_cast<mod_beam_forming::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_beam_forming_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_beam_forming::VariableStruct* vars = reinterpret_cast<mod_beam_forming::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_beam_forming_( csParamDef* pdef ) {
  params_mod_beam_forming_( pdef );
}
extern "C" void _init_mod_beam_forming_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_beam_forming_( param, env, writer );
}
extern "C" bool _start_exec_mod_beam_forming_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_beam_forming_( env, writer );
}
extern "C" void _exec_mod_beam_forming_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_beam_forming_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_beam_forming_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_beam_forming_( env, writer );
}

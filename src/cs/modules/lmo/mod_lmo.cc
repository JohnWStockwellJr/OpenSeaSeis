/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csInterpolation.h"
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_lmo {
  struct VariableStruct {
    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_sou_z;
    int hdrId_rec_x;
    int hdrId_rec_y;
    int hdrId_rec_z;
    int hdrId_lmo;
    int hdrId_velocity;
    int mode;
    int dim;
    float bulk_shift;
    float velocity;
    cseis_geolib::csInterpolation* interpol;
    float* buffer;
  };

  static int const MODE_APPLY   = 1;
  static int const MODE_COMPUTE = 2;
  static int const MODE_REMOVE  = 3;

  static int const DIM_2D   = 11;
  static int const DIM_3D   = 12;
}
using mod_lmo::VariableStruct;

//*************************************************************************************************
// Init phase
//
void init_mod_lmo_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  vars->hdrId_lmo   = -1;
  vars->hdrId_sou_x = -1;
  vars->hdrId_sou_y = -1;
  vars->hdrId_sou_z = -1;
  vars->hdrId_sou_x = -1;
  vars->hdrId_rec_x = -1;
  vars->hdrId_rec_y = -1;
  vars->hdrId_rec_z = -1;

  vars->mode = mod_lmo::MODE_APPLY;
  vars->dim  = mod_lmo::DIM_3D;
  vars->bulk_shift = 0;
  vars->velocity   = 1500;
  vars->interpol = NULL;
  vars->buffer = NULL;
  
  std::string text;
  if( param->exists("dimension") ) {
    param->getString("dimension", &text);
    if( !text.compare("2d") ) {
      vars->dim = mod_lmo::DIM_2D;
    }
    else if( !text.compare("3d") ) {
      vars->dim = mod_lmo::DIM_3D;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
      env->addError();
    }
  }
  if( param->exists("mode") ) {
    param->getString("mode", &text);
    if( !text.compare("apply") ) {
      vars->mode = mod_lmo::MODE_APPLY;
    }
    else if( !text.compare("compute") ) {
      vars->mode = mod_lmo::MODE_COMPUTE;
    }
    else if( !text.compare("remove") ) {
      vars->mode = mod_lmo::MODE_REMOVE;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
      env->addError();
    }
  }

  std::string hdrname_lmo("stat_lmo");
  std::string hdrname_sou_x("sou_x");
  std::string hdrname_sou_y("sou_y");
  std::string hdrname_sou_z("sou_z");
  std::string hdrname_rec_x("rec_x");
  std::string hdrname_rec_y("rec_y");
  std::string hdrname_rec_z("rec_z");

  if( param->exists("hdr_lmo")  ) {
    param->getString("hdr_lmo", &hdrname_lmo);
  }
  if( param->exists("hdr_sou")  ) {
    param->getString("hdr_sou", &hdrname_sou_x, 0);
    param->getString("hdr_sou", &hdrname_sou_y, 1);
    if( vars->dim == mod_lmo::DIM_3D ) param->getString("hdr_sou", &hdrname_sou_z, 2);
  }
  if( param->exists("hdr_rec")  ) {
    param->getString("hdr_rec", &hdrname_rec_x, 0);
    param->getString("hdr_rec", &hdrname_rec_y, 1);
    if( vars->dim == mod_lmo::DIM_3D ) param->getString("hdr_rec", &hdrname_rec_z, 2);
  }
  
  param->getString( "velocity", &text );
  csFlexNumber number;
  if( !number.convertToNumber( text ) ) {
    vars->hdrId_velocity = hdef->headerIndex(text);
  }
  else {
    vars->velocity = number.floatValue();
  }
  
  if( param->exists("bulk_shift")  ) {
    param->getFloat( "bulk_shift", &vars->bulk_shift );
  }

  vars->hdrId_sou_x = hdef->headerIndex(hdrname_sou_x);
  vars->hdrId_sou_y = hdef->headerIndex(hdrname_sou_y);
  vars->hdrId_rec_x = hdef->headerIndex(hdrname_rec_x);
  vars->hdrId_rec_y = hdef->headerIndex(hdrname_rec_y);
  if( vars->dim == mod_lmo::DIM_3D ) {
    vars->hdrId_sou_z = hdef->headerIndex(hdrname_sou_z);
    vars->hdrId_rec_z = hdef->headerIndex(hdrname_rec_z);
  }

  if( vars->mode != mod_lmo::MODE_COMPUTE ) {
    int numCoefficients = 8;
    vars->interpol = new csInterpolation( shdr->numSamples, shdr->sampleInt, numCoefficients );
    vars->buffer = new float[shdr->numSamples];
  }
  if( vars->mode != mod_lmo::MODE_REMOVE && !hdef->headerExists(hdrname_lmo) ) {
    hdef->addHeader(cseis_geolib::TYPE_FLOAT,hdrname_lmo,"LMO correction [ms]");
  }
  vars->hdrId_lmo = hdef->headerIndex(hdrname_lmo);
}

//*************************************************************************************************
// Exec phase
//
void exec_mod_lmo_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const*    shdr = env->superHeader;
  
  int numTraces = traceGather->numTraces();
  for( int itrc = 0; itrc < numTraces; itrc++ ) {
    csTrace* trace = traceGather->trace( itrc );
    csTraceHeader* trcHdr = trace->getTraceHeader();
    double sou_x = trcHdr->doubleValue( vars->hdrId_sou_x );
    double sou_y = trcHdr->doubleValue( vars->hdrId_sou_y );
    double rec_x = trcHdr->doubleValue( vars->hdrId_rec_x );
    double rec_y = trcHdr->doubleValue( vars->hdrId_rec_y );
    double distance;
    if( vars->dim == mod_lmo::DIM_2D ) {
      distance = sqrt( pow(sou_x-rec_x,2) + pow(sou_y-rec_y,2) );
    }
    else { //if( 
      double sou_z = trcHdr->floatValue( vars->hdrId_sou_z );
      double rec_z = trcHdr->floatValue( vars->hdrId_rec_z );
      distance = sqrt( pow(sou_x-rec_x,2) + pow(sou_y-rec_y,2) + pow(sou_z-rec_z,2) );
    }
    if( vars->hdrId_velocity >= 0 ) vars->velocity = trcHdr->floatValue( vars->hdrId_velocity );

    float statLMO   = 1000.0f * distance / vars->velocity; // 1000x to convert to [ms]
    trcHdr->setFloatValue( vars->hdrId_lmo, statLMO );

    if( vars->mode != mod_lmo::MODE_COMPUTE ) {
      float statShift_ms = -statLMO + vars->bulk_shift;
      if( vars->mode == mod_lmo::MODE_REMOVE ) {
        statShift_ms = statLMO - vars->bulk_shift;
      }
      // Apply static shift to data
      float* samples = trace->getTraceSamples();
      memcpy( vars->buffer, samples, shdr->numSamples*sizeof(float) );
      vars->interpol->static_shift( statShift_ms, vars->buffer, samples );
    }
  }
  
  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_lmo_( csParamDef* pdef ) {
  pdef->setModule( "LMO", "Linear moveout correction" );

  pdef->addParam( "dimension", "Dimension of LMO application", NUM_VALUES_FIXED );
  pdef->addValue( "3d", VALTYPE_OPTION );
  pdef->addOption( "3d", "Apply 3D LMO, use XYZ coordinates" );
  pdef->addOption( "2d", "Apply 2D LMO, use XY coordinates" );
  
  pdef->addParam( "mode", "Mode of application", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "compute", "(Re-)compute LMO correction term" );
  pdef->addOption( "apply", "(Re-)compute LMO correction term & apply LMO correction" );
  pdef->addOption( "remove", "Remove LMO correction stored in trace header specified in user parameter 'hdr_lmo', honoring the provided 'bulk_shift'" );

  pdef->addParam( "velocity", "Velocity to usefor LMO [m/s]", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER );

  pdef->addParam( "hdr_lmo", "Trace header name where computed & applied", NUM_VALUES_FIXED );
  pdef->addValue( "stat_lmo", VALTYPE_STRING );

  pdef->addParam( "hdr_sou", "Trace header names containing source XYZ coordinates", NUM_VALUES_FIXED );
  pdef->addValue( "sou_x", VALTYPE_STRING );
  pdef->addValue( "sou_y", VALTYPE_STRING );
  pdef->addValue( "sou_z", VALTYPE_STRING );

  pdef->addParam( "hdr_rec", "Trace header names containing receiver XYZ coordinates", NUM_VALUES_FIXED );
  pdef->addValue( "rec_x", VALTYPE_STRING );
  pdef->addValue( "rec_y", VALTYPE_STRING );
  pdef->addValue( "rec_z", VALTYPE_STRING );

  pdef->addParam( "bulk_shift", "Apply static bulk shift to all traces", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Static shift [ms]. Positive value shifts samples downwards" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_lmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_lmo::VariableStruct* vars = reinterpret_cast<mod_lmo::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_lmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_lmo::VariableStruct* vars = reinterpret_cast<mod_lmo::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->interpol != NULL ) {
    delete vars->interpol;
    vars->interpol = NULL;
  }
  if( vars->buffer != NULL ) {
    delete [] vars->buffer;
    vars->buffer = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_lmo_( csParamDef* pdef ) {
  params_mod_lmo_( pdef );
}
extern "C" void _init_mod_lmo_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_lmo_( param, env, writer );
}
extern "C" bool _start_exec_mod_lmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_lmo_( env, writer );
}
extern "C" void _exec_mod_lmo_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_lmo_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_lmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_lmo_( env, writer );
}

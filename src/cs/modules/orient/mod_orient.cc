/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_methods.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ORIENT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_orient {
  struct VariableStruct {
    int hdrId_rec_x;
    int hdrId_rec_y;
    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_rec_z;
    int hdrId_sou_z;
    int hdrId_vec_x;
    int hdrId_vec_y;
    int hdrId_vec_z;
    int hdrId_sensor;

    int hdrId_azim_in;
    int hdrId_azim_out;
    int hdrId_roll;
    int hdrId_tilt;
    int hdrId_tiltx;
    int hdrId_tilty;

    int input;
    float angle_inline_dip;
    float angle_xline_dip;
    float v1v2_ratio;
    bool isOBCMode;
    int nTraces;
    int output_option;
    bool compute_azim;

    bool compute_tilt_xy;
    bool compute_tilt_roll;
  };

  static const double EPSILON = 0.000001;

  static const int INPUT_XYZ  = 20;
  static const int INPUT_XYZP = 21;
  static const int INPUT_ENSEMBLE = 22;

  static const int SENSOR_INDEX[] = { 3, 4, 5, 1 };
  static const int ID_X = 0;
  static const int ID_Y = 1;
  static const int ID_Z = 2;
  static const int ID_P = 3;

  static const int SENSOR_INDEX_Z2 = 2;
  static const int SENSOR_INDEX_Z5 = 5;

  static int const OUTPUT_ALL     = 1;
  static int const OUTPUT_FIRST   = 2;
  static int const OUTPUT_LAST    = 3;

  void checkHeader( std::string headerName, csTraceHeaderDef* hdef, csLogWriter* writer ) {
    if( !hdef->headerExists(headerName) ){
      writer->error("Trace header '%s' not found.", headerName.c_str());
    }
    int type = hdef->headerType( headerName );
    if( type != TYPE_FLOAT &&  type != TYPE_DOUBLE ) {
      writer->error("Trace header '%s' has the wrong type. Should be FLOAT or DOUBLE", headerName.c_str());
    }
  }
}
using namespace mod_orient;

/*
void orientation_compute_roll_tilt (
  float const* vec_polarisation, // Polarisation vector in as-laid coordinate system
  double rcvx,     // Receiver X position [m]
  double rcvy,     // Receiver Y position [m]
  double rcvz,     // Receiver Z position [m]
  double srcx,     // Source X position [m]
  double srcy,     // Source Y position [m]
  double srcz,     // Source Z position [m]
  float v1v2_ratio,         // Velocity ratio at seabed
  int   isOBCMode,          // true(1) if source is above receiver, false(0) if source is beneath receiver
  float angle_azim,         // Receiver azimuth [rad]
  float* angle_tilt,        // Tilt angle [rad]
  float* angle_roll );       // Roll angle [rad]

void orientation_compute_azimuth(
  float const* vec_polarisation,
  double rcvx,     // Receiver X position [m]
  double rcvy,     // Receiver Y position [m]
  double srcx,     // Source X position [m]
  double srcy,     // Source Y position [m]
  float* angle_azim // Computed azimuth [rad]
);
*/

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_orient_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );



  vars->hdrId_rec_x   = -1;
  vars->hdrId_rec_y   = -1;
  vars->hdrId_sou_x   = -1;
  vars->hdrId_sou_y   = -1;
  vars->hdrId_rec_z = -1;
  vars->hdrId_sou_z   = -1;
  vars->hdrId_vec_x   = -1;
  vars->hdrId_vec_y   = -1;
  vars->hdrId_vec_z   = -1;
  vars->hdrId_azim_in = -1;
  vars->hdrId_azim_out= -1;
  vars->hdrId_sensor  = -1;
  vars->hdrId_roll    = -1;
  vars->hdrId_tilt    = -1;
  vars->hdrId_tiltx   = -1;
  vars->hdrId_tilty   = -1;
  vars->input      = 0;
  vars->angle_inline_dip = 0;
  vars->angle_xline_dip  = 0;
  vars->v1v2_ratio       = 1.0;
  vars->isOBCMode        = true;
  vars->nTraces          = 0;
  vars->output_option    = OUTPUT_ALL;
  vars->compute_azim     = false;
  vars->compute_tilt_xy  = false;
  vars->compute_tilt_roll= false;

  //-----------------------------------
  std::string text;
  if( param->exists("input") ) {
    param->getString("input",&text);
    if( !text.compare("xyz") ) {
      vars->input = INPUT_XYZ;
      edef->setTraceSelectionMode( TRCMODE_FIXED, 3 );
    }
    else if( !text.compare("xyzp") ) {
      vars->input = INPUT_XYZP;
      edef->setTraceSelectionMode( TRCMODE_FIXED, 4 );
    }
    else if( !text.compare("ensemble") ) {
      vars->input = INPUT_ENSEMBLE;
      edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  else {  // Default setting
    vars->input = INPUT_XYZ;
    edef->setTraceSelectionMode( TRCMODE_FIXED, 3 );
  }

  //-------------------------------------
  if( param->exists("v1v2_ratio") ) {
    param->getFloat( "v1v2_ratio", &vars->v1v2_ratio );
  }
  else {
    vars->v1v2_ratio = 1.0;
  }

  //-------------------------------------
  param->getString( "output", &text );
  if( !text.compare("all") ) {
    vars->output_option = OUTPUT_ALL;
  }
  else if( !text.compare("last") ) {
    vars->output_option = OUTPUT_LAST;
  }
  else {
    writer->error("Output option not recognized: %s.", text.c_str());
  }

  //-------------------------------------
  if( param->exists( "inline_dip" ) ) {
    param->getFloat( "inline_dip", &vars->angle_inline_dip );
    vars->angle_inline_dip = vars->angle_inline_dip * M_PI / 180.0;
  }
  else {
    vars->angle_inline_dip = 0.0;
  }
  if( param->exists( "xline_dip" ) ) {
    param->getFloat( "xline_dip", &vars->angle_xline_dip );
    vars->angle_xline_dip = vars->angle_xline_dip * M_PI / 180.0;
  }
  else {
    vars->angle_xline_dip = 0.0;
  }

  // Set default options
  vars->isOBCMode = true;
  if( param->exists("advanced") ) {
    param->getString( "advanced", &text );
    if( !text.compare("yes") ) {
      param->getString( "acquisition_mode", &text );
      if( !text.compare("obc") ) {
        vars->isOBCMode = true;
      }
      else if( !text.compare("land") ) {
        vars->isOBCMode = false;
      }
      else {
        writer->error("Unknown option: '%s'", text.c_str());
      }
    }
    else if( !text.compare("no") ) {
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  //----------------------------------------------------------------------
  //
  if( param->exists("compute_azim") ) {
    param->getString( "compute_azim", &text );
    if( !text.compare("yes") ) {
      vars->compute_azim = true;
      int numValues = param->getNumValues( "compute_azim" );
      if( numValues > 2 ) writer->error("Parameter 'azim' expects one or two user supplied values. Found %d", numValues );
      if( numValues > 1 ) {
        param->getString( "compute_azim", &text, 1 );
        vars->hdrId_azim_out = hdef->headerIndex( text );
      }
      else {
        // Default trace header
        vars->hdrId_azim_out = hdef->headerIndex( "an_azim" );
      }
    }
    else if( !text.compare("no") ) {
      vars->compute_azim = false;
    }
    else {
      writer->line("Option not recognized: '%s'", text.c_str());
    }
  }

  if( param->exists("compute_tilt_xy") ) {
    param->getString( "compute_tilt_xy", &text );
    if( !text.compare("yes") ) {
      vars->compute_tilt_xy = true;
      int numValues = param->getNumValues( "compute_tilt_xy" );
      if( numValues > 3 ) writer->error("Parameter 'compute_tilt_xy' expects up to three user supplied values. Found %d", numValues );
      if( numValues > 2 ) {
        param->getString( "compute_tilt_xy", &text, 2 );
        vars->hdrId_tilty = hdef->headerIndex( text );
      }
      else { // Default trace header
        vars->hdrId_tilty = hdef->headerIndex( "an_tilty" );
      }
      if( numValues > 1 ) {
        param->getString( "compute_tilt_xy", &text, 1 );
        vars->hdrId_tiltx = hdef->headerIndex( text );
      }
      else { // Default trace header
        vars->hdrId_tiltx = hdef->headerIndex( "an_tiltx" );
      }
    }
    else if( !text.compare("no") ) {
      vars->compute_tilt_xy = false;
    }
    else {
      writer->line("Option not recognized: '%s'", text.c_str());
    }
  }

  if( param->exists("compute_tilt_roll") ) {
    param->getString( "compute_tilt_roll", &text );
    if( !text.compare("yes") ) {
      vars->compute_tilt_roll = true;
      int numValues = param->getNumValues( "compute_tilt_roll" );
      if( numValues > 3 ) writer->error("Parameter 'compute_tilt_roll' expects up to three user supplied values. Found %d", numValues );
      if( numValues > 2 ) {
        param->getString( "compute_tilt_roll", &text, 2 );
        vars->hdrId_roll = hdef->headerIndex( text );
      }
      else { // Default trace header
        vars->hdrId_roll = hdef->headerIndex( "an_roll" );
      }
      if( numValues > 1 ) {
        param->getString( "compute_tilt_roll", &text, 1 );
        vars->hdrId_tilt = hdef->headerIndex( text );
      }
      else { // Default trace header
        vars->hdrId_tilt = hdef->headerIndex( "an_tilt" );
      }
    }
    else if( !text.compare("no") ) {
      vars->compute_tilt_roll = false;
    }
    else {
      writer->line("Option not recognized: '%s'", text.c_str());
    }
  }

  if( !vars->compute_tilt_roll && !vars->compute_tilt_xy && !vars->compute_azim ) {
    writer->error("No angles selected for computation: Select either parameter 'tilt_roll', 'tilt_xy' or 'azim', or any combination of these");
  }
  if( vars->compute_tilt_xy || vars->compute_tilt_roll ) {
    if( param->exists("hdr_azim") ) {
      param->getString("hdr_azim",&text);
      vars->hdrId_azim_in = hdef->headerIndex(text);
    }
    else {
      //      vars->hdrId_azim_in = hdef->headerIndex("an_azim");
      writer->error("Parameter 'hdr_azim' not specified, required for the chosen options. This gives the input trace header name containing the sensor azimuth.");
    }
  }

  //----------------------------------------------------------------------
  //

  if( !hdef->headerExists("sensor") ){
    writer->error("Trace header 'sensor' not found.");
  }

  vars->hdrId_sou_x    = hdef->headerIndex( "sou_x" );
  vars->hdrId_sou_y    = hdef->headerIndex( "sou_y" );
  vars->hdrId_rec_x    = hdef->headerIndex( "rec_x" );
  vars->hdrId_rec_y    = hdef->headerIndex( "rec_y" );
  vars->hdrId_sou_z    = hdef->headerIndex( "sou_z" );
  vars->hdrId_rec_z = hdef->headerIndex( "rec_z" );
  vars->hdrId_vec_x    = hdef->headerIndex( "vec_x" );
  vars->hdrId_vec_y    = hdef->headerIndex( "vec_y" );
  vars->hdrId_vec_z    = hdef->headerIndex( "vec_z" );
  vars->hdrId_sensor   = hdef->headerIndex( "sensor" );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_orient_(
                      csTraceGather* traceGather,
                      int* port,
                      int* numTrcToKeep,
                      csExecPhaseEnv* env,
                      csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
//  csTraceHeaderDef const* hdef = env->headerDef;


  //---------------------------------------------
  int nTraces = traceGather->numTraces();

  if( vars->input == INPUT_XYZP ) {
    if( traceGather->numTraces() != 4 ) {
      writer->error("Incorrect number of traces in input gather. Expected: 4 (XYZP), found: %d", nTraces );
    }
  }
  else if( vars->input == INPUT_XYZ ) {
    if( traceGather->numTraces() != 3 ) {
      writer->error("Incorrect number of traces in input gather. Expected: 3 (XYZ), found: %d", nTraces );
    }
  }
  int no_value = -999;
  int trace_index[] = {no_value,no_value,no_value};
  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    int sensor = traceGather->trace(itrc)->getTraceHeader()->intValue(vars->hdrId_sensor);
    if( sensor == SENSOR_INDEX_Z2 ) sensor = SENSOR_INDEX_Z5;  // Convert sensor index 2 to 5 (both are Z component)
    if( sensor >= SENSOR_INDEX[ID_X] && sensor <= SENSOR_INDEX[ID_Z] ) {
      if( trace_index[sensor-3] == no_value ) {
        trace_index[sensor-3] = itrc;
      }
      else {
        writer->error("Input gather contains more than one trace for sensor %d.", sensor);
        return;
      }
    }
  }

  //
  // Read in sensor numbers for input traces
  // Check that all sensors 3, 4 and 5 are present
  //

  csTraceHeader* trcHdrSensor[3];

  for( int idSensor = 0; idSensor < 3; idSensor++ ) {
    if( trace_index[idSensor] == no_value ) {
      writer->error("Input gather is missing a sensor %d trace.", SENSOR_INDEX[idSensor]);
    }
    else {
      trcHdrSensor[idSensor] = traceGather->trace(trace_index[idSensor])->getTraceHeader();
    }
  }

  csTraceHeader* trcHdr = trcHdrSensor[0];
  
  double srcx = trcHdr->doubleValue( vars->hdrId_sou_x );
  double srcy = trcHdr->doubleValue( vars->hdrId_sou_y );
  float  srcz = trcHdr->floatValue( vars->hdrId_sou_z );

  double rcvx = trcHdr->doubleValue( vars->hdrId_rec_x );
  double rcvy = trcHdr->doubleValue( vars->hdrId_rec_y );
  float  rcvz = trcHdr->floatValue( vars->hdrId_rec_z );

  float vec_polarisation[3];
  vec_polarisation[0] = trcHdr->floatValue( vars->hdrId_vec_x );
  vec_polarisation[1] = trcHdr->floatValue( vars->hdrId_vec_y );
  vec_polarisation[2] = trcHdr->floatValue( vars->hdrId_vec_z );

  //--------------------------------------------------------------------
  //
  
  if( vars->compute_tilt_roll || vars->compute_tilt_xy ) {
    float angle_azim = DEG2RAD( trcHdr->floatValue( vars->hdrId_azim_in ) );
    float angle_tilt = 0;
    float angle_roll = 0;

    // Derive recorded unit vector. Make sure vector is normalised
    float norm = sqrt( CS_SQR(vec_polarisation[0]) + CS_SQR(vec_polarisation[1]) + CS_SQR(vec_polarisation[2]));
    if( norm != 1.0 ) {
      if( norm < EPSILON ) {
        norm = EPSILON;
      }
      for( int i = 0; i < 3; i++ ) {
        vec_polarisation[i] /= norm;
      }
    }
    int isOBCMode;
    if( vars->isOBCMode ) {
      isOBCMode = 1;
    }
    else {
      isOBCMode = 0;
    }

    orientation_compute_roll_tilt (
      vec_polarisation, // Polarisation vector in as-laid coordinate system
      rcvx,     // Receiver X position [m]
      rcvy,     // Receiver Y position [m]
      rcvz,     // Receiver Z position [m]
      srcx,     // Source X position [m]
      srcy,     // Source Y position [m]
      srcz,     // Source Z position [m]
      vars->v1v2_ratio,         // Velocity ratio at seabed
      isOBCMode,          // true(1) if source is above receiver, false(0) if source is beneath receiver
      angle_azim,         // Receiver azimuth [rad]
      &angle_tilt,        // Tilt angle [rad]
      &angle_roll );       // Roll angle [rad]

    angle_roll -= vars->angle_xline_dip;
    angle_tilt -= vars->angle_inline_dip;

    angle_roll = RAD2DEG( angle_roll );
    angle_tilt = RAD2DEG( angle_tilt );

    if( vars->compute_tilt_roll ) {
      for( int i = 0; i < nTraces; i++ ) {
        trcHdr = traceGather->trace(i)->getTraceHeader();
        trcHdr->setFloatValue( vars->hdrId_roll, angle_roll );
        trcHdr->setFloatValue( vars->hdrId_tilt, angle_tilt );
      }
    }

    if( vars->compute_tilt_xy ) {
      angle_roll -= 90;  // Convert roll angle to tilty
      for( int i = 0; i < nTraces; i++ ) {
        trcHdr = traceGather->trace(i)->getTraceHeader();
        trcHdr->setFloatValue( vars->hdrId_tilty, angle_roll );
        trcHdr->setFloatValue( vars->hdrId_tiltx, angle_tilt );
      }
    }
  
  }

  if( vars->compute_azim ) {
    float angle_azim = 0;

    // Make sure 2D vector is normalised
    float norm = sqrt( CS_SQR(vec_polarisation[0]) + CS_SQR(vec_polarisation[1]) );
    if( norm < EPSILON ) {
      norm = EPSILON;
    }
    vec_polarisation[0] /= norm;  // Inline
    vec_polarisation[1] /= norm;  // Crossline

    orientation_compute_azimuth( vec_polarisation, rcvx, rcvy, srcx, srcy, &angle_azim );

    angle_azim = RAD2DEG( angle_azim );
    if( angle_azim < 0.0 ) angle_azim += 360.0;
    for( int i = 0; i < nTraces; i++ ) {
      trcHdr = traceGather->trace(i)->getTraceHeader();
      trcHdr->setFloatValue( vars->hdrId_azim_out, angle_azim );
    }
  }
//  else {
//    angle_tilt = 0.0;
//    angle_roll = M_PI_2;
//    angle_azim = 0.0;
//  } // END if
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_orient_( csParamDef* pdef ) {
  pdef->setModule( "ORIENT", "Solve sensor orientation", "Derive sensor orientation angles from input XYZ vector, computed by module $HODOGRAM" );

  pdef->addParam( "input", "Input data", NUM_VALUES_FIXED, "Specifies what type of input data is expected in this module" );
  pdef->addValue( "xyz", VALTYPE_OPTION );
  pdef->addOption( "xyz", "3 sequential input traces, containing sensors 3(X), 4(Y) and 5(Z). Sorting of these three may be arbitrary." );
  pdef->addOption( "xyzp", "4 sequential input traces, containing sensors 3(X), 4(Y) and 5(Z), and 4th trace with sensor 1(hydrophone). Sorting of these three may be arbitrary." );
  pdef->addOption( "ensemble", "Input whole ensemble (user sorted). Each input ensemble must contain one trace each for sensor 3(X), 4(Y) and 5(Z). Sorting of these three may be arbitrary." );

  pdef->addParam( "compute_tilt_roll", "Compute tilt & roll angles?", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not compute tilt & roll angles." );
  pdef->addOption( "yes", "Compute roll & tilt angles." );
  pdef->addValue( "an_tilt", VALTYPE_STRING, "Output trace header name for tilt angle [deg]" );
  pdef->addValue( "an_roll", VALTYPE_STRING, "Output trace header name for roll angle [deg]" );
  
  pdef->addParam( "compute_tilt_xy", "Compute tilt X & tilt Y angles?", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not compute tilt X and tilt Y angles." );
  pdef->addOption( "yes", "Compute tilt X and tilt Y angles." );
  pdef->addValue( "an_tiltx", VALTYPE_STRING, "Output trace header name for tilt X angle [deg]" );
  pdef->addValue( "an_tilty", VALTYPE_STRING, "Output trace header name for tilt Y angle [deg]" );
  
  pdef->addParam( "compute_azim", "Compute azimuth angle?", NUM_VALUES_VARIABLE, "Usually, this is the positive X sensor direction" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not compute azimuth angle." );
  pdef->addOption( "yes", "Compute azimuth angle." );
  pdef->addValue( "an_azim", VALTYPE_STRING, "Output trace header name for azimuth angle [deg]" );

  pdef->addParam( "hdr_azim", "Trace header containing azimuth (= positive X direction)", NUM_VALUES_FIXED );
  pdef->addValue( "an_azim", VALTYPE_STRING, "Trace header name containing azimuth / sensor X direction [deg]" );

  pdef->addParam( "output", "Specify what traces shall be output by this module", NUM_VALUES_FIXED );
  pdef->addValue( "all", VALTYPE_OPTION );
  pdef->addOption( "all", "Output all traces that were input to this module." );
  pdef->addOption( "last", "Only output last trace." );

  pdef->addParam( "inline_dip", "Seabed dip in inline direction", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Seabed inline dip [deg]" );
  pdef->addParam( "xline_dip", "Seabed dip in xline direction", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Seabed xline dip [deg]" );

  pdef->addParam( "advanced", "Advanced options?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Use advanced options." );
  pdef->addOption( "no", "Do not use advanced options." );

  pdef->addParam( "acquisition_mode", "Acquisition mode", NUM_VALUES_FIXED );
  pdef->addValue( "obc", VALTYPE_OPTION );
  pdef->addOption( "obc", "OBC acquisition mode. Sources are above receivers." );
  pdef->addOption( "land", "Land acquisition mode. Sources are beneath receivers." );

  pdef->addParam( "v1v2_ratio", "V1/V2 ratio at seabed (v1: water vel, v2: seabed vel)", NUM_VALUES_FIXED );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "V1/V2 velocity ratio at seabed" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_orient_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_orient::VariableStruct* vars = reinterpret_cast<mod_orient::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_orient_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_orient::VariableStruct* vars = reinterpret_cast<mod_orient::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_orient_( csParamDef* pdef ) {
  params_mod_orient_( pdef );
}
extern "C" void _init_mod_orient_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_orient_( param, env, writer );
}
extern "C" bool _start_exec_mod_orient_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_orient_( env, writer );
}
extern "C" void _exec_mod_orient_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_orient_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_orient_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_orient_( env, writer );
}

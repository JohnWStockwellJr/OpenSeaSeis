/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csRotation.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ROTATE
 *
 * This module has only been -partially- converted from using C-style methods to using C++ class csRotation
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_rotate {
  struct VariableStruct {
    int nTraces;
    int mode;
    int modeEffective;
    int hdrId_an_from;
    int hdrId_an_to;
    int orient_xyz;
    int orient_z;
    int hdrId_azim;
    int hdrId_incl_i;
    int hdrId_incl_c;
    int hdrId_incl_v;
    int hdrId_roll;
    int hdrId_tilt;
    int hdrId_tiltx;
    int hdrId_tilty;
    int hdrId_recx;
    int hdrId_recy;
    int hdrId_soux;
    int hdrId_souy;
    int hdrId_sensor;
    int rot_method;
    int input;
    int numInputTraces;
    float max_valid_tilt;
    float tiltx_prev;
    float tilty_prev;
    bool log_tiltxy_warnings;
    cseis_geolib::csRotation* rotation;
  };

  static const int INCLINOMETER  = 31;
  static const int SEISMIC       = 32;

  static const int METHOD_2D_AZIM   = 41;
  static const int METHOD_TILT_ROLL = 42;
  static const int METHOD_ICV       = 43;
  static const int METHOD_GALPERIN  = 44;
  static const int METHOD_ARMSS     = 45;
  static const int METHOD_ARMSS_TILT_ROLL = 46;
  static const int METHOD_TRILOBIT  = 47;
  static const int METHOD_TILT_XY   = 48;
  static const int METHOD_2D_RADIAL = 49;
  static const int METHOD_FROM_TO   = 50;

  static const int ROLLTILT_2_ICV = 10;
  static const int ICV_2_ROLLTILT = 11;
  static const int INPUT_XY   = 20;
  static const int INPUT_XYZ  = 21;
  static const int INPUT_XYZP = 22;
  static const int INPUT_ENSEMBLE = 23;
  static const int SENSOR_INDEX[] = { 3, 4, 5, 1 };
  static const int ID_X = 0;
  static const int ID_Y = 1;
  static const int ID_Z = 2;
  static const int ID_P = 3;

  static const int ORIENT_RIGHTHAND = 71;
  static const int ORIENT_LEFTHAND  = 72;
  static const int ORIENT_Z_DOWN    = 73;
  static const int ORIENT_Z_UP      = 74;

  static const int SENSOR_INDEX_Z2 = 2;
  static const int SENSOR_INDEX_Z5 = 5;

  static const int MODE_APPLY  = 11;
  static const int MODE_REMOVE = 12;

  static const float NO_TILT = 999;

  void checkHeaderFloat( std::string headerName, csTraceHeaderDef* hdef, csLogWriter* writer ) {
    if( !hdef->headerExists(headerName) ){
      writer->error("Trace header '%s' not found.", headerName.c_str());
    }
    int type = hdef->headerType( headerName );
    if( type != TYPE_FLOAT &&  type != TYPE_DOUBLE ) {
      writer->error("Trace header '%s' has the wrong type. Should be FLOAT or DOUBLE", headerName.c_str());
    }
  }
  void checkHeader( std::string headerName, csTraceHeaderDef* hdef, csLogWriter* writer ) {
    if( !hdef->headerExists(headerName) ){
      writer->error("Trace header '%s' not found.", headerName.c_str());
    }
  }

}
using mod_rotate::VariableStruct;

void rotate_xyz_icv( float* xTrace,
                     float* yTrace,
                     float* zTrace,
                     float* weight_i,
                     float* weight_c,
                     float* weight_v,
                     int nSamples,
                     bool apply );

void rotate_xy_azim( float* xTrace,
                     float* yTrace,
                     float azim,
                     int nSamples,
                     bool apply );

void rotate_xy( float* xTrace,
                float* yTrace,
                float azim,
                int nSamples );

void rotate_xyz_roll_tilt( float* xTrace,
                           float* yTrace,
                           float* zTrace,
                           float roll,
                           float tilt,
                           int nSamples,
                           bool apply );

void rotate_xyz_tiltxy( float* xTrace,
                         float* yTrace,
                         float* zTrace,
                         float roll,
                         float tilt,
                         int nSamples,
                         bool apply );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_rotate_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );


  // Defaults
  vars->nTraces = 0;
  vars->mode    = mod_rotate::MODE_APPLY;
  vars->modeEffective = mod_rotate::MODE_APPLY;
  vars->hdrId_azim   = -1;
  vars->hdrId_incl_i = -1;
  vars->hdrId_incl_c = -1;
  vars->hdrId_incl_v = -1;
  vars->hdrId_roll   = -1;
  vars->hdrId_tilt   = -1;
  vars->hdrId_tiltx  = -1;
  vars->hdrId_tilty  = -1;
  vars->hdrId_sensor = -1;
  vars->hdrId_recx   = -1;
  vars->hdrId_recy   = -1;
  vars->hdrId_soux   = -1;
  vars->hdrId_an_from = -1;
  vars->hdrId_an_to   = -1;

  vars->orient_xyz = mod_rotate::ORIENT_RIGHTHAND;
  vars->orient_z   = mod_rotate::ORIENT_Z_DOWN;

  vars->rot_method   = 0;
  vars->input        = 0;
  vars->numInputTraces = 0;
  vars->rotation     = NULL;

  vars->tiltx_prev = mod_rotate::NO_TILT;
  vars->tilty_prev = mod_rotate::NO_TILT;
  vars->max_valid_tilt = mod_rotate::NO_TILT;
  //--------------------------------------------
  
  std::string text;
  if( param->exists("mode") ) {
    param->getString( "mode", &text );
    if( !text.compare("apply") ) {
      vars->mode = mod_rotate::MODE_APPLY;
    }
    else if( !text.compare("remove") ) {
      vars->mode = mod_rotate::MODE_REMOVE;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  if( param->exists("orient") ) {
    param->getString( "orient", &text, 0 );
    if( !text.compare("right") ) {
      vars->orient_xyz = mod_rotate::ORIENT_RIGHTHAND;
    }
    else if( !text.compare("left") ) {
      vars->orient_xyz = mod_rotate::ORIENT_LEFTHAND;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
    param->getString( "orient", &text, 1 );
    if( !text.compare("down") ) {
      vars->orient_z = mod_rotate::ORIENT_Z_DOWN;
    }
    else if( !text.compare("up") ) {
      vars->orient_z = mod_rotate::ORIENT_Z_UP;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  std::string hdrName_roll("an_roll");
  std::string hdrName_tilt("an_tilt");
  std::string hdrName_azim("an_azim");
  std::string hdrName_tiltx("an_tiltx");
  std::string hdrName_tilty("an_tilty");
  std::string hdrName_from("none");
  std::string hdrName_to("none");

  param->getString( "method", &text );
  if( !text.compare("azim") ) {
    vars->rot_method = mod_rotate::METHOD_2D_AZIM;
    if( param->exists("hdr_azim") ) {
      param->getString("hdr_azim",&hdrName_azim);
    }
  }
  else if( !text.compare("radial") ) {
    vars->rot_method = mod_rotate::METHOD_2D_RADIAL;
    if( param->exists("hdr_azim") ) {
      param->getString("hdr_azim",&hdrName_azim);
    }
//    vars->mode = mod_rotate::MODE_APPLY;  // Fixed to mod_rotate::MODE_APPLY
  }
  else if( !text.compare("from_to") ) {
    vars->rot_method = mod_rotate::METHOD_FROM_TO;
    param->getString("angle_from",&hdrName_from);
    param->getString("angle_to",&hdrName_to);

    if( (vars->orient_xyz == mod_rotate::ORIENT_RIGHTHAND && vars->orient_z == mod_rotate::ORIENT_Z_UP) ||
        (vars->orient_xyz == mod_rotate::ORIENT_LEFTHAND  && vars->orient_z == mod_rotate::ORIENT_Z_DOWN) ) {
      vars->modeEffective = mod_rotate::MODE_APPLY;
      writer->line("Info: Effective mode for 'from-to' rotation has been automatically set to 'REMOVE'");
    }
    else {
      vars->modeEffective = mod_rotate::MODE_REMOVE;  // --> -1
      writer->line("Info: Effective mode for 'from-to' rotation has been automatically set to 'APPLY'");
    }
  }
  else if( !text.compare("tilt_roll") ) {
    vars->rot_method = mod_rotate::METHOD_TILT_ROLL;
    if( param->exists("hdr_roll") ) {
      param->getString("hdr_roll",&hdrName_roll);
    }
    if( param->exists("hdr_tilt") ) {
      param->getString("hdr_tilt",&hdrName_tilt);
   }
  }
  else if( !text.compare("tilt_xy") ) {
    vars->rot_method = mod_rotate::METHOD_TILT_XY;
    param->getFloat("max_tiltxy", &vars->max_valid_tilt, 0 );
    string text;
    param->getString("max_tiltxy", &text, 1 );
    if( !text.compare("yes") ) {
      vars->log_tiltxy_warnings = true;
    }
    else if( !text.compare("no") ) {
      vars->log_tiltxy_warnings = true;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
    if( param->exists("hdr_tiltx") ) {
      param->getString("hdr_tiltx",&hdrName_tiltx);
    }
    if( param->exists("hdr_tilty") ) {
      param->getString("hdr_tilty",&hdrName_tilty);
    }
  }
  else if( !text.compare("icv") ) {
    vars->rot_method = mod_rotate::METHOD_ICV;
  }
  else if( !text.compare("armss_tilt_roll") ) {
    vars->rot_method = mod_rotate::METHOD_ARMSS_TILT_ROLL;
    if( param->exists("hdr_roll") ) {
      param->getString("hdr_roll",&hdrName_roll);
    }
    if( param->exists("hdr_tilt") ) {
      param->getString("hdr_tilt",&hdrName_tilt);
    }
  }
  else if( !text.compare("armss") ) {
    vars->rot_method = mod_rotate::METHOD_ARMSS;
  }
  else if( !text.compare("galperin") ) {
    vars->rot_method = mod_rotate::METHOD_GALPERIN;
  }
  else if( !text.compare("trilobit") ) {
    vars->rot_method = mod_rotate::METHOD_TRILOBIT;
  }
  else {
    writer->error("Unknown method: %s", text.c_str());
  }

  switch( vars->rot_method ) {
    case mod_rotate::METHOD_2D_AZIM:
      mod_rotate::checkHeader( hdrName_azim, hdef, writer );
      vars->hdrId_azim = hdef->headerIndex( hdrName_azim );
      break;
    case mod_rotate::METHOD_FROM_TO:
      mod_rotate::checkHeader( hdrName_from, hdef, writer );
      vars->hdrId_an_from = hdef->headerIndex( hdrName_from );
      mod_rotate::checkHeader( hdrName_to, hdef, writer );
      vars->hdrId_an_to   = hdef->headerIndex( hdrName_to );
      break;
    case mod_rotate::METHOD_2D_RADIAL:
      mod_rotate::checkHeader( hdrName_azim, hdef, writer );
      vars->hdrId_azim = hdef->headerIndex( hdrName_azim );
      mod_rotate::checkHeader( "sou_x", hdef, writer );
      mod_rotate::checkHeader( "sou_y", hdef, writer );
      mod_rotate::checkHeader( "rec_x", hdef, writer );
      mod_rotate::checkHeader( "rec_y", hdef, writer );
      vars->hdrId_soux = hdef->headerIndex("sou_x");
      vars->hdrId_souy = hdef->headerIndex("sou_y");
      vars->hdrId_recx = hdef->headerIndex("rec_x");
      vars->hdrId_recy = hdef->headerIndex("rec_y");
      break;
    case mod_rotate::METHOD_ARMSS_TILT_ROLL:
      mod_rotate::checkHeader( hdrName_roll, hdef, writer );
      mod_rotate::checkHeader( hdrName_tilt, hdef, writer );
      vars->hdrId_roll = hdef->headerIndex( hdrName_roll );
      vars->hdrId_tilt = hdef->headerIndex( hdrName_tilt );
      break;
    case mod_rotate::METHOD_TILT_ROLL:
      mod_rotate::checkHeader( hdrName_roll, hdef, writer );
      mod_rotate::checkHeader( hdrName_tilt, hdef, writer );
      vars->hdrId_roll = hdef->headerIndex( hdrName_roll );
      vars->hdrId_tilt = hdef->headerIndex( hdrName_tilt );
      break;
    case mod_rotate::METHOD_TILT_XY:
      mod_rotate::checkHeader( hdrName_tiltx, hdef, writer );
      mod_rotate::checkHeader( hdrName_tilty, hdef, writer );
      vars->hdrId_tiltx = hdef->headerIndex( hdrName_tiltx );
      vars->hdrId_tilty = hdef->headerIndex( hdrName_tilty );
      break;
    case mod_rotate::METHOD_ICV:
      mod_rotate::checkHeaderFloat("incl_i",hdef,writer);
      mod_rotate::checkHeaderFloat("incl_c",hdef,writer);
      mod_rotate::checkHeaderFloat("incl_v",hdef,writer);
      vars->hdrId_incl_i = hdef->headerIndex( "incl_i" );
      vars->hdrId_incl_c = hdef->headerIndex( "incl_c" );
      vars->hdrId_incl_v = hdef->headerIndex( "incl_v" );
      break;
    default:
      break;
  }

  mod_rotate::checkHeader("sensor",hdef,writer);
  vars->hdrId_sensor = hdef->headerIndex( "sensor" );

//  vars->nTraces = gr->maxdtr;

  //-----------------------------------
//  if( param->exists("input") ) {
    param->getString("input",&text);
    if( !text.compare("xy") ) {
      vars->input = mod_rotate::INPUT_XY;
      edef->setTraceSelectionMode( TRCMODE_FIXED, 2 );
    }
    else if( !text.compare("xyz") ) {
      vars->input = mod_rotate::INPUT_XYZ;
      edef->setTraceSelectionMode( TRCMODE_FIXED, 3 );
    }
    else if( !text.compare("xyzp") ) {
      vars->input = mod_rotate::INPUT_XYZP;
      edef->setTraceSelectionMode( TRCMODE_FIXED, 4 );
    }
    else if( !text.compare("ensemble") ) {
      vars->input = mod_rotate::INPUT_ENSEMBLE;
      edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
//  }
//  else {  // Default setting
//    vars->input = INPUT_XYZ;
//    edef->setTraceSelectionMode( TRCMODE_FIXED, 3 );
//  }

  int rotMode;
  if( vars->mode == mod_rotate::MODE_APPLY ) {
    rotMode = csRotation::ROT_MODE_APPLY;
  }
  else {
    rotMode = csRotation::ROT_MODE_REMOVE;
  }

  vars->rotation = new csRotation( shdr->numSamples, rotMode );
  if( vars->rot_method == mod_rotate::METHOD_TRILOBIT ) {
    vars->rotation->setRotationType( csRotation::ROT_TYPE_TRILOBIT );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_rotate_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;


//---------------------------------------------
  int nTraces = traceGather->numTraces();
  
  int numRotTraces = 3;
  if( vars->input == mod_rotate::INPUT_XYZP ) {
    if( traceGather->numTraces() != 4 ) {
      writer->error(edef->moduleName().c_str(), 0, "Incorrect number of traces in input gather. Expected: 4 (XYZP), found: %d", nTraces );
    }
  }
  else if( vars->input == mod_rotate::INPUT_XYZ ) {
    if( traceGather->numTraces() != 3 ) {
      writer->error( edef->moduleName().c_str(), 0, "Incorrect number of traces in input gather. Expected: 3 (XYZ), found: %d", nTraces );
    }
  }
  else if( vars->input == mod_rotate::INPUT_XY ) {
    if( traceGather->numTraces() != 2 ) {
      writer->error( edef->moduleName().c_str( ), 0, "Incorrect number of traces in input gather. Expected: 2 (XY), found: %d", nTraces );
    }
    numRotTraces = 2;
  }
  int no_value = -999;
  int trace_index[] = {no_value,no_value,no_value};
  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    int sensor = traceGather->trace(itrc)->getTraceHeader()->intValue(vars->hdrId_sensor);
    if( sensor == mod_rotate::SENSOR_INDEX_Z2 ) sensor = mod_rotate::SENSOR_INDEX_Z5;  // Convert sensor index 2 to 5 (both are Z component)
    if( sensor >= mod_rotate::SENSOR_INDEX[mod_rotate::ID_X] && sensor <= mod_rotate::SENSOR_INDEX[mod_rotate::ID_Z] ) {
      if( trace_index[sensor-3] == no_value ) {
        trace_index[sensor-3] = itrc;
      }
      else {
        writer->line("\nTrace header dump of current trace:");
        traceGather->trace(itrc)->getTraceHeader()->dump( writer->getFile() );
        if( sensor == mod_rotate::SENSOR_INDEX[mod_rotate::ID_Z] ) {
          writer->error( edef->moduleName().c_str(), 0, "Input gather contains more than one trace for sensor 2/5 (= Z).");
        }
        else {
          writer->error( edef->moduleName().c_str(), 0, "Input gather contains more than one trace for sensor %d.", sensor);
        }
        return;
      }
    }
  }

  //
  // Read in sensor numbers for input traces
  // Check that all sensors 3, 4 and 5/2 are present
  //
  csTraceHeader** trcHdr = new csTraceHeader*[numRotTraces];
  for( int idSensor = 0; idSensor < numRotTraces; idSensor++ ) {
    if( trace_index[idSensor] == no_value ) {
      writer->line("\nTrace header dump of current trace:");
      traceGather->trace(0)->getTraceHeader()->dump( writer->getFile() );
      if( idSensor == mod_rotate::ID_Z ) {
        writer->error( edef->moduleName().c_str(), 0, "Input gather is missing a sensor 2/5 trace (=Z).");
      }
      else {
        writer->error( edef->moduleName().c_str(), 0, "Input gather is missing a sensor %d trace.", mod_rotate::SENSOR_INDEX[idSensor]);
      }
    }
    else {
      trcHdr[idSensor] = traceGather->trace(trace_index[idSensor])->getTraceHeader();
    }
  }

  if( vars->rot_method == mod_rotate::METHOD_2D_AZIM ) {
    float azim = DEG2RAD( trcHdr[0]->floatValue(vars->hdrId_azim) );
    rotate_xy_azim(
                   traceGather->trace(trace_index[mod_rotate::ID_X])->getTraceSamples(),
                   traceGather->trace(trace_index[mod_rotate::ID_Y])->getTraceSamples(),
                   azim, shdr->numSamples, vars->mode == mod_rotate::MODE_APPLY );
  }
  else if( vars->rot_method == mod_rotate::METHOD_2D_RADIAL ) {
    float geophone_azim = DEG2RAD( trcHdr[mod_rotate::ID_X]->floatValue(vars->hdrId_azim) );
    double sou_x = trcHdr[trace_index[mod_rotate::ID_X]]->doubleValue(vars->hdrId_soux);
    double sou_y = trcHdr[trace_index[mod_rotate::ID_X]]->doubleValue(vars->hdrId_souy);
    double rec_x = trcHdr[trace_index[mod_rotate::ID_X]]->doubleValue(vars->hdrId_recx);
    double rec_y = trcHdr[trace_index[mod_rotate::ID_X]]->doubleValue(vars->hdrId_recy);
    float sr_azim = fmod( (float)atan2( rec_x-sou_x, rec_y-sou_y ) + 2*M_PI, 2*M_PI );  // Positive, clockwise from North
    float rotation_angle = sr_azim - geophone_azim;
    if( vars->mode == mod_rotate::MODE_APPLY ) {
      rotation_angle *= -1.0;
    }
//    fprintf(stderr,"Rotation angle: %f   sr_azim: %f    %f %f\n", rotation_angle*180/M_PI, sr_azim*180/M_PI, sou_x, sou_y);
    // Remember that the following rotation method rotates clock-wise rotation to COORDINATE SYSTEM
    // ...this means, to rotate a vector clock-wise, the angle has to be negative
    rotate_xy_azim(
      traceGather->trace(trace_index[mod_rotate::ID_X])->getTraceSamples(),
      traceGather->trace(trace_index[mod_rotate::ID_Y])->getTraceSamples(),
      rotation_angle, shdr->numSamples, true);
//    trcHdr[trace_index[ID_X]]->setFloatValue( vars->hdrId_azim, RAD2DEG(sr_azim) );
//    trcHdr[trace_index[ID_Y]]->setFloatValue( vars->hdrId_azim, RAD2DEG(sr_azim) + 90.0 );
  }
  else if( vars->rot_method == mod_rotate::METHOD_FROM_TO ) {
    float angle_from_rad = DEG2RAD( trcHdr[0]->floatValue(vars->hdrId_an_from) );
    float angle_to_rad   = DEG2RAD( trcHdr[0]->floatValue(vars->hdrId_an_to) );
    float rotation_angle = angle_to_rad - angle_from_rad;
    if( vars->modeEffective == mod_rotate::MODE_REMOVE ) {
      rotation_angle *= -1.0;
    }
    //    if( vars->mode == mod_rotate::MODE_REMOVE ) {
    //  rotation_angle *= -1.0;
    // }
    rotate_xy(
              traceGather->trace(trace_index[mod_rotate::ID_X])->getTraceSamples(),
              traceGather->trace(trace_index[mod_rotate::ID_Y])->getTraceSamples(),
              rotation_angle, shdr->numSamples );
  }
  else {  // 3D rotation
    float* xTrace = traceGather->trace(trace_index[mod_rotate::ID_X])->getTraceSamples();
    float* yTrace = traceGather->trace(trace_index[mod_rotate::ID_Y])->getTraceSamples();
    float* zTrace = traceGather->trace(trace_index[mod_rotate::ID_Z])->getTraceSamples();

    if( vars->rot_method == mod_rotate::METHOD_TILT_ROLL ) {
      float roll = DEG2RAD( trcHdr[0]->floatValue(vars->hdrId_roll) );
      float tilt = DEG2RAD( trcHdr[0]->floatValue(vars->hdrId_tilt) );
      rotate_xyz_roll_tilt( xTrace, yTrace, zTrace, roll, tilt, shdr->numSamples, vars->mode == mod_rotate::MODE_APPLY);
    }
    else if( vars->rot_method == mod_rotate::METHOD_TILT_XY ) {
      float tiltx = trcHdr[trace_index[mod_rotate::ID_X]]->floatValue(vars->hdrId_tiltx);
      float tilty = trcHdr[trace_index[mod_rotate::ID_Y]]->floatValue(vars->hdrId_tilty);
      if( fabs(tiltx) > vars->max_valid_tilt ) {
        if( vars->log_tiltxy_warnings ) {
          writer->line("Warning: Tilt X exceeds valid maximum: |%f| > %f", tiltx, vars->max_valid_tilt);
        }
        if( vars->tiltx_prev == mod_rotate::NO_TILT ) {
          writer->error("Tilt X exceeds valid maximum: |%f| > %f.\nNo valid tilt X found on previous trace. Repair tilt values first before rotating data", tiltx, vars->max_valid_tilt);
        }
        tiltx = vars->tiltx_prev;
      }
      if( fabs(tilty) > vars->max_valid_tilt ) {
        if( vars->log_tiltxy_warnings ) {
          writer->line("Warning: Tilt Y exceeds valid maximum: |%f| > %f", tilty, vars->max_valid_tilt);
        }
        if( vars->tilty_prev == mod_rotate::NO_TILT ) {
          writer->error("Tilt Y exceeds valid maximum: |%f| > %f.\nNo valid tilt Y found on previous trace. Repair tilt values first before rotating data", tilty, vars->max_valid_tilt);
        }
        tilty = vars->tilty_prev;
      }
      rotate_xyz_tiltxy( xTrace, yTrace, zTrace, DEG2RAD(tiltx), DEG2RAD(tilty), shdr->numSamples, vars->mode == mod_rotate::MODE_APPLY);
      vars->tiltx_prev = tiltx;
      vars->tilty_prev = tilty;
    }
    else if( vars->rot_method == mod_rotate::METHOD_ICV ) {
      float incl_i[3], incl_v[3], incl_c[3];
      for( int idSensor = 0; idSensor < 3; idSensor++ ) {
        incl_i[idSensor] = trcHdr[idSensor]->floatValue(vars->hdrId_incl_i);
        incl_c[idSensor] = trcHdr[idSensor]->floatValue(vars->hdrId_incl_c);
        incl_v[idSensor] = trcHdr[idSensor]->floatValue(vars->hdrId_incl_v);
      }
      rotate_xyz_icv( xTrace, yTrace, zTrace, incl_i, incl_c, incl_v, shdr->numSamples, vars->mode);
    }
    else if( vars->rot_method == mod_rotate::METHOD_GALPERIN ) {
      vars->rotation->rotate3d_galperin( xTrace, yTrace, zTrace );
    }
    else if( vars->rot_method == mod_rotate::METHOD_ARMSS ) {
      vars->rotation->rotate3d_armss( xTrace, yTrace, zTrace );
    }
    else if( vars->rot_method == mod_rotate::METHOD_ARMSS_TILT_ROLL ) {
      float armss_roll = DEG2RAD( trcHdr[0]->floatValue(vars->hdrId_roll) );
      float armss_tilt = DEG2RAD( trcHdr[0]->floatValue(vars->hdrId_tilt) );
//      printf("ROll tilt: %f %f\n", trcHdr[0]->floatValue(vars->hdrId_roll), trcHdr[0]->floatValue(vars->hdrId_tilt) );
      vars->rotation->rotate3d_armss( xTrace, yTrace, zTrace, armss_roll, armss_tilt );
    }
    else if( vars->rot_method == mod_rotate::METHOD_TRILOBIT ) {
      vars->rotation->rotate3d_trilobit( xTrace, yTrace, zTrace );
    }
  } // END 3D rotation
  
  delete [] trcHdr;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_rotate_( csParamDef* pdef ) {
  pdef->setModule( "ROTATE", "Perform 2D/3D rotation to input traces", "For 2D rotation (i.e. azimuth), the coordinate system is assumed to be X/East, Y/North" );

  pdef->addDoc("The values found in the input headers must be in degrees. Rotation is done by this angle.");
  pdef->addDoc("Note that the COORDINATE system is rotated clock-wise by the given rotation angle, not the data VECTOR.");

  pdef->addParam( "input", "Input data", NUM_VALUES_FIXED, "Specifies what type of input data is expected in this module" );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "xy", "2 sequential input traces, containing sensors 3(X) and 4(Y). Sorting of these two may be arbitrary." );
  pdef->addOption( "xyz", "3 sequential input traces, containing sensors 3(X), 4(Y) and 5/2(Z). Sorting of these three may be arbitrary." );
  pdef->addOption( "xyzp", "4 sequential input traces, containing sensors 3(X), 4(Y) and 5/2(Z), and 4th trace with sensor 1(hydrophone). Sorting of these four may be arbitrary." );
  pdef->addOption( "ensemble", "Input whole ensemble (user sorted). Each input ensemble must contain one trace each for sensor 3(X), 4(Y) and 5/2(Z). Sorting of these three may be arbitrary." );

  pdef->addParam( "method", "Rotation method", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "azim", "Rotate input trace pair XY by azimuth angle.", "Required input header: an_azim" );
  pdef->addOption( "radial", "Rotate input trace pair XY to radial/transverse.", "Required input headers: an_azim, sou_x, sou_y, rec_x, rec_y" );
  pdef->addOption( "icv", "Rotate input trace triplet XYZ by inclinometer ICV values.", "Required input headers: incl_i, incl_c, incl_v" );
  pdef->addOption( "tilt_roll", "Rotate input trace triplet XYZ by roll/tilt angles.", "Required input headers: an_roll, an_tilt" );
  pdef->addOption( "tilt_xy", "Rotate input trace triplet XYZ by tilt X/Y angles.", "Required input headers: an_tiltx, an_tilty" );
  pdef->addOption( "galperin", "Rotate to galperin configuration", "Geophones point outwards down, first geophone's azimuth is unchanged by rotation" );
  pdef->addOption( "armss", "Rotate ARMSS node data: Galperin configuration and 162deg ARMSS azimuth" );
  pdef->addOption( "armss_tilt_roll", "Rotate ARMSS node data: Galperin, azimuth and ARMSS roll & tilt", "Required input header: an_roll, an_tilt" );
  pdef->addOption( "trilobit", "Rotate Trilobit node data from Galperin g123 configuration to XYZ" );
  pdef->addOption( "from_to", "Rotate XY trace pair (=XY coordinate system) from angle 'angle_from' to 'angle_to'. Angles usually give X orientation" );

  pdef->addParam( "orient", "Orientation of XYZ/XY components", NUM_VALUES_FIXED );
  pdef->addValue( "right", VALTYPE_OPTION );
  pdef->addOption( "right", "XYZ components describe a right hand system" );
  pdef->addOption( "left", "XYZ components describe a left hand system" );
  pdef->addValue( "up", VALTYPE_OPTION );
  pdef->addOption( "up", "Z component is pointing downwards", "For XY input only, the orientation of a 'virtual' Z component defines the relative orientation of XY in the same way" );
  pdef->addOption( "down", "Z component is pointing upwards", "For XY input only, the orientation of a 'virtual' Z component defines the relative orientation of XY in the same way" );

  pdef->addParam( "hdr_roll", "Trace header containing roll angle", NUM_VALUES_FIXED );
  pdef->addValue( "an_roll", VALTYPE_STRING, "Trace header name containing roll angle [deg]" );

  pdef->addParam( "hdr_tilt", "Trace header containing tilt angle", NUM_VALUES_FIXED );
  pdef->addValue( "an_tilt", VALTYPE_STRING, "Trace header name containing tilt angle [deg]" );

  pdef->addParam( "hdr_tiltx", "Trace header containing tilt X angle", NUM_VALUES_FIXED );
  pdef->addValue( "an_tiltx", VALTYPE_STRING, "Trace header name containing tilt X angle [deg]" );

  pdef->addParam( "hdr_tilty", "Trace header containing tilt Y angle", NUM_VALUES_FIXED );
  pdef->addValue( "an_tilty", VALTYPE_STRING, "Trace header name containing tilt Y angle [deg]" );

  pdef->addParam( "hdr_azim", "Trace header containing azimuth", NUM_VALUES_FIXED );
  pdef->addValue( "an_azim", VALTYPE_STRING, "Trace header name containing azimuth [deg]" );

  pdef->addParam( "angle_from", "Angle that input data is rotated to (usually direction of X component)", NUM_VALUES_FIXED, "Azimuth angle is defined clock-wise from North" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name containing angle [deg]" );
  pdef->addParam( "angle_to", "Angle that data shall be rotated to (usually output direction of X component)", NUM_VALUES_FIXED, "Azimuth angle is defined clock-wise from North" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name containing angle [deg]" );

  pdef->addParam( "max_tiltxy", "Maximum valid tilt value for tilt XY rotation", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum valid tilt value", "If tilt exceeds the maximum valid tilt value, it is replaced by the value from previous trace, if available" );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Write warnings to job log" );
  pdef->addOption( "no", "Do not write warnings to job log" );

  pdef->addParam( "mode", "Mode of operation", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply rotation.", "For 2D rotation, apply clockwise rotation to coordinate system (X/East and Y/North)." );
  pdef->addOption( "remove", "Remove rotation. For field data, specify 'remove' to rotate from as-laid to final coordinate system.", "For 2D rotation, apply ANTI-clockwise rotation to coordinate system (X/East and Y/North)." );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_rotate_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_rotate::VariableStruct* vars = reinterpret_cast<mod_rotate::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_rotate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_rotate::VariableStruct* vars = reinterpret_cast<mod_rotate::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->rotation != NULL ) {
    delete vars->rotation;
    vars->rotation = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_rotate_( csParamDef* pdef ) {
  params_mod_rotate_( pdef );
}
extern "C" void _init_mod_rotate_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_rotate_( param, env, writer );
}
extern "C" bool _start_exec_mod_rotate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_rotate_( env, writer );
}
extern "C" void _exec_mod_rotate_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_rotate_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_rotate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_rotate_( env, writer );
}

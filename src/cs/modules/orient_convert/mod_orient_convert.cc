/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ORIENT_CONVERT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_orient_convert {
  struct VariableStruct {
    int hdrId_incl_i;
    int hdrId_incl_c;
    int hdrId_incl_v;
    int hdrId_roll;
    int hdrId_tilt;
    int hdrId_sensor;
    int hdrId_ain1;
    int hdrId_ain2;
    int hdrId_ain3;
    int method;
    int input;
  };

  static const int TILTROLL_2_ICV = 10;
  static const int ICV_2_TILTROLL = 11;
  static const int ARMSS = 12;
  static const int INPUT_XYZ  = 20;
  static const int INPUT_XYZP = 21;
  static const int INPUT_ENSEMBLE = 22;
  static const int INPUT_SINGLE_TRACE = 30;
  static const int SENSOR_INDEX[] = { 3, 4, 5, 1 };
  static const int ID_X = 0;
  static const int ID_Y = 1;
  static const int ID_Z = 2;
  static const int ID_P = 3;

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
using namespace mod_orient_convert;


//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_orient_convert_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 3 );

  // Set defaults
  vars->hdrId_incl_i = -1;
  vars->hdrId_incl_c = -1;
  vars->hdrId_incl_v = -1;
  vars->hdrId_roll   = -1;
  vars->hdrId_tilt   = -1;
  vars->hdrId_sensor = -1;
  vars->hdrId_ain1   = -1;
  vars->hdrId_ain2   = -1;
  vars->hdrId_ain3   = -1;
  vars->method       = -1;
  vars->input        = INPUT_XYZ;

//----------------------------------------------------
  std::string text;

  param->getString( "method", &text );
  text = toLowerCase(text);
  if( !text.compare("tiltroll_2_icv") ) {
    vars->method = TILTROLL_2_ICV;
  }
  else if( !text.compare("icv_2_tiltroll") ) {
    vars->method = ICV_2_TILTROLL;
  }
  else if( !text.compare("armss") ) {
    vars->method = ARMSS;
    // Override input method
    vars->input = INPUT_SINGLE_TRACE;
    edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }
  else {
    writer->error("Method option not recognized: %s.", text.c_str());
  }

  //-----------------------------------
  if( param->exists("input") ) {
    if( vars->method == ARMSS ) {
      writer->warning("'input' parameter ignored for ARMSS accelerometer conversion. Single trace input.");
    }
    else {
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
  }

//----------------------------------------------------
  if( vars->method == ICV_2_TILTROLL ) {
    checkHeader( "incl_i", hdef, writer );
    checkHeader( "incl_c", hdef, writer );
    checkHeader( "incl_v", hdef, writer );

    vars->hdrId_incl_i = hdef->headerIndex( "incl_i" );
    vars->hdrId_incl_c = hdef->headerIndex( "incl_c" );
    vars->hdrId_incl_v = hdef->headerIndex( "incl_v" );

    if( hdef->headerExists ("an_roll") ) {
      checkHeader( "an_roll", hdef, writer );
      vars->hdrId_roll = hdef->headerIndex("an_roll" );
    }
    else {
      vars->hdrId_roll = hdef->addStandardHeader( "an_roll" );
    }
    if( hdef->headerExists ("an_tilt") ) {
      vars->hdrId_tilt = hdef->headerIndex("an_tilt" );
    }
    else {
      vars->hdrId_tilt = hdef->addStandardHeader( "an_tilt" );
    }
  }
  else if( vars->method == TILTROLL_2_ICV ) {
    checkHeader( "an_roll", hdef, writer );
    checkHeader( "an_tilt", hdef, writer );
    vars->hdrId_roll = hdef->headerIndex( "an_roll" );
    vars->hdrId_tilt = hdef->headerIndex( "an_tilt" );

    if( !hdef->headerExists ("incl_i") ) {
      hdef->addStandardHeader( "incl_i" );
    }
    if( !hdef->headerExists ("incl_c") ) {
      hdef->addStandardHeader( "incl_c" );
    }
    if( hdef->headerExists ("incl_v") ) {
      hdef->addStandardHeader( "incl_v" );
    }
    vars->hdrId_incl_i = hdef->headerIndex("incl_i" );
    vars->hdrId_incl_c = hdef->headerIndex("incl_c" );
    vars->hdrId_incl_v = hdef->headerIndex("incl_v" );

  }
  else {
    checkHeader( "ain1", hdef, writer );
    checkHeader( "ain2", hdef, writer );
    checkHeader( "ain3", hdef, writer );

    vars->hdrId_ain1 = hdef->headerIndex( "ain1" );
    vars->hdrId_ain2 = hdef->headerIndex( "ain2" );
    vars->hdrId_ain3 = hdef->headerIndex( "ain3" );

    if( !hdef->headerExists ("an_roll") ) {
      hdef->addStandardHeader( "an_roll" );
    }
    if( !hdef->headerExists ("an_tilt") ) {
      hdef->addStandardHeader( "an_tilt" );
    }
    vars->hdrId_roll = hdef->headerIndex("an_roll" );
    vars->hdrId_tilt = hdef->headerIndex("an_tilt" );
  }

//-------------------------------------------------------------
  if( !hdef->headerExists("sensor") ){
    writer->error("Trace header 'sensor' was not found in the trace headers.");
  }
  vars->hdrId_sensor = hdef->headerIndex( "sensor" );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_orient_convert_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;


//---------------------------------------------
  double roll, tilt;
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
  else if( vars->input == INPUT_SINGLE_TRACE ) {
    if( traceGather->numTraces() != 1 ) {
      writer->error("Incorrect number of traces in input gather. Expected: 1, found: %d", nTraces );
    }
    csTrace* trace = traceGather->trace(0);
    csTraceHeader* trcHeader = trace->getTraceHeader();
    double ain1 = trcHeader->doubleValue(vars->hdrId_ain1);
    double ain2 = trcHeader->doubleValue(vars->hdrId_ain2);
    double ain3 = trcHeader->doubleValue(vars->hdrId_ain3);

    if( edef->isDebug() ) writer->line( "ain1 ain2 ain3: %10.4f %10.4f %10.f\n", ain1, ain2, ain3 );

    // Rotate from system fixed coordinate system to XYZ coordinate system
    // --> Rotate clock-wise by 42deg, and flip some polarities
    double a_rad = DEG2RAD( 42.0 );
    double ay =  cos(a_rad)*ain2 - sin(a_rad)*ain3;
    double ax =  sin(a_rad)*ain2 + cos(a_rad)*ain3;
    double az = -ain1;

    double roll   = atan2(ay,-az);
    double az_rot = -sin(roll)*ay + cos(roll)*az;
    double tilt   = RAD2DEG( atan2(ax,-az_rot) );
    roll = RAD2DEG( roll );
    roll = -roll;  // Follow same convention as Sercel
    tilt = -tilt;  // Follow same convention as Sercel

    if( edef->isDebug() ) writer->line("roll tilt: %10.4f %10.4f   %5d  %5d\n", roll, tilt, NINT(roll*100), NINT(tilt*100) );
    trcHeader->setFloatValue( vars->hdrId_roll, (float)roll );
    trcHeader->setFloatValue( vars->hdrId_tilt, (float)tilt );

    return;
  }

  int no_value = -999;
  int trace_index[] = {no_value,no_value,no_value};
  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    int sensor = traceGather->trace(itrc)->getTraceHeader()->intValue(vars->hdrId_sensor);
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

  csTraceHeader* trcHdr[3];

  for( int idSensor = 0; idSensor < 3; idSensor++ ) {
    if( trace_index[idSensor] == no_value ) {
      writer->error("Input gather is missing a sensor %d trace.", SENSOR_INDEX[idSensor]);
    }
    else {
      trcHdr[idSensor] = traceGather->trace(trace_index[idSensor])->getTraceHeader();
    }
  }

  // Read in header values and perform conversion
  if( vars->method == TILTROLL_2_ICV ) {
//    for( int idSensor = 0; idSensor < 3; idSensor++ ) {
      roll = DEG2RAD( trcHdr[0]->doubleValue(vars->hdrId_roll) );
      tilt = DEG2RAD( trcHdr[0]->doubleValue(vars->hdrId_tilt) );
//    }
    trcHdr[ID_X]->setDoubleValue( vars->hdrId_incl_i, cos(tilt) );
    trcHdr[ID_X]->setDoubleValue( vars->hdrId_incl_c, 0.0 );
    trcHdr[ID_X]->setDoubleValue( vars->hdrId_incl_v, -sin(tilt) );

    trcHdr[ID_Y]->setDoubleValue( vars->hdrId_incl_i, cos(roll)*sin(tilt) );
    trcHdr[ID_Y]->setDoubleValue( vars->hdrId_incl_c, sin(roll) );
    trcHdr[ID_Y]->setDoubleValue( vars->hdrId_incl_v, cos(roll)*cos(tilt) );

    trcHdr[ID_Z]->setDoubleValue( vars->hdrId_incl_i, sin(roll)*sin(tilt) );
    trcHdr[ID_Z]->setDoubleValue( vars->hdrId_incl_c, -cos(roll) );
    trcHdr[ID_Z]->setDoubleValue( vars->hdrId_incl_v, sin(roll)*cos(tilt) );
  }
  else if( vars->method == ICV_2_TILTROLL ) {
    double value_i[3], value_c[3], value_v[3];
    for( int idSensor = 0; idSensor < 3; idSensor++ ) {
      value_i[idSensor] = DEG2RAD( trcHdr[idSensor]->doubleValue(vars->hdrId_incl_i) );
      value_c[idSensor] = DEG2RAD( trcHdr[idSensor]->doubleValue(vars->hdrId_incl_c) );
      value_v[idSensor] = DEG2RAD( trcHdr[idSensor]->doubleValue(vars->hdrId_incl_v) );
    }
    tilt = RAD2DEG(atan2( -value_v[ID_X], value_i[ID_X] ));
    roll = RAD2DEG(atan2( value_c[ID_Y], -value_c[ID_Z] )) + 360.0;
    roll -= 360.0*(double)( (int)(roll/360.0) );
    for( int itrc = 0; itrc < nTraces; itrc++ ) {
      csTraceHeader* trcHdrTmp = traceGather->trace(itrc)->getTraceHeader();
      trcHdrTmp->setDoubleValue( vars->hdrId_tilt, tilt );
      trcHdrTmp->setDoubleValue( vars->hdrId_roll, roll );
    }
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_orient_convert_( csParamDef* pdef ) {
  pdef->setModule( "ORIENT_CONVERT", "Convert sensor orientation parameters", "Applies some process to trace headers dealing with sensor orientation" );

  pdef->addParam( "input", "Input data", NUM_VALUES_FIXED, "Specifies what type of input data is expected in this module" );
  pdef->addValue( "xyz", VALTYPE_OPTION );
  pdef->addOption( "xyz", "3 sequential input traces, containing sensors 3(X), 4(Y) and 5(Z). Sorting of these three may be arbitrary." );
  pdef->addOption( "xyzp", "4 sequential input traces, containing sensors 3(X), 4(Y) and 5(Z), and 4th trace with sensor 1(hydrophone). Sorting of these three may be arbitrary." );
  pdef->addOption( "ensemble", "Input whole ensemble (user sorted). Each input ensemble must contain one trace each for sensor 3(X), 4(Y) and 5(Z). Sorting of these three may be arbitrary." );

  pdef->addParam( "method", "", NUM_VALUES_VARIABLE, "Time knee points at which specified scalars apply" );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "tiltroll_2_icv", "Convert roll and tilt angles into ICV inclinometer values.", "Required input headers: an_roll, an_tilt" );
  pdef->addOption( "icv_2_tiltroll", "Convert ICV inclinometer values into roll and tilt angles.", "Required input headers: incl_i, incl_c, incl_v" );
  pdef->addOption( "armss", "Convert ARMSS accelerometer values into (ARMSS) roll and tilt angles.", "Required input headers: ain1, ain2, ain3" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_orient_convert_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_orient_convert::VariableStruct* vars = reinterpret_cast<mod_orient_convert::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_orient_convert_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_orient_convert::VariableStruct* vars = reinterpret_cast<mod_orient_convert::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_orient_convert_( csParamDef* pdef ) {
  params_mod_orient_convert_( pdef );
}
extern "C" void _init_mod_orient_convert_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_orient_convert_( param, env, writer );
}
extern "C" bool _start_exec_mod_orient_convert_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_orient_convert_( env, writer );
}
extern "C" void _exec_mod_orient_convert_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_orient_convert_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_orient_convert_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_orient_convert_( env, writer );
}

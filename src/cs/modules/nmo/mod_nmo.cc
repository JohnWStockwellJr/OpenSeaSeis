/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csNMOCorrection.h"
#include "csTimeFunction.h"
#include "csTableManagerNew.h"
#include "csTableAll.h"
#include "csTableNew.h"
#include "csTableValueList.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: NMO
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_nmo {
  struct VariableStruct {
    float* time_sec;  // [s]
    float* velocities;  // [m/s]
    csNMOCorrection* nmo;
    int numTimes;
    int hdrId_offset;
    int hdrId_velocity;
    int hdrId_diff_offset;
    int mode;
    float percentVel; // Percent velocity change to be applied to input velocity
    double offset_diffNMO;
    bool isDiffNMO;
    bool dump;
    csTableNew* table;
    csTableManagerNew* oldTableManager;
    int* hdrId_keys;
    float timeSample1_ms;
  };
}
using mod_nmo::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_nmo_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


//---------------------------------------------
//
  vars->numTimes       = 0;
  vars->time_sec          = NULL;
  vars->velocities     = NULL;
  vars->nmo            = NULL;
  vars->hdrId_offset   = -1;
  vars->hdrId_velocity = -1;
  vars->hdrId_diff_offset = -1;
  vars->isDiffNMO = false;
  vars->offset_diffNMO = 0;
  vars->table          = NULL;
  vars->oldTableManager   = NULL;
  vars->hdrId_keys     = NULL;
  vars->dump           = false; 
  vars->percentVel = 0.0;
  vars->timeSample1_ms = 0.0;

  csVector<std::string> valueList;
  
  //-------------------------------------
  std::string text;
  int mode_nmo = csNMOCorrection::PP_NMO;
  float offsetApex_m = 4000;
  float zeroOffsetDamping = 5;
  if( param->exists("wave_mode") ) {
    param->getString( "wave_mode", &text );
    if( !text.compare("pp_iso") ) {
      mode_nmo = csNMOCorrection::PP_NMO;
    }
    else if( !text.compare("ps_iso") ) {
      mode_nmo = csNMOCorrection::PS_NMO;
    }
    else if( !text.compare("pp_vti") ) {
      mode_nmo = csNMOCorrection::PP_NMO_VTI;
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }
  }
  if( param->exists("empirical") ) {
    mode_nmo = csNMOCorrection::EMPIRICAL_NMO;
    param->getFloat("empirical", &offsetApex_m, 0 );
    param->getFloat("empirical", &zeroOffsetDamping, 1 );
  }
  if( param->exists("output_vel") ) {
    param->getString( "output_vel", &text );
    if( !text.compare("no") ) {
      // Do nothing, keep NMO mode set earlier
    }
    else if( !text.compare("yes") ) {
      mode_nmo = csNMOCorrection::OUTPUT_VEL;
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }
  }
  vars->mode = csNMOCorrection::NMO_APPLY;
  if( param->exists("mode") ) {
    param->getString( "mode", &text );
    if( !text.compare("apply") ) {
      vars->mode = csNMOCorrection::NMO_APPLY;
    }
    else if( !text.compare("remove") ) {
      vars->mode = csNMOCorrection::NMO_REMOVE;
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }
  }

  if( param->exists("dump") ) {
    param->getString( "dump", &text );
    if( !text.compare("yes") ) {
      vars->dump = true;
    }
    else if( !text.compare("no") ) {
      vars->dump = false;
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }

  }
//---------------------------------------------
  if( param->exists("diff_nmo_offset") ) {
    param->getString("diff_nmo_offset",&text);
    csFlexNumber number;
    if( !number.convertToNumber( text ) ) {
      vars->hdrId_diff_offset = hdef->headerIndex(text);
    }
    else {
      vars->offset_diffNMO = number.doubleValue();
    }
    vars->isDiffNMO = true;
  }
//---------------------------------------------
  if( param->exists("percent") ) {
    param->getFloat("percent",&vars->percentVel);
    if( vars->percentVel <= -100 || vars->percentVel >=100 ) writer->error("Specified percent velocity change (=%.2f%%) out of valid range (-99%% - +99%%)");
  }
//---------------------------------------------
  if( param->exists("time_samp1") ) {
    param->getFloat("time_samp1",&vars->timeSample1_ms);
    if( vars->timeSample1_ms > 0 ) writer->error("Specified time of first sample must be smaller than or equal to 0. Specified: %f", vars->timeSample1_ms);
    if( vars->timeSample1_ms != 0 && vars->isDiffNMO ) writer->error("Differential NMO is currently not supported in the case where the time of first sample is not zero");
  }
//---------------------------------------------
// Retrieve velocity table
//
  if( param->exists("table_old") ) {
    param->getString("table_old", &text );
    try {
      vars->oldTableManager = new csTableManagerNew( text, csTableAll::TABLE_TYPE_TIME_FUNCTION, hdef );
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s': %s\n", text.c_str(), exc.getMessage() );
    }
    if( edef->isDebug() ) {
      vars->oldTableManager->dump();
    }
  }
  if( param->exists("table") ) {
    if( vars->oldTableManager != NULL ) writer->error("Specify either user parameter 'table' or 'table_old', not both");
    std::string tableFilename;
    int colTime = 0;
    int colVel  = 1;
    int colEta  = -1;
    param->getString("table", &tableFilename );
    param->getInt("table_col",&colTime,0);
    param->getInt("table_col",&colVel,1);
    if( colTime < 1 || colVel < 1 ) writer->error("Column numbers in table (user parameter 'table_col') must be larger than 0");
    if( param->getNumValues("table_col") > 2 ) {
      param->getInt("table_col",&colEta,2);
      if( colEta < 1 ) writer->error("Column number for eta (user parameter 'table_col') must be larger than 0");
    }
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
        vars->table->addKey( col-1, interpolate );  // -1 to convert from 'user' column to 'C++' column
        if( !hdef->headerExists( headerName ) ) {
          writer->error("No matching trace header found for table key '%s'", headerName.c_str() );
        }
        vars->hdrId_keys[ikey] = hdef->headerIndex( headerName );
      } // END for ikey
    }

    // Add velocity 'value' column to input table
    vars->table->addValue( colVel-1 );  // -1 to convert from 'user' column to 'C++' column
    if( colEta > 0 ) vars->table->addValue( colEta-1 );

    bool sortTable = false;
    try {
      vars->table->initialize( tableFilename, sortTable );
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s': %s\n", text.c_str(), exc.getMessage() );
    }
  }
  
//---------------------------------------------
// Retrieve velocity-time pairs
//
  if( vars->table == NULL && vars->oldTableManager == NULL ) {
    if( param->exists( "time" ) ) {
      param->getAll( "time", &valueList );
      if( valueList.size() < 1 ){
        writer->error("No times specified in user parameter 'time'!");
      }
    }
    vars->numTimes = valueList.size();

    if( vars->numTimes > 0 ) {
      vars->time_sec = new float[vars->numTimes];
      for( int i = 0; i < vars->numTimes; i++ ) {
        vars->time_sec[i] = atof( valueList.at(i).c_str() ) / 1000.0;  // Convert to seconds
      }
    }
    else {
      vars->numTimes = 1;
      vars->time_sec = new float[vars->numTimes];
      vars->time_sec[0] = 0;
    }

    vars->velocities = new float[vars->numTimes];

    valueList.clear();
    param->getAll( "velocity", &valueList );

    if( valueList.size() < 1 ){
      writer->error("Missing user parameter 'velocity'!");
    }
    else if( valueList.size() != vars->numTimes ) {
      writer->error("Unequal number of velocities(%d) and times(%d)", valueList.size(), vars->numTimes );
    }
    csFlexNumber number;
    for( int i = 0; i < vars->numTimes; i++ ) {
      if( !number.convertToNumber( valueList.at(i) ) ) {
        text = valueList.at(i);
        if( vars->numTimes > 1 ) writer->error("Specified velocity is not a valid number: '%s'", text.c_str() );
        if( !hdef->headerExists(text) ) {
          writer->error("Specified trace header name containing velocity does not exist: '%s'", text.c_str());
        }
        vars->hdrId_velocity = hdef->headerIndex(text);
        break;
      }
      vars->velocities[i] = number.floatValue();
      if( vars->percentVel != 0.0 ) vars->velocities[i] *= (100.0+vars->percentVel)/100.0;
      if( edef->isDebug() ) writer->line("Velocity #%d: '%s' --> %f", i, valueList.at(i).c_str(), vars->velocities[i] );
    }
  }
  
  //---------------------------------------------
  // Set NMO object
  //
  vars->nmo = new csNMOCorrection( shdr->sampleInt, shdr->numSamples, mode_nmo );
  vars->nmo->setModeOfApplication( vars->mode );
  if( mode_nmo == csNMOCorrection::EMPIRICAL_NMO ) {
    vars->nmo->setEmpiricalNMO( offsetApex_m, zeroOffsetDamping );
  }
  if( vars->timeSample1_ms != 0 ) {
    vars->nmo->setTimeSample1( vars->timeSample1_ms );
  }

  if( param->exists("horizon_nmo") ) {
    param->getString( "horizon_nmo", &text );
    if( !text.compare("yes") ) {
      int horMethod = csNMOCorrection::HORIZON_METHOD_LINEAR;
      if( param->getNumValues("horizon_nmo") > 1 ) {
        param->getString("horizon_nmo", &text, 1);
        if( !text.compare("linear") ) {
          horMethod = csNMOCorrection::HORIZON_METHOD_LINEAR;
        }
        else if( !text.compare("quadratic") ) {
          horMethod = csNMOCorrection::HORIZON_METHOD_QUAD;
        }
        else {
          writer->error("Option not recognized: %s.", text.c_str());
        }
      }
      vars->nmo->setHorizonBasedNMO( true, horMethod );
    }
    else if( !text.compare("no") ) {
      vars->nmo->setHorizonBasedNMO( false );
    }
    else {
      writer->error("Option not recognized: %s.", text.c_str());
    }
  }
  
  if( !hdef->headerExists( "offset" ) ) {
    writer->error("Trace header 'offset' does not exist. Cannot perform NMO correction.");
  }
  else if( hdef->headerType( "offset" ) != TYPE_FLOAT && hdef->headerType( "offset" ) != TYPE_DOUBLE ) {
    writer->error("Trace header 'offset' exists but has the wrong number type. Should be FLOAT.");
  }

  vars->hdrId_offset = hdef->headerIndex( "offset" );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_nmo_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csTrace* trace = traceGather->trace(0);


  float* samples = trace->getTraceSamples();
  double offset  = trace->getTraceHeader()->doubleValue( vars->hdrId_offset );

  if( vars->hdrId_velocity >= 0 ) {
    vars->velocities[0] = trace->getTraceHeader()->floatValue(vars->hdrId_velocity);
    if( vars->percentVel != 0.0 ) vars->velocities[0] *= (100.0+vars->percentVel)/100.0;
  }
  if( vars->hdrId_diff_offset >= 0 ) {
    vars->offset_diffNMO  = trace->getTraceHeader()->doubleValue( vars->hdrId_diff_offset );
  }

  if( vars->table != NULL ) {
    csTimeFunction<double> const* timeFunc;
    if( vars->table->numKeys() > 0 ) {
      double* keyValueBuffer = new double[vars->table->numKeys()];
      for( int ikey = 0; ikey < vars->table->numKeys(); ikey++ ) {
        keyValueBuffer[ikey] = trace->getTraceHeader()->doubleValue( vars->hdrId_keys[ikey] );
      }
      timeFunc = vars->table->getFunction( keyValueBuffer, vars->dump );
      delete [] keyValueBuffer;
    }
    else {
      timeFunc = vars->table->getFunction( NULL, vars->dump );
    }
    if( !vars->isDiffNMO ) {
      vars->nmo->perform_nmo( timeFunc, offset, samples );
    }
    else {
      vars->nmo->perform_differential_nmo( timeFunc, offset, vars->offset_diffNMO, samples );
    }
  }
  else if( vars->oldTableManager != NULL ) {
    csTimeFunction<double> const* timeFunc = vars->oldTableManager->getFunction( trace->getTraceHeader() );
    if( !vars->isDiffNMO ) {
      vars->nmo->perform_nmo( timeFunc, offset, samples );
    }
    else {
      vars->nmo->perform_differential_nmo( timeFunc, offset, vars->offset_diffNMO, samples );
    }
  }
  else if( !vars->isDiffNMO ) {
    vars->nmo->perform_nmo( vars->numTimes, vars->time_sec, vars->velocities, offset, samples );
  }
  else {
    vars->nmo->perform_differential_nmo( vars->numTimes, vars->time_sec, vars->velocities, offset, vars->offset_diffNMO, samples );
  }
  
  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_nmo_( csParamDef* pdef ) {
  pdef->setModule( "NMO", "Perform NMO (normal moveout) correction" );

  pdef->addDoc("When reading in the velocity field from an ASCII table, velocity-time functions are interpolated at the trace locations based on the");
  pdef->addDoc("specified 'table_key' key trace headers. The interpolation assumes a regular cartesian grid, oriented W-E/N-S.");
  pdef->addDoc("This will give incorrect results if the grid has an oblique orientation and if XY coordinates are used as keys. If inline/crossline numbers");
  pdef->addDoc("are used as keys, the interpolation will work correctly regardless of the grid's orientation.");
  pdef->addDoc("If XY coordinates must be used as keys, the work-around is to re-grid the input velocity field to a W-E/N-S grid.");

  pdef->addParam( "wave_mode", "Mode of wave to model", NUM_VALUES_FIXED );
  pdef->addValue( "pp_iso", VALTYPE_OPTION );
  pdef->addOption( "pp_iso", "Isotropic PP mode (downgoing P, upgoing P)" );
  pdef->addOption( "ps_iso", "Isotropic PS mode (downgoing P, upgoing S)" );
  pdef->addOption( "pp_vti", "VTI anisotropic PP mode, using Alkhalifah and Tsvankin (1995)", "Supply an eta field along with the (zero-offset) NMO velocity field. Provide ASCII table with at least thee columns: time vel eta" );

  pdef->addParam( "time", "List of time value [ms]", NUM_VALUES_VARIABLE, "Time knee points at which specified velocities apply" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of time values [ms]." );

  pdef->addParam( "velocity", "List of NMO velocities [m/s]", NUM_VALUES_VARIABLE,
                  "Velocities values at specified time knee points. In between time knee points, velocities will be linearly interpolated." );
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "List of velocities [m/s]..." );

  pdef->addParam( "table", "Velocity table", NUM_VALUES_FIXED,
                  "Velocity table format: The ASCII file should contain only numbers, no text. Required are two columns containing time and velocity pairs. Optional: Up to 2 key values. Lines starting with '#' are considered comment lines." );
  pdef->addValue( "", VALTYPE_STRING, "Filename containing velocity table.");

  pdef->addParam( "table_col", "Table columns containing time in [ms] and NMO velocity in [m/s]", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing time [ms]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table containing NMO velocity [m/s]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Optional: Column number in input table containing additional parameter", "For example, for VTI the eta parameter is required" );

  pdef->addParam( "table_key", "Key trace header used to match values found in specified table columns", NUM_VALUES_VARIABLE,
                  "Specify the 'key' parameter for each key in the velocity file" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name of key header" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table" );
  pdef->addValue( "yes", VALTYPE_OPTION, "Interpolate based to this key?" );
  pdef->addOption( "yes", "Use this key for interpolation of value" );
  pdef->addOption( "no", "Do not use this key for interpolation", "The input table is expected to contain the exact key values for this trace header" );

  pdef->addParam( "mode", "Mode of NMO application", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply NMO." );
  pdef->addOption( "remove", "Remove NMO." );

  pdef->addParam( "empirical", "Empirical 'NMO' function", NUM_VALUES_FIXED,
                  "Empirical NMO is a residual NMO, to be applied after initial NMO or LMO. Best scan for the two parameters to find out their effect on the data" );
  pdef->addValue( "4000", VALTYPE_NUMBER, "Offset apex [m]" );
  pdef->addValue( "5", VALTYPE_NUMBER, "Zero-offset damping (High value: high damping, 0: No damping)" );

  pdef->addParam( "horizon_nmo", "Apply horizon based NMO?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not apply horizon based NMO" );
  pdef->addOption( "yes", "Apply horizon based NMO" );
  pdef->addValue( "linear", VALTYPE_OPTION );
  pdef->addOption( "linear", "Use linear interpolation method" );
  pdef->addOption( "quadratic", "Use quadratic interpolation method" );

  pdef->addParam( "diff_nmo_offset", "Offset value for differential NMO", NUM_VALUES_FIXED, "Perform differential NMO to a non-zero offset value" );
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "Offset value (or header name) for differential NMO correction", "The output trace will be NMO corrected to this offset value, not zero offset" );

  pdef->addParam( "dump", "Dump velocity functions to standard output?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not dump velocity functions" );
  pdef->addOption( "yes", "Dump velocity functons retrieved at every location to standard output" );

  pdef->addParam( "output_vel", "Output velocity field instead of NMO'd data?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Output input data with NMO applied/removed" );
  pdef->addOption( "yes", "Output velocity field read in from input velocity table" );

  pdef->addParam( "table_old", "Velocity table (OLD format)", NUM_VALUES_FIXED, "This option is provided for backward compatibility" );
  pdef->addValue( "", VALTYPE_STRING, "Velocity table file name.",
                  "Velocity table format: The first line must contain a header title giving the names of all columns (space separated). The remaining lines give the data values. Key headers are specified in the front with a preceding '@' sign. The column specifying the time in milliseconds must be named 'time'. The last column must contain a velocity in m/s. It doesn't matter what name is chosen for the velocity column, unit is meters/second [m/s]. Example:  '@source  time  velocity'. The table does not need to contain a key trace header name; Example: 'time  velocity'." );

  pdef->addParam( "time_samp1", "Time of first sample (in [ms])", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_HEADER_NUMBER, "Time of first sample [ms]", "Only negative (and zero) times are supported" );

  //  pdef->addParam( "percent", "Percent change to be applied to input velocity", NUM_VALUES_FIXED );
  //  pdef->addValue( "0", VALTYPE_HEADER_NUMBER, "Percent [%]" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_nmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_nmo::VariableStruct* vars = reinterpret_cast<mod_nmo::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_nmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_nmo::VariableStruct* vars = reinterpret_cast<mod_nmo::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->time_sec != NULL ) {
    delete [] vars->time_sec;
    vars->time_sec = NULL;
  }
  if( vars->velocities != NULL ) {
    delete [] vars->velocities;
    vars->velocities = NULL;
  }
  if( vars->nmo != NULL ) {
    delete vars->nmo;
    vars->nmo = NULL;
  }
  if( vars->table != NULL ) {
    delete vars->table;
    vars->table = NULL;
  }
  if( vars->oldTableManager != NULL ) {
    delete vars->oldTableManager;
    vars->oldTableManager = NULL;
  }
  if( vars->hdrId_keys != NULL ) {
    delete [] vars->hdrId_keys;
    vars->hdrId_keys = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_nmo_( csParamDef* pdef ) {
  params_mod_nmo_( pdef );
}
extern "C" void _init_mod_nmo_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_nmo_( param, env, writer );
}
extern "C" bool _start_exec_mod_nmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_nmo_( env, writer );
}
extern "C" void _exec_mod_nmo_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_nmo_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_nmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_nmo_( env, writer );
}

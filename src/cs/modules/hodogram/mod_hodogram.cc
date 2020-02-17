/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"
#include "csFlexNumber.h"
#include "geolib_methods.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: HODOGRAM
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_hodogram {
  struct VariableStruct {
    int input;
    int first_sample;
    int last_sample;
    int pol_first_sample;
    int pol_last_sample;
    int method_linefit;
    int method_pol;
    int match_comp;
    float minMatchValue;
    float val_x;
    float val_y;
    float val_z;
    bool solve_polarity;
    bool pol_isMinPhase;
    bool force_origin;
    bool outputMajorAxis;
    int hdrId_vec_x;
    int hdrId_vec_y;
    int hdrId_vec_z;
    int hdrId_sensor;
    int hdrId_axis_maj;
    int hdrId_axis_med;
    int hdrId_axis_min;
    int hdrId_match_x;
    int hdrId_match_y;
    int hdrId_match_z;
    int hdrType_match_x;
    int hdrType_match_y;
    int hdrType_match_z;
  };
#define MATCH_X 1
#define MATCH_Y 2
#define MATCH_Z 3

#define LINEFIT_3D    11
#define LINEFIT_LS    12
#define LINEFIT_SVD   13
  static const int POL_PEAK_TROUGH  = 1;
  static const int POL_HYDROPHONE   = 2;
  static const int POL_MATCH_HEADER = 3;

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

#define NO_VALUE -999

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
using namespace mod_hodogram;

void polarity_correction_match( VariableStruct* vars, csTraceHeader* trcHdr, float* vec_unit );


//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_hodogram_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );


  //-----------------------------------
  vars->input        = 0;
  vars->first_sample = 0;
  vars->last_sample  = 0;
  vars->pol_first_sample = 0;
  vars->pol_last_sample  = 0;
  vars->method_linefit   = 0;
  vars->minMatchValue   = 0;
  vars->method_pol = 0;
  vars->match_comp = 0;
  vars->val_x      = 0;
  vars->val_y      = 0;
  vars->val_z      = 0;
  vars->solve_polarity  = false;
  vars->pol_isMinPhase  = false;
  vars->force_origin    = false;
  vars->outputMajorAxis = true;
  vars->hdrId_vec_x = -1;
  vars->hdrId_vec_y = -1;
  vars->hdrId_vec_z = -1;
  vars->hdrId_sensor = -1;
  vars->hdrId_axis_maj = -1;
  vars->hdrId_axis_med = -1;
  vars->hdrId_axis_min = -1;
  vars->hdrId_match_x = -1;
  vars->hdrId_match_y = -1;
  vars->hdrId_match_z = -1;
  vars->hdrType_match_x = -1;
  vars->hdrType_match_y = -1;
  vars->hdrType_match_z = -1;

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

//-----------------------------------------------------------------
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
    vars->first_sample = (int)round(startTime / shdr->sampleInt);
    vars->last_sample  = (int)round(endTime / shdr->sampleInt);
  }

  if( vars->first_sample < 0 || vars->last_sample < vars->first_sample || vars->last_sample >= shdr->numSamples) {
    writer->error("Specified start/end times are invalid. Please check input parameters.");
  }

  vars->outputMajorAxis = true;
  if( param->exists("output_axis") ) {
    param->getString( "output_axis", &text );
    if( !text.compare("major") ) {
      vars->outputMajorAxis = true;
    }
    else if( text.compare("minor") ) {
      vars->outputMajorAxis = false;
    }
    else {
      writer->error("Unknown option '%s'", text.c_str());
    }
  }

//--------------------------------------------------------------
  vars->method_pol = POL_PEAK_TROUGH;
  param->getString( "solve_pol", &text );
  if( !text.compare("yes") ) {
    vars->solve_polarity = true;
  }
  else if( !text.compare("no") ) {
    vars->solve_polarity = false;
  }
  else {
    writer->error("Unknown option '%s'", text.c_str());
  }
  if( vars->solve_polarity ) {
    startTime = -1.0;
    endTime   = -2.0;
    param->getFloat( "pol_start_time", &startTime );
    param->getFloat( "pol_end_time", &endTime );

    if( startTime == 0.0 && endTime == 0.0 ) {
      vars->pol_first_sample = 0;
      vars->pol_last_sample  = shdr->numSamples-1;
    }
    else if( startTime > endTime ) {
      writer->error("Specified start/end times for polarity analysis window are invalid. Please check input parameters.");
    }
    else {
      vars->pol_first_sample = (int)round(startTime / shdr->sampleInt);
      vars->pol_last_sample  = (int)round(endTime / shdr->sampleInt);
    }
    if( vars->pol_first_sample < 0 || vars->pol_last_sample < vars->pol_first_sample || vars->pol_last_sample > shdr->numSamples ) {
      writer->error("Specified start/end times of polarity analysis window are invalid. Please check input parameters.");
    }
    param->getString( "pol_method", &text );
    if( !text.compare("hydrophone") ) {
      vars->method_pol = POL_HYDROPHONE;
      if( vars->input != INPUT_XYZP ) {
        writer->error("The chosen polarity method '%s' requires input method 'xyzp', which means there 4 input traces: XYZ and P (hydrophone)",
          text.c_str());
      }
    }
    else if( !text.compare("min_phase") ) {
      vars->pol_isMinPhase = true;
      vars->method_pol = POL_PEAK_TROUGH;
    }
    else if( !text.compare("zero_phase") ) {
      vars->pol_isMinPhase = false;
      vars->method_pol = POL_PEAK_TROUGH;
    }
    else if( !text.compare("match_hdr") ) {
      vars->method_pol = POL_MATCH_HEADER;
    }
    else {
      writer->error("Unknown option '%s'", text.c_str());
    }

    if( vars->method_pol == POL_MATCH_HEADER ) {
      csVector<string> valueList(3);
      param->getAll( "hdr_pol_match", &valueList );
      if( valueList.size() != 4 ) {
        writer->error("Incorrect number of values for parameter HDR_POL_MATCH. Expected 4, found %d. Example: HDR_POL_MATCH  DX  DY  DZ  10.0",
          valueList.size() );
      }
      std::string match_hdr_x = valueList.at(0);
      std::string match_hdr_y = valueList.at(1);
      std::string match_hdr_z = valueList.at(2);
      std::string textNumber = valueList.at(3);
      csFlexNumber number;
      if( !number.convertToNumber( textNumber ) ) {
        writer->error("Last value for user parameter HDR_POL_MATCH not a number: '%s'", textNumber.c_str());
      }
      vars->minMatchValue = number.floatValue();

      vars->val_x = NO_VALUE;
      vars->val_y = NO_VALUE;
      vars->val_z = NO_VALUE;

      if( !hdef->headerExists(match_hdr_x.c_str()) ) {
        if( number.convertToNumber( match_hdr_x ) ) {
          if( number.doubleValue() >= 0 ) {
            vars->val_x = number.doubleValue() > 0 ? vars->minMatchValue : 0.0;
          }
          else {
            vars->val_x = -vars->minMatchValue;
          }
        }
        else {
          writer->error("Parameter HDR_POL_MATCH: Trace header %s not found.", match_hdr_x.c_str());
        }
      }
      else {
        vars->hdrId_match_x   = hdef->headerIndex(match_hdr_x.c_str());
        vars->hdrType_match_x = hdef->headerType(match_hdr_x.c_str());
      }
      if( !hdef->headerExists(match_hdr_y.c_str()) ) {
        if( number.convertToNumber( match_hdr_y ) ) {
          if( number.doubleValue() >= 0 ) {
            vars->val_y = number.doubleValue() > 0 ? vars->minMatchValue : 0.0;
          }
          else {
            vars->val_y = -vars->minMatchValue;
          }
        }
        else {
          writer->error("Option MATCH_HDR: Trace header %s not found.", match_hdr_y.c_str());
        }
      }
      else {
        vars->hdrId_match_y   = hdef->headerIndex(match_hdr_y.c_str());
        vars->hdrType_match_y = hdef->headerType(match_hdr_y.c_str());
      }
      if( !hdef->headerExists(match_hdr_z.c_str()) ) {
        if( number.convertToNumber( match_hdr_z ) ) {
          if( number.doubleValue() >= 0 ) {
            vars->val_z = number.doubleValue() > 0 ? vars->minMatchValue : 0.0;
          }
          else {
            vars->val_z = -vars->minMatchValue;
          }
        }
        else {
          writer->error("Option MATCH_HDR: Trace header %s not found.", match_hdr_z.c_str());
        }
      }
      else {
        vars->hdrId_match_z   = hdef->headerIndex(match_hdr_z.c_str());
        vars->hdrType_match_z = hdef->headerType(match_hdr_z.c_str());
      }
    }
  }
//--------------------------------------------------------------------
//
  if( param->exists("advanced_opt") ) {
    param->getString( "advanced_opt", &text );
    if( !text.compare("yes") ) {
      param->getString( "method_linefit", &text );
      if( !text.compare("3d") ) {
        vars->method_linefit = LINEFIT_3D;
      }
      else if( !text.compare("ls") ) {
        vars->method_linefit = LINEFIT_LS;
      }
      else if( !text.compare("svd") ) {
        vars->method_linefit = LINEFIT_SVD;
      }
      else {
        writer->error("Unknown option '%s'", text.c_str());
      }
      
      param->getString( "force_origin", &text );
      if( !text.compare("yes") ) {
        vars->force_origin = true;
      }
      else if( !text.compare("no") ) {
        vars->force_origin = false;
      }
      else {
        writer->error("Unknown option '%s'", text.c_str());
      }
    }
    else if( !text.compare("no") ) {  // Use default options
      vars->method_linefit = LINEFIT_3D;
      vars->force_origin   = 0;  // NO
    }
    else {
      writer->error("Unknown option '%s'", text.c_str());
    }
  }
  else {
    vars->method_linefit = LINEFIT_3D;
    vars->force_origin   = 0;  // NO
  }
//-------------------------------------------------------------------------
  // Add further necessary headers
  if( !hdef->headerExists("axis_maj") ) {
    vars->hdrId_axis_maj = hdef->addHeader( TYPE_FLOAT, "axis_maj", "Size of major axis of best fit spheroid" );
  }
  if( !hdef->headerExists("axis_med") ) {
    vars->hdrId_axis_med = hdef->addHeader( TYPE_FLOAT, "axis_med", "Size of medium (=first minor) axis of best fit spheroid" );
  }
  if( !hdef->headerExists("axis_min") ) {
    vars->hdrId_axis_min = hdef->addHeader( TYPE_FLOAT, "axis_min", "Size of minor (smallest) axis of best fit spheroid" );
  }

  if( hdef->headerType( "sensor" ) != TYPE_INT ) {
    writer->error("Trace header 'sensor' has the wrong type.");
  }
  vars->hdrId_sensor  = hdef->headerIndex( "sensor" );



  if( hdef->headerExists ("vec_x") ) {
    checkHeader( "vec_x", hdef, writer );
    vars->hdrId_vec_x = hdef->headerIndex("vec_x" );
  }
  else {
    vars->hdrId_vec_x = hdef->addHeader( TYPE_FLOAT, "vec_x", "Unit direction vector - X component" );
  }
  if( hdef->headerExists ("vec_y") ) {
    checkHeader( "vec_y", hdef, writer );
    vars->hdrId_vec_y = hdef->headerIndex("vec_y" );
  }
  else {
    vars->hdrId_vec_y = hdef->addHeader( TYPE_FLOAT, "vec_y", "Unit direction vector - Y component" );
  }
  if( hdef->headerExists ("vec_z") ) {
    checkHeader( "vec_z", hdef, writer );
    vars->hdrId_vec_z = hdef->headerIndex("vec_z" );
  }
  else {
    vars->hdrId_vec_z = hdef->addHeader( TYPE_FLOAT, "vec_z", "Unit direction vector - Z component" );
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_hodogram_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;


  float vec_unit[3];
  float vec_axes[3];
  float vec_minor[3];

//  if(  != vars->nTraces ) {
//    writer->error("Wrong number of traces in ensemble.\nNumber of traces found: %d.\nNumber of traces expected: %d.\nUse 'ENSEMBLE REDEFINE' or similar tool to set correct ensemble sorting.",
//      traceGather->numTraces(), vars->nTraces);
 //   return;
 // }

  int  trace_index[4];
  int nTraces = traceGather->numTraces();
  for( int i = 0; i < 4; i++ ) {
    trace_index[i] = NO_VALUE;
  }

  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    int sensor = traceGather->trace(itrc)->getTraceHeader()->intValue(vars->hdrId_sensor);
    if( sensor == SENSOR_INDEX_Z2 ) sensor = SENSOR_INDEX_Z5;  // Convert sensor index 2 to 5 (both are Z component)
    if( sensor >= SENSOR_INDEX[ID_X] && sensor <= SENSOR_INDEX[ID_Z] ) {
      if( trace_index[sensor-3] == NO_VALUE ) {
        trace_index[sensor-3] = itrc;
      }
      else {
        writer->error("Input gather contains more than one trace for sensor %d.", sensor);
        return;
      }
    }
    else if( sensor == SENSOR_INDEX[ID_P] ) {
      if( trace_index[3] == NO_VALUE ) {
        trace_index[3] = itrc;
      }
      else {
        writer->error("Input gather contains more than one trace for sensor %d.", sensor);
        return;
      }
    }
  }

//---------------------------------------------------
  csTraceHeader* trcHdr[3];
  for( int idSensor = 0; idSensor < 3; idSensor++ ) {
    if( trace_index[idSensor] == NO_VALUE ) {
      writer->error("Input gather is missing a sensor %d trace.", SENSOR_INDEX[idSensor]);
    }
    else {
      trcHdr[idSensor] = traceGather->trace(trace_index[idSensor])->getTraceHeader();
    }
  }
  if( vars->method_pol == POL_HYDROPHONE && trace_index[3] == NO_VALUE ) {
    writer->error("Ensemble is missing a sensor 1 trace (= hydrophone).");
    return;
  }

//---------------------------------------------------
// Compute unit vector in specified data window
//

  bool ret  = false;
  if( vars->method_linefit == LINEFIT_3D ) {
    ret = linefit_3d_all(
              traceGather->trace(trace_index[ID_X])->getTraceSamples(),
              traceGather->trace(trace_index[ID_Y])->getTraceSamples(),
              traceGather->trace(trace_index[ID_Z])->getTraceSamples(),
              vars->first_sample,
              vars->last_sample,
              vars->force_origin,
              vec_unit,
              vec_axes,
              vec_minor );
  }
  else if( vars->method_linefit == LINEFIT_LS || vars->method_linefit == LINEFIT_SVD ) {
    ret = linefit_3d_2step(
              traceGather->trace(trace_index[ID_X])->getTraceSamples() + vars->first_sample,
              traceGather->trace(trace_index[ID_Y])->getTraceSamples() + vars->first_sample,
              traceGather->trace(trace_index[ID_Z])->getTraceSamples() + vars->first_sample,
              vars->last_sample - vars->first_sample + 1,
              vars->method_linefit,
              vars->force_origin,
              vec_unit );
  }

  if( !ret ) {
    vec_unit[0] = 0.0;
    vec_unit[1] = 0.0;
    vec_unit[2] = 1.0;
    vec_minor[0] = 0.0;
    vec_minor[1] = 0.0;
    vec_minor[2] = 1.0;
    writer->warning("Hodogram analysis failed.. Bad/noisy/zero input data?");
  }

  if( edef->isDebug() ) {
    writer->line("Unit vector = %12.6f %12.6f %12.6f\n", vec_unit[0], vec_unit[1], vec_unit[2] );
  }

//
// Adjust polarity of unit vector
//

  if( vars->solve_polarity ) {
    switch( vars->method_pol ) {
      case POL_PEAK_TROUGH:
        polarity_correction(
          traceGather->trace(trace_index[ID_X])->getTraceSamples(),
          traceGather->trace(trace_index[ID_Y])->getTraceSamples(),
          traceGather->trace(trace_index[ID_Z])->getTraceSamples(),
          vars->pol_first_sample,
          vars->pol_last_sample,
          vars->pol_isMinPhase,
          vec_unit );
        break;
      case POL_HYDROPHONE:
        polarity_correction_p(
          traceGather->trace(trace_index[ID_X])->getTraceSamples(),
          traceGather->trace(trace_index[ID_Y])->getTraceSamples(),
          traceGather->trace(trace_index[ID_Z])->getTraceSamples(),
          traceGather->trace(trace_index[ID_P])->getTraceSamples(),
          vars->pol_first_sample,
          vars->pol_last_sample,
          vec_unit );
        break;
      case POL_MATCH_HEADER:
        {
          polarity_correction_match( vars, trcHdr[ID_X], vec_unit );
        }
        break;
    }
  }

//
// Set output headers
//
  if( vars->outputMajorAxis ) {
    for( int idSensor = 0; idSensor < 3; idSensor++ ) {
      trcHdr[idSensor]->setDoubleValue( vars->hdrId_vec_x, vec_unit[0] );
      trcHdr[idSensor]->setDoubleValue( vars->hdrId_vec_y, vec_unit[1] );
      trcHdr[idSensor]->setDoubleValue( vars->hdrId_vec_z, vec_unit[2] );
    }
  }
  else {
    for( int idSensor = 0; idSensor < 3; idSensor++ ) {
      trcHdr[idSensor]->setDoubleValue( vars->hdrId_vec_x, vec_minor[0] );
      trcHdr[idSensor]->setDoubleValue( vars->hdrId_vec_y, vec_minor[1] );
      trcHdr[idSensor]->setDoubleValue( vars->hdrId_vec_z, vec_minor[2] );
    }
  }
  for( int idSensor = 0; idSensor < 3; idSensor++ ) {
    trcHdr[idSensor]->setDoubleValue( vars->hdrId_axis_maj, vec_axes[0] );
    trcHdr[idSensor]->setDoubleValue( vars->hdrId_axis_med, vec_axes[1] );
    trcHdr[idSensor]->setDoubleValue( vars->hdrId_axis_min, vec_axes[2] );
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_hodogram_( csParamDef* pdef ) {
  pdef->setModule( "HODOGRAM", "Hodogram analysis", "Compute 3D polarisation vector from hodogram." );

  pdef->addParam( "input", "Input data", NUM_VALUES_FIXED, "Specifies how many input traces shall be processed at once" );
  pdef->addValue( "xyz", VALTYPE_OPTION );
  pdef->addOption( "xyz", "3 sequential input traces, containing sensors 3(X), 4(Y) and 5(Z). Sorting of these three may be arbitrary." );
  pdef->addOption( "xyzp", "4 sequential input traces, containing sensors 3(X), 4(Y) and 5(Z), and 4th trace with sensor 1(hydrophone). Sorting of these three may be arbitrary." );
  pdef->addOption( "ensemble", "Input whole ensemble (user sorted). Each input ensemble must contain one trace each for sensor 3(X), 4(Y) and 5(Z). Sorting of these three may be arbitrary." );

  pdef->addParam( "start_time", "Start time of analysis window", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Start time [ms]" );
  pdef->addParam( "end_time", "End time of analysis window", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "End time [ms]" );

  pdef->addParam( "output_axis", "Output axis", NUM_VALUES_FIXED );
  pdef->addValue( "major", VALTYPE_OPTION );
  pdef->addOption( "major", "Output 3d vector along major axis" );
  pdef->addOption( "minor", "Output 3d vector along minor axis" );

  pdef->addParam( "solve_pol", "Solve for output vector polarity", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Yes" );
  pdef->addOption( "no", "No", "In this case, the polarity of the output vector is unchanged" );

  pdef->addParam( "pol_method", "Method to establish polarity of polarisation vector", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "min_phase", "Determine polarity from the first peak and trough (minimum phase wavelet).",
    "This specified time window should include the first two lobes of the wavelet, of which the first one determines the polarity. The first lobes (peak/trough or trough/peak) are determined as the minimum and maximum values in the specified time window." );
  pdef->addOption( "zero_phase", "Determine polarity from the largest peak or trough in the specified time window.", "This method works best for a zero phase wavelet." );
  pdef->addOption( "hydrophone", "Polarity is established by cross-correlation with hydrophone trace.", "Requires input of 4th trace = hydrophone" );
  pdef->addOption( "match_hdr", "Match polarity of polarisation vector with the values in the specified trace headers. Use for example dx, dy, dz derived from the geometry to guide the polarity estimation of the polarisation vector.", "Specify trace headers, or values, in parameter 'hdr_pol_match'" );

  pdef->addParam( "pol_start_time", "Start time of polarity analysis window", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Polarity start time [ms]" );
  pdef->addParam( "pol_end_time", "End time of polarity analysis window", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Polarity end time [ms]" );

  pdef->addParam( "hdr_pol_match", "Header names (or fixed values) for polarity method 'match_hdr'", NUM_VALUES_FIXED );
  pdef->addValue( "<hdr_x>", VALTYPE_STRING, "Header giving approximate polarity of X polarisaton" );
  pdef->addValue( "<hdr_y>", VALTYPE_STRING, "Header giving approximate polarity of Y polarisation" );
  pdef->addValue( "<hdr_z>", VALTYPE_STRING, "Header giving approximate polarity of Z polarisation" );
  pdef->addValue( "0", VALTYPE_STRING, "Minimum value to match. If all of the three input headers/values are lower than this value, no attempt will be made to change the polarity" );

  pdef->addParam( "advanced_opt", "Specify advanced options.", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Use advanced options" );
  pdef->addOption( "no", "Do not use advanced options." );

  pdef->addParam( "force_origin", "Force vector through origin.", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Force output vector through origin" );
  pdef->addOption( "no", "Do not force output vector through origin." );

  pdef->addParam( "method_linefit", "Method for line fitting.", NUM_VALUES_FIXED );
  pdef->addValue( "3d", VALTYPE_OPTION );
  pdef->addOption( "3d", "Perform full 3d line fitting through data points" );
}

void polarity_correction_match( VariableStruct* vars, csTraceHeader* trcHdr, float* vec_unit ) {
  float sign = 0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  if( vars->val_x != NO_VALUE ) {
    dx = vars->val_x;
  }
  else if( vars->hdrType_match_x == TYPE_FLOAT || vars->hdrType_match_x == TYPE_DOUBLE ) {
    dx = trcHdr->doubleValue( vars->hdrId_match_x );
  }
  else if( vars->hdrType_match_x == TYPE_INT ) {
    dx = (double)trcHdr->intValue( vars->hdrId_match_x );
  }

  if( vars->val_y != NO_VALUE ) {
    dy = vars->val_y;
  }
  else if( vars->hdrType_match_y == TYPE_FLOAT || vars->hdrType_match_y == TYPE_DOUBLE ) {
    dy = trcHdr->doubleValue( vars->hdrId_match_y );
  }
  else if( vars->hdrType_match_y == TYPE_INT ) {
    dy = (double)trcHdr->intValue( vars->hdrId_match_y );
  }

  if( vars->val_z != NO_VALUE ) {
    dz = vars->val_z;
  }
  else if( vars->hdrType_match_z == TYPE_FLOAT || vars->hdrType_match_z == TYPE_DOUBLE ) {
    dz = trcHdr->doubleValue( vars->hdrId_match_z );
  }
  else if( vars->hdrType_match_z == TYPE_INT ) {
    dz = (double)trcHdr->intValue( vars->hdrId_match_z );
  }

  double dx_abs = abs(dx);
  double dy_abs = abs(dy);
  double dz_abs = abs(dz);

  if( vars->minMatchValue > 0 ) {
    if( dx_abs < vars->minMatchValue && dy_abs < vars->minMatchValue && dz_abs < vars->minMatchValue ) {
      return;
    }
  }

  if( dx_abs > dy_abs && dx_abs > dz_abs ) {
    sign = (dx >= 0) ? 1 : -1;
    sign *= vec_unit[0];
  }
  else if( dy_abs > dz_abs ) {
    sign = (dy >= 0) ? 1 : -1;
    sign *= vec_unit[1];
  }
  else {
    sign = (dz >= 0) ? 1 : -1;
    sign *= vec_unit[2];
  }
  if( sign < 0) {
    vec_unit[0] *= -1.0;
    vec_unit[1] *= -1.0;
    vec_unit[2] *= -1.0;
  }
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_hodogram_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_hodogram::VariableStruct* vars = reinterpret_cast<mod_hodogram::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_hodogram_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_hodogram::VariableStruct* vars = reinterpret_cast<mod_hodogram::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_hodogram_( csParamDef* pdef ) {
  params_mod_hodogram_( pdef );
}
extern "C" void _init_mod_hodogram_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_hodogram_( param, env, writer );
}
extern "C" bool _start_exec_mod_hodogram_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_hodogram_( env, writer );
}
extern "C" void _exec_mod_hodogram_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_hodogram_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_hodogram_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_hodogram_( env, writer );
}

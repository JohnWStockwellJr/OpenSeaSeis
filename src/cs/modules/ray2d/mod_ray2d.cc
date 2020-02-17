/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"
#include "geolib_math.h"
#include "geolib_methods.h"
#include "csMatrixFStyle.h"
#include "wfront_sub.h"
#include "csTimeStretch.h"
#include "csSort.h"
#include "csSortManager.h"
#include "csFileUtils.h"
#include <cmath>
#include <cstring>
#include "csSegyWriter.h"
#include "csSegyDefines.h"
#include "csSegyHdrMap.h"
#include "csSegyBinHeader.h"
#include "csSegyTraceHeader.h"

using namespace cseis_system;
using namespace cseis_geolib;using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: RAY2D
 *
 * @author Bjorn Olofsson
 * @date   2011
 */
namespace mod_ray2d {
  struct VariableStruct {
    int hdrId_rec_x;
    int hdrId_rec_z;
    int hdrId_sou_x;
    int hdrId_sou_z;
    int* hdrId_rayTime;
    int* hdrId_rayAngle;
    int* hdrId_ccpx;
    int* hdrId_ccpz;
    float rec_x;
    float rec_z;
    float sou_x;
    float sou_z;
    int gatherType;
    int nmoMethod;
    bool applyBinning;
    int modelType;

    // Model definition
    int numInterfaces; // Number of model interfaces
    int numInterfacesAll; // Number of model interfaces & boreholes
    int* numPointsInt; // Number of points defining each interface
    int maxPointsInt;  // Maximum number of points defining layer interfaces = max(numPointsInt)
    csMatrixFStyle<int> iii;
    csMatrixFStyle<float> x_int;
    csMatrixFStyle<float> z_int;
    csMatrixFStyle<int> iii2;
    csMatrixFStyle<float> b_int;
    csMatrixFStyle<float> c_int;
    csMatrixFStyle<float> d_int;
    float left_border;
    float right_border;

    float* rho1;      // Density
    float* rho2;      // Density linear coefficient
    int* nx_grid;
    int* nz_grid;
    csMatrixFStyle<float> x_grid;
    csMatrixFStyle<float> z_grid;
    csMatrixFStyle<float> v_grid;
    int maxPointsXGrid;
    int maxPointsZGrid;
    float* vpvs_ratio;

    int rec_interface;
    int rec_layer;
    int sou_interface;
    int sou_layer;
    float sou_dzdx;

    int numCodes;
    int maxNumCodeSteps;
    int* numCodeSteps;
    int* codeStepPlot;
    csMatrixFStyle<int> rayCode;
    csMatrixFStyle<int> rayHold;

    float dtray;  // Time interval between ray-tracing steps [s]
    float spread_dist; // Maximum spreading
    float spread_angle_rad;  // Maximum angle between ray before new ray is constructed
    int numRays;
    float anglein_rad;
    float angleout_rad;

    int seismogramType;
    float rickerPeakFrequency;

    int flag_stack;
    int flag_twopoint;
    int flag_compounds;
    int f_rayouttmp;
    int f_rayout;
    int f_xrays;
    int f_wfronts;
    int step_wfronts;
    int f_timeout;
    FILE* file_timeout;

    std::string filename_timeout;
    std::string filename_intout;
    bool isFirstCall;
 };
  static int const MODEL_1D        = 1;
  static int const MODEL_2D_SIMPLE = 2;
  static int const MODEL_2D_GRID   = 3;
  static int const RECEIVER_GATHER = 11;
  static int const SOURCE_GATHER   = 12;
  static int const SHEAR_PARAM_VPVS = 21;
  static int const SHEAR_PARAM_VS   = 22;

  static int const NMO_NONE  = 31;
  static int const NMO_TIME  = 32;
  static int const NMO_DEPTH = 33;

  static int const SEISMOGRAM_NONE      = 44;
  static int const SEISMOGRAM_DYNAMIC   = 45;
  static int const SEISMOGRAM_KINEMATIC = 46;

  static int const VELOUT_VP      = 51;
  static int const VELOUT_VS      = 52;
  static int const VELOUT_VPVS    = 53;
  static int const VELOUT_RHO     = 54;
  static int const VELOUT_INT     = 55;
}

void compute_z( mod_ray2d::VariableStruct* vars, float xvalue, int intIndex, float* zout, int* pointIndex_out );
void exportVelocityModel( mod_ray2d::VariableStruct* vars, cseis_geolib::csSegyWriter* segyWriter, double dx, double dz, float xMin, float xMax, float zMin, float zMax, int velout_flag );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_ray2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  mod_ray2d::VariableStruct* vars = new mod_ray2d::VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->hdrId_rec_x = -1;
  vars->hdrId_rec_z = -1;
  vars->hdrId_sou_x = -1;
  vars->hdrId_sou_z = -1;
  vars->rec_x = 0;
  vars->rec_z = 0;
  vars->sou_x = 0;
  vars->sou_z = 0;
  vars->gatherType = 0;
  vars->nmoMethod = mod_ray2d::NMO_NONE;
  vars->applyBinning = false;

  vars->numInterfaces = 0;
  vars->numInterfacesAll = 0;
  vars->numPointsInt  = NULL;
  vars->maxPointsInt  = 0;
  vars->left_border = 0;
  vars->right_border = 0;

  vars->rho1 = NULL;
  vars->rho2 = NULL;
  vars->nx_grid = NULL;
  vars->nz_grid = NULL;
  //vars->x_grid;
  //vars->z_grid;
  vars->maxPointsXGrid = 0;
  vars->maxPointsZGrid = 0;
  vars->vpvs_ratio = NULL;

  // vars->rec_interface    = 2;
  // vars->rec_layer        = 2;
  // vars->sou_interface    = 0;
  // vars->sou_layer        = 1;
  // CHANGE: Set layer/interface more intelligently especially for 1D models
  // The following setup works for OBN receiver gather geometry. It swaps receivers and sources since the ray tracer assumes shot geometry.
  vars->rec_interface    = 1;
  vars->rec_layer        = 1;
  vars->sou_interface    = 0;
  vars->sou_layer        = 2;

  vars->sou_dzdx = 0.0;

  vars->numCodes = 0;
  vars->maxNumCodeSteps = 10;
  vars->numCodeSteps = NULL;
  vars->codeStepPlot = NULL;

  // Set up spreading & wavefront construction criteria
  vars->dtray      = 0.04f;  // Time interval between ray-tracing steps [s]
  vars->spread_dist = 0.06f; // Maximum spreading
  vars->spread_angle_rad  = 2.3f;  // Maximum angle between ray before new ray is constructed. Convert from [deg] to [rad] later on
  vars->numRays = 400;
  vars->anglein_rad  = (float)(0.0   * M_PI / 180.0);
  vars->angleout_rad = (float)(180.0 * M_PI / 180.0);

  vars->seismogramType = mod_ray2d::SEISMOGRAM_NONE;
  vars->rickerPeakFrequency = 30.0;

  vars->flag_stack = 0;
  vars->flag_twopoint = 0;
  vars->flag_compounds = 0;
  vars->f_rayouttmp = 0;
  vars->f_rayout   = 0;
  vars->f_xrays    = 0;
  vars->f_wfronts  = 0;
  vars->step_wfronts = 4;
  vars->f_timeout = 0;
  vars->file_timeout = NULL;
  vars->filename_timeout = "";
  vars->filename_intout = "";
  vars->isFirstCall = true;

  vars->modelType = -1;


  float EPSILON = 0.00001f;


  //********************************************************************************
  //
  // Read in misc parameters
  //

  std::string text;
  int flag_duplicate_raycodes = 0;
  if( param->exists( "duplicate_raycodes" ) ) {
    param->getString("duplicate_raycodes", &text);
    if( !text.compare("no") ) {
      flag_duplicate_raycodes = 0;
    }
    else if( !text.compare("yes") ) {
      flag_duplicate_raycodes = 1;
    }
    else  {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  if( param->exists( "flags" ) ) {
    int numFlags = param->getNumValues("flags");
    if( numFlags > 0 ) param->getInt( "flags", &vars->flag_twopoint, 0 );
    if( numFlags > 1 ) param->getInt( "flags", &vars->flag_stack, 1 );
    if( numFlags > 2 ) param->getInt( "flags", &vars->flag_compounds, 2 );
  }
  if( param->exists( "dump_rays" ) ) {
    param->getInt( "dump_rays", &vars->f_rayout, 0 );
    param->getInt( "dump_rays", &vars->f_rayouttmp, 1 );
    param->getInt( "dump_rays", &vars->f_xrays, 2 );
  }
  if( param->exists( "dump_wfronts" ) ) {
    param->getInt( "dump_wfronts", &vars->f_wfronts, 0 );
    param->getInt( "dump_wfronts", &vars->step_wfronts, 1 );
  }
  if( param->exists( "dump_times" ) ) {
    param->getString( "dump_times", &vars->filename_timeout );
    if( !csFileUtils::createDoNotOverwrite( vars->filename_timeout ) ) {
      writer->error("Error opening output file %s\n", vars->filename_timeout.c_str() );
    }
  }

  if( param->exists( "dt_ray" ) ) {
    param->getFloat( "dt_ray", &vars->dtray );
    vars->dtray /= 1000.0f; // Convert from [ms] to [s]
  }
  if( param->exists( "spread_dist" ) ) {
    param->getFloat( "spread_dist", &vars->spread_dist );
    vars->spread_dist /= 1000.0f; // Convert from [m] to [km]
  }
  if( param->exists( "spread_angle" ) ) {
    param->getFloat( "spread_angle", &vars->spread_angle_rad );
  }
  vars->spread_angle_rad *= (float)( M_PI / 180.0 );  // Convert from deg to rad

  if( param->exists( "nrays" ) ) {
    param->getInt( "nrays", &vars->numRays );
  }
  if( param->exists( "beam" ) ) {
    param->getFloat( "beam", &vars->anglein_rad, 0 );
    param->getFloat( "beam", &vars->angleout_rad, 1 );
    vars->anglein_rad  *= (float)(M_PI / 180.0);
    vars->angleout_rad *= (float)(M_PI / 180.0);
  }

  //********************************************************************************
  //
  // Read in gather type
  //

  if( param->exists( "gather" ) ) {
    param->getString( "gather", &text );
    if( !text.compare("receiver") ) {
      vars->gatherType = mod_ray2d::RECEIVER_GATHER;
    }
    else if( !text.compare("source") ) {
      vars->gatherType = mod_ray2d::SOURCE_GATHER;
    }
    else  {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  //********************************************************************************
  //********************************************************************************
  //
  // Read in source/receiver coordinates
  //
  csFlexNumber number;

  param->getString( "rec_x", &text );
  if( !number.convertToNumber( text ) ) {
    if( !hdef->headerExists(text) ) {
      writer->error("Specified trace header name containing receiver X coordinate does not exist: '%s'", text.c_str());
    }
    vars->hdrId_rec_x = hdef->headerIndex(text);
  }
  else { // User specified a constant value
    vars->rec_x = number.floatValue();
    if( vars->gatherType == mod_ray2d::SOURCE_GATHER ) {
      writer->error("For source gathers, the receiver X coordinate cannot be constant (%f). Please specify a trace header name containing the receiver X coordinate", vars->rec_x);
    }
    writer->line("Constant receiver X coordinate = %.2f", vars->rec_x);
  }

  param->getString( "rec_z", &text );
  if( !number.convertToNumber( text ) ) {
    if( !hdef->headerExists(text) ) {
      writer->error("Specified trace header name containing receiver Z coordinate does not exist: '%s'", text.c_str());
    }
    vars->hdrId_rec_z = hdef->headerIndex(text);
  }
  else { // User specified a constant value
    vars->rec_z = number.floatValue();
    writer->line("Constant receiver Z coordinate = %.2f", vars->rec_z);
  }

  param->getString( "sou_x", &text );
  if( !number.convertToNumber( text ) ) {
    if( !hdef->headerExists(text) ) {
      writer->error("Specified trace header name containing source X coordinate does not exist: '%s'", text.c_str());
    }
    vars->hdrId_sou_x = hdef->headerIndex(text);
  }
  else { // User specified a constant value
    vars->sou_x = number.floatValue();
    if( vars->gatherType == mod_ray2d::RECEIVER_GATHER ) {
      writer->error("For receiver gathers, the source X coordinate cannot be constant (%f). Please specify a trace header name containing the source X coordinate", vars->sou_x);
    }
    writer->line("Constant source X coordinate = %.2f", vars->sou_x);
  }

  param->getString( "sou_z", &text );
  if( !number.convertToNumber( text ) ) {
    if( !hdef->headerExists(text) ) {
      writer->error("Specified trace header name containing source Z coordinate does not exist: '%s'", text.c_str());
    }
    vars->hdrId_sou_z = hdef->headerIndex(text);
  }
  else { // User specified a constant value
    vars->sou_z = number.floatValue();
    writer->line("Constant source Z coordinate = %.2f", vars->sou_z);
  }

  vars->rec_x /= 1000.0f; // Convert from [m] to [km]
  vars->rec_z /= 1000.0f; // Convert from [m] to [km]
  vars->sou_x /= 1000.0f; // Convert from [m] to [km]
  vars->sou_z /= 1000.0f; // Convert from [m] to [km]

  //********************************************************************************
  //********************************************************************************
  //
  // Read in other parameters
  //
  if( param->exists( "apply_nmo" ) ) {
    param->getString( "apply_nmo", &text );
    if( !text.compare("no") ) {
      vars->nmoMethod = mod_ray2d::NMO_NONE;
    }
    else if( !text.compare("time") ) {
      vars->nmoMethod = mod_ray2d::NMO_TIME;
    }
    else if( !text.compare("depth") ) {
      vars->nmoMethod = mod_ray2d::NMO_DEPTH;
      shdr->domain = cseis_geolib::DOMAIN_XD;
    }
    else  {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  if( param->exists( "apply_binning" ) ) {
    param->getString( "apply_binning", &text );
    if( !text.compare("yes") ) {
      vars->applyBinning = true;
    }
    else if( !text.compare("no") ) {
      vars->applyBinning = false;
    }
    else  {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  if( param->exists( "seismogram" ) ) {
    param->getString( "seismogram", &text, 1 );
    if( !text.compare("dynamic") ) {
      vars->seismogramType = mod_ray2d::SEISMOGRAM_DYNAMIC;
    }
    else if( !text.compare("kinematic") ) {
      vars->seismogramType = mod_ray2d::SEISMOGRAM_KINEMATIC;
    }
    else  {
      writer->error("Unknown option: %s", text.c_str() );
    }
    param->getFloat( "seismogram", &vars->rickerPeakFrequency, 0 );
  }

  //********************************************************************************
  //********************************************************************************
  //
  // Read in source/receiver location
  //

  int layer;
  int interface = 0;
  param->getInt( "rec_loc", &layer, 0 );      
  if( param->getNumValues("rec_loc") > 1 ) param->getInt( "rec_loc", &interface, 1 );
  if( vars->gatherType == mod_ray2d::RECEIVER_GATHER ) {
    vars->sou_layer     = layer;
    vars->sou_interface = interface;
  }
  else {
    vars->rec_layer     = layer;
    vars->rec_interface = interface;
  }

  interface = 0;
  param->getInt( "sou_loc", &layer, 0 );      
  if( param->getNumValues("sou_loc") > 1 ) param->getInt( "sou_loc", &interface, 1 );
  if( vars->gatherType == mod_ray2d::RECEIVER_GATHER ) {
    vars->rec_layer     = layer;
    vars->rec_interface = interface;
  }
  else {
    vars->sou_layer     = layer;
    vars->sou_interface = interface;
  }

  //********************************************************************************
  //********************************************************************************
  //
  // Read in model
  //
  int numBoreholes = 0; // Number of boreholes (not supported)

  int shearParam = mod_ray2d::SHEAR_PARAM_VPVS;
  float maxOffset1DModel = 10.000; // [km]
  if( param->exists("max_offset") ) {
    param->getFloat("max_offset", &maxOffset1DModel);
    maxOffset1DModel /= 1000.0f; // Convert to [km]
  }


  if( param->exists("model_type") ) {
    param->getString("model_type", &text, 0);
    if( !text.compare("1d") ) {
      vars->modelType = mod_ray2d::MODEL_1D;
    }
    else if( !text.compare("2d") ) {
      vars->modelType = mod_ray2d::MODEL_2D_SIMPLE;
    }
    else if( !text.compare("2d_grid") ) {
      vars->modelType = mod_ray2d::MODEL_2D_GRID;
      writer->error("Model type not supported yet: 2D grid");
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
    param->getString("model_type", &text, 1);
    if( !text.compare("vpvs") ) {
      shearParam = mod_ray2d::SHEAR_PARAM_VPVS;
    }
    else if( !text.compare("vs") ) {
      shearParam = mod_ray2d::SHEAR_PARAM_VS;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  //  else {
  //    writer->error("No 1D model defined (user parameter 'model1d'), and no 2D model defined (user parameters 'int2d' and 'layer')");
  //  }

  int numLines = 0;
  if( vars->modelType == mod_ray2d::MODEL_1D ) {
    if( !param->exists("int1d") ) {
      writer->error("No 1D model interface defined (user parameter 'int1d')");
    }
    numLines = param->getNumLines("int1d");
    vars->numInterfaces = numLines;
  }
  else if( vars->modelType == mod_ray2d::MODEL_2D_SIMPLE ) {
    if( !param->exists("int2d") ) {
      writer->error("No 2D model interface defined (user parameter 'int2d')");
    }
    else if( !param->exists("layer") ) {
      writer->error("No local 1D layer properties defined for 2D model (user parameter 'layer')");
    }
    numLines = param->getNumLines("int2d");
    param->getIntAtLine( "int2d", &vars->numInterfaces, numLines-1 );
  }
  else { // 2D grid model
    writer->error("Error");
  }
  if( vars->numInterfaces < 2 ) {
    writer->error("Please specify at least two model interfaces. Only %d has been defined", vars->numInterfaces);
  }

  vars->numInterfacesAll = vars->numInterfaces + numBoreholes;  // Number of all model "interfaces"
  vars->numPointsInt = new int[vars->numInterfacesAll];

  csVector<float> xList;
  csVector<float> zList;
  csVector<int> iiiList;

  if( vars->modelType == mod_ray2d::MODEL_1D ) {
    for( int inter = 0; inter < vars->numInterfaces; inter++ ) {
      vars->numPointsInt[inter] = 2;
    }
    float zval;
    int id;
    for( int iline = 0; iline < numLines; iline++ ) {
      param->getIntAtLine("int1d", &id, iline, 0);
      param->getFloatAtLine("int1d", &zval, iline, 1);
      zList.insertEnd(zval/1000.0f);
      zList.insertEnd(zval/1000.0f);
      xList.insertEnd(-maxOffset1DModel);
      xList.insertEnd(maxOffset1DModel);
      iiiList.insertEnd(-1);
      iiiList.insertEnd(-1);
    }
  }
  else { // 2D model
    int interCurrent = 1;
    int pointCounter = 0;
    for( int iline = 0; iline < numLines; iline++ ) {
      int inter;
      param->getIntAtLine("int2d", &inter, iline, 0);
      if( inter < 1 ) {
        writer->error("Inconsistent 2D model definition: Interface number = %d (<1)", inter);
      }
      else if( inter < interCurrent ) {
        writer->error("Inconsistent 2D model definition: Point #%d in model specified for previous interface (#%d)", iline, inter);
      }
      else if( inter > interCurrent+1 ) {
      writer->error("Inconsistent 2D model definition: Interface #%d follows interface #%d. Model interface numbers must be consecutive.", interCurrent, inter);
      }
      if( inter != interCurrent ) {
        vars->numPointsInt[interCurrent-1] = pointCounter;
        pointCounter = 0;
        interCurrent = inter;
      }
      float xval;
      float zval;
      int iii = 0;
      param->getFloatAtLine("int2d", &xval, iline, 1);
      param->getFloatAtLine("int2d", &zval, iline, 2);
      if( param->getNumValues("int2d", iline ) > 3 ) {
        param->getIntAtLine("int2d", &iii, iline, 3);
      }
      xList.insertEnd(xval/1000.0f);
      zList.insertEnd(zval/1000.0f);
      iiiList.insertEnd(iii);
      pointCounter += 1;
    }
    vars->numPointsInt[interCurrent-1] = pointCounter;
    int numCheckLayers = param->getNumLines("layer");
    if( numCheckLayers < vars->numInterfaces-1 ) {
      writer->error("Too few layers defined for 2D model, user parameter 'layer'. Specified layers: %d, actual number of layers: %d (=#interfaces-1)",
                 numCheckLayers, vars->numInterfaces-1);
    }
    else if( numCheckLayers > vars->numInterfaces-1 ) {
      writer->error("Too many layers defined for 2D model, user parameter 'layer'. Specified layers: %d, actual number of layers: %d (=#interfaces-1)",
                 numCheckLayers, vars->numInterfaces-1);
    }
  } // END: Read in 2D model

  //----------------------------------------------------------------------

  for( int i = 0; i < vars->numInterfaces; i++ ) {
    if( vars->numPointsInt[i] > vars->maxPointsInt ) vars->maxPointsInt = vars->numPointsInt[i];
  }
  vars->iii.resize(vars->maxPointsInt,vars->numInterfacesAll);
  vars->x_int.resize(vars->maxPointsInt,vars->numInterfacesAll);
  vars->z_int.resize(vars->maxPointsInt,vars->numInterfacesAll);

  int counter = 0;
  float xMin = xList.at(0);
  float xMax = xMin;
  float zMin = zList.at(0);
  float zMax = zMin;

  writer->line("Number of interfaces: %d, max number of points per interface: %d\n", vars->numInterfaces, vars->maxPointsInt);
  for( int ip = 0; ip < xList.size(); ip++ ) {
    float xtmp = xList.at(ip);
    float ztmp = zList.at(ip);
    if( xtmp < xMin ) xMin = xtmp;
    if( xtmp > xMax ) xMax = xtmp;
    if( ztmp < zMin ) zMin = ztmp;
    if( ztmp > zMax ) zMax = ztmp;
    writer->line("Model point #%-2d :   %12.3f %12.3f", ip, xList.at(ip), zList.at(ip));
  }

  for( int kInt = 0; kInt < vars->numInterfaces; kInt++ ) {
    for( int ip = 0; ip < vars->numPointsInt[kInt]; ip++ ) {
      vars->x_int(ip,kInt) = xList.at(counter);
      vars->z_int(ip,kInt) = zList.at(counter);
      vars->iii(ip,kInt)   = iiiList.at(counter);
      //      fprintf(stdout,"%d %d %d\n", ip, kInt, vars->iii(ip,kInt));
      counter += 1;
    }
    if( vars->iii(0,kInt) != -1 ) {
      vars->iii(0,kInt) = -1;
    }
    if( vars->iii(vars->numPointsInt[kInt]-1,kInt) != -1 ) {
      vars->iii(vars->numPointsInt[kInt]-1,kInt) = -1;
    }
    if( vars->x_int(0,kInt) != xMin ) {
      writer->error("First point X coordinate (=%f) of interface #%d not at left model boundary (=%f)",
                 vars->x_int(0,kInt), kInt+1, xMin);
    }
    else if( vars->x_int(vars->numPointsInt[kInt]-1,kInt) != xMax ) {
      writer->error("Last point X coordinate (=%f) of interface #%d not at right model boundary (=%f)",
                 vars->x_int(vars->numPointsInt[kInt]-1,kInt), kInt+1, xMax);
    }
  }

  vars->iii2.resize(vars->maxPointsInt,vars->numInterfaces);
  vars->b_int.resize(vars->maxPointsInt,vars->numInterfaces);
  vars->c_int.resize(vars->maxPointsInt,vars->numInterfaces);
  vars->d_int.resize(vars->maxPointsInt,vars->numInterfaces);

  //
  // END: Read in model
  //********************************************************************************
  //********************************************************************************
  // Read in velocity model
  //
  float* vp_top    = NULL;   // Velocity at top of layer
  float* vp_bottom = NULL;   // Velocity at bottom of layer
  vars->rho1 = new float[vars->numInterfaces];      // Density
  vars->rho2 = new float[vars->numInterfaces];      // Density linear coefficient
  vars->nx_grid = NULL;
  vars->nz_grid = NULL;

  int* index_grid = new int[vars->numInterfaces]; // Flag specifying whether this layer is defined by top/bottom velocity (=1) or a grid (=0)
  vars->vpvs_ratio = new float[vars->numInterfaces]; // Vp/Vs ratio in each layer

  vars->maxPointsXGrid = 2;
  vars->maxPointsZGrid = 2;
  vars->nx_grid = new int[vars->numInterfaces];
  vars->nz_grid = new int[vars->numInterfaces];

  vars->x_grid.resize( vars->maxPointsXGrid, vars->numInterfaces );
  vars->z_grid.resize( vars->maxPointsZGrid, vars->numInterfaces );
  vars->v_grid.resize( 4, vars->maxPointsZGrid, vars->maxPointsXGrid, vars->numInterfaces );

  if( edef->isDebug() ) {
    writer->line("Size: %d %d %d\n", vars->x_grid.size(), vars->z_grid.size(), vars->v_grid.size());
    writer->line("Dimensions: %d %d %d %d\n", 4, vars->maxPointsZGrid, vars->maxPointsXGrid, vars->numInterfaces );
  }

  /*
  bool isVelocity1D = true;
  if( param->exists("velocity") ) {
    param->getString("velocity", &text, 0);
    if( !text.compare("1d") ) {
      isVelocity1D = true;
    }
    else if( !text.compare("2d") ) {
      isVelocity1D = false;
      writer->error("Option not supported yet: 2D velocity field");
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  */

  int numLayers = vars->numInterfaces - 1;
  if( vars->modelType == mod_ray2d::MODEL_1D || vars->modelType == mod_ray2d::MODEL_2D_SIMPLE ) {
    vp_top      = new float[vars->numInterfaces];  // Velocity at top of layer
    vp_bottom   = new float[vars->numInterfaces];  // Velocity at bottom of layer
  }
  /*
  if( vars->modelType == mod_ray2d::MODEL_1D ) {
    for( int ilay = 0; ilay < numLayers; ilay++ ) {
      int numValues = param->getNumValues("model1d", ilay);
      param->getFloatAtLine("model1d", &vp_top[ilay], ilay, 1);
      param->getFloatAtLine("model1d", &vars->vpvs_ratio[ilay], ilay, 2);
      if( shearParam == mod_ray2d::SHEAR_PARAM_VS ) {
        if( vars->vpvs_ratio[ilay] < EPSILON ) vars->vpvs_ratio[ilay] = EPSILON;
        vars->vpvs_ratio[ilay] = vp_top[ilay]/vars->vpvs_ratio[ilay];  // Convert Vs to Vp/Vs
      }
      vars->nx_grid[ilay] = 2;
      vars->nz_grid[ilay] = 2;
      vp_bottom[ilay] = vp_top[ilay];
      vars->rho1[ilay] = 1.0;
      vars->rho2[ilay] = 0.0;
      if( numValues > 3 ) {
        param->getFloatAtLine("model1d", &vars->rho1[ilay], ilay, 3);
        if( numValues > 4 ) {
          param->getFloatAtLine("model1d", &vars->rho2[ilay], ilay, 4);
        }
        if( numValues > 5 ) {
          param->getFloatAtLine("model1d", &vp_bottom[ilay], ilay, 5);
        }
      }
      index_grid[ilay] = 1;  // Flag specifying that this layer is defined by top/bottom velocity, not a grid
    }
  }
  else {
  */
  for( int ilay = 0; ilay < numLayers; ilay++ ) {
    int numValues = param->getNumValues("layer", ilay);
    int layerIndex;
    param->getIntAtLine("layer", &layerIndex, ilay, 0);
    if( layerIndex != ilay+1 ) {
      writer->error("Inconsistent layer numbering for layer #%d (user parameter 'layer'). User specified index: #%d",ilay+1,layerIndex);
    }
    param->getFloatAtLine("layer", &vp_top[ilay], ilay, 1);
    param->getFloatAtLine("layer", &vars->vpvs_ratio[ilay], ilay, 2);
    if( shearParam == mod_ray2d::SHEAR_PARAM_VS ) {
      if( vars->vpvs_ratio[ilay] < EPSILON ) vars->vpvs_ratio[ilay] = EPSILON;
      vars->vpvs_ratio[ilay] = vp_top[ilay]/vars->vpvs_ratio[ilay];  // Convert Vs to Vp/Vs
    }
    vars->nx_grid[ilay] = 2;
    vars->nz_grid[ilay] = 2;
    vp_bottom[ilay] = vp_top[ilay];
    vars->rho1[ilay] = 1.0;
    vars->rho2[ilay] = 0.0;
    if( numValues > 3 ) {
      param->getFloatAtLine("layer", &vars->rho1[ilay], ilay, 3);
      if( numValues > 4 ) {
        param->getFloatAtLine("layer", &vars->rho2[ilay], ilay, 4);
        if( numValues > 5 ) {
          param->getFloatAtLine("layer", &vp_bottom[ilay], ilay, 5);
        }
      }
    }
    index_grid[ilay] = 1;  // Flag specifying that this layer is defined by top/bottom velocity, not a grid
  }

  if( vars->modelType == mod_ray2d::MODEL_1D || vars->modelType == mod_ray2d::MODEL_2D_SIMPLE ) {
    for( int ilay = 0; ilay < numLayers; ilay++ ) {
      vp_top[ilay]    /= 1000.0f; // Convert from [m/s] to [km/s]
      vp_bottom[ilay] /= 1000.0f;
    }
  }

  //  for( int ilay = 0; ilay < numLayers; ilay++ ) {
  //   fprintf(stdout,"VEL: %d  vp: %f %f  vpvs: %f  rho: %f %f\n", ilay, vp_top[ilay], vp_bottom[ilay], vars->vpvs_ratio[ilay], vars->rho1[ilay], vars->rho2[ilay]);
  //  }

  //  float* v_grid = new float[4*maxPointsZGrid*maxPointsXGrid*vars->numInterfaces];
  int flag_smooth = 0;
  int error = 0;

  int f_out = 0;

  for( int kInt = 0; kInt < vars->numInterfaces; kInt++ ) {
    for( int i = 0; i < vars->numPointsInt[kInt]; i++ ) {
      vars->b_int(i,kInt) = -999;
      vars->c_int(i,kInt) = -999;
      vars->d_int(i,kInt) = -999;
    }
  }

  model_( vars->numInterfaces, numBoreholes, vars->numInterfacesAll, vars->numPointsInt, vars->x_int.getPointer(), vars->z_int.getPointer(), vars->iii.getPointer(),
          vars->nx_grid, vars->nz_grid, vars->x_grid.getPointer(), vars->z_grid.getPointer(), vars->v_grid.getPointer(),
          index_grid, vp_top, vp_bottom, vars->vpvs_ratio, f_out, flag_smooth,
          vars->b_int.getPointer(), vars->c_int.getPointer(), vars->d_int.getPointer(),
          vars->left_border, vars->right_border,
          vars->iii2.getPointer(), error,
          vars->maxPointsInt, vars->maxPointsXGrid, vars->maxPointsZGrid);

  if( vp_top != NULL ) {
    delete [] vp_top;
  }
  if( vp_bottom != NULL ) {
    delete [] vp_bottom;
  }

  /*
  for( int kInt = 0; kInt < vars->numInterfaces; kInt++ ) {
    for( int i = 0; i < vars->numPointsInt[kInt]; i++ ) {
      fprintf(stdout,"%2d %2d   %9.2f %9.2f",kInt,i, vars->x_int(i,kInt), vars->z_int(i,kInt));
      fprintf(stdout," %9.2f %9.2f %9.2f", vars->b_int(i,kInt), vars->c_int(i,kInt), vars->d_int(i,kInt));
      fprintf(stdout," %d %d\n", kInt*vars->numInterfacesAll+i, vars->numPointsInt[kInt] );
    }
  }
  */
  
  //--------------------------------------------------
  // Dump model
  //
  if( param->exists("dump_model") ) {
    param->getString("dump_model",&vars->filename_intout);
    if( !csFileUtils::createDoNotOverwrite( vars->filename_intout ) ) {
      writer->error("Error opening output file %s\n", vars->filename_intout.c_str() );
    }
  }

  if( param->exists( "dump_vel" ) ) {
    double velout_dx;
    double velout_dz;
    param->getString( "dump_vel", &text, 0 );
    param->getDouble( "dump_vel", &velout_dx, 1 );
    param->getDouble( "dump_vel", &velout_dz, 2 );
    std::string text2;
    param->getString( "dump_vel", &text2, 3 );
    int velout_flag = mod_ray2d::VELOUT_VP;
    if( !text2.compare("vp") ) {
      velout_flag = mod_ray2d::VELOUT_VP;
    }
    else if( !text2.compare("vs") ) {
      velout_flag = mod_ray2d::VELOUT_VS;
    }
    else if( !text2.compare("vpvs") ) {
      velout_flag = mod_ray2d::VELOUT_VPVS;
    }
    else if( !text2.compare("rho") ) {
      velout_flag = mod_ray2d::VELOUT_RHO;
    }
    else if( !text2.compare("int") ) {
      velout_flag = mod_ray2d::VELOUT_INT;
    }

    csSegyWriter* segyWriter;
    int numTracesBuffer = 20;
    bool rev_byte_order = false;
    bool isSUFormat = false;
    bool autoscale_hdrs = true;
    char* charHdr = new char[cseis_segy::SIZE_CHARHDR+1];
    csSegyHdrMap hdrMap( csSegyHdrMap::SEGY_STANDARD, false );

    try {
      int nz_velout = (int)round( 1000*(zMax-zMin)/(velout_dz) ) + 1;

      segyWriter = new csSegyWriter( text, numTracesBuffer, rev_byte_order, autoscale_hdrs, isSUFormat );
      cseis_geolib::csSegyBinHeader* binHdr = segyWriter->binHdr();
      binHdr->numAuxTraces     = 0;
      binHdr->sampleIntUS     = (unsigned short)(velout_dz*1000);
      binHdr->sampleIntOrigUS = (unsigned short)(velout_dz*1000);
      binHdr->numSamples       = nz_velout;
      binHdr->numSamplesOrig   = nz_velout;
      binHdr->dataSampleFormat = cseis_segy::DATA_FORMAT_IEEE;

      segyWriter->initialize( &hdrMap, charHdr );
      segyWriter->openFile();
      segyWriter->freeCharBinHdr();  // Don't need these headers anymore -> free memory
    }
    catch( csException& e ) {
      segyWriter = NULL;
      writer->error("Error when generating SEGY writer object.\nSystem message: %s", e.getMessage() );
    }
    delete [] charHdr;
    exportVelocityModel( vars, segyWriter, velout_dx/1000, velout_dz/1000, xMin, xMax, zMin, zMax, velout_flag );
    delete segyWriter;
  }

  //********************************************************************************
  //********************************************************************************
  //
  // Read in ray codes
  //
  vars->numCodes = param->getNumLines("ray_code");
  int numCodeStepPlots = 0;
  if( param->exists("plot_code_step") ) numCodeStepPlots = param->getNumValues("plot_code_step");
  if( vars->numCodes == 0 ) {
    writer->line("Error: No ray code specified (user parameter 'ray_code')");
    env->addError();
  }
  vars->numCodeSteps = new int[vars->numCodes];
  vars->codeStepPlot = new int[vars->numCodes];
  vars->maxNumCodeSteps = 0;
  for( int icode = 0; icode < vars->numCodes; icode++ ) {
    vars->numCodeSteps[icode] = param->getNumValues("ray_code",icode) - 1;
    if( vars->numCodeSteps[icode] <= 0 ) {
      writer->line("Error: User parameter 'ray_code' (#%d) specified without value/ray code", icode+1);
      env->addError();
    }
    if( icode < numCodeStepPlots ) {
      param->getInt( "plot_code_step", &vars->codeStepPlot[icode], icode );
    }
    else {
      vars->codeStepPlot[icode] = 0;
    }
    if( vars->numCodeSteps[icode] > vars->maxNumCodeSteps ) {
      vars->maxNumCodeSteps = vars->numCodeSteps[icode];
    }
  }
  vars->rayCode.resize(vars->maxNumCodeSteps+1,vars->numCodes);
  vars->rayHold.resize(vars->maxNumCodeSteps+1,vars->numCodes);

  for( int icode = 0; icode < vars->numCodes; icode++ ) {
    param->getStringAtLine( "ray_code", &text, icode, 0 );
    if( !text.compare("down") ) {
      vars->rayCode(0,icode) = 1;  // Set up/down flag: First ray goes downwards from source
    }
    else if( !text.compare("up") ) {
      vars->rayCode(0,icode) = -1;  // Set up/down flag: First ray goes upwards from source
    }
    else {
      writer->error("Unrecognized option for ray code direction: %s. Specify either 'up' or 'down'", text.c_str());
    }
    if( vars->gatherType == mod_ray2d::SOURCE_GATHER ) {
      for( int istep = 0; istep < vars->numCodeSteps[icode]; istep++ ) {
        param->getIntAtLine( "ray_code", &vars->rayCode(istep+1,icode), icode, istep+1 );
      }
    }
    else { // For receiver gather, reverse order of ray code...
      for( int istep = 0; istep < vars->numCodeSteps[icode]; istep++ ) {
        param->getIntAtLine( "ray_code", &vars->rayCode(istep+1,icode), icode, vars->numCodeSteps[icode]-(istep+1)-1 );
      }
      // ...also, check whether ray code arrives at receiver from above --> set up/down flag 'upwards' = -1
      if( vars->numCodeSteps[icode] > 1 && vars->rayCode(2,icode) < vars->rayCode(1,icode) ) {
        vars->rayCode(0,icode) = -1;
      }
    }
  }


  //********************************************************************************
  // Set up remaining parameters & dummies
  //
  int nreclines = 1;
  
  int retValue = wfront_check_params(
                                     vars->numInterfaces,
                                     vars->x_int,
                                     vars->numPointsInt,
                                     index_grid,
                                     vars->x_grid,
                                     vars->z_grid,
                                     vars->nx_grid,
                                     vars->nz_grid,
                                     vars->sou_layer, vars->sou_interface,
                                     vars->anglein_rad, vars->angleout_rad,
                                     vars->numRays,
                                     nreclines,
                                     &vars->rec_layer,
                                     &vars->rec_interface,
                                     vars->f_wfronts,
                                     stdout,
                                     vars->step_wfronts,
                                     vars->numCodes,
                                     vars->rayCode,
                                     vars->numCodeSteps,
                                     vars->flag_twopoint,
                                     vars->flag_stack );

  delete [] index_grid;
  if( retValue < 0 ) {
    exit(-1);
  }

  error = 0;

  if( edef->isDebug() ) {
    writer->write("\n *** Ray codes BEFORE sorting ***");
    for( int i = 0; i < vars->numCodes; i++ ) {
      writer->write("\nRay code #%d: ", i+1);
      for( int istep = 0; istep <= vars->numCodeSteps[i]; istep++ ) {
        writer->write("%d  ", vars->rayCode(istep,i) );
      }
    }
    writer->line("\n *** Ray codes AFTER sorting ***");
  }

  // Sort ray codes
  sortraycode_( vars->numCodeSteps, vars->numCodes, vars->codeStepPlot, f_out,
                vars->sou_interface, vars->sou_layer, vars->rayCode.getPointer(), vars->rayHold.getPointer(),
                flag_duplicate_raycodes, error, vars->maxNumCodeSteps );

  if( error != 0 ) {
    writer->error("Unknown error occurred in subroutine sortraycode. Error code: %d\n", error);
  }

  char name[11];
  char name2[7];
  char name3[12];
  vars->hdrId_rayTime  = new int[vars->numCodes];
  vars->hdrId_rayAngle = new int[vars->numCodes];
  vars->hdrId_ccpx = new int[vars->numCodes];
  vars->hdrId_ccpz = new int[vars->numCodes];
  for( int i = 0; i < vars->numCodes; i++ ) {
    sprintf(name,"ray_time%02d",i+1);
    name[10] = '\0';
    vars->hdrId_rayTime[i] = hdef->addHeader(TYPE_FLOAT,name);

    sprintf(name3,"ray_angle%02d",i+1);
    name3[11] = '\0';
    vars->hdrId_rayAngle[i] = hdef->addHeader(TYPE_FLOAT,name3);

    sprintf(name2,"ccpx%02d",i+1);
    name[6] = '\0';
    vars->hdrId_ccpx[i] = hdef->addHeader(TYPE_FLOAT,name2);
    sprintf(name2,"ccpz%02d",i+1);
    name[6] = '\0';
    vars->hdrId_ccpz[i] = hdef->addHeader(TYPE_FLOAT,name2);
    if( edef->isDebug() ) {
      writer->write("Ray code #%d: ", i+1);
      for( int istep = 0; istep <= vars->numCodeSteps[i]; istep++ ) {
        writer->write("%d  ", vars->rayCode(istep,i) );
      }
      writer->write("\n");
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_ray2d_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  mod_ray2d::VariableStruct* vars = reinterpret_cast<mod_ray2d::VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;

  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    if( vars->filename_timeout.size() > 0 ) {
      vars->file_timeout = fopen(vars->filename_timeout.c_str(),"w");
      if( vars->file_timeout == NULL ) {
        writer->error("Error opening output file %s\n", vars->filename_timeout.c_str() );
      }
    }
    if( vars->filename_intout.size() > 0 ) {
      FILE* file_intout = fopen(vars->filename_intout.c_str(),"w");
      if( file_intout == NULL ) {
        writer->error("Error opening model output file %s\n", vars->filename_intout.c_str() );
      }
      int npoints_intout = 50;
      float dx = (vars->x_int(vars->numPointsInt[0]-1,0) - vars->x_int(0,0)) / float(npoints_intout-1);

      for( int kInt = 0; kInt < vars->numInterfaces; kInt++ ) {
        int j = 0;
        for( int ip = 0; ip < npoints_intout; ip++ ) {
          float x_output2 = vars->x_int(0,0) + (float)(ip) * dx;
          while( x_output2 > vars->x_int(j+1,kInt) ) {
            j = j + 1;
            fprintf(file_intout,"%f %f\n", vars->x_int(j,kInt),vars->z_int(j,kInt));
          }
          float tmpx = x_output2 - vars->x_int(j,kInt);
          float zout = ((vars->d_int(j,kInt)*tmpx + vars->c_int(j,kInt))*tmpx + vars->b_int(j,kInt))*tmpx + vars->z_int(j,kInt);
          fprintf(file_intout,"%f %f\n", x_output2, zout);
        }
        fprintf(file_intout,"\n");
      }
      fclose(file_intout);
    }
  }

  int numTraces = traceGather->numTraces();

  // Set source & receiver geometry
  // Add one trace in order to coompute zero offset (in case of option NMO_TIME)
  /*  int numTracesInternal = numTraces;
  if( vars->nmoMethod == mod_ray2d::NMO_TIME ) {
    numTracesInternal = numTraces + 1;
    } */
  csMatrixFStyle<float> xz_rec(numTraces,2);
  float sou_x;
  float sou_z;
  if( vars->gatherType == mod_ray2d::RECEIVER_GATHER ) {
    // Input is receiver gather
    // --> Flip receiver/source coordinates to ray-trace the whole gather in one go (ray tracer assumes shot gather geometry)
    if( vars->hdrId_sou_x == -1 ) {
      writer->error("Program bug: Source X coordinate cannot be constant for receiver gather");
    }
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      xz_rec(itrc,0) = (float)traceGather->trace(itrc)->getTraceHeader()->doubleValue(vars->hdrId_sou_x) / 1000.0f;
    }
    if( vars->hdrId_sou_z == -1 ) {
      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        xz_rec(itrc,1) = vars->sou_z;
      }
    }
    else {
      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        xz_rec(itrc,1) = (float)traceGather->trace(itrc)->getTraceHeader()->doubleValue(vars->hdrId_sou_z) / 1000.0f;
      }
    }
    if( vars->hdrId_rec_x != -1 ) {
      sou_x = (float)traceGather->trace(0)->getTraceHeader()->doubleValue(vars->hdrId_rec_x) / 1000.0f;
    }
    else {
      sou_x = vars->rec_x;
    }
    if( vars->hdrId_rec_z != -1 ) {
      sou_z = (float)traceGather->trace(0)->getTraceHeader()->doubleValue(vars->hdrId_rec_z) / 1000.0f;
    }
    else {
      sou_z = vars->rec_z;
    }
  }
  else { // Source gather
    if( vars->hdrId_rec_x == -1 ) {
      writer->error("Program bug: Receiver X coordinate cannot be constant for source gather");
    }
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      xz_rec(itrc,0) = (float)traceGather->trace(itrc)->getTraceHeader()->doubleValue(vars->hdrId_rec_x) / 1000.0f;
    }
    if( vars->hdrId_rec_z == -1 ) {
      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        xz_rec(itrc,1) = vars->rec_z;
      }
    }
    else {
      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        xz_rec(itrc,1) = (float)traceGather->trace(itrc)->getTraceHeader()->doubleValue(vars->hdrId_rec_z) / 1000.0f;
      }
    }
    if( vars->hdrId_sou_x != -1 ) {
      sou_x = (float)traceGather->trace(0)->getTraceHeader()->doubleValue(vars->hdrId_sou_x) / 1000.0f;
    }
    else {
      sou_x = vars->sou_x;
    }
    if( vars->hdrId_sou_z != -1 ) {
      sou_z = (float)traceGather->trace(0)->getTraceHeader()->doubleValue(vars->hdrId_sou_z) / 1000.0f;
    }
    else {
      sou_z = vars->sou_z;
    }
  }

  // Check source position, set sou_dzdx (derivative at source)
  float sou_dzdx;
  int error = 0;
  int f_out = 6; // stdout
  check_sou_loc_( vars->x_int.getPointer(), vars->z_int.getPointer(), vars->d_int.getPointer(),
                  vars->c_int.getPointer(), vars->b_int.getPointer(), vars->numPointsInt,
                  sou_x, sou_z, sou_dzdx,
                  vars->sou_interface, vars->sou_layer,
                  error, f_out, vars->numInterfaces, vars->maxPointsInt );

  if( error != 0 ) {
    writer->error("Unknown error occurred when checking source location (X= %f, Z= %f).", sou_x, sou_z);
  }

  //--------------------------------------------------------------------------------
  // Dummy parameters needed for wfront
  int numBoreholes = 0;
  int maxAngles2Point = 3*numTraces;
  // Set up output
  int step_spread   = 5;
  int step1_wfronts = 2;
  int last_wfronts  = 0;
  if( last_wfronts < step1_wfronts ) {
    last_wfronts = 100000;
  }
  int numReclines = 1;
  int numReceivers = numTraces;
  writer->line("numTraces: %d, nrec: %d %d  %f %f %f\n",
            numTraces, numReceivers, vars->maxNumCodeSteps, vars->dtray, vars->spread_dist, vars->spread_angle_rad);

  if( vars->nmoMethod != mod_ray2d::NMO_NONE ) {
    numReceivers = numTraces + 1;
    csMatrixFStyle<float> xz_rec_new(numReceivers,2);
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      xz_rec_new(itrc,0) = xz_rec(itrc,0);
      xz_rec_new(itrc,1) = xz_rec(itrc,1);
    }
    xz_rec.resize(numReceivers,2);
    for( int itrc = 0; itrc < numReceivers; itrc++ ) {
      xz_rec(itrc,0) = xz_rec_new(itrc,0);
      xz_rec(itrc,1) = xz_rec_new(itrc,1);
    }
  }
  int maxNumReceivers = numReceivers;
  float sou_time  = 0.0;

  int flag_smooth    = 0;
  int flag_surface   = 0;

  csSortManager sortManager( 1, csSortManager::TREE_SORT );
  sortManager.resetValues( numReceivers );
  for( int itrc = 0; itrc < numReceivers; itrc++ ) {
    sortManager.setValue( itrc, 0, xz_rec(itrc,0) );
  }
  sortManager.sort();

  csMatrixFStyle<float> xz_rec_sorted(numReceivers,2);
  for( int itrc = 0; itrc < numReceivers; itrc++ ) {
    int index = sortManager.sortedIndex(itrc);
    xz_rec_sorted(itrc,0) = xz_rec(index,0);
    xz_rec_sorted(itrc,1) = xz_rec(index,1);
  }

  /*
  writer->line("SOU_XZ  %f %f", sou_x, sou_z);
  for( int itrc = 0; itrc < numTraces; itrc++ ) {
    writer->line("REC_XZ  %d %f %f", itrc, xz_rec_sorted(itrc,0), xz_rec_sorted(itrc,1));
  }
  */

  int maxNumArrivals = 20;
  int numTimeCodes = 5 * vars->numCodes + 1;  // "Allocate" enough space for potential additional arrivals (second arrivals, direct arrival etc)
  // !CHANGE! Problem: How to record direct arrival, or multiples? (numCodes is only the number of reflection codes. Direct arrival is in addition to that)
  csMatrixFStyle<float> time;
  csMatrixFStyle<float> amplitude;
  csMatrixFStyle<float> phase;
  csMatrixFStyle<float> angle;
  time.resize( maxNumReceivers, maxNumArrivals, numTimeCodes );
  amplitude.resize( maxNumReceivers, maxNumArrivals, numTimeCodes );
  phase.resize( maxNumReceivers, maxNumArrivals, numTimeCodes );
  angle.resize( maxNumReceivers, maxNumArrivals, numTimeCodes );
  for( int i = 0; i < maxNumReceivers; i++ ) {
    for( int j = 0; j < maxNumArrivals; j++ ) {
      for( int k = 0; k < numTimeCodes; k++ ) {
        time(i,j,k) = 0.0;
        phase(i,j,k) = 0.0;
        angle(i,j,k) = 0.0;
      }
    }
  }

  csMatrixFStyle<float> ccp;
  ccp.resize( 2, maxNumReceivers, numTimeCodes );
  for( int i = 0; i < maxNumReceivers; i++ ) {
    for( int k = 0; k < numTimeCodes; k++ ) {
      ccp(0,i,k) = 0.0;
      ccp(1,i,k) = 0.0;
    }
  }

  //  int numRays = vars->numRays;
  //  fprintf(stdout,"%d %d %d %d %d\n", maxNumReceivers, vars->maxNumCodeSteps, vars->maxPointsXGrid, vars->maxPointsZGrid, maxAngles2Point);
  //  fprintf(stdout,"%d %d %d %d\n", vars->numInterfacesAll, vars->maxPointsInt, maxNumArrivals, numTimeCodes );

  if( edef->isDebug() ) {
    writer->line("RAY2D: Number of input traces: %d, Number of receivers/sources: %d\n", numTraces, numReceivers );
  }

  int componentOut = 3; // 3: Z, 2: X component


  wfsub_( vars->numInterfaces, vars->numPointsInt, vars->x_int.getPointer(), vars->z_int.getPointer(), vars->iii.getPointer(),
          vars->rho1, vars->rho2, vars->vpvs_ratio,
          vars->nx_grid, vars->nz_grid, vars->x_grid.getPointer(), vars->z_grid.getPointer(), vars->v_grid.getPointer(),
          vars->b_int.getPointer(), vars->c_int.getPointer(), vars->d_int.getPointer(),
          vars->left_border, vars->right_border, vars->iii2.getPointer(),
          vars->numRays, vars->sou_layer, vars->sou_interface, sou_x, sou_z, sou_time, vars->sou_dzdx,
          vars->anglein_rad, vars->angleout_rad,
          vars->dtray, vars->spread_dist, vars->spread_angle_rad,
          flag_smooth, vars->flag_compounds, flag_surface, vars->flag_stack, vars->flag_twopoint,
          numReclines, &vars->rec_interface, &vars->rec_layer, &numReceivers, xz_rec_sorted.getPointer(), numBoreholes,
          step_spread, vars->step_wfronts, step1_wfronts, last_wfronts,
          vars->rayCode.getPointer(), vars->numCodeSteps, vars->numCodes, vars->codeStepPlot, vars->rayHold.getPointer(),
          f_out, vars->f_wfronts, vars->f_xrays, vars->f_rayout, vars->f_rayouttmp, vars->f_timeout,
          maxNumReceivers, vars->maxNumCodeSteps,
          vars->maxPointsXGrid, vars->maxPointsZGrid, maxAngles2Point,
          vars->numInterfacesAll, vars->maxPointsInt, maxNumArrivals,
          time.getPointer(), amplitude.getPointer(), phase.getPointer(), angle.getPointer(), componentOut, numTimeCodes, ccp.getPointer() );
  /*
  for( int icode = 0; icode < numTimeCodes; icode++ ) {
    for( int iarr = 0; iarr < maxNumArrivals; iarr++ ) {
      for( int irec = 0; irec < maxNumReceivers; irec++ ) {
        angle(irec,iarr,icode) = fmod( 360 + angle(irec,iarr,icode)*180.0/M_PI - 270, 180 );
      }
    }
  } */
  if( vars->file_timeout != 0 ) {
    for( int icode = 0; icode < numTimeCodes; icode++ ) {
      for( int iarr = 0; iarr < maxNumArrivals; iarr++ ) {
        for( int irec = 0; irec < maxNumReceivers; irec++ ) {
          if( time(irec,iarr,icode) != 0.0 ) {
            fprintf( vars->file_timeout, "%f %f %f %f %f %d\n",
                     xz_rec_sorted(irec,0), time(irec,iarr,icode), amplitude(irec,iarr,icode), phase(irec,iarr,icode), angle(irec,iarr,icode), iarr );
          }
        }
      }
    }
  }

  for( int icode = 0; icode < numTimeCodes; icode++ ) {
    bool found = false;
    for( int irec = 0; irec < maxNumReceivers; irec++ ) {
      if( time(irec,0,icode) != 0.0 ) {
        found = true;
        break;
      }
    }
    if( !found ) { // Skip remaining time codes since they are empty and contain no arrivals
      numTimeCodes = icode;
      break;
    }
  }

  writer->line("Number of computed ray path arrivals: %d   (may include direct arrival path and/or secondary arrivals)", numTimeCodes);


  if( vars->seismogramType != mod_ray2d::SEISMOGRAM_NONE ) {
    float dampingFactor = 0.0;
    float phaseShift_rad = 0.0;
    float sampleInt_s = shdr->sampleInt / 1000.0f;
    float timeSamp1_s = 0.0;
    for( int itrc = 0; itrc < numTraces; itrc++ ) {
      // Make sure to match sorted time array with unsorted input tracegather
      int sortedIndex = sortManager.sortedIndex(itrc);
      float* samples  = traceGather->trace(sortedIndex)->getTraceSamples();

      for( int icode = 0; icode < numTimeCodes; icode++ ) {
        for( int iarr = 0; iarr < maxNumArrivals; iarr++ ) {      
          float eventTime_s    = time(itrc,iarr,icode);
          if( eventTime_s == 0.0 ) continue;
          //      fprintf(stderr,"Event time: %f   %d\n", eventTime_s, iarr);
          float eventAmplitude = 1.0f;
          float eventPhase_rad = 0.0f;
          if( vars->seismogramType == mod_ray2d::SEISMOGRAM_DYNAMIC ) {
            eventPhase_rad = phase(itrc,0,icode);
            eventAmplitude = amplitude(itrc,0,icode);
          }
          bool success = addRickerWavelet( vars->rickerPeakFrequency, dampingFactor, phaseShift_rad, sampleInt_s, timeSamp1_s,
                                           eventTime_s, eventPhase_rad, eventAmplitude, samples, shdr->numSamples );
          if( !success ) {
            // Nothing..
          }
        } // for iarr
      }
    }
  }

  for( int icode = 0; icode < numTimeCodes; icode++ ) {
    int irec = 0;
    while( irec < maxNumReceivers ) {
     //   TEMP   if( time(irec,iarr,icode) != 0.0 ) {
      if( time(irec,0,icode) == 0.0 ) {
        float dx = 0.0;
        float dt = 0.0;
        float da = 0.0;
        float dp = 0.0;
        float dangle = 0.0;
        if( irec > 1 ) {
          dx = xz_rec_sorted(irec-1,0) - xz_rec_sorted(irec-2,0);
          dt = time(irec-1,0,icode) - time(irec-2,0,icode);
          da = amplitude(irec-1,0,icode) - amplitude(irec-2,0,icode);
          dp = phase(irec-1,0,icode) - phase(irec-2,0,icode);
          dangle = angle(irec-1,0,icode) - angle(irec-2,0,icode);
          float ratio = 0.0;
          if( abs(dx) > 0.00001 ) {
            ratio = (xz_rec_sorted(irec,0) - xz_rec_sorted(irec-1,0) ) / dx;
          }
          time(irec,0,icode) = time(irec-1,0,icode) + dt*ratio;
          amplitude(irec,0,icode) = amplitude(irec-1,0,icode) + da*ratio;
          phase(irec,0,icode) = phase(irec-1,0,icode) + dp*ratio;

          angle(irec,0,icode) = angle(irec-1,0,icode) + dangle*ratio;

          float dccpx = ccp(0,irec-1,icode) - ccp(0,irec-2,icode);
          float dccpz = ccp(1,irec-1,icode) - ccp(1,irec-2,icode);
          ccp(0,irec,icode) = ccp(0,irec-1,icode) + dccpx * ratio;
          ccp(1,irec,icode) = ccp(1,irec-1,icode) + dccpz * ratio;
        }
        else {
          // Don't know what to do
        }
      }
      irec += 1;
    } // END while
    irec = numReceivers-1;
    while( irec >= 0 ) {
      if( time(irec,0,icode) == 0.0 ) {
        float dx = 0.0;
        float dt = 0.0;
        float da = 0.0;
        float dp = 0.0;
        float dangle = 0.0;
        if( irec < numReceivers-2 ) {
          dx = xz_rec_sorted(irec+1,0) - xz_rec_sorted(irec+2,0);
          dt = time(irec+1,0,icode) - time(irec+2,0,icode);
          da = amplitude(irec+1,0,icode) - amplitude(irec+2,0,icode);
          dp = phase(irec+1,0,icode) - phase(irec+2,0,icode);
          dangle = angle(irec+1,0,icode) - angle(irec+2,0,icode);
          float ratio = 0.0;
          if( abs(dx) > 0.00001 ) {
            ratio = (xz_rec_sorted(irec,0) - xz_rec_sorted(irec+1,0) ) / dx;
          }
          time(irec,0,icode) = time(irec+1,0,icode) + dt * ratio;
          amplitude(irec,0,icode) = amplitude(irec+1,0,icode) + da * ratio;
          phase(irec,0,icode) = phase(irec+1,0,icode) + dp * ratio;
          angle(irec,0,icode) = angle(irec+1,0,icode) + dangle * ratio;

          float dccpx = ccp(0,irec+1,icode) - ccp(0,irec+2,icode);
          float dccpz = ccp(1,irec+1,icode) - ccp(1,irec+2,icode);
          ccp(0,irec,icode) = ccp(0,irec+1,icode) + dccpx * ratio;
          ccp(1,irec,icode) = ccp(1,irec+1,icode) + dccpz * ratio;
        }
        else {
          // Don't know what to do
        }
      }
      irec -= 1;
    } // END while
  }


  if( vars->nmoMethod != mod_ray2d::NMO_NONE ) {
    csTimeStretch timeStretch( shdr->sampleInt, shdr->numSamples, csTimeStretch::SAMPLE_INTERPOL_SINC );
    timeStretch.setBottomStretch( true );
    int numTimes    = vars->numCodes+1;
    float* buffer   = new float[shdr->numSamples];
    float* timesIn  = new float[numTimes];
    float* timesOut = new float[numTimes];
    timesIn[0]  = 0.0f;
    timesOut[0] = 0.0f;
    for( int itrc = 0; itrc < numReceivers; itrc++ ) {
      // Make sure to match sorted time array with unsorted input tracegather
      int sortedIndex = sortManager.sortedIndex(itrc);
      if( sortedIndex >= numTraces ) continue;
      float* samples = traceGather->trace(sortedIndex)->getTraceSamples();

      memcpy( buffer, samples, shdr->numSamples*sizeof(float) );
      for( int icode = 0; icode < vars->numCodes; icode++ ) {
        //        int index = vars->numCodes-icode-1;
        int index = icode;
        timesIn[index+1]  = time(itrc,0,icode) * 1000.0f;  // Convert from [s] to [ms]
        if( vars->nmoMethod == mod_ray2d::NMO_TIME ) {
          timesOut[index+1] = time(0,0,icode) * 1000.0f;  // Convert from [s] to [ms]
        }
        else {
          timesOut[index+1] = ccp(1,itrc,icode) * 1000.0f;  // Convert from [km] to [m]
        }
      }

      if( edef->isDebug() ) {
        for( int icode = 0; icode <= vars->numCodes; icode++ ) {
          writer->line(" NMO time IN/OUT:  %d   %.2f  %.2f   (%.2f)", icode,  timesIn[icode], timesOut[icode], time(0,0,icode)*1000.0f);
        }
        writer->line("");
      }
      //      float offset =  traceGather->trace(sortedIndex)->getTraceHeader()->floatValue( hdef->headerIndex("offset") );
      //  for( int icode = 0; icode <= vars->numCodes; icode++ ) {
      //    fprintf(stdout,"STRETCH %d %d  %f %f   %f\n", itrc, icode, timesIn[icode], timesOut[icode], offset);
      //  }
      timeStretch.applyTimeInterval( buffer, timesIn, timesOut, numTimes, samples );
    }
    delete [] buffer;
    delete [] timesIn;
    delete [] timesOut;
  }

  for( int itrc = 0; itrc < numReceivers; itrc++ ) {
    // Sort back computed values to original sorting
    int sortedIndex = sortManager.sortedIndex(itrc);
    if( sortedIndex >= numTraces ) continue;
    csTraceHeader* trcHdr = traceGather->trace(sortedIndex)->getTraceHeader();
    for( int icode = 0; icode < vars->numCodes; icode++ ) {
      trcHdr->setFloatValue( vars->hdrId_rayTime[icode], time(itrc,0,icode)*1000.0f );
      trcHdr->setFloatValue( vars->hdrId_rayAngle[icode], angle(itrc,0,icode) );
      trcHdr->setFloatValue( vars->hdrId_ccpx[icode], ccp(0,itrc,icode)*1000.0f );
      trcHdr->setFloatValue( vars->hdrId_ccpz[icode], ccp(1,itrc,icode)*1000.0f );
    }
    if( vars->nmoMethod != mod_ray2d::NMO_NONE && vars->modelType == mod_ray2d::MODEL_1D ) {
      if( vars->nmoMethod == mod_ray2d::NMO_DEPTH ) {
        for( int icode = 0; icode < vars->numCodes; icode++ ) {
          int maxCodeLayer = abs(vars->rayCode(0,icode));
          for( int istep = 1; istep < vars->numCodeSteps[icode]; istep++ ) {
            int rayCode = abs(vars->rayCode(istep,icode));
            if( rayCode > maxCodeLayer ) maxCodeLayer = rayCode;
          }
          if( vars->rayCode(0,icode) < 0 && vars->numCodeSteps[icode] == 1 ) maxCodeLayer += 1; // Direct arrival
          trcHdr->setFloatValue( vars->hdrId_rayTime[icode], vars->z_int(0,maxCodeLayer)*1000 ); // x1000 to convert from [km] to [m]
        }
      }
    }
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_ray2d_( csParamDef* pdef ) {
  pdef->setModule( "RAY2D", "2D isotropic ray tracer" );

  pdef->addParam( "gather", "Type of input trace gather/ensemble", NUM_VALUES_FIXED);
  pdef->addValue( "receiver", VALTYPE_OPTION );
  pdef->addOption( "receiver", "Input ensemble is a receiver gather" );
  pdef->addOption( "source", "Input ensemble is a source gather" );

  pdef->addParam( "model_type", "Type of model", NUM_VALUES_FIXED);
  pdef->addValue( "1d", VALTYPE_OPTION );
  pdef->addOption( "1d", "Specify 1D model", "Each interface is defined by its depth and layer properties (user parameter 'model1d').");
  pdef->addOption( "2d", "Specify 2D model", "Each interface is defined by at least one 2D point (XZ coordinate, user parameter 'int2d'), and local 1D layer properties");
  pdef->addValue( "vs", VALTYPE_OPTION );
  pdef->addOption( "vs", "Specify model S-wave velocity at the top of the layer (will be converted internally into vp/vs ratio that is constant for this layer)");
  pdef->addOption( "vpvs", "Specify model S-wave velocity indirectly as a constant vp/vs ratio");

  pdef->addParam( "rec_x", "Receiver X coordinate", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "Value or trace header name for receiver X coordinate [m]" );

  pdef->addParam( "rec_z", "Receiver Z coordinate", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "Value or trace header name for receiver Z coordinate [m]" );

  pdef->addParam( "sou_x", "Source X coordinate", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "Value or trace header name for source X coordinate [m]" );

  pdef->addParam( "sou_z", "Source Z coordinate", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "Value or trace header name for source Z coordinate [m]" );

  pdef->addParam( "sou_loc", "Source location", NUM_VALUES_VARIABLE);
  pdef->addValue( "1", VALTYPE_HEADER_NUMBER, "Layer in which source is located" );
  pdef->addValue( "0", VALTYPE_HEADER_NUMBER, "Interface on which source is located (if any)" );

  pdef->addParam( "rec_loc", "Receiver location", NUM_VALUES_VARIABLE);
  pdef->addValue( "2", VALTYPE_HEADER_NUMBER, "Layer in which receiver is located" );
  pdef->addValue( "2", VALTYPE_HEADER_NUMBER, "Interface on which receiver is located (if any)" );

  pdef->addParam( "apply_nmo", "Apply ray-traced normal moveout correction?", NUM_VALUES_FIXED);
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not apply NMO" );
  pdef->addOption( "time", "Apply NMO, output zero-offset time [ms]" );
  pdef->addOption( "depth", "Apply NMO, output depth [m]" );

  //  pdef->addParam( "apply_binning", "Apply ray-traced binning correction?", NUM_VALUES_FIXED);
  //  pdef->addValue( "yes", VALTYPE_OPTION );
  //  pdef->addOption( "yes", "Apply binning" );
  //  pdef->addOption( "no", "Do not apply binning" );

  pdef->addParam( "seismogram", "Generate seismogram by adding Ricker wavelet to input traces at ray-traced travel times", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_NUMBER, "Peak frequency [Hz] of Ricker wavelet" );
  pdef->addValue( "dynamic", VALTYPE_OPTION );
  pdef->addOption( "dynamic", "Generate seismogram including amplitude/phase correction" );
  pdef->addOption( "kinematic", "Generate seismogram with constant wavelet" );

  pdef->addParam( "ray_code", "Ray code definition", NUM_VALUES_VARIABLE, "Each ray code consists of a list of layer numbers through which the ray travels");
  pdef->addValue( "down", VALTYPE_OPTION );
  pdef->addOption( "down", "Ray travels downwards from source" );
  pdef->addOption( "up", "Ray travels upwards from source" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of layer numbers defining one ray path from source to receiver, for example '1 2 2 1' for a ray traveling from the source through layer 1, then through layer 2, then reflecting back through layer 2, then through layer 1 to the receiver", "Specify negative layer numbers for S-waves" );

  pdef->addParam( "plot_code_step", "If applicable: First ray code step to plot (wave fronts, rays...)", NUM_VALUES_VARIABLE, "Used for wave front plots and ray plots, if specified" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of values: For each used-defined ray code ('ray_code'), specify the first code step which shall be plotted, starting with 1 for the first code step" );

  pdef->addParam( "max_offset", "Define maximum source-receiver offset for 1D model", NUM_VALUES_FIXED );
  pdef->addValue( "10000", VALTYPE_NUMBER, "Maximum source-receiver offset [m] in case of 1D model" );

  pdef->addParam( "int1d", "Define 1D model interface", NUM_VALUES_VARIABLE, "Repeat this parameter for each model interface.");
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of interface. First interface = 1" );
  pdef->addValue( "", VALTYPE_NUMBER, "Depth of interface [m]" );

  pdef->addParam( "int2d", "Define one point of one 2D model interface", NUM_VALUES_VARIABLE );
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of interface. First interface = 1" );
  pdef->addValue( "", VALTYPE_NUMBER, "X coordinate [m]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Z coordinate [m]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Optional: Point flag. 0 for continuous interface point, -1 for discontinuous",
                  "Corner points are automatically set to -1" );

  pdef->addParam( "layer", "Define local 1D layer properties", NUM_VALUES_VARIABLE,
                  "Repeat this parameter for each model layer. Layer numbers refer to top interface. The model bottom interface does not require layer properties");
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of layer. First layer = 1" );
  pdef->addValue( "", VALTYPE_NUMBER, "P-wave velocity at top of layer [m/s]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Vp/Vs velocity ratio, or S-wave velocity at top of layer [m/s], see user parameter 'model_type'" );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Density (=rho1)" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Density coefficient (=rho2)", "Density is computed as  rho = rho1 + vp*rho2" );
  pdef->addValue( "", VALTYPE_NUMBER, "Optional: P-wave velocity at bottom of layer [m/s]",
                  "If specified, P-wave velocity function is a vertical gradient. If not specified, layer P-wave velocity is set constant" );

  pdef->addParam( "dt_ray", "Ray propagation time interval", NUM_VALUES_FIXED );
  pdef->addValue( "40", VALTYPE_NUMBER, "Time interval [ms]" );

  pdef->addParam( "spread_dist", "Maximum spreading distance between rays before new ray is interpolated", NUM_VALUES_FIXED );
  pdef->addValue( "60", VALTYPE_NUMBER, "Maximum spreading distance [m]" );

  pdef->addParam( "spread_angle", "Maximum spreading angle between rays before new ray is interpolated", NUM_VALUES_FIXED );
  pdef->addValue( "2.3", VALTYPE_NUMBER, "Maximum spreading angle [deg]" );

  pdef->addParam( "nrays", "Number of rays starting from the source", NUM_VALUES_FIXED );
  pdef->addValue( "400", VALTYPE_NUMBER );

  pdef->addParam( "beam", "Definition of beam of rays departing from source location", NUM_VALUES_FIXED, "Starting rays are equally spaced inside the beam" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Start angle [deg]" );
  pdef->addValue( "180.0", VALTYPE_NUMBER, "End angle [deg]" );

  //  pdef->addParam( "layer2d", "Define local 2D layer properties (for 2D model)", NUM_VALUES_VARIABLE );
  //  grid based definition of velocities etc

  //  pdef->addParam( "rho2d", "Define layer density", NUM_VALUES_VARIABLE );
  //  pdef->addValue( "1", VALTYPE_NUMBER, "Number of layer. First layer = 1" );

  //  pdef->addOption( "2d", "Use gridded velocity model with arbitrary knee points",
  //                 "Each layer velocity is defined by a 2D grid (user parameter 'vel_grid')");

  //  pdef->addParam( "filename_int2d", "Read in 2D model from ASCII file", NUM_VALUES_FIXED,
  //                "Input format is   X Z FLAG.  Use empty lines to separate interfaces" );
  // pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "dump_model", "Dump model into ASCII file", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output ASCII file" );

  pdef->addParam( "flags", "Various flags", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Two-point ray-tracing flag. 1: Only compute rays for two-point" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Stack flag. 1: Simulate zero-offset geometry" );
  pdef->addValue( "0", VALTYPE_NUMBER, "'Compound' ray flag. 1: Do not trace 'compound' rays = turning/diving waves" );

  pdef->addParam( "dump_rays", "Dump rays", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Fortran output file number for rays (temporary)", "Output file will be named 'fort.<num>'" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Fortran output file number for rays (final)", "Output file will be named 'fort.<num>'" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Ray step (only output evey n'th ray)" );

  pdef->addParam( "dump_wfronts", "Dump wave fronts", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Fortran output file number for wave fronts", "Output file will be named 'fort.<num>'" );
  pdef->addValue( "4", VALTYPE_NUMBER, "Wave front 'step' (only output every n'th wave front)" );

  pdef->addParam( "dump_vel", "Dump veocity model", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Output cell interval in horizontal direction [m]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Output cell interval in vertical direction [m]" );
  pdef->addValue( "vp", VALTYPE_OPTION );
  pdef->addOption( "vp", "Output Vp" );
  pdef->addOption( "vs", "Output Vs" );
  pdef->addOption( "int", "Output interface number" );
  //  pdef->addOption( "vpvs", "Output Vp/Vs ratio" );
  //  pdef->addOption( "rho", "Output density" );

  pdef->addParam( "dump_times", "Dump travel times", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Filename of output ASCII file" );

  pdef->addParam( "duplicate_raycodes", "Allow duplicate ray codes?", NUM_VALUES_FIXED, "If 'yes', each user-defined ray code will be ray-traced. Otherwise, ray codes contained in other (longer) ray codes will not be ray-traced again");
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not allow duplicate ray codes" );
  pdef->addOption( "yes", "Allow duplicate ray codes", "This is the preferred mode when applying ray-traced NMO, since only one arrival per user-defined ray code is used for NMO" );


  pdef->addDoc("This module computes travel times using a 2D isotropic, dynamic ray tracer, and applies");
  pdef->addDoc("travel time corrections to the data, equivalent to a horizon based normal moveout correction.");
  pdef->addDoc("The ray tracer can be supplied with an almost arbitrarily complex 2D model. It is based on Psencik's SEIS88, with 2d wavefront construction on top of it. The module is rather difficult to set up since the code was barely made to work from old Fortran 77 code. The following example should help:");

  pdef->addDoc("<br><br><kbd>");
  pdef->addDoc("$INPUT_CREATE<br>");
  pdef->addDoc(" ntraces      120<br>");
  pdef->addDoc(" length       5000<br>");
  pdef->addDoc(" sample_int   2<br>");
  pdef->addDoc(" value        0.0<br>");
  pdef->addDoc("<br>");
  pdef->addDoc("$HDR_MATH<br>");
  pdef->addDoc(" new rec_x<br>");
  pdef->addDoc(" new sou_x<br>");
  pdef->addDoc(" equation rec_x \"trcno*50.0\"<br>");
  pdef->addDoc(" equation sou_x \"0.0\"<br>");
  pdef->addDoc("<br>");
  pdef->addDoc("$RAY2D<br>");
  pdef->addDoc(" seismogram   20 dynamic<br>");
  pdef->addDoc(" gather       source<br>");
  pdef->addDoc(" model_type   1d  vs<br>");
  pdef->addDoc(" max_offset   8000<br>");
  pdef->addDoc(" dt_ray       40<br>");
  pdef->addDoc(" spread_dist  20<br>");
  pdef->addDoc(" spread_angle 2.3<br>");
  pdef->addDoc(" nrays        400<br>");
  pdef->addDoc(" rec_x    rec_x<br>");
  pdef->addDoc(" rec_z    1620<br>");
  pdef->addDoc(" rec_loc  2 2<br>");
  pdef->addDoc(" sou_x    sou_x<br>");
  pdef->addDoc(" sou_z    10<br>");
  pdef->addDoc(" sou_loc  1<br>");
  pdef->addDoc(" int1d 1 0<br>");
  pdef->addDoc(" int1d 2 1620<br>");
  pdef->addDoc(" int1d 3 1780<br>");
  pdef->addDoc(" int1d 4 2100<br>");
  pdef->addDoc(" int1d 5 3300<br>");
  pdef->addDoc(" layer 1 1475 1<br>");
  pdef->addDoc(" layer 2 1500 200<br>");
  pdef->addDoc(" layer 3 1670 280<br>");
  pdef->addDoc(" layer 4 2040 600 <br>");
  pdef->addDoc("layer 5 3500 1400<br>");
  pdef->addDoc(" # P-wave ray paths (P-wave: Positive layer numbers, S-wave: Negative layer numbers)<br>");
  pdef->addDoc(" ray_code 1 2 2<br>");
  pdef->addDoc(" ray_code 1 2 3 3 2<br>");
  pdef->addDoc(" <br>");
  pdef->addDoc("$OUTPUT<br>");
  pdef->addDoc("  filename test_ray2d.cseis<br></kbd>");
}

//--------------------------------------------------------------------------------
//
void compute_z( mod_ray2d::VariableStruct* vars, float xvalue, int intIndex, float* zout, int* pointIndex_out ) {
  int pointIndex = 0;
  do {
    if( xvalue < vars->x_int(pointIndex+1,intIndex) ) break;
    pointIndex += 1;
  } while( pointIndex < vars->numPointsInt[intIndex]-1 );
  pointIndex = std::min( pointIndex, vars->numPointsInt[intIndex]-2 );
  *pointIndex_out = pointIndex;
  float dx = xvalue - vars->x_int(pointIndex,intIndex);
  *zout = ((vars->d_int(pointIndex,intIndex)*dx + vars->c_int(pointIndex,intIndex))*dx + vars->b_int(pointIndex,intIndex))*dx + vars->z_int(pointIndex,intIndex);
}
//--------------------------------------------------------------------------------
//
void exportVelocityModel( mod_ray2d::VariableStruct* vars, cseis_geolib::csSegyWriter* segyWriter, double dx, double dz, float xMin, float xMax, float zMin, float zMax, int velout_flag ) {
  int nx_velout = (int)round( (xMax-xMin)/dx ) + 1;
  int nz_velout = (int)round( (zMax-zMin)/dz ) + 1;

  if( nx_velout < 2 ) {
    nx_velout = 2;
    dx = xMax-xMin;
  }
  if( nz_velout < 2 ) {
    nz_velout = 2;
    dz = zMax-zMin;
  }

  //  fprintf(stderr,"min/max:   X: %f %f   Z: %f %f\n", xMin, xMax, zMin, zMax);
  //  fprintf(stderr,"nx/nz   %d %d\n", nx_velout, nz_velout);

  float* samples = new float[nz_velout];
  float* interface_depths = new float[vars->numInterfaces];

  int veltype = 1;
  if( velout_flag == mod_ray2d::VELOUT_VS ) {
    veltype = -1;
  }
  //  float ptos = 2.0;

  cseis_geolib::csSegyTraceHeader* segyTrcHdr = segyWriter->getTraceHeader(); 
  int hdrId_chan = segyTrcHdr->headerIndex( cseis_geolib::HDR_CHAN.name );
  int hdrId_recx = segyTrcHdr->headerIndex( cseis_geolib::HDR_REC_X.name );
  int hdrId_soux = segyTrcHdr->headerIndex( cseis_geolib::HDR_SOU_X.name );

  for( int ix = 0; ix < nx_velout; ix++ ) {
    float x_output = (float)( xMin + (double)ix * dx );

    segyTrcHdr->setIntValue( hdrId_chan, ix+1 );
    segyTrcHdr->setDoubleValue( hdrId_soux, x_output*1000 ); // x1000 to convert from [km] to [m]
    segyTrcHdr->setDoubleValue( hdrId_recx, x_output*1000 ); // x1000 to convert from [km] to [m]

    for( int intIndex = 0; intIndex < vars->numInterfaces; intIndex++ ) {
      int pointIndex;
      compute_z( vars, x_output, intIndex, &interface_depths[intIndex], &pointIndex );
      //      fprintf(stdout,"DEPTH %d %.4f %.4f\n", intIndex, x_output, interface_depths[intIndex]);
    }

    for( int iz = 0; iz < nz_velout; iz++ ) {
      float z_output = (float)( zMin + (double)iz * dz );
      
      int intIndex = 0;
      do {
        if( interface_depths[intIndex] > z_output ) {
          if( intIndex > 0 ) intIndex -= 1;
          else break;
        }
        else if( interface_depths[intIndex+1] < z_output ) {
          if( intIndex < vars->numInterfaces-2 ) intIndex += 1;
          else break;
        }
        else {
          break;
        }
      } while( intIndex < vars->numInterfaces-1 );

      int actlayer = intIndex+1;
      float velocity;
      if( velout_flag == mod_ray2d::VELOUT_INT ) {
        velocity = intIndex;
      }
      else {
        velonly_( x_output, z_output, vars->nx_grid[intIndex], vars->nz_grid[intIndex],
                  &vars->x_grid(0,intIndex), &vars->z_grid(0,intIndex), &vars->v_grid( 0, 0, 0, intIndex ),
                  veltype, vars->vpvs_ratio[intIndex], actlayer, velocity,
                  vars->maxPointsXGrid, vars->maxPointsZGrid );
        //        if( velout_flag == mod_ray2d::VELOUT_VP ) {
        // }
        // else if( velout_flag == mod_ray2d::VELOUT_VS ) {
        //  }
      }

      samples[iz] = velocity;
      //      fprintf(stdout,"%7.4f %7.4f  %f   %d\n", x_output, z_output, velocity, intIndex );
    } // END for iz
    try {
      segyWriter->writeNextTrace( (byte_t*)samples, nz_velout );
    }
    catch( csException& exc ) {
      fprintf(stderr, "System error message: %s" ,exc.getMessage() );
    }
  } // END for ix

  delete [] samples;
  delete [] interface_depths;
}

//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_ray2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_ray2d::VariableStruct* vars = reinterpret_cast<mod_ray2d::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_ray2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_ray2d::VariableStruct* vars = reinterpret_cast<mod_ray2d::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->file_timeout != 0 ) {
    fclose( vars->file_timeout );
    vars->file_timeout = NULL;
  }
  if( vars->hdrId_rayTime != NULL ) {
    delete [] vars->hdrId_rayTime;
    vars->hdrId_rayTime = NULL;
  }
  if( vars->hdrId_rayAngle != NULL ) {
    delete [] vars->hdrId_rayAngle;
    vars->hdrId_rayAngle = NULL;
  }
  if( vars->hdrId_ccpx != NULL ) {
    delete [] vars->hdrId_ccpx;
    vars->hdrId_ccpx = NULL;
  }
  if( vars->hdrId_ccpz != NULL ) {
    delete [] vars->hdrId_ccpz;
    vars->hdrId_ccpz = NULL;
  }
  if( vars->numPointsInt != NULL ) {
    delete [] vars->numPointsInt;
    vars->numPointsInt = NULL;
  }
  if( vars->rho1 != NULL ) {
    delete [] vars->rho1;
    vars->rho1 = NULL;
  }
  if( vars->rho2 != NULL ) {
    delete [] vars->rho2;
    vars->rho2 = NULL;
  }
  if( vars->nx_grid != NULL ) {
    delete [] vars->nx_grid;
    vars->nx_grid = NULL;
  }
  if( vars->nz_grid != NULL ) {
    delete [] vars->nz_grid;
    vars->nz_grid = NULL;
  }
  if( vars->codeStepPlot != NULL ) {
    delete [] vars->codeStepPlot;
    vars->codeStepPlot = NULL;
  }
  if( vars->numCodeSteps != NULL ) {
    delete [] vars->numCodeSteps;
    vars->numCodeSteps = NULL;
  }
  if( vars->vpvs_ratio != NULL ) {
    delete [] vars->vpvs_ratio;
    vars->vpvs_ratio = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_ray2d_( csParamDef* pdef ) {
  params_mod_ray2d_( pdef );
}
extern "C" void _init_mod_ray2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_ray2d_( param, env, writer );
}
extern "C" bool _start_exec_mod_ray2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_ray2d_( env, writer );
}
extern "C" void _exec_mod_ray2d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_ray2d_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_ray2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_ray2d_( env, writer );
}

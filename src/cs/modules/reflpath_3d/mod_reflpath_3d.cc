/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csGeolibUtils.h"
#include "csPlane3D.h"
#include "csPoint3D.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: REFLPATH_3D
 *
 * @author Bjorn Olofsson
 * @date   2010
 */
namespace mod_reflpath_3d {
  struct VariableStruct {
    int hdrId_p1;
    int hdrId_p2;
    int hdrId_p3;
    int hdrId_rayCMP;
    int hdrId_rayDist;

    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_sou_z;

    int hdrId_rec_x;
    int hdrId_rec_y;
    int hdrId_rec_z;
    
    cseis_geolib::csPlane3D* planeRefl;

    std::string filename_plane;
    FILE* fout_plane;
    double out_xmin;
    double out_xmax;
    double out_ymin;
    double out_ymax;
  };
}
using mod_reflpath_3d::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_reflpath_3d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  //  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->hdrId_p1 = -1;
  vars->hdrId_p2 = -1;
  vars->hdrId_p3 = -1;
  vars->hdrId_rayCMP  = -1;
  vars->hdrId_rayDist = -1;
  vars->planeRefl = new csPlane3D();

  vars->fout_plane = NULL;
  vars->filename_plane = "";
  vars->out_xmin = 0;
  vars->out_xmax = 0;
  vars->out_ymin = 0;
  vars->out_ymax = 0;

  std::string name1("plane_p1");
  std::string name2("plane_p2");
  std::string name3("plane_p3");
  std::string name;

  if( param->exists("hdr_in") ) {
    param->getString("hdr_in", &name1, 0);
    param->getString("hdr_in", &name2, 1);
    param->getString("hdr_in", &name3, 2);
  }
  vars->hdrId_p1 = hdef->headerIndex(name1);
  vars->hdrId_p2 = hdef->headerIndex(name2);
  vars->hdrId_p3 = hdef->headerIndex(name3);

  std::string nameRayCMP("ray_cmp");
  std::string nameRayDist("ray_dist");
  if( param->exists("hdr_out") ) {
    param->getString("hdr_out", &nameRayCMP, 0);
    param->getString("hdr_out", &nameRayDist, 1);
  }
  if( !hdef->headerExists(nameRayCMP) ) {
    hdef->addHeader( TYPE_VECTOR, nameRayCMP, "Ray CMP point XYZ coordinates" );
    /*    int pos = nameRayCMP.size() + 1;
    name = nameRayCMP;
    name += ".x";
    hdef->addHeader( TYPE_DOUBLE, name, "Ray CMP point X coordinate" );
    name.replace( pos, 1, "y" );
    hdef->addHeader( TYPE_DOUBLE, name, "Ray CMP point Y coordinate" );
    name.replace( pos, 1, "z" );
    hdef->addHeader( TYPE_DOUBLE, name, "Ray CMP point Z coordinate" );
    */
  }
  if( !hdef->headerExists(nameRayDist) ) {
    hdef->addHeader( TYPE_DOUBLE, nameRayDist );
  }
  //  else if( !hdef->headerType(nameOut) != cseis_geolib::TYPE_VECTOR_DOUBLE ) {
  //    writer->error("Output trace header '%s' must be of vector_double type, i.e. must be three headers name.x name.y and name.z", nameOut.c_str());
  //  }
  vars->hdrId_rayCMP = hdef->headerIndex(nameRayCMP);
  vars->hdrId_rayDist = hdef->headerIndex(nameRayDist);

  if( param->exists("hdr_sou_xyz") ) {
    param->getString("hdr_sou_xyz", &name, 0);
    vars->hdrId_sou_x = hdef->headerIndex(name);
    param->getString("hdr_sou_xyz", &name, 1);
    vars->hdrId_sou_y = hdef->headerIndex(name);
    param->getString("hdr_sou_xyz", &name, 2);
    vars->hdrId_sou_z = hdef->headerIndex(name);
  }
  else {
    vars->hdrId_sou_x = hdef->headerIndex(cseis_geolib::HDR_SOU_X.name);
    vars->hdrId_sou_y = hdef->headerIndex(cseis_geolib::HDR_SOU_Y.name);
    vars->hdrId_sou_z = hdef->headerIndex(cseis_geolib::HDR_SOU_Z.name);
  }
  if( param->exists("hdr_rec_xyz") ) {
    param->getString("hdr_rec_xyz", &name, 0);
    vars->hdrId_rec_x = hdef->headerIndex(name);
    param->getString("hdr_rec_xyz", &name, 1);
    vars->hdrId_rec_y = hdef->headerIndex(name);
    param->getString("hdr_rec_xyz", &name, 2);
    vars->hdrId_rec_z = hdef->headerIndex(name);
  }
  else {
    vars->hdrId_rec_x = hdef->headerIndex(cseis_geolib::HDR_REC_X.name);
    vars->hdrId_rec_y = hdef->headerIndex(cseis_geolib::HDR_REC_Y.name);
    vars->hdrId_rec_z = hdef->headerIndex(cseis_geolib::HDR_REC_Z.name);
  }

  if( param->exists("output_plane") ) {
    std::string text;
    param->getString("output_plane", &vars->filename_plane, 0);
    param->getDouble("output_plane", &vars->out_xmin, 1);
    param->getDouble("output_plane", &vars->out_xmax, 2);
    param->getDouble("output_plane", &vars->out_ymin, 3);
    param->getDouble("output_plane", &vars->out_ymax, 4);
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_reflpath_3d_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;

  if( vars->filename_plane.size() > 0 ) {
    vars->fout_plane = fopen(vars->filename_plane.c_str(),"w");
    if( vars->fout_plane == NULL ) writer->error("Cannot open ASCII output file %s", vars->filename_plane.c_str() );
  }

  int numTraces = traceGather->numTraces();
  if( numTraces == 0 ) return;

  FILE* debugStream = edef->isDebug() ? writer->getFile() : NULL;

  csTraceHeader* trcHdr = traceGather->trace(0)->getTraceHeader();
  csPoint3D p1 = trcHdr->vectorValue( vars->hdrId_p1 );
  csPoint3D p2 = trcHdr->vectorValue( vars->hdrId_p2 );
  csPoint3D p3 = trcHdr->vectorValue( vars->hdrId_p3 );
  vars->planeRefl->setFromPoints( p1, p2, p3 );
  writer->line("Reflector plane dip/azim angles:  %.3f  %.3f", vars->planeRefl->dipAngle(), vars->planeRefl->azimuthAngle() );
  writer->line("Reflector plane normal:           %.5f  %.5f  %.5f", vars->planeRefl->normalVector().x, vars->planeRefl->normalVector().y, vars->planeRefl->normalVector().z );

  if( vars->fout_plane != NULL ) {
    int num = 20;
    double dx = (vars->out_xmax - vars->out_xmin) / (num-1);
    double dy = (vars->out_ymax - vars->out_ymin) / (num-1);
    for( int ix = 0; ix < num; ix++ ) {
      double x = ix*dx + vars->out_xmin;
      for( int iy = 0; iy < num; iy++ ) {
        double y = iy*dy + vars->out_ymin;
        fprintf(vars->fout_plane,"%.2f %.2f %.4f  %.5f %.5f  x/y/z[m] dip/azim[deg]\n", x, y, vars->planeRefl->computeZ( x, y ), vars->planeRefl->dipAngle(), vars->planeRefl->azimuthAngle() );
      }
    }
  }

  for( int itrc = 0; itrc < numTraces; itrc++ ) {
    trcHdr = traceGather->trace(itrc)->getTraceHeader();
    csPoint3D pSrc( trcHdr->doubleValue(vars->hdrId_sou_x), trcHdr->doubleValue(vars->hdrId_sou_y), trcHdr->doubleValue(vars->hdrId_sou_z) );
    csPoint3D pRcv( trcHdr->doubleValue(vars->hdrId_rec_x), trcHdr->doubleValue(vars->hdrId_rec_y), trcHdr->doubleValue(vars->hdrId_rec_z) );

    cseis_geolib::rayPathParam param = vars->planeRefl->computeReflectionPoint( pSrc, pRcv, debugStream );
    trcHdr->setVectorValue( vars->hdrId_rayCMP, param.cmp );
    trcHdr->setDoubleValue( vars->hdrId_rayDist, param.distanceSrcCMP + param.distanceRcvCMP );
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_reflpath_3d_( csParamDef* pdef ) {
  pdef->setModule( "REFLPATH_3D", "Compute 3D ray path assuming straight rays for a reflection from a 3D plane reflector, assuming constant velocity etc." );

  pdef->addParam( "hdr_in", "Input trace headers describing reflector plane: Three (X,Y,Z) vector points defining a 3D plane", NUM_VALUES_FIXED, "Run HDR_MATH_ENS, method 'fit_plane', to create these trace headers" );
  pdef->addValue( "plane_p1", VALTYPE_STRING, "" );
  pdef->addValue( "plane_p2", VALTYPE_STRING, "" );
  pdef->addValue( "plane_p3", VALTYPE_STRING, "" );

  pdef->addParam( "hdr_out", "Output trace header to store sub-surface reflection point", NUM_VALUES_FIXED );
  pdef->addValue( "ray_cmp", VALTYPE_STRING, "3D vector point / CMP position (ray_cmp.x, ray_cmp.y, ray_cmp.z)" );
  pdef->addValue( "ray_dist", VALTYPE_STRING, "Ray path distance from source via CMP to receiver" );

  pdef->addParam( "hdr_sou_xyz", "Input trace headers containing source XYZ coordinates", NUM_VALUES_FIXED );
  pdef->addValue( "sou_x", VALTYPE_STRING );
  pdef->addValue( "sou_y", VALTYPE_STRING );
  pdef->addValue( "sou_z", VALTYPE_STRING );

  pdef->addParam( "hdr_rec_xyz", "Input trace headers containing receiver XYZ coordinates", NUM_VALUES_FIXED );
  pdef->addValue( "rec_x", VALTYPE_STRING );
  pdef->addValue( "rec_y", VALTYPE_STRING );
  pdef->addValue( "rec_z", VALTYPE_STRING );

  pdef->addParam( "output_plane", "Specify to write 3D reflector plane to ASCII file", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Name of ASCII output file" );
  pdef->addValue( "", VALTYPE_NUMBER, "Minimum X coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum X coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Minimum Y coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum Y coordinate" );

  // pdef->addParam( "option", "...", NUM_VALUES_FIXED );
  // pdef->addValue( "no", VALTYPE_OPTION );
  // pdef->addOption( "yes", "..." );
  // pdef->addOption( "no", "..." );
}

//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_reflpath_3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_reflpath_3d::VariableStruct* vars = reinterpret_cast<mod_reflpath_3d::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_reflpath_3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_reflpath_3d::VariableStruct* vars = reinterpret_cast<mod_reflpath_3d::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->planeRefl != NULL ) {
    delete vars->planeRefl;
    vars->planeRefl = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_reflpath_3d_( csParamDef* pdef ) {
  params_mod_reflpath_3d_( pdef );
}
extern "C" void _init_mod_reflpath_3d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_reflpath_3d_( param, env, writer );
}
extern "C" bool _start_exec_mod_reflpath_3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_reflpath_3d_( env, writer );
}
extern "C" void _exec_mod_reflpath_3d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_reflpath_3d_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_reflpath_3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_reflpath_3d_( env, writer );
}

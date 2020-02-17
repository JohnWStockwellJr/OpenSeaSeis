/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_methods.h"
#include <cmath>
#include "csPoint3D.h"
#include "csProfile.h"
#include "csTableAll.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: MIRROR
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_mirror {
  struct VariableStruct {
    csTableAll* table;
    int hdrId_bin_x;
    int hdrId_bin_y;
    int hdrId_bin_z;
    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_sou_z;
    int hdrId_rec_x;
    int hdrId_rec_y;
    int hdrId_rec_z;
    double constantDepth;
    double resolution;   // Distance between knee points in [m]
    double threshold;

    int method;
  };
  static int const METHOD_OBS_MIRROR = 1;
}
using namespace mod_mirror;

double tryLamda( double lamda, csPoint3D const& p_src, csProfile const* profile, bool doPrint = false );


//*************************************************************************************************
// Init phase
//
//
//*************************************************************************************************
void init_mod_mirror_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->table        = NULL;
  vars->hdrId_bin_x  = -1;
  vars->hdrId_bin_y  = -1;
  vars->hdrId_bin_z  = -1;
  vars->constantDepth = 0;
  vars->resolution   = 500;
  vars->threshold    = -1;

  vars->hdrId_sou_x    = -1;
  vars->hdrId_sou_y    = -1;
  vars->hdrId_sou_z    = -1;
  vars->hdrId_rec_x    = -1;
  vars->hdrId_rec_y    = -1;
  vars->hdrId_rec_z    = -1;

//---------------------------------------------
// Retrieve table
//
  std::string text;
  
  if( param->exists("table") ) {
    param->getString("table", &text );
    try {
      vars->table = new csTableAll( csTableAll::TABLE_TYPE_UNIQUE_KEYS );
      vars->table->initialize( text );
      vars->table->readTableContents( true );
    }
    catch( csException& exc ) {
      writer->error("Error when initializing input table '%s':\n %s\n", text.c_str(), exc.getMessage() );
    }
    if( vars->table->numKeys() != 2 ) {
      writer->error("Number of table keys in input table = %d. Required number of keys = 2. Specify table key by placing the character '%c' in front of the key name. Example: %cx %cy depth  (x and y are table keys)", csTableAll::KEY_CHAR, csTableAll::KEY_CHAR, csTableAll::KEY_CHAR );
    }
    if( edef->isDebug() ) {
      vars->table->dump();
    }
    if( param->exists("resolution") ) {
      param->getDouble("resolution", &vars->resolution );
    }
    if( param->exists("target_diff") ) {
      param->getDouble("target_diff", &vars->threshold );
    }
  }
  else {
    param->getDouble("const_depth",&vars->constantDepth);
  }

  //-------------------------------------------------------------
  //
  /*
  if( param->exists("method") ) {
    param->getString("method",&text);
    if( !text.compare("normal") ) {
      vars->method = METHOD_NORMAL;
    }
    else if( !text.compare("obc_asym") ) {
      vars->method = METHOD_OBC_ASYM;
      param->getFloat( "depth", &vars->depth );
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  */
  //-------------------------------------------------------------
  // Headers
  //

  vars->hdrId_sou_x = hdef->headerIndex( "sou_x" );
  vars->hdrId_sou_y = hdef->headerIndex( "sou_y" );
  vars->hdrId_sou_z = hdef->headerIndex( "sou_z" );

  vars->hdrId_rec_x = hdef->headerIndex( "rec_x" );
  vars->hdrId_rec_y = hdef->headerIndex( "rec_y" );
  vars->hdrId_rec_z = hdef->headerIndex( "rec_z" );

  if( !hdef->headerExists( "bin_x" ) ) {
    hdef->addHeader( TYPE_DOUBLE, "bin_x", "Bin X coordinate [m]");
  }
  else if( hdef->headerType( "bin_x" ) != TYPE_DOUBLE ) {
    writer->error("Trace header 'bin_x' exists but has the wrong number type. Should be DOUBLE.");
  }
  if( !hdef->headerExists( "bin_y" ) ) {
    hdef->addHeader( TYPE_DOUBLE, "bin_y", "Bin Y coordinate [m]");
  }
  else if( hdef->headerType( "bin_y" ) != TYPE_DOUBLE ) {
    writer->error("Trace header 'bin_y' exists but has the wrong number type. Should be DOUBLE.");
  }
  if( !hdef->headerExists( "bin_z" ) ) {
    hdef->addHeader( TYPE_DOUBLE, "bin_z", "Bin Z coordinate [m]");
  }
  else if( hdef->headerType( "bin_z" ) != TYPE_DOUBLE ) {
    writer->error("Trace header 'bin_z' exists but has the wrong number type. Should be DOUBLE.");
  }
  vars->hdrId_bin_x = hdef->headerIndex("bin_x");
  vars->hdrId_bin_y = hdef->headerIndex("bin_y");
  vars->hdrId_bin_z = hdef->headerIndex("bin_z");

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_mirror_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;

  csTrace* trace = traceGather->trace(0);


  csPoint3D p_src;
  csPoint3D p_rcv;
  csTraceHeader* trcHdr = trace->getTraceHeader();
  p_src.x = trcHdr->doubleValue( vars->hdrId_sou_x );
  p_src.y = trcHdr->doubleValue( vars->hdrId_sou_y );
  p_src.z = trcHdr->doubleValue( vars->hdrId_sou_z );
  p_rcv.x = trcHdr->doubleValue( vars->hdrId_rec_x );
  p_rcv.y = trcHdr->doubleValue( vars->hdrId_rec_y );
  p_rcv.z = trcHdr->doubleValue( vars->hdrId_rec_z );
  csProfile profile(p_src,p_rcv, csProfile::EXTRAPOLATION_CONSTANT);

  //--------------------------------------------------
  // Extract bathymetry values
  //
  if( vars->table != NULL ) {
    int numDepths = (int)( profile.offset() / vars->resolution ) + 1;
    if( numDepths < 2 ) numDepths = 2;
    double keys[2];
    profile.initDepths( numDepths );
    for( int i = 0; i < numDepths; i++ ) {
      double lamda = (double)i / (double)(numDepths-1);
      csPoint3D p = profile.pointAt( lamda );
      keys[0] = p.x;
      keys[1] = p.y;
      double depth = vars->table->getValue( keys );
      profile.setDepth( i, depth );
      //      fprintf(stdout,"%f %f %f\n", p.x, p.y, depth);
    }
  }
  else {
    profile.initDepths( 1, &vars->constantDepth );
  }

  //------------------------------------------------
  // Two-point 'ray-tracing'
  //
  double lamdaLeft  = 0.1;
  double lamdaRight = 0.5;
  double lamdaMid   = 1.0/3.0;
  double depthMid   = tryLamda( lamdaMid, p_src, &profile );
  double minDepthDiff= 2*p_src.z;
  if( vars->threshold > 0.0 ) minDepthDiff = vars->threshold;
  double depthDiff = depthMid-p_rcv.z;

  int counter = 1;
  while( fabs(depthDiff) > minDepthDiff && counter < 30 ) {
    if( depthDiff < 0.0 ) {
      lamdaRight = lamdaMid;
    }
    else {
      lamdaLeft = lamdaMid;
    }
    lamdaMid = 0.5 * (lamdaLeft + lamdaRight);
    depthMid = tryLamda( lamdaMid, p_src, &profile );
    depthDiff = depthMid-p_rcv.z;
    counter += 1;
  }

  csPoint3D p = profile.pointAt( lamdaMid );
  trcHdr->setDoubleValue( vars->hdrId_bin_x, p.x );
  trcHdr->setDoubleValue( vars->hdrId_bin_y, p.y );
  trcHdr->setDoubleValue( vars->hdrId_bin_z, p.z );

  //double keys[2];
  //keys[0] = p_src.x;
  //keys[1] = p_src.y;

  // int source = trcHdr->intValue( hdef->headerIndex("source") );
  // tryLamda( lamdaMid, p_src, &profile, true );
  // fprintf(stderr,"%d %f %f %f  0.0   %f\n", source, p_src.x, p_src.y, p_src.z, vars->table->getValue(keys) );
  // fprintf(stderr,"%d %f %f %f  %f\n", source, p.x, p.y, p.z, lamdaMid);
  // fprintf(stderr,"%d %f %f %f  1.0\n\n", source, p_rcv.x, p_rcv.y, p_rcv.z);
  //  for( int i = 0; i <= 100; i++ ) {
  //   double lamda = (double)i / 100;
  //  fprintf(stderr,"%f %f\n", lamda, profile.depthAt(lamda) );
  // }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_mirror_( csParamDef* pdef ) {
  pdef->setModule( "MIRROR", "Perform mirror image binning", "Sets trace headers bin_x, bin_y, and bin_z" );

  pdef->addDoc("NOTE: The input bathymetry table, if specified, must contain three columns giving XYZ coordinates.");
  pdef->addDoc("The bathymetry table must be given as a regular 2D grid, in Seaseis table format: First line must define @x @y depth (two table 'keys' x and y, and one 'value', i.e. depth)");

  pdef->addParam( "method", "Binning method", NUM_VALUES_FIXED );
  pdef->addValue( "obs_mirror", VALTYPE_OPTION );
  pdef->addOption( "obs_mirror", "OBS mirror image binning: Set bin centre at seabed reflection point of first water bottom multiple." );

  pdef->addParam( "table", "Table containing water depths (bathymetry)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name of Seaseis table",
                  "The table must have at least three columns, the first two giving X and Y coordinates [m], and the third giving the depth [m]. The bathymetry table must be binned to a regular grid in North-South direction" );

  pdef->addParam( "const_depth", "Use constant water depth", NUM_VALUES_FIXED, "If a constant water depth is specified, all parameters relating to an input bathymetry table are not used" );
  pdef->addValue( "", VALTYPE_NUMBER, "Constant water depth [m]");

  pdef->addParam( "resolution", "Resolution of 2D bathymetry profile used in 'ray tracing' algorithm", NUM_VALUES_FIXED );
  pdef->addValue( "500", VALTYPE_NUMBER, "Distance between knee points extracted from bathymetry table [m]");

  pdef->addParam( "target_diff", "Termination threshold for 'ray tracing'. Iterative search terminates when a ray is found that hits the target point within the specified accuracy",
                  NUM_VALUES_FIXED, "When not specified, 2x source depth is used as the termination threshold" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum distance to target point [m]",
                  "Iteration ends when ray is found that hits the target point within the specified distance");
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_mirror_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_mirror::VariableStruct* vars = reinterpret_cast<mod_mirror::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_mirror_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_mirror::VariableStruct* vars = reinterpret_cast<mod_mirror::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->table != NULL ) {
    delete vars->table;
    vars->table = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_mirror_( csParamDef* pdef ) {
  params_mod_mirror_( pdef );
}
extern "C" void _init_mod_mirror_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_mirror_( param, env, writer );
}
extern "C" bool _start_exec_mod_mirror_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_mirror_( env, writer );
}
extern "C" void _exec_mod_mirror_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_mirror_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_mirror_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_mirror_( env, writer );
}

//--------------------------------------------------------------------------------
//

double tryLamda( double lamda, csPoint3D const& p_src, csProfile const* profile, bool doPrint ) {
  double depth     = profile->depthAt( lamda );
  double offset    = profile->offsetAt( lamda );
  double an_depart = atan2( offset, (depth-p_src.z) );

  // Gradient:
  // double an_gradient = 0.0;
  // double an_refl   = 2*an_gradient + an_depart;
  // double tan_refl  = tan(an_refl);
  //   Gradient = 0  -->  an_refl = an_depart

  double tan_refl  = tan(an_depart);
  double lamda2    = lamda + (tan_refl * depth )/profile->offset();
  if( lamda2 >= 1.0 ) return -profile->depthAt( 1.0 );
  double delta_s   = (1.0 - lamda2) * profile->offset();
  double depth_trial = delta_s / tan_refl;

  if( doPrint ) {
    fprintf(stdout,"0.0 %f\n", p_src.z);
    fprintf(stdout,"%f %f\n", lamda, depth);
    fprintf(stdout,"%f 0.0\n", lamda2);
    fprintf(stdout,"1.0 %f\n\n", depth_trial);
  }
  
  return depth_trial;
}

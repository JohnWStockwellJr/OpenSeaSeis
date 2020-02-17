/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: GEOTOOLS
 * NOTE: This module is a prototype used for testing only
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_geotools {
  struct VariableStruct {
    int tool;
    int method;
    int hdrId_out;
    type_t hdrType_out;
    int hdrId_an_inci;

    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_sou_z;
    int hdrId_rec_x;
    int hdrId_rec_y;
    int hdrId_rec_elev;

    int nvalues;
    int* stations;
    double* xvalues;
    double* yvalues;
    double* zvalues;
    double mean_value;

    float accuracy;
    float velocity;
  };
}
using mod_geotools::VariableStruct;

static int const TOOL_OBC_MULT1 = 1;
static int const SEABED_FLAT         = 230;
static int const SEABED_LINEAR_SLOPE = 231;
static int const MAX_ITERATIONS      = 1000;

double compute_dist_mult1( double angle, double aoffset, double srcz, double rcvz, double* total_dist );
double compute_dist_mult1_slope( double angle, double aoffset, double srcz, double rcvz, double wdep_src, double* total_dist, double* angle_of_incidence );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_geotools_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


//---------------------------------------------
//
//
  std::string text;

  //  vars->method = SEABED_LINEAR_SLOPE;
  vars->method = SEABED_FLAT;

  param->getString( "method", &text );
  if( !text.compare("seabed_flat") ) {
    vars->method = SEABED_FLAT;
  }
  else if( !text.compare("seabed_linear_slope") ) {
    vars->method = SEABED_LINEAR_SLOPE;
  }
  else {
    writer->error("Unknown tool method: '%s'", text.c_str() );
  }

  //----------------------------------------------------------------
  
  param->getString( "tool", &text );
  if( !text.compare("comp_obc_mult1") ) {
    vars->tool = TOOL_OBC_MULT1;
  }
  else {
    writer->error("Unknown tool option: '%s'", text.c_str() );
  }

  param->getFloat( "velocity", &vars->velocity );

  if( param->exists("accuracy") ) {
    param->getFloat( "accuracy", &vars->accuracy );
  }
  else {
    vars->accuracy = 10.0;
  }

  //----------------------------------------------------------------

  param->getString( "header", &text );
  if( hdef->headerExists( text ) ) {
    vars->hdrId_out   = hdef->headerIndex( text );
    vars->hdrType_out = hdef->headerType( text );
  }
  else {
    writer->error("Trace header does not exist: '%s'", text.c_str() );
  }

  if( !hdef->headerExists( "an_inci" ) ) {
    hdef->addStandardHeader( "an_inci" );
  }
  vars->hdrId_an_inci = hdef->headerIndex( "an_inci" );

  //---------------------------------------------------------------
  // Read in ASCII file
  //
  FILE* f_in;
  csVector<std::string> valueList;

  param->getString( "filename", &text );
  if( (f_in = fopen( text.c_str(), "r" )) == (FILE*) NULL ) {
    writer->error("Error occurred when opening file '%s'", text.c_str() );
  }

  char buffer[1024];
  int counterLines = 0;
  while( fgets( buffer, 1024, f_in ) != NULL ) {
    counterLines++;
  }
  vars->nvalues = counterLines;
  rewind( f_in );

  counterLines = 0;
  vars->xvalues  = new double[vars->nvalues];
  vars->yvalues  = new double[vars->nvalues];
  vars->zvalues  = new double[vars->nvalues];
  vars->stations = new int[vars->nvalues];
  double sum = 0.0;
  
  while( fgets( buffer, 1024, f_in ) != NULL ) {
    valueList.clear();
    tokenize( buffer, valueList );
    if( valueList.size() < 4 ) {
      writer->error("Input file contains too few columns...");
    }
    vars->stations[counterLines]  = atoi(valueList.at(0).c_str());
    vars->xvalues[counterLines]   = atof(valueList.at(1).c_str());
    vars->yvalues[counterLines]   = atof(valueList.at(2).c_str());
    vars->zvalues[counterLines]   = atof(valueList.at(3).c_str());
    sum += vars->zvalues[counterLines];
    //    fprintf(stdout,"--- %d %f\n", vars->stations[counterLines], vars->values[counterLines]);
    counterLines++;
  }

  fclose( f_in );

  vars->mean_value = sum/vars->nvalues;

//----------------------------------------------
  vars->hdrId_sou_x = hdef->headerIndex("sou_x");
  vars->hdrId_sou_y = hdef->headerIndex("sou_y");
  vars->hdrId_sou_z = hdef->headerIndex("sou_z");
  vars->hdrId_rec_x = hdef->headerIndex("rec_x");
  vars->hdrId_rec_y = hdef->headerIndex("rec_y");
  vars->hdrId_rec_elev = hdef->headerIndex("rec_elev");
  
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_geotools_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
//  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);

  if( vars->stations ) {
    delete [] vars->stations;
    vars->stations = NULL;
  }

  csTraceHeader* trcHdr = trace->getTraceHeader();
  
  if( vars->tool == TOOL_OBC_MULT1 ) {
    float time = 0.0;
    double total_dist;
    
    double srcx = trcHdr->doubleValue( vars->hdrId_sou_x );
    double srcy = trcHdr->doubleValue( vars->hdrId_sou_y );
    double srcz = trcHdr->doubleValue( vars->hdrId_sou_z );
    double rcvx = trcHdr->doubleValue( vars->hdrId_rec_x );
    double rcvy = trcHdr->doubleValue( vars->hdrId_rec_y );
    double rcvz = trcHdr->doubleValue( vars->hdrId_rec_elev );

    double dx = rcvx-srcx;
    double dy = rcvy-srcy;
    double dz = abs(rcvz-srcz);
    double aoffset  = sqrt( dx*dx + dy*dy );
    
    double angle_left  = 0.0;
    double angle_right = atan2( aoffset, dz );
    double angle_mid = 0.0;

    double delta_left = 0.0;
    //    double delta_right = 0.0;
    double delta_mid = 0.0;
    double an_inci = 0.0;

    int counterIteration = 0;
    //------------------------------------
    // Assume flat seabed
    //
    if( vars->method == SEABED_FLAT ) {
      delta_left  = compute_dist_mult1( angle_left, aoffset, srcz, rcvz, &total_dist );
      //      delta_right = compute_dist_mult1( angle_right, aoffset, srcz, rcvz, &total_dist );    
      time = (float)total_dist / vars->velocity;
      do {
        angle_mid = 0.5 * ( angle_left + angle_right );
        delta_mid   = compute_dist_mult1( angle_mid, aoffset, srcz, rcvz, &total_dist );
        
        time = (float)total_dist / vars->velocity;
        //      writer->line("angle, delta, time:  %9.5f  %9.5f  %9.5f", RAD2DEG(angle_mid), delta_mid, 1000*time );
        //       writer->line("angle, delta, time:  %9.5f  %9.5f  %9.5f", delta_left, delta_mid, delta_right );
        //   printf("angle, delta, time:  %9.5f  %9.5f  %9.5f\n", delta_left, delta_mid, delta_right );
        
        if( delta_mid*delta_left > 0.0 ) {
          angle_left = angle_mid;
          delta_left   = delta_mid;
        }
        else {
          angle_right = angle_mid;
          //          delta_right  = delta_mid;
        }
        counterIteration += 1;
      } while( abs(delta_mid) > vars->accuracy && counterIteration < MAX_ITERATIONS );
      an_inci = angle_mid;
    }
    //------------------------------------
    // Assume seabed has linear slope
    //
    else if( vars->method == SEABED_LINEAR_SLOPE ) {
      int minIndex = 0;
      double minOffset = 9999999.9;
      for( int i = 0; i < vars->nvalues; i++ ) {
        float tmpoff = sqrt( pow(srcx-vars->xvalues[i],2) + pow(srcy-vars->yvalues[i],2) );
        if( tmpoff < minOffset ) {
          minOffset = tmpoff;
          minIndex = i;
        }
      }
      double wdep_src = vars->zvalues[minIndex];
      delta_left  = compute_dist_mult1_slope( angle_left, aoffset, srcz, rcvz, wdep_src, &total_dist, &an_inci );
      //      delta_right = compute_dist_mult1_slope( angle_right, aoffset, srcz, rcvz, wdep_src, &total_dist, &an_inci );
      time = (float)total_dist / vars->velocity;
      do {
        angle_mid = 0.5 * ( angle_left + angle_right );
        delta_mid = compute_dist_mult1_slope( angle_mid, aoffset, srcz, rcvz, wdep_src, &total_dist, &an_inci );
        
        time = (float)total_dist / vars->velocity;
        //      writer->line("an_inci, delta, time:  %9.5f  %9.5f  %9.5f", RAD2DEG(angle_mid), delta_mid, 1000*time );
        //       writer->line("angle, delta, time:  %9.5f  %9.5f  %9.5f", delta_left, delta_mid, delta_right );
        //   printf("angle, delta, time:  %9.5f  %9.5f  %9.5f\n", delta_left, delta_mid, delta_right );
        
        if( delta_mid*delta_left > 0.0 ) {
          angle_left = angle_mid;
          delta_left   = delta_mid;
        }
        else {
          angle_right = angle_mid;
          //          delta_right  = delta_mid;
        }
        counterIteration += 1;
      } while( abs(delta_mid) > vars->accuracy && counterIteration < MAX_ITERATIONS);
    }
    if( counterIteration == MAX_ITERATIONS ) {
      writer->warning("Maximum number of iterations reached!");
    }
    writer->line("an_inci, delta, time, iterations:  %7.3f %7.3f   %9.5f  %9.5f  %6d", RAD2DEG(an_inci), RAD2DEG(angle_mid), delta_mid, 1000*time, counterIteration );
    //    printf("an_inci, delta, time, iterations:  %7.3f %7.3f   %9.5f  %9.5f  %6d\n", RAD2DEG(an_inci), RAD2DEG(angle_mid-an_inci)/2, delta_mid, 1000*time, counterIteration );

    trcHdr->setFloatValue( vars->hdrId_out, (float)time*1000 );
    trcHdr->setFloatValue( vars->hdrId_an_inci, (float)RAD2DEG(an_inci) );
  }
  
  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_geotools_( csParamDef* pdef ) {
  pdef->setModule( "GEOTOOLS", "Various geophysical tools", "This module is still a prototype for testing only." );

  pdef->addParam( "tool", "Tool", NUM_VALUES_FIXED );
  pdef->addValue( "comp_obc_mult1", VALTYPE_OPTION );
  pdef->addOption( "comp_obc_mult1", "Compute time of first water bottom multiple for OBC geometry" );

  pdef->addParam( "method", "Tool method", NUM_VALUES_FIXED );
  pdef->addValue( "seabed_flat", VALTYPE_OPTION );
  pdef->addOption( "seabed_flat", "Assume flat seabed" );
  pdef->addOption( "seabed_linear_slope", "Assume seabed has linear slope" );

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name", "...containing receiver depths" );

  pdef->addParam( "header", "Trace header name to store result", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name", "...to store computed travel time" );

  pdef->addParam( "velocity", "Velocity [m/s]", NUM_VALUES_FIXED );
  pdef->addValue( "1500.0", VALTYPE_NUMBER, "Velocity [m/s]" );

  pdef->addParam( "accuracy", "Accuracy [m]", NUM_VALUES_FIXED );
  pdef->addValue( "10.0", VALTYPE_NUMBER, "Accuracy [m]" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_geotools_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_geotools::VariableStruct* vars = reinterpret_cast<mod_geotools::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_geotools_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_geotools::VariableStruct* vars = reinterpret_cast<mod_geotools::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars ) {
    if( vars->xvalues ) delete [] vars->xvalues;
    if( vars->yvalues ) delete [] vars->yvalues;
    if( vars->zvalues ) delete [] vars->zvalues;
    if( vars->stations ) delete [] vars->stations;
    delete vars; vars = NULL;
  }
}

extern "C" void _params_mod_geotools_( csParamDef* pdef ) {
  params_mod_geotools_( pdef );
}
extern "C" void _init_mod_geotools_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_geotools_( param, env, writer );
}
extern "C" bool _start_exec_mod_geotools_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_geotools_( env, writer );
}
extern "C" void _exec_mod_geotools_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_geotools_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_geotools_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_geotools_( env, writer );
}

//*************************************************************************************************
//
//
double compute_dist_mult1( double angle, double aoffset, double srcz, double rcvz, double* total_dist ) {
  double dz = abs( rcvz - srcz);
  double x1 = tan(angle) * dz;
  double x2 = tan(angle) * rcvz;
  
  double total_off = x1 + 2*x2;
  *total_dist = sqrt(x1*x1+dz*dz) + 2*sqrt(rcvz*rcvz + x2*x2);
  
  return ( aoffset - total_off );
}
//*************************************************************************************************
//
// 
//
//
//
//
double compute_dist_mult1_slope(
                              double angle,                 // Angle of emergence of ray from source
                              double aoffset,               // Absolute offset between true source and receiver position
                              double srcz,                  // Source depth
                              double rcvz,                  // Receiver depth (=seabed depth at receiver station)
                              double wdep_src,              // Water depth vertically beneath source position [m]
                              double* total_dist,           // Total distance of multiple ray path
                              double* angle_of_incidence )  // Angle of incidence at receiver, relative to true vertical
{
  double an_inci;

  // Damping: Needed to damp down high frequency jumps in water depth table
  // Damping term flattens seabed for near offset traces, and leaves measured slope untouched for far offset traces
  double damping = aoffset/wdep_src;
  if( damping > 1.0 ) damping = 1.0;
  
  double slope = atan2( wdep_src-rcvz, aoffset ) * damping;
  double lambda = (wdep_src - srcz) / ( cos(angle) + sin(angle)*tan(slope) );

  // p: Direction vector of ray
  double px = -sin(angle);
  double pz = -cos(angle);

  // Q: Reflection point at seabed
  double Qx = 0.0  + lambda * sin(angle);
  double Qz = srcz + lambda * cos(angle);
  // Incidence angle at seabed
  an_inci = angle - slope;

  // p2: Direction vector of ray
  double p2x = px * cos(2*an_inci) - pz * sin(2*an_inci);
  double p2z = px * sin(2*an_inci) + pz * cos(2*an_inci);

  // R: Reflection point at sea surface
  double Rx = Qx - Qz * (p2x/p2z);
  double Rz = 0.0;

  // Incidence angle at sea surface, at point R
  an_inci = 2*an_inci - angle;
  //  an_inci = angle - 2*slope;

  // tmpz: Water depth beneath surface reflection point R
  double tmpz = rcvz + (aoffset-Rx) * tan(slope);
  lambda = tmpz / ( cos(an_inci) + sin(an_inci)*tan(slope) );

  // S: Point of impact at seabed
  double Sx = Rx + lambda * sin(an_inci);
  double Sz = Rz + lambda * cos(an_inci);

  *total_dist =
    sqrt( pow(Qx-0.0,2.0) + pow(Qz-srcz,2) ) +
    sqrt( pow(Rx-Qx,2) + pow(Rz-Qz,2) ) +
    sqrt( pow(Sx-Rx,2) + pow(Sz-Rz,2) );
  *angle_of_incidence = an_inci;
  /*  
  fprintf(stdout,"0.0 0.0  --O\n");
  fprintf(stdout,"%f %f  --Q\n", Qx, Qz);
  fprintf(stdout,"%f %f  --R, an_inci: %8.4f\n", Rx, Rz, an_inci*180.0/M_PI);
  //  fprintf(stdout,"%f %f  --D\n", Rx, tmpz);
  //  fprintf(stdout,"%f %f  --D\n", Rx, Rz);
  fprintf(stdout,"%f %f  --S, an_inci: %8.4f\n\n", Sx, Sz, an_inci*180.0/M_PI);
  */
  //  fprintf(stdout,"--slope: %7.3f\n", RAD2DEG(slope) );
  return( aoffset - Sx );
}


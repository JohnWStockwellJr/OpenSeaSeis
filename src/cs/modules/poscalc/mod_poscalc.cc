/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <assert.h>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: POSCALC
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_poscalc {
  struct VariableStruct {
    float timeDelay;
    float velocity;
    float sou_z;
    bool verbose;
    int hdrId;
    int hdrId_srcx;
    int hdrId_srcy;
    int hdrId_srcz;
    int hdrId_rcvx;
    int hdrId_rcvy;
    int hdrId_rcvdep;
    int hdrId_rcv;
    int hdrId_dx;
    int hdrId_dy;
    int hdrId_dz;
    int hdrId_dt;
    int hdrId_dv;
    int hdrId_velocity;
    double** xSrc;
    double* xDelta;   // Result: Delta correction for all unknowns
    double* xStddev;  // Result: Standard deviation for all unknowns
    double* fbpTime;  // [ms]
    int maxObs;
    int nUnknowns;
    bool solveReceiver;
    bool isConstantSrcDepth;
  };
  static int const METHOD_PEAK_TROUGH   = 11;
  static int const METHOD_ZERO_CROSSING = 22;

  void sub_poscalc_indv( double timeDelay, // Time delay to apply to FB picks in advance, default=0.0
                         bool verbose,
                         int nObs, // Number of observations
                         double velocity,
                         double rcvx,
                         double rcvy,
                         double rcvdep,
                         double** xSrc,
                         double* fbpTime,
                         int nUnknowns,
                         double* xDelta,
                         double* xStddev );
  void LUDecomposition( double** AA, int nIn, int* indx, double* dd );
  void LUBackSub( double** AA, int nIn, int* indx, double* bb );
}
using namespace mod_poscalc;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_poscalc_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );


    vars->timeDelay   = 0;
    vars->velocity    = 0;
    vars->sou_z      = 0;
    vars->verbose     = false;
    vars->hdrId        = -1;
    vars->hdrId_srcx   = -1;
    vars->hdrId_srcy   = -1;
    vars->hdrId_srcz   = -1;
    vars->hdrId_rcvx   = -1;
    vars->hdrId_rcvy   = -1;
    vars->hdrId_rcvdep = -1;
    vars->hdrId_rcv    = -1;
    vars->hdrId_dx     = -1;
    vars->hdrId_dy     = -1;
    vars->hdrId_dz     = -1;;
    vars->hdrId_dt     = -1;
    vars->hdrId_dv     = -1;
    vars->hdrId_velocity = -1;

    vars->xSrc      = NULL;
    vars->xDelta = NULL;
    vars->xStddev   = NULL;
    vars->fbpTime   = NULL;
    vars->maxObs    = 0;
    vars->nUnknowns = 0;
    vars->solveReceiver = true;
    vars->isConstantSrcDepth = false;

  //---------------------------------------------------------
  std::string text;

  if( param->exists("verbose") ) {
    param->getString( "verbose", &text );
    if( !text.compare("yes") ) {
      vars->verbose = true;
    }
    else if( !text.compare("no") ) {
      vars->verbose = false;
    }
    else {
      writer->line("Unknown option '%s'", text.c_str() );
    }
  }

  if( param->exists("domain") ) {
    param->getString( "domain", &text );
    if( !text.compare("receiver") ) {
      vars->solveReceiver = true;
    }
    else if( !text.compare("source") ) {
      vars->solveReceiver = false;
    }
    else {
      writer->line("Unknown option '%s'", text.c_str() );
    }
  }

  if( param->exists("timedelay") ) param->getFloat( "timedelay", &vars->timeDelay );

  //------------------------
  param->getString( "velocity", &text );
  csFlexNumber number;
  if( !number.convertToNumber( text ) ) {
    if( !hdef->headerExists(text) ) {
      writer->error("Specified velocity value is not an existing trace header name, nor a number: '%s'", text.c_str());
    }
    vars->hdrId_velocity = hdef->headerIndex(text);
    writer->line("Read in water velocity from trace header '%s' (read in from first trace of input ensemble)", text.c_str());
  }
  else {
    vars->velocity = number.floatValue();
    writer->line("Using constant water velocity = %f", vars->velocity);
  }
  //------------------------

  if( param->exists("sou_z") ) {
    param->getFloat( "sou_z", &vars->sou_z );
    vars->isConstantSrcDepth = true;
  }
  else if( hdef->headerExists("sou_z") ) {
    vars->hdrId_srcz = hdef->headerIndex("sou_z");
    vars->isConstantSrcDepth = false;
  }
  else {
    writer->error("Source depth must be given either in trace header sou_z, or input parameter 'sou_z'");
  }
  param->getInt( "num_unknowns", &vars->nUnknowns );

  if( param->exists( "maxobs" ) ) {
    param->getInt( "maxobs", &vars->maxObs );
  }
  std::string ensHeaderName;
  param->getString( "hdr_ens", &ensHeaderName );
  if( !hdef->headerExists( ensHeaderName ) ) {
    writer->error("Ensemble trace header '%s' does not exist", ensHeaderName.c_str() );
  }

  param->getString( "header", &text );
  if( !hdef->headerExists( text ) ) {
    writer->error("Trace header '%s' does not exist", text.c_str() );
  }
  vars->hdrId        = hdef->headerIndex( text );

  // Allocate memory
  vars->xDelta = new double[vars->nUnknowns];
  vars->xStddev   = new double[vars->nUnknowns];
  if( vars->maxObs > 0 ) {
    vars->xSrc = new double*[vars->maxObs];
    for( int i = 0; i < vars->maxObs; i++ ) {
      vars->xSrc[i] = new double[5];
    }
    vars->fbpTime = new double[vars->maxObs];
  }

  vars->hdrId_srcx   = hdef->headerIndex( "sou_x" );
  vars->hdrId_srcy   = hdef->headerIndex( "sou_y" );
  vars->hdrId_rcvx   = hdef->headerIndex( "rec_x" );
  vars->hdrId_rcvy   = hdef->headerIndex( "rec_y" );
  vars->hdrId_rcvdep = hdef->headerIndex( "rec_z" );
  vars->hdrId_rcv    = hdef->headerIndex( ensHeaderName );

  // If source side shall be solved, simply exchange source and receiver headers:
  if( !vars->solveReceiver ) {
    int tmp;

    tmp = vars->hdrId_srcx;
    vars->hdrId_srcx = vars->hdrId_rcvx;
    vars->hdrId_rcvx = tmp;

    tmp = vars->hdrId_srcy;
    vars->hdrId_srcy = vars->hdrId_rcvy;
    vars->hdrId_rcvy = tmp;

    if( vars->isConstantSrcDepth ) {
      //      tmp = vars->hdrId_srcz;
      //      vars->hdrId_rcvdep = tmp;
      vars->hdrId_srcz = vars->hdrId_rcvdep;
    }
    else {
      tmp = vars->hdrId_srcz;
      vars->hdrId_srcz = vars->hdrId_rcvdep;
      vars->hdrId_rcvdep = tmp;
    }
  }

  // Output headers
  if( !hdef->headerExists("dx") ) hdef->addHeader( TYPE_FLOAT, "dx" );
  if( !hdef->headerExists("dy") ) hdef->addHeader( TYPE_FLOAT, "dy" );
  if( !hdef->headerExists("dz") ) hdef->addHeader( TYPE_FLOAT, "dz" );
  vars->hdrId_dx   = hdef->headerIndex( "dx" );
  vars->hdrId_dy   = hdef->headerIndex( "dy" );
  vars->hdrId_dz   = hdef->headerIndex( "dz" );
  if( hdef->headerType("dx") != TYPE_FLOAT && hdef->headerType("dx") != TYPE_DOUBLE ) {
    writer->error("Trace header 'dx' exists but has wrong header type. Should be FLOAT.");
  }
  if( hdef->headerType("dy") != TYPE_FLOAT && hdef->headerType("dy") != TYPE_DOUBLE ) {
    writer->error("Trace header 'dy' exists but has wrong header type. Should be FLOAT.");
  }
  if( hdef->headerType("dz") != TYPE_FLOAT && hdef->headerType("dz") != TYPE_DOUBLE ) {
    writer->error("Trace header 'dz' exists but has wrong header type. Should be FLOAT.");
  }

  if( vars->nUnknowns > 3 ) {
    if( !hdef->headerExists("dt") ) {
      hdef->addHeader( TYPE_FLOAT, "dt" );
    }
    else if( hdef->headerType("dt") != TYPE_FLOAT && hdef->headerType("dt") != TYPE_DOUBLE ) {
      writer->error("Trace header 'dt' exists but has wrong header type. Should be FLOAT.");
    }
    vars->hdrId_dt   = hdef->headerIndex( "dt" );
  }
  if( vars->nUnknowns > 4 ) {
    if( !hdef->headerExists("dv") ) {
      hdef->addHeader( TYPE_FLOAT, "dv" );
    }
    else if( hdef->headerType("dv") != TYPE_FLOAT && hdef->headerType("dv") != TYPE_DOUBLE ) {
      writer->error("Trace header 'dv' exists but has wrong header type. Should be FLOAT.");
    }
    vars->hdrId_dv = hdef->headerIndex( "dv" );
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_poscalc_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;


  int nTraces = traceGather->numTraces();

  if( nTraces >= vars->maxObs ) {
    if( vars->maxObs > 0 )
      writer->line("Info: More input traces (=%d) than specified maxobs (= %d). Enlarge arrays....", nTraces, vars->maxObs);
    if( vars->xSrc != NULL ) {
      for( int iobs = 0; iobs < vars->maxObs; iobs++ ) {
        if( vars->xSrc[iobs] != NULL ) {
          delete [] vars->xSrc[iobs];
        }
      }
      delete [] vars->xSrc;
    }
    if( vars->fbpTime != NULL ) {
      delete [] vars->fbpTime;
    }

    vars->maxObs = nTraces;

    vars->xSrc = new double*[vars->maxObs];
    for( int i = 0; i < vars->maxObs; i++ ) {
      vars->xSrc[i] = new double[5];
    }
    vars->fbpTime = new double[vars->maxObs];
  }


  for( int iobs = 0; iobs < nTraces; iobs++ ) {
    csTraceHeader* trcHdr = traceGather->trace(iobs)->getTraceHeader();
    vars->fbpTime[iobs] = trcHdr->doubleValue( vars->hdrId );  // convert to [ms]
    vars->xSrc[iobs][0] = trcHdr->doubleValue( vars->hdrId_srcx );
    vars->xSrc[iobs][1] = trcHdr->doubleValue( vars->hdrId_srcy );
    if( vars->solveReceiver && vars->isConstantSrcDepth ) {
      vars->xSrc[iobs][2] = vars->sou_z;
    }
    else {
      vars->xSrc[iobs][2] = trcHdr->doubleValue( vars->hdrId_srcz );
    }
    vars->xSrc[iobs][3] = 0.0;
  }
  csTraceHeader* trcHdr = traceGather->trace(0)->getTraceHeader();
  double rcvx   = trcHdr->doubleValue( vars->hdrId_rcvx );
  double rcvy   = trcHdr->doubleValue( vars->hdrId_rcvy );
  double rcvdep;
  if( vars->solveReceiver || !vars->isConstantSrcDepth ) {
    rcvdep = trcHdr->doubleValue( vars->hdrId_rcvdep );
  }
  else {
    rcvdep = vars->sou_z;
  }

  int rcv = trcHdr->intValue( vars->hdrId_rcv );
  if( edef->isDebug() ) writer->line("Ensemble %d XYZ:  %f %f %f", rcv, rcvx, rcvy, rcvdep );

  // In case velocity trace header was given: Read velocity from first input trace
  if( vars->hdrId_velocity >= 0 ) {
    vars->velocity = traceGather->trace(0)->getTraceHeader()->floatValue( vars->hdrId_velocity );
    writer->line("POSCALC: Read in water velocity = %f [m/s] from trace header of first trace, ensemble header = %d",
              vars->velocity, rcv);
  }
  /*
  for( int iobs = 0; iobs < nTraces; iobs++ ) {
    double dx = rcvx - vars->xSrc[iobs][0];
    double dy = rcvy - vars->xSrc[iobs][1];
    double dz = rcvdep - vars->xSrc[iobs][2];
    double dist = sqrt( dx*dx + dy*dy + dz*dz );
    vars->fbpTime[iobs] = dist/vars->velocity*1000 + 10;
  }
  */
  sub_poscalc_indv(
                   vars->timeDelay,
                   vars->verbose,
                   nTraces, // Number of observations
                   vars->velocity,
                   rcvx,
                   rcvy,
                   rcvdep,
                   vars->xSrc,
                   vars->fbpTime,
                   vars->nUnknowns,
                   vars->xDelta,
                   vars->xStddev );

  if( edef->isDebug() ) {
    for( int i = 0; i < nTraces; i++ ) {
      writer->line( "%d  %f  %f %f %f\n", i, vars->fbpTime[i], vars->xSrc[i][0], vars->xSrc[i][1], vars->xSrc[i][2] );
    }
  }


  // Temp: Set array elements for printout
  vars->xSrc[0][0] = rcvx;
  vars->xSrc[1][0] = rcvy;
  if( vars->nUnknowns > 2 ) {
    vars->xSrc[2][0] = rcvdep;
    if( vars->nUnknowns > 3 ) {
      vars->xSrc[3][0] = vars->timeDelay;
      if( vars->nUnknowns > 4 ) {
        vars->xSrc[4][0] = vars->velocity;
      }
    }
  }
  

  for( int k = 0; k < vars->nUnknowns; k++) {
    writer->write("%d  %d  %-11.3f %-11.3f %-11.3f (Solution for unknown %d (x,dx,stdev))\n",
               rcv,k,vars->xSrc[k][0],vars->xDelta[k],vars->xStddev[k], k);
    //    fprintf(stdout,"%d  %d  %-11.3f %-11.3f %-11.3f (Solution for unknown %d (x,dx,stdev))\n",
    //        rcv,k,vars->xSrc[k][0],vars->xDelta[k],vars->xStddev[k], k);
  }

  // Set output headers
  for( int iobs = 0; iobs < nTraces; iobs++ ) {
    csTraceHeader* trcHdr = traceGather->trace(iobs)->getTraceHeader();
    trcHdr->setFloatValue( vars->hdrId_dx, vars->xDelta[0]);
    trcHdr->setFloatValue( vars->hdrId_dy, vars->xDelta[1]);
    trcHdr->setFloatValue( vars->hdrId_dz, vars->xDelta[2]);
    if( vars->nUnknowns > 3 ) {
      trcHdr->setFloatValue( vars->hdrId_dt, vars->xDelta[3] );
      if( vars->nUnknowns > 4 ) {
        trcHdr->setFloatValue( vars->hdrId_dv, vars->xDelta[4] );
      }
    }
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_poscalc_( csParamDef* pdef ) {
  pdef->setModule( "POSCALC", "Compute receiver position from time picks" );

  pdef->addParam( "verbose", "Print diagnostic output", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "yes" );
  pdef->addOption( "no", "No" );

  pdef->addParam( "timedelay", "Known time delay on input traces", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Constant time delay" );

  pdef->addParam( "num_unknowns", "Number of unknowns to solve for", NUM_VALUES_FIXED );
  pdef->addValue( "3", VALTYPE_NUMBER, "Number of unknowns in inversion: X(1), Y(2), Z(3), time shift(4), velocity(5)" );

  pdef->addParam( "velocity", "Water velocity", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Water velocity [m/s]", "Alternatively, specify trace header name containing water velocity. Velocity will be read in from first input trace" );

  pdef->addParam( "maxobs", "Maximum number of observations (=traces) in ensemble", NUM_VALUES_FIXED, "Required for pre-allocation purposes" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum number of observations (=traces) in ensemble" );

  pdef->addParam( "sou_z", "Source depth [m]", NUM_VALUES_FIXED, "If source depth is not specified it will be read in from trace header 'sou_z'" );
  pdef->addValue( "", VALTYPE_NUMBER, "Source depth [m]" );

  pdef->addParam( "header", "Header where picked time is stored", NUM_VALUES_FIXED );
  pdef->addValue( "time_pick", VALTYPE_STRING, "Name of trace header where picked time is stored" );

  pdef->addParam( "hdr_ens", "Ensemble header (to be used for printout)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Ensemble header (to be used for printout)" );

  pdef->addParam( "domain", "Specify which side shall be solved: Receiver or source domain", NUM_VALUES_FIXED );
  pdef->addValue( "receiver", VALTYPE_OPTION );
  pdef->addOption( "receiver", "Solve for receiver XYZ" );
  pdef->addOption( "source", "Solve for source XYZ" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_poscalc_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_poscalc::VariableStruct* vars = reinterpret_cast<mod_poscalc::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_poscalc_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_poscalc::VariableStruct* vars = reinterpret_cast<mod_poscalc::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->xSrc != NULL ) {
    for( int iobs = 0; iobs < vars->maxObs; iobs++ ) {
      if( vars->xSrc[iobs] != NULL ) {
        delete [] vars->xSrc[iobs];
      }
    }
    delete [] vars->xSrc;
  }
  if( vars->xDelta != NULL ) {
    delete [] vars->xDelta;
    vars->xDelta = NULL;
  }
  if( vars->xStddev != NULL ) {
    delete [] vars->xStddev;
    vars->xStddev = NULL;
  }
  if( vars->fbpTime != NULL ) {
    delete [] vars->fbpTime;
    vars->fbpTime = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_poscalc_( csParamDef* pdef ) {
  params_mod_poscalc_( pdef );
}
extern "C" void _init_mod_poscalc_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_poscalc_( param, env, writer );
}
extern "C" bool _start_exec_mod_poscalc_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_poscalc_( env, writer );
}
extern "C" void _exec_mod_poscalc_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_poscalc_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_poscalc_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_poscalc_( env, writer );
}



/*
 * Program poscalc_indv.cpp
 *
 * Author: Bjorn Olofsson
 * Compile: mload lang/6.1
 *          CC poscalc_indv.cpp lubksbBO.c ludcmpBO.c -o poscalc_indv
 * Date of last update: 27.09.2001
 *
 *
 * Calculates new positions for either receivers or shots using
 * first break picks
 *
 * The position of each individual receiver is calculated independently.
 *
 *
 * Methode: "Ausgleichung nach vermittelnden Beobachtungen"
 *
 ****  L_i + v_i = f_i(xHat_j)
 ****            = f_i(xStar_j + delta_xHat_j)
 ****
 **** Taylor expansion, linearising:
 ****  v_i = (df_i / dxHat_j)_0 * delta_xHat_j - l_i  /  l_i = L_i - f_i(x_Hat_j)
 *
 * L_i       : Observations (supposed to be correct in a least-square's sense)
 * v_i       : Corrections (will be minimised)
 *  (f_i       : functional relationship)
 * xStar_j   : 'Estimated' parameters (original values for the parameters)
 * xHat_j    : 'Compensated' parameters (newly calculated values)
 * delta_xHat_j : 'Reduced' vectors (corrections to the original parameter values)
 *
 *
 *----------------------------------------------------------------------*
 * Implementation of the above problem
 *
 * L_i     : First break pick time  --> fbpTime[i]
 * xStar_j : Original receiver (x,y,z) position and time delay --> xStar[j], (j=1..4)
 * xHat_j  : Newly calculated receiver position and time delay --> xHat[j],  (j=1..4)
 * l_i --> ll[i]
 * (df_i / dxHat_j)_0 --> AA[i][j]
 *
 *
 *----------------------------------------------------------------------*
 *
 * Input file format:

 #  1117                               (Unique receiver station number)
 360.0 1485.0 8.0  546980.3 6887101.0  (rcv depth, water velocity, src depth, rcvx, rcvy)
 546703.6    6887542.0    459.500      (srcx, srcy, first break time) (1)
 546703.7    6887266.5    353.500      (srcx, srcy, first break time) (2)
 546703.7    6887291.5    360.500      (srcx, srcy, first break time) (3)
 546704.1    6887492.0    437.000      (srcx, srcy, first break time) (4)
 ...
 #  1118                               (Unique receiver station number)
 359.0 1485.0 8.0  547005.3 6887101.0  (rcv depth, water velocity, src depth, rcvx, rcvy)
 ...

 *
 */

namespace mod_poscalc {
void sub_poscalc_indv( double timeDelay, // Time delay to apply to FB picks in advance, default=0.0
                       bool verbose,
                       int nObs, // Number of observations
                       double velocity,
                       double rcvx,
                       double rcvy,
                       double rcvdep,
                       double** xSrc,
                       double* fbpTime,
                       int nUnknowns,
                       double* xDelta,
                       double* xStddev )
{

  //  const unsigned short int nUnknowns = 5; // Unknowns : rcvx, rcvy, rcvdep, timeDelay/velocity
  const unsigned short int nStations = 1; // Stations : 1 receiver only
  
  double ** AA; // Matrix A
  double * ll; // Vector l_i
  //  double * ssStar; // Vectors S*i, approximation / distances between shots and receiver, by positions
  //  double ssx, ssy, ssz; // Help variables...
  //  double x, y, z, t; // Help variables
  //  double distance; // Maximum distance to process, distance calc. from approximate positions
  double ** NN; // Coefficient matrix N
  double ** NNinv; // Inverse of NN
  double * nn; // normal vector...
  double * xStar; // Approximate unknowns / original receiver position and time delay
  double * xHat; // Compensated position / new receiver position
  double * delta_xHat; // 'Reduced' positions
  double * sigma_xx; // Covariance matrix --> only diagonal elements!
  double sigma2; // Variance
  double * vv; // Correction vector

  double * srcRcvDist; // Distance between rcv & src (offset and depth accounted for)
  double left_side, right_side, ss;

  int * indxTmp; // Only needed for LU programs
  double dTmp, * colTmp; // Only needed for LU programs
  double ** LUmatrix;

  AA = NULL;
  ll = NULL;
  NN = NNinv = NULL;

  //--------------------------------------------------------------------------------
  // Command line entries
  xStar = new double [nUnknowns];
  assert(xStar != NULL);
  srcRcvDist = new double[nObs];
  assert( srcRcvDist != NULL );

  //  velocity / 1000.0 in order to convert from [m/s] into [m/ms]
  velocity /= 1000.0;
  
  xStar[0] = rcvx; // X
  xStar[1] = rcvy; // Y
  if( nUnknowns > 2 ) {
    xStar[2] = rcvdep; // Z, dept
    if( nUnknowns > 3 ) {
      xStar[3] = timeDelay;  // use variable 'timeDelay' instead of xStar[3] where possible
      if( nUnknowns > 4 ) {
        xStar[4] = velocity; // [m/ms]]
      }
    }
  }
  
  double sum;
  for( int iobs = 0; iobs < nObs; iobs++ ) {
    sum = 0.0;
    for( int j = 0; j < nUnknowns; j++ ) {
      ss = xStar[j] - xSrc[iobs][j];
      sum += ss * ss;
    }
    srcRcvDist[iobs] = sqrt(sum); // Approximate distance
  }

  //--------------------------------------------------------------------------------
  // Initialize some fields

  nn = new double [nUnknowns];
  assert(nn != NULL);
  delta_xHat = new double [nUnknowns];
  assert(delta_xHat != NULL);
  xHat = new double [nUnknowns];
  assert(xHat != NULL);
  sigma_xx = new double [nUnknowns];
  assert(sigma_xx != NULL);

  colTmp = new double [nUnknowns];
  assert(colTmp != NULL);
  indxTmp = new int [nUnknowns];
  assert(indxTmp != NULL);

  //--------------------------------------------------
  // Initialize vectors and matrizes:

  ll = new double [nObs];
  assert(ll != NULL);
  vv = new double [nObs];
  assert(vv != NULL);

  AA = new double * [nObs];
  assert(AA != NULL);
  for( int i = 0; i < nObs; i++ ) {
    AA[i] = new double [nUnknowns]; // 4 unknowns: rcvx, rcvy, rcvdep, timeDelay
    assert(AA[i] != NULL);
  }

  NN = new double * [nUnknowns];
  assert(NN != NULL);
  NNinv = new double * [nUnknowns];
  assert(NNinv != NULL);
  LUmatrix = new double * [nUnknowns];
  assert(LUmatrix != NULL);
  for( int j = 0; j < nUnknowns; j++ ) {
    NN[j] = new double [nUnknowns]; // 3 unknowns: rcvx, rcvy, rcvdep
    assert(NN[j] != NULL);
    NNinv[j] = new double [nUnknowns];
    assert(NNinv[j] != NULL);
    LUmatrix[j] = new double [nUnknowns];
    assert(LUmatrix[j] != NULL);
  }

  //--------------------------------------------------
  // Step 1: Setting initial matrizes and vectors

  for( int i = 0; i < nObs; i++ ) {
    ll[i] = fbpTime[i] - (srcRcvDist[i]/velocity + timeDelay);
    AA[i][0] = (xStar[0] - xSrc[i][0]) / (srcRcvDist[i] * velocity);
    AA[i][1] = (xStar[1] - xSrc[i][1]) / (srcRcvDist[i] * velocity);
    if( nUnknowns > 2 )
      AA[i][2] = (xStar[2] - xSrc[i][2]) / (srcRcvDist[i] * velocity);
    if( nUnknowns > 3 )
      AA[i][3] = 1.0;
    if( nUnknowns > 4 )
      AA[i][4] = srcRcvDist[i] / (velocity*velocity);

    if (verbose) {
      fprintf(stdout,"ll[%d] = %f srcrcvdist= %f\n",i, ll[i], srcRcvDist[i]);
      fprintf(stdout,"  [%d] = xStar=%f yStar=%f zStar=%f\n",i, xStar[0], xStar[1], xStar[2]);
      fprintf(stdout,"  [%d] = xSrc[0]=%f xSrc[1]=%f zSrc=%f\n",i, xSrc[i][0], xSrc[i][1], xSrc[i][2]);
    }
  }

  //  if (1) { exit(-1); }

  //--------------------------------------------------
  // Step 2: Calculate coefficient matrix NN and normal vector nn
  // NN = AA(T) * AA    / AA(T): Transposed AA
  // nn = AA(T) * ll

  for( int j = 0; j < nUnknowns; j++) {
    for( int k = 0; k < nUnknowns; k++) {
      sum = 0;
      for( int i = 0; i < nObs; i++) {
        sum += AA[i][j] * AA[i][k];
      }
      NN[j][k] = sum;
    }
  }

  for( int j = 0; j < nUnknowns; j++) {
    sum = 0;
    for( int i = 0; i < nObs; i++) {
      sum += AA[i][j] * ll[i];
    }
    nn[j] = sum;
  }

  //--------------------------------------------------
  // Step 3: Invert matrix NN --> NNinv
  //

  // Copy matrix NN into temporary matrix LUmatrix
  for( int j = 0; j < nUnknowns; j++) {
    for( int k = 0; k < nUnknowns; k++) {
      LUmatrix[j][k] = NN[j][k];
    }
  }

  LUDecomposition(LUmatrix,nUnknowns,indxTmp,&dTmp);

  for( int j = 0; j < nUnknowns; j++) {
    for( int k = 0; k < nUnknowns; k++) colTmp[k] = 0.0;
    colTmp[j] = 1.0;
    LUBackSub(LUmatrix,nUnknowns,indxTmp,colTmp);
    for( int k = 0; k < nUnknowns; k++) NNinv[k][j] = colTmp[k];
  }

  //--------------------------------------------------
  // Test paragraph only:
  // Print NN and NNinv matrix

  if( verbose ) {
    for( int i = 0; i < nUnknowns; i++) {
      for( int j = 0; j < nUnknowns; j++) {
        fprintf(stdout,"%-10.4f ",NN[i][j]);
      }
      cout << endl;
    }
    cout << endl;
    for( int i = 0; i < nUnknowns; i++) {
      for( int j = 0; j < nUnknowns; j++) {
        fprintf(stdout,"%-10.4f ",NNinv[i][j]);
      }
      cout << endl;
    }
    cout << endl;

    fprintf(stdout,"Check matrix inversion: NN * NNinv =\n");
    // Check matrix inversion
    for( int i = 0; i < nUnknowns; i++) {
      for( int j = 0; j < nUnknowns; j++) {
        sum = 0.0;
        for( int k = 0; k < nUnknowns; k++) {
          sum += NN[i][k]*NNinv[k][j];
        }
        fprintf(stdout,"%-10.4f ",sum);
        //          LUmatrix[i][j] = sum;
      }
      cout << endl;
    }
  }

  //--------------------------------------------------
  // Step 4: Calculate reduced vector delta_xHat and
  // compensated parameters (vector xHat)
  //

  for( int j = 0; j < nUnknowns; j++) {
    sum = 0.0;
    for( int k = 0; k < nUnknowns; k++) {
      sum += NNinv[j][k]*nn[k];
    }
    delta_xHat[j] = sum;
    xHat[j] = xStar[j] + delta_xHat[j];
  }

  //--------------------------------------------------
  // Step 5: Calculate correction vector vv
  //

  for( int i = 0; i < nObs; i++) {
    sum = 0.0;
    for( int k = 0; k < nUnknowns; k++) {
      sum += AA[i][k]*delta_xHat[k];
    }
    vv[i] = sum - ll[i];
  }

  //--------------------------------------------------
  // Step 6: Calculate variance matrix sigma2 and
  // covariance matrix sigma_xx
  // Don't need full matrix sigma_xx --> only diagonal elements!

  sum = 0.0;
  for( int i = 0; i < nObs; i++) {
    sum += vv[i]*vv[i];
  }
  sigma2 = sum / (nUnknowns-nStations);
  for( int k = 0; k < nUnknowns; k++) {
    sigma_xx[k] = sigma2 * NNinv[k][k];
  }

  //--------------------------------------------------
  // Step 7: Cross check results
  //  Method 1: v(T)*v = -v(T)*ll
  //  Method 2: A(T)*v = 0

  left_side = sigma2 * (nUnknowns-nStations);
  sum = 0.0;
  for( int i = 0; i < nObs; i++) {
    sum += vv[i]*ll[i];
  }
  right_side = -sum;

  if( verbose ) {
    fprintf(stdout,"Error check 1:  %-8.4f-%-8.4f = %-14.10f\n",left_side,right_side,left_side-right_side);
    for( int i = 0; i < nObs; i++) {
      sum = 0.0;
      for( int k = 0; k < nUnknowns; k++) {
        sum += AA[i][k]*vv[k];
      }
      fprintf(stdout,"Error check 2(%-2d):  0 = %-14.10f\n",i,sum);
    }
  }

  //--------------------------------------------------
  // Output solutions and standard deviation for each Unknown

  for( int k = 0; k < nUnknowns; k++) {
    xDelta[k]  = delta_xHat[k];
    xStddev[k] = sqrt(fabs(sigma_xx[k]));
  }
  // Convert velocity back to [m/s], and reverse polarity (why? only when reversed the sign is OK)
  if( nUnknowns > 4 ) {
    xDelta[4]  *= -1000;
    xStddev[4] *= 1000;
  }
}

#define TINY 1.0e-10;
#define FREE_ARG char*

void LUDecomposition( double** AA, int nIn, int* indx, double* dd )
{
  float dum, sum, temp;
  unsigned int imax = 0;
  unsigned int nIn_unsigned = (unsigned int)nIn;

  float* vv = new float[ nIn_unsigned ];
  if( !vv ) fprintf(stderr,"Allocation failure in ludcmpBO\n");

  *dd = 1.0;
  for( unsigned int i = 0; i < nIn_unsigned; i++ ) {
    float big = 0.0;
    for( unsigned int j = 0; j < nIn_unsigned; j++ ) {
      if ((temp = fabs(AA[i][j])) > big) big = temp;
    }
    if (big == 0.0) fprintf(stderr,"Singular matrix in routine ludcmpBO\n");
    vv[i] = 1.0 / big;
  }

  for( unsigned int j = 0; j < nIn_unsigned; j++ ) {
    for( unsigned int i = 0; i < j; i++ ) {
      sum = AA[i][j];
      for( unsigned int k = 0; k < i; k++) sum -= AA[i][k]*AA[k][j];
      AA[i][j] = sum;
    }
    float big = 0.0;
    for( unsigned int i = j; i < nIn_unsigned; i++ ) {
      sum = AA[i][j];
      for( unsigned int k = 0; k < j; k++ ) {
        sum -= AA[i][k]*AA[k][j];
      }
      AA[i][j] = sum;
      if( (dum = vv[i]*fabs(sum)) >= big ) {
        big = dum;
        imax = i;
      }
    }
    if( j != imax ) {
      for( unsigned int k = 0; k < nIn_unsigned; k++ ) {
        dum = AA[imax][k];
        AA[imax][k] = AA[j][k];
        AA[j][k] = dum;
      }
      *dd = -(*dd);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if( AA[j][j] == 0.0 ) AA[j][j] = 1.0e-10;
    if( j != nIn_unsigned ) {
      dum = 1.0 / AA[j][j];
      for( unsigned int i = j+1; i < nIn_unsigned; i++ ) AA[i][j] *= dum;
    }
  }

  delete [] vv;
}

void LUBackSub( double** AA, int nIn, int* indx, double* bb )
{
  int ii = -1;

  for( int i = 0; i < nIn; i++ ) {
    int ip = indx[i];
    float sum = bb[ip];
    bb[ip] = bb[i];
    if( ii != -1 ) {
      for( int j = ii; j <= i-1; j++ ) sum -= AA[i][j]*bb[j];
    }
    else if( sum ) {
      ii = i;
    }
    bb[i] = sum;
  }
  for( int i = nIn-1; i >= 0; i-- ) {
    float sum = bb[i];
    for( int j = i+1; j < nIn; j++ ) sum -= AA[i][j]*bb[j];
    bb[i] = sum/AA[i][i];
  }
}
} // END namespace

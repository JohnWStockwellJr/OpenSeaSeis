/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csGeolibUtils.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: OFF2ANGLE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_off2angle {
  struct VariableStruct {
    float angle_first;
    float angle_inc;
    float angle_last;
    float angle_width;
    float v1v2_ratio;
    int hdrId_offset;
    int hdrId_rec_x;
    int hdrId_rec_y;
    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_rec_elev;
    int hdrId_depth;
    int hdrId_angle_incidence;
    int maxTracesOut;
    float* angle_incidence;
  };
  void checkHeader( std::string headerName, csTraceHeaderDef* hdef, csLogWriter* writer, type_t checkType = TYPE_UNKNOWN ) {
    if( !hdef->headerExists(headerName) ){
      writer->error("Trace header '%s' not found.", headerName.c_str());
    }
    if( checkType != TYPE_UNKNOWN ) {
      int type = hdef->headerType(headerName);
      if( TYPE_FLOAT_DOUBLE ) {
        if( type != TYPE_FLOAT &&  type != TYPE_DOUBLE ) {
          writer->error("Trace header '%s' has the wrong type. Should be %s or %s",
            headerName.c_str() );
        }
      }
      else {
        if( type != checkType ) {
          writer->error("Trace header '%s' has the wrong type. Should be %s",
            headerName.c_str(), csGeolibUtils::typeText(type) );
        }
      }
    } // END if checkType != TYPE_UNKNOWN
  }
}
using namespace mod_off2angle;

//*************************************************************************************************
// Init phase
//
// Transforms trace gather from offset to angle of incidence
// - Only method currently available transforms into angle of incidence at seabed, according to src-rcv geometry
//
// Necessary improvement:
// - Implement angle band width != 0
//
//*************************************************************************************************
void init_mod_off2angle_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  //------------------------------------------------------
  vars->angle_first   = 0;
  vars->angle_inc     = 0;
  vars->angle_last    = 0;
  vars->angle_width   = 0;
  vars->v1v2_ratio    = 0;
  vars->hdrId_offset  = -1;
  vars->hdrId_rec_x   = -1;
  vars->hdrId_rec_y   = -1;
  vars->hdrId_sou_x   = -1;
  vars->hdrId_sou_y   = -1;
  vars->hdrId_rec_elev = -1;
  vars->hdrId_depth    = -1;
  vars->hdrId_angle_incidence = -1;
  vars->maxTracesOut   = 0;
  vars->angle_incidence = NULL;


  cseis_geolib::csVector<std::string> valueList;
  
  param->getAll( "angle", &valueList );
  if( valueList.size() != 4 && valueList.size() != 3 ) {
    writer->error("User parameter 'angle': Wrong number of input parameters. Expected: 3 or 4, found: %d", valueList.size() );
  }
  vars->angle_first = atof( valueList.at(0).c_str() );
  vars->angle_last  = atof( valueList.at(1).c_str() );
  vars->angle_inc   = atof( valueList.at(2).c_str() );
  vars->angle_width = atof( valueList.at(3).c_str() );

  if( vars->angle_last < vars->angle_first ) {
    writer->error(" Last angle must be smaller or equal to first angle.");
  }
  if( vars->angle_first == vars->angle_last ) {
    vars->maxTracesOut = 1;
  }
  else if( vars->angle_inc <= 0.0 ) {
    writer->error(" Angle increment must be greater or equal to zero.");
  }
  else {
    vars->maxTracesOut = (int)((vars->angle_last - vars->angle_first) / vars->angle_inc) + 1;
  }
  vars->angle_first = DEG2RAD( vars->angle_first );
  vars->angle_last  = DEG2RAD( vars->angle_last );
  vars->angle_inc   = DEG2RAD( vars->angle_inc );
  vars->angle_width = DEG2RAD( vars->angle_width );

  vars->angle_incidence = NULL; //allocate_vector( gr->maxdtr );
//  vars->headerOut = allocate_vector( vars->maxTracesOut * gr->nth_size );
//  vars->tracesOut = allocate_matrix( vars->maxTracesOut, gr->numsmp );

  if( !hdef->headerExists("an_inci") ) {
    hdef->addStandardHeader( "an_inci" );
  }

  vars->hdrId_sou_x    = hdef->headerIndex( "sou_x" );
  vars->hdrId_sou_y    = hdef->headerIndex( "sou_y" );
  vars->hdrId_rec_x    = hdef->headerIndex( "rec_x" );
  vars->hdrId_rec_y    = hdef->headerIndex( "rec_y" );
  vars->hdrId_depth    = hdef->headerIndex( "sou_z" );
  vars->hdrId_offset   = hdef->headerIndex( "offset" );
  vars->hdrId_rec_elev = hdef->headerIndex( "rec_elev" );
  vars->hdrId_angle_incidence = hdef->headerIndex("an_inci" );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_off2angle_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;


  //---------------------------------------------------------------
  int nSamples = shdr->numSamples;
  int nTraces  = traceGather->numTraces();

  if( nTraces < 2 ) {
    writer->warning("%s: Gather with less than two traces.. This cannot be processed", edef->moduleName().c_str() );
    return;
  }

  float rcvz, srcz;
  double offset, daoffset, dz;
  float* an_inci = new float[nTraces];
//  int traceLeft, traceRight;
//  float weightLeft, weightRight;
//  float angleFirst, angleLast, angleCurrent, sin_angleLeft, sin_angleRight;
//  float direction, diff, min_angle, max_angle;
  // Read in header values
  // Compute angle of incidence
  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    csTraceHeader* trcHdr = traceGather->trace( itrc )->getTraceHeader();
    offset = trcHdr->doubleValue( vars->hdrId_offset );
    rcvz   = trcHdr->floatValue( vars->hdrId_rec_elev );
    srcz   = trcHdr->floatValue( vars->hdrId_depth );
    dz = rcvz - srcz;
    daoffset = sqrt( offset*offset + dz*dz );
    daoffset = std::max( daoffset, 0.0001 );

    an_inci[itrc] = acos(dz/daoffset);  // Always positive, assumes dz > 0
    if( offset < 0.0 ) {
      an_inci[itrc] *= -1.0;
    }
    if( edef->isDebug() ) writer->line("Trace: %3d, Offset: %f, Angle of incidence: %f\n", itrc, offset, RAD2DEG(an_inci[itrc]) );
  }

  // Check that input data is sorted correctly
  float direction, diff, min_angle, max_angle;
  direction = an_inci[nTraces-1] > an_inci[0] ? 1.0 : -1.0;
  for( int itrc = 1; itrc < nTraces; itrc++ ) {
    diff = an_inci[itrc] - an_inci[itrc-1];
    if( diff * direction < 0.0 ) {
      writer->error ("Input gather is not sorted correctly. Expected trace sorting is increasing or decreasing offsets.");
    }
  }
  if( direction > 0.0 ) {
    min_angle = an_inci[0];
    max_angle = an_inci[nTraces-1];
  }
  else {
    min_angle = an_inci[nTraces-1];
    max_angle = an_inci[0];
  }

  //***********************************************************************
  // Set available angle range for this gather
  //
  float angleFirst, angleLast;
  if( min_angle <= vars->angle_first ) {
    angleFirst = vars->angle_first;
  }
  else {
    angleFirst = (float)( (int)((min_angle - vars->angle_first)/vars->angle_inc) + 1 ) * vars->angle_inc + vars->angle_first;
  }
  if( max_angle >= vars->angle_last ) {
    angleLast = vars->angle_last;
  }
  else {
    angleLast = (float)(int)( (max_angle - vars->angle_first)/vars->angle_inc ) * vars->angle_inc + vars->angle_first;
  }

  int nTracesOut = (int)((angleLast - angleFirst) / vars->angle_inc + 0.5) + 1;
  if( nTracesOut > nTraces ) {
    traceGather->createTraces( nTraces, nTracesOut-nTraces, hdef, shdr->numSamples );
  }

  if( edef->isDebug() ) {
   writer->line(" in: %d, out: %d,  aFirst: %f, aLast: %f\n", nTraces, nTracesOut, RAD2DEG(angleFirst), RAD2DEG(angleLast));
   writer->line(" aFirst: %f, aLast: %f, nSamples: %d\n", RAD2DEG(vars->angle_first), RAD2DEG(vars->angle_last), nSamples);
  }

//***********************************************************************************
/*
  float angleCurrent, sin_angleRight, sin_angleLeft;
  int traceLeft, traceRight;
  float weightRight;
//  float weightLeft;

  if( direction > 0.0 ) {
    traceRight = 0;
    traceLeft = 0;
    for( int itrc = 0; itrc < nTracesOut; itrc++ ) {
      angleCurrent = angleFirst + (float)(itrc) * vars->angle_inc;
      while( traceRight < nTraces && an_inci[traceRight] < angleCurrent ) traceRight++;
      if( traceRight > 0 && traceRight < nTraces ) {
        traceLeft = traceRight-1;
        sin_angleLeft  = sin( an_inci[traceLeft] );
        sin_angleRight = sin( an_inci[traceRight] );
        if( sin_angleRight > sin_angleLeft ) {
          weightRight = (sin(angleCurrent) - sin_angleLeft) / (sin_angleRight - sin_angleLeft);
        }
        else {
          weightRight = 0.0;
        }
      }
      else {
        weightRight = 0.0;
      }
    //  weightLeft = 1.0 - weightRight;

      // Average trace headers:
      //headerStack( &rthdrs[traceLeft*gr->nth_size], weightLeft, &rthdrs[traceRight*gr->nth_size], weightRight, &vars->headerOut[itrc*gr->nth_size]);
      //  vars->headerOut[itrc*gr->nth_size+vars->hdrId_an_inci] = RAD2DEG( angleCurrent );

    //  for( int isamp = 0; isamp < nSamples; isamp++ ) {
    //    tracesOut[itrc][isamp] = weightLeft * tracesIn[traceLeft][isamp] + weightRight * tracesIn[traceRight][isamp];
    //  }
    }
  }
  else {  // Sort direction < 0.0
    traceRight = nTraces-1;
    traceLeft  = nTraces-1;
    for( int itrc = 0; itrc < nTracesOut; itrc++ ) {
      angleCurrent = angleFirst + (float)(itrc) * vars->angle_inc;
      while( traceRight >= 0 && an_inci[traceRight] < angleCurrent ) traceRight--;
      if( traceRight >= 0 && traceRight < nTraces-1 ) {
        traceLeft = traceRight+1;
        sin_angleLeft  = sin( an_inci[traceLeft] );
        sin_angleRight = sin( an_inci[traceRight] );
        if( sin_angleRight > sin_angleLeft ) {
          weightRight = (sin(angleCurrent) - sin_angleLeft) / (sin_angleRight - sin_angleLeft);
        }
        else {
          weightRight = 0.0;
        }
      }
      else {
        weightRight = 0.0;
      }
      //      weightLeft = 1.0 - weightRight;

      // Average trace headers:
      //headerStack( &rthdrs[traceLeft*gr->nth_size], weightLeft, &rthdrs[traceRight*gr->nth_size], weightRight, &vars->headerOut[itrc*gr->nth_size]);
      //  vars->headerOut[itrc*gr->nth_size+vars->hdrId_an_inci] = RAD2DEG( angleCurrent );

      //      for( int isamp = 0; isamp < nSamples; isamp++ ) {
      //  tracesOut[itrc][isamp] = weightLeft * tracesIn[traceLeft][isamp] + weightRight * tracesIn[traceRight][isamp];
      //  }
    }
  }

  // Copy output traces to input/output buffer
  for( int itrc = 0; itrc < nTracesOut; itrc++ ) {
 //   memcpy( &rthdrs[itrc*gr->nth_size], &vars->headerOut[itrc*gr->nth_size], gr->nth_size*4 );
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
//      tracesIn[itrc][isamp] = tracesOut[itrc][isamp];
    }
//    ithdrs[itrc*gr->nth_size+hdrId_traceno] = itrc + 1;
  }
*/

}
//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_off2angle_( csParamDef* pdef ) {
  pdef->setModule( "OFF2ANGLE", "Offset to angle transform (DEFUNCT MODULE)", "Transform seismic gather to regular trace interval of constant angle of incidence at seabed" );

  pdef->addParam( "angle", "Setup of output angle bands", NUM_VALUES_FIXED, "" );
  pdef->addValue( "0", VALTYPE_NUMBER, "First angle band [deg]" );
  pdef->addValue( "50", VALTYPE_NUMBER, "Last angle band [deg]" );
  pdef->addValue( "10", VALTYPE_NUMBER, "Angle band increment [deg]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Width of each angle band [deg]" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_off2angle_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_off2angle::VariableStruct* vars = reinterpret_cast<mod_off2angle::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_off2angle_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_off2angle::VariableStruct* vars = reinterpret_cast<mod_off2angle::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_off2angle_( csParamDef* pdef ) {
  params_mod_off2angle_( pdef );
}
extern "C" void _init_mod_off2angle_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_off2angle_( param, env, writer );
}
extern "C" bool _start_exec_mod_off2angle_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_off2angle_( env, writer );
}
extern "C" void _exec_mod_off2angle_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_off2angle_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_off2angle_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_off2angle_( env, writer );
}

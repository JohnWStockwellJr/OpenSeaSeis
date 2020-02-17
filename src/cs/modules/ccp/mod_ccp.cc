/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_methods.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: CCP
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_ccp {
  struct VariableStruct {
    float* times;  // [s]
    float* velocities;  // [m/s]
    float vpvs;
    float depth;
    float bin_scalar;
    int mode;
    int numTimes;
    int hdrId_offset;
    int hdrId_ccp;
    int hdrId_rcv;
    int hdrId_source;

    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_rec_x;
    int hdrId_rec_y;

    int hdrId_ccp_x;
    int hdrId_ccp_y;

    int* hdrId_ray2D_ccpDist;
    int* hdrId_ray2D_ccpDepth;
    int* hdrId_ray2D_rayTime;
    double* ccpDist;
    double* ccpDepth;
    double* rayTime;
    int numRay2D;
    int method;
  };
  static int const MODE_NONE   = 0;
  static int const MODE_NUMBER = 1;
  static int const MODE_COORD  = 2;
  static int const MODE_BOTH   = 3;
  static int const METHOD_CCP_ISO = 11;
  static int const METHOD_RAY2D   = 12;
}
using mod_ccp::VariableStruct;

//*************************************************************************************************
// Init phase//
//
//*************************************************************************************************
void init_mod_ccp_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->times        = NULL;
  vars->velocities   = NULL;
  vars->vpvs         = 0;
  vars->depth        = 0;
  vars->bin_scalar   = 2.0;
  vars->numTimes     = 0;
  vars->hdrId_offset = -1;
  vars->hdrId_ccp    = -1;
  vars->hdrId_rcv    = -1;
  vars->hdrId_source = -1;

  vars->mode          = mod_ccp::MODE_BOTH;
  vars->hdrId_rec_x    = -1;
  vars->hdrId_rec_y    = -1;
  vars->hdrId_sou_x    = -1;
  vars->hdrId_sou_y    = -1;
  vars->hdrId_ccp_x    = -1;
  vars->hdrId_ccp_y    = -1;

  vars->method       = mod_ccp::METHOD_CCP_ISO;
  vars->hdrId_ray2D_ccpDist = NULL;
  vars->hdrId_ray2D_ccpDepth = NULL;
  vars->hdrId_ray2D_rayTime    = NULL;
  vars->numRay2D            = 0;
  vars->ccpDist  = NULL;
  vars->ccpDepth = NULL;
  vars->rayTime  = NULL;

//---------------------------------------------
// Retrieve times and velocities
//
  csVector<std::string> valueList;

  vars->times = NULL;
  vars->velocities = NULL;

  //-------------------------------------------------------------
  //
  std::string text;

  if( param->exists("method") ) {
     param->getString("method",&text);
    if( !text.compare("ray2d") ) {
      vars->method = mod_ccp::METHOD_RAY2D;
    }
    else if( !text.compare("ccp_iso") ) {
      vars->method = mod_ccp::METHOD_CCP_ISO;
    }
    else {
      writer->error("Uknown option: %s", text.c_str());
    }
  }

  if( vars->method == mod_ccp::METHOD_CCP_ISO ) {
    param->getFloat( "depth", &vars->depth );
    param->getFloat( "vpvs", &vars->vpvs );
    
    if( param->exists("compute_xy") ) {
      param->getString("compute_xy",&text);
      if( !text.compare("yes") ) {
        vars->mode = mod_ccp::MODE_BOTH;
      }
      else if( !text.compare("no") ) {
        vars->mode = mod_ccp::MODE_NUMBER;
      }
      else {
        writer->error("Uknown option: %s", text.c_str());
      }
    }
    
    if( param->exists("compute_bin") ) {
      param->getString("compute_bin",&text);
      if( !text.compare("yes") ) {
        // Keep current setting
      }
      else if( !text.compare("no") ) {
        if( vars->mode == mod_ccp::MODE_BOTH ) {
          vars->mode = mod_ccp::MODE_COORD;
        }
        else {
          vars->mode = mod_ccp::MODE_NONE;
        }
      }
      else {
        writer->error("Uknown option: %s", text.c_str());
      }
    }
  }
  //----------------------------------------------------------------
  else {
    vars->mode = mod_ccp::MODE_BOTH;
    csVector<int> ccpxHdrIdList;
    csVector<int> ccpzHdrIdList;
    csVector<int> rayTimeHdrIdList;
    bool found = false;
    int counter = 1;
    do {
      char name_ccpx[7], name_ccpz[7], name_time[11];
      sprintf(name_ccpx,"ccpx%02d",counter);
      sprintf(name_ccpz,"ccpz%02d",counter);
      sprintf(name_time,"ray_time%02d",counter);
      name_ccpx[6] = '\0';
      name_ccpz[6] = '\0';
      name_time[10] = '\0';
      found = hdef->headerExists( name_ccpx ) && hdef->headerExists( name_ccpz ) && hdef->headerExists( name_time );
      if( found ) {
        ccpxHdrIdList.insertEnd( hdef->headerIndex(name_ccpx) );
        ccpzHdrIdList.insertEnd( hdef->headerIndex(name_ccpz) );
        rayTimeHdrIdList.insertEnd( hdef->headerIndex(name_time) );
        if( edef->isDebug() ) writer->line("Ray2D trace header names for CCP binning (#%d): %s %s %s",counter,name_ccpx,name_ccpz,name_time);
        counter += 1;
      }
    } while( found );
    vars->numRay2D = ccpxHdrIdList.size();
    //    fprintf(stderr,"Number of CCP trace headers: %d\n", vars->numRay2D);
    vars->hdrId_ray2D_ccpDist  = new int[vars->numRay2D];
    vars->hdrId_ray2D_ccpDepth = new int[vars->numRay2D];
    vars->hdrId_ray2D_rayTime  = new int[vars->numRay2D];
    vars->ccpDist  = new double[vars->numRay2D];
    vars->ccpDepth = new double[vars->numRay2D];
    vars->rayTime  = new double[vars->numRay2D];
    for( int i = 0; i < vars->numRay2D; i++ ) {
      vars->hdrId_ray2D_ccpDist[i]  = ccpxHdrIdList.at(i);
      vars->hdrId_ray2D_ccpDepth[i] = ccpzHdrIdList.at(i);
      vars->hdrId_ray2D_rayTime[i]  = rayTimeHdrIdList.at(i);
    }
  }

  if( vars->method == mod_ccp::METHOD_CCP_ISO && vars->mode == mod_ccp::MODE_NONE ) {
    writer->error("Computation of both the CCP bin number and CCP bin XY coordinates is turned off. At least one of these need to be comptued.");
  }

  //-------------------------------------------------------------
  // Velocities needed for time-depth conversion
  /*
  if( param->exists( "time" ) ) {
    param->getAll( "time", &valueList );
    if( valueList.size() < 1 ){
      writer->error("No times specified in user parameter 'time'!");
    }
  }
  vars->numTimes = valueList.size();

  if( vars->numTimes > 0 ) {
    vars->times = new float[vars->numTimes];
    for( int i = 0; i < vars->numTimes; i++ ) {
      vars->times[i] = atof( valueList.at(i).c_str() ) / 1000.0;  // Convert to seconds
    }
  }
  else {
    vars->numTimes = 1;
    vars->times = new float[vars->numTimes];
    vars->times[0] = 0;
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
//    vars->scalar[i] = atof( valueList.at(i).c_str() );
    if( !number.convertToNumber( valueList.at(i) ) ) {
      writer->error("Specified velocity is not a valid number: '%s'", valueList.at(i).c_str() );
    }
    vars->velocities[i] = number.floatValue();
    if( edef->isDebug() ) writer->line("Velocity #%d: '%s' --> %f", i, valueList.at(i).c_str(), vars->velocities[i] );
  }

*/

//  if( vars->mode == mod_ccp::MODE_BOTH || vars->mode == mod_ccp::MODE_NUMBER ) {
  if( (vars->mode & mod_ccp::MODE_NUMBER) != 0 ) {
    if( param->exists("bin_scalar") ) {
      param->getFloat("bin_scalar",&vars->bin_scalar);
    }
    if( !hdef->headerExists( "rcv" ) ) {
      writer->error("Trace header 'rcv' does not exist.");
    }
    if( !hdef->headerExists( "source" ) ) {
      writer->error("Trace header 'source' does not exist.");
    }
    if( !hdef->headerExists( "ccp" ) ) {
      hdef->addStandardHeader( "ccp" );
    }
    vars->hdrId_rcv    = hdef->headerIndex( "rcv" );
    vars->hdrId_source = hdef->headerIndex( "source" );
    vars->hdrId_ccp = hdef->headerIndex( "ccp" );
  }
//  if( vars->mode == mod_ccp::MODE_BOTH || vars->mode == mod_ccp::MODE_COORD ) {
  if( (vars->mode & mod_ccp::MODE_COORD) != 0 ) {
    vars->hdrId_sou_x = hdef->headerIndex( "sou_x" );
    vars->hdrId_sou_y = hdef->headerIndex( "sou_y" );
    vars->hdrId_rec_x = hdef->headerIndex( "rec_x" );
    vars->hdrId_rec_y = hdef->headerIndex( "rec_y" );
    if( !hdef->headerExists( "ccp_x" ) ) {
      hdef->addHeader( TYPE_DOUBLE, "ccp_x" );
    }
    if( !hdef->headerExists( "ccp_y" ) ) {
      hdef->addHeader( TYPE_DOUBLE, "ccp_y" );
    }
    vars->hdrId_ccp_x = hdef->headerIndex( "ccp_x" );
    vars->hdrId_ccp_y = hdef->headerIndex( "ccp_y" );
  }
  //-------------------------------------------------------------
  // Headers
  //
  if( !hdef->headerExists( "offset" ) ) {
    writer->error("Trace header 'offset' does not exist.");
  }
  vars->hdrId_offset = hdef->headerIndex( "offset" );

  // CCP scalar = 1: CCP numbering follows that of source and receiver stations
  // CCP scalar = 2: CCP numbering is twice that of source and receiver stations
  // ...
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_ccp_(
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


  csTrace* trace = traceGather->trace(0);
  csTraceHeader* trcHdr = trace->getTraceHeader();
  double sou_x  = trcHdr->doubleValue( vars->hdrId_sou_x );
  double sou_y  = trcHdr->doubleValue( vars->hdrId_sou_y );
  double rec_x  = trcHdr->doubleValue( vars->hdrId_rec_x );
  double rec_y  = trcHdr->doubleValue( vars->hdrId_rec_y );
  //  fprintf(stderr,"RS  %f %f   %f %f\n", rec_x, rec_y, sou_x, sou_y);

  if( vars->method == mod_ccp::METHOD_CCP_ISO ) {
    double offset  = trcHdr->doubleValue( vars->hdrId_offset );
    double offset_pside;
    double min_offset = 0.0001;
    double ratio = 1.0;

    double rcv     = (double)trcHdr->intValue( vars->hdrId_rcv );
    double source  = (double)trcHdr->intValue( vars->hdrId_source );

    double ccp_x = rec_x;
    double ccp_y = rec_y;

    if( fabs(offset) >= min_offset ) {
      ccp_offset_iso( vars->vpvs, offset, vars->depth, &offset_pside );
      ratio = (offset_pside/offset);
      if( edef->isDebug() ) writer->line("CCP ratio: %f, offset: %f, offset_pside: %f", ratio, offset, offset_pside);
    }
    if( (vars->mode & mod_ccp::MODE_COORD) != 0 ) {
      if( fabs(offset) >= min_offset ) {
        ccp_x = sou_x + ratio * (rec_x - sou_x);
        ccp_y = sou_y + ratio * (rec_y - sou_y);
      }
      //    else {
      //   ccp_x = rec_x;
      //   ccp_y = rec_y;
      // }
      trcHdr->setDoubleValue( vars->hdrId_ccp_x, ccp_x );
      trcHdr->setDoubleValue( vars->hdrId_ccp_y, ccp_y );
    }
    if( (vars->mode & mod_ccp::MODE_NUMBER) != 0 ) {
      int ccp;
      if( fabs(offset) >= min_offset ) {
        ccp = (int)round( vars->bin_scalar * ( source + ratio*( rcv - source ) ));
      }
      else {
        ccp = (int)round( vars->bin_scalar*source );
      }
      trcHdr->setIntValue( vars->hdrId_ccp, ccp );
    }
  }
  // RAY2D:
  else {
    for( int iseg = 0; iseg < vars->numRay2D; iseg++ ) {
      vars->rayTime[iseg] = trcHdr->floatValue( vars->hdrId_ray2D_rayTime[iseg] );
      vars->ccpDist[iseg] = fabs( trcHdr->floatValue( vars->hdrId_ray2D_ccpDist[iseg] ) );
      vars->ccpDepth[iseg] = trcHdr->floatValue( vars->hdrId_ray2D_ccpDepth[iseg] );
      //      fprintf(stderr,"#%d %f %f %f\n", iseg, vars->rayTime[iseg], vars->ccpDist[iseg], vars->ccpDepth[iseg]);
    }

    //    csVector
    double bin_size = 12.5; // [m]
    double bin_size_half = 0.5*bin_size; // [m]
    //    double posBorderBin1 = rec_x - bin_size_half;
    //    double currentPos = posBorderBin1;  // (x) or really xy
    csVector<double> timeCutList;  // List of exact times where cut shall be made (+/- halftaper)
    csVector<double> distCutList;   // List of CCP X locations for each newly cut segment (first position = rcvpos, so no need to store in list)
    int bin = 0;
    float maxTime = (float)shdr->numSamples * shdr->sampleInt;
    for( int iseg = vars->numRay2D-1; iseg >= 0; iseg-- ) {
      //      double topTime  = 0.0;
      double topDepth = 0.0;
      double topDist  = 0.0;
      if( iseg < vars->numRay2D-1 ) {
        //        topTime  = vars->rayTime[iseg+1];
        topDepth = vars->ccpDepth[iseg+1];
        topDist  = vars->ccpDist[iseg+1];
      }
      double dx = vars->ccpDist[iseg] - topDist;
      if( dx < 0.1 ) dx = 0.1;
      double dz = vars->ccpDepth[iseg] - topDepth;
      //      double dt = vars->rayTime[iseg] - topTime;
      int bin1 = (int)( topDist / bin_size ) + 1;
      int bin2 = (int)( vars->ccpDist[iseg] / bin_size ) + 1;
      int numCuts = bin2 - bin1;
      bin = bin1;
      //      fprintf(stdout,"SEG  %d %d %d   %f    %f   topDepth:  %f   %f %f\n", bin, bin1, bin2, dz, dx, topDist, topDepth, vars->ccpDist[iseg] );
      for( int icut = 0; icut < numCuts; icut++ ) {
        double distCut = ((double)bin * bin_size) - bin_size_half;
        double time = (dz/dx) * (distCut-topDist) + topDepth;
        //        fprintf(stdout,"%d  %d    %f  %f    %f  %f    time/topDepth:  %f  %f\n", bin, icut, dz, dx, distCut, topDist, time, topDepth);
        timeCutList.insertEnd(time);
        distCutList.insertEnd(distCut);
        bin += 1;
        if( timeCutList.last() > maxTime ) break;
      }
      if( timeCutList.size() > 0 && timeCutList.last() > maxTime ) break;
    }

    if( timeCutList.size() == 0 || timeCutList.last() < maxTime ) {
      timeCutList.insertEnd( maxTime );
      if( distCutList.size() != 0 ) {
        //        float cutTmp = distCutList.at( distCutList.size()-1 );
        float cutTmp = distCutList.last();
        distCutList.insertEnd( cutTmp );
      }
      else {
        distCutList.insertEnd(0);
      }
    }

    int numCuts = timeCutList.size();
    /*
    fprintf(stderr,"Number of cuts: %d\n", numCuts);
    for( int i = 0; i < numCuts; i++ ) {
      double diff = 0.0;
      if( i > 0 ) diff = timeCutList.at(i) - timeCutList.at(i-1);
            fprintf(stdout,"%d %f %f %f  %f\n", i, distCutList.at(i), timeCutList.at(i),diff,trcHdr->floatValue(hdef->headerIndex("offset")) );
    }
    */
    if( numCuts > 0 ) {
      int numTracesIn = traceGather->numTraces();
      if( numCuts > 1 ) traceGather->createTraces( numTracesIn, numCuts-1, hdef, shdr->numSamples );  // numCuts-1: Existing trace serves as first output trace
      double rs_azim_rad = atan2( sou_x-rec_x, sou_y-rec_y );

      float* samplesIn = traceGather->trace(0)->getTraceSamples();
      int currentTopSample = 0;
      for( int icut = 0; icut < numCuts; icut++ ) {
        float* samplesNew = traceGather->trace(icut)->getTraceSamples();
        double dist = distCutList.at(icut) - bin_size_half;
        if( dist < 0.0 ) dist = 0.0;  // Should be corrected further up. Avoid negative 'cut' distance
        double ccp_x = rec_x + sin(rs_azim_rad) * dist;
        double ccp_y = rec_y + cos(rs_azim_rad) * dist;
        //        fprintf(stdout,"BINXY  %f %f   %f %f   %f %f  %f  --  %d\n", rec_x, rec_y, ccp_x, ccp_y, sou_x, sou_y, dist,
        //        trcHdr->intValue(hdef->headerIndex("source")) );
        int sampleBottom = (int)round( timeCutList.at(icut) / shdr->sampleInt );
        if( icut == numCuts-1 || sampleBottom > shdr->numSamples ) sampleBottom = shdr->numSamples;
        int numSamples2Copy = sampleBottom - currentTopSample;
        csTraceHeader* trcHdrNew = traceGather->trace(icut)->getTraceHeader();
        if( icut > 0 ) {
          trcHdrNew->copyFrom( trcHdr );
          //          fprintf(stdout,"SAMPLES: %d %d   %d\n", currentTopSample, sampleBottom, numSamples2Copy);
          if( currentTopSample < 0 ) {
            //            fprintf(stderr,"SAMPLES: %d %d   %d\n", currentTopSample, sampleBottom, numSamples2Copy);
            writer->error("CCP  ERROR: Sample index < 0: %d", currentTopSample);
          }
          else if( currentTopSample >= shdr->numSamples ) {
            //            fprintf(stderr,"SAMPLES: %d %d   %d\n", currentTopSample, sampleBottom, numSamples2Copy);
            writer->error("CCP  ERROR: Sample index > numSamples: %d", currentTopSample);
          }
          else if( sampleBottom > shdr->numSamples ) {
            //            fprintf(stderr,"SAMPLES: %d %d   %d\n", currentTopSample, sampleBottom, numSamples2Copy);
            writer->error("CCP  ERROR: Bottom Sample index > numSamples: %d", sampleBottom );
          }

          if( currentTopSample+numSamples2Copy > shdr->numSamples ) {
            //            fprintf(stderr,"SAMPLES: %d %d   %d\n", currentTopSample, sampleBottom, numSamples2Copy);
            writer->error("CCP  ERROR: Sample index++ > numSamples: %d", currentTopSample+numSamples2Copy,numSamples2Copy);
          }

          for( int isamp = 0; isamp < currentTopSample; isamp++ ) {
            samplesNew[isamp] = 0.0;
          }
          float* dummy = new float[numSamples2Copy];
          memcpy( dummy, &samplesIn[currentTopSample], numSamples2Copy*sizeof(float) );
          memcpy( &samplesNew[currentTopSample], dummy, numSamples2Copy*sizeof(float) );
          memcpy( &samplesNew[currentTopSample], &samplesIn[currentTopSample], numSamples2Copy*sizeof(float) );
          delete [] dummy;
          for( int isamp = sampleBottom; isamp < shdr->numSamples; isamp++ ) {
            samplesNew[isamp] = 0.0;
          }
        }
        trcHdrNew->setDoubleValue( vars->hdrId_ccp_x, ccp_x );
        trcHdrNew->setDoubleValue( vars->hdrId_ccp_y, ccp_y );
        currentTopSample = sampleBottom;
      }
      // Set bottom samples in first trace to zero
      int sampleBottom = (int)round( timeCutList.at(0) / shdr->sampleInt );
      if( sampleBottom < shdr->numSamples ) {
        for( int isamp = std::max(0,sampleBottom); isamp < shdr->numSamples; isamp++ ) {
          samplesIn[isamp] = 0.0;
        }
      }
    } // END if numCuts != 0
    else {
      trcHdr->setDoubleValue( vars->hdrId_ccp_x, rec_x );
      trcHdr->setDoubleValue( vars->hdrId_ccp_y, rec_y );
    }
    /*
- User input: Minimum output trace interval in meters: bin_size
- Check rec/sou XY's. Divide source-receiver path into N equal 'bins' of size ~bin_size m
- Create CCP depth function, interpolate linearly
   time  ccpx_dist  (= offset from source to CCP)
- Start at time = 0 & ccpx_dist = receiver location, and move down
- for each function segment:
  currentTop = Keep segment i top position 
  if top position and bottom position are within distance <= bin_size/2, continue with next segment
  if( bottom position of segment is beyond >= bin_size/2

bin_size = 12.5; // [m]
bin_size_half = 0.5*bin_size; // [m]
double posBorderBin1 = rcvPos + bin_size_half;
double currentPos = posBorderBin1;  // (x) or really xy
csVector<double> timeCutList();  // List of exact times where cut shall be made (+/- halftaper)
csVector<double> distCutList();   // List of CCP X locations for each newly cut segment (first position = rcvpos, so no need to store in list)
int bin = 1;
for( int iseg = 0; iseg < numSegments; iseg++ ) {
  double dx = segment[iseg].bottomPos - ((double)bin * bin_size + posBorderBin1);
  int numCuts = (int)( dx / bin_size ) + 1;
  for( int icut = 0; icut < numCuts; icut++ ) {
    double distCut = ((double)bin * bin_size + posBorderBin1);
    double dt = segment[iseg].bottomTime - segment[iseg].topTime;
    double dx = segment[iseg].bottomPos - segment[iseg].topPos;
    double time = (dt/dx) * (distCut-segment[iseg].topPos) + segment[iseg].topTime;
    timeCutList.insertEnd(time);
    distCutList.insertEnd(distCut-bin_size_half);
    bin += 1;
  }
}
    */

  }


  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_ccp_( csParamDef* pdef ) {
  pdef->setModule( "CCP", "Perform ccp binning" );
  pdef->addDoc("Current implementation is for isotropic assumption only");

  pdef->addParam( "method", "Method for CCP binning", NUM_VALUES_FIXED );
  pdef->addValue( "ccp_iso", VALTYPE_OPTION );
  pdef->addOption( "ccp_iso", "Use isotropic CCP equation, assuming source & receiver are at same datum" );
  pdef->addOption( "ray2d", "Perform CCP binning using output from module 'RAY2D'", "No further user parameter input is required" );

  pdef->addParam( "compute_xy", "Compute bin XY location?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Compute CCP bin XY location" );

  pdef->addParam( "compute_bin", "Compute CCP bin number?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Compute CCP bin number" );
  pdef->addOption( "no", "Do not compute CCP bin number" );

  pdef->addParam( "bin_scalar", "Numbering of CCP: Apply scalar to 'normal' CCP number", NUM_VALUES_FIXED );
  pdef->addValue( "2", VALTYPE_NUMBER, "CCP number is computed as follows: ccp = bin_scalar * ratio * (source + rcv)",
    "...the ratio is the ratio between the offset between the source and CCP location and the source-receiver offset" );

  pdef->addParam( "depth", "Target depth at which CCP binning shall be performed", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Target depth [m]" );

  pdef->addParam( "vpvs", "Effective Vp/Vs ratio", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Effective Vp/Vs ratio at target level" );

  pdef->addParam( "time", "List of time value [ms]", NUM_VALUES_VARIABLE, "Time knee points at which specified velocities apply" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of time values [ms]..." );

  pdef->addParam( "velocity", "List of average vertical velocities [m/s]", NUM_VALUES_VARIABLE,
      "Velocity values at specified time knee points. In between time knee points, velocities are linearly interpolated." );
  pdef->addValue( "", VALTYPE_NUMBER, "List of average vertical velocities [m/s]..." );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_ccp_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_ccp::VariableStruct* vars = reinterpret_cast<mod_ccp::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_ccp_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_ccp::VariableStruct* vars = reinterpret_cast<mod_ccp::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->times != NULL ) {
    delete [] vars->times;
    vars->times = NULL;
  }
  if( vars->velocities != NULL ) {
    delete [] vars->velocities;
    vars->velocities = NULL;
  }
  if( vars->hdrId_ray2D_ccpDist != NULL ) {
    delete [] vars->hdrId_ray2D_ccpDist;
    vars->hdrId_ray2D_ccpDist = NULL;
  }
  if( vars->hdrId_ray2D_ccpDepth != NULL ) {
    delete [] vars->hdrId_ray2D_ccpDepth;
    vars->hdrId_ray2D_ccpDepth = NULL;
  }
  if( vars->hdrId_ray2D_rayTime != NULL ) {
    delete [] vars->hdrId_ray2D_rayTime;
    vars->hdrId_ray2D_rayTime    = NULL;
  }
  if( vars->rayTime != NULL ) {
    delete [] vars->rayTime;
    vars->rayTime    = NULL;
  }
  if( vars->ccpDist != NULL ) {
    delete [] vars->ccpDist;
    vars->ccpDist    = NULL;
  }
  if( vars->ccpDepth != NULL ) {
    delete [] vars->ccpDepth;
    vars->ccpDepth    = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_ccp_( csParamDef* pdef ) {
  params_mod_ccp_( pdef );
}
extern "C" void _init_mod_ccp_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_ccp_( param, env, writer );
}
extern "C" bool _start_exec_mod_ccp_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_ccp_( env, writer );
}
extern "C" void _exec_mod_ccp_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_ccp_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_ccp_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_ccp_( env, writer );
}

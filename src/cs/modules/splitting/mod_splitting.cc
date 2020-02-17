/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"
#include "geolib_math.h"
#include "geolib_methods.h"
#include "csTimeStretch.h"
#include "cseis_curveFitting.h"
#include "csFileUtils.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPLITTING
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_splitting {
  struct VariableStruct {
    int outputHdrOption;
    int outputOption;
    
    int nAngles;
    float angleInc_rad;
    float** s1;  // Trace sample buffer holding S1 amplitudes
    float** s2;  // Trace sample buffer holding S2 amplitudes
    FILE* fout;

    float angleWidthOmit;  // Width of angle band to omit from analysis
    bool doLayerStripping;
    int ensembleCounter;

    int normMethod;
    int s2Scaling;
    float azim_rad;
    int corrMaxlag_samples;
    float minOffset;
    float maxOffset;
    int numWindows;
    float* winStart;
    float* winEnd;
    float* corr_lags;
    int rotationMode;
    bool stretchTopHalf;

    int hdrId_offset;
    int hdrId_rec_x;
    int hdrId_rec_y;
    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_sensor;
    int hdrId_s1az;
    int hdrId_window;
    int hdrId_identifier;

    std::string filename_info;
    bool isFirstCall;
 };
  static int const OUTPUT_S1S2_STACKS     = 4;
  static int const OUTPUT_CORRECTED_DATA  = 8;
  static int const OUTPUT_LAST_S1S2_STACKS = 12;

  static int const OUTPUT_HDR_FIRST   = 1;
  static int const OUTPUT_HDR_LAST    = 2;
  static int const OUTPUT_HDR_AVERAGE = 0;

  static int const NORM_METHOD_LS = 11;
  static int const NORM_METHOD_NTRACES = 12;

  static int const S2_SCALING_ISOTROPIC   = 21;
  static int const S2_SCALING_ANISOTROPIC = 22;

  static const int MODE_APPLY  = 11;
  static const int MODE_REMOVE = 12;
}

void rotate_xy( float* xTrace,
                float* yTrace,
                float azim,
                int nSamples );

void compute_s1s2_splitting_stacks( float** samples, float const* sr_azim, int numSamples, int numPairs, float azim,
    bool isDebug, FILE* logFile,
    float angleWidthOmit, int normMethod, int s2Scaling, int nAngles, float angleInc,
    float** s1, float** s2 );

  void compute_twosided_correlation2( float const* samplesLeft, float const* samplesRight,
                                   int nSampIn, float* corr, int maxlag_in_num_samples );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_splitting_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  mod_splitting::VariableStruct* vars = new mod_splitting::VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->outputHdrOption = mod_splitting::OUTPUT_HDR_AVERAGE;
  vars->outputOption    = mod_splitting::OUTPUT_S1S2_STACKS;
  vars->nAngles    = 0;
  vars->angleInc_rad = 0.0;
  vars->s1         = NULL;
  vars->s2         = NULL;
  vars->normMethod = mod_splitting::NORM_METHOD_LS;
  vars->azim_rad     = 0;
  vars->hdrId_rec_x  = -1;
  vars->hdrId_rec_y  = -1;
  vars->hdrId_sou_x  = -1;
  vars->hdrId_sou_y  = -1;
  vars->hdrId_sensor = -1;
  vars->hdrId_s1az   = -1;
  vars->hdrId_window = -1;
  vars->hdrId_identifier = -1;
  vars->angleWidthOmit = DEG2RAD( 5.0f );
  vars->s2Scaling    = mod_splitting::S2_SCALING_ANISOTROPIC;
  vars->stretchTopHalf = false;

  vars->fout = NULL;

  vars->minOffset  = 0.0;
  vars->maxOffset  = 0.0;
  vars->corrMaxlag_samples = 0;
  vars->numWindows = 0;
  vars->winStart   = NULL;
  vars->winEnd     = NULL;
  vars->doLayerStripping = false;
  vars->corr_lags = NULL;
  vars->rotationMode = -1;
  vars->ensembleCounter = 0;

  vars->filename_info = "";
  vars->isFirstCall = true;

  std::string text;

  if( param->exists("layer_strip") ) {
    param->getString("layer_strip", &text);
    if( !text.compare("yes") ) {
      vars->doLayerStripping = true;
    }
    else if( !text.compare("no") ) {
      vars->doLayerStripping = false;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
  }
  param->getFloat("azim", &vars->azim_rad);
  vars->azim_rad = DEG2RAD( vars->azim_rad );
  param->getInt("num_angles", &vars->nAngles);
  float angleInc_deg = 180.0f / (float)vars->nAngles;
  vars->angleInc_rad = DEG2RAD( angleInc_deg );
  writer->line("Number of angles = %d, Angle increment = %.2fdeg", vars->nAngles, angleInc_deg);

  if( vars->doLayerStripping ) {
    vars->corr_lags = new float[vars->nAngles];
    float maxLag_ms;
    param->getFloat("corr_maxlag", &maxLag_ms);
    vars->corrMaxlag_samples = (int)(maxLag_ms/shdr->sampleInt);
    writer->line("Perform layer stripping. Maximum lag time for S1/S2 cross-correlation = %.2fms (= Minimum expected wavelet length)",
      maxLag_ms);

    param->getString( "orient_xy", &text, 0 );
    if( !text.compare("right") ) {
      vars->rotationMode = mod_splitting::MODE_REMOVE;
    }
    else if( !text.compare("left") ) {
      vars->rotationMode = mod_splitting::MODE_APPLY;
    }
  }

  if( param->exists("stretch_top_half") ) {
    param->getString("stretch_top_half", &text);
    if( !text.compare("yes") ) {
      vars->stretchTopHalf = true;
    }
    else if( !text.compare("no") ) {
      vars->stretchTopHalf = false;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
  }

  if( param->exists("output_hdr") ) {
    param->getString("output_hdr", &text);
    if( !text.compare("first") ) {
      vars->outputHdrOption = mod_splitting::OUTPUT_HDR_FIRST;
    }
    else if( !text.compare("last") ) {
      vars->outputHdrOption = mod_splitting::OUTPUT_HDR_LAST;
    }
    else if( !text.compare("average") ) {
      writer->error("Option 'average' is not supported yet...");
      vars->outputHdrOption = mod_splitting::OUTPUT_HDR_AVERAGE;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
  }

  if( param->exists("output") ) {
    param->getString("output", &text);
    if( !text.compare("s1s2_stacks") ) {
      vars->outputOption = mod_splitting::OUTPUT_S1S2_STACKS;
    }
    else if( !text.compare("corrected_data") ) {
      vars->outputOption = mod_splitting::OUTPUT_CORRECTED_DATA;
      if( !vars->doLayerStripping ) {
        writer->error("Output option 'corrected_data' is only applicable when performing layer stripping");
      }
    }
    else if( !text.compare("last_s1s2_stacks") ) {
      vars->outputOption = mod_splitting::OUTPUT_LAST_S1S2_STACKS;
      if( !vars->doLayerStripping ) {
        vars->outputOption = mod_splitting::OUTPUT_S1S2_STACKS;
      }
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
  }

  if( param->exists("norm_method") ) {
    param->getString("norm_method", &text);
    if( !text.compare("ls") ) {
      vars->normMethod = mod_splitting::NORM_METHOD_LS;
    }
    else if( !text.compare("ntraces") ) {
      vars->normMethod = mod_splitting::NORM_METHOD_NTRACES;
    }
    else if( !text.compare("isotropic") ) {
      vars->normMethod = mod_splitting::NORM_METHOD_NTRACES;
      vars->s2Scaling  = mod_splitting::S2_SCALING_ISOTROPIC;
    }
    else {
      writer->line("Unknown option: '%s'", text.c_str());
    }
  }

  if( param->exists("angle_omit") ) {
    param->getFloat("angle_omit", &vars->angleWidthOmit);
    vars->angleWidthOmit = DEG2RAD( vars->angleWidthOmit );
  }

  if( param->exists("offset_range") ) {
    param->getFloat("offset_range", &vars->minOffset,0);
    if( param->getNumValues("offset_range") > 1 ) {
      param->getFloat("offset_range", &vars->maxOffset,1);
    }
  }

  if( param->exists("win_start") ) {
    vars->numWindows = param->getNumValues("win_start");
    if( vars->numWindows != param->getNumValues("win_end") ) {
      writer->error("Unequal window start and end times: %d != %d", vars->numWindows, param->getNumValues("win_end") );
    }
    vars->winStart = new float[vars->numWindows];
    vars->winEnd   = new float[vars->numWindows];
    for( int iwin = 0; iwin < vars->numWindows; iwin += 1 ) {
      param->getFloat("win_start", &vars->winStart[iwin], iwin);
      param->getFloat("win_end", &vars->winEnd[iwin], iwin);
    }
  }
  if( param->exists("write_info") ) {
    param->getString("write_info", &text, 0);
    vars->hdrId_identifier = hdef->headerIndex(text);
    param->getString("write_info", &vars->filename_info, 1);
    if( !csFileUtils::createDoNotOverwrite( vars->filename_info ) ) {
      writer->error("Error opening output file %s\n", vars->filename_info.c_str() );
    }
  }

  vars->s1 = new float*[vars->nAngles];
  vars->s2 = new float*[vars->nAngles];
  for( int iangle = 0; iangle < vars->nAngles; iangle++ ) {
    vars->s1[iangle] = new float[shdr->numSamples];
    vars->s2[iangle] = new float[shdr->numSamples];
  }

  vars->hdrId_rec_x  = hdef->headerIndex("rec_x");
  vars->hdrId_rec_y  = hdef->headerIndex("rec_y");
  vars->hdrId_sou_x  = hdef->headerIndex("sou_x");
  vars->hdrId_sou_y  = hdef->headerIndex("sou_y");
  vars->hdrId_sensor = hdef->headerIndex("sensor");
  vars->hdrId_offset = hdef->headerIndex("offset");
  if( !hdef->headerExists("s1az") ) {
    hdef->addHeader( TYPE_FLOAT, "s1az", "S1 angle [deg]" );
  }
  vars->hdrId_s1az  = hdef->headerIndex("s1az");
  if( !hdef->headerExists("window") ) {
    hdef->addHeader( TYPE_INT, "window", "Analysis window" );
  }
  vars->hdrId_window = hdef->headerIndex("window");

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_splitting_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  mod_splitting::VariableStruct* vars = reinterpret_cast<mod_splitting::VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;


  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    if( vars->filename_info.size() > 0 ) {
      vars->fout = fopen(vars->filename_info.c_str(),"w");
      if( vars->fout == NULL ) {
        writer->error("Error when opening ASCII output file '%s'", vars->filename_info.c_str());
      }
    }
  }

  //---------------------------------------------------------------------
  // Ensemble mode
  // Current assumptions:
  //   - Input data is one consistent ensemble that can be stacked together (NMO etc has been applied)
  //   - Data is sorted into adjacent XY traces
  //
  
  int numTotalTraces  = traceGather->numTraces();
  if( numTotalTraces % 2 != 0 ) {
    writer->error("SPLITTING: Wrong input data sorting. Expected pairs of XY traces. Found uneven number (=%d) of input traces", numTotalTraces);
  }
  int numSamples = shdr->numSamples;
  int numTotalPairs   = numTotalTraces/2;

  csVector<int> traceIndexList;

  for( int ip = 0; ip < numTotalPairs; ip++ ) {
    int traceIndex = ip*2; // 'X' trace
    csTraceHeader* trcHdr = traceGather->trace(traceIndex)->getTraceHeader();
    float offset = fabs(trcHdr->floatValue(vars->hdrId_offset));
    if( offset >= vars->minOffset && (vars->maxOffset == 0.0 || offset <= vars->maxOffset) ) {
      traceIndexList.insertEnd(traceIndex);
    }
  }

  int numUsedPairs = traceIndexList.size();
  int numUsedTraces = numUsedPairs*2;
  float* sr_azim = new float[numUsedPairs];
  float** samples = new float*[numUsedTraces];
  for( int ipair = 0; ipair < numUsedPairs; ipair++ ) {
    int traceIndexOut = ipair*2;
    int trueTraceIndex = traceIndexList.at(ipair);
    csTraceHeader* trcHdr = traceGather->trace(trueTraceIndex)->getTraceHeader();
    double rec_x = trcHdr->doubleValue(vars->hdrId_rec_x);
    double rec_y = trcHdr->doubleValue(vars->hdrId_rec_y);
    double sou_x = trcHdr->doubleValue(vars->hdrId_sou_x);
    double sou_y = trcHdr->doubleValue(vars->hdrId_sou_y);
    int sensor1  = trcHdr->intValue(vars->hdrId_sensor);
    int sensor2  = traceGather->trace(trueTraceIndex+1)->getTraceHeader()->intValue(vars->hdrId_sensor);
    if( sensor1 != 3 || sensor2 != 4 ) {
      writer->error("SPLITTING: Wrong input sensor code or sorting found: Trace #%d: %d, Trace #%d: %d. Expected values: 3 and 4 for X and Y",
        trueTraceIndex+1, trueTraceIndex+2, sensor1, sensor2);
    }
    double dx = rec_x - sou_x;
    double dy = rec_y - sou_y;
    sr_azim[ipair] = (float)fmod( atan2(dx,dy) + 2.0*M_PI, 2.0*M_PI );
    if( sensor1 == 3 ) { // 'X' trace
      if( sensor2 != 4 ) {
        writer->error("SPLITTING: Wrong input sensor code found: %d. Expected 3 and 4 for X and Y", sensor2);
      }
      samples[traceIndexOut]   = traceGather->trace(trueTraceIndex)->getTraceSamples();
      samples[traceIndexOut+1] = traceGather->trace(trueTraceIndex+1)->getTraceSamples();
    }
    else if( sensor1 == 4 ) { // 'Y' trace
      if( sensor2 != 3 ) {
        writer->error("SPLITTING: Wrong input sensor code found: %d. Expected 3 and 4 for X and Y", sensor2);
      }
      samples[traceIndexOut]   = traceGather->trace(trueTraceIndex+1)->getTraceSamples();
      samples[traceIndexOut+1] = traceGather->trace(trueTraceIndex)->getTraceSamples();
    }
    else {
      writer->error("SPLITTING: Wrong input sensor code found: %d. Expected 3 and 4 for X and Y", sensor1);
    }
  }
  // The following assumes a 2D right-handed XY coordinate system
  // Problem is that for some OBC, XY are defined anti-symmetric, so Y polarity has to be flipped

  if( edef->isDebug() ) {
    for( int iangle = 0; iangle < vars->nAngles; iangle++ ) {
      float s1az = (float)(iangle)*vars->angleInc_rad;
      float angle_xy_to_s1s2 = s1az - vars->azim_rad;
      writer->line("-------------------------------\nTESTING S1 ANGLE   %9.2f  %9.2f\n", RAD2DEG(s1az), RAD2DEG(angle_xy_to_s1s2) );
      for( int ip = 0; ip < numUsedPairs; ip++ ) {
        float angle_s1s2_to_srazim = sr_azim[ip] - s1az;
        writer->line("Trace XY pair #%3d, sr_azim: %9.2f, rotation angle: %9.2f  (s1az: %9.2f)", ip, RAD2DEG(sr_azim[ip]), RAD2DEG(angle_s1s2_to_srazim), RAD2DEG(s1az) );
      }
    }
  }

  //******************************************************************************************
  // Perform splitting analysis

  // Generate S1/S2 stacks. To be done regardless of layer stripping
  compute_s1s2_splitting_stacks( samples, sr_azim, numSamples, numUsedPairs, vars->azim_rad,
    edef->isDebug(), writer->getFile(),
    vars->angleWidthOmit, vars->normMethod, vars->s2Scaling, vars->nAngles/2, vars->angleInc_rad,
    vars->s1, vars->s2 );
  
   // Trick: Only a quarter of a circle is actually needed to uniquely sample the directions (see nAngles/2 in function call above)
   // --> Copy redundant traces to full 180deg half circle for nicer output display and easier handling regarding cross-correlation etc
   for( int iangle = 0; iangle < vars->nAngles/2; iangle++ ) {
     int indexTo = iangle + vars->nAngles/2;
     memcpy( &vars->s1[indexTo][0], &vars->s2[iangle][0], sizeof(float)*shdr->numSamples );
     memcpy( &vars->s2[indexTo][0], &vars->s1[iangle][0], sizeof(float)*shdr->numSamples );
   }
   
  //******************************************************************************************
  // If S1S2 stacks are requested to be written to output:
  // - Create new traces to hold S1S2 stacks and copy trace headers
  // -
  // Write S1S2 stacks to output if requested
  //
  int startIndexOutputTrace = numTotalTraces;

  if( vars->outputOption == mod_splitting::OUTPUT_S1S2_STACKS || vars->outputOption == mod_splitting::OUTPUT_LAST_S1S2_STACKS ) {
    // a) Create new traces to hold S1S2 stacks
    int numS1S2Stacks = (vars->outputOption == mod_splitting::OUTPUT_S1S2_STACKS) ? (vars->numWindows+1) : 1;
    traceGather->createTraces( numTotalTraces, numS1S2Stacks*2*vars->nAngles, hdef, shdr->numSamples );
    // b) Set trace headers for output S1S2 stack traces
    csTraceHeader* trcHdrRef = (vars->outputHdrOption == mod_splitting::OUTPUT_HDR_FIRST) ?
        traceGather->trace(0)->getTraceHeader():traceGather->trace(numTotalTraces-1)->getTraceHeader();
    for( int istack = 0; istack < numS1S2Stacks; istack++ ) {
      for( int iangle = 0; iangle < vars->nAngles; iangle++ ) {
        float s1az = RAD2DEG( (float)(iangle)*vars->angleInc_rad );
        int traceIndex = iangle*2 + istack*2*vars->nAngles + startIndexOutputTrace;
        csTraceHeader* trcHdr = traceGather->trace(traceIndex)->getTraceHeader();
        trcHdr->copyFrom(trcHdrRef);
        trcHdr->setIntValue( vars->hdrId_sensor, 3 );
        trcHdr->setFloatValue( vars->hdrId_s1az, s1az );
        trcHdr->setIntValue( vars->hdrId_window, istack );  // Window '0' = raw input data, no layer stripping
        trcHdr = traceGather->trace(traceIndex+1)->getTraceHeader();
        trcHdr->copyFrom(trcHdrRef);
        trcHdr->setIntValue( vars->hdrId_sensor, 4 );
        trcHdr->setFloatValue( vars->hdrId_s1az, s1az );
        trcHdr->setIntValue( vars->hdrId_window, istack );
      }
    }
    // c) Copy first (raw) S1S2 stack to output trace if required
    if( vars->outputOption == mod_splitting::OUTPUT_S1S2_STACKS ) {
      for( int iangle = 0; iangle < vars->nAngles; iangle++ ) {
        int traceIndex = iangle*2 + startIndexOutputTrace;
        memcpy( traceGather->trace(traceIndex)->getTraceSamples(), vars->s1[iangle], shdr->numSamples*sizeof(float) );
        memcpy( traceGather->trace(traceIndex+1)->getTraceSamples(), vars->s2[iangle], shdr->numSamples*sizeof(float) );
      }
      startIndexOutputTrace += 2*vars->nAngles;
    }
  }


  // Perform layer stripping, if requested
  if( vars->doLayerStripping ) {
    int nSampCorr  = 2 * vars->corrMaxlag_samples + 1;
    float* corrBuffer  = new float[nSampCorr];
    float* traceBuffer = new float[shdr->numSamples];
    int winIndex = 0;

    // Loop over all layer-stripping windows
    do {
      int startIndex = (int)round( vars->winStart[winIndex] / shdr->sampleInt );
      int nSampIn    = (int)round( vars->winEnd[winIndex] / shdr->sampleInt ) - startIndex + 1;
      // Cross-correlate S1/S2 stacks in current window
      float maxLagTime = 0.0;
      //      int angleIndex_maxTime = 0;
      int nAnglesOK = 0;
      double* val11 = new double[vars->nAngles];
      double* val22 = new double[vars->nAngles];
      for( int iangle = 0; iangle < vars->nAngles; iangle++ ) {
        // Cross-correlate S1 & S2 trace for each tested S1 angle
        compute_twosided_correlation2(
                          &vars->s1[iangle][startIndex],
                          &vars->s2[iangle][startIndex],
                          nSampIn,
                          corrBuffer,
                          vars->corrMaxlag_samples );
        int sampleIndex_maxAmp = 0;
        float maxAmp = corrBuffer[sampleIndex_maxAmp];
        // Determine maximum cross-correlation lag time & amplitude
        for( int isamp = 0; isamp < nSampCorr; isamp++ ) {
          if( corrBuffer[isamp] > maxAmp ) {
            maxAmp = corrBuffer[isamp];
            sampleIndex_maxAmp = isamp;
          }
        }
        // Interpolate maximum cross-correlation lag
        float sampleIndex_maxAmp_float = getQuadMaxSample( corrBuffer, sampleIndex_maxAmp, nSampCorr, &maxAmp );
        vars->corr_lags[iangle] = (sampleIndex_maxAmp_float-(float)vars->corrMaxlag_samples)*shdr->sampleInt;
        // Set
        if( vars->corr_lags[iangle] > maxLagTime ) {
          maxLagTime = vars->corr_lags[iangle];
          //          angleIndex_maxTime = iangle;
        }
        if( sampleIndex_maxAmp_float != 0.0f && sampleIndex_maxAmp_float != (float)(nSampCorr-1) ) {
          val11[nAnglesOK] = (float)(iangle)*vars->angleInc_rad;
          val22[nAnglesOK] = vars->corr_lags[iangle];
//          fprintf(stderr,"OK  %d %f %f   %f %d\n", iangle, RAD2DEG(val11[nAnglesOK]), val22[nAnglesOK], sampleIndex_maxAmp_float, vars->corrMaxlag_samples);
          nAnglesOK += 1;
        }
//        fprintf(stderr,"%d %d %f %f\n", winIndex, iangle, sampleIndex_maxAmp_float, vars->corr_lags[iangle] );
      } // END for( iangle )
      int periodicity = 2;
      double angleInc_deg = 1.0;
      double s1az_double, s2lag_double, stddev;
      computeXcorCos( val11, val22, nAnglesOK, periodicity, angleInc_deg, false, s1az_double, stddev, s2lag_double );
      float s1az  = (float)s1az_double;
      float s2lag = (float)s2lag_double;
      delete [] val11;
      delete [] val22;
//      fprintf(stderr,"\nMax lag time: %f %d   %f %lf\n", maxLagTime, angleIndex_maxTime, result, result2 );
/*
// Other method, seeking maximum amplitude. Less robust, for example in case cross-correlation max is cycle skipped
       float val1[3];
      float val2[3];
      for( int ival = 0; ival < 3; ival++ ) {
        int iangle = angleIndex_maxTime-1+ival;
        if( iangle < 0 ) {
          iangle += vars->nAngles;
        }
        else if( iangle >= vars->nAngles ) {
          iangle -= vars->nAngles;
        }
        val1[ival] = RAD2DEG( (float)(iangle)*vars->angleInc_rad );
        val2[ival] = vars->corr_lags[iangle];
      }

      float s2lag;
      float sampleIndex_maxAmp_float = getQuadMaxSample( val2, 1, vars->nAngles, &s2lag );
      int sampleIndex_int = (int)sampleIndex_maxAmp_float;
      if( sampleIndex_int == 2 ) sampleIndex_int = 1;
      float s1az  = (sampleIndex_maxAmp_float-(float)sampleIndex_int) * ( val1[sampleIndex_int+1]-val1[sampleIndex_int] ) + val1[sampleIndex_int];
*/
      if( vars->fout != NULL ) {
        int ensemble_id = traceGather->trace(0)->getTraceHeader()->intValue( vars->hdrId_identifier );
        fprintf(vars->fout, "%10d %12.2f %12.2f  %10.3f  %10.3f\n", ensemble_id, vars->winStart[winIndex], vars->winEnd[winIndex], s1az, s2lag);
      }
      writer->line("%4d %12.2f %12.2f  %10.3f  %10.3f", vars->ensembleCounter+1, vars->winStart[winIndex], vars->winEnd[winIndex], s1az, s2lag);
      if( edef->isDebug() ) fprintf(stderr,"Window  %f - %f  s1az %f  s2lag %f\n", vars->winStart[winIndex], vars->winEnd[winIndex], s1az, s2lag);
      // 3) Correct splitting in current window
      //  a) Rotate from azim to s1az
      //  b) Time stretch S2 back to S1
      //  c) Rotate back from s1az to azim
      float rotation_angle_rad = (float)(s1az * M_PI / 180.0) - vars->azim_rad;
      if( vars->rotationMode == mod_splitting::MODE_REMOVE ) {
        rotation_angle_rad *= -1.0f;
      }

      float times[2];
      times[0] = vars->winStart[winIndex];
      if( vars->stretchTopHalf ) {
        times[1] = 0.5f * ( vars->winStart[winIndex] + vars->winEnd[winIndex] );
      }
      else {
        times[1] = vars->winEnd[winIndex];
      }
      csTimeStretch timeStretchObj( shdr->sampleInt, shdr->numSamples );
      float stretch = -s2lag;
      for( int ipair = 0; ipair < numTotalPairs; ipair++ ) {
        int trueTraceIndex = ipair*2;
        float* traceXPtr = traceGather->trace(trueTraceIndex)->getTraceSamples();
        float* traceYPtr = traceGather->trace(trueTraceIndex+1)->getTraceSamples();
        // Rotate from azim to s1az
        rotate_xy( traceXPtr, traceYPtr, rotation_angle_rad, shdr->numSamples );
        // Time stretch/squeeze S2
        timeStretchObj.applyStretchFunction( traceYPtr, times, &stretch, 2, traceBuffer );
        memcpy( traceYPtr, traceBuffer, sizeof(float)*shdr->numSamples );
        // Rotate back from s1az to azim
        rotate_xy( traceXPtr, traceYPtr, -rotation_angle_rad, shdr->numSamples );
      }

      winIndex += 1;

      // Regenerate S1/S2 stacks after layer-stripping.
      compute_s1s2_splitting_stacks( samples, sr_azim, numSamples, numUsedPairs, vars->azim_rad,
        edef->isDebug(), writer->getFile(),
        vars->angleWidthOmit, vars->normMethod, vars->s2Scaling, vars->nAngles/2, vars->angleInc_rad,
        vars->s1, vars->s2 );

      // Trick: Only a quarter of a circle is actually needed to uniquely sample the directions (see nAngles/2 in function call above)
      // --> Copy redundant traces to full 180deg half circle for nicer output display and easier handling regarding cross-correlation etc
      for( int iangle = 0; iangle < vars->nAngles/2; iangle++ ) {
        int indexTo = iangle + vars->nAngles/2;
        memcpy( &vars->s1[indexTo][0], &vars->s2[iangle][0], sizeof(float)*shdr->numSamples );
        memcpy( &vars->s2[indexTo][0], &vars->s1[iangle][0], sizeof(float)*shdr->numSamples );
      }

      //******************************************************************************************
      // Write S1S2 stacks to output if requested
      //
      if( vars->outputOption == mod_splitting::OUTPUT_S1S2_STACKS || (vars->outputOption == mod_splitting::OUTPUT_LAST_S1S2_STACKS && (winIndex == vars->numWindows)) ) {
        for( int iangle = 0; iangle < vars->nAngles; iangle++ ) {
          int traceIndex = iangle*2 + startIndexOutputTrace;
          memcpy( traceGather->trace(traceIndex)->getTraceSamples(), vars->s1[iangle], shdr->numSamples*sizeof(float) );
          memcpy( traceGather->trace(traceIndex+1)->getTraceSamples(), vars->s2[iangle], shdr->numSamples*sizeof(float) );
        }
        startIndexOutputTrace += vars->nAngles*2;
      }
    } while( winIndex < vars->numWindows );
    if( corrBuffer != NULL ) delete [] corrBuffer;
    if( traceBuffer != NULL ) delete [] traceBuffer;
  } // END if layer_stripping

  // Remove original seismic traces
  if( vars->outputOption == mod_splitting::OUTPUT_LAST_S1S2_STACKS || vars->outputOption == mod_splitting::OUTPUT_S1S2_STACKS ) {
    traceGather->freeTraces( 0, numTotalTraces );
  }

  //******************************************************************************************

  delete [] samples;
  delete [] sr_azim;
  vars->ensembleCounter += 1;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_splitting_( csParamDef* pdef ) {
  pdef->setModule( "SPLITTING", "Shear-wave splitting analysis" );

  pdef->addDoc("This module generates stacked S1 and S2 traces, for each input ensemble and each test angle specified by the parameters 'num_angles' and 'angle_inc'.");
  pdef->addDoc("");
  pdef->addDoc("Input data must be one ensemble with pairs of XY traces. NMO etc has been applied, so that the result can be stacked");
  pdef->addDoc("A 2D right-handed XY coordinate system is assumed (e.g. X pointing North, and Y pointing East)");

  pdef->addParam( "layer_strip", "Perform layer stripping?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Generate S1/S2 stacks, no layer stripping" );
  pdef->addOption( "yes", "Perform layer-stripping", "Specify windows for layer stripping, data output option, and cross-correlation parameters" );

  pdef->addParam( "azim", "Sensor azimuth [deg]", NUM_VALUES_FIXED, "Azimuth of X sensor. Y sensor must be 90deg clock-wise (right-hand system)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Sensor azimuth [deg]" );

  pdef->addParam( "num_angles", "Number of angles to test in a 180deg half-circle", NUM_VALUES_FIXED,
    "The angle increment for S1/S2 analysis is computed as 180deg/num_angles" );
  pdef->addValue( "", VALTYPE_NUMBER, "" );

//  pdef->addParam( "angle_inc", "Angle increment", NUM_VALUES_FIXED );
//  pdef->addValue( "", VALTYPE_NUMBER, "[deg]" );

  pdef->addParam( "angle_omit", "Width of angle band to omit from S1/S2 stack", NUM_VALUES_FIXED, "Traces with source-receiver azimuth of 90deg +/-angle are omitted from S1 stack, and 180deg +/-angle are omitted from S2 stack" );
  pdef->addValue( "5", VALTYPE_NUMBER, "[deg]" );

  pdef->addParam( "norm_method", "Normalisation method for S1/S2 stacks", NUM_VALUES_FIXED );
  pdef->addValue( "ls", VALTYPE_OPTION );
  pdef->addOption( "ls", "Least-square normalisation" );
  pdef->addOption( "ntraces", "Normalise by number of traces" );
  pdef->addOption( "isotropic", "Normalise by number of traces, isotropic S2 scaling", "This normalisation enhances isotropic side scatterers on 'S2' stack" );

  pdef->addParam( "offset_range", "Min/max absolute offset range to use for splitting analysis/layer stripping", NUM_VALUES_VARIABLE );
  pdef->addValue( "0", VALTYPE_NUMBER, "Minimum absolute offset [m]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum absolute offset [m]", "Not specified = No maximum offset limitation" );

  pdef->addParam( "orient_xy", "Orientation of XY components", NUM_VALUES_FIXED, "Required for layer stripping" );
  pdef->addValue( "right", VALTYPE_OPTION );
  pdef->addOption( "right", "XY components describe a right hand system" );
  pdef->addOption( "left", "XY components describe a left hand system" );

  pdef->addParam( "corr_maxlag", "Parameters for cross-correlation when layer stripping", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Maximum time lag for cross-correlation" );

  pdef->addParam( "win_start", "List of start times defining layer stripping windows", NUM_VALUES_VARIABLE );
  pdef->addValue( "0", VALTYPE_NUMBER, "List of times in [ms]" );
  pdef->addParam( "win_end", "List of end times defining layer stripping windows", NUM_VALUES_VARIABLE );
  pdef->addValue( "0", VALTYPE_NUMBER, "List of times in [ms]", "0 = Use full time window of input data" );

  pdef->addParam( "output", "What data shall be output?", NUM_VALUES_FIXED );
  pdef->addValue( "s1s2_stacks", VALTYPE_OPTION );
  pdef->addOption( "s1s2_stacks", "Write S1/S2 analysis stacks to output traces" );
  pdef->addOption( "corrected_data", "Write layer-stripped, corrected input data to output traces" );
  pdef->addOption( "last_s1s2_stacks", "Write only last S1/S2 analysis stacks (layer-stripped, corrected) to output traces" );
//  pdef->addOption( "all_s1s2_stacks", "Write both S1/S2 analysis stacks AND layer-stripped, corrected input data to output traces",
//     "Trace header 'splitting_ens' for the two different outputs is set to 1 and 2, respectively" );

  pdef->addParam( "output_hdr", "How shall trace header of output trace be determined?", NUM_VALUES_FIXED );
  pdef->addValue( "first", VALTYPE_OPTION );
  pdef->addOption( "first", "Write trace header values of first input trace to output trace" );
  pdef->addOption( "last", "Write trace header values of last input trace to output trace" );
  pdef->addOption( "average", "Take average over all input trace headers for output trace" );

  pdef->addParam( "write_info", "Write layer stripping results to ASCII output file", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header that serves as unique identifier for each analysed gather" );
  pdef->addValue( "", VALTYPE_STRING, "File name of output ASCII file", "Output format:  <br>id_value win_start[ms] win_end[ms] s1az[deg] s2lag[ms]" );

  pdef->addParam( "stretch_top_half", "Apply time stretch from top to centre of analysis window?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "No, apply time stretch over while analysis window" );
  pdef->addOption( "yes", "Yes, only apply time stretch over top half of analysis window" );
}



//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_splitting_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_splitting::VariableStruct* vars = reinterpret_cast<mod_splitting::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_splitting_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_splitting::VariableStruct* vars = reinterpret_cast<mod_splitting::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->s1 != NULL && vars->s2 != NULL ) {
    for( int iangle = 0; iangle < vars->nAngles; iangle++ ) {
      delete [] vars->s1[iangle];
      delete [] vars->s2[iangle];
    }
    delete [] vars->s1;
    delete [] vars->s2;
    vars->s1 = NULL;
    vars->s2 = NULL;
  }
  if( vars->winStart != NULL && vars->winEnd != NULL ) {
    delete [] vars->winStart;
    vars->winStart = NULL;
    delete [] vars->winEnd;
    vars->winEnd = NULL;
  }
  if( vars->corr_lags != NULL ) {
    delete [] vars->corr_lags;
    vars->corr_lags = NULL;
  }
  if( vars->fout != NULL ) {
    fclose(vars->fout);
    vars->fout = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_splitting_( csParamDef* pdef ) {
  params_mod_splitting_( pdef );
}
extern "C" void _init_mod_splitting_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_splitting_( param, env, writer );
}
extern "C" bool _start_exec_mod_splitting_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_splitting_( env, writer );
}
extern "C" void _exec_mod_splitting_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_splitting_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_splitting_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_splitting_( env, writer );
}

void compute_s1s2_splitting_stacks( float** samples, float const* sr_azim, int numSamples, int numPairs, float azim,
    bool isDebug, FILE* logFile,
    float angleWidthOmit, int normMethod, int s2Scaling, int nAngles, float angleInc,
    float** s1, float** s2 ) {

  //  fprintf(stdout,"BB %d %d  %d  %f  %f\n", numSamples, numPairs, isDebug, angleWidthOmit, angleInc );
  float rad_90deg  = (float)( 0.5 * M_PI );
  float rad_180deg = (float)M_PI;
  float rad_270deg = (float)( 1.5 * M_PI );
  float rad_360deg = (float)( 2.0 * M_PI );

  for( int iangle = 0; iangle < nAngles; iangle++ ) {
    float s1az = (float)(iangle)*angleInc;
    float angle_xy_to_s1s2 = fmod( s1az - azim + rad_360deg, rad_360deg );

    for( int isamp = 0; isamp < numSamples; isamp++ ) {
      s1[iangle][isamp] = 0;
      s2[iangle][isamp] = 0;
      double sum_cos  = 0.0;
      double sum_cos2 = 0.0;
      double sum_sin  = 0.0;
      double sum_sin2 = 0.0;
      // TEMP
      double sumS1 = 0.0;
      double sumS2 = 0.0;
      double sum_sin1 = 0.0;
      double sum_cos1 = 0.0;

      int counterPairsS1 = 0;
      int counterPairsS2 = 0;

      for( int ip = 0; ip < numPairs; ip++ ) {
        int trcX = ip*2; // 'X' trace index
        int trcY = trcX+1; // 'Y' trace index
//        fprintf(stdout,"%d %d    %d\n", trcX, trcY, numPairs);
//        fflush(stdout);
        float ampX = samples[trcX][isamp];
        float ampY = -samples[trcY][isamp];   // Take negative Y amplitude due to OBC definition of XYZ coordinate system!
//        fprintf(stdout,"        %f %f\n", ampX, ampY);
//        fflush(stdout);
        if( ampX == 0.0 && ampY == 0.0 ) continue;
        float ampS1  = ampX*cos(angle_xy_to_s1s2) - ampY*sin(angle_xy_to_s1s2);
        float ampS2  = ampX*sin(angle_xy_to_s1s2) + ampY*cos(angle_xy_to_s1s2);

        float angle_srazim_to_s1s2 = s1az - sr_azim[ip];

        double cos_angle = cos(angle_srazim_to_s1s2);
        double sin_angle = sin(angle_srazim_to_s1s2);

        angle_srazim_to_s1s2 = fmod( angle_srazim_to_s1s2 + rad_360deg, rad_360deg );
        float sign_s1 = 1.0;
        float sign_s2 = 1.0;
        if( angle_srazim_to_s1s2 > rad_90deg && angle_srazim_to_s1s2 <= rad_270deg ) {
          sign_s1 = -1.0;
        }
        if( angle_srazim_to_s1s2 > rad_180deg ) {
          sign_s2 = -1.0;
        }
        if( s2Scaling == mod_splitting::S2_SCALING_ISOTROPIC ) {
          sign_s2 = sign_s1;
        }
        if( (angle_srazim_to_s1s2 < rad_90deg-angleWidthOmit) || (angle_srazim_to_s1s2 > rad_270deg+angleWidthOmit) ||
            (angle_srazim_to_s1s2 > rad_90deg+angleWidthOmit && angle_srazim_to_s1s2 < rad_270deg-angleWidthOmit) ) {
          sum_cos   += ampS1 * cos_angle;
          sum_cos2  += cos_angle*cos_angle;

          sumS1     += sign_s1*ampS1;
          sum_cos1  += cos_angle;
          counterPairsS1 += 1;
        }

        if( (angle_srazim_to_s1s2 > angleWidthOmit && angle_srazim_to_s1s2 < rad_180deg-angleWidthOmit) ||
            (angle_srazim_to_s1s2 > rad_180deg+angleWidthOmit && angle_srazim_to_s1s2 < rad_360deg-angleWidthOmit) ) {
          sum_sin   += ampS2 * sin_angle;
          sum_sin2  += sin_angle*sin_angle;

          sumS2     += sign_s2*ampS2;
          sum_sin1  += sin_angle;
          counterPairsS2 += 1;
        }
      }
      if( normMethod == mod_splitting::NORM_METHOD_NTRACES ) {
        if( counterPairsS1 != 0 && counterPairsS2 != 0 ) {
          s1[iangle][isamp] = (float)sumS1/(float)counterPairsS1;
          s2[iangle][isamp] = (float)sumS2/(float)counterPairsS2;
        }
      }
      else {
        if( sum_cos2 != 0.0 && sum_sin2 != 0.0 ) {
          s1[iangle][isamp] = (float)( sum_cos/sum_cos2 );
          s2[iangle][isamp] = (float)( sum_sin/sum_sin2 );
        }
      }
    }
  }

}

void compute_twosided_correlation2( float const* samplesLeft, float const* samplesRight, int nSampIn, float* corr, int maxlag_in_num_samples ) {

  int sampStart;
  int sampEnd;

  //---------------------------------------
  // Compute negative lags
  //
  sampEnd    = nSampIn;
  for( int ilag = -maxlag_in_num_samples; ilag < 0; ilag++ ) {
    sampStart  = -ilag;
    float sum = 0;
    for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
      sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
    }
    corr[ilag+maxlag_in_num_samples] = sum;
  }

  //---------------------------------------
  // Compute positive lags
  //
  sampStart  = 0;
  for( int ilag = 0; ilag <= maxlag_in_num_samples; ilag++ ) {
    sampEnd    = nSampIn-ilag;
    float sum = 0;
    for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
      sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
    }
    corr[ilag+maxlag_in_num_samples] = sum;
  }
}



/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>
#include <cstring>
#include "hdr_math_ens_defines.h"
#include "csEquationSolver.h"
#include "csSort.h"
#include "csInterpolation.h"
#include "geolib_math.h"
#include "geolib_methods.h"

using std::string;
using namespace cseis_geolib;
using namespace cseis_system;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TRC_MATH_ENS
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_trc_math_ens {
  struct VariableStruct {
    int  method;
    int numEns;
    int hdrId_ens;
    int hdrId_dim2;
    int num_dim2;
    mod_trc_math_ens::Ens* ens;
    cseis_system::csTraceGather* gather;
    bool isFirstCall;
    cseis_geolib::csEquationSolver* solver;
    double* userVarValues;
    int  numVariables;
    int* varIndexList;
    int medianFilt_numTraces;
    int medianFilt_numTracesALL;
    int medianFilt_output;
    bool medianFilt_applyAtEdges;
    int medianFilt_maxSample;
    cseis_geolib::csSort<float>* medianFilt_sortObj;
    int despike_winSpike_traces;
    int despike_winRef_traces;
    float despike_ratio;
    bool despike_apply;
  };
  static const int METHOD_DEBIAS = 1;
  static const int METHOD_ROLL_MEAN   = 2;
  static const int METHOD_EQUATION = 3;
  static const int METHOD_MIN   = 4;
  static const int METHOD_MAX   = 5;
  static const int METHOD_MEAN  = 6;
  static const int METHOD_DIFF  = 7;
  static const int METHOD_MEDIAN_FILT  = 8;
  static const int METHOD_DESPIKE      = 9;
  static const int METHOD_DIFF_ADAPTIVE  = 10;
  static const int METHOD_ROLL_MEDIAN   = 11;

  static int const OUTPUT_FILT = 41;
  static int const OUTPUT_DIFF = 42;
}
using namespace mod_trc_math_ens;

void computeEnsemble( int indexEns, cseis_system::csTraceGather* gatherIn, mod_trc_math_ens::Ens* ens, int numSamples, int numDim2, csTraceGather* gatherOut );

// ERROR: Wrong firstIndex or nElements passed to vector object, remove(firstIndex,nElements) function
// 

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_trc_math_ens_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->method   = METHOD_DEBIAS;
  vars->hdrId_ens  = -1;
  vars->hdrId_dim2 = -1;
  vars->numEns = 0;
  vars->num_dim2 = 0;
  vars->ens = NULL;
  vars->gather = NULL;
  vars->isFirstCall = true;
  vars->solver = NULL;
  vars->numVariables = 0;
  vars->userVarValues = NULL;
  vars->varIndexList = NULL;
  vars->medianFilt_numTracesALL = 0;
  vars->medianFilt_numTraces    = 0;
  vars->medianFilt_sortObj      = NULL;
  vars->medianFilt_applyAtEdges = false;
  vars->medianFilt_maxSample = -1;

  vars->despike_winSpike_traces  = 1;
  vars->despike_winRef_traces    = 5;
  vars->despike_ratio = 3.0;
  vars->despike_apply = true;

  int numTracesFixed = 0;

  //----------------------------------------------------------------
  std::string text;
  if( param->exists("mode") ) {
    param->getString("mode", &text);
    if( !text.compare("ensemble") ) {
      edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
    }
    else if( !text.compare("fixed") ) {
      param->getInt("ntraces", &numTracesFixed);
      if( numTracesFixed < 1 ) writer->error("Incorrect number of traces specified in user parameter 'ntraces' : %d", numTracesFixed);
      edef->setTraceSelectionMode( TRCMODE_FIXED, numTracesFixed );
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  else {
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }

  param->getString( "method", &text );
  text = toLowerCase( text );
  if( !text.compare("debias") ) {
    vars->method = METHOD_DEBIAS;
    //    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE, 2 );
  }
  else if( !text.compare("mean") ) {
    vars->method = METHOD_MEAN;
  }
  else if( !text.compare("median_filt") ) {
    vars->method = METHOD_MEDIAN_FILT;
  }
  else if( !text.compare("roll_mean") ) {
    vars->method = METHOD_ROLL_MEAN;
  }
  else if( !text.compare("roll_median") ) {
    vars->method = METHOD_ROLL_MEDIAN;
  }
  else if( !text.compare("min") ) {
    vars->method = METHOD_MIN;
  }
  else if( !text.compare("max") ) {
    vars->method = METHOD_MAX;
  }
  else if( !text.compare("diff") ) {
    vars->method = METHOD_DIFF;
  }
  else if( !text.compare("diff_adapt") ) {
    vars->method = METHOD_DIFF_ADAPTIVE;
  }
  else if( !text.compare("equation") ) {
    vars->method = METHOD_EQUATION;
  }
  else if( !text.compare("despike") ) {
    vars->method = METHOD_DESPIKE;
    param->getInt("despike",&vars->despike_winSpike_traces,0);
    param->getInt("despike",&vars->despike_winRef_traces,1);
    vars->despike_winRef_traces   = (int)( vars->despike_winRef_traces/2 )*2 + 1;
    vars->despike_winSpike_traces = (int)( vars->despike_winSpike_traces/2 )*2 + 1;
    param->getFloat("despike",&vars->despike_ratio,2);
    param->getString("despike",&text,3);
    if( !text.compare("yes") ) {
      vars->despike_apply = true;
    }
    else if( !text.compare("ratio") ) {
      vars->despike_apply = false;
    }
    else {
      writer->line("Option not recognized: %s.", text.c_str());
      env->addError();
    }
  }
  else {
    writer->line("Method not recognized: %s.", text.c_str());
    env->addError();
  }

  if( vars->method == mod_trc_math_ens::METHOD_ROLL_MEAN || vars->method == mod_trc_math_ens::METHOD_ROLL_MEDIAN ) {
    param->getString("roll_ens",&text);
    param->getInt("roll_ens",&vars->numEns,1);
    vars->hdrId_ens = hdef->headerIndex(text);
    vars->numEns = 2*vars->numEns + 1;

    param->getString("roll_dim2",&text);
    param->getInt("roll_dim2",&vars->num_dim2,1);
    vars->hdrId_dim2 = hdef->headerIndex(text);
    vars->num_dim2 = vars->num_dim2;

    vars->ens = new mod_trc_math_ens::Ens(vars->numEns);
    vars->gather = new cseis_system::csTraceGather( hdef );

    writer->line("# ensembles in filter: %d", vars->numEns);
    writer->line("# traces in filter:    %d", 2*vars->num_dim2+1);
  }
  else if( vars->method == mod_trc_math_ens::METHOD_EQUATION ) {
    std::string equationText;
    param->getString( "equation", &equationText );

    char* name = new char[4];
    int numTraces = 0;
    for( int i = 1; i < 10; i++ ) {
      sprintf(name,"x%1d",i);
      //      fprintf(stderr,"Testing '%s' in '%s'\n", name, equationText.c_str() );
      if( equationText.find(name) == string::npos ) {
        break;
      }
      numTraces += 1;
    }
    writer->line("Number of trace variables (xN) found in equation: %d", numTraces);
    if( numTraces < 2 ) {
      writer->error("Did not find at least two trace variables (x1 & x2) in equation");
    }
    else if( numTraces > 9 ) {
      writer->error("Maximum number of traces in gather that can be referenced in equation (xN) is 9");
    }
    std::string* allTraceNames = new std::string[numTraces];
    for( int i = 0; i < numTraces; i++ ) {
      sprintf(name,"x%1d",i+1);
      allTraceNames[numTraces-i-1] = name;
    }
    vars->solver = new csEquationSolver();
    if( !vars->solver->prepare( equationText, allTraceNames, numTraces ) ) {
      writer->error("Error occurred: %s", vars->solver->getErrorMessage().c_str() );
    }
    delete [] allTraceNames;

    csVector<string> constList;
    vars->solver->prepareUserConstants( &constList );

    vars->numVariables = constList.size();
    vars->userVarValues = new double[vars->numVariables];
    vars->varIndexList = new int[vars->numVariables];

    if( vars->numVariables != numTraces ) {
      writer->error("Inconsistent number of input traces. Expected %d, found %d.\n", vars->numVariables, numTraces);
    }
    else if ( numTraces > numTracesFixed ) {
      writer->error("Number of traces accessed in equation (=%d) exceeds number of traces specified as input ensemble (=%d)", numTraces, numTracesFixed);
    }
    if( numTracesFixed == 0 ) edef->setTraceSelectionMode( TRCMODE_FIXED, numTraces );

    for( int ivar = 0; ivar < vars->numVariables; ivar++ ) {
      writer->line("Variable #%d: %s", ivar+1, constList.at(ivar).c_str());
      if( constList.at(ivar).at(0) != 'x' ) {
        writer->error("Unknown variable name in equation: '%s'. Use variable 'xN' to reference trace N.", constList.at(ivar).c_str() );
        env->addError();
      }
      if( edef->isDebug() ) writer->line("Variable #%d: %s", ivar, constList.at(ivar).c_str() );
      char numChar = constList.at(ivar).at(1);
      int traceNum = numChar - 48;
      if( traceNum < 1 || traceNum > numTraces ) {
        writer->error("Unknown variable name in equation: '%s'. Use variable 'xN' to reference trace N.", constList.at(ivar).c_str() );
        env->addError();
      }
      vars->varIndexList[ivar] = traceNum-1;
    }
  }

  if( vars->method == mod_trc_math_ens::METHOD_MEDIAN_FILT || vars->method == mod_trc_math_ens::METHOD_ROLL_MEDIAN ) {
    vars->medianFilt_numTracesALL = 11;
    if( param->exists("median_filt") ) {
      param->getInt( "median_filt", &vars->medianFilt_numTracesALL );
      if( vars->medianFilt_numTracesALL < 3 ) {
        vars->medianFilt_numTracesALL = 3;
      }
      else if( (int)(vars->medianFilt_numTracesALL/2)*2 == vars->medianFilt_numTracesALL ) {
        vars->medianFilt_numTracesALL += 1;
      }
      param->getString("median_filt", &text, 1 );
      if( !text.compare("filt") ) {
        vars->medianFilt_output = OUTPUT_FILT;
      }
      else if( !text.compare("diff") ) {
        vars->medianFilt_output = OUTPUT_DIFF;
      }
      else {
        writer->error("Unknown option: %s", text.c_str() );
      }
      if( param->getNumValues( "median_filt" ) > 2 ) {
        param->getString("median_filt", &text, 2 );
        if( !text.compare("no") ) {
          vars->medianFilt_applyAtEdges = false;
        }
        else if( !text.compare("yes") ) {
          vars->medianFilt_applyAtEdges = true;
        }
        else {
          writer->error("Unknown option: %s", text.c_str() );
        }
      }
      if( param->getNumValues( "median_filt" ) > 3 ) {
        float maxTime;
        param->getFloat("median_filt", &maxTime, 3 );
        vars->medianFilt_maxSample = (maxTime <= 0) ? -1 : std::min( shdr->numSamples-1, (int)(maxTime / shdr->sampleInt) );
      }
    }
    vars->medianFilt_numTraces = (vars->medianFilt_numTracesALL - 1 )/2;
    vars->medianFilt_sortObj   = new csSort<float>();
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_trc_math_ens_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;

  int nTraces = traceGather->numTraces();

  //  fprintf(stdout,"TRC_MATH_ENS %d: Num traces %d, lastcall %d\n", vars->method, nTraces, edef->isLastCall());

  if( vars->method == mod_trc_math_ens::METHOD_DIFF ) {
    // Compute difference between first trace and all following traces
    float* samples1 = traceGather->trace(0)->getTraceSamples();
    for( int itrc = 1; itrc < nTraces; itrc++ ) {
      float* samples2 = traceGather->trace(itrc)->getTraceSamples();
      for( int isamp = 0; isamp < shdr->numSamples; isamp++) {
        samples1[isamp] -= samples2[isamp];
      }
    }
    traceGather->freeTraces( 1, nTraces-1 );
  }
  else if( vars->method == mod_trc_math_ens::METHOD_DIFF_ADAPTIVE ) {
    float* samples1 = traceGather->trace(0)->getTraceSamples();
    float* samples2 = traceGather->trace(1)->getTraceSamples();
    int maxLag_samples = 25;
    int startSamp = 250;
    int endSamp =   1000;
    int numSamp = endSamp - startSamp + 1;
    int numSampBuffer = 2*maxLag_samples + 1;
    float* buffer = new float[shdr->numSamples];
    cseis_geolib::csInterpolation interpol( shdr->numSamples, shdr->sampleInt, 8 );

    compute_twosided_correlation( &samples1[startSamp], &samples2[startSamp], numSamp, buffer, maxLag_samples, false );
    float maxAmp = buffer[0];
    int sampleIndex_maxAmp = 0;
    for( int isamp = 0; isamp < numSampBuffer; isamp++ ) {
      if( buffer[isamp] > maxAmp ) {
        maxAmp = buffer[isamp];
        sampleIndex_maxAmp = isamp;
      }
    }

    float rms1 = 0;
    float rms2 = 0;
    for( int isamp = startSamp; isamp <= endSamp; isamp++ ) {
      rms1 += samples1[isamp]*samples1[isamp];
      rms2 += samples2[isamp]*samples2[isamp];
    }
    rms1 = sqrt(rms1/(float)numSamp);
    rms2 = sqrt(rms2/(float)numSamp);
    float sampleIndex_maxAmp_float = getQuadMaxSample( buffer, sampleIndex_maxAmp, numSampBuffer, &maxAmp );
    float timeShift_ms = ( maxLag_samples - sampleIndex_maxAmp_float ) * shdr->sampleInt;

    //    fprintf(stderr,"Time shift: %.2f  %f   %f %f\n", timeShift_ms, rms1/rms2, rms1, rms2);
    interpol.static_shift( timeShift_ms, samples2, buffer );
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples1[isamp] -= buffer[isamp];
    }

    delete [] buffer;
    traceGather->freeTraces( 1, nTraces-1 );
  } 
  else if( vars->method == mod_trc_math_ens::METHOD_DEBIAS || vars->method == mod_trc_math_ens::METHOD_MIN ||
           vars->method == mod_trc_math_ens::METHOD_MAX || vars->method == mod_trc_math_ens::METHOD_MEAN ) {
    float** tracePtr = new float*[nTraces];
    for( int itrc = 0; itrc < nTraces; itrc++) {
      tracePtr[itrc] = traceGather->trace(itrc)->getTraceSamples();
    }
    if( vars->method == mod_trc_math_ens::METHOD_DEBIAS || vars->method == mod_trc_math_ens::METHOD_MEAN ) {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++) {
        float mean = 0.0;
        for( int itrc = 0; itrc < nTraces; itrc++) {
          mean += tracePtr[itrc][isamp];
        }
        mean /= (float)nTraces;
        if( vars->method == mod_trc_math_ens::METHOD_DEBIAS ) {
          for( int itrc = 0; itrc < nTraces; itrc++) {
            tracePtr[itrc][isamp] -= mean;
          }
        }
        else { // vars->method == mod_trc_math_ens::METHOD_MEAN ) {
          tracePtr[0][isamp] = mean;
        }
      }  
    }
    else if( vars->method == mod_trc_math_ens::METHOD_MIN ) {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++) {
        float valmin = tracePtr[0][isamp];
        for( int itrc = 1; itrc < nTraces; itrc++) {
          if( tracePtr[itrc][isamp] < valmin ) valmin = tracePtr[itrc][isamp];
        }
        tracePtr[0][isamp] = valmin;
      }
    }
    else if( vars->method == mod_trc_math_ens::METHOD_MAX ) {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++) {
        float valmax = tracePtr[0][isamp];
        for( int itrc = 1; itrc < nTraces; itrc++) {
          if( tracePtr[itrc][isamp] > valmax ) valmax = tracePtr[itrc][isamp];
        }
        tracePtr[0][isamp] = valmax;
      }
    }
    delete [] tracePtr;
    if( vars->method == mod_trc_math_ens::METHOD_MIN || vars->method == mod_trc_math_ens::METHOD_MAX || vars->method == mod_trc_math_ens::METHOD_MEAN ) {
      traceGather->freeTraces( 1, nTraces-1 );
    }
  }
  //--------------------------------------------------------------------------------
  else if( vars->method == mod_trc_math_ens::METHOD_ROLL_MEAN ) {
    if( nTraces != 0 ) {
      double ensValue = traceGather->trace( nTraces-1 )->getTraceHeader()->doubleValue( vars->hdrId_ens );
      vars->ens->addEns( nTraces, ensValue );
      traceGather->moveTracesTo( 0, nTraces, vars->gather );
    }

    if( vars->ens->numEns() == vars->numEns || edef->isLastCall() ) {
      int indexCenterEns = (int)( vars->ens->numEns()/2 );
      if( vars->isFirstCall ) {
        // Compute first half of collected traces (done only once for first call to mean filter)
        for( int iens = 0; iens < indexCenterEns; iens++ ) {
          computeEnsemble( iens, vars->gather, vars->ens, shdr->numSamples, vars->num_dim2, traceGather );
        }
      }
      if( !edef->isLastCall() ) {
        computeEnsemble( indexCenterEns, vars->gather, vars->ens, shdr->numSamples, vars->num_dim2, traceGather );
        int nTracesFirstEns = vars->ens->numTraces(0);
        vars->ens->releaseFirstEns(); // Release first ensemble, not needed anymore
        vars->gather->freeTraces(0,nTracesFirstEns);
      }
      else { //if( edef->isLastCall() ) {
        int numEns = vars->ens->numEns();
        for( int iens = indexCenterEns; iens < numEns; iens++ ) {
          computeEnsemble( iens, vars->gather, vars->ens, shdr->numSamples, vars->num_dim2, traceGather );
        }
        vars->gather->freeAllTraces();

        for( int iens = 0; iens < numEns; iens++ ) {
          vars->ens->releaseFirstEns(); // Release first ensemble, not needed anymore
        }
      }
      vars->isFirstCall = false;
    }
    if( vars->gather->numTraces() > 0 ) edef->setTracesAreWaiting();
  }
  else if( vars->method == mod_trc_math_ens::METHOD_MEDIAN_FILT ) {
    if( nTraces < vars->medianFilt_numTracesALL ) return;
    float* bufferIn      = new float[nTraces];
    float* bufferOut     = new float[nTraces];
    float* buffer_filter = new float[nTraces];
    float** samplesPtr   = new float*[nTraces];
    int trc1_out = 0;
    int trc2_out = nTraces-1;
    if( !vars->medianFilt_applyAtEdges ) {
      trc1_out = vars->medianFilt_numTraces;
      trc2_out = nTraces - vars->medianFilt_numTraces - 1;
    }
    for( int itrc = 0; itrc < nTraces; itrc++ ) {
      samplesPtr[itrc] = traceGather->trace(itrc)->getTraceSamples();
    }
    //--------------------------------------------------------------------------------
    int lastSample = (vars->medianFilt_maxSample <= 0) ? shdr->numSamples-1 : vars->medianFilt_maxSample;

    for( int isamp = 0; isamp <= lastSample; isamp++ ) {
      int median_numTraces = vars->medianFilt_numTraces;
      int median_numTracesALL = vars->medianFilt_numTracesALL;
      if( vars->medianFilt_maxSample > 0 ) {
        median_numTracesALL = (int)( ( (float)vars->medianFilt_numTracesALL - 3.0 ) * ( 1.0 - (float)isamp/(float)vars->medianFilt_maxSample ) ) + 3;
        if( median_numTracesALL < 3 ) {
          median_numTracesALL = 3;
        }
        median_numTraces = ( median_numTracesALL - 1 )/2;
      }
      for( int itrc = 0; itrc < nTraces; itrc++ ) {
        bufferIn[itrc] = samplesPtr[itrc][isamp];
      }
      for( int itrc = median_numTraces; itrc < nTraces-median_numTraces; itrc++ ) {
        int index1 = itrc - median_numTraces;
        memcpy( buffer_filter, &bufferIn[index1], median_numTracesALL*sizeof(float) );
        vars->medianFilt_sortObj->treeSort( buffer_filter, median_numTracesALL );
        bufferOut[itrc] = buffer_filter[median_numTraces];
      }
      // Edges:
      if( vars->medianFilt_applyAtEdges ) {
        for( int itrc = 0; itrc < median_numTraces; itrc++ ) {
          int ntrc2 = std::max( itrc + median_numTraces, 3 );
          memcpy( buffer_filter, bufferIn, ntrc2*sizeof(float) );
          vars->medianFilt_sortObj->treeSort( buffer_filter, ntrc2 );
          bufferOut[itrc] = buffer_filter[ntrc2/2];
        }
        for( int itrc = std::max(0,nTraces-median_numTraces); itrc < nTraces; itrc++ ) {
          int index1 = std::max( 0, itrc - median_numTraces );
          int ntrc2  = nTraces - index1 - 1;
          memcpy( buffer_filter, &bufferIn[index1], ntrc2*sizeof(float) );
          vars->medianFilt_sortObj->treeSort( buffer_filter, ntrc2 );
          bufferOut[itrc] = buffer_filter[ntrc2/2];
        }
      }
      if( vars->medianFilt_output == mod_trc_math_ens::OUTPUT_FILT ) {
        for( int itrc = trc1_out; itrc <= trc2_out; itrc++ ) {
          samplesPtr[itrc][isamp] = bufferOut[itrc];
        }
      }
      else {
        for( int itrc = trc1_out; itrc <= trc2_out; itrc++ ) {
          samplesPtr[itrc][isamp] -= bufferOut[itrc];
        }
      }
    } // END for isamp
    if( !vars->medianFilt_applyAtEdges && vars->medianFilt_output == mod_trc_math_ens::OUTPUT_DIFF ) {
      for( int itrc = 0; itrc < trc1_out; itrc++ ) {
        float* sPtr_trace = samplesPtr[itrc];
        for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
          sPtr_trace[isamp] = 0;
        }
      }
      for( int itrc = trc2_out+1; itrc < nTraces; itrc++ ) {
        float* sPtr_trace = samplesPtr[itrc];
        for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
          sPtr_trace[isamp] = 0;
        }
      }
    }
    delete [] bufferIn;
    delete [] bufferOut;
    delete [] buffer_filter;
    delete [] samplesPtr;
    //    fprintf(stdout,"EXIT MEDIAN numTraces out: %d\n", traceGather->numTraces());
  }
  else if( vars->method == METHOD_EQUATION ) {
    if( vars->numVariables != traceGather->numTraces() ) {
      writer->warning("Inconsistent number of input traces (at end of file). Expected %d, found %d. Killed traces.\n",
                      vars->numVariables, traceGather->numTraces());
      traceGather->freeAllTraces();
      return;
    }
    float* samplesOut = traceGather->trace(0)->getTraceSamples();
    int nSamples = shdr->numSamples;
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      for( int ivar = 0; ivar < vars->numVariables; ivar++ ) {
        float* samples = traceGather->trace(ivar)->getTraceSamples();
        vars->userVarValues[ vars->varIndexList[ivar] ] = (double)samples[isamp];
      }
      vars->solver->setUserConstants( vars->userVarValues, vars->numVariables );
      samplesOut[isamp] = (float)vars->solver->solve();
    }
    traceGather->freeTraces( 1, traceGather->numTraces()-1 );
  }
  else if( vars->method == mod_trc_math_ens::METHOD_DESPIKE ) {
    float* buffer = new float[nTraces];

    float ratioNumTraces = (float)vars->despike_winRef_traces / (float)vars->despike_winSpike_traces;
    int trcWinMid   = (int)( vars->despike_winRef_traces/2 );
    int trc1WinSpike = trcWinMid - (int)( vars->despike_winSpike_traces/2 );
    int trcIndexLast   = nTraces - 1 - vars->despike_winRef_traces;

    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      for( int itrc = 0; itrc < nTraces; itrc++ ) {
        buffer[itrc] = traceGather->trace(itrc)->getTraceSamples()[isamp] + 1e-11;
      }

      float sumWinRef   = 0;
      float sumWinSpike = 0;
      for( int itrc = 0; itrc < vars->despike_winRef_traces; itrc++ ) {
        sumWinRef += buffer[itrc];
      }
      for( int itrc = 0; itrc < vars->despike_winSpike_traces; itrc++ ) {
        sumWinSpike += buffer[itrc+trc1WinSpike];
      }
      int traceIndex = 0;
      while( traceIndex <= trcIndexLast ) {
        float ratio = ( sumWinSpike / sumWinRef ) * ratioNumTraces;
        if( vars->despike_apply ) {
          if( ratio > vars->despike_ratio ) traceGather->trace(traceIndex+trcWinMid)->getTraceSamples()[isamp] = 0;
        }
        else {
          traceGather->trace(traceIndex+trcWinMid)->getTraceSamples()[isamp] = ratio;
        }
        sumWinRef   -= buffer[traceIndex];
        sumWinSpike -= buffer[traceIndex+trc1WinSpike];
        sumWinRef   += buffer[traceIndex+vars->despike_winRef_traces];
        sumWinSpike += buffer[traceIndex+trc1WinSpike+vars->despike_winSpike_traces];
        traceIndex += 1;
      }
    } // END for isamp
    delete [] buffer;
  }
}

//********************************************************************************
//
//
void params_mod_trc_math_ens_( csParamDef* pdef ) {
  pdef->setModule( "TRC_MATH_ENS", "Multi-trace sample computation", "Apply multi-trace mathematical equation to sample values" );

  pdef->addParam( "method", "", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "debias", "Debias trace samples: For each sample time, compute mean value across all traces and remove it from trace sample values" );
  pdef->addOption( "roll_mean", "Rolling mean filter, applied on a sample-by-sample basis" );
  pdef->addOption( "equation", "General math equation. Specify under parameter 'equation'" );
  pdef->addOption( "min", "Output minimum value, applied on a sample-by-sample basis" );
  pdef->addOption( "max", "Output maximum value, applied on a sample-by-sample basis" );
  pdef->addOption( "mean", "Output mean value, applied on a sample-by-sample basis" );
  pdef->addOption( "diff", "Take difference between consecutive traces. For each input ensemble, output 1 trace", "If the ensemble contains more than 2 traces, the difference is computed as follows: diff = trc1 - (trc2 + trc3 + ...)" );
  //  pdef->addOption( "diff_adapt", "Take difference between consecutive traces using adaptive subtraction", "Specify 'ntraces' for number of samples to use in window" );
  pdef->addOption( "median_filt", "Apply median filter, applied on a sample-by-sample basis across the input ensemble" );
  pdef->addOption( "despike", "Apply despike in spatial direction" );
  // pdef->addOption( "roll_median", "Rolling median filter, applied on a sample-by-sample basis" );

  pdef->addParam( "roll_ens", "", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Ensemble trace header" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of ensembles to roll. Specify 0 to work on a single ensemble at a time (no smearing between ensembles)" );

  pdef->addParam( "roll_dim2", "", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header" );
  pdef->addValue( "3", VALTYPE_NUMBER, "Number of traces to roll. Specify 0 to work on a single trace (=no filtering within current ensemble)" );

  pdef->addParam( "equation", "Mathematical equation on N consecutive traces", NUM_VALUES_FIXED, "One trace will be output" );
  pdef->addValue( "", VALTYPE_STRING, "Arbitrary mathematical equation to apply to each individual trace sample. Use 'xN' to reference sample value of trace N.",
                  "Constants: pi,e. Functions: abs,acos,asin,atan,atan2,ceil,cos,cosh,exp,floor,log,log10,max,min,mod,pow,int,round,sin,sinh,sqrt,tan,tanh,todegrees,toradians,sign");

  pdef->addParam( "mode", "Mode of operation: Input data gathering", NUM_VALUES_FIXED );
  pdef->addValue( "ensemble", VALTYPE_OPTION );
  pdef->addOption( "ensemble", "Input one ensemble at a time" );
  pdef->addOption( "fixed", "Input a fixed number of traces at a time (specified in 'ntraces')" );

  pdef->addParam( "ntraces", "Number of consecutive traces to input at once", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Number of traces" );

  pdef->addParam( "median_filt", "Parameters for median filter", NUM_VALUES_VARIABLE );
  pdef->addValue( "11", VALTYPE_NUMBER, "Total length of median filter in number of traces", "This must be an odd number, minimum 3" );
  pdef->addValue( "filt", VALTYPE_OPTION );
  pdef->addOption( "filt", "Output filtered data" );
  pdef->addOption( "diff", "Output difference between input and filtered data" );
  pdef->addValue( "no", VALTYPE_OPTION, "Apply filter at edges?" );
  pdef->addOption( "no", "Do not apply filter at edges" );
  pdef->addOption( "yes", "Apply filter at edges" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Taper application: Maximum time where median filter is applied. Filter length will be gradually tapered down", "Specify 0 to apply full median filter over full trace" );

  pdef->addParam( "despike", "Parameters for despike process", NUM_VALUES_FIXED, "Sets samples to 0" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Width of despike window [#traces]" );
  pdef->addValue( "5", VALTYPE_NUMBER, "Width of reference window [#traces]" );
  pdef->addValue( "3.0", VALTYPE_NUMBER, "Maximum ratio between despike & reference window" );
  pdef->addValue( "yes", VALTYPE_OPTION, "Apply despike?" );
  pdef->addOption( "yes", "Output data with despike applied" );
  pdef->addOption( "ratio", "Output RMS ratio between spike and reference window" );
  //  pdef->addOption( "none", "Apply additional processing to input trace?" );
  //  pdef->addValue( "none", VALTYPE_OPTION, "Use input trace as-is" );
  //  pdef->addValue( "envelope", VALTYPE_OPTION, "Compute envelope of input trace" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_trc_math_ens_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_trc_math_ens::VariableStruct* vars = reinterpret_cast<mod_trc_math_ens::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_trc_math_ens_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_trc_math_ens::VariableStruct* vars = reinterpret_cast<mod_trc_math_ens::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;

  if( vars->solver != NULL ) {
    delete vars->solver; vars->solver = NULL;
  }
  if( vars->userVarValues != NULL ) {
    delete [] vars->userVarValues;
    vars->userVarValues = NULL;
  }
  if( vars->varIndexList != NULL ) {
    delete [] vars->varIndexList;
    vars->varIndexList = NULL;
  }
  if( vars->ens != NULL ) {
    delete vars->ens;
    vars->ens = NULL;
  }
  if( vars->gather != NULL ) {
    delete vars->gather;
    vars->gather = NULL;
  }
  if( vars->medianFilt_sortObj != NULL ) {
    delete vars->medianFilt_sortObj;
    vars->medianFilt_sortObj = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_trc_math_ens_( csParamDef* pdef ) {
  params_mod_trc_math_ens_( pdef );
}
extern "C" void _init_mod_trc_math_ens_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_trc_math_ens_( param, env, writer );
}
extern "C" bool _start_exec_mod_trc_math_ens_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_trc_math_ens_( env, writer );
}
extern "C" void _exec_mod_trc_math_ens_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_trc_math_ens_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_trc_math_ens_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_trc_math_ens_( env, writer );
}


//********************************************************************************
//
//
void computeEnsemble( int indexEns, csTraceGather* gatherIn, mod_trc_math_ens::Ens* ens, int numSamples, int numDim2, csTraceGather* gatherOut ) {
  int nTracesOutCurrent = gatherOut->numTraces();

  int numEnsHalf = (int)( ens->numEns()/2 );
  int ensFirst = std::max(0, indexEns - numEnsHalf );
  int ensLast  = std::max( indexEns - numEnsHalf, ens->numEns()-1 );

  int nTracesCurrentEnsemble = ens->numTraces(indexEns);

  int traceIndex1CentreEns = 0;
  //  for( int iensIn = ensFirst; iensIn < indexEns; iensIn++ ) {
  for( int iensIn = 0; iensIn < indexEns; iensIn++ ) {
    traceIndex1CentreEns += ens->numTraces(iensIn);
  }

  // ...for every trace in the ensemble that shall be computed:
  for( int itrc = 0; itrc < nTracesCurrentEnsemble; itrc++ ) {
    // Copy currently processed trace from centre ensemble to output gather       
    int traceIndexIn = traceIndex1CentreEns + itrc;
    //    fprintf(stdout,"ntracesIn: %d, ntracesOut: %d, traceIndex: %d/%d\n", gatherIn->numTraces(), gatherOut->numTraces(), itrc, traceIndexIn );
    gatherIn->copyTraceTo( traceIndexIn, gatherOut );
    // Retrieve pointer to data samples
    float* samplesCurrentTrace = gatherOut->trace( nTracesOutCurrent+itrc )->getTraceSamples();
    // Zero out current trace samples
    for( int isamp = 0; isamp < numSamples; isamp++ ) {
      samplesCurrentTrace[isamp] = 0.0;
    }
    int traceIndex1ProcessedEns = 0;
    int traceCounter = 0;
    // ..for every ensemble that shall be averaged for this trace
    for( int iensIn = ensFirst; iensIn <= ensLast; iensIn++ ) {
      int nTracesProcessedEnsemble = ens->numTraces(iensIn);
      int trcFirst = std::max( 0, itrc - numDim2 );
      int trcLast  = std::min( itrc + numDim2, nTracesProcessedEnsemble-1 );
      // ..for every trace that shall be averaged for this trace
      for( int itrcIn = trcFirst; itrcIn <= trcLast; itrcIn++ ) {
        float* samplesIn = gatherIn->trace( traceIndex1ProcessedEns + itrcIn )->getTraceSamples(); 
        for( int isamp = 0; isamp < numSamples; isamp++ ) {
          samplesCurrentTrace[isamp] += samplesIn[isamp];
        }
        traceCounter += 1;
      } // END: for itrcIn
      traceIndex1ProcessedEns += nTracesProcessedEnsemble;
    } // END: for iensIn
    // Normalization:
    for( int isamp = 0; isamp < numSamples; isamp++ ) {
      samplesCurrentTrace[isamp] /= traceCounter;
    }
  } // END: for itrc
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>
#include <cstring>
#include "csToken.h"
#include "csVector.h"
#include "csEquationSolver.h"
#include "csSort.h"
#include "csFFTTools.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: trc_math
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_trc_math {
  struct VariableStruct {
    float value;
    bool isAdd;
    cseis_geolib::csEquationSolver* solver;
    double* userConstantValues;
    int  numVariables;
    int option;
    bool isGradient;
    float gradient;
    float gradientScalar;
    double dbScalar;
    float* buffer;
    bool isMedianFilt;
    bool isMeanFilt;
    int numSampMedianFilt;
    int numSampMedianFiltALL;
    cseis_geolib::csSort<float>* sortObj;
    float* medianBuffer;
    bool isRMSThreshFilter;
    float rmsThresh_amp;
    int rmsThresh_nsampWin;
    int rmsThresh_nsampTaper;
    int hdrId_rmsThreshold;
    int startSample;
    int endSample;
    cseis_geolib::csFFTTools* fftTool;
    int despike_nsampTaper;
    int minNumZeroSamples;
    bool nanReplace;
    float nanReplacementValue;
    int hdrId_val;
    int hdrValFlag;
  };
  static int const OPTION_NONE     = 0;
  static int const OPTION_EQUATION = 1;
  static int const OPTION_ADD      = 2;
  static int const OPTION_DB       = 4;
  static int const OPTION_FLIP     = 8;
  static int const OPTION_RMS_THRESH = 16;
  static int const OPTION_ENVELOPE = 32;
  static int const OPTION_TAPER_ZEROS = 64;
  static int const OPTION_INST_FREQ  = 128;
  static int const OPTION_INST_PHASE = 256;
  static int const HDR_VAL_NONE = 0;
  static int const HDR_VAL_ADD  = 1;
  static int const HDR_VAL_MULT = 2;
}
using namespace mod_trc_math;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_trc_math_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->value           = 0;
  vars->option          = OPTION_NONE;
  vars->isGradient      = false;
  vars->gradient = 0;
  vars->gradientScalar = 0.0;
  vars->solver          = NULL;
  vars->numVariables    = 0;
  vars->userConstantValues = NULL;
  vars->dbScalar = 10.0;
  vars->buffer = NULL;
  vars->medianBuffer = NULL;
  
  vars->isMedianFilt = false;
  vars->isMeanFilt = false;
  vars->numSampMedianFilt = 0;
  vars->numSampMedianFiltALL = 0;
  vars->sortObj = NULL;
  vars->minNumZeroSamples = 5;
  
  vars->rmsThresh_amp   = 0;
  vars->rmsThresh_nsampWin  = 0;
  vars->rmsThresh_nsampTaper = 0;
  vars->isRMSThreshFilter = false;
  vars->hdrId_rmsThreshold = -1;
  vars->hdrId_val = -1;
  vars->hdrValFlag = mod_trc_math::HDR_VAL_NONE;
  
  vars->startSample = 0;
  vars->endSample   = shdr->numSamples-1;

  vars->fftTool        = NULL;
  vars->despike_nsampTaper = 25;

  vars->nanReplace = false;
  vars->nanReplacementValue = 0;

  if( param->exists("replace_nan") ) {
    vars->nanReplace = true;
    param->getFloat( "replace_nan", &vars->nanReplacementValue );
  }

  if( param->exists("window") ) {
    double startTime;
    double endTime;
    param->getDouble( "window", &startTime, 0 );
    param->getDouble( "window", &endTime, 1 );
    vars->startSample = (int)round( startTime / shdr->sampleInt );
    vars->endSample   = (int)round( endTime / shdr->sampleInt );
    if( vars->startSample < 0 ) vars->startSample = 0;
    if( vars->endSample >= shdr->numSamples ) vars->endSample = shdr->numSamples-1;
  }
  if( param->exists("gradient") ) {
    vars->isGradient = true;
    param->getFloat( "gradient", &vars->gradient );
    vars->gradientScalar = vars->isGradient ? (vars->gradient*shdr->sampleInt/1000.0) : 0.0;
  }
  if( param->exists("median_filt") ) {
    vars->isMedianFilt = true;
    float length;
    param->getFloat( "median_filt", &length );
    vars->numSampMedianFiltALL = (int)round( length / shdr->sampleInt );
    if( vars->numSampMedianFiltALL < 3 ) {
      vars->numSampMedianFiltALL = 3;
    }
    else if( (int)(vars->numSampMedianFiltALL/2)*2 == vars->numSampMedianFiltALL ) {
      vars->numSampMedianFiltALL += 1;
    }
    vars->numSampMedianFilt = ( vars->numSampMedianFiltALL - 1 ) / 2;
    writer->line("Length of median filter = %d samples  (=%f ms/Hz..)", vars->numSampMedianFiltALL, vars->numSampMedianFiltALL * shdr->sampleInt  );
    vars->sortObj = new csSort<float>();
    vars->buffer = new float[shdr->numSamples];
    vars->medianBuffer = new float[shdr->numSamples];
  }
  if( param->exists("mean_filt") ) {
    if( vars->isMedianFilt ) writer->error("Cannot specify both median and mean filter");
    vars->isMeanFilt = true;
    float length;
    param->getFloat( "mean_filt", &length );
    vars->numSampMedianFilt = (int)( 0.5 * length / shdr->sampleInt );
    vars->numSampMedianFiltALL = 2*vars->numSampMedianFilt + 1;
    writer->line("Length of mean filter = %d samples  (%f)", 2*vars->numSampMedianFilt+1, length );
    vars->buffer = new float[shdr->numSamples];
    vars->medianBuffer = new float[shdr->numSamples];
  }
  if( param->exists("add") ) {
    vars->option = OPTION_ADD;
    param->getFloat( "add", &vars->value );
  }
  if( param->exists("flip") ) {
    string text;
    param->getString( "flip", &text );
    if( !text.compare("yes") ) {
      vars->option |= OPTION_FLIP;
      if( vars->buffer == NULL ) vars->buffer = new float[shdr->numSamples];
    }
    else if( text.compare("no") ) {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  if( param->exists("taper_zeros") ) {
    vars->option = OPTION_TAPER_ZEROS;
    float taperLen_ms;
    param->getFloat("taper_zeros",&taperLen_ms);
    vars->despike_nsampTaper = (int)round(taperLen_ms / shdr->sampleInt);
    if( param->getNumValues("taper_zeros") > 1 ) {
      param->getInt("taper_zeros",&vars->minNumZeroSamples,1);
    }
  }

  if( param->exists("rms_thresh") ) {
    vars->isRMSThreshFilter = true;
    string text;
    csFlexNumber number_rms;
    param->getString( "rms_thresh", &text, 0 );
    if( !number_rms.convertToNumber( text ) ) {
      if( !hdef->headerExists(text) ) {
        writer->error("Specified trace header name containing RMS threshold amplitude does not exist: '%s'", text.c_str());
      }
      vars->hdrId_rmsThreshold = hdef->headerIndex(text);
    }
    else {
      vars->rmsThresh_amp = number_rms.floatValue();
    }
    float length;
    param->getFloat( "rms_thresh", &length, 1 );
    vars->rmsThresh_nsampWin = (int)round(length/shdr->sampleInt);
    if( vars->rmsThresh_nsampWin < 1 ) vars->rmsThresh_nsampWin = 1;
    if( param->getNumValues("rms_thresh") > 2 ) {
      param->getFloat( "rms_thresh", &length, 2 );
      vars->rmsThresh_nsampTaper = (int)round(length/shdr->sampleInt);
    }
    if( vars->rmsThresh_nsampTaper < 0 ) vars->rmsThresh_nsampTaper = 0;
    if( vars->rmsThresh_nsampTaper > vars->rmsThresh_nsampWin ) writer->error("RMS threshold filter: Taper length must be smaller than window length");
    if( vars->buffer == NULL ) vars->buffer = new float[shdr->numSamples];
  }
  if( param->exists("db") ) {
    //    if( vars->option != OPTION_NONE ) writer->error("Only one trc_math option can be specified at once.");
    vars->option |= OPTION_DB;
    param->getFloat( "db", &vars->value );
    if( param->getNumValues("db") > 1 ) {
      string text;
      param->getString( "db", &text, 1 );
      if( !text.compare("power") ) {
        vars->dbScalar = 10.0;
      }
      else if( !text.compare("amp") ) {
        vars->dbScalar = 20.0;  
      }
      else {
        writer->error("Unknown option: %s", text.c_str());
      }
    }
  }
  if( param->exists("equation") ) {
    //    if( vars->option != OPTION_NONE ) writer->error("Only one trc_math option can be specified at once.");
    vars->option |= OPTION_EQUATION;
    std::string equationText;
    param->getString( "equation", &equationText );
    vars->solver = new csEquationSolver();
    string variableName = "x";
    if( !vars->solver->prepare( equationText, &variableName, 1 ) ) {
      writer->error("Error occurred: %s", vars->solver->getErrorMessage().c_str() );
    }
    csVector<string> constList;
    vars->solver->prepareUserConstants( &constList );

    vars->numVariables = constList.size();
    vars->userConstantValues = new double[vars->numVariables];
    for( int ivar = 0; ivar < vars->numVariables; ivar++ ) {
      if( constList.at(ivar).compare("x") ) {
        writer->error("Unknown variable name in equation: '%s'. Use variable 'x' to reference sample value.", constList.at(ivar).c_str() );
        env->addError();
      }
      if( edef->isDebug() ) writer->line("Variable #%d: %s", ivar, constList.at(ivar).c_str() );
    }
  }
  if( param->exists("header") ) {
    string text;
    param->getString( "header", &text );
    if( !hdef->headerExists(text) ) writer->error("Trace header not found: %s", text.c_str() );
    vars->hdrId_val = hdef->headerIndex(text);
    param->getString( "header", &text, 1 );
    if( !text.compare("add") ) {
      vars->hdrValFlag = mod_trc_math::HDR_VAL_ADD;
    }
    else if( !text.compare("multiply") ) {
      vars->hdrValFlag = mod_trc_math::HDR_VAL_MULT;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }
  if( param->exists("transform") ) {
    //    if( vars->option != OPTION_NONE ) writer->error("Transform can only be specified by itself. Remove all other parameters from module.");
    string text;
    param->getString( "transform", &text );
    if( !text.compare("inst_freq") ) {
      vars->option = OPTION_INST_FREQ;
    }
    else if( !text.compare("inst_phase") ) {
      vars->option = OPTION_INST_PHASE;
    }
    else if( !text.compare("envelope") ) {
      vars->option = OPTION_ENVELOPE;
    }
    vars->fftTool = new csFFTTools( shdr->numSamples, shdr->sampleInt );
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_trc_math_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);

  float* samplesALL = trace->getTraceSamples();
  float* samples = &samplesALL[vars->startSample];
  int nSamples = vars->endSample - vars->startSample + 1;


  if( vars->solver != NULL ) {
    if( vars->numVariables == 1 ) {
      for( int isamp = 0; isamp < nSamples; isamp++ ) {
        vars->userConstantValues[0] = (double)samples[isamp];
        vars->solver->setUserConstants( vars->userConstantValues, 1 );
        samples[isamp] = (float)vars->solver->solve();
      }
    }
    else if( vars->numVariables > 1 ) {
      for( int isamp = 0; isamp < nSamples; isamp++ ) {
        for( int ivar = 0; ivar < vars->numVariables; ivar++ ) {
          vars->userConstantValues[ivar] = (double)samples[isamp];
        }
        vars->solver->setUserConstants( vars->userConstantValues, vars->numVariables );
        samples[isamp] = (float)vars->solver->solve();
      }
    }
    else {
      for( int isamp = 0; isamp < nSamples; isamp++ ) {
        samples[isamp] = (float)vars->solver->solve();
      }
    }
  }

  if( (vars->option & OPTION_ADD) != 0 ) {
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      samples[isamp] += vars->value;
    }
  }
  
  if( (vars->option & OPTION_DB) != 0 ) {
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      samples[isamp] = vars->dbScalar * log10( samples[isamp] + vars->value );
    }
  }

  if( vars->isGradient ) {
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      samples[isamp] += vars->gradientScalar * (float)isamp;
    }
  }

  if( vars->isMedianFilt ) {
    for( int isamp = vars->numSampMedianFilt; isamp < nSamples-vars->numSampMedianFilt; isamp++ ) {
      int index1 = isamp-vars->numSampMedianFilt;
      memcpy( vars->buffer, &samples[index1], vars->numSampMedianFiltALL * sizeof(float) );
      vars->sortObj->treeSort( vars->buffer, vars->numSampMedianFiltALL );
      vars->medianBuffer[isamp] = vars->buffer[vars->numSampMedianFilt];
    }
    memcpy( &samples[vars->numSampMedianFilt], &vars->medianBuffer[vars->numSampMedianFilt],  (nSamples-vars->numSampMedianFilt) * sizeof(float) );
  }
  else if( vars->isMeanFilt ) {
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      int startSamp = std::max( isamp-vars->numSampMedianFilt, 0 );
      int endSamp   = std::min( isamp+vars->numSampMedianFilt, nSamples-1 );
      int numValues = endSamp - startSamp + 1;
      float sum = 0;
      for( int i = 0; i < numValues; i++ ) {
        sum += samples[i+startSamp];
      }
      vars->buffer[isamp] = sum / (float)numValues;
    }
    memcpy( samples, vars->buffer, nSamples * sizeof(float) );
  }

  if( vars->isRMSThreshFilter ) {
    if( vars->hdrId_rmsThreshold >= 0 ) {
      vars->rmsThresh_amp = trace->getTraceHeader()->floatValue( vars->hdrId_rmsThreshold );
    }
    int half = vars->rmsThresh_nsampWin/2;
    int startSamp = 0;
    int endSamp = nSamples - vars->rmsThresh_nsampWin;
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      vars->buffer[isamp] = 1;
    }
    for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
      float rms = 0.0;
      for( int irms = 0; irms < vars->rmsThresh_nsampWin; irms++ ) {
        rms += pow( samples[isamp+irms], 2 );
      }
      rms = sqrt( rms ) / (float)vars->rmsThresh_nsampWin;
      if( rms < vars->rmsThresh_amp ) vars->buffer[isamp+half] = 0;
    }
    // Edges:
    for( int isamp = 0; isamp < half; isamp++ ) {
      float rms = 0.0;
      int endSamp2 = std::min( nSamples-1, isamp+half );
      for( int irms = 0; irms <= endSamp2; irms++ ) {
        rms += pow( samples[irms], 2 );
      }
      rms = sqrt( rms ) / (float)(endSamp2 + 1);
      if( rms < vars->rmsThresh_amp ) vars->buffer[isamp] = 0;
    }
    
    for( int isamp = nSamples-half; isamp < nSamples; isamp++ ) {
      float rms = 0.0;
      int startSamp2 = std::max( 0, isamp-half );
      for( int irms = startSamp2; irms < nSamples; irms++ ) {
        rms += pow( samples[irms], 2 );
      }
      rms = sqrt( rms ) / (float)( nSamples - startSamp2 + 1 );
      if( rms < vars->rmsThresh_amp ) vars->buffer[isamp] = 0;
    }
    
    if( vars->rmsThresh_nsampTaper > 0 ) {
      // Apply taper to zeroed sections:
      float prevSamp = samples[0];
      int isamp = 1;
      while( isamp < nSamples ) {
        float currentSamp = vars->buffer[isamp];
        if( prevSamp * currentSamp == 0 && currentSamp > 0 ) {
          int samp1 = isamp;
          int samp2 = std::min( nSamples, samp1+vars->rmsThresh_nsampTaper-1 );
          int ns = samp2 - samp1 + 1;
          for( int is = samp1; is <= samp2; is++ ) {
            vars->buffer[is] *= 1.0 - 0.5 * ( 1.0 + cos( M_PI * (float)(is-samp1)/(float)(ns-1) ) );
          }
          isamp += 1;
          currentSamp = 1;
        }
        else {
          isamp += 1;
        }
        prevSamp = currentSamp;
      }
      prevSamp = samples[0];
      isamp = 1;
      while( isamp < nSamples ) {
        float currentSamp = vars->buffer[isamp];
        if( prevSamp * currentSamp == 0 && prevSamp > 0 ) {
          int samp1 = std::max( 0, isamp-vars->rmsThresh_nsampTaper );
          int samp2 = std::min( nSamples, samp1+vars->rmsThresh_nsampTaper-1 );
          int ns = samp2 - samp1 + 1;
          for( int is = samp1; is <= samp2; is++ ) {
            vars->buffer[is] *= 0.5 * ( 1.0 + cos( M_PI * (float)(is-samp1)/(float)(ns-1) ) );
          }
          isamp += 1;
          currentSamp = 0;
        }
        else {
          isamp += 1;
        }
        prevSamp = currentSamp;
      } // END while
      
    }
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      samples[isamp] *= vars->buffer[isamp];
    }
  }
  if( vars->option == OPTION_TAPER_ZEROS ) {

    // Final step: Apply cosine tapers to edges of areas with zero values

    // Apply taper to zeroed sections:
    float* samples = trace->getTraceSamples();
    int isamp = vars->startSample;
    float prevSamp = samples[isamp++];
    int zeroSampleCounter = 0;
    while( isamp <= vars->endSample ) {
      if( prevSamp != 0.0f ) {
        zeroSampleCounter = 0;
      }
      else {
        zeroSampleCounter += 1;
      }
      float currentSamp = samples[isamp];
      if( zeroSampleCounter >= vars->minNumZeroSamples && currentSamp != 0 ) {
        int samp1 = isamp;
        int samp2 = std::min( vars->endSample, samp1+vars->despike_nsampTaper-1 );
        int ns = samp2 - samp1 + 1;
        for( int is = samp1; is <= samp2; is++ ) {
          samples[is] *= 1.0 - 0.5 * ( 1.0 + cos( M_PI * (float)(is-samp1)/(float)(ns-1) ) );
        }
        zeroSampleCounter = 0;
      }
      isamp += 1;
      prevSamp = currentSamp;
    }
    isamp = vars->endSample;
    prevSamp = samples[isamp--];
    zeroSampleCounter = 0;
    while( isamp >= 0 ) {
      if( prevSamp != 0.0f ) {
        zeroSampleCounter = 0;
      }
      else {
        zeroSampleCounter += 1;
      }
      float currentSamp = samples[isamp];
      if( zeroSampleCounter >= vars->minNumZeroSamples && currentSamp != 0 ) {
        int samp1 = std::max( 0, isamp-vars->despike_nsampTaper+1 );
        int samp2 = std::min( vars->endSample, samp1+vars->despike_nsampTaper-1 );
        int ns = samp2 - samp1 + 1;
        for( int is = samp1; is <= samp2; is++ ) {
          samples[is] *= 0.5 * ( 1.0 + cos( M_PI * (float)(is-samp1)/(float)(ns-1) ) );
        }
        zeroSampleCounter = 0;
      }
      isamp -= 1;
      prevSamp = currentSamp;
    } // END while
  }
  
  if( (vars->option & OPTION_FLIP) != 0 ) {
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      vars->buffer[nSamples-isamp-1] = samples[isamp];
    }
    memcpy(samples,vars->buffer,nSamples*sizeof(float));
  }
  

  if( (vars->option & OPTION_ENVELOPE) != 0 ) {
    vars->fftTool->envelope( samples );
  }
  else if( (vars->option & OPTION_INST_FREQ) != 0 ) {
    vars->fftTool->inst_freq( samples );
  }
  else if( (vars->option & OPTION_INST_PHASE) != 0 ) {
    vars->fftTool->inst_phase( samples );
  }

  if( vars->nanReplace ) {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      if( cseis_geolib::is_nan_inf(samples[isamp]) ) {
        samples[isamp] = vars->nanReplacementValue;
      }
    }
  }

  if( vars->hdrValFlag != mod_trc_math::HDR_VAL_NONE ) {
    double hdrValue = trace->getTraceHeader()->doubleValue( vars->hdrId_val );
    if( vars->hdrValFlag == mod_trc_math::HDR_VAL_ADD ) {
      for( int is = vars->startSample; is <= vars->endSample; is++ ) {
        samples[is] += hdrValue;
      }
    }
    else {
      for( int is = vars->startSample; is <= vars->endSample; is++ ) {
        samples[is] *= hdrValue;
      }
    }
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_trc_math_( csParamDef* pdef ) {
  pdef->setModule( "trc_math", "Trace sample computation", "Perform mathematical computations on trace sample values." );

  pdef->addDoc("If more than one user parameter is specified, operations are applied in the order they are defined in this help document, NOT in the order they are supplied in the flow file");

  pdef->addParam( "equation", "Mathematical equation", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Arbitrary mathematical equation to apply to each individual trace sample. Use 'x' to reference original sample value.",
                  "Constants: pi,e. Functions: abs,acos,asin,atan,atan2,ceil,cos,cosh,exp,floor,log,log10,max,min,mod,pow,int,round,sin,sinh,sqrt,tan,tanh,todegrees,toradians,sign");

  pdef->addParam( "add", "Add constant value to sample value", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Constant value to add to trace sample value" );

  pdef->addParam( "gradient", "Add time-variant gradient to trace sample value ", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Gradient. Time-variant gradient is computed as follows: time_variant_gradient_value = gradient * time_in_s" );

  pdef->addParam( "db", "Convert to dB", NUM_VALUES_VARIABLE );
  pdef->addValue( "0", VALTYPE_NUMBER, "Added noise. Set to other than zero to prevent taking logarithm of zero" );
  pdef->addValue( "amp", VALTYPE_OPTION );
  pdef->addOption( "amp", "Apply amplitude dB equation: (20*log10(x+noise))" );
  pdef->addOption( "power", "Apply power dB equation: (10*log10(x+noise))" );

  pdef->addParam( "flip", "Flip trace upside down?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not flip trace" );
  pdef->addOption( "yes", "Flip trace upside down: Reverse the order of the sample values" );

  pdef->addParam( "median_filt", "Apply median filter to input trace", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Total length of median filter in unit of trace", "If the specified filter length translates into an even number of samples, one sample is added" );

  pdef->addParam( "mean_filt", "Apply mean filter to input trace", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Length of mean filter in unit of trace" );

  pdef->addParam( "rms_thresh", "Apply RMS amplitude threshold sliding window filter to input trace", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "RMS amplitude threshold: Amplitudes above this are passed through the filter" );
  pdef->addValue( "", VALTYPE_NUMBER, "Window length [ms]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Cosine taper length to avoid steep edges [ms]" );

  pdef->addParam( "taper_zeros", "Taper areas with zero values", NUM_VALUES_VARIABLE );
  pdef->addValue( "50", VALTYPE_NUMBER, "Taper length [ms]" );
  pdef->addValue( "5", VALTYPE_NUMBER, "Number of consecutive samples that need to be zero for taper to be applied" );

  pdef->addParam( "replace_nan", "Replace NAN and INF values", NUM_VALUES_FIXED, "If this user parameter is specified, NAN and INF values are replaced with the given value" );
  pdef->addValue( "", VALTYPE_NUMBER, "Replacement value" );

  pdef->addParam( "header", "Add or multiply trace header value", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "add", VALTYPE_OPTION );
  pdef->addOption( "add", "Add trace header value" );
  pdef->addOption( "multiply", "Multiply trace header value" );
  
  pdef->addParam( "transform", "Transform data trace", NUM_VALUES_FIXED );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "none", "Do not apply any data transform" );
  pdef->addOption( "inst_freq", "Compute instanteous frequency" );
  pdef->addOption( "inst_phase", "Compute instanteous phase" );
  pdef->addOption( "envelope", "Compute envelope" );

  pdef->addParam( "window", "Application window - Limit application of specified math to specific time window", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Start time/frequency [ms] or [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "End time/frequency [ms] or [Hz]" );  
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_trc_math_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_trc_math::VariableStruct* vars = reinterpret_cast<mod_trc_math::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_trc_math_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_trc_math::VariableStruct* vars = reinterpret_cast<mod_trc_math::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->solver != NULL ) {
    delete vars->solver; vars->solver = NULL;
  }
  if( vars->userConstantValues != NULL ) {
    delete [] vars->userConstantValues;
    vars->userConstantValues = NULL;
  }
  if( vars->buffer != NULL ) {
    delete [] vars->buffer; vars->buffer = NULL;
  }
  if( vars->medianBuffer != NULL ) {
    delete [] vars->medianBuffer; vars->medianBuffer = NULL;
  }
  if( vars->sortObj != NULL ) {
    delete vars->sortObj;
    vars->sortObj = NULL;
  }
  if( vars->fftTool != NULL ) {
    delete vars->fftTool;
    vars->fftTool = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_trc_math_( csParamDef* pdef ) {
  params_mod_trc_math_( pdef );
}
extern "C" void _init_mod_trc_math_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_trc_math_( param, env, writer );
}
extern "C" bool _start_exec_mod_trc_math_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_trc_math_( env, writer );
}
extern "C" void _exec_mod_trc_math_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_trc_math_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_trc_math_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_trc_math_( env, writer );
}


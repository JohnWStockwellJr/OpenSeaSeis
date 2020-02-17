/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cstddef>
#include <cmath>
#include <cstring>

extern "C" {
 #include <fftw3.h>
 #include <limits.h>
}

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: FFT_2D
 *
 * @author Bernie Pinsonnault
 * @date   2011
 */
namespace mod_fft_2d {
  struct VariableStruct {

    bool firstCall;         // First call to exec phase?

    // Control logic in exec phase
    int ntraceLastCall;     // Number of traces in last exec call.
    int nChangeNtr;         // Count the number of times the trace count changes between calls.

    // User information
    int    direction;       // Direction of transform: 'forward' or 'inverse'
    string direction_str;
    int fftDataType;        // Datatype of the transformed data after FORWARD. 
    int fixedNtr;           // Fixed number of traces on output

    int taperType;            // Type of taper tp apply to input 
    int taperLengthInSamples; // Taper length in number of samples (from 0 to 1)

    bool normalize;      // Normalize data: Yes or no?
    double normScalar;   // Normalization scalar to be applied before inverse FFT.
                         // This scalar accounts for zeros that may have been padded to input trace

    // Data details
    int nSamplesInput;      // Number of samples input (4-byte float)
    float sampleRateInput;  // Sample rate of input data (ms or hz).

    float timeSampleRate;   // Sample rate in time domain

    int nSamplesTime;       // Number of samples in time domain *before* FOWARD
    int nSamplesToFFT;      // Number of samples after padding *before* FOWARD
    int nFreqToFFT;         // Number of samples in frequency domain (complex == 2*float)

    float freqSampleRate;   // Sample rate in frequency domain
    float nyquist;          // Nyquist frequency

    int nSamplesToOutput;   // Number of samples output as equivalent # of floats (ex: 1 complex = 2 floats)

    // Related headers
    string KHeader;         // Name of header to store delta K
    int    KID;             // Header index of KHeader
    int hdrID_knorm; // Header index of k_norm header
    float  deltaK;          // wave number increment. 
    string padHeader;       // Name of header to flag padded traces (1=pad, 0=orig)
    int    padID;           // Header index of padHeader

    // FFTW stuff
    unsigned fftw_flag;           // Flag for FFTW plan routine
    float* realBuffer;            // Float buffer 
    fftwf_complex* complexBuffer; // Complex buffer
    fftwf_plan  myPlan;           // FFTW plan (presumably a structure?)

  };
  static const string MY_NAME = "fft_2d";

  static const int FORWARD = FFTW_FORWARD;
  static const int INVERSE = FFTW_BACKWARD; 

  static const int TAPER_NONE     = -1;
  static const int TAPER_COSINE   = 1;
  static const int TAPER_HANNING  = 2;
  static const int TAPER_BLACKMAN = 3;

}
using namespace mod_fft_2d;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_fft_2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  // Initialize
  std::string KHeader_OLD   = "k-wavenumber";
  std::string padHeader_OLD = "fft2d-padtrace";
  vars->KHeader        = "k_wavenumber";
  vars->padHeader      = "fft2d_padtrace";
  vars->ntraceLastCall = 0;
  vars->fixedNtr       = 0;
  vars->direction      = 0;
  vars->direction_str  = "UNKNOWN";
  vars->realBuffer     = NULL;
  vars->complexBuffer  = NULL;
  vars->myPlan         = NULL;
  vars->nChangeNtr     = 0;
  vars->firstCall      = true;
  vars->normalize      = false;
  vars->taperType      = TAPER_NONE;
  vars->taperLengthInSamples = 20 / shdr->sampleInt;
  vars->hdrID_knorm = -1;

  vars->nSamplesInput   = shdr->numSamples;
  vars->sampleRateInput = shdr->sampleInt;

  // User parameters:
  std::string text;

  // Normalize the output.
  if( param->exists("norm") ) {
    param->getString("norm", &text);
    text = toLowerCase( text );
    if( text.compare("yes") == 0 ) {
      vars->normalize = true;
    } else if( text.compare("no") == 0 ) {
      vars->normalize = false;
    }
    else {
      writer->line("ERROR:Unknown 'norm' option: %s", text.c_str());
      env->addError();
    }
  }

  //--------------------------------------------------
  param->getString( "direction", &text ); 
  text = toLowerCase( text );
  // FORWARD
  if( text.compare("forward") == 0 ) {
    vars->direction = FORWARD;
    vars->direction_str = "FORWARD";
  }
  else if( text.compare("inverse") == 0 ) {
    vars->direction = INVERSE;
    vars->direction_str = "INVERSE";
  }
  else {
    writer->line("ERROR: Unknown 'direction' parameter specified: %s", text.c_str() );
    env->addError();
  }

  // Override domain stored in the superheader if it conflicts with the user parameters.
  if( param->exists("override") ) {
    param->getString("override", &text);
    text = toLowerCase( text );
    bool overrideFK = false;
    if( text.compare("fk") == 0 ) {
      shdr->domain = DOMAIN_FK;
      overrideFK = true;
    }
    //      else if( text.compare("fx") == 0 ) {
    //   shdr->domain = DOMAIN_FX;
    //  }
    else if( text.compare("xt") == 0 ) {
      shdr->domain = DOMAIN_XT;
    }
    else if( text.compare("xd") == 0 ) {
      shdr->domain = DOMAIN_XD;
    }
    else if( text.compare("no") == 0 ) {
      // Nothing
    }
    else {
      writer->line("ERROR: Unknown option '%s'", text.c_str() );
      env->addError();
    }
    if( overrideFK ) {
      // bool isFXDomain = (shdr->domain == DOMAIN);
      param->getString("override", &text, 1);
      if( text.compare("amp_phase") == 0 ) {
        //shdr->fftDataType = isFXDomain ? FX_AMP_PHASE : FK_AMP_PHASE;
        shdr->fftDataType = FK_AMP_PHASE;
      }
      else if( text.compare("complex") == 0 ) {
        shdr->fftDataType = FK_COMPLEX;
      }
      else if( text.compare("real_imag") == 0 ) {
        shdr->fftDataType = FK_REAL_IMAG;
      }
      else if( text.compare("amp") == 0 ) {
        shdr->fftDataType = FK_AMP;
      }
      else if( text.compare("psd") == 0 ) {
        shdr->fftDataType = FK_PSD;
      }
      else {
        writer->line("ERROR: Unknown option '%s'", text.c_str() );
        env->addError();
      }
      if( vars->direction == INVERSE ) {
        param->getInt("override", &shdr->numSamplesXT, 2);
        param->getFloat("override", &shdr->sampleIntXT, 3);
      }
    }
  }

  // FORWARD
  if( vars->direction == FORWARD ) {
    writer->line("FORWARD transform.");

    // Taper the input
    if( param->exists("taper_type") ) {
      param->getString("taper_type", &text);
      text = toLowerCase( text );      
      if( text.compare("none") == 0 ) {
        vars->taperType = TAPER_NONE;
        
      }else if( text.compare("cos") == 0 ) {
        vars->taperType = TAPER_COSINE;

      } else if( text.compare("hanning") == 0 ) {
        vars->taperType = TAPER_HANNING;
        vars->taperLengthInSamples = shdr->numSamples/2;

      } else if( text.compare("blackman") == 0 ) {
        vars->taperType = TAPER_BLACKMAN;
        vars->taperLengthInSamples = shdr->numSamples/2;

      } else {
        writer->line("ERROR:Unknown 'taper_type' option: %s", text.c_str());
        env->addError();
      }
    }
    if( param->exists("taper_len") ) {
      if( vars->taperType == TAPER_HANNING ) {
        writer->line("ERROR:Cannot specify 'taper_len' for Hanning taper, taper length is fixed to half the trace length");
        env->addError();
      } else {
        float taperLength;
        param->getFloat("taper_len", &taperLength);
        vars->taperLengthInSamples = (int)round( taperLength / shdr->sampleInt );
      }
    }

    // Optimize the trace length for the FFT.
    bool optimize_length = true;
    if( param->exists("length_opt") ) {
      param->getString("length_opt", &text);
      text = toLowerCase( text );
      if( text.compare("yes") == 0 ) {
        optimize_length = true;
      } else if( text.compare("no")  == 0 ) {
        optimize_length = false;
      } else {
        writer->line("ERROR:Unknown 'optimize_length' option: %s", text.c_str());
        env->addError();
      }
    }

    // Check the input domain and set the output domain.
    if( shdr->domain != DOMAIN_XT && shdr->domain != DOMAIN_XD ) {
      writer->line("ERROR:Input traces are not in XT (or XD) domain. FFT forward transform not possible." );
      env->addError();
    }

    // Optimize trace length for FFTW
    vars->nSamplesTime  = vars->nSamplesInput;    
    int nin = vars->nSamplesTime;
    int nout = 0;
    if ( optimize_length ) {
      if( cseis_geolib::factor_2357( nin, nout ) ){
        writer->line("Optimizing number of real samples to transform to %d (input is %d).", nout, nin );
      }
      else {
        writer->line("ERROR:Failed optimizing FFT length %d.", nin );
        env->addError();              
      }    
    } else {
      nout = nin;
      writer->line("Trace length NOT optimized for FFT: %d.", nin );
    }
    vars->nSamplesToFFT = (int)((nout+1)/2)*2;  // Make sure number of FFT samples is even number
    vars->nFreqToFFT    = vars->nSamplesToFFT/2+1;

    // Output option
    text = "amp_phase";
    vars->nSamplesToOutput = 2*vars->nFreqToFFT; 
    vars->fftDataType      = FK_AMP_PHASE;
    if( param->exists("output") ) { 
      param->getString( "output", &text ); 
      vars->fftDataType = 0;
    }
    text = toLowerCase( text );
    if ( text.compare("complex") == 0 ){
      vars->nSamplesToOutput = 2*vars->nFreqToFFT; 
      vars->fftDataType      = FK_COMPLEX;

    } else if ( text.compare("real_imag") == 0 ){
      vars->nSamplesToOutput = 2*vars->nFreqToFFT; 
      vars->fftDataType      = FK_REAL_IMAG;

    } else if ( text.compare("amp_phase") == 0  ){
      vars->nSamplesToOutput = 2*vars->nFreqToFFT; 
      vars->fftDataType      = FK_AMP_PHASE;
        
    } else if ( text.compare("amp") == 0 ){
      vars->nSamplesToOutput = vars->nFreqToFFT; 
      vars->fftDataType      = FK_AMP;

    } else if ( text.compare("psd") == 0 ){
      vars->nSamplesToOutput = vars->nFreqToFFT; 
      vars->fftDataType      = FK_PSD;

    } else {
      writer->line("ERROR: Unknown 'output' parameter specified: %s", text.c_str() );
      env->addError();
    }
    writer->line("Transformed data will be output in format: %s", text.c_str() );
    if ( vars->fftDataType != FK_AMP_PHASE 
         && vars->fftDataType != FK_REAL_IMAG
         && vars->fftDataType != FK_COMPLEX ){
      writer->warning("Will not be able to inverse transform data if only %s data output.", text.c_str() );
    }

    writer->line("Number of output samples is %d.",vars->nSamplesToOutput );

    // Fixed number of output traces, optionally optimize.
    if( param->exists("fix_ntr") ){
      int fixntr;
      param->getInt( "fix_ntr", &fixntr, 0 );

      bool do_opt = true;
      if ( param->getNumValues("fix_ntr") == 2 ){
        param->getString("fix_ntr", &text, 1 );
        text = toLowerCase( text );
        if ( text.compare("yes") == 0 ){
          do_opt = true;
        } else if ( text.compare("no") == 0 ){
          do_opt = false;
        } else {
          writer->line("ERROR:Unknown value for 'fix_ntr': %s.", text.c_str());
          env->addError();    
        }
      }

      if ( do_opt ){
        int out;
        if( cseis_geolib::factor_2357( fixntr, out ) ){
          writer->line("Fixing optimum number of traces at %d (specified was %d).", out, fixntr);
          fixntr = out;
        } else {
          writer->line("ERROR:Cannot optimize 'fix_ntr': %d.",  fixntr );
          env->addError();              
        }
      } else {
        writer->line("Fixing number of output traces at %d. NOT optimized for FFT.", fixntr);
      }

      vars->fixedNtr = fixntr;
    } else {
      writer->line("Number of traces in output ensemble will be same as number of input traces. Not optimized");
    }

    // At the start, always measure; reset automatically in exec if neeeded.
    vars->fftw_flag = FFTW_MEASURE;

    // Use unit wave numbers. Need spacing to determine exact delta-K.
    vars->deltaK = 1.0;

    // Required headers
    if( !hdef->headerExists( vars->padHeader ) ) {
      if( hdef->headerExists( padHeader_OLD ) ) {
        vars->padHeader = padHeader_OLD;
      }
      else {
        hdef->addHeader( TYPE_INT, vars->padHeader );
      }
    }
    vars->padID = hdef->headerIndex( vars->padHeader );
    if( !hdef->headerExists( vars->KHeader ) ) {
      if( hdef->headerExists( KHeader_OLD ) ) {
        vars->KHeader = KHeader_OLD;
      }
      else {
        hdef->addHeader( TYPE_FLOAT, vars->KHeader );
      }
    }
    vars->KID = hdef->headerIndex( vars->KHeader );    
    if( !hdef->headerExists( "k_norm" ) ) {
      hdef->addHeader( TYPE_FLOAT, "k_norm" );
    }
    vars->hdrID_knorm = hdef->headerIndex( "k_norm" );

    // Evaluate frequency info
    vars->timeSampleRate = vars->sampleRateInput;
    vars->nyquist = 1./(2.*(vars->timeSampleRate/1000.0));
    vars->freqSampleRate = vars->nyquist/((float)vars->nFreqToFFT-1);
    writer->line("Input sample rate = %f -> Nyquist is %f @ delta f = %f.",
              vars->timeSampleRate, vars->nyquist, vars->freqSampleRate );
    
    // Save the time-domain information now.
    shdr->numSamplesXT = vars->nSamplesTime;
    shdr->sampleIntXT  = vars->sampleRateInput;

    // Reset samples & sample rate for output
    shdr->numSamples = vars->nSamplesToOutput;
    shdr->sampleInt  = vars->freqSampleRate;

    shdr->domain      = DOMAIN_FK;
    shdr->fftDataType = vars->fftDataType;

    // INVERSE
  }
  else if( vars->direction == INVERSE ) {
    writer->line("INVERSE transform.");

    if( param->exists("output") ) {writer->warning("Option 'ouput' ignored for INVERSE." ); }
    if( param->exists("fix_ntr") ) {writer->warning("Option 'fix_ntr' ignored for INVERSE." ); }
    if( param->exists("length_opt") ) {writer->warning("Option 'length_opt' ignored for INVERSE." ); }
    if( param->exists("taper_type") ) {writer->warning("Option 'taper_type' ignored for INVERSE." ); }
    if( param->exists("taper_len") ) {writer->warning("Option 'taper_len' ignored for INVERSE." ); }

    // Check the input domain and set the output domain.
    if( shdr->domain != DOMAIN_FK ) {
      writer->line("ERROR:Input traces are not in FK domain. FFT inverse transform not possible" );
      env->addError();
    }

    // Evaluate input data type.
    vars->fftDataType = shdr->fftDataType;
    if ( vars->fftDataType == FK_AMP_PHASE ){
      text = "FK_AMP_PHASE";
    } else if ( vars->fftDataType == FK_REAL_IMAG ){
      text = "FK_REAL_IMAG";
    } else if ( vars->fftDataType == FK_COMPLEX ){
      text = "FK_COMPLEX";
    } else if ( vars->fftDataType == FK_AMP ){
      text = "FK_AMP";
    } else if ( vars->fftDataType == FK_PSD ){
      text = "FK_PSD";
    } else {
      text = "UNKNOWN";
      writer->line("ERROR:Cannot inverse transform data of type %d.", vars->fftDataType );
      env->addError();      
    }

    writer->line("Determined input data is in format: %s", text.c_str() );
    if ( vars->fftDataType != FK_AMP_PHASE 
         && vars->fftDataType != FK_REAL_IMAG
         && vars->fftDataType != FK_COMPLEX ){
      writer->line("ERROR: Cannot inverse transform data if only %s data present.", text.c_str() );
      env->addError();
    }

    vars->nSamplesTime  = shdr->numSamplesXT; 
    // vars->nSamplesToFFT = vars->nSamplesTime;
    // vars->nFreqToFFT    = vars->nSamplesToFFT/2+1;
    vars->nFreqToFFT = vars->nSamplesInput/2;
    vars->nSamplesToFFT = (vars->nFreqToFFT-1)*2;

    // NOTE: This is same as FFT (1D) module but doesn't appear to be used anywhere?
    if( vars->normalize ) vars->normScalar = (double)shdr->numSamplesXT/(double)shdr->numSamples;

    // At the start, always measure; reset automatically in exec if neeeded.
    vars->fftw_flag = FFTW_MEASURE; 

    // Required headers
    vars->padID = 0;
    if( hdef->headerExists( vars->padHeader ) ){
      vars->padID = hdef->headerIndex( vars->padHeader );
      hdef->deleteHeader( vars->padHeader );
    } else {
      writer->warning("Missing %s header. Cannot determine if traces were padded. All traces will be output.", vars->padHeader.c_str() );
    }

    // FORWARD headers no longer needed
    if( hdef->headerExists( vars->KHeader ) ) hdef->deleteHeader( vars->KHeader );    

    // Number of samples on output & reset sample rate for time.
    vars->nSamplesToOutput = shdr->numSamplesXT;
    shdr->numSamples   = vars->nSamplesToOutput;
    shdr->sampleInt    = shdr->sampleIntXT;
    shdr->numSamplesXT = 0;
    shdr->sampleIntXT  = 0;

    shdr->domain      = DOMAIN_XT;    
    shdr->fftDataType = FK_NONE;

  }
}

//*******************************************************************************************
// Exec phase
//*******************************************************************************************

void exec_mod_fft_2d_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader    const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;


  int numTracesIn = traceGather->numTraces();

  // If the number of traces out is fixed, we should never recalculate 
  // or reallocate, except on the first call. 
  bool recalc = false;
  if ( vars->fixedNtr > 0 ){

    //Make sure the number of traces doesn't exceed the fixed value.
    if ( traceGather->numTraces() >  vars->fixedNtr ) {
      writer->error("Number of traces input %d exceeds user specified limit % d", numTracesIn, vars->fixedNtr );
    }

  } else {

    // If the number of traces changes, recalculate plan & realloc buffers.
    // Always true on the first exec call.
    if ( traceGather->numTraces() != vars->ntraceLastCall ){
      recalc = true;
      vars->nChangeNtr++;

      // If the size changes too often, switch to ESTIMATE as the plan calculation
      // takes longer for MEASURE so we don't want to do it a lot. 
      // Note this may result in a slower performance of the actual transform. 
      if ( vars->nChangeNtr == 5 ) vars->fftw_flag = FFTW_ESTIMATE;
    }
    vars->ntraceLastCall = traceGather->numTraces();
  }

  int isamp, itrc, index;
  float           *fftR;
  fftwf_complex   *fftC;
  fftwf_plan  fftPlan2D;

  // FORWARD  
  if ( vars->direction == FORWARD ){

    // Apply taper to input traces (time direction only!)
    if( vars->taperType == TAPER_COSINE || vars->taperType == TAPER_HANNING ) {
      int nSteps = vars->taperLengthInSamples;

      for(itrc=0; itrc<numTracesIn; itrc++ ) {
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        for( int i = 0; i < nSteps; i++ ) {
          float scalar = (float)cos( M_PI_2 * (float)(nSteps-i)/(float)nSteps );
          samples[i] *= scalar;
        }
        for( int i = vars->nSamplesInput-nSteps; i < vars->nSamplesInput; i++ ) {
          float scalar = (float)cos( M_PI_2 * (float)(nSteps-vars->nSamplesInput+i+1)/(float)nSteps );
          samples[i] *= scalar;
        }
      }

    } else if( vars->taperType == TAPER_BLACKMAN ) {
      float alpha = 0.16;
      float a0 = 0.5 * (1.0 - alpha);
      float a1 = 0.5;
      float a2 = 0.5 * alpha;
      for(itrc=0; itrc<numTracesIn; itrc++ ) {
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        for( int i = 0; i < vars->nSamplesInput; i++ ) {
          float piFactor = (2.0 * M_PI) * (float)i / (float)(vars->nSamplesInput - 1);
          float weight = a0 - a1*cos( piFactor ) + a2*cos( 2 * piFactor );
          samples[i] *= weight;
        }
      }
    }

    int nx  = numTracesIn;
    int nxp = nx;
    if ( vars->fixedNtr > 0 ) nxp = vars->fixedNtr;

    int nt    = vars->nSamplesTime;
    int ntp   = vars->nSamplesToFFT;
    int nfreq = vars->nFreqToFFT;

    // Allocate buffers & calculate the plan
    if ( recalc || vars->firstCall ){
      if ( vars->firstCall ) vars->firstCall = false;

      // Deallocate previous
      if( vars->realBuffer != NULL ) {
        fftwf_free(vars->realBuffer); vars->realBuffer = NULL;
      }
      if( vars->complexBuffer != NULL ) {
        fftwf_free(vars->complexBuffer); vars->complexBuffer = NULL;
      }   
      if( vars->myPlan != NULL ) {
        fftwf_destroy_plan (vars->myPlan); vars->myPlan = NULL;
      }

      // Buffer allocations
      fftR = (float*)fftwf_malloc(sizeof(float) * ntp*nxp);
      fftC = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) *nxp*nfreq);
      fftPlan2D = fftwf_plan_dft_r2c_2d(nxp, ntp, fftR, fftC, vars->fftw_flag);

      vars->realBuffer    = fftR;
      vars->complexBuffer = fftC;
      vars->myPlan        = fftPlan2D;
    }
    fftR = vars->realBuffer;
    fftC = vars->complexBuffer;
    fftPlan2D = vars->myPlan;

    // Copy input to real array
    memset (fftR, 0, ntp*nxp*sizeof(float));
    //    for(itrc=0; itrc<nx; itrc++ ) {
    for(itrc=0; itrc<numTracesIn; itrc++ ) {
      float *samples = traceGather->trace(itrc)->getTraceSamples();
      //BAP: Why is this mulitplied by -1?
      if (itrc%2 !=0)
        for(isamp=0; isamp<nt; isamp++) 
          samples[isamp] *= -1;      
      memcpy (&fftR[itrc*ntp], samples, nt*sizeof(float));
    }
  
    // FFT
    memset (fftC, 0, sizeof(fftwf_complex)*nxp*nfreq);
    fftwf_execute(fftPlan2D);

    // Create padded traces, copy headers from last input trace to pad traces.
    if ( vars->fixedNtr-numTracesIn > 0 ){
      traceGather->createTraces( numTracesIn, vars->fixedNtr-numTracesIn, hdef, shdr->numSamples );
      csTraceHeader* trcHdrOrig = traceGather->trace(numTracesIn-1)->getTraceHeader();
      for (int itrc = numTracesIn; itrc < nxp; itrc++){
        traceGather->trace(itrc)->getTraceHeader()->copyFrom( trcHdrOrig );
      }      
    }

    // Set pad header
    int KID   = vars->KID;
    int padID = vars->padID;
    int ipad  = 0; 
    float knormValue = 2.0 / (float)(numTracesIn-1);
    for(itrc=0; itrc<nxp; itrc++ ) {
      if ( itrc == numTracesIn  ) ipad=1;
      traceGather->trace(itrc)->getTraceHeader()->setIntValue( padID, ipad );
    }

    // Set K number header
    int n2_1 = nxp/2+1;
    float myK;
    for(itrc=0; itrc<n2_1; itrc++ ) {
      myK = -1.0 * (float)(n2_1-itrc-1) * vars->deltaK; 
      csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
      trcHdr->setIntValue( KID, (int)myK );
      trcHdr->setFloatValue( vars->hdrID_knorm, myK*knormValue );
    }
    for(itrc=n2_1; itrc<nxp; itrc++ ) {
      myK = (float)(itrc-n2_1+1) * vars->deltaK; 
      csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
      trcHdr->setIntValue( KID, (int)myK );
      trcHdr->setFloatValue( vars->hdrID_knorm, myK*knormValue );
    }

    double normSpectralDensity = 1.0;
    int output = vars->fftDataType;

    // Copy COMPLEX output to output traces
    if ( output == FK_COMPLEX ){
      int nOut = nfreq*2;
      float* out = (float*)(fftC);
      for(itrc=0; itrc<nxp; itrc++ ){
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        memcpy (samples, &out[itrc*nOut], nOut*sizeof(float));
      }

      // BAP: Is this the correct calculation of normSpectralDensity?
      if( vars->normalize ) {
        normSpectralDensity = 1.0 / sqrt( (double)vars->nSamplesInput * (double)numTracesIn );
        for(itrc=0; itrc<nxp; itrc++ ){
          float *samples = traceGather->trace(itrc)->getTraceSamples();
          for (isamp=0; isamp<nfreq; isamp++) {
            samples[isamp]         *= normSpectralDensity;
            samples[isamp + nfreq] *= normSpectralDensity;
          }        
        }      
      }

      // Copy REAL/IMAG values to first/second halves of output traces
    } else if ( output == FK_REAL_IMAG ){
      
      // BAP: Is this the correct calculation of normSpectralDensity?
      if( vars->normalize ) {
        normSpectralDensity = 1.0 / sqrt( (double)vars->nSamplesInput * (double)numTracesIn );
      }

      for(itrc=0; itrc<nxp; itrc++ ){
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        fftwf_complex *bufC = fftC + itrc*nfreq;
        for (isamp=0; isamp<nfreq; isamp++) {
          samples[isamp]         = (float)(bufC[isamp][0]*normSpectralDensity);
          samples[isamp + nfreq] = (float)(bufC[isamp][1]*normSpectralDensity);
        }        
      }      
          
      // Calculate PSD and store in output traces
    } else if ( output == FK_PSD ){
      
      // BAP: Is this the correct calculation of normSpectralDensity?
      normSpectralDensity = 1.0;
      if( vars->normalize ) {
        int numSamplesNonZero = 0;
        for(itrc=0; itrc<numTracesIn; itrc++ ){
          float *samples = traceGather->trace(itrc)->getTraceSamples();
          for( int i = 0; i < nt; i++ ) {
            if( samples[i] != 0.0 ) numSamplesNonZero += 1;
          }
        }
        if( numSamplesNonZero == 0 ) numSamplesNonZero = 1;
        normSpectralDensity = 1.0 / (numSamplesNonZero * (1000.0/vars->sampleRateInput));
      }

      for(itrc=0; itrc<nxp; itrc++ ){
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        fftwf_complex *bufC = fftC + itrc*nfreq;
        for (isamp=0; isamp<nfreq; isamp++) {
          samples[isamp] = (bufC[isamp][0]*bufC[isamp][0] 
                            + bufC[isamp][1]*bufC[isamp][1])* normSpectralDensity;
        }
      }

      // Calculate amplitude and store in output traces
    } else if ( output == FK_AMP ){
      float amplitude;
      for(itrc=0; itrc<nxp; itrc++ ){
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        fftwf_complex *bufC = fftC + itrc*nfreq;
        for (isamp=0; isamp<nfreq; isamp++) {
          amplitude = sqrt (bufC[isamp][0]*bufC[isamp][0] + 
                            bufC[isamp][1]*bufC[isamp][1]);
          samples[isamp] = amplitude;
        }
      }

      // BAP: Is this the correct calculation of normSpectralDensity?
      if( vars->normalize ) {
        normSpectralDensity = 1.0 / sqrt( (double)vars->nSamplesInput * (double)numTracesIn );
        int numSamplesNormalise = shdr->numSamples;
        for(itrc=0; itrc<nxp; itrc++ ){
          float *samples = traceGather->trace(itrc)->getTraceSamples();
          for( int is = 0; is < numSamplesNormalise; is++ ) {
            samples[is] *= normSpectralDensity;
          }
        }
      }
          
      // Calculate amplitude & phase and store in output traces
    } else if ( output == FK_AMP_PHASE ){
      float amplitude, phase;
      for(itrc=0; itrc<nxp; itrc++ ){
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        fftwf_complex *bufC = fftC + itrc*nfreq;
        for (isamp=0; isamp<nfreq; isamp++) {
          amplitude = sqrt (bufC[isamp][0]*bufC[isamp][0] + 
                            bufC[isamp][1]*bufC[isamp][1]);
          phase     = atan2(bufC[isamp][1], bufC[isamp][0]);
          samples[isamp]         = amplitude;
          samples[isamp + nfreq] = phase;
        }
      }

      // Normalize amp only (?)
      // BAP: Is this the correct calculation of normSpectralDensity?
      if( vars->normalize ) {
        normSpectralDensity = 1.0 / sqrt( (double)vars->nSamplesInput * (double)numTracesIn );
        int numSamplesNormalise = vars->nSamplesToFFT/2;
        for(itrc=0; itrc<nxp; itrc++ ){
          float *samples = traceGather->trace(itrc)->getTraceSamples();
          for( int is = 0; is < numSamplesNormalise; is++ ) {
            samples[is] *= normSpectralDensity;
          }
        }
      }
     
    }
  }  

  // INVERSE
  if ( vars->direction == INVERSE ){
    int nx    = numTracesIn;
    int nxp   = nx;

    int nt    = vars->nSamplesTime;
    int ntp   = vars->nSamplesToFFT;
    int nfreq = vars->nFreqToFFT;

    float scale = 1./(nxp*ntp);

    // Allocate buffers
    if ( recalc ){
      // Deallocate previous
      if( vars->realBuffer != NULL ) {
        fftwf_free(vars->realBuffer); vars->realBuffer = NULL;
      }
      if( vars->complexBuffer != NULL ) {
        fftwf_free(vars->complexBuffer); vars->complexBuffer = NULL;
      }   
      if( vars->myPlan != NULL ) {
        fftwf_destroy_plan (vars->myPlan); vars->myPlan = NULL;
      }
    
      // Buffer allocations
      fftR = (float*)fftwf_malloc(sizeof(float) * ntp*nxp);
      fftC = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) *nxp*nfreq);    
      fftPlan2D = fftwf_plan_dft_c2r_2d(nxp, ntp, fftC, fftR, vars->fftw_flag);

      vars->realBuffer    = fftR;
      vars->complexBuffer = fftC;
      vars->myPlan        = fftPlan2D;      
    }
    fftR = vars->realBuffer;
    fftC = vars->complexBuffer;
    fftPlan2D = vars->myPlan;

    // Copy input complex data to complex buffer
    int output = vars->fftDataType;
    memset (fftC, 0, sizeof(fftwf_complex)*nxp*nfreq);
    if ( output == FK_COMPLEX ){
      int nOut = nfreq*2;
      float* out = (float*)(fftC);
      for(itrc=0; itrc<nx; itrc++ ){
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        memcpy (&out[itrc*nOut], samples, nOut*sizeof(float));
      }

      // 
    } else if ( output == FK_REAL_IMAG ){
      for(itrc=0; itrc<nx; itrc++ ){
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        fftwf_complex *bufC = fftC + itrc*nfreq;
        for (isamp=0; isamp<nfreq; isamp++) {
          bufC[isamp][0] = samples[isamp];
          bufC[isamp][1] = samples[isamp + nfreq];
        }        
      }
      // Convert input amp & phase data to complex numbers and store in buffer.
    } else if ( output == FK_AMP_PHASE ){
      float amplitude, phase;
      for(itrc=0; itrc<nx; itrc++ ){
        float *samples = traceGather->trace(itrc)->getTraceSamples();
        fftwf_complex *bufC = fftC + itrc*nfreq;
        for (isamp=0; isamp<nfreq; isamp++) {
          amplitude = samples[isamp];
          phase     = samples[isamp + nfreq];
          bufC[isamp][0] = amplitude * cos(phase);
          bufC[isamp][1] = amplitude * sin(phase);
        }        
      }

    }

    // FFT
    memset (fftR, 0, ntp*nxp*sizeof(float));
    fftwf_execute(fftPlan2D);

    // Copy result back into seismic trace (discard padded traces if exist).
    int padID = vars->padID;
    for(itrc=0; itrc<nx; itrc++ ) {

      // Skip padded traces
      if ( traceGather->trace(itrc)->getTraceHeader()->intValue( padID ) == 1 ) break;

      float *samples = traceGather->trace(itrc)->getTraceSamples(); 
      for (isamp=0; isamp<nt; isamp++) {
        index = itrc*ntp + isamp;
        //BAP: Why is this mulitplied by -1?
        if (itrc%2==0)
          samples[isamp] = fftR[index]*scale;
        else 
          samples[isamp] = -fftR[index]*scale;
      }
    }

    // Remove padded traces
    traceGather->freeTraces( itrc, nx-itrc );
  }

}

//*************************************************************************************************
// Parameter definition
//*************************************************************************************************
void params_mod_fft_2d_( csParamDef* pdef ) {
  pdef->setModule( "FFT_2D", "2D FFT", "2D Fast Fourier Transform to/from XT - FK domains");

  pdef->addParam( "direction", "Direction of transform.", NUM_VALUES_FIXED );
  pdef->addValue( "forward", VALTYPE_OPTION );
  pdef->addOption( "forward", "Forward transform from x-t to k-w. Input must be in XT domain" );
  pdef->addOption( "inverse", "Inverse transform from k-w to x-t. Input must be in FK domain" );

  pdef->addParam("output", "Type of output after FORWARD transform.", NUM_VALUES_FIXED );
  pdef->addValue("amp_phase", VALTYPE_OPTION );
  pdef->addOption("amp", "output amplitude only");
  pdef->addOption("amp_phase", "output both amplitude and phase");
  pdef->addOption("psd", "Output PSD spectrum" );
  pdef->addOption("real_imag", "Output complex data such that the real part of the complex data is stored contiguously in the first half of the trace and the imaginary part of the complex data is stored contiguously in the second half of the trace");
  pdef->addOption("complex", "Output complex data such that the real and imaginary parts of the first complex element are stored, respectively, in the first and second real float words of the trace; the real and imanginary parts of the second complex element are stored, respecitvely, in the third and fourth real float words of the trace, etc");

  pdef->addParam( "fix_ntr", "For FORWARD transform, fix the number of traces output traces at the specified value.", NUM_VALUES_VARIABLE,
                  "Ensembles with less traces will be padded. This ensures consistent number of wave numbers after trasform. If not specified, the number of traces in each ensemble will be used 'as is' for the transform. This can be inefficent and may result in an inconsistent number of wave numbers per ensemble after the transform.");
  pdef->addValue( "", VALTYPE_NUMBER, "Number of traces to fix the size of the output ensemble. The job will fail if an input ensemble contains more traces than this value" );
  pdef->addValue("yes", VALTYPE_OPTION, "Optionally reset 'fix_ntr' to be optimum for FFT performace. This may slightly increase the supplied value" );
  pdef->addOption("yes", "Reset 'fix_ntr' to improve performance");
  pdef->addOption("no", "Use 'fix_ntr' as specified by the user");

  pdef->addParam( "length_opt", "For FORWARD transform, reset the input traces length to be optimized for the FFT.", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Optimize trace length" );
  pdef->addOption( "no", "Do not optimize trace length, use the input length" );

  pdef->addParam( "taper_type", "For FORWARD transform, taper type of taper to apply to input trace in the time direction. THIS FEATURE IS UNTESTED!", NUM_VALUES_FIXED );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "none", "Do not apply any taper to input trace" );
  pdef->addOption( "cos", "Apply cosine taper to input trace" );
  pdef->addOption( "hanning", "Apply 'Hanning' cosine taper to input trace. Taper length is 1/2 trace" );
  pdef->addOption( "blackman", "Apply 'Blackman' taper (alpha=0.16) to input trace" );

  pdef->addParam( "taper_len", "For FORWARD transform, taper length [ms]. THIS FEATURE IS UNTESTED!", NUM_VALUES_FIXED );
  pdef->addValue( "20", VALTYPE_NUMBER, "Taper length [ms]" );

  pdef->addParam( "norm", "Normalize output. THIS FEATURE IS UNTESTED!", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Normalize output values", "Example: Using the same input data but with different amount of added zeros, the amplitude spectrum will look exactly the same for the same output frequency" );
  pdef->addOption( "no", "Do not normalize output values" );

  pdef->addParam( "override", "Override domain specified in superheader.", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Acknowledge domain found in super header" );
  pdef->addOption( "fk", "Override domain found in super header. Set to FK." );
  pdef->addOption( "xt", "Override domain found in super header. Set to XT (time)." );
  pdef->addOption( "xd", "Override domain found in super header. Set to XD (depth)." );
  pdef->addValue( "amp_phase", VALTYPE_OPTION, "Override data type" );
  pdef->addOption( "amp_phase", "Set FFT data type to amplitude/phase spectrum" );
  pdef->addOption( "amp", "Set FFT data type to amplitude spectrum" );
  pdef->addOption( "psd", "Set FFT data type to psd spectrum" );
  pdef->addOption( "complex", "Set FFT data type to complex spectrum" );
  pdef->addOption( "real_imag", "Set FFT data type to real/imaginary spectrum" );
  pdef->addValue( "", VALTYPE_NUMBER, "Output number of samples (inverse transform only)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Output sample interval [ms] (inverse transform only)" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_fft_2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_fft_2d::VariableStruct* vars = reinterpret_cast<mod_fft_2d::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_fft_2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_fft_2d::VariableStruct* vars = reinterpret_cast<mod_fft_2d::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->realBuffer != NULL ) {
    fftwf_free(vars->realBuffer); vars->realBuffer = NULL;
  }
  if( vars->complexBuffer != NULL ) {
    fftwf_free(vars->complexBuffer); vars->complexBuffer = NULL;
  }   
  if( vars->myPlan != NULL ) {
    fftwf_destroy_plan (vars->myPlan); vars->myPlan = NULL;
  }
  delete vars; 
  vars = NULL;
}

extern "C" void _params_mod_fft_2d_( csParamDef* pdef ) {
  params_mod_fft_2d_( pdef );
}
extern "C" void _init_mod_fft_2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_fft_2d_( param, env, writer );
}
extern "C" bool _start_exec_mod_fft_2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_fft_2d_( env, writer );
}
extern "C" void _exec_mod_fft_2d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_fft_2d_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_fft_2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_fft_2d_( env, writer );
}

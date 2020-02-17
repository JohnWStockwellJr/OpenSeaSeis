/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csFFTDesignature.h"
#include "csASCIIFileReader.h"
#include "csFileUtils.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: DESIGNATURE
 *
 * @author Bjorn Olofsson
 * @date   2011
 */
namespace mod_designature {
  struct VariableStruct {
    csFFTDesignature** fftDesigList;
    int asciiFormat;
    std::string filenameOut;
    bool isFileOutWavelet;
    float filterTimeShift_s;

    int hdrId_input;
    std::string hdrName_input; 
    int* hdrValueList;
    int numInputWavelets;
    int currentIndexHdrValue;
  };
  static int const UNIT_MS = 3;
  static int const UNIT_S  = 4;
}
using mod_designature::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_designature_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->asciiFormat     = -1;
  vars->fftDesigList    = NULL;
  vars->filenameOut     = "";
  vars->isFileOutWavelet = false;
  vars->hdrId_input      = -1;
  vars->hdrValueList     = NULL;
  vars->numInputWavelets = 0;
  vars->currentIndexHdrValue = 0;
  vars->hdrName_input = "";
  vars->filterTimeShift_s = 0;

//---------------------------------------------
//

  std::string filename;
  std::string text;
  int unitTime_ascii = mod_designature::UNIT_MS;

  bool overrideSampleInt = false;
  if( param->exists("override_sample_int") ) {
    param->getString("override_sample_int", &text);
    if( !text.compare("no") ) {
      overrideSampleInt = false;
    }
    else if( !text.compare("yes") ) {
      overrideSampleInt = true;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }

  param->getString("input_wavelet", &filename);
  if( param->getNumValues("input_wavelet") > 1 ) {
    param->getString("input_wavelet", &text, 1);
    if( !text.compare("ms") ) {
      unitTime_ascii = mod_designature::UNIT_MS;
    }
    else if( !text.compare("s") ) {
      unitTime_ascii = mod_designature::UNIT_S;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }
  param->getString("format",&text,0);
  if( !text.compare("columns") ) {
    vars->asciiFormat = cseis_io::csASCIIFileReader::FORMAT_COLUMNS;
  }
  else if( !text.compare("signature") ) {
    vars->asciiFormat = cseis_io::csASCIIFileReader::FORMAT_NUCLEUS_SIGNATURE;
  }
  else {
    writer->error("Unknown option: %s", text.c_str() );
  }

  int colIndexTime  = 0;
  int colIndexValue = 1;
  int colIndexTrace = -1;
  if( vars->asciiFormat == cseis_io::csASCIIFileReader::FORMAT_COLUMNS ) {
    if( param->exists("format") && param->getNumValues("format") > 1 ) {
      param->getInt( "format", &colIndexTime, 1 );
      colIndexTime  -= 1;
      if( param->getNumValues("format") > 2 ) {
        param->getInt( "format", &colIndexValue, 2 );
        colIndexValue -= 1;
        if( param->getNumValues("format") > 3 ) {
          param->getInt( "format", &colIndexTrace, 3 );
          colIndexTrace -= 1;
        }
      }
    }
  }

  if( param->exists("input_hdr") ) {
    param->getString("input_hdr", &vars->hdrName_input);
    if( colIndexTrace < 0 ) {
      writer->line("Input trace header '%s' specified in user parameter 'input_hdr' but identifying column in input wavelet file is undefined (=%d)", vars->hdrName_input.c_str(), colIndexTrace+1);
      writer->line("Specify column number of input wavelet file (user parameter 'input_wavelet') which identifies the wavelet to be used for which input trace");
      writer->line("Identification is done via the specified trace header '%s'", vars->hdrName_input.c_str());
      env->addError();
    }
    if( !hdef->headerExists(vars->hdrName_input) ) {
      writer->error("Trace header does not exist: '%s'", vars->hdrName_input.c_str() );
    }
    vars->hdrId_input   = hdef->headerIndex(vars->hdrName_input);
    int type = hdef->headerType(vars->hdrName_input);
    if( type != cseis_geolib::TYPE_INT ) {
      writer->error("Trace header '%s': Only integer type supported. Use integer trace header as wavelet identifier", vars->hdrName_input.c_str() );
    }
  }


  float percWhiteNoise = 0.01f;
  if( param->exists("white_noise") ) {
    param->getFloat("white_noise",&percWhiteNoise);
  }

  int filterType = csFFTDesignature::AMP_PHASE;
  if( param->exists("option") ) {
    param->getString("option", &text );
    if( !text.compare("phase_only") ) {
      filterType = csFFTDesignature::PHASE_ONLY;
    }
    else if( !text.compare("amp_only") ) {
      filterType = csFFTDesignature::AMP_ONLY;
    }
    else if( !text.compare("amp_phase") ) {
      filterType = csFFTDesignature::AMP_PHASE;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  //----------------------------------------------------------------------
  // Read in first input signature wavelet from external ASCII file
  // ...more wavelets are read in later, if applicable
  //
  cseis_io::ASCIIParam asciiParam;
  cseis_io::csASCIIFileReader asciiFileReader( filename, vars->asciiFormat );
  try {
    bool success = asciiFileReader.initialize( &asciiParam, colIndexTime, colIndexValue, colIndexTrace );
    if( !success ) writer->error("Unknown error occurred during initialization of signature input file. Incorrect or unsupported format?");
    success = asciiFileReader.readNextTrace( &asciiParam );
    if( !success ) writer->error("Unknown error occurred when reading in samples from signature input file. Incorrect or unsupported format?");
  }
  catch( csException& e ) {
    writer->error("Error occurred when initializing input ASCII file: %s", e.getMessage() );
  }

  //
  //----------------------------------------------------------------------

  float sampleInt_ms = asciiParam.sampleInt;
  float timeZero_ms  = -(float)asciiParam.timeFirstSamp;
  // ^^ Minus sign in case time of first sample is negative...or so. See Nucleus format
  // ..this is required because the csDesignature class simply assumes time of first sample = 0

  if( unitTime_ascii == mod_designature::UNIT_S ) {
    sampleInt_ms *= 1000.0;
    timeZero_ms  *= 1000.0;
  }
  if( sampleInt_ms != shdr->sampleInt ) {
    if( !overrideSampleInt ) {
      writer->error("ASCII input file has different sample interval (=%f ms) than input data (=%f ms). Unsupported case.", sampleInt_ms, shdr->sampleInt);
    }
    else {
      writer->warning("ASCII input file has different sample interval (=%f ms) than input data (=%f ms). Ignored.", sampleInt_ms, shdr->sampleInt);
    }
  }
  if( asciiParam.numSamples() > shdr->numSamples ) {
    writer->error("ASCII input file has more samples (=%d) than input data (=%d). Unsupported case.", asciiParam.numSamples(), shdr->numSamples);
  }
  if( asciiParam.numSamples() <= 0 ) {
    writer->error("Could not read in any sample values from input file. Unsupported or incorrect file format?");
  }

  cseis_geolib::csVector<float*> bufferList;
  cseis_geolib::csVector<float> hdrValueList;
  bool readAnotherWavelet = false;
  do {
    float* buffer = new float[shdr->numSamples];
    for( int isamp = 0; isamp < asciiParam.numSamples(); isamp++ ) {
      buffer[isamp] = asciiParam.sample(isamp);
    }
    // Pad signature file with zeros...
    for( int isamp = asciiParam.numSamples(); isamp < shdr->numSamples; isamp++ ) {
      buffer[isamp] = 0.0;
    }
    bufferList.insertEnd( buffer );
    hdrValueList.insertEnd( asciiParam.traceNumber );
    if( edef->isDebug() ) {
      //    writer->line("---------------ASCII file dump-------------------");
      for( int isamp = 0; isamp < asciiParam.numSamples(); isamp++ ) {
        //      writer->line("%d %f", isamp, asciiParam.sample(isamp) );
        fprintf( stdout, "%d %f %d\n", isamp, asciiParam.sample(isamp), asciiParam.traceNumber );
      }
      if( hdrValueList.size() > 1 ) fprintf( stdout, "\n" );
    }
    if( vars->hdrId_input != 0 ) {
      readAnotherWavelet = asciiFileReader.readNextTrace( &asciiParam );
    }
  } while( readAnotherWavelet );

  vars->numInputWavelets = bufferList.size();
  //
  //----------------------------------------------------------------------
  

  //----------------------------------------------------------------------
  // Read in output wavelet from external ASCII file
  //
  std::string filenameOut;
  int unitTimeOut_ascii = unitTime_ascii;
  float* bufferOut = NULL;
  if( param->exists("output_wavelet") ) {
    param->getString("output_wavelet", &filenameOut);
    if( param->getNumValues("output_wavelet") > 1 ) {
      param->getString("output_wavelet", &text, 1);
      if( !text.compare("ms") ) {
        unitTimeOut_ascii = mod_designature::UNIT_MS;
      }
      else if( !text.compare("s") ) {
        unitTimeOut_ascii = mod_designature::UNIT_S;
      }
      else {
        writer->error("Unknown option: '%s'", text.c_str());
      }
    }

    cseis_io::ASCIIParam asciiParamOut;
    try {
      cseis_io::csASCIIFileReader asciiFileReader( filenameOut, vars->asciiFormat );
      bool success = asciiFileReader.initialize( &asciiParamOut, colIndexTime, colIndexValue, colIndexTrace );
      if( !success ) writer->error("Unknown error occurred during initialization of output wavelet file. Incorrect or unsupported format?");
      success = asciiFileReader.readNextTrace( &asciiParamOut );
      if( !success ) writer->error("Unknown error occurred when reading in samples from output wavelet file. Incorrect or unsupported format?");
    }
    catch( csException& e ) {
      writer->error("Error occurred when initializing input ASCII file: %s", e.getMessage() );
    }
    if( asciiParamOut.numSamples() != asciiParam.numSamples() ) {
      writer->error("Output wavelet contains %d number of samples, not matching the number of samples in input wavelet = %d.\nThis is not supported",
                 asciiParamOut.numSamples() != asciiParam.numSamples() );
    }
    bufferOut = new float[shdr->numSamples];
    for( int isamp = 0; isamp < asciiParamOut.numSamples(); isamp++ ) {
      bufferOut[isamp] = asciiParamOut.sample(isamp);
    }
    // Pad signature file with zeros...
    for( int isamp = asciiParamOut.numSamples(); isamp < shdr->numSamples; isamp++ ) {
      bufferOut[isamp] = 0.0;
    }
  }
  if( unitTimeOut_ascii != unitTime_ascii ) {
    writer->error("Different time units specified for input & output wavelet. This is not supported.");
  }
  //--------------------------------------------------

  if( param->exists("zero_time") ) {
    param->getFloat("zero_time",&timeZero_ms);
    if( param->getNumValues("zero_time") > 1 ) {
      string unit;
      param->getString( "zero_time", &unit, 1 );
      if( !unit.compare("samples") ) {
        timeZero_ms *= asciiParam.sampleInt;
      }
      else if( !unit.compare("s") ) {
        timeZero_ms *= 1000.0;
      }
      else if( !unit.compare("ms") ) {
        // Nothing
      }
    }
  }
  if( edef->isDebug() ) {
    writer->line("Zero time = %f ms", timeZero_ms);
  }


  // Create all desiignature objects: One per input signature wavelet
  vars->fftDesigList = new csFFTDesignature*[vars->numInputWavelets];
  vars->hdrValueList = new int[vars->numInputWavelets];
  for( int iwavelet = 0; iwavelet < vars->numInputWavelets; iwavelet++ ) {
    float* bufferPtr = bufferList.at(iwavelet);
    vars->fftDesigList[iwavelet] = new csFFTDesignature( shdr->numSamples, shdr->sampleInt, bufferPtr, timeZero_ms/1000.0, percWhiteNoise, bufferOut );
    vars->hdrValueList[iwavelet] = hdrValueList.at(iwavelet);
    vars->fftDesigList[iwavelet]->setDesigFilterType( filterType );
    delete [] bufferPtr;
  }
  bufferList.clear();
  hdrValueList.clear();
  //  vars->fftDesig = new csFFTDesignature( shdr->numSamples, shdr->sampleInt, buffer, timeZero_ms/1000.0, percWhiteNoise, bufferOut );

  //  delete [] buffer;
  if( bufferOut != NULL ) delete [] bufferOut;

  float freqNy = 500.0f/asciiParam.sampleInt;
  float slopeDefault = 30;
  if( param->exists("highpass") ) {
    float cutOff;
    float slope = slopeDefault;
    param->getFloat("highpass", &cutOff );
    if( cutOff < 0 || cutOff > freqNy ) {
      writer->error("Filter frequency exceeds valid range (0-%fHz): %fHz", freqNy, cutOff );
    }
    if( param->getNumValues("highpass") > 1 ) {
      param->getFloat("highpass", &slope, 1 );
    }
    for( int iwavelet = 0; iwavelet < vars->numInputWavelets; iwavelet++ ) {
      vars->fftDesigList[iwavelet]->setDesigHighPass( cutOff, slope );
    }
  }
  if( param->exists("lowpass") ) {
    float cutOff;
    float slope = slopeDefault;
    param->getFloat("lowpass", &cutOff );
    if( param->getNumValues("lowpass") > 1 ) {
      param->getFloat("lowpass", &slope, 1 );
    }
    if( cutOff < 0 || cutOff > freqNy ) {
      writer->error("Filter frequency exceeds valid range (0-%fHz): %fHz", freqNy, cutOff );
    }
    for( int iwavelet = 0; iwavelet < vars->numInputWavelets; iwavelet++ ) {
      vars->fftDesigList[iwavelet]->setDesigLowPass( cutOff, slope );
    }
  }

  if( param->exists("high_set") ) {
    float freqHighSet = 0;
    param->getFloat("high_set", &freqHighSet);
    for( int iwavelet = 0; iwavelet < vars->numInputWavelets; iwavelet++ ) {
      vars->fftDesigList[iwavelet]->setDesigHighEnd( freqHighSet );
    }
  }

  //--------------------------------------------------
  if( param->exists("notch_suppression") ) {
    float notchFreq;
    float notchWidth;
    int numLines = param->getNumLines("notch_suppression");
    for( int iline = 0; iline < numLines; iline++ ) {
      param->getFloatAtLine("notch_suppression", &notchFreq, iline, 0);
      param->getFloatAtLine("notch_suppression", &notchWidth, iline, 1);
      for( int iwavelet = 0; iwavelet < vars->numInputWavelets; iwavelet++ ) {
        vars->fftDesigList[iwavelet]->setNotchSuppression( notchFreq, notchWidth );
      }
    }
  }


  // Write filter to output file
  if( param->exists("filename_output") ) {
    param->getString("filename_output", &vars->filenameOut);
    if( !csFileUtils::createDoNotOverwrite( vars->filenameOut ) ) {
      writer->error("Unable to open filter output file %s. Wrong path name..?", vars->filenameOut.c_str() );
    }
    if( param->getNumValues("filename_output") > 1 ) {
      param->getString("filename_output", &text, 1);
      if( !text.compare("wavelet") ) {
        vars->isFileOutWavelet = true;
      }
      else if( !text.compare("spectrum") ) {
        vars->isFileOutWavelet = false;
      }
      else {
        writer->error("Unknown option: '%s'", text.c_str());
      }
      if( param->getNumValues("filename_output") > 2 ) {
        param->getFloat("filename_output", &vars->filterTimeShift_s, 2);
        vars->filterTimeShift_s /= 1000.0; // Convert from [s] to [ms]
      }
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_designature_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  if( vars->filenameOut.size() > 0 ) {
    FILE* fout = fopen(vars->filenameOut.c_str(), "w");
    //    for( int iwavelet = 0; iwavelet < vars->fftDesigList->size(); iwavelet++ ) {
    if( vars->isFileOutWavelet ) {
      vars->fftDesigList[0]->dump_wavelet(fout,true,vars->filterTimeShift_s);
    }
    else {
      vars->fftDesigList[0]->dump_spectrum(fout);
    }
    fclose(fout);
    fout = NULL;
    vars->filenameOut = ""; // Only dump output file once
  }

  float* samples = trace->getTraceSamples();
  // More than one input wavelet
  if( vars->numInputWavelets > 1 ) {
    // Check if trace header value matches wavelet identifier number:
    int hdrValue = trace->getTraceHeader()->intValue( vars->hdrId_input );
    if( hdrValue != vars->hdrValueList[vars->currentIndexHdrValue] ) {
      bool found = false;
      for( int iwavelet = vars->currentIndexHdrValue+1; iwavelet < vars->numInputWavelets; iwavelet++ ) {
        if( vars->hdrValueList[iwavelet] == hdrValue ) {
          vars->currentIndexHdrValue = iwavelet;
          found = true;
          break;
        }
      }
      if( !found ) {
        for( int iwavelet = 0; iwavelet < vars->currentIndexHdrValue;iwavelet++ ) {
          if( vars->hdrValueList[iwavelet] == hdrValue ) {
            vars->currentIndexHdrValue = iwavelet;
            found = true;
            break;
          }
        } // END for
        if( !found ) {
          writer->error("No wavelet found for trace with trace header value = %d   (trace header name: '%s')", hdrValue, vars->hdrName_input.c_str() );
        }
      } // END !found
    } // END: reset currentIndexHdrValue
  } // END: More than one wavelet

  vars->fftDesigList[vars->currentIndexHdrValue]->applyFilter( samples, shdr->numSamples );

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_designature_( csParamDef* pdef ) {
  pdef->setModule( "DESIGNATURE", "Designature filter operation" );
  pdef->addDoc("This module generates a filter (in the frequency domain) from a specified signature wavelet, and applies this filter to all input traces.");
  pdef->addDoc("Currently supported options: Create zero-phasing filter that removes the specified signature; in other words, apply the transfer function that converts the specified signature wavelet into a zero-phase spike.");
  pdef->addDoc("Add specified white noise to signature spectrum before computing transfer function. Finally, the designature filter spectrum can be further bandpass filtered.");

  pdef->addDoc("Note that this is a naive implementation, using spectral division. The user has the choice to limit the frequency range, add white noise, and to fill notches to alleviate artefacts.");


  pdef->addParam( "input_wavelet", "Name of file containing input signature/response wavelet", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit of time values found in input file, if any" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );

  pdef->addParam( "format", "Input wavelet ASCII file format", NUM_VALUES_VARIABLE);
  pdef->addValue( "columns", VALTYPE_OPTION );
  pdef->addOption( "signature", "Read in source signature from Nucleus ASCII file" );
  pdef->addOption( "columns", "Simple file format with 2 or 3 columns:  Time[ms]  Amplitude  (Trace number)",
    "Trace number column is optional." );
  pdef->addValue( "1", VALTYPE_NUMBER, "For 'column' format: Column number specifying time (or other sample unit)" );
  pdef->addValue( "2", VALTYPE_NUMBER, "For 'column' format: Column number specifying sample value" );
  pdef->addValue( "-1", VALTYPE_NUMBER, "Optional: Column number specifying trace header value identifying wavelet trace. -1: Not applicable" );

  pdef->addParam( "input_hdr", "Name of trace header that is used to identify wavelet trace in input file", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING );

  pdef->addParam( "output_wavelet", "Name of file containing output wavelet", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit of time values found in input file, if any" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );

  pdef->addParam( "zero_time", "Zero time in input wavelet", NUM_VALUES_VARIABLE,
                  "Specifying this parameter overrides the zero time that is found in the ASCII signature file" );
  pdef->addValue( "0", VALTYPE_NUMBER );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );
  pdef->addOption( "samples", "Samples" );

  pdef->addParam( "white_noise", "Amount of white noise (in percent of maximum amplitude) to add to the signature/response spectrum",
                  NUM_VALUES_FIXED );
  pdef->addValue( "0.01", VALTYPE_NUMBER );

  pdef->addParam( "notch_suppression", "Suppress notches in signature wavelet", NUM_VALUES_FIXED,
                  "Apply cosine taper around notch frequency to filter" );
  pdef->addValue( "", VALTYPE_NUMBER, "Notch frequency [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Width of suppression filter [Hz]" );

  pdef->addParam( "option", "Filter type option", NUM_VALUES_FIXED);
  pdef->addValue( "amp_phase", VALTYPE_OPTION );
  pdef->addOption( "amp_phase", "Create amplitude & phase designature filter" );
  pdef->addOption( "amp_only", "Create amplitude only designature filter" );
  pdef->addOption( "phase_only", "Create phase only designature filter" );

  pdef->addParam( "lowpass", "Lowpass filter to apply to designature filter before application", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Cutoff frequency for low-pass filter [Hz]",
     "The cutoff frequency will be damped by -3db" );
  pdef->addValue( "30", VALTYPE_NUMBER, "Filter slope (db/oct)" );

  pdef->addParam( "highpass", "Highpass filter to apply to designature filter before application", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Cutoff frequency for highpass filter [Hz]",
     "The cutoff frequency will be damped by -3db" );
  pdef->addValue( "30", VALTYPE_NUMBER, "Filter slope (db/oct)" );

  pdef->addParam( "high_set", "Set high-end of filter constant", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Frequency [Hz]", "Above this frequency, keep amplitude and phase constant");

  pdef->addParam( "filename_output", "Name of file where designature filter shall be written to", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addValue( "spectrum", VALTYPE_OPTION );
  pdef->addOption( "spectrum", "Output filter as amplitude/phase spectrum" );
  pdef->addOption( "wavelet", "Output filter as wavelet" );
  pdef->addValue( "0", VALTYPE_STRING, "Time shift [ms] to apply to filter wavelet" );
  //  pdef->addValue( "no", VALTYPE_OPTION );
  //  pdef->addOption( "yes", "Apply phase shift to filter wavelet to place main energy at center" );
  //  pdef->addOption( "no", "Do not center-shift wavelet, output as-is" );

  pdef->addParam( "override_sample_int", "Override sample interval?", NUM_VALUES_FIXED);
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "" );
  pdef->addOption( "yes", "Ignore sample interval of input wavelet. Assume it is the same as the input data" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_designature_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_designature::VariableStruct* vars = reinterpret_cast<mod_designature::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_designature_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_designature::VariableStruct* vars = reinterpret_cast<mod_designature::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->fftDesigList != NULL ) {
    for( int iwavelet = 0; iwavelet < vars->numInputWavelets; iwavelet++ ) {
      delete vars->fftDesigList[iwavelet];
    }
    delete [] vars->fftDesigList;
    vars->fftDesigList = NULL;
  }
  if( vars->hdrValueList != NULL ) {
    delete [] vars->hdrValueList;
    vars->hdrValueList = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_designature_( csParamDef* pdef ) {
  params_mod_designature_( pdef );
}
extern "C" void _init_mod_designature_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_designature_( param, env, writer );
}
extern "C" bool _start_exec_mod_designature_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_designature_( env, writer );
}
extern "C" void _exec_mod_designature_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_designature_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_designature_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_designature_( env, writer );
}

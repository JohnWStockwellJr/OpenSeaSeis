/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csTableNew.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: MUTE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_mute {
  struct VariableStruct {

    bool debug; 

    int domain;         // Domain of the input trace.
    int fftDataType;    // FFT data type if input data transformed.

    int inLength;       // Input length. Actual number of mutable samples.

    string tableName;   // Name of user table.
    csTableNew* table;  // User's Table containing keys & mute value(s).
    int* hdrId_keys;    // Key's header IDs 

    int tableValueIndex1;   // Column in table with first mute value
    int tableValueIndex2;   // Column in table with (optional) second mute value

    int mode;           // Mute mode: FRONT or END
    int nValued;        // Number of mute values per trace. For now, only 1 or 2

    float mute_time;    // Current mute value in units of trace domain
    float mute_time2;   // Current mute value in units of trace domain (optional)
    int mute_start_samp; // Starting sample of mute
    int mute_end_samp;   // Ending sample of mute

    bool killZeroTraces; // Kill trace if the mute will completely zero the trace.

    int taperType;          // Type of taper
    string taperType_str;   // Type of taper
    int taperLengthSamples; // Length of taper
    int taperApply;         // where is the taper applied
    string taperApply_str;  // where is the taper applied

    bool  indicate;             // Flag to delineate the mute area with spikes rather than actually mute
    float indicateValue;        // Value of spike
    int   indicateWidthSamples; // Width of spike

    int windowStartSample;
    int windowEndSample;
  };
  static const string MY_NAME = "mute";

  static int const TAPER_APPLY_MIN    = 1;
  static int const TAPER_APPLY_CENTER = 2;
  static int const TAPER_APPLY_MAX    = 3;

  static int const TAPER_LINEAR = 1;
  static int const TAPER_COSINE = 2;

  static int const APPLY_INSIDE  = 1;
  static int const APPLY_OUTSIDE = 2;
  static int const APPLY_FRONT   = 3;
  static int const APPLY_END     = 4;

  // Utility functions
  void  do_mute( int const startSamp, int const endSamp, float*samples ){
    for( int isamp = startSamp; isamp <= endSamp; isamp++ ) {samples[isamp] = 0.0; }
  }

  void do_taper_front( int const type, int const length,  
                       int const startSamp, int const endSamp, float*samples, 
                       int const zero_sample, int const final_sample ){

    // LINEAR Taper
    if( type == TAPER_LINEAR ) {
      for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
        if ( isamp < zero_sample ) continue;
        if ( isamp > final_sample ) break;
        float weight = (float)(isamp-startSamp+1)/(float)length;
        samples[isamp] *= weight;
      }

      // COSINE taper
    } else if( type == TAPER_COSINE ) {
      for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
        if ( isamp < zero_sample ) continue;
        if ( isamp > final_sample ) break;
        float phase  = ( ( (float)(isamp-startSamp) / (float)length ) - 1.0 ) * M_PI;
        float weight = 0.5 * (cos(phase) + 1.0);
        samples[isamp] *= weight;
      }

    }
    return;
  }

  void do_taper_end( int const type, int const length,  
                     int const startSamp, int const endSamp, float*samples,
                     int const zero_sample, int const final_sample ){
    // LINEAR Taper
    if( type == TAPER_LINEAR ) {
      for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
        if ( isamp < zero_sample ) continue;
        if ( isamp > final_sample ) break;
        float weight = (float)(endSamp-isamp)/(float)length;
        samples[isamp] *= weight;
      }      
      
      // COSINE taper
    } else if( type == TAPER_COSINE ) {
      for( int isamp = startSamp; isamp < endSamp; isamp++ ) {
        if ( isamp < zero_sample ) continue;
        if ( isamp > final_sample ) break;
        float phase  = ( ( (float)(endSamp-isamp) / (float)length ) - 1.0 ) * M_PI;
        float weight = 0.5 * (cos(phase) + 1.0);
        samples[isamp] *= weight;
      }
      
    }
    return;
  }

} 

using namespace mod_mute ;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_mute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->inLength         = shdr->numSamples;
  vars->table            = NULL;
  vars->hdrId_keys       = NULL;
  vars->tableValueIndex1 = -1;
  vars->tableValueIndex2 = -1;
  vars->mute_start_samp  = 0;
  vars->mute_end_samp    = 0;
  vars->killZeroTraces   = false; 
  vars->indicate         = false;
  vars->indicateValue    = 0.0;
  vars->indicateWidthSamples = 1;
  vars->nValued          = 0;
  float taperLength      = 0; // See below
  vars->taperType        = TAPER_COSINE;
  vars->taperType_str    = "COSINE";
  vars->taperApply       = TAPER_APPLY_MIN;
  vars->taperApply_str   = "MINIMUM";
  vars->windowStartSample = 0; // By default, apply mute to full trace length
  vars->windowEndSample   = shdr->numSamples-1; // By default, apply mute to full trace length

  std::string text;

  // Check for valid domains.
  vars->domain      = shdr->domain;
  vars->fftDataType = shdr->fftDataType;
  if ( shdr->domain == DOMAIN_XT || shdr->domain == DOMAIN_KT || shdr->domain == DOMAIN_XD ) {
    writer->line("Input data is in TIME/DEPTH domain.");
  }
  else if ( shdr->domain == DOMAIN_FX ){
    writer->line("Input data is in FX domain data.");
    
    if ( shdr->fftDataType == FX_AMP_PHASE ){
      writer->line("Muting amplitudes only of FX domain AMPLITUDE & PHASE data.");
      vars->inLength = vars->inLength/2;
      if ( 2*vars->inLength != shdr->numSamples ) {
        writer->warning("For FX data, typically expect number of input samples to be even but that is not the case.", text.c_str());
      }

    } else if ( shdr->fftDataType == FX_AMP ){
      writer->line("Muting amplitudes of FX domain AMPLITUDE-ONLY data.");

    } else if ( shdr->fftDataType == FX_REAL_IMAG ){
      writer->warning("Applying mute to FX domain REAL-IMAGINARY data.");

    } else if ( shdr->fftDataType == FX_PSD ){
      writer->warning("Applying mute to FX domain PSD data.");

    } else {
      writer->line("WARN:Unknown FX domain type. Will process as if time. Proceed at your own risk!");
    }      

  } else if ( shdr->domain == DOMAIN_FK ){
    writer->line("Input data is in FX domain data.");

    if ( shdr->fftDataType == FK_AMP_PHASE ){
      writer->line("Muting amplitudes only of FK domain AMPLITUDE & PHASE data.");
      vars->inLength = vars->inLength/2;
      if ( 2*vars->inLength != shdr->numSamples ) {
        writer->warning("For FX data, typically expect number of input samples to be even but that is not the case.", text.c_str());
      }

    } else if ( shdr->fftDataType == FK_AMP ){
      writer->line("Muting amplitudes of FK domain AMPLITUDE-ONLY data.");

    } else if ( shdr->fftDataType == FK_REAL_IMAG ){
      writer->warning("Applying mute to FK domain REAL-IMAGINARY data.");

    } else if ( shdr->fftDataType == FX_PSD ){
      writer->warning("Applying mute to FK domain PSD data.");

    } else {
      writer->warning("Unknown FK domain type. Will process as if time. Proceed at your own risk!");
    }      

  } else {
    writer->warning("Unknown domain type. Will process as if time. Proceed at your own risk!");
  }

  // Retrieve mute table
  if( param->exists("table") ) {
    param->getString("table", &text );
    vars->tableName = text;

    vars->table = new csTableNew( csTableNew::TABLE_TYPE_UNIQUE_KEYS );

    // Require "table_key" info for table option
    if ( !param->exists("table_key") ){
      writer->line("ERROR:Parameter 'table_key' required with 'table' option.");
      env->addError();
    } else {
      int numKeys = param->getNumLines("table_key");
      if( numKeys == 0 ) {
        writer->line("ERROR:No table_key(s) specified.");
        env->addError(); 
      }      
      vars->hdrId_keys = new int[numKeys];
      bool interpolate = true;
      for( int ikey = 0; ikey < numKeys; ikey++ ) {
        std::string headerName;
        int col;
        param->getStringAtLine( "table_key", &headerName, ikey, 0 );
        param->getIntAtLine( "table_key", &col, ikey, 1 );
        vars->table->addKey( col-1, interpolate );  // -1 to convert from 'user' column to 'C' column
        if( !hdef->headerExists( headerName ) ) {
          writer->line("ERROR:No matching trace header found for table_key = '%s'", headerName.c_str() );
          env->addError();    
        } else {
          vars->hdrId_keys[ikey] = hdef->headerIndex( headerName );
        }
      } // END for ikey
    }

    // Require mute value info for table option
    int numValues = param->getNumValues("table_mute");
    if( numValues == 0 ) {
      writer->line("ERROR:No 'table_mute' specified.");
      env->addError();    
    } else if( numValues > 2 ) {
      writer->line("ERROR:Only two 'table_mute' columns supported right now.");
      env->addError();    
    } else {
      vars->nValued = numValues;
      int col;
      param->getInt( "table_mute", &col );
      vars->table->addValue( col-1 );
      vars->tableValueIndex1 = 0;
      if ( numValues == 2 ){
        param->getInt( "table_mute", &col, 1 );      
        vars->table->addValue( col-1 );
        vars->tableValueIndex2 = 1;
      }
    }

    // Init the table if no errors previously
    if ( env->errorCount() == 0 ){

      bool sortTable = false;  // Skip sorting for now. 
      param->getString("table", &text );
      try {
        vars->table->initialize( text, sortTable );
      }
      catch( csException& exc ) {
        writer->line("ERROR:Initializing input table '%s': %s\n", text.c_str(), exc.getMessage() );
        env->addError();    
      }

    } else {
      writer->line("ERROR:Skipping initialization table due to previous errors\n", text.c_str());
      env->addError();          
    }

  } else {
    if ( param->exists("table_key") || param->exists("table_mute") ){
      writer->line("Parameters 'table_key'/'table_mute' ignored since 'table' not specified.");
    }
  }

  // User must specify constant mute if no table specified.
  if( !vars->table ) {
    if ( param->exists("time") ){

      float recordLength = (float)(vars->inLength-1) * shdr->sampleInt;

      vars->nValued++;
      param->getFloat( "time", &vars->mute_time );
      if( vars->mute_time < 0 || vars->mute_time > recordLength ) {
        writer->line("ERROR:Specified 'time' (%.0f) is outside of valid range (0-%.0f).",
                  vars->mute_time, vars->inLength*shdr->sampleInt );
        env->addError();    
      }

      if ( param->exists("time2") ){
        vars->nValued++;
        param->getFloat( "time2", &vars->mute_time2 );
        if( vars->mute_time2 < 0 || vars->mute_time2 > recordLength ) {
          writer->line("ERROR:Specified 'time2' (%.0f) is outside of valid range (0-%.0f).",
                    vars->mute_time2, recordLength );
          env->addError();    
        }
      }
    } else if ( param->exists("time2") ){
      writer->line("ERROR:Cannot specify 'time2' without 'time'.");
      env->addError();
    } else {
      writer->line("ERROR:Must specify only one of either 'table' or 'time' parameter.");
      env->addError();
    }

  } else if ( param->exists("time") ||  param->exists("time2") ){
    writer->line("Parameters 'time'/'time2' ignored since 'table' already specified.");
  }

  // For now only 1 or 2 mute values per trace
  if ( vars->nValued < 0 || vars->nValued > 2 ){
    writer->line("ERROR:Currently can only specify one or two mute values per trace.");
    env->addError();
  }

  // Fetch mute mode
  if( param->exists("mode") ) {
    param->getString("mode",&text);
    text = toLowerCase( text );    
    if( text.compare("front") == 0 ) {
      vars->mode = APPLY_FRONT;
    } else if( text.compare("end") == 0 ) {
      vars->mode = APPLY_END;
    } else if( text.compare("inside") == 0 ) {
      vars->mode = APPLY_INSIDE;
    } else if( text.compare("outside") == 0 ) {
      vars->mode = APPLY_OUTSIDE;
    } else {
      vars->mode = 0;
      writer->line("ERROR:Unknown 'mode' option: '%s'", text.c_str());
      env->addError();    
    }

    // Default mode
  } else {
    if ( vars->nValued == 1 ){
      vars->mode = APPLY_FRONT;
    } else {
      vars->mode = APPLY_OUTSIDE;
    }
  }

  // Evaluate parameters so far. 
  if ( vars->mode == APPLY_FRONT || vars->mode == APPLY_END ) {
    if ( vars->nValued != 1 ){
      writer->line("ERROR:Cannot specify 'mode = %s' if more than one mute value per trace.", text.c_str() );
      env->addError();          
    }
    int sampleIndex = (int)( vars->mute_time / shdr->sampleInt );
    if( vars->mode == APPLY_FRONT ) {
      vars->mute_start_samp = 0;
      vars->mute_end_samp   = sampleIndex-1;
    } else if( vars->mode == APPLY_END ) {
      vars->mute_start_samp = sampleIndex+1;
      vars->mute_end_samp   = vars->inLength-1;
    }

  } else if ( vars->mode == APPLY_INSIDE || vars->mode == APPLY_OUTSIDE ) {
    if ( vars->nValued == 1 ){
      writer->line("ERROR:Cannot specify 'mode = %s' if only one mute value per trace.", text.c_str() );
      env->addError();          
    }
    vars->mute_start_samp = (int)( vars->mute_time / shdr->sampleInt );
    vars->mute_end_samp   = (int)( vars->mute_time2 / shdr->sampleInt );
  }

  // Taper information
  if( param->exists("taper_len") ) {param->getFloat("taper_len", &taperLength );}
  vars->taperLengthSamples = (int)( taperLength/shdr->sampleInt + 0.5 );

  if( param->exists("taper_type") ) {
    param->getString("taper_type", &text);
    text = toLowerCase( text );    
    if( text.compare("cos") == 0 ) {
      vars->taperType = TAPER_COSINE;
    } else if( text.compare("linear") == 0 ) {
      vars->taperType = TAPER_LINEAR;
    } else {
      writer->line("ERROR:Unknown 'taper_type': %s", text.c_str());
      env->addError();    
    }
    vars->taperType_str = toUpperCase( text );
  }

  if( param->exists("taper_apply") ) {
    param->getString("taper_apply", &text);
    text = toLowerCase( text );    
    if( text.compare("minimum") == 0 ) {
      vars->taperApply = TAPER_APPLY_MIN;
    } else if( text.compare("center") == 0 ) {
      vars->taperApply = TAPER_APPLY_CENTER;
    } else if( text.compare("maximum") == 0 ) {
      vars->taperApply = TAPER_APPLY_MAX;
    } else {
      writer->line("ERROR:Unknown 'taper_apply': %s", text.c_str());
      env->addError();    
    }
    vars->taperApply_str = toUpperCase( text );
  }

  if ( vars->taperLengthSamples > 0 ){
    writer->line("A %s taper of length %d samples will be applied at the %s of the mute.", 
              vars->taperType_str.c_str(), 
              vars->taperLengthSamples,
              vars->taperApply_str.c_str() );
  }

  if( param->exists("window") ) {
    param->getString("window",&text);
    if( !text.compare("no") ) {
      // Nothing
    }
    else if( !text.compare("yes") ) {
      float startTime;
      float endTime;
      param->getFloat( "window", &startTime, 1 );
      param->getFloat( "window", &endTime, 2 );
      vars->windowStartSample = std::max( 0, (int)round( startTime / (float)shdr->sampleInt ) );
      vars->windowEndSample   = std::min( shdr->numSamples-1, (int)round( endTime / (float)shdr->sampleInt ) );
      if( vars->windowStartSample >= shdr->numSamples ) writer->error("Window start time (=%f) exceeds trace length (=%f)", startTime, shdr->numSamples );
      if( vars->windowEndSample <= vars->windowStartSample ) writer->error("Window start time (=%f) exceeds or equals window end time (=%f)", startTime, endTime );
    }
    else {
      writer->error("Unknown mode option: '%s'", text.c_str());
    }
  }

  // Kill traces if they are completely muted?
  if( param->exists("kill") ) {
    param->getString("kill",&text);
    text = toLowerCase( text );    
    if( text.compare("yes") == 0 ) {
      vars->killZeroTraces = true;
    } else if( text.compare("no") == 0 ) {
      vars->killZeroTraces = false;
    } else {
      writer->line("ERROR:Unknown 'kill' option: '%s'", text.c_str());
      env->addError();    
    }
  }

  // Delineate the mute area with a "spike" rather than actually mute it. 
  if( param->exists("indicate") ) {
    vars->indicate = true;
    param->getFloat( "indicate", &vars->indicateValue );
    if( param->getNumValues("indicate") > 1 ) {
      param->getInt( "indicate", &vars->indicateWidthSamples, 1 );
    }
  }

}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_mute_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  float* samples = trace->getTraceSamples();

  // Fetch mutes from table is specified
  if( vars->table != NULL ) {
    double* keyValueBuffer = new double[vars->table->numKeys()];
    for( int ikey = 0; ikey < vars->table->numKeys(); ikey++ ) {
      keyValueBuffer[ikey] = trace->getTraceHeader()->doubleValue( vars->hdrId_keys[ikey] );
    }
    double value, value2=-1;
    try {
      value = vars->table->getValue( keyValueBuffer, vars->tableValueIndex1 );
    }
    catch( csException& e ) {
      delete [] keyValueBuffer;
      writer->error("MUTE:ERROR:Retrieving mute value from table %s", e.getMessage());
      throw(e);
    }
    if ( vars->nValued == 2 ){
      try {
        value2 = vars->table->getValue( keyValueBuffer, vars->tableValueIndex2 );
      }
      catch( csException& e ) {
        delete [] keyValueBuffer;
        writer->error("MUTE:ERROR:Retrieving second mute value from table %s", e.getMessage());
        throw(e);
      }
    }
    delete [] keyValueBuffer;  
    vars->mute_time  = (float)value;
    vars->mute_time2 = (float)value2;
  }

  // If requested, apply indicate option & return
  if( vars->indicate ) { 
    int sampleIndex = (int)( vars->mute_time / shdr->sampleInt );
    int endSamp = std::min( sampleIndex+vars->indicateWidthSamples, vars->inLength-1 );
    for( int isamp = std::max( sampleIndex-vars->indicateWidthSamples, 0); isamp < endSamp; isamp++ ) {
      samples[isamp] = vars->indicateValue;
    }
    if ( vars->nValued == 2 ){
      int sampleIndex = (int)( vars->mute_time2 / shdr->sampleInt );
      int endSamp = std::min( sampleIndex+vars->indicateWidthSamples, vars->inLength-1 );
      for( int isamp = std::max( sampleIndex-vars->indicateWidthSamples, 0); isamp < endSamp; isamp++ ) {
        samples[isamp] = vars->indicateValue;
      }      
    }
    return;
  }

  // Note: Taper_start/_end are samples, so taper_start is always < taper_end. 
  // The apply routines will figure it out. 
  int mute_start, mute_end;
  int zero = 0;
  mute_start = mute_end = zero;

  int taperlength = vars->taperLengthSamples;

  // Single mute value: FRONT/END muting
  if (vars->nValued == 1 ){
    int taper_start, taper_end; 

    int sampleIndex = (int)( vars->mute_time / shdr->sampleInt );

    taper_start = taper_end = sampleIndex;

    // FRONT mute
    if( vars->mode == APPLY_FRONT ) {
      vars->mute_start_samp = 0;
      vars->mute_end_samp   = sampleIndex-1;

      mute_start = vars->mute_start_samp;
      mute_end = vars->mute_end_samp;

      // Taper?
      if ( taperlength > 0 ){
        mute_start = vars->mute_start_samp;

        if ( vars->taperApply == TAPER_APPLY_MAX ){
          taper_start = mute_end - taperlength + 1;
          
        } else if ( vars->taperApply == TAPER_APPLY_CENTER ){ 
          taper_start = mute_end - taperlength/2 + 1;
          
        } else if ( vars->taperApply == TAPER_APPLY_MIN ){ 
          taper_start = mute_end+1;
        }
        taper_end   = taper_start + taperlength;
        mute_end    = taper_start - 1; 

        // nothing muted
        if ( taper_end < 0 ){ return; }

        // Do taper here
        do_taper_front( vars->taperType, vars->taperLengthSamples, 
                        taper_start, taper_end, samples, 
                        zero, vars->inLength-1 );
      }

      // END mute
    } else if( vars->mode == APPLY_END ) {
      vars->mute_start_samp = sampleIndex;
      vars->mute_end_samp   = vars->inLength-1;

      mute_start = vars->mute_start_samp;
      mute_end = vars->mute_end_samp;

      // Taper?
      if ( taperlength > 0 ){
        mute_end = vars->mute_end_samp;

        if ( vars->taperApply == TAPER_APPLY_MAX ){
          taper_start = mute_start;
          
        } else if ( vars->taperApply == TAPER_APPLY_CENTER ){ 
          taper_start = mute_start - taperlength/2;
          
        } else if ( vars->taperApply == TAPER_APPLY_MIN ){ 
          taper_start = mute_start - taperlength;
          
        }
        taper_end   = taper_start + taperlength;
        mute_start  = taper_end;

        // nothing muted
        if ( taper_start > vars->inLength-1 ){return ;}

        // Do taper here
        do_taper_end( vars->taperType, vars->taperLengthSamples, 
                      taper_start, taper_end, samples, 
                      zero, vars->inLength-1 );

      }
    }

    // Kill trace?
    if( vars->killZeroTraces && mute_end >= vars->inLength-1 && mute_start <= 0 ) {
      traceGather->freeAllTraces();
      return;
    }

    // Do hard FRONT/END mute here
    do_mute( mute_start, mute_end, samples );

    // Two-valued mute: SURGICAL muting
  } else if (vars->nValued == 2 ){
    int taper_start1, taper_end1, taper_start2, taper_end2;

    int sampleIndex = (int)( vars->mute_time / shdr->sampleInt );
    vars->mute_start_samp = std::max( vars->windowStartSample, sampleIndex );

    taper_start1 = taper_end1 = sampleIndex;

    sampleIndex = (int)( vars->mute_time2 / shdr->sampleInt );
    vars->mute_end_samp   = std::min( vars->windowEndSample+1, std::min( vars->mute_end_samp+vars->taperLengthSamples, shdr->numSamples ) );

    taper_start2 = taper_end2 = sampleIndex;

    // INSIDE mute
    if( vars->mode == APPLY_INSIDE ) {
      mute_start = vars->mute_start_samp;
      mute_end   = vars->mute_end_samp;

      // Taper?
      if ( taperlength > 0 ) {
        if ( vars->taperApply == TAPER_APPLY_MAX ){
          taper_start1 = mute_start;
          taper_start2 = mute_end - taperlength+1;

        } else if ( vars->taperApply == TAPER_APPLY_CENTER ){ 
          taper_start1 = mute_start - taperlength/2;
          taper_start2 = mute_end - taperlength/2+1;

        } else if ( vars->taperApply == TAPER_APPLY_MIN ){ 
          taper_start1 = mute_start - taperlength;
          taper_start2 = mute_end + 1;

        }
        taper_end1   = taper_start1 + taperlength;
        mute_start   = taper_end1;
        taper_end2   = taper_start2 + taperlength;
        mute_end     = taper_start2 - 1;

        // Nothing muted        
        if ( taper_start1 >= taper_end2 ){return;}

        // Deal with special case of crossing tapers!
        int taper_end_final   = vars->inLength-1;
        int  taper_front_zero = zero;
        if ( taper_end1 > taper_start2 ) {
          int cross = (taper_end1 + taper_start2)/2;
          taper_end_final = taper_front_zero = cross;
        }

        // Apply END taper at first mute
        do_taper_end( vars->taperType, vars->taperLengthSamples, 
                      taper_start1, taper_end1, samples, 
                      zero, taper_end_final );
      
        // Apply FRONT taper at second mute
        do_taper_front( vars->taperType, vars->taperLengthSamples, 
                        taper_start2, taper_end2, samples, 
                        taper_front_zero, vars->inLength-1 );
      }

      if( vars->killZeroTraces && mute_end >= vars->inLength-1 && mute_start <= 0 ) {
        traceGather->freeAllTraces();
        return;
      }

      // Hard mute between the two mutes
      do_mute( mute_start, mute_end, samples );
 
      // OUTSIDE mute
    } else if( vars->mode == APPLY_OUTSIDE ) {
      mute_start = vars->mute_start_samp - 1;
      mute_end   = vars->mute_end_samp + 1;

      // Taper?
      if ( taperlength > 0 ) {
        if ( vars->taperApply == TAPER_APPLY_MAX ){
          taper_start1 = mute_start - taperlength+1;
          taper_start2 = mute_end;

        } else if ( vars->taperApply == TAPER_APPLY_CENTER ){ 
          taper_start1 = mute_start - taperlength/2+1;          
          taper_start2 = mute_end - taperlength/2;

        } else if ( vars->taperApply == TAPER_APPLY_MIN ){ 
          taper_start1 = mute_start+1;          
          taper_start2 = mute_end - taperlength;

        }
        taper_end1   = taper_start1 + taperlength;
        mute_start   = taper_start1 - 1;
        taper_end2   = taper_start2 + taperlength;
        mute_end     = taper_end2;

        // Nothing muted
        if ( taper_end1 < 0 && taper_start2 > vars->inLength-1 ){return;}

        // Do taper here
        do_taper_front( vars->taperType, vars->taperLengthSamples, 
                        taper_start1, taper_end1, samples, 
                        zero, vars->inLength-1 );
        do_taper_end( vars->taperType, vars->taperLengthSamples, 
                      taper_start2, taper_end2, samples, 
                      zero, vars->inLength-1 );        
      }

      if( vars->killZeroTraces && mute_start >= mute_end ) {
        traceGather->freeAllTraces();
        return;
      }

      // Hard mute outside the two times
      do_mute( 0, mute_start, samples );
      do_mute( mute_end, vars->inLength-1, samples );
    }
  }  
  
  return;
}

//*************************************************************************************************
// Parameter definition
//*************************************************************************************************
void params_mod_mute_( csParamDef* pdef ) {
  pdef->setModule( "MUTE", "Mute trace data" );
  pdef->setVersion( 2, 0 );

  pdef->addDoc("Performs simple trace muting (single value per trace) or simple 'surgical' muting (two values per trace.)<br>");
  pdef->addDoc("While typically this is done in the time domain, it can be used to zero the trace samples regardless of the domain of the trace.");
  pdef->addDoc("So while the parameters and descriptions often refer to 'time' and 'milliseconds (ms)', be sure to take care that values are provided in terms of the untis appropraite to the domain of the trace.");
  pdef->addDoc("For example, for FK or FX data, the values should be provided in terms of hertz (Hz.)<br>");

  pdef->addParam( "table", "Mute table.", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Mute table file, full path. Use instead of (constant) 'time' parameters.",
                  "The mute table must have at least two columns, one giving a key and the second giving the mute in units of the current domain of the trace" );

  pdef->addParam( "table_key", "Required with 'table' option. Key trace header(s) used to match values found in specified table columns.", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name of key header" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number of the header values in input table" );
  // pdef->addValue( "no", VALTYPE_OPTION, "Interpolate based to this key?" );
  // pdef->addOption( "yes", "Use this key for interpolation of value" );
  // pdef->addOption( "no", "Do not use this key for interpolation", "The input table is expected to contain the exact key values for this trace header" );
  pdef->addParam( "table_mute", "Required with 'table' option. Column numbers for mute values.", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number of first mute value. Required. If only one, use with FRONT or END mode" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number of second mute value. Optional. Use with surgical muting: INSIDE or OUTSIDE mode" );

  pdef->addParam( "time", "Constant mute at specified value in units of the current domain of the trace. Use instead of 'table'.", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Constant mute value. Use alone with FRONT or END mode or with 'time2' for INSIDE or OUTSIDE" );

  pdef->addParam( "time2", "Second constant mute at specified value in units of the current domain of the trace. Use instead of 'table'.", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Second constant mute value. Use with surgical muting: INSIDE or OUTSIDE mode" );

  pdef->addParam( "mode", "Supply 'front' or 'end' for single-value mute, 'inside' or 'outside' for two-value mute.", NUM_VALUES_FIXED );
  pdef->addValue( "front", VALTYPE_OPTION );
  pdef->addOption( "front", "Front mute (default). If only one mute value per key, specifies that 'time' is the value of first unmuted sample", "Samples are muted from first sample to the sample at the mute value. The mute value is the first unmuted sample" );
  pdef->addOption( "end", "End mute. If only one mute value per key, specifies that 'time' is the value of last unmuted sample", "Samples are muted from sample at the mute value to the end of the trace. The mute value is the first muted sample" );
  pdef->addOption( "inside", "Inside mute. If two mute values per key, specify to mute between the two values", "Samples above and below those values are un-muted" );
  pdef->addOption( "outside", "Outside mute. If two mute values per key, specify to mute above and below the two values respectively", "Samples between those two values are un-muted" );

  pdef->addParam( "taper_len", "Taper length.", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Mute taper length" );

  pdef->addParam( "taper_type", "Type of mute taper.", NUM_VALUES_FIXED );
  pdef->addValue( "cos", VALTYPE_OPTION );
  pdef->addOption( "linear", "Apply linear taper" );
  pdef->addOption( "cos", "Apply cosine taper" );

  pdef->addParam( "taper_apply", "How to apply the taper.", NUM_VALUES_FIXED );
  pdef->addValue( "minimum", VALTYPE_OPTION );
  pdef->addOption( "minimum", "The minimim value of the taper (0) is at the mute time" );
  pdef->addOption( "center", "Taper is centered at the mute time" );
  pdef->addOption( "maximum", "The maximim value of the taper (1) is at the mute time" );


  // not fully implemented yet:
  //  pdef->addParam( "window", "Restrict mute to certain window?", NUM_VALUES_FIXED, "If specified, mute will not be applied outside of the specified window, including no tapering" );
  //  pdef->addValue( "no", VALTYPE_OPTION);
  //  pdef->addOption( "no", "Do not restrict mute to a window: Apply mute over full trace length" );
  //  pdef->addOption( "yes", "Apply mute only inside  specified window" );
  //  pdef->addValue( "", VALTYPE_NUMBER, "Start time" );
  //  pdef->addValue( "", VALTYPE_NUMBER, "End time" );

  pdef->addParam( "kill", "Kill zero traces?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Kill zero traces" );
  pdef->addOption( "yes", "Do not kill zero traces" );

  pdef->addParam( "indicate", "Do not mute, rather delineate the mute region by setting the mute samples to the given value.", NUM_VALUES_VARIABLE,
                  "Input data will not be muted. Instead, spikes are placed at the mute values, with the given amplitude." );
  pdef->addValue( "0", VALTYPE_NUMBER, "Value that mute samples are set to" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Width in samples to indicate with given value" );

}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_mute::VariableStruct* vars = reinterpret_cast<mod_mute::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_mute::VariableStruct* vars = reinterpret_cast<mod_mute::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->table != NULL ) {
    delete vars->table;
    vars->table = NULL;
  }
  if( vars->hdrId_keys != NULL ) {
    delete [] vars->hdrId_keys;
    vars->hdrId_keys = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_mute_( csParamDef* pdef ) {
  params_mod_mute_( pdef );
}
extern "C" void _init_mod_mute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_mute_( param, env, writer );
}
extern "C" bool _start_exec_mod_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_mute_( env, writer );
}
extern "C" void _exec_mod_mute_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_mute_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_mute_( env, writer );
}

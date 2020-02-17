/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstring>
#include "cseis_includes.h"
#include "geolib_endian.h"
#include "geolib_string_utils.h"
#include "csCompareVector.h"

#include "csSegyTraceHeader.h"
#include "csSegyHdrMap.h"
#include "csSUTraceManager.h"
#include "csSUArguments.h"
#include "su_modules.h"
extern "C" {
  #include <pthread.h>
}

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SU
 *
 * @author Bjorn Olofsson
 */
namespace mod_su {
  struct VariableStruct {
    int traceCounter;
    cseis_su::csSUTraceManager* suPush;
    cseis_su::csSUTraceManager* suPull;
    cseis_su::csSUArguments* args;
    pthread_t thread;
    bool isThreadAlive;
    std::string suModuleName;
    bool ignoreSUHdr;

    int numSamplesIn;
    float sampleIntIn;
    csSegyTraceHeader* suTrcHdrWrite;
    csSegyTraceHeader* suTrcHdrRead;
    csSegyHdrMap* suHdrMap;
    int*  hdrIndexSU;
    type_t* hdrTypeSU;

    int counter;
  };
  static cseis_geolib::csCompareVector<std::string>* moduleList = NULL;
  static int CALL_COUNTER = 0;
  // Idle sleep time while pushing traces from Seaseis to SU, in [microseconds]
  static int PUSH_IDLE_SLEEP_TIME_US = 1000;
}
using mod_su::VariableStruct;

void setHeaders( csSegyTraceHeader* segyTrcHdr, csSegyHdrMap* segyHdrMap, int* hdrIndexSegy, type_t* hdrTypeSegy,
                 int numSamples, float sampleInt, csTraceHeaderDef* hdef, int scalarPolarity );
void performPull( mod_su::VariableStruct* vars,
                  csTraceGather* traceGather,
                  csTraceHeaderDef const* hdef,
                  csSuperHeader const* shdr,
                  bool isDebug );

//*************************************************************************************************
// Init phase
//
void init_mod_su_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  mod_su::CALL_COUNTER += 1;
  vars->counter = mod_su::CALL_COUNTER;

  vars->traceCounter    = 0;
  edef->setExecType( EXEC_TYPE_MULTITRACE );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->suPush = new cseis_su::csSUTraceManager();
  vars->suPull = new cseis_su::csSUTraceManager();
  vars->suPull->setLogFile( log->getFile() );
  vars->args   = new cseis_su::csSUArguments();  // MUST be allocated on heap!! Memory on stack cannot be used by other threads
  vars->ignoreSUHdr = false;

  vars->numSamplesIn = shdr->numSamples;
  vars->sampleIntIn = shdr->sampleInt;
  vars->suTrcHdrWrite = NULL;
  vars->suTrcHdrRead  = NULL;
  vars->suHdrMap      = NULL;
  vars->hdrIndexSU = NULL;
  vars->hdrTypeSU  = NULL;
  vars->isThreadAlive = false;
  vars->suModuleName = "";

  std::string text;
  int suModuleIndex = -1;

  cseis_geolib::csVector<std::string> argvList;

  param->getString( "name", &text );
  vars->suModuleName = cseis_geolib::toLowerCase( text );
  for( int i = 0; i < NUM_SU_MODULES; i++ ) {
    if( !vars->suModuleName.compare(SU_MODULE_NAMES[i]) ) {
      suModuleIndex = i;
      break;
    }
  }
  if( suModuleIndex < 0 ) {
    log->error("SU module not found: %s.\nSee help for available modules", vars->suModuleName.c_str());
  }

  if( param->exists("ignore_su_hdr") ) {
    param->getString("ignore_su_hdr", &text );
    if( !text.compare("yes") ) {
      vars->ignoreSUHdr = true;
    }
    else if( !text.compare("no") ) {
      vars->ignoreSUHdr = false;
    }
    else {
      log->error("Unknown option: %s", text.c_str() );
    }
  }

  if( param->exists("param") ) {
    int numLines = param->getNumLines("param");
    cseis_geolib::csVector<std::string> tempArgvList;
    for( int iline = 0; iline < numLines; iline++ ) {
      param->getAll( "param", &tempArgvList, iline );
      // Split up each argument in tempArgvList into finer arguments.
      // This is only necessary if the user specified several separate arguments within double quotes
      for( int iarg = 0; iarg < tempArgvList.size(); iarg++ ) {
        csVector<std::string> tokenList;
        tokenList.clear();
        tokenize( tempArgvList.at(iarg).c_str(), tokenList );
        for( int itoken = 0; itoken < tokenList.size(); itoken++ ) {
          argvList.insertEnd( tokenList.at(itoken) );
        }
      } // END for iarg
    } // END for iline
  }
  argvList.insertStart( vars->suModuleName );
  log->write("Full SU command line:\n   '");
  for( int iarg = 0; iarg < argvList.size(); iarg++ ) {
    if( iarg > 0 ) log->write(" ");
    log->write("%s", argvList.at(iarg).c_str() );
  }
  log->write("'\n");

  vars->args->cs2su = vars->suPush;
  vars->args->su2cs = vars->suPull;
  vars->args->setArgv( &argvList );
  if( edef->isDebug() ) vars->args->debugFlag = 1;

  if( param->exists("nsamples") ) {
    param->getInt( "nsamples", &shdr->numSamples );
  }
  if( param->exists("sample_int") ) {
    param->getFloat( "sample_int", &shdr->sampleInt );
  }

  // Create separate thread where SU module/process will run
  // SU module (csTraceManager) will wait in idle loop until new trace is received
  // SU thread will terminate when SeaSeis has indicated that no more traces will be pushed (csSUTraceManager::setEOF)
  int rc = pthread_create( &vars->thread, NULL, SU_MAIN_METHODS[suModuleIndex], (void*)vars->args );
  if( rc ) {
    log->error("Could not launch thread for module '%s': Return code from pthread_create() is %d\n", vars->suModuleName.c_str(), rc);
  }

  log->line("Successfully launched thread for module '%s': Return code from pthread_create() is %d\n", vars->suModuleName.c_str(), rc);
  vars->isThreadAlive = true;

  //--------------------------------------------------
  if( mod_su::moduleList == NULL ) {
    mod_su::moduleList = new cseis_geolib::csCompareVector<std::string>();
  }
  if( mod_su::moduleList->contains( vars->suModuleName ) ) {
    log->error("This flow already contains a call to SU module '%s'.\nThis is not supported yet due to 'static variables' used in the source code of many SU modules.\n", vars->suModuleName.c_str() );
  }
  mod_su::moduleList->insertEnd( vars->suModuleName );
  //--------------------------------------------------  


  //--------------------------------------------------------------------------------
  // Prepare CSEIS -> SU header mapping, and vice versa
  //

  vars->suHdrMap      = new csSegyHdrMap(csSegyHdrMap::SEGY_SU,true);
  vars->suTrcHdrWrite = new cseis_geolib::csSegyTraceHeader(vars->suHdrMap);
  vars->suTrcHdrRead  = new cseis_geolib::csSegyTraceHeader(vars->suHdrMap);
  int numHeaders = vars->suTrcHdrWrite->numHeaders();
  vars->hdrIndexSU = new int[numHeaders];
  vars->hdrTypeSU  = new type_t[numHeaders];
  setHeaders( vars->suTrcHdrWrite, vars->suHdrMap, vars->hdrIndexSU, vars->hdrTypeSU, vars->numSamplesIn, vars->sampleIntIn, hdef, 1 );
  setHeaders( vars->suTrcHdrRead, vars->suHdrMap, vars->hdrIndexSU, vars->hdrTypeSU, shdr->numSamples, shdr->sampleInt, hdef, -1 );

  // Check if SU module has already failed
  if( vars->suPull->isError() ) {
    log->error("SU module '%s' returned an error message.", vars->suModuleName.c_str() );
  }
}

//*************************************************************************************************
// Exec phase
//

void exec_mod_su_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>(env->execPhaseDef->variables());
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;


  if( edef->isCleanup() ) {
    // Important to keep this order:
    // 1) Set EOF flag for pull and push objects
    if( vars->suPull ) {
      vars->suPull->setEOF();
    }
    if( vars->suPush ) {
      vars->suPush->setEOF();
    }
    // 2) Free thread
    if( vars->isThreadAlive ) {
      void *retval;
      pthread_join( vars->thread, &retval );
    }
    // 3) Free memory
    if( vars->suPull != NULL ) {
      delete vars->suPull;
      vars->suPull = NULL;
    }
    if( vars->suPush != NULL ) {
      delete vars->suPush;
      vars->suPush = NULL;
    }
    if( vars->suHdrMap != NULL ) {
      delete vars->suHdrMap;
      vars->suHdrMap = NULL;
    }
    if( vars->suTrcHdrRead != NULL ) {
      delete vars->suTrcHdrRead;
      vars->suTrcHdrRead = NULL;
    }
    if( vars->suTrcHdrWrite != NULL ) {
      delete vars->suTrcHdrWrite;
      vars->suTrcHdrWrite = NULL;
    }
    if( vars->hdrIndexSU != NULL ) {
      delete [] vars->hdrIndexSU;
      vars->hdrIndexSU = NULL;
    }
    if( vars->hdrTypeSU != NULL ) {
      delete [] vars->hdrTypeSU;
      vars->hdrTypeSU = NULL;
    }
    if( vars->args != NULL ) {
      delete vars->args;
      vars->args = NULL;
    }
    delete vars; vars = NULL;
    return;
  }

  // Check if SU module has already failed
  if( vars->suPull->isError() ) {
    log->error("SU module '%s' returned an error message.", vars->suModuleName.c_str() );
  }

  if( !edef->isLastCall() ) {
    // If last call has not been reached, indicate to CSEIS base system that there are potentially traces still being processed
    edef->setTracesAreWaiting();
  }

  // --------------------------------------------------
  // If there is at least one input trace, push this to the SU process
  if( traceGather->numTraces() > 0 ) {
    vars->traceCounter += 1;

    // Before pushing, make sure the 'push' object can receive a new trace (=it must be 'empty')
    while( !vars->suPush->isEmpty() ) {
      if( vars->suPull->isEmpty() ) cseis_su::csSUTraceManager::sleep_ms( (float)mod_su::PUSH_IDLE_SLEEP_TIME_US/1000.0f );
      // Meanwhile, if the 'pull' object has a trace waiting, retrieve it
      while( !vars->suPull->isEmpty() ) {
        performPull( vars, traceGather, hdef, shdr, edef->isDebug() );
      }
    }

    //**************************************************
    // Push trace to SU process
    csTrace* trace = traceGather->trace(0);
    // Prepare output trace header to be passed on to SU command
    csTraceHeader* trcHdr = trace->getTraceHeader();
    int nHeadersSU = vars->suTrcHdrWrite->numHeaders();
    for( int ihdr = 0; ihdr < nHeadersSU; ihdr++ ) {
      int hdrIdOut = vars->hdrIndexSU[ihdr];
      if( hdrIdOut < 0 ) continue;
      if( vars->hdrTypeSU[ihdr] == TYPE_FLOAT ) {
        vars->suTrcHdrWrite->setFloatValue( ihdr, trcHdr->floatValue( hdrIdOut )  );
      }
      else if( vars->hdrTypeSU[ihdr] == TYPE_DOUBLE ) {
        vars->suTrcHdrWrite->setDoubleValue( ihdr, trcHdr->doubleValue(hdrIdOut) );
      }
      else if( vars->hdrTypeSU[ihdr] == TYPE_INT ) {
        vars->suTrcHdrWrite->setIntValue( ihdr, trcHdr->intValue(hdrIdOut) );
      }
      else if( vars->hdrTypeSU[ihdr] == TYPE_INT64 ) {
        vars->suTrcHdrWrite->setIntValue( ihdr, (int)trcHdr->int64Value(hdrIdOut) );
      }
      else {
        vars->suTrcHdrWrite->setStringValue( ihdr, trcHdr->stringValue(hdrIdOut) );
      }
    }

    if( vars->suPush->putTrace( vars->suTrcHdrWrite, trace->getTraceSamples(), vars->numSamplesIn ) ) {
      traceGather->freeTrace(0);  // Remove trace
      if( edef->isDebug() ) fprintf(stdout,"CSEIS: Push trace %d\n", vars->traceCounter);
    }
    else {
      // Shouldn't happen
      fprintf(stderr,"CSEIS: Couldn't push trace %d\n", vars->traceCounter);
      log->error("CSEIS: Couldn't push trace %d\n", vars->traceCounter);
    }

  } // END: if( traceGather->numTraces() > 0 )
  //**************************************************


  // Retrieve any traces from the 'pull' object
  while( !vars->suPull->isEmpty() ) {
    performPull( vars, traceGather, hdef, shdr, edef->isDebug() );
  }

  // --------------------------------------------------
  // Last call: Pull all remaining traces from SU process
  if( edef->isLastCall() ) {
    if( edef->isDebug() ) fprintf(stdout,"CSEIS: AT EOF, last call!\n");
    vars->suPush->setEOF();  // Indicate to SU process that no more traces will arrive
    while( !vars->suPull->isEOF() ) {
      if( vars->suPull->isEmpty() ) cseis_su::csSUTraceManager::sleep_ms(1);
      performPull( vars, traceGather, hdef, shdr, edef->isDebug() );
    }
    if( edef->isDebug() ) fprintf(stdout,"CSEIS: END\n");
  }

  return;
}


//*************************************************************************************************
// Parameter definition
//
void params_mod_su_( csParamDef* pdef ) {
  pdef->setModule( "SU", "Seismic Unix module" );
  pdef->addDoc("Run SU module within Seaseis. The following SU modules are available:");

  int num10 = (int)(NUM_SU_MODULES/10) + 1;
  int counter = 0;
  for( int i = 0; i < num10; i++ ) {
    std::string text = "";
    for( int k = 0; k < 10; k++ ) {
      if( counter == NUM_SU_MODULES ) break;
      if( k > 0 ) text.append(", ");
      text.append(SU_MODULE_NAMES[counter]);
      counter += 1;
    }
    pdef->addDoc( text.c_str() );
  }

  pdef->addParam( "name", "Name of SU module to run", NUM_VALUES_FIXED, "Example: sugain" );
  pdef->addValue( "", VALTYPE_STRING );

  pdef->addParam( "param", "SU command line parameters", NUM_VALUES_FIXED, "Example (for sugain): tpow=2 epow=0.5. An arbitrary number of parameters can be specified. Specify a new line starting with 'param' as many times as needed, to specify further command line parameters" );
  pdef->addValue( "", VALTYPE_STRING );

  pdef->addParam( "sample_int", "Sample interval in output trace", NUM_VALUES_FIXED, "This parameter needs to specified only if the sample interval output by the SU module differs from the input. In this case it is required in order to inform Seaseis about the correct output sample interval." );
  pdef->addValue( "", VALTYPE_NUMBER, "Sample interval in output trace" );

  pdef->addParam( "nsamples", "Number of samples in output trace", NUM_VALUES_FIXED, "This parameter needs to specified only if the number of samples output by the SU module differs from the input. In this case it is required in order to inform Seaseis about the correct number of output samples." );
  pdef->addValue( "", VALTYPE_NUMBER, "Number of samples in output trace" );

  pdef->addParam( "ignore_su_hdr", "Ignore SU headers for number of samples and sample interval?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Number of samples and sample interval in Seaseis output trace must match output from SU" );
  pdef->addOption( "yes", "Ignore number of samples and sample interval specified in SU header" );
}

extern "C" void _params_mod_su_( csParamDef* pdef ) {
  params_mod_su_( pdef );
}
extern "C" void _init_mod_su_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_su_( param, env, log );
}
extern "C" void _exec_mod_su_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_su_( traceGather, port, numTrcToKeep, env, log );
}


void setHeaders( csSegyTraceHeader* segyTrcHdr, csSegyHdrMap* segyHdrMap, int* hdrIndexSegy, type_t* hdrTypeSegy,
                 int numSamples, float sampleInt, csTraceHeaderDef* hdef, int scalarPolarity ) {

  // Add mandatory trace headers, if defined in user specified, pre-defined, trace header map
  int num_id = 0;
  int id_nsamp        = num_id++;
  int id_sampint_us   = num_id++;
  int id_scalar_coord = num_id++;
  int id_scalar_elev  = num_id++;
  int id_scalar_stat  = num_id++;
  int id_fold         = num_id++;
  int id_fold_vert    = num_id++;
  int id_data_type    = num_id++;
  int id_trc_type     = num_id++;
  int id_unit_coord   = num_id++;
  int id_time_code    = num_id++;;
  csVector<string> mandatoryHdrNames(num_id);
  csVector<int>    mandatoryHdrValues(num_id);
  mandatoryHdrNames.allocate( num_id, "" );
  mandatoryHdrValues.allocate( num_id, 0 );
  // Set default values for mandatory SEGY output headers:
  mandatoryHdrNames.set("nsamp",       id_nsamp);        mandatoryHdrValues.set( numSamples, id_nsamp );
  mandatoryHdrNames.set("sampint_us",  id_sampint_us);   mandatoryHdrValues.set( (int)(sampleInt*1000), id_sampint_us );
  mandatoryHdrNames.set("scalar_coord", id_scalar_coord); mandatoryHdrValues.set( scalarPolarity * 100, id_scalar_coord );
  mandatoryHdrNames.set("scalar_elev", id_scalar_elev);  mandatoryHdrValues.set( scalarPolarity * 10, id_scalar_elev );
  mandatoryHdrNames.set("scalar_stat", id_scalar_stat);  mandatoryHdrValues.set( scalarPolarity * 1, id_scalar_stat );
  mandatoryHdrNames.set("fold",        id_fold);         mandatoryHdrValues.set( 1, id_fold );
  mandatoryHdrNames.set("fold_vert",   id_fold_vert);    mandatoryHdrValues.set( 1, id_fold_vert );
  mandatoryHdrNames.set("data_type",   id_data_type);    mandatoryHdrValues.set( 1, id_data_type );
  mandatoryHdrNames.set("trc_type",    id_trc_type);     mandatoryHdrValues.set( 1, id_trc_type );
  mandatoryHdrNames.set("unit_coord",  id_unit_coord);   mandatoryHdrValues.set( 1, id_unit_coord );
  mandatoryHdrNames.set("time_code",   id_time_code);    mandatoryHdrValues.set( 4, id_time_code );
  // Check if mandatory headers are defined in user selected pre-set header map
  // Each header that is not defined will be removed from the list of mandatory headers
  int hdrIndex = -1;
  for( int ihdr = mandatoryHdrNames.size()-1; ihdr >= 0; ihdr-- ) {
    string name = mandatoryHdrNames.at(ihdr);
    if( !segyHdrMap->contains(name,&hdrIndex) ) {
      mandatoryHdrNames.remove( ihdr );
      mandatoryHdrValues.remove( ihdr );
    }
  }

  int numHeaders = segyTrcHdr->numHeaders();
  for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
    string const& name = segyTrcHdr->headerName(ihdr);
    //    log->line(" #%2d  %s", ihdr+1, name.c_str());
    if( hdef->headerExists(name) ) {
      hdrIndexSegy[ihdr] = hdef->headerIndex( name );
      hdrTypeSegy[ihdr]  = hdef->headerType( name );
    }
    else {  // Must be mandatory header...
      hdrIndexSegy[ihdr] = -1;
    }
  }

  for( int ihdr = 0; ihdr < mandatoryHdrNames.size(); ihdr++ ) {
    string name = mandatoryHdrNames.at(ihdr);
    int hdrIndex = segyTrcHdr->headerIndex(name);
    //    if( hdrIndex < 0 ) log->error("Program bug! Wrong mandatory header index....");
    if( hdef->headerExists(name) ) {
      hdrIndexSegy[hdrIndex] = -1;
    }
    segyTrcHdr->setIntValue( hdrIndex, mandatoryHdrValues.at(ihdr) );
  }

}

//--------------------------------------------------------------------------------
//
void performPull( mod_su::VariableStruct* vars,
                  csTraceGather* traceGather,
                  csTraceHeaderDef const* hdef,
                  csSuperHeader const* shdr,
                  bool isDebug )
{
  unsigned char const* bufferPtr = NULL;
  if( vars->suPull->getTracePtr( &bufferPtr ) ) {
    // First, check consistency of sample interval and number of samples in trace from SU:
    if(  !vars->ignoreSUHdr ) {
      if( vars->suPull->numSamples() != shdr->numSamples ) {
        throw( csException("Number of samples output by SU module '%s' (=%d) does not match number of samples specified for SeaSeis (=%d). Specify ' num_samples  %d' in Seaseis flow",
                           vars->suModuleName.c_str(), vars->suPull->numSamples(), shdr->numSamples, vars->suPull->numSamples())  );
      }
      if( fabs(vars->suPull->sampleInt()-shdr->sampleInt) > 0.0001 ) {
        throw( csException("Sample interval output by SU module '%s' (=%.4fms) does not match sample interval specified for SeaSeis (=%.4fms). Specify ' sample_int  %f' or ' ignore_su_hdr yes' in Seaseis flow",
                           vars->suModuleName.c_str(), vars->suPull->sampleInt(), shdr->sampleInt, vars->suPull->sampleInt())  );
      }
    }

    csTrace* newTrace = traceGather->createTrace( hdef, shdr->numSamples );
    csTraceHeader* newTrcHdr = newTrace->getTraceHeader();
    //    fprintf(stdout,"CSEIS: Successful pull of trace %d   %x\n", vars->traceCounter, bufferPtr );
    //  fflush(stdout);
    float* newSamples = newTrace->getTraceSamples();

    memcpy(newSamples,&bufferPtr[240],shdr->numSamples*sizeof(float));

    vars->suTrcHdrRead->readHeaderValues( bufferPtr, false, true );  // swapEndian=false, isAutoScaleHeaders=true
    vars->suPull->freeTrace(); // Crucial to free trace after calling getTracePtr()
          
    int nHeaders = vars->suTrcHdrRead->numHeaders();
    for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
      int hdrIdOut = vars->hdrIndexSU[ihdr];
      if( hdrIdOut < 0 ) continue;
      switch( vars->hdrTypeSU[ihdr] ) {
      case TYPE_FLOAT:
        newTrcHdr->setFloatValue( hdrIdOut, vars->suTrcHdrRead->floatValue(ihdr) );
        break;
      case TYPE_DOUBLE:
        newTrcHdr->setDoubleValue( hdrIdOut, vars->suTrcHdrRead->doubleValue(ihdr) );
        break;
      case TYPE_INT:
        newTrcHdr->setIntValue( hdrIdOut, vars->suTrcHdrRead->intValue(ihdr) );
        break;
      case TYPE_INT64:
        newTrcHdr->setInt64Value( hdrIdOut, (csInt64_t)vars->suTrcHdrRead->intValue(ihdr) );
        break;
      case TYPE_STRING:
        newTrcHdr->setStringValue( hdrIdOut, vars->suTrcHdrRead->stringValue(ihdr) );
        break;
      }
    } // END for ihdr
    
    if( isDebug ) fprintf(stdout,"CSEIS: Successful pull of trace %d\n", vars->traceCounter);
  }
  else {
    if( isDebug ) fprintf(stdout,"CSEIS: Couldn't pull trace %d\n", vars->traceCounter);
  }
}

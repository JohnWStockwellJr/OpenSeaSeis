/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csModule.h"
#include "cseis_includes.h"

#include "csTraceGather.h"
#include "csTraceHeaderDef.h"
#include "csTrace.h"
#include "csTraceHeader.h"
#include "csUserParam.h"
#include "csParamDef.h"
#include "csMemoryPoolManager.h"
#include "csMethodRetriever.h"
#include "csExecPhaseDef.h"
#include "csInitExecEnv.h"

#include "csException.h"
#include "csVector.h"
#include "csQueue.h"
#include "csFlexNumber.h"
#include "csTimer.h"
#include "csTable.h"
#include <cstdlib>
#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace cseis_system;

csModule::csModule( std::string const& name, int id_number, csMemoryPoolManager* memManager, int mpiProcID, int mpiNumProc ) : myUniqueID(id_number) {
  myMemoryPoolManager = memManager;
  myHeaderDef = NULL;
  myName = name;
  myMPIProcID  = mpiProcID;
  myMPINumProc = mpiNumProc;
  init();
}
void csModule::init() {
  myExecPhaseDef = new csExecPhaseDef( myName, myMPIProcID );
  myTraceGather  = new csTraceGather( myMemoryPoolManager );
  myTraceQueue   = new cseis_geolib::csQueue<csTrace*>();

  myMethodParam     = NULL;
  myMethodInit      = NULL;
  myMethodExecStart = NULL;
  myMethodExec      = NULL;
  myMethodCleanup   = NULL;

  myNumFixedTracesToCollect = 0;
  myNumTracesToBePassed     = 0;
  myTotalNumProcessedTraces = 0;
  myTotalNumIncomingTraces  = 0;
  myExecEnvPtr              = NULL;
  mySuperHeader             = new csSuperHeader();

  myEnsembleKeyHeaderIndex = NULL;
  myEnsembleKeyValue       = NULL;
  myNextEnsembleKeyValue   = NULL;
  myHelperHdrValues = NULL;
  myIsEnsembleFull = false;

  myNumOutputPorts         = 1;
  myNumInputPorts          = 1;
  myTimeExecPhaseCPU       = 0.0;

  myIsFinishedProcessing = false;
  
  myVersion[MAJOR] = 1;
  myVersion[MINOR] = 0;

  // Special modules
  if( !myName.compare("IF") ) {
    myModuleType = MODTYPE_IF;
    myNumOutputPorts = 2;
  }
  else if( !myName.compare("ELSEIF") ) {
    myModuleType = MODTYPE_ELSEIF;
    myNumOutputPorts = 2;
  }
  else if( !myName.compare("ELSE") ) {
    myModuleType = MODTYPE_ELSE;
  }
  else if( !myName.compare("ENDIF") ) {
    myModuleType = MODTYPE_ENDIF;
  }
  else if( !myName.compare("SPLIT") ) {
    myModuleType = MODTYPE_SPLIT;
    myNumOutputPorts = 2;
  }
  else if( !myName.compare("ENDSPLIT") ) {
    myModuleType = MODTYPE_ENDSPLIT;
  }
  else {
    myModuleType = MODTYPE_UNKNOWN;
  }
}
//---------------------------------------------------------------------
//
csModule::~csModule() {
  if( myExecPhaseDef != NULL ) {
    delete myExecPhaseDef;
    myExecPhaseDef = NULL;
  }
  if( myExecEnvPtr != NULL ) {
    delete myExecEnvPtr;
    myExecEnvPtr = NULL;
  }
  if( myHeaderDef != NULL ) {
    delete myHeaderDef;
    myHeaderDef = NULL;
  }
  if( myTraceGather != NULL ) {
    delete myTraceGather;
    myTraceGather = NULL;
  }
  if( myTraceQueue != NULL ) {
    delete myTraceQueue;
    myTraceQueue = NULL;
  }
  if( mySuperHeader ) {
    delete mySuperHeader;
    mySuperHeader = NULL;
  }
  if( myEnsembleKeyHeaderIndex ) {
    delete [] myEnsembleKeyHeaderIndex;
    myEnsembleKeyHeaderIndex = NULL;
  }
  if( myEnsembleKeyValue ) {
    delete [] myEnsembleKeyValue;
    myEnsembleKeyValue = NULL;
  }
  if( myNextEnsembleKeyValue ) {
    delete [] myNextEnsembleKeyValue;
    myNextEnsembleKeyValue = NULL;
  }
  if( myHelperHdrValues != NULL ) {
    delete [] myHelperHdrValues; myHelperHdrValues = NULL;
  }
}
//*********************************************************************
//
// Module submission methods
//
//*********************************************************************

void csModule::retrieveParamInitMethods() {
  csMethodRetriever::getParamInitMethod( myName, myVersion[MAJOR], myVersion[MINOR], myMethodParam, myMethodInit );
  if( myMethodParam == NULL ) {
    throw( cseis_geolib::csException("Program bug: Module's Parameter definition method not found.") );
  }
  else if( myMethodInit == NULL ) {
    throw( cseis_geolib::csException("Program bug: Module's Init phase method not found.") );
  }
}

void csModule::getParamDef( csParamDef* paramDef ) {
  if( myMethodParam == NULL ) {
    retrieveParamInitMethods();
  }
  (*myMethodParam)( paramDef );
}
//-------------------------------------------------------------------
void csModule::submitInitPhase( csParamManager* paramManager, csLogWriter* writer, cseis_geolib::csTable const** tables, int numTables ) {
  if( myHeaderDef == NULL ) throw("csModule::submitInitPhase: Program bug: header definition object has not been initialized yet.");
  csInitPhaseEnv initEnv( myHeaderDef, myExecPhaseDef, mySuperHeader, tables, numTables );
  
  if( myMethodInit == NULL ) {
    retrieveParamInitMethods();
  }

  // Run init phase:
  (*myMethodInit)( paramManager, &initEnv, writer );

  // Set exec phase method pointer
  if( myExecPhaseDef->traceMode == TRCMODE_UNKNOWN ) {
    throw( cseis_geolib::csException("Program bug: Trace mode not set in module init method. Should be set to FIXED, ENSEMBLE...") );
  }

  csMethodRetriever::getExecMethod( myName, myVersion[MAJOR], myVersion[MINOR], myMethodExecStart, myMethodExec, myMethodCleanup );
  if( myMethodExecStart == NULL ) {
    throw( cseis_geolib::csException("Program bug: Module's StartExec phase method not found.") );
  }
  if( myMethodExec == NULL ) {
    throw( cseis_geolib::csException("Program bug: Module's Exec phase method not found.") );
  }
  if( myMethodCleanup == NULL ) {
    throw( cseis_geolib::csException("Program bug: Module's Clean-up phase method not found.") );
  }

  if( initEnv.errorCount() > 0 ) {
    throw( cseis_geolib::csException("Error(s) occurred in init phase. See log file for details.") );
  }
  //--------------------------------------------------------------------------------
  // Consistency check of certain settings (basically checking that module init phase module is correct)
  //
  // 1) Check that all ensemble key headers exist. If not, remove this key from superheader, issue warning message
  int keyIndex = 0;
  while( keyIndex < mySuperHeader->numEnsembleKeys() ) {
    std::string const* name = mySuperHeader->ensembleKey( keyIndex );
    if( myHeaderDef->headerExists( name->c_str() ) ) {
      keyIndex += 1;
    }
    else {
      writer->warning("Super header: Ensemble key header '%s' does not exist in trace header. Ensemble key has been removed from super header.", name->c_str());
      mySuperHeader->removeEnsembleKey( keyIndex );
    }
  }
  // Store trace header index of ensemble keys, for quick access
  int numKeys = mySuperHeader->numEnsembleKeys();
  if( numKeys > 0 ) {
    myEnsembleKeyHeaderIndex = new int [numKeys];
    myEnsembleKeyValue       = new cseis_geolib::csFlexNumber[numKeys];
    myNextEnsembleKeyValue   = new cseis_geolib::csFlexNumber[numKeys];
    for( int ikey = 0; ikey < numKeys; ikey++ ) {
      myEnsembleKeyHeaderIndex[ikey] = myHeaderDef->headerIndex( *mySuperHeader->ensembleKey(ikey) );
    }
    myHelperHdrValues = new cseis_geolib::csFlexNumber[numKeys];
  }
  myHeaderDef->resetByteLocation();

  if( myExecPhaseDef->traceMode == TRCMODE_FIXED ) {
    myNumFixedTracesToCollect = myExecPhaseDef->numTraces;
#ifdef USE_MPI
    if( isMPISupported() ) myNumFixedTracesToCollect *= myMPINumProc;
#endif
  }
}

bool csModule::isMPISupported() const {
  return myExecPhaseDef->isMPISupported();
}
//-------------------------------------------------------------------
//
//
//
bool csModule::submitExecStartPhase( csLogWriter* writer ) {
  if( !myMethodExecStart ) {
    writer->line( "Cannot find Start Exec Phase function for module '%s', required since v3.00.", getName() );
    return false;
  }
  if( myMPIProcID != 0 && !isMPISupported() ) return true; // MPI not supported AND not the main process  --> Skip start exec, pretend it was OK. StartExec phase is only run for main process
  
  // Run Start Exec Phase
  return( (*myMethodExecStart)( myExecEnvPtr, writer ) );
}
//-------------------------------------------------------------------
//
//
//
bool csModule::isReadyToSubmitExecPhase( bool forceToRun ) const {
  if( myMPIProcID != 0 ) {
    throw( cseis_geolib::csException("Program bug in csModule::isReadyToExecPhase(): Should only be called for main MPI process\n") );
  }
  if( myNumTracesToBePassed != 0 ) {
    throw( cseis_geolib::csException("csModule::isReadyToSubmitExecPhase(): Program bug? Processed traces in trace gather have not been passed yet.." ) );
  }
  //---------------------------------------------------------------------------
  if( myExecPhaseDef->execType() != EXEC_TYPE_INPUT ) {  // Any module that is not an INPUT module 
    int totalNumTraces = myTraceGather->numTraces()+myTraceQueue->size();
    if( myExecPhaseDef->traceMode == TRCMODE_FIXED ) {
      if( totalNumTraces >= myNumFixedTracesToCollect ) return true;
    }
    else { // if( myExecPhaseDef->traceMode == TRCMODE_ENSEMBLE ) {
      if( mySuperHeader->numEnsembleKeys() > 0 && myIsEnsembleFull ) return true;
      // else: No ensemble key set --> process entire data set only when forced (when all traces have been read in)
    }
    return( forceToRun && (totalNumTraces > 0 || myExecPhaseDef->tracesAreWaiting()) );
  }
  //---------------------------------------------------------------------------
  else { // if( myExecPhaseDef->execType() == EXEC_TYPE_INPUT ) {  // Input module
    return true;
  }
  return false;
}

//********************************************************************************
//
bool csModule::submitExecPhase( bool forceToProcess, csLogWriter* writer, int& outPort ) {
  if( myMPIProcID != 0 ) {
    throw( cseis_geolib::csException("Program bug in csModule::submitExecPhase(): Should only be called for main MPI process\n") );
  }
  cseis_geolib::csTimer timer;
  timer.start();
  int nProcessedTraces = 0;
  myExecPhaseDef->myIsLastCall = forceToProcess;
  
  // Case (1) Input module
  if( myExecPhaseDef->execType() == EXEC_TYPE_INPUT ) {
    // Input: Retrieve/allocate new trace, then submit exec phase to read in this one trace
    csTrace* trace = myMemoryPoolManager->getNewTrace( myHeaderDef, mySuperHeader->numSamples );
    myTraceGather->addTrace( trace );
    outPort = 0;  // Output port should be 0, because input module only has one output port...
    int numTrcToKeep = 0;
    (*myMethodExec)( myTraceGather, &outPort, &numTrcToKeep, myExecEnvPtr, writer );
 
    if( myTraceGather->numTraces() > 0 ) {
      nProcessedTraces = myTraceGather->numTraces(); // Should be 1 trace
    }
    else {   // No trace has been read in
      myTraceGather->freeAllTraces();
    }
    myIsFinishedProcessing = true;
  }
  //----------------------------------------------------------------------------------------
  // Case (2) Not an Input module
  else {
    // MPIDebug writer->line("csModule:  Submit fixed trace module with %d traces. %d MPI processes\n", myExecPhaseDef->numTraces, myMPINumProc);
    if( myExecPhaseDef->traceMode == TRCMODE_FIXED ) {
      while( myTraceGather->numTraces() < myNumFixedTracesToCollect && !myTraceQueue->isEmpty() ) {
        myTraceGather->addTrace( myTraceQueue->pop() );
      }
    }
    myExecPhaseDef->myIsLastCall = (forceToProcess && myTraceQueue->size() == 0); // Do not force to run when there are still traces in "queue"    
    myTotalNumIncomingTraces += myTraceGather->numTraces();  // Accumulate number of traces input to this module so far

    int numTrcToKeep = 0;  // Number of traces to keep, i.e. traces to roll into next pass of this exec method
    myExecPhaseDef->myTracesAreWaiting = false;

#ifdef USE_MPI
    if( isMPISupported() ) mpi_distribute_traces_before_exec_phase( writer );
#endif
    //********************************************************************************
    (*myMethodExec)( myTraceGather, &outPort, &numTrcToKeep, myExecEnvPtr, writer );
    //********************************************************************************
#ifdef USE_MPI
    if( isMPISupported() ) {
      if( numTrcToKeep != 0 || myExecPhaseDef->myTracesAreWaiting ) {
        throw( cseis_geolib::csException("Program bug in module %s: \nMPI is supported but module exec phase routine\n (a) ...asks to keep 'numTrcToKeep'(=%d) traces, or (b) ...reports that traces are 'waiting' (=%s)\n", myName.c_str(), numTrcToKeep, myExecPhaseDef->myTracesAreWaiting ? "true" : "false") );
      }
      mpi_distribute_traces_after_exec_phase( writer );
    }
#endif
    myIsEnsembleFull = false;
    myTotalNumIncomingTraces -= numTrcToKeep;  // Correct number of traces input to this module by traces that are kept

    if( numTrcToKeep > myTraceGather->numTraces() ) {
      writer->warning("Module %s: Supposed number of traces to keep is larger than number of traces in gather. This is probably due to a program bug in the module method.", myName.c_str());
      numTrcToKeep = myTraceGather->numTraces();
    }
    else if( myExecPhaseDef->myIsLastCall && numTrcToKeep ) {
      // Do nothing. Allow modules to continue processing even if forced to run. This may lead to infinite loops if modules are not handled correctly
    }

    if( myExecPhaseDef->traceMode == TRCMODE_FIXED ) {
      nProcessedTraces = myTraceGather->numTraces() - numTrcToKeep;
      myIsFinishedProcessing = true;
      int totalNumTracesLeft = numTrcToKeep + myTraceQueue->size();

      if( totalNumTracesLeft >= myNumFixedTracesToCollect || (forceToProcess && totalNumTracesLeft != 0) ) {
        myIsFinishedProcessing = false;
      }
    }
    else { //if( myExecPhaseDef->traceMode == TRCMODE_ENSEMBLE ) {  // Case 2
      nProcessedTraces = myTraceGather->numTraces() - numTrcToKeep;
      if( myTraceGather->numTraces() == 0 ) { // Case 2.1
        updateTracesEnsembleModule();
      }
      else if( numTrcToKeep == 0 ) {  // Case 2.2
        updateTracesEnsembleModule();
      }
      else { // if( nProcessedTraces == 0 ) { // Case 2.3 & 2.4
        myIsEnsembleFull = true;
        myIsFinishedProcessing = false;
      }
      if( forceToProcess && (myTraceGather->numTraces()-nProcessedTraces + myTraceQueue->size()) != 0 ) {
        myIsEnsembleFull = true;
        myIsFinishedProcessing = false;
      }
    }

  }

  // ******** DELETE HEADERS (if any) ON TRACES OUTPUT TO NEXT MODULE *********
  //
  if( myHeaderDef->getIndexOfHeadersToDel()->size() > 0 ) {
    for( int itrc = 0; itrc < nProcessedTraces; itrc++ ) {
      myTraceGather->trace(itrc)->getTraceHeader()->deleteHeaders( myHeaderDef );
    }
  }
  // ******** SET NUMBER OF SAMPLES FOR ALL TRACES MOVED TO FURTHER MODULES ********
  for( int itrc = 0; itrc < nProcessedTraces; itrc++ ) {
    myTraceGather->trace(itrc)->getTraceDataObject()->set( mySuperHeader->numSamples );
  }
  // **************************************************************************
  if( outPort > myNumOutputPorts ) {
    throw( cseis_geolib::csException("csModule::submitExecPhase(): Output port number too large. This is most likely due to a program bug in the module method...") );
  }
  myNumTracesToBePassed     += nProcessedTraces; // Add processed traces to number of traces to be passed to next module
  myTotalNumProcessedTraces += nProcessedTraces; // Accumulate number of traces processed by this module
  myTimeExecPhaseCPU += timer.getElapsedTime();  // Accumulate exec phase CPU time
  //  fprintf(stderr,"MPI proc %d: END of csModule %s submit %d %d\n", myMPIProcID, getName(), nProcessedTraces, myNumTracesToBePassed);
  return( nProcessedTraces > 0 );
}
//-------------------------------------------------------------------
//
//
void csModule::submitCleanupPhase( csLogWriter* writer, bool reportMissing ) {
  myExecPhaseDef->myIsCleanup = true;
  if( myMethodCleanup ) {
    (*myMethodCleanup)( myExecEnvPtr, writer );
  }
  else if( reportMissing ) {
    writer->warning( "Cannot find Cleanup Phase function for module '%s', required since v3.00.", getName() );
  }
}

//-------------------------------------------------------------------
//
#ifdef USE_MPI
// START of MPI code
//
//
bool csModule::mpi_submitExecPhase( csLogWriter* writer, int& outPort ) {
  if( myMPIProcID == 0 || !isMPISupported() ) {
    throw( cseis_geolib::csException("Program bug in csModule::mpi_submitExecPhase(): Should only be called for non-main MPI processes for modules with MPI suport\n") );
  }
  cseis_geolib::csTimer timer;
  timer.start();

  int numTrcToKeep = 0;  // Number of traces to keep, i.e. traces to roll into next pass of this exec method
  myExecPhaseDef->myTracesAreWaiting = false;

  mpi_distribute_traces_before_exec_phase( writer );
  (*myMethodExec)( myTraceGather, &outPort, &numTrcToKeep, myExecEnvPtr, writer );
  if( numTrcToKeep != 0 || myExecPhaseDef->myTracesAreWaiting ) {
    throw( cseis_geolib::csException("Program bug in module %s: \nMPI is supposedly supported but module exec phase asks\n (a) to keep 'numTrcToKeep'(=%d) traces...or (b) that traces are 'waiting' (=%s)\n", myName.c_str(), numTrcToKeep, myExecPhaseDef->myTracesAreWaiting ? "true" : "false") );
  }
  mpi_distribute_traces_after_exec_phase( writer );
  
  myTimeExecPhaseCPU += timer.getElapsedTime();  // Accumulate exec phase CPU time for other MPI processes
  return true;
}
void csModule::mpi_decompressInitPhaseObjects( char const* dataSuperHdr, char const* dataHdrDef ) {
  mySuperHeader->mpi_decompress( dataSuperHdr );
  myHeaderDef->mpi_decompress( dataHdrDef );
}
void csModule::mpi_distribute_traces_before_exec_phase( csLogWriter* writer ) {
  if( !isMPISupported() || myMPINumProc == 1 ) return;
  int numTracesAll = myTraceGather->numTraces();

  //  MPIDebug  writer->line("csModule::mpi_distribute_before(): MPI %d: csModule %s, distribute %d traces\n", myMPIProcID, getName(), numTracesAll);
  
  int tag = 99;
  MPI_Status stat;
  if( myMPIProcID == 0 ) {
    int numSamples       = myTraceGather->trace(0)->numSamples();
    int numBytesHdrBlock = myHeaderDef->getTotalNumBytes();
    int numTracesProc0   = std::min( myExecPhaseDef->numTraces, numTracesAll );
    int traceCounter     = numTracesProc0;
    // MPIDebug writer->line("mpi_distribute_before: MPI %d: csModule %s, keep %d traces & %d blocks", myMPIProcID, getName(), numTracesProc0, numBytesHdrBlock);
    for( int mpiTargetProc = 1; mpiTargetProc < myMPINumProc; mpiTargetProc++ ) {
      int numTracesToSend = std::min( myExecPhaseDef->numTraces, numTracesAll-traceCounter );
      MPI_Send( &numTracesToSend, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
      // MPIDebug writer->line("mpi_distribute_before: MPI %d to %d: csModule %s, sent %d traces (ntrc_gather: %d, ntrc_queue: %d)", myMPIProcID, mpiTargetProc, getName(), numTracesToSend, myTraceGather->numTraces(), myTraceQueue->size() );
      if( numTracesToSend > 0 ) {
        MPI_Send( &numSamples, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
        MPI_Send( &numBytesHdrBlock, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
        // MPIDebug writer->line("mpi_distribute_before: MPI %d: csModule %s, sent %d numSamples, %d numBytes", myMPIProcID, getName(), numSamples, numBytesHdrBlock);
        for( int itrc = 0; itrc < numTracesToSend; itrc++ ) {
          int traceIndex = itrc + traceCounter;
          csTrace* trace = myTraceGather->trace( traceIndex );

          float* dataSamples = trace->getTraceSamples();
          char* hdrBlockHandle = trace->getTraceHeader()->getTraceHeaderValueBlockHandle();
          MPI_Send( dataSamples, numSamples, MPI::FLOAT, mpiTargetProc, tag, MPI_COMM_WORLD );
          MPI_Send( hdrBlockHandle, numBytesHdrBlock, MPI::CHAR, mpiTargetProc, tag, MPI_COMM_WORLD );
          // MPIDebug writer->line("mpi_distribute_before: MPI %d to %d: csModule %s, sent trace index %d", myMPIProcID, mpiTargetProc, getName(), traceIndex );
        }
        traceCounter += numTracesToSend;
      }
    }
    // MPIDebug writer->line("mpi_distribute_before: MPI %d: csModule %s, Completed sending traces", myMPIProcID, getName() );
    myTraceGather->reduceNumTracesTo( numTracesProc0 );
  }
  else { // ...not the main MPI process
    int numTracesToReceive;
    MPI_Recv( &numTracesToReceive, 1, MPI::INT, 0, tag, MPI_COMM_WORLD, &stat );
    // MPIDebug writer->line("mpi_distribute_before: MPI %d: csModule %s, receive %d traces\n", myMPIProcID, getName(), numTracesToReceive );
    if( numTracesToReceive > 0 ) {
      int numSamples;
      int numBytesHdrBlock;
      MPI_Recv( &numSamples, 1, MPI::INT, 0, tag, MPI_COMM_WORLD, &stat );
      MPI_Recv( &numBytesHdrBlock, 1, MPI::INT, 0, tag, MPI_COMM_WORLD, &stat );
      // MPIDebug writer->line("mpi_distribute_before: MPI %d: csModule %s, received %d numSamples, %d numBytes", myMPIProcID, getName(), numSamples, numBytesHdrBlock);
      myTraceGather->createTraces( 0, numTracesToReceive, myHeaderDef, numSamples );
      for( int itrc = 0; itrc < numTracesToReceive; itrc++ ) {
        csTrace* trace = myTraceGather->trace( itrc );
        float* dataSamples = trace->getTraceSamples();
        char* hdrBlockHandle = trace->getTraceHeader()->getTraceHeaderValueBlockHandle();
        MPI_Recv( dataSamples, numSamples, MPI::FLOAT, 0, tag, MPI_COMM_WORLD, &stat );
        MPI_Recv( hdrBlockHandle, numBytesHdrBlock, MPI::CHAR, 0, tag, MPI_COMM_WORLD, &stat );
        // MPIDebug writer->line("mpi_distribute_before: MPI %d: csModule %s, received trace index %d", myMPIProcID, getName(), itrc );
      }
    }
  }
  //  MPIDebug writer->line("mpi_distribute_before(): MPI %d: csModule %s, finish distributing traces. Num traces now: %d\n", myMPIProcID, getName(), myTraceGather->numTraces());
}
//------------------------------------------------------------------------------------------------------------------------
// Distribute traces from MPI processes back to main process
//
void csModule::mpi_distribute_traces_after_exec_phase( csLogWriter* writer ) {
  if( !isMPISupported() || myMPINumProc == 1 ) return;
  int numTracesIn = myTraceGather->numTraces();
  // MPIDEBUG  fprintf(stderr,"MPI %d: csModule %s, back-distribute %d traces\n", myMPIProcID, getName(), numTracesIn);
  int tag = 99;
  MPI_Status stat;
  if( myMPIProcID != 0 ) {
    int mpiTargetProc = 0;
    int numTracesToSend = numTracesIn;
    MPI_Send( &numTracesToSend, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
    if( numTracesToSend > 0 ) {
      int numSamples       = myTraceGather->trace(0)->numSamples();
      int numBytesHdrBlock = myHeaderDef->getTotalNumBytes();
      // MPIDEBUG      fprintf(stderr,"MPI %d: csModule %s, send %d traces & %d blocks\n", myMPIProcID, getName(), numTracesToSend, numBytesHdrBlock);
      MPI_Send( &numSamples, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
      MPI_Send( &numBytesHdrBlock, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
      for( int itrc = 0; itrc < numTracesToSend; itrc++ ) {
        csTrace* trace = myTraceGather->trace( itrc );
        float* dataSamples = trace->getTraceSamples();
        char* hdrBlockHandle = trace->getTraceHeader()->getTraceHeaderValueBlockHandle();
        MPI_Send( dataSamples, numSamples, MPI::FLOAT, mpiTargetProc, tag, MPI_COMM_WORLD );
        MPI_Send( hdrBlockHandle, numBytesHdrBlock, MPI::CHAR, mpiTargetProc, tag, MPI_COMM_WORLD );
      }
    }
    myTraceGather->freeAllTraces();
  }
  else { // ...This is the main MPI process
    for( int mpiTargetProc = 1; mpiTargetProc < myMPINumProc; mpiTargetProc++ ) {
      // MPIDEBUG      fprintf(stderr,"MPI %d: csModule %s, back-distribute from mpi process #%d. Number of traces in gather: %d\n", myMPIProcID, getName(), mpiTargetProc, myTraceGather->numTraces() );
      int numTracesToReceive;
      MPI_Recv( &numTracesToReceive, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD, &stat );
      if( numTracesToReceive > 0 ) {
        int numSamples;
        int numBytesHdrBlock;
        MPI_Recv( &numSamples, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD, &stat );
        MPI_Recv( &numBytesHdrBlock, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD, &stat );
        // MPIDEBUG     fprintf(stderr,"MPI %d: csModule %s, receive %d traces & %d blocks\n", myMPIProcID, getName(), numTracesToReceive, numBytesHdrBlock);
        int numTracesCurrent = myTraceGather->numTraces();
        myTraceGather->createTraces( numTracesCurrent, numTracesToReceive, myHeaderDef, numSamples );
        for( int itrc = 0; itrc < numTracesToReceive; itrc++ ) {
          csTrace* trace = myTraceGather->trace( itrc + numTracesCurrent );
          float* dataSamples = trace->getTraceSamples();
          char* hdrBlockHandle = trace->getTraceHeader()->getTraceHeaderValueBlockHandle();
          MPI_Recv( dataSamples, numSamples, MPI::FLOAT, mpiTargetProc, tag, MPI_COMM_WORLD, &stat );
          MPI_Recv( hdrBlockHandle, numBytesHdrBlock, MPI::CHAR, mpiTargetProc, tag, MPI_COMM_WORLD, &stat );
          // MPIDEBUG       fprintf(stderr,"MPI %d: csModule %s, back-distribute\n", myMPIProcID, getName() );
        }
      }
    }
    // MPIDEBUG    fprintf(stderr,"MPI %d: csModule %s, back-distribute finished, received all traces, number %d\n", myMPIProcID, getName(), myTraceGather->numTraces() );
  }
  // MPIDebug writer->line("mpi_distribute_after(): MPI %d: csModule %s, finish distributing traces. Num traces now: %d\n", myMPIProcID, getName(), myTraceGather->numTraces());
}
// END of MPI code
#endif

//*********************************************************************
//
// Trace management methods
//
//*********************************************************************

//------------------------------------------------------
void csModule::addNewTraceToGather( csTrace* trace, int inPort ) {
  trace->getTraceHeader()->setHeaders( myHeaderDef, inPort );
  trace->getTraceDataObject()->setMax( mySuperHeader->numSamples );
  myTraceGather->addTrace( trace );
}
void csModule::addNewTraceToQueue( csTrace* trace, int inPort ) {
  trace->getTraceHeader()->setHeaders( myHeaderDef, inPort );
  trace->getTraceDataObject()->setMax( mySuperHeader->numSamples );
  myTraceQueue->push( trace );
}
//------------------------------------------------------
//
bool csModule::traceIsPartOfCurrentEnsemble( csTrace* trace ) {
  csTraceHeader* trcHdrPtr = trace->getTraceHeader();
  int numKeys = mySuperHeader->numEnsembleKeys();
  for( int ikey = 0; ikey < numKeys; ikey++ ) {
    int hdrIndex = myEnsembleKeyHeaderIndex[ikey];
    char hdrType = myHeaderDef->headerType( hdrIndex );
    if( hdrType == cseis_geolib::TYPE_INT ) {
      myHelperHdrValues[ikey].setIntValue( trcHdrPtr->intValue(hdrIndex) );
    }
    else if( hdrType == cseis_geolib::TYPE_FLOAT ) {
      myHelperHdrValues[ikey].setFloatValue( trcHdrPtr->floatValue(hdrIndex) );
    }
    else if( hdrType == cseis_geolib::TYPE_DOUBLE ) {
      myHelperHdrValues[ikey].setDoubleValue( trcHdrPtr->doubleValue(hdrIndex) );
    }
    else if( hdrType == cseis_geolib::TYPE_STRING ) {
      throw( cseis_geolib::csException("Encountered ensemble key of type string. This is currently not supported.") );
    }
  }

  if( myTraceGather->numTraces() > 0 ) {  // Current module already has traces in gather -> Check if ensemble header values agree
    for( int ikey = 0; ikey < numKeys; ikey++ ) {
      if( myEnsembleKeyValue[ikey] != myHelperHdrValues[ikey] ) {
        for( int ikeyNew = 0; ikeyNew < numKeys; ikeyNew++ ) {
          myNextEnsembleKeyValue[ikeyNew] = myHelperHdrValues[ikeyNew];
        }
        myIsEnsembleFull = true;
        return false;
      }
    }
  }
  else {  // Current module does not have any traces yet in gather -> set ensemble keys
    for( int ikey = 0; ikey < numKeys; ikey++ ) {
      myEnsembleKeyValue[ikey] = myHelperHdrValues[ikey];
    }
  }
  return true;
}

//------------------------------------------------------------------
//
void csModule::updateTracesEnsembleModule() {
  for( int i = 0; i < mySuperHeader->numEnsembleKeys(); i++ ) {
    myEnsembleKeyValue[i] = myNextEnsembleKeyValue[i];
  }
  myIsEnsembleFull = false;
  myIsFinishedProcessing = true;
  if( !myTraceQueue->isEmpty() ) {
    // Move first trace from queue to gather, this doesn't need check
    myTraceGather->addTrace( myTraceQueue->pop() );
    while( !myTraceQueue->isEmpty() ) {
      csTrace* trace = myTraceQueue->peek();
      if( traceIsPartOfCurrentEnsemble( trace ) ) {
        myTraceGather->addTrace( myTraceQueue->pop() );
      }
      else {  // Trace is NOT part of current ensemble
        myIsEnsembleFull = true;
        myIsFinishedProcessing = false;
        break;
      }
    }  // end while
  }
}
//------------------------------------------------------
// Move traces. Start with first trace (trace index 0)
void csModule::moveTracesFrom( csModule* module, int inPort ) {
  if( myMPIProcID != 0 ) return; // Only main process is allowed to pass on traces
  if( module->myNumTracesToBePassed == 0 ) {  // Does this ever occur?
    throw( cseis_geolib::csException("csModule::moveTracesFrom: No traces passed to be moved from module %s to this module. Is that good or bad? (MPI proc %d)", module->getName(), myMPIProcID) );
  }

  //  fprintf(stdout,"Move %d traces from %s to %s\n", module->myNumTracesToBePassed, module->getName(), getName() );
  // Case A) Single trace or fixed trace module
  if( myExecPhaseDef->execType() == EXEC_TYPE_INPUT || myExecPhaseDef->traceMode == TRCMODE_FIXED ) {
    int firstTraceIndex = 0;
    int numFixedTraces = myExecPhaseDef->numTraces;
    // A1. Move as many traces as possible from the other 'module' to the 'trace gather'.
    while( myTraceGather->numTraces() < numFixedTraces && firstTraceIndex < module->myNumTracesToBePassed ) {
      //      fprintf(stdout,"  %s toTraceGather chan %d\n", getName(), module->myTraceGather->trace( firstTraceIndex )->getTraceHeader()->intValue(5) );
      addNewTraceToGather( module->myTraceGather->trace( firstTraceIndex++ ), inPort );
      //      printf("A1. PreMod gather: %d, ThisMod gather: %d   --- trace1: %d  ntrc2BeTraced:(%d)\n", module->myTraceGather->numTraces(), myTraceGather->numTraces(), firstTraceIndex, module->myNumTracesToBePassed );
    }
    // A2. Move remaining traces to the 'trace queue'
    for( int itrc = firstTraceIndex; itrc < module->myNumTracesToBePassed; itrc++ ) {
      addNewTraceToQueue( module->myTraceGather->trace( itrc ), inPort );
      //      fprintf(stdout,"  %s toTraceQueue  chan %d\n", getName(), module->myTraceGather->trace( itrc )->getTraceHeader()->intValue(5) );
      //      printf("A2. PreMod gather: %d, ThisMod queue: %d   ---  trace1: %d\n", module->myTraceGather->numTraces(), myTraceQueue->size(), firstTraceIndex);
    }
  }
  // Case B) Ensemble trace module or module with variable number of traces
  else { // else if( myExecPhaseDef->traceMode == TRCMODE_ENSEMBLE ) {
    // First attempt: - Move all traces to the 'trace queue'.
    // - Next, try to move as many traces as possible from trace queue to trace gather

    // B1. The main case for multi-trace ensemble modules:
    // Check whether each trace is part of the current 'trace gather'.
    // If yes, add trace to trace gather. Otherwise, move this and all remaining traces into trace queue.
    if( mySuperHeader->numEnsembleKeys() > 0 ) {
      int firstTraceIndex = 0;
      // 1. If next ensemble has not been found yet, keep moving traces to trace gather
      while( firstTraceIndex < module->myNumTracesToBePassed ) {
        csTrace* trace = module->myTraceGather->trace( firstTraceIndex );
        //printf("moveTracesFrom: Trace is FULL: %s\n", myIsEnsembleFull ? "yes" : "no" );
        if( myIsEnsembleFull || !traceIsPartOfCurrentEnsemble( trace ) ) {
          //printf("moveTracesFrom: Trace is FULL or NOT part of current ensemble.  %s\n",  myName.c_str() );
          addNewTraceToQueue( trace, inPort );
        }
        else { //if( traceIsPartOfCurrentEnsemble( trace ) ) {
          // printf("moveTracesFrom: Trace is part of current ensemble.  %s\n", myName.c_str());
          addNewTraceToGather( trace, inPort );
        }
        firstTraceIndex += 1;
      }  // end while
    }
    // B2. No ensemble key set --> buffer entire data set directly into the trace gather
    else { //if( mySuperHeader->numEnsembleKeys() == 0 ) {
      for( int itrc = 0; itrc < module->myNumTracesToBePassed; itrc++ ) {
        addNewTraceToGather( module->myTraceGather->trace( itrc ), inPort );
      }
    }
  }
  // Remove traces from trace gather, but DO NOT FREE TRACES. Traces are freed when last module is reached.
  // Traces can not be freed yet because they are still in use by the remaining modules.
  module->myTraceGather->deleteTraces( 0, module->myNumTracesToBePassed );
  module->myNumTracesToBePassed = 0;
}

//------------------------------------------------------
//
void csModule::lastModuleTraceCleanup() {
  // Free processed traces from trace gather
  myTraceGather->freeTraces( 0, myNumTracesToBePassed );
  myNumTracesToBePassed = 0;
}
//-------------------------------------------------------------------
//
//

bool csModule::finishedProcessing() const {
  //  printf("finishedProcessing '%s': gather ntraces %d, numTracesPassed: %d, queue ntraces: %d, finished: %d\n", getName(), myTraceGather->numTraces(), myNumTracesToBePassed,  myTraceQueue->size(), (myTraceGather->numTraces() == myNumTracesToBePassed && myTraceQueue->isEmpty()) );
  return( myIsFinishedProcessing );
}
void csModule::setExecType( int execType ) {
  myExecPhaseDef->setExecType(execType);
}
int csModule::getExecType() const {
  return myExecPhaseDef->execType();
}
csTraceHeaderDef const* csModule::getHeaderDef() const {
  return myHeaderDef;
}
csExecPhaseDef const* csModule::getExecPhaseDef() const {
  return myExecPhaseDef;
}

//*********************************************************************
//
// Input port methods
//
//*********************************************************************
//
void csModule::setInputPorts( cseis_geolib::csVector<csModule const*> const* moduleList ) {
  myNumInputPorts = moduleList->size();

  if( myHeaderDef != NULL ) {
    delete myHeaderDef;
  }
  csTraceHeaderDef const** hdefPrev = new csTraceHeaderDef const*[myNumInputPorts];
  for( int inPort = 0; inPort < myNumInputPorts; inPort++ ) {
    hdefPrev[inPort] = moduleList->at(inPort)->myHeaderDef;
    if( inPort == 0 ) {
      mySuperHeader->set( moduleList->at(inPort)->mySuperHeader );
    }
    else {
      csSuperHeader* shdr = moduleList->at(inPort)->mySuperHeader;
      if( mySuperHeader->sampleInt != shdr->sampleInt ) {
        throw( cseis_geolib::csException("Module #%d (%s): Data from different input ports to this module have unequal sample interval", myUniqueID+1, getName()) );
      }
      else if( mySuperHeader->numSamples != shdr->numSamples ) {
        throw( cseis_geolib::csException("Module #%d (%s): Data from different input ports to this module have unequal number of samples", myUniqueID+1, getName()) );
      }
    }
  }
  myHeaderDef = new csTraceHeaderDef( myNumInputPorts, hdefPrev, myMemoryPoolManager );
  delete [] hdefPrev;

  if( myExecEnvPtr != NULL ) {
    delete myExecEnvPtr;
  }
  myExecEnvPtr = new csExecPhaseEnv( myHeaderDef, myExecPhaseDef, mySuperHeader );
}
//------------------------------------------------------------
//
void csModule::setDebugFlag( bool doDebug ) {
  myExecPhaseDef->myIsDebug = doDebug;
}
//------------------------------------------------------------
//
void csModule::mpiSetNumProc( int mpiNumProc ) {
  myMPINumProc = mpiNumProc;
}
int csModule::tempNumTraces() { return myTraceGather->numTraces(); }  // debug only, remove!

//--------------------------------------------------------------
std::string csModule::versionString() const {
  char ver[40];
  sprintf( ver, "%d.%d%c", myVersion[MAJOR], myVersion[MINOR], '\0' );
  return( std::string(ver) );
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csRSFHeader.h"
#include "csRSFReader.h"
#include <cstring>
#include <cmath>

#include "csStandardHeaders.h"
#include "csFileUtils.h"
#include "csGeolibUtils.h"

#include "csSort.h"
#include "csIOSelection.h"
#include "csSortManager.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace cseis_io;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: INPUT_RSF
 *
 */

namespace mod_input_rsf {
  struct VariableStruct {
    long traceCounter;
    long tracePerFileCounter;
    int numTracesIn;
    csRSFReader** rsfReaders;
    csRSFHeader** rsfHdr;
    int nTracesToRead;


    int dim1Counter;
    int dim2Counter;

    int hdrId_dim2;
    int hdrId_dim3;
    int hdrId_dim4;
    int hdrId_dim5;
    int hdrId_trcno;
    int hdrId_fileno;
    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_sou_z;

    int hdrId_time_samp1;
    int hdrId_time_samp1_us;
    int hdrId_delay_time;
    double delayTime;
    bool atEOF;

    int count_dim2;
    int count_dim3;
    int count_dim4;
    int count_dim5;
    int numDim;

    int numTracesBuffer;
    bool rev_byte_order;
    int mergeOption;

    int numFiles;
    std::string* filenameList;
    int currentFileIndex;

    bool isHdrSelection;
  };
  static int const MERGE_ALL    = 1;
  static int const MERGE_TRACE  = 2;
}

using mod_input_rsf::VariableStruct;

void initialize(  VariableStruct* vars, csLogWriter* writer );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_rsf_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  //----------------------------------------------
  // Initialise
  
  vars->mergeOption   = mod_input_rsf::MERGE_ALL;
  vars->traceCounter = 0;
  vars->rsfReaders   = NULL;
  vars->rsfHdr       = NULL;
  vars->nTracesToRead = 0;
  vars->numDim = 0;

  vars->hdrId_time_samp1    = -1;
  vars->hdrId_time_samp1_us = -1;
  vars->hdrId_dim2 = -1;
  vars->hdrId_dim3 = -1;
  vars->hdrId_dim4 = -1;
  vars->hdrId_dim5 = -1;
  vars->hdrId_trcno = -1;
  vars->hdrId_delay_time = -1;
  vars->hdrId_fileno = -1;
  vars->hdrId_sou_x = -1;
  vars->hdrId_sou_y = -1;
  vars->hdrId_sou_z = -1;
  vars->delayTime = 0;
  vars->atEOF = false;
  vars->count_dim2 = 0;
  vars->count_dim3 = 0;
  vars->count_dim4 = 0;
  vars->count_dim5 = 0;
  vars->numTracesBuffer = 0;
  vars->rev_byte_order = false;
  vars->isHdrSelection = false;

  //------------------------------------------------
  
  std::string headerName;
  std::string yesno;
  std::string text;

  //------------------------------------
  if( param->exists("merge") ) {
    param->getString( "merge", &text );
    if( !text.compare("all") ) {
      vars->mergeOption = mod_input_rsf::MERGE_ALL;
    }
    else if( !text.compare("trace") ) {
      vars->mergeOption = mod_input_rsf::MERGE_TRACE;
    }
    else {
      writer->error("Option for user parameter 'merge' not recognised: %s", text.c_str());
    }
  }

  vars->numFiles = param->getNumLines("filename");
  if( vars->numFiles == 1 ) vars->mergeOption = mod_input_rsf::MERGE_ALL;
  if( vars->numFiles == 0 ) {
    if( !param->exists("directory") ) writer->error("Required input parameter missing: 'filename'");
    std::string directory;
    param->getString("directory", &directory);
    std::string extension = "rsf";
    csVector<string> fileList;
    bool searchSubDirs = true;
    csFileUtils::retrieveFiles( directory, extension, &fileList, searchSubDirs, writer->getFile() );
    vars->numFiles = fileList.size();
    if( vars->numFiles == 0 ) writer->error("No RSF file found in directory %s, including its sub-directories", directory.c_str() );

    vars->filenameList = new std::string[vars->numFiles];
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      vars->filenameList[ifile] = fileList.at(ifile);
    }
    // Sort file names
    csSort<std::string>().simpleSort( vars->filenameList, vars->numFiles );
  }
  else {
    vars->filenameList = new std::string[vars->numFiles];
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      param->getStringAtLine("filename",&vars->filenameList[ifile],ifile);
    }
  }
  vars->currentFileIndex = 0;
  
  //----------------------------------------------------

  if( param->exists("ntraces") ) {
    param->getInt( "ntraces", &vars->nTracesToRead );
    if( vars->nTracesToRead < 0 ) {
      vars->nTracesToRead = 0;
    }
  }
  if( vars->nTracesToRead <= 0 ) vars->nTracesToRead = -1;  // Do not bother how many traces, read in all

  if( param->exists( "ntraces_buffer" ) ) {
    param->getInt( "ntraces_buffer", &vars->numTracesBuffer );
    if( vars->numTracesBuffer < 0 || vars->numTracesBuffer > 9999999 ) {
      writer->warning("Number of buffered traces out of range (=%d). Changed to default.", vars->numTracesBuffer);
      vars->numTracesBuffer = 0;
    }
  }

  if( param->exists("reverse_byte_order") ) {
    param->getString( "reverse_byte_order", &text );
    if( !text.compare("yes") ) {
      vars->rev_byte_order = true;
    }
    else if( !text.compare("no") ) {
      vars->rev_byte_order = false;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  int numSamplesOut = 0;
  if( param->exists("nsamples") ) {
    param->getInt( "nsamples", &numSamplesOut );
  }

  std::string selectionText = "";
  std::string selectionHdrName = "";
  if( param->exists( "header" ) ) {
    vars->isHdrSelection = true;
    vars->numTracesBuffer = 1;  // Force single trace buffer.
    csVector<std::string> valueList;
    std::string text;
    param->getAll( "header", &valueList );
    if( valueList.size() > 1 ) {
      writer->error("Currently, trace selection is only supported for one trace header. Number of supplied header names: %d", valueList.size());
    }
    selectionHdrName = valueList.at(0);
    valueList.clear();
    param->getString( "select", &selectionText );
  }

  //-------------------------------------------------------
  //
  vars->rsfReaders = new csRSFReader*[vars->numFiles];
  vars->rsfHdr     = new csRSFHeader*[vars->numFiles];
  int numSamplesPrev  = 0;
  float sampleIntPrev = 0;
  int numTracesPrev = 0;
  for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
    std::string filename = vars->filenameList[ifile];
    try {
      vars->rsfReaders[ifile] = new csRSFReader( filename, vars->numTracesBuffer, vars->rev_byte_order );
      vars->rsfHdr[ifile] = new csRSFHeader();
      vars->rsfReaders[ifile]->initialize( vars->rsfHdr[ifile] );
    }
    catch( csException& e ) {
      vars->rsfReaders = NULL;
      writer->error("Error when opening RSF file '%s'.\nSystem message: %s", filename.c_str(), e.getMessage() );
    }
    if( ifile > 0 ) {
      if( vars->rsfReaders[ifile]->numSamples() != numSamplesPrev ) {
        writer->error("Input file '%s' has different number of samples (=%d) than previous input file (=%d)", filename.c_str(), vars->rsfReaders[ifile]->numSamples(), numSamplesPrev);
      }
      if( vars->rsfReaders[ifile]->sampleInt() != sampleIntPrev ) {
        writer->error("Input file '%s' has different number of sample interval (=%f) than first input file (=%f)", filename.c_str(), vars->rsfReaders[ifile]->sampleInt(), sampleIntPrev);
      }
      if( vars->mergeOption == mod_input_rsf::MERGE_TRACE && vars->rsfReaders[ifile]->numTraces() != numTracesPrev ) {
        writer->error("Input file '%s' has different number of traces (=%d) than first input file (=%d)", filename.c_str(), vars->rsfReaders[ifile]->numTraces(), numTracesPrev);
      }
    }
    numSamplesPrev = vars->rsfReaders[ifile]->numSamples();
    sampleIntPrev  = vars->rsfReaders[ifile]->sampleInt();
    numTracesPrev  = vars->rsfReaders[ifile]->numTraces();
    if( numTracesPrev == 0 ) writer->error("Input file '%s' does not contain any traces", filename.c_str() );
  }

  vars->count_dim2   = 1;
  vars->count_dim3   = 1;
  vars->count_dim4   = 1;
  vars->count_dim5   = 1;

  // OLD  initialize( vars, log );

  if( numSamplesOut == 0 ) {
    shdr->numSamples = vars->rsfReaders[0]->numSamples();
  }
  else {
    shdr->numSamples = numSamplesOut;
  }
  shdr->sampleInt = vars->rsfReaders[0]->sampleInt();
  vars->delayTime = vars->rsfHdr[0]->o1;

  writer->line("");
  for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
    writer->line( "  File name(%2d):        %s    (# traces: %d)", ifile+1, vars->filenameList[ifile].c_str(), vars->rsfReaders[ifile]->numTraces() );
  }
  writer->line( "  Sample interval [ms]: %f", shdr->sampleInt );
  writer->line( "  Number of samples:    %d", shdr->numSamples );
  writer->line("");
  writer->line(" RSF ASCII file dump:");
  vars->rsfReaders[0]->dump( writer->getFile() );

  vars->numDim = vars->rsfReaders[0]->numDim();
  vars->hdrId_dim2 = hdef->addHeader( cseis_geolib::TYPE_DOUBLE, "dim2", "RSF data dimension 2" );
  if( vars->numDim >= 3 ) vars->hdrId_dim3 = hdef->addHeader( cseis_geolib::TYPE_DOUBLE, "dim3", "RSF data dimension 3" );
  if( vars->numDim >= 4 ) vars->hdrId_dim4 = hdef->addHeader( cseis_geolib::TYPE_DOUBLE, "dim4", "RSF data dimension 4" );
  if( vars->numDim >= 5 ) vars->hdrId_dim5 = hdef->addHeader( cseis_geolib::TYPE_DOUBLE, "dim5", "RSF data dimension 4" );
  vars->hdrId_trcno = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrId_delay_time = hdef->addStandardHeader( HDR_DELAY_TIME.name );
  hdef->addStandardHeader( HDR_FILENO.name );
  vars->hdrId_fileno = hdef->headerIndex( HDR_FILENO.name );

  if( vars->rsfHdr[0]->sou_x != 0 || vars->rsfHdr[0]->sou_y != 0 || vars->rsfHdr[0]->sou_z != 0 ) {
    vars->hdrId_sou_x = hdef->addStandardHeader( HDR_SOU_X.name );
    vars->hdrId_sou_y = hdef->addStandardHeader( HDR_SOU_Y.name );
    vars->hdrId_sou_z = hdef->addStandardHeader( HDR_SOU_Z.name );
  }

  //--------------------------------------------------------------------------------
  // ...important to do the following after all headers have been set for hdef
  //
  if( vars->isHdrSelection ) {
    hdef->resetByteLocation();  // This will otherwise be done by base system AFTER init phase
    if( !hdef->headerExists( selectionHdrName ) ) {
      writer->error("Selection trace header '%s' is not defined in input file '%s'", selectionHdrName.c_str(), vars->filenameList[0].c_str());
    }
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      bool success = vars->rsfReaders[ifile]->setSelection( selectionText, selectionHdrName, cseis_geolib::csIOSelection::SORT_NONE, cseis_geolib::csSortManager::SIMPLE_SORT );
      if( !success ) {
        writer->warning("Error occurred when intializing header selection for RSF file '%s'.\n --> No input traces found that match specified selection '%s' for header '%s'.\n",
                     vars->filenameList[ifile].c_str(), selectionText.c_str(), selectionHdrName.c_str() );
      }
    }
  }

  vars->traceCounter = 0;
  vars->tracePerFileCounter = 0;
  vars->numTracesIn = vars->rsfReaders[0]->numTraces();

  shdr->domain = vars->rsfReaders[0]->getCSEISDomain();
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_input_rsf_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csTrace* trace = traceGather->trace(0);


  if( vars->atEOF ) {
    traceGather->freeAllTraces();
    return;
  }
  if( vars->nTracesToRead > 0 && vars->traceCounter == vars->nTracesToRead ) {
    vars->atEOF = true;
    traceGather->freeAllTraces();
    return;
  }

  csSuperHeader const*    shdr = env->superHeader;
  csTraceHeader* trcHdr = trace->getTraceHeader();
  float* samples = trace->getTraceSamples();

  csRSFReader* readerPtr = vars->rsfReaders[vars->currentFileIndex];
  bool success = readerPtr->getNextTrace( (byte_t*)samples, shdr->numSamples );

  //  if( vars->tracePerFileCounter == vars->numTracesIn ) {  // All traces of current file have been read in
  if( !success ) {  // All traces of current file have been read in
    if( vars->mergeOption == mod_input_rsf::MERGE_ALL ) {
      // Close current file, set pointer to next input file
      delete vars->rsfReaders[vars->currentFileIndex];
      vars->rsfReaders[vars->currentFileIndex] = NULL;
      vars->currentFileIndex += 1;

      while( !success && vars->currentFileIndex < vars->numFiles ) {
        success = vars->rsfReaders[vars->currentFileIndex]->getNextTrace( (byte_t*)samples, shdr->numSamples );
        if( !success ) {
          delete vars->rsfReaders[vars->currentFileIndex];
          vars->rsfReaders[vars->currentFileIndex] = NULL;
          vars->currentFileIndex += 1;
        }
      }
      if( vars->currentFileIndex == vars->numFiles ) {
        vars->atEOF = true;
        traceGather->freeAllTraces();
        return;
      }
      readerPtr = vars->rsfReaders[vars->currentFileIndex];
      vars->numTracesIn = readerPtr->numTraces();
      vars->tracePerFileCounter = 0;
    }
    else { // MERGE_TRACE
      vars->atEOF = true;
      traceGather->freeAllTraces();
      return;
    }
  }

  int traceIndexInFile = vars->rsfReaders[vars->currentFileIndex]->getCurrentTraceIndex();

  // (int)vars->tracePerFileCounter
  trcHdr->setDoubleValue( vars->hdrId_dim2, readerPtr->computeDim2( traceIndexInFile) );
  if( vars->numDim >= 3 ) {
    trcHdr->setDoubleValue( vars->hdrId_dim3, readerPtr->computeDim3( traceIndexInFile ) );
    if( vars->numDim >= 4 ) {
      trcHdr->setDoubleValue( vars->hdrId_dim4, readerPtr->computeDim4( traceIndexInFile ) );
      if( vars->numDim >= 5 ) {
        trcHdr->setDoubleValue( vars->hdrId_dim5, readerPtr->computeDim5( traceIndexInFile ) );
      }
    }
  }
  trcHdr->setDoubleValue( vars->hdrId_delay_time, vars->delayTime );
  trcHdr->setIntValue( vars->hdrId_fileno, vars->currentFileIndex+1 );
  trcHdr->setIntValue( vars->hdrId_trcno, (int)vars->traceCounter+1 );
  if( vars->hdrId_sou_x >= 0 ) {
    trcHdr->setDoubleValue( vars->hdrId_sou_x, vars->rsfHdr[vars->currentFileIndex]->sou_x );
    trcHdr->setDoubleValue( vars->hdrId_sou_y, vars->rsfHdr[vars->currentFileIndex]->sou_y );
    trcHdr->setDoubleValue( vars->hdrId_sou_z, vars->rsfHdr[vars->currentFileIndex]->sou_z );
  }

  vars->traceCounter += 1;  // Count all traces read in so far
  if( vars->mergeOption == mod_input_rsf::MERGE_ALL ) {
    vars->tracePerFileCounter += 1; // Count traces read in from each file
  }
  else { // MERGE_TRACE
    vars->currentFileIndex += 1; // Switch to next input file
    if( vars->currentFileIndex == vars->numFiles ) {
      vars->currentFileIndex = 0; // Switch back to first input file
      vars->tracePerFileCounter += 1;
    }
  }

  return;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_rsf_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_RSF", "Input data in RSF format" );

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "directory", "Name of directory to search for input files", NUM_VALUES_FIXED, "All RSF files in this directory and its subdirectories will be read" );
  pdef->addValue( "", VALTYPE_STRING );

  pdef->addParam( "nsamples", "Number of samples to read in", NUM_VALUES_FIXED,
                  "If number of samples in input data set is smaller, traces will be filled with zeros. Set 0 to set number of samples from input data set.");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of samples to read in" );

  pdef->addParam( "ntraces", "Number of traces to read in", NUM_VALUES_FIXED,
                  "Input of traces will stop when all traces have been read in, or if the number of traces specified has been reached. Traces will not be filled up to the specified range");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to read in" );

  pdef->addParam( "merge", "Method to merge traces from multiple input files", NUM_VALUES_FIXED);
  pdef->addValue( "all", VALTYPE_OPTION );
  pdef->addOption( "all", "Read in all traces of file 1, then all traces of file 2 etc" );
  pdef->addOption( "trace", "Read in one trace per input file, then repeat until first file has been fully read in. Stop then." );

  pdef->addParam( "ntraces_buffer", "Number of traces to read into buffer at once", NUM_VALUES_FIXED,
                  "Reading in a large number of traces at once may enhance performance, but requires more memory" );
  pdef->addValue( "20", VALTYPE_NUMBER, "Number of traces to buffer" );

  pdef->addParam( "reverse_byte_order", "Reverse byte order of input file (endian byte swapping)", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Reverse byte order of input file" );
  pdef->addOption( "no", "Do not reverse byte order" );


  pdef->addParam( "header", "Name of trace header used for trace selection", NUM_VALUES_FIXED, "Use in combination with user parameter 'select'. NOTE: With the current Seaseis disk data format, selecting traces on input is typically slower than reading in all traces and making the trace selection later on, e.g. by using module 'SELECT'" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of header values", NUM_VALUES_FIXED, "Only traces which fit the trace value selection will be read in. Use in combination with user parameter 'header'" );
  pdef->addValue( "", VALTYPE_STRING, "Selection string. See documentation for more detailed description of selection syntax" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_input_rsf_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_input_rsf::VariableStruct* vars = reinterpret_cast<mod_input_rsf::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_input_rsf_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_input_rsf::VariableStruct* vars = reinterpret_cast<mod_input_rsf::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->rsfReaders != NULL ) {
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      if( vars->rsfReaders[ifile] != NULL ) {
        delete vars->rsfReaders[ifile];
        vars->rsfReaders[ifile] = NULL;
      }
    }
    delete [] vars->rsfReaders;
    vars->rsfReaders = NULL;
  }
  if( vars->rsfHdr != NULL ) {
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      if( vars->rsfHdr[ifile] != NULL ) {
        delete vars->rsfHdr[ifile];
        vars->rsfHdr[ifile] = NULL;
      }
    }
    delete [] vars->rsfHdr;
    vars->rsfHdr = NULL;
  }
  if( vars->filenameList != NULL ) {
    delete [] vars->filenameList;
    vars->filenameList = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_input_rsf_( csParamDef* pdef ) {
  params_mod_input_rsf_( pdef );
}
extern "C" void _init_mod_input_rsf_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_rsf_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_rsf_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_rsf_( env, writer );
}
extern "C" void _exec_mod_input_rsf_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_rsf_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_rsf_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_rsf_( env, writer );
}

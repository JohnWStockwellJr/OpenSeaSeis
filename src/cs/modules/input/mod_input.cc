/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSeismicReader.h"
#include "csStandardHeaders.h"
#include "csFlexHeader.h"
#include "csIOSelection.h"
#include "csSortManager.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: INPUT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_input {
  struct VariableStruct {
    cseis_system::csSeismicReader** readers;
    char* hdrValueBlock;
    char* hdrValueBlock_dummy;
    int hdrValueBlockSize;
    int numFiles;
    int nTracesToRead;
    int currentFile;
    std::string* filenames;
    int numTracesBuffer;
    bool atEOF;
    int mergeOption;
    long traceCounter;
    int hdrId_fileno;
    int forcedNumSamples;
    float sampleInt_override;

    // Merging of traces from different input files with same header value
    type_t hdrType_merge;     // ...
    int    mergeOrder;        // Read in traces by increasing or decreasing order of merge header value
    csFlexHeader*  mergeHdrValues;    // Merge header value of current trace, for each input file
    csFlexHeader currentMergeHdrValue;// Current value of merge header of the trace that is read in next
    int currentFilePointerIndex;      // Index of current file where next trace shall be read in
    csVector<int>* fileIndexList;     // List of index pointers to remaining files (that have not been read in fully yet)

    cseis_system::csTraceHeaderDef** trcHdrDef;
    cseis_system::csTraceHeader**    trcHdr;

    bool* hdrIsEqual; // flag that indicates if the trace header block of file N equals the trace header block from file 1
    int** hdrId;
    type_t** hdrType;

    bool isHdrSelection;
    int sortOrder;
    int sortMethod;
  };
  static int const MERGE_ALL    = 1;
  static int const MERGE_TRACE  = 2;
  static int const MERGE_HEADER = 3;
  static int const MERGE_HDR_COPY = 4;

  static int const MERGE_INCREASING = 11;
  static int const MERGE_DECREASING = 22;

  static int const NOT_AVAILABLE = -1;
}
using mod_input::VariableStruct;

//*********************************************************************************bool****************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->mergeOption   = mod_input::MERGE_ALL;
  vars->readers       = NULL;
  vars->hdrValueBlock = NULL;
  vars->hdrValueBlock_dummy = NULL;
  vars->numFiles      = 0;
  vars->filenames     = NULL;
  vars->atEOF         = false;
  vars->nTracesToRead = 0;
  vars->trcHdrDef  = NULL;
  vars->trcHdr     = NULL;
  vars->hdrIsEqual = NULL;
  vars->hdrId      = NULL;
  vars->hdrType    = NULL;
  vars->fileIndexList = NULL;
  vars->forcedNumSamples = -1;

  vars->hdrType_merge  = TYPE_UNKNOWN;
  vars->mergeHdrValues = NULL;
  vars->mergeOrder     = mod_input::MERGE_INCREASING;
  vars->currentFilePointerIndex = 0;

  vars->hdrValueBlockSize = 0;
  vars->currentFile       = 0;
  vars->numTracesBuffer   = 0;
  vars->traceCounter      = 0;
  vars->hdrId_fileno      = -1;
  vars->sampleInt_override = -1;

  vars->isHdrSelection = false;
  vars->sortOrder = cseis_geolib::csIOSelection::SORT_NONE;
  vars->sortMethod = cseis_geolib::csSortManager::SIMPLE_SORT;

//------------------------------------------------------------
  vars->numFiles = param->getNumLines( "filename" );
  if( vars->numFiles <= 0 ) {
    writer->error("No input file specified. No parameter 'filename' found.");
  }
  vars->filenames = new std::string[vars->numFiles];
  for( int i = 0; i < vars->numFiles; i++ ) {
    string filename;
    param->getStringAtLine( "filename", &filename, i );
    int length = (int)filename.length();
    if( length < 6 || (filename.substr( length-6, 6 ).compare(".oseis") && filename.substr( length-6, 6 ).compare(".cseis")) ) {
      writer->warning("Input file name does not have standard Cseis extension '.cseis' (or '.oseis'): '%s'\n", filename.c_str());
    }
    vars->filenames[i] = filename;
  }

  if( param->exists( "override_sampleint" ) ) {
    param->getFloat( "override_sampleint", &vars->sampleInt_override );
  }

  if( param->exists( "ntraces_buffer" ) ) {
    param->getInt( "ntraces_buffer", &vars->numTracesBuffer );
    if( vars->numTracesBuffer <= 0 || vars->numTracesBuffer > 999999 ) {
      writer->warning("Number of buffered traces out of range (=%d). Changed to 0.", vars->numTracesBuffer);
      vars->numTracesBuffer = 0;
    }
  }

  if( param->exists( "ntraces" ) ) {
    param->getInt( "ntraces", &vars->nTracesToRead );
    if( vars->nTracesToRead < 0 ) {
      vars->nTracesToRead = 0;
    }
  }

  if( param->exists( "nsamples" ) ) {
    param->getInt( "nsamples", &vars->forcedNumSamples );
    if( vars->forcedNumSamples <= 0 ) {
      vars->forcedNumSamples = -1;
    }
    else {
      shdr->numSamples = vars->forcedNumSamples;
    }
  }

  //----------------------------------------------------
  string mergeHeaderName = ""; 
  bool enableRandomAccess = false;

  if( param->exists("merge") ) {
    string text;
    param->getString( "merge", &text );
    if( !text.compare("all") ) {
      vars->mergeOption = mod_input::MERGE_ALL;
    }
    else if( !text.compare("trace") ) {
      vars->mergeOption = mod_input::MERGE_TRACE;
    }
    else if( !text.compare("hdr_copy_file1") ) {
      vars->mergeOption = mod_input::MERGE_HDR_COPY;
    }
    else if( !text.compare("header") ) {
      vars->mergeOption = mod_input::MERGE_HEADER;
      param->getString( "header_merge", &mergeHeaderName );
      enableRandomAccess = true;
    }
    else {
      writer->error("Option for user parameter 'merge' not recognised: %s", text.c_str());
    }
  }

  //-------------------------------------------------------------------------

  /*
   * Trace selection based on header value
   * Note that with the current version of the input data format there is no performance increase, rather a decrease compared
   * to reading in all traces and selecting later on using the module $SELECT.
   * Best to use this functionality only for sorting of input data.
   */
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
    if( param->exists( "sort" ) ) {
      string text;
      param->getString( "sort", &text );
      if( !text.compare("no") ) {
        vars->sortOrder = cseis_geolib::csIOSelection::SORT_NONE;
      }
      else if( !text.compare("increasing") ) {
        vars->sortOrder = cseis_geolib::csIOSelection::SORT_INCREASING;
      }
      else if( !text.compare("decreasing") ) {
        vars->sortOrder = cseis_geolib::csIOSelection::SORT_DECREASING;
      }
      else {
        writer->error("Unknow option: %s", text.c_str());
      }
      if( param->getNumValues( "sort" ) > 1 ) {
        param->getString( "sort", &text, 1 );
        if( !text.compare("simple") ) {
          vars->sortMethod = cseis_geolib::csSortManager::SIMPLE_SORT;
        }
        else if( !text.compare("tree") ) {
          vars->sortMethod = cseis_geolib::csSortManager::TREE_SORT;
        }
        else {
          writer->error("Unknown option: %s", text.c_str());
        }
      }
    }
    if( vars->sortOrder == cseis_geolib::csIOSelection::SORT_NONE ) {
      writer->warning("With the current Seaseis disk data format, selecting traces on input using user parameters 'header' & 'select' is typically slower than reading in all traces and performing the selection afterwards, e.g. by using module 'SELECT'. It is recommended to use input selection only when traces shall be sorted on input");
    }
  }

  //-------------------------------------------------------------------------
  // Open all input files
  //
  // (1) Allocate fields that are required in any case
  //
  vars->readers       = new csSeismicReader*[vars->numFiles];
  vars->hdrIsEqual    = new bool[vars->numFiles];
  vars->readers[0]    = NULL;
  vars->hdrIsEqual[0] = true;

  // (2) Allocate and init fields that are only required when more than one input file is being read in
  if( vars->numFiles > 1 ) {
    vars->trcHdrDef  = new csTraceHeaderDef*[vars->numFiles];
    vars->trcHdr     = new csTraceHeader*[vars->numFiles];
    vars->hdrId      = new int*[vars->numFiles];
    vars->hdrType    = new type_t*[vars->numFiles];
    vars->fileIndexList = new csVector<int>(vars->numFiles);
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      vars->fileIndexList->insertEnd(ifile);
      vars->hdrIsEqual[ifile] = true;
      vars->readers[ifile]    = NULL;
      vars->trcHdrDef[ifile]  = NULL;
      vars->trcHdr[ifile]     = NULL;
      vars->hdrId[ifile]      = NULL;
      vars->hdrType[ifile]    = NULL;
    }
  }
  else {  // Only single input file --> input order doesn't matter. Use simplest one..
    vars->mergeOption = mod_input::MERGE_ALL;
  }

  if( vars->mergeOption == mod_input::MERGE_HEADER ) {
    vars->mergeHdrValues = new csFlexHeader[vars->numFiles];
  }

  if( param->exists("header") ) {
    enableRandomAccess = true;
  }

  try {
    // (3) Open first input file
    // First file determines trace headers. All following files have to either match or are matched to trace headers in first file
    vars->readers[0] = new csSeismicReader( vars->filenames[0], enableRandomAccess, vars->numTracesBuffer );
    bool success = vars->readers[0]->readFileHeader( shdr, hdef, &vars->hdrValueBlockSize, writer->getFile() );
    if( vars->forcedNumSamples > 0 ) shdr->numSamples = vars->forcedNumSamples;
    if( !success ) {
      writer->error("Unknown error occurred when reading file header from SeaSeis file '%s'.\n", vars->filenames[0].c_str() );
    }
    if( vars->mergeOption == mod_input::MERGE_HEADER ) {
      if( !hdef->headerExists( mergeHeaderName ) ) {
        writer->error("Required merge header '%s' is not set in input file '%s'", mergeHeaderName.c_str(), vars->filenames[0].c_str());
      }
      vars->hdrType_merge = hdef->headerType( mergeHeaderName );
    }
//    writer->line("");
    // NOTE: vars->trcHdrDef[0] is not set because this is the same as hdef defined in this module
    int maxHdrValueBlockSize = vars->hdrValueBlockSize;

    // (4) Open additional files
    for( int ifile = 1; ifile < vars->numFiles; ifile++ ) {
      vars->trcHdrDef[ifile] = new csTraceHeaderDef( hdef );
      int newBlockSize;

      vars->readers[ifile] = new csSeismicReader( vars->filenames[ifile], enableRandomAccess, vars->numTracesBuffer );
      csSuperHeader newShdr;
      success = vars->readers[ifile]->readFileHeader( &newShdr, vars->trcHdrDef[ifile], &newBlockSize );
      vars->trcHdrDef[ifile]->resetByteLocation(); // Don't forget to reset byte location!
      // ^..Resetting byte loc is usually done by SeaSeis system after init phase but has to be done manually for the additional trcHdrDef objects

      if( newBlockSize > maxHdrValueBlockSize ) maxHdrValueBlockSize = newBlockSize;
      if( !success ) {
        writer->error("Unknown error occurred when reading file header from SeaSeis file '%s'.\n", vars->filenames[ifile].c_str() );
      }
      else if( newShdr.numSamples != shdr->numSamples || newShdr.sampleInt != shdr->sampleInt ) {

        if( vars->sampleInt_override <= 0 && newShdr.sampleInt != shdr->sampleInt && vars->forcedNumSamples < 0 ) {
          writer->error("SeaSeis file '%s' has different properties (num samples = %d, sampleInt = %.5f) than previous file(s) (num samples = %d, sampleInt = %.5f.\nYou can force the number of samples by setting the user parameter 'nsamples', but different sample intervals are not supported.\nUse user parameter 'override_sampleint' to override samnple interval of input data sets.\n",
                        vars->filenames[ifile].c_str(), newShdr.numSamples, newShdr.sampleInt, shdr->numSamples, shdr->sampleInt );
        }
      }

      if( newBlockSize != vars->hdrValueBlockSize || !vars->trcHdrDef[ifile]->equals(hdef) ) {
        if( newBlockSize != vars->hdrValueBlockSize ) {
          writer->warning("SeaSeis file '%s' has a different header block size (%d != %d) than first file. This is supported but will slow down the merge process.\n",  vars->filenames[ifile].c_str(), newBlockSize, vars->hdrValueBlockSize );
        }
        else {
          writer->warning("SeaSeis file '%s' has a different header definition than first file. This is supported but will slow down the merge process.\n",  vars->filenames[ifile].c_str() );
        }
        vars->hdrIsEqual[ifile] = false;
        vars->trcHdr[ifile]     = new csTraceHeader();
        vars->trcHdr[ifile]->setHeaders( vars->trcHdrDef[ifile] );
        int numHeaders = hdef->numHeaders();
        vars->hdrId[ifile]      = new int[hdef->numHeaders()];
        vars->hdrType[ifile]    = new type_t[hdef->numHeaders()];
        for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
          vars->hdrId[ifile][ihdr]   = mod_input::NOT_AVAILABLE;
          vars->hdrType[ifile][ihdr] = hdef->headerType( ihdr );
          std::string name           = hdef->headerName( ihdr );
          if( vars->trcHdrDef[ifile]->headerExists( name ) ) {
            vars->hdrId[ifile][ihdr] = vars->trcHdrDef[ifile]->headerIndex( name );
          }
        }
        if( edef->isDebug() ) {
          for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
            writer->line("Input trace header #%2d (%d): %d  '%s'",ihdr, ifile, vars->hdrId[ifile][ihdr], hdef->headerName(ihdr).c_str());
          }
        }
      }  // end: else if...
      else if( vars->mergeOption == mod_input::MERGE_HEADER &&
               (!vars->trcHdrDef[ifile]->headerExists( mergeHeaderName ) || (vars->trcHdrDef[ifile]->headerType( mergeHeaderName ) != vars->hdrType_merge ) ) ) {
        writer->error("Required merge header '%s' is not set in input file '%s', or has different type", mergeHeaderName.c_str(), vars->filenames[ifile].c_str());
      }
    }
    // Determine maximum numSamples
    //    int maxNumSamples = shdr->numSamples;
    //  bool found = false;
    //  for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
    //    int numSamples = vars->readers[ifile]->numSamples();
    //    found = ( numSamples != maxNumSamples );
    //    if( numSamples > maxNumSamples ) maxNumSamples = numSamples;
    //  }
    vars->hdrValueBlock = new char[maxHdrValueBlockSize];
    if( vars->mergeOption == mod_input::MERGE_HDR_COPY ) vars->hdrValueBlock_dummy = new char[maxHdrValueBlockSize];
  }
  catch( csException& exc ) {
    writer->error("Error occurred when opening SeaSeis file. System message:\n%s", exc.getMessage() );
  }

  if( !hdef->headerExists( HDR_FILENO.name ) ) {
    hdef->addStandardHeader( HDR_FILENO.name );
  }
  vars->hdrId_fileno = hdef->headerIndex( HDR_FILENO.name );

  writer->line("");
  for( int i = 0; i < vars->numFiles; i++ ) {
    writer->line("  File name #%-2d:      %s", i+1, vars->filenames[i].c_str() );
  }
  writer->line("  Header block size:             %d", vars->hdrValueBlockSize );
  writer->line("  Max number of buffered traces: %d", vars->readers[0]->numTracesCapacity() );
  shdr->dump( writer->getFile() );
  writer->line("");
  writer->flush();

  //--------------------------------------------------------------------------------
  // ...important to do the following after all headers have been set for hdef and all other input files
  //
  if( vars->isHdrSelection ) {
    hdef->resetByteLocation();  // This will otherwise be done by base system AFTER init phase
    if( !hdef->headerExists( selectionHdrName ) ) {
      writer->error("Selection trace header '%s' is not defined in input file '%s'", selectionHdrName.c_str(), vars->filenames[0].c_str());
    }
    type_t selectionHdrType = hdef->headerType( selectionHdrName );
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      if( ifile > 0 ) {
        if( !vars->trcHdrDef[ifile]->headerExists( selectionHdrName ) || vars->trcHdrDef[ifile]->headerType( selectionHdrName ) != selectionHdrType ) {
          writer->error("Selection trace header '%s' is not defined in input file '%s', or has different type",
                     selectionHdrName.c_str(), vars->filenames[ifile].c_str());
        }
      }

      //      bool success = vars->readers[ifile]->setSelection( selectionText, selectionHdrName, selectionHdrType, vars->sortOrder );
      bool success = vars->readers[ifile]->setSelection( selectionText, selectionHdrName, vars->sortOrder, vars->sortMethod );
      if( !success ) {
        writer->error("Error occurred when intializing header selection for SeaSeis file '%s'.\n --> No input traces found that match specified selection '%s' for header '%s'.\n",
                   vars->filenames[ifile].c_str(), selectionText.c_str(), selectionHdrName.c_str() );
      }
    }
  }

  // For header merge, retrieve first header value
  if( vars->mergeOption == mod_input::MERGE_HEADER ) {
    hdef->resetByteLocation();  // This will otherwise be done by base system AFTER init phase
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      vars->readers[ifile]->setHeaderToPeek( mergeHeaderName );
      bool success = vars->readers[ifile]->peekHeaderValue( &vars->mergeHdrValues[ifile] );
      if( !success ) {
        writer->error("Peek into of first header value of header '%s' for file '%s' failed.\nUnknown error cause...",
                   mergeHeaderName.c_str(), vars->filenames[ifile].c_str() );
      }
      else if( edef->isDebug() ) {
        writer->line("Header value of first trace, file #%-2d:  %s", ifile+1, vars->mergeHdrValues[ifile].toString().c_str() );
      }
    }
    csFlexHeader minValue = vars->mergeHdrValues[0];
    vars->currentFilePointerIndex = 0;
    for( int ifile = 1; ifile < vars->numFiles; ifile++ ) {
      if( vars->mergeHdrValues[ifile] < minValue ) {
        minValue = vars->mergeHdrValues[ifile];
        vars->currentFilePointerIndex = ifile;
      }
    }
    vars->currentMergeHdrValue = minValue;
    vars->currentFile = vars->fileIndexList->at(vars->currentFilePointerIndex);
    if( edef->isDebug() ) {
      writer->line("---------------------------------\nHeader value of first trace, file #%-2d:  %s\n",
                vars->currentFile+1, vars->mergeHdrValues[vars->currentFile].toString().c_str() );
    }
  }


  if( edef->isDebug() ) {
    hdef->dump( stdout );
  }
  vars->traceCounter = 0;
  if( vars->sampleInt_override > 0 ) shdr->sampleInt = vars->sampleInt_override;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_input_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const*  shdr = env->superHeader;
  csTraceHeaderDef const*  hdef = env->headerDef;

  csTrace* trace = traceGather->trace(0);

  if( vars->atEOF ) {
    traceGather->freeAllTraces();
    return;
  }
  
  csTraceHeader* trcHdr = trace->getTraceHeader();
  float* samples = trace->getTraceSamples();

  if( vars->nTracesToRead > 0 && vars->nTracesToRead == vars->traceCounter ) {
    vars->atEOF = true;
    traceGather->freeAllTraces();
    return;
  }

  //------------------------------------------------------------------------------------
  int currentFileNumberToSet = vars->currentFile;

  if( (vars->mergeOption == mod_input::MERGE_TRACE || vars->mergeOption == mod_input::MERGE_HDR_COPY) && vars->currentFile == vars->numFiles ) vars->currentFile = 0;

  try {
    char* hdrValueBlockPtr = vars->hdrValueBlock;
    if( vars->mergeOption == mod_input::MERGE_HDR_COPY && vars->currentFile > 0 ) {
      hdrValueBlockPtr = vars->hdrValueBlock_dummy; // For all but first input file, read trace header into dummy block that is going to be discarded.
    }
    //    bool success = vars->readers[vars->currentFile]->readTrace( samples, vars->hdrValueBlock, shdr->numSamples );
    bool success = vars->readers[vars->currentFile]->readTrace( samples, hdrValueBlockPtr, shdr->numSamples );
    
    if( !success ) {
      if( vars->mergeOption == mod_input::MERGE_ALL ) {
        delete vars->readers[vars->currentFile];
        vars->readers[vars->currentFile] = NULL;
        vars->currentFile += 1;
        while( !success && vars->currentFile < vars->numFiles ) {
          success = vars->readers[vars->currentFile]->readTrace( samples, vars->hdrValueBlock, shdr->numSamples );
//          if( edef->isDebug() ) fprintf(stdout,"Read in next file(%d, success = %s): %s\n", vars->currentFile, success ? "true" : "false", vars->filenames[vars->currentFile].c_str());
          if( !success ) {
            delete vars->readers[vars->currentFile];
            vars->readers[vars->currentFile] = NULL;
            vars->currentFile += 1;
          }
        }
        if( vars->currentFile == vars->numFiles ) {
          vars->atEOF = true;
          traceGather->freeAllTraces();
          return;
        }
        if( !success ) {
          traceGather->freeAllTraces();
          return;
        }
        currentFileNumberToSet = vars->currentFile;
      }
      else if( vars->mergeOption == mod_input::MERGE_TRACE || vars->mergeOption == mod_input::MERGE_HDR_COPY ) {
        for( int i = 0; i < vars->numFiles; i++ ) {
          delete vars->readers[i];
          vars->readers[i] = NULL;
        }
        delete [] vars->readers;
        vars->readers = NULL;
        vars->atEOF = true;
        traceGather->freeAllTraces();
        return;
      }
      else if( vars->mergeOption == mod_input::MERGE_HEADER ) {
        if( edef->isDebug() ) {
          writer->line("Unsuccessful reading from file %d (pointer %d)", vars->currentFile+1, vars->currentFilePointerIndex);
        }
        delete vars->readers[vars->currentFile];
        vars->readers[vars->currentFile] = NULL;
        vars->fileIndexList->remove( vars->currentFilePointerIndex );  // Remove current file from file list
        if( vars->fileIndexList->size() == 0 ) {  // This was the last trace, from the last file. Exit now.
          vars->atEOF = true;
          traceGather->freeAllTraces();
          return;
        }
        csFlexHeader minValue = vars->mergeHdrValues[vars->fileIndexList->at(0)];
        vars->currentFilePointerIndex = 0;
        for( int index = 1; index < vars->fileIndexList->size(); index++ ) {
          int fileIndex = vars->fileIndexList->at(index);
          if( vars->mergeHdrValues[fileIndex] < minValue ) {
            minValue = vars->mergeHdrValues[fileIndex];
            vars->currentFilePointerIndex = index;
          }
        }
        vars->currentMergeHdrValue = minValue;
        vars->currentFile = vars->fileIndexList->at(vars->currentFilePointerIndex);
        success = vars->readers[vars->currentFile]->readTrace( samples, vars->hdrValueBlock, shdr->numSamples );
        if( !success ) writer->error("Unexpected end of file encountered for file '%s'... File corrupted..?", vars->filenames[vars->currentFile].c_str());
        if( edef->isDebug() ) {
          writer->line("Header value of first trace, file #%-2d:  %s\n", vars->currentFile+1, vars->mergeHdrValues[vars->currentFile].toString().c_str() );
        }
        currentFileNumberToSet = vars->currentFile;
      }
    } // END: Read operation was unsuccessful
    // Read operation was a 'success':
    else if( vars->mergeOption == mod_input::MERGE_HEADER ) {
      success = vars->readers[vars->currentFile]->peekHeaderValue( &vars->mergeHdrValues[vars->currentFile] );
      if( !success ) {
        if( edef->isDebug() ) {
          writer->line("EOF, file #%-2d, size of index list: %d\n", vars->currentFile+1, vars->fileIndexList->size() );
        }
        // Nothing. Let next read statement fail, and deal with result then
      }
      else if( vars->mergeHdrValues[vars->currentFile] != vars->currentMergeHdrValue ) {
        csFlexHeader minValue = vars->mergeHdrValues[vars->fileIndexList->at(0)];
        vars->currentFilePointerIndex = 0;
        for( int index = 1; index < vars->fileIndexList->size(); index++ ) {
          int fileIndex = vars->fileIndexList->at(index);
          if( vars->mergeHdrValues[fileIndex] < minValue ) {
            minValue = vars->mergeHdrValues[fileIndex];
            vars->currentFilePointerIndex = index;
          }
        }
        vars->currentMergeHdrValue = minValue;
        vars->currentFile = vars->fileIndexList->at(vars->currentFilePointerIndex);
        if( edef->isDebug() ) {
          writer->line("Header value of first trace, file #%-2d:  %s\n", vars->currentFile+1, vars->mergeHdrValues[vars->currentFile].toString().c_str() );
        }
      }
    } // END: read operation was a success, and option MERGE_HEADER
  }
  catch( csException& exc ) {
    writer->error("Error occurred when reading SeaSeis file. System message:\n%s", exc.getMessage() );
  }
  
  // Set header value block
  if( vars->hdrIsEqual[currentFileNumberToSet] || vars->mergeOption == mod_input::MERGE_HDR_COPY ) {
    //    writer->line("Num headers/bytes %d %d, current file: %d", trcHdr->numHeaders(), hdef->getTotalNumBytes(), currentFileNumberToSet );
    trcHdr->setTraceHeaderValueBlock( vars->hdrValueBlock, vars->hdrValueBlockSize );
  }
  // Header def of input file was different than that of first file --> Set each trace header value manually
  else {
    // writer->line("Num headers: %d, bytes %d %d, current file: %d", vars->trcHdr[currentFileNumberToSet]->numHeaders(), hdef->getTotalNumBytes(), vars->trcHdrDef[currentFileNumberToSet]->getTotalNumBytes(), currentFileNumberToSet );
    int numHeaders     = hdef->numHeaders();
    int* hdrIdPtr      = vars->hdrId[currentFileNumberToSet];
    type_t* hdrTypePtr = vars->hdrType[currentFileNumberToSet];
    csTraceHeader* trcHdrInPtr = vars->trcHdr[currentFileNumberToSet];
    trcHdrInPtr->setTraceHeaderValueBlock( vars->hdrValueBlock, vars->trcHdrDef[currentFileNumberToSet]->getTotalNumBytes() );
    for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
      // writer->line("Input trace header #%2d: %d  '%s'",ihdr,vars->hdrId[currentFileNumberToSet][ihdr], hdef->headerName(ihdr).c_str());
      if( hdrIdPtr[ihdr] != mod_input::NOT_AVAILABLE ) {
        switch( hdrTypePtr[ihdr] ) {
          case TYPE_INT:
            trcHdr->setIntValue( ihdr, trcHdrInPtr->intValue( hdrIdPtr[ihdr] ) );
            break;
          case TYPE_FLOAT:
            trcHdr->setFloatValue( ihdr, trcHdrInPtr->floatValue( hdrIdPtr[ihdr] ) );
            break;
          case TYPE_DOUBLE:
            trcHdr->setDoubleValue( ihdr, trcHdrInPtr->doubleValue( hdrIdPtr[ihdr] ) );
            break;
          case TYPE_INT64:
            trcHdr->setInt64Value( ihdr, trcHdrInPtr->int64Value( hdrIdPtr[ihdr] ) );
            break;
          case TYPE_STRING:
            trcHdr->setStringValue( ihdr, trcHdrInPtr->stringValue( hdrIdPtr[ihdr] ) );
            break;
        }
      }
    }
  }

  trcHdr->setIntValue( vars->hdrId_fileno, currentFileNumberToSet+1 );

  if( vars->mergeOption == mod_input::MERGE_TRACE || vars->mergeOption == mod_input::MERGE_HDR_COPY ) {
    vars->currentFile += 1;
    if( vars->currentFile == vars->numFiles ) vars->currentFile = 0;
  }

  vars->traceCounter += 1;
  return;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_( csParamDef* pdef ) {
  pdef->setModule( "INPUT", "Input SeaSeis data", "Reads SeaSeis formatted data from disk file, file extension '.cseis'" );

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "ntraces", "Number of traces to read in", NUM_VALUES_FIXED,
                  "Input of traces will stop when all traces have been read in, or if the number of traces specified has been reached. Traces will not be filled up to the specified range");
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of traces to read in" );

  pdef->addParam( "merge", "Method to merge traces from multiple input files", NUM_VALUES_FIXED);
  pdef->addValue( "all", VALTYPE_OPTION );
  pdef->addOption( "all", "Read in all traces of file 1, then all traces of file 2 etc" );
  pdef->addOption( "trace", "Read in one trace per input file, then repeat until first file has been fully read in. Stop then." );
  pdef->addOption( "header", "Merge traces from input files by header value. Assumes data has been previously sorted.", "Requires header name to be specified" );
  pdef->addOption( "hdr_copy_file1", "Same as option 'trace', except that the trace header of each trace is copied as-is from the first input file to all other input files." );

  pdef->addParam( "header_merge", "Trace header name for merged input (parameter 'order', option 'header').", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "nsamples", "Number of samples to read in", NUM_VALUES_FIXED,
                  "If number of samples in input data set is smaller, traces will be filled with zeros.");
  pdef->addValue( "", VALTYPE_NUMBER, "Number of samples to read in" );

  pdef->addParam( "override_sampleint", "Override sample interval", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Sample interval ([ms] [Hz]...). Specify '0' to NOT override sample interval" );
  
  pdef->addParam( "header", "Name of trace header used for trace selection", NUM_VALUES_FIXED, "Use in combination with user parameter 'select'. NOTE: With the current Seaseis disk data format, selecting traces on input is typically slower than reading in all traces and making the trace selection later on, e.g. by using module 'SELECT'" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addParam( "select", "Selection of header values", NUM_VALUES_FIXED, "Only traces which fit the trace value selection will be read in. Use in combination with user parameter 'header'" );
  pdef->addValue( "", VALTYPE_STRING, "Selection string. See documentation for more detailed description of selection syntax" );

  pdef->addParam( "sort", "Sort input data traces specified in user parameter 'header' and 'select'?", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not sort input data traces" );
  pdef->addOption( "increasing", "Sort input data traces by specified trace header, in increasing order" );
  pdef->addOption( "decreasing", "Sort input data traces by specified trace header, in decreasing order" );
  pdef->addValue( "simple", VALTYPE_OPTION );
  pdef->addOption( "simple", "Simplest sort method. Fastest for small and partially pre-sorted data sets" );
  pdef->addOption( "tree", "Tree sorting method. Most efficient for large, totally un-sorted data sets" );

  pdef->addParam( "ntraces_buffer", "Number of traces to buffer", NUM_VALUES_FIXED,
                  "Reading a large number of traces at once may enhance performance, but requires more memory" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to buffer when reading" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_input_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_input::VariableStruct* vars = reinterpret_cast<mod_input::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_input_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_input::VariableStruct* vars = reinterpret_cast<mod_input::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->filenames != NULL ) {
    delete [] vars->filenames;
    vars->filenames = NULL;
  }
  if( vars->readers != NULL ) {
    for( int i = 0; i < vars->numFiles; i++ ) {
      delete vars->readers[i];
    }
    delete [] vars->readers;
    vars->readers = NULL;
  }
  if( vars->mergeHdrValues != NULL ) {
    delete [] vars->mergeHdrValues;
    vars->mergeHdrValues = NULL;
  }
  if( vars->hdrValueBlock != NULL ) {
    delete [] vars->hdrValueBlock;
    vars->hdrValueBlock = NULL;
  }
  if( vars->hdrValueBlock_dummy != NULL ) {
    delete [] vars->hdrValueBlock_dummy;
    vars->hdrValueBlock_dummy = NULL;
  }
  if( vars->trcHdrDef != NULL ) {
    for( int i = 0; i < vars->numFiles; i++ ) {
      if( vars->trcHdrDef[i] != NULL ) {
        delete vars->trcHdrDef[i];
      }
    }
    delete [] vars->trcHdrDef;
    vars->trcHdrDef = NULL;
  }
  if( vars->trcHdr != NULL ) {
    for( int i = 0; i < vars->numFiles; i++ ) {
      if( vars->trcHdr[i] != NULL ) {
        delete vars->trcHdr[i];
      }
    }
    delete [] vars->trcHdr;
    vars->trcHdr = NULL;
  }
  if( vars->hdrId != NULL && vars->hdrType != NULL ) {
    for( int i = 0; i < vars->numFiles; i++ ) {
      if( vars->hdrId[i] != NULL ) {
        delete [] vars->hdrId[i];
      }
      if( vars->hdrType[i] != NULL ) {
        delete [] vars->hdrType[i];
      }
    }
    delete [] vars->hdrId;
    delete [] vars->hdrType;
    vars->hdrId   = NULL;
    vars->hdrType = NULL;
  }
  if( vars->hdrIsEqual != NULL ) {
    delete [] vars->hdrIsEqual;
    vars->hdrIsEqual = NULL;
  }
  if( vars->fileIndexList != NULL ) {
    delete vars->fileIndexList;
    vars->fileIndexList = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_input_( csParamDef* pdef ) {
  params_mod_input_( pdef );
}
extern "C" void _init_mod_input_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_( env, writer );
}
extern "C" void _exec_mod_input_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_( env, writer );
}

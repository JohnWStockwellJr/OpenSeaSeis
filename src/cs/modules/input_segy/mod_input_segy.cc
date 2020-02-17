/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstring>
#include <cmath>
 
#include "cseis_includes.h"
#include "csSegyReader.h"
#include "csSegyTraceHeader.h"
#include "csSegyHeaderInfo.h"
#include "csSegyDefines.h"
#include "csSegyBinHeader.h"
#include "csStandardHeaders.h"
#include "csFileUtils.h"
#include "csSort.h"
#include "csGeolibUtils.h"
#include "csIOSelection.h"
#include "csSortManager.h"
 
using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;
 
namespace mod_input_segy {
  struct CharInfo {
    int line;
    int byte1;
    int byte2;
    int numBytes;
    int hdrId;
    int hdrType;
    int intValue;
    int doubleValue;
  };
  struct VariableStruct {
    long traceCounter;
    double startTimeUNIXsec;
    csSegyReader* segyReader;
    int*  hdrIndexSegy;
    type_t* hdrTypeSegy;
    int   nsampMAX;
    int   nTracesToRead;
    int hdr_mapping;
    bool dumpTrchdr;
    int numFiles;
    int currentFile;
    char** filenames;
    int jobID;    // Job ID from binary header
    int reelNum;   // Reel number from binary header
    int lineNum;   // Line number from binary header
    int hdrId_jobID;
    int hdrId_lineNum;
    int hdrId_reelNum;
    int hdrId_fileno;
    int hdrId_time_samp1;
    int hdrId_time_samp1_us;
    int hdrId_time_year;
    int hdrId_time_day;
    int hdrId_time_hour;
    int hdrId_time_min;
    int hdrId_time_sec;
    int hdrId_time_msec;
    bool atEOF;
    csSegyReader::SegyReaderConfig config;
    csSegyHdrMap* hdrMap;
    bool override_year;
    int year;
 
    bool isHdrSelection;
    int sortOrder;
    int sortMethod;
    std::string selectionText;
    std::string selectionHdrName;
    bool readFirstLast;
    CharInfo* charInfo;
  };
  void dumpFileHeaders( csLogWriter* writer, csSegyReader* reader, int hdr_mapping );
  void extractInfoFromCharHdr( char const* charHdr, CharInfo* charInfo );
}
using mod_input_segy::VariableStruct;
 
//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_segy_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
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
 
  vars->startTimeUNIXsec = 0.0;
  vars->traceCounter = 0;
  vars->segyReader   = NULL;
  vars->hdrIndexSegy = NULL;
  vars->hdrTypeSegy = NULL;
  vars->nsampMAX = 0;
  vars->nTracesToRead = 0;
  vars->hdr_mapping = csSegyHdrMap::SEGY_STANDARD;
  vars->hdrMap = NULL;
  vars->dumpTrchdr = false;
  vars->numFiles = 0;
  vars->currentFile = 0;
  vars->filenames = NULL;
  vars->jobID     = 0;    // Job ID from binary header
  vars->reelNum   = 0;   // Reel number from binary header
  vars->lineNum   = 0;   // Line number from binary header
  vars->hdrId_jobID   = -1;
  vars->hdrId_lineNum = -1;
  vars->hdrId_reelNum = -1;
  vars->hdrId_fileno  = -1;
 
  vars->override_year = false;
  vars->year = 1970;
 
  vars->hdrId_time_samp1  = -1;
  vars->hdrId_time_samp1_us = -1;
  vars->hdrId_time_year   = -1;
  vars->hdrId_time_day    = -1;
  vars->hdrId_time_hour   = -1;
  vars->hdrId_time_min    = -1;
  vars->hdrId_time_sec    = -1;
  vars->hdrId_time_msec   = -1;
 
  vars->atEOF = false;
  vars->isHdrSelection = false;
  vars->sortOrder = cseis_geolib::csIOSelection::SORT_NONE;
  vars->sortMethod = cseis_geolib::csSortManager::SIMPLE_SORT;
  vars->selectionText = "";
  vars->selectionHdrName = "";

  vars->readFirstLast = false;
  vars->charInfo = NULL;

  //------------------------------------------------
 
  std::string headerName;
  std::string yesno;
  bool printEbcdic = true;
  int numTracesBuffer = 0;
  bool searchSubDirs  = false;
 
  if( edef->isDebug() ) writer->line("Starting init phase of INPUT_SEGY...");
 
  /*
   * Trace selection based on header value
   * Note that with the current version of the input data format there is no performance increase, rather a decrease compared
   * to reading in all traces and selecting later on using the module $SELECT.
   * Best to use this functionality only for sorting of input data.
   */
  if( param->exists( "header_select" ) ) {
    vars->isHdrSelection = true;
    numTracesBuffer = 1;  // Force single trace buffer.
    csVector<std::string> valueList;
    std::string text;
    param->getAll( "header_select", &valueList );
    if( valueList.size() > 1 ) {
      writer->error("Currently, trace selection is only supported for one trace header. Number of supplied header names: %d", valueList.size());
    }
    vars->selectionHdrName = valueList.at(0);
    valueList.clear();
    param->getString( "select", &vars->selectionText );
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
      writer->warning("Selecting traces on input using user parameters 'header' & 'select' is typically slower than reading in all traces and performing the selection afterwards, e.g. by using module 'SELECT'. It is recommended to use input selection only when traces shall be sorted on input");
    }
  }
 
  //------------------------------------
  if( param->exists( "year" ) ) {
    param->getInt("year",&vars->year);
    vars->override_year = true;
  }
 
  //------------------------------------
  int numSingleFiles = param->getNumLines( "filename" );
  csVector<string> fileList;
  if( param->exists("directory") ) {
    int numValues = param->getNumValues("directory");
    string directory;
    param->getString("directory", &directory, 0);
    string extension = "segy";
    if( numValues > 1 ) {
      param->getString("directory", &extension, 1);
      if( numValues > 2 ) {
        string text;
        param->getString("directory", &text, 2);
        if( !text.compare("yes") ) {
          searchSubDirs = true;
        }
        else if( !text.compare("no") ) {
          searchSubDirs = false;
        }
        else {
          writer->error("Unknown option: %s", text.c_str() );
        }
      }
    }
    csFileUtils::retrieveFiles( directory, extension, &fileList, searchSubDirs, writer->getFile() );
 
    // Sort file names
    string* names = new string[fileList.size()];
    for( int i = 0; i < fileList.size(); i++ ) {
      names[i] = fileList.at(i);
    }
    csSort<std::string>().simpleSort( names, fileList.size() );
    for( int i = 0; i < fileList.size(); i++ ) {
      fileList.set( names[i], i );
    }
    delete [] names;
  }
 
  if( numSingleFiles <= 0 && fileList.size() == 0 ) {
    writer->error("No filename specified.");
  }
  vars->numFiles = numSingleFiles + fileList.size();
  vars->filenames = new char*[vars->numFiles];
 
  for( int i = numSingleFiles; i < vars->numFiles; i++ ) {
    string s1 = fileList.at(i-numSingleFiles);
    int length = (int)s1.size();
    vars->filenames[i] = new char[length+1];
    memcpy(vars->filenames[i],s1.c_str(),length);
    vars->filenames[i][length] = '\0';
  }
 
  for( int i = 0; i < numSingleFiles; i++ ) {
    string filename;
    param->getStringAtLine( "filename", &filename, i );
    int length = (int)filename.length();
    vars->filenames[i] = new char[length+1];
    memcpy(vars->filenames[i],filename.c_str(),length);
    vars->filenames[i][length] = '\0';
  }
 
  for( int i = 0; i < vars->numFiles; i++ ) {
    writer->line("Input file #%4d: %s", i, vars->filenames[i]);
  }
  writer->line("");
 
  //----------------------------------------------------
 
  int numSamplesOut = 0;
  bool overrideNumSamples = false;
  if( param->exists("nsamples" ) ) {
    param->getInt( "nsamples", &numSamplesOut, 0 );
    if( param->getNumValues("nsamples") > 1 ) {
      string text;
      param->getString( "nsamples", &text, 1 );
      if( !text.compare("yes") ) {
        overrideNumSamples = true;
      }
      else if( !text.compare("no") ) {
        overrideNumSamples = false;
      }
      else {
        writer->error("Unknown option: %s", text.c_str() );
      }
    }
  }
 
  bool overrideSampleInt = false;
  float sampleInt = 0.0;
  if( param->exists("sample_int") ) {
    param->getFloat( "sample_int", &sampleInt );
    overrideSampleInt = true;
    if( sampleInt <= 0 ) {
      writer->error("Wrong sample interval specified: %f", sampleInt);
    }
  }
 
  if( param->exists("ntraces") ) {
    param->getInt( "ntraces", &vars->nTracesToRead );
    if( vars->nTracesToRead < 0 ) {
      vars->nTracesToRead = 0;
    }
  }
 
  //---------------------------------------------------------------------
  bool isSUFormat = false;
  if( param->exists("su_format") ) {
    std::string text;
    param->getString( "su_format", &text );
    if( !text.compare("yes") ) {
      isSUFormat = true;
    }
    else if( !text.compare("no") ) {
      isSUFormat = false;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  else {
    int length = (int)strlen(vars->filenames[0]);
    if( length >= 3 ) {
      char* ptr = &vars->filenames[0][length-3];
      if( !strcmp(ptr,".su") || !strcmp(ptr,".SU") || !strcmp(ptr,".Su") || !strcmp(ptr,".sU") ) {
        isSUFormat = true;
      }
      else {
        isSUFormat = false;
      }
    }
  }
 
  if( param->exists("hdr_map") ) {
    string hdr_mappingName;
    param->getString( "hdr_map", &hdr_mappingName );
    if( !hdr_mappingName.compare("obc") ) {
      vars->hdr_mapping = csSegyHdrMap::SEGY_OBC;
    }
    else if( !hdr_mappingName.compare("standard") ) {
      vars->hdr_mapping = csSegyHdrMap::SEGY_STANDARD;
    }
    else if( !hdr_mappingName.compare("node") ) {
      vars->hdr_mapping = csSegyHdrMap::SEGY_NODE;
    }
    else if( !hdr_mappingName.compare("das") ) {
      vars->hdr_mapping = csSegyHdrMap::SEGY_DAS;
    }
    else if( !hdr_mappingName.compare("node_old") ) {
      vars->hdr_mapping = csSegyHdrMap::SEGY_NODE_OLD;
    }
    else if( !hdr_mappingName.compare("send") ) {
      vars->hdr_mapping = csSegyHdrMap::SEGY_SEND;
    }
    else if( !hdr_mappingName.compare("armss") ) {
      vars->hdr_mapping = csSegyHdrMap::SEGY_ARMSS;
    }
    else if( !hdr_mappingName.compare("psegy") ) {
      vars->hdr_mapping   = csSegyHdrMap::SEGY_PSEGY;
      numTracesBuffer = 1;
    }
    else if( !hdr_mappingName.compare("su") ) {
      vars->hdr_mapping   = csSegyHdrMap::SEGY_SU;
    }
    else if( !hdr_mappingName.compare("su_standard") ) {
      vars->hdr_mapping   = csSegyHdrMap::SEGY_SU_ONLY;
    }
    else if( !hdr_mappingName.compare("su_both") ) {
      vars->hdr_mapping   = csSegyHdrMap::SEGY_SU_BOTH;
    }
    else if( !hdr_mappingName.compare("none") ) {
      vars->hdr_mapping   = csSegyHdrMap::NONE;
    }
    else {
      writer->error("Option for user parameter 'hdr_map' not recognised: %s", hdr_mappingName.c_str());
    }
  }
  else if( isSUFormat ) { // Set default map for SU input data
    vars->hdr_mapping   = csSegyHdrMap::SEGY_SU;
  }
 
  if( param->exists("filename_hdrmap") ) {
    string filename_hdrmap;
    param->getString("filename_hdrmap", &filename_hdrmap );
    vars->hdrMap = new csSegyHdrMap( vars->hdr_mapping, false, filename_hdrmap );
  }
  else {
    vars->hdrMap = new csSegyHdrMap( vars->hdr_mapping, false );
  }
 
  //---------------------------------------------------------------------
 
  if( param->exists("print") ) {
    param->getString( "print", &yesno );
    yesno = toLowerCase( yesno );
    if( !yesno.compare( "yes" ) ) {
      printEbcdic = true;
    }
    else if( !yesno.compare( "no" ) ) {
      printEbcdic = false;
    }
    else {
      writer->line("%s: Parameter value not recognized: %s", edef->moduleName().c_str(), yesno.c_str() );
      env->addError();
    }
  }
  else {
    printEbcdic = true;
  }
 
  if( param->exists("dump_trchdr") ) {
    param->getString( "dump_trchdr", &yesno );
    yesno = toLowerCase( yesno );
    if( !yesno.compare( "yes" ) ) {
      vars->dumpTrchdr = true;
    }
    else if( !yesno.compare( "no" ) ) {
      vars->dumpTrchdr = false;
    }
    else {
      writer->line("%s: Parameter value not recognized: %s", edef->moduleName().c_str(), yesno.c_str() );
      env->addError();
    }
  }
  else {
    vars->dumpTrchdr = false;
  }
 
  if( param->exists( "ntraces_buffer" ) ) {
    param->getInt( "ntraces_buffer", &numTracesBuffer );
    if( numTracesBuffer < 0 || numTracesBuffer > 9999999 ) {
      writer->warning("Number of buffered traces out of range (=%d). Changed to default.", numTracesBuffer);
      numTracesBuffer = 0;
    }
  }
  bool autoscale_hdrs = true;
  if( param->exists("auto_scale") ) {
    string text;
    param->getString( "auto_scale", &text );
    if( !text.compare("no") ) {
      autoscale_hdrs = false;
    }
  }
 
  bool overrideHdrs = false;
  if( param->exists("hdr_duplicate") ) {
    string text;
    param->getString( "hdr_duplicate", &text );
    if( !text.compare("override") ) {
      overrideHdrs = true;
    }
    else if( !text.compare("abort") ) {
      overrideHdrs = false;
    }
    else if( !text.compare("both") ) {
      overrideHdrs = false;
    }
    else {
      writer->error("Unknown option '%s'", text.c_str() );
    }
  }
 
  bool isEBCDIC = true;
  if( param->exists("charhdr_format") ) {
    string text;
    param->getString( "charhdr_format", &text );
    if( !text.compare("ascii") ) {
      isEBCDIC = false;
    }
    else if( !text.compare("ebcdic") ) {
      isEBCDIC = true;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  bool rev_byte_order_hdr  = false;
  bool rev_byte_order_data = false;
  if( param->exists("reverse_byte_order") ) {
    string text;
    param->getString( "reverse_byte_order", &text, 0 );
    if( !text.compare("yes") ) {
      rev_byte_order_hdr  = true;
      rev_byte_order_data = true; // By default, reverse byte order for both header & data
    }
    if( param->getNumValues("reverse_byte_order") > 1 ) {
      param->getString( "reverse_byte_order", &text, 1 );
      if( !text.compare("yes") ) {
        rev_byte_order_data = true;
      }
    }
  }
 
  int data_format = cseis_segy::AUTO;
  if( param->exists("data_format") ) {
    string text;
    param->getString( "data_format", &text );
    if( !text.compare("auto") ) {
      data_format = cseis_segy::AUTO;
    }
    else if( !text.compare("ieee") ) {
      data_format = cseis_segy::DATA_FORMAT_IEEE;
    }
    else if( !text.compare("ibm") ) {
      data_format = cseis_segy::DATA_FORMAT_IBM;
    }
    else {
      writer->error("Option for user parameter 'data_format' not recognised: %s", text.c_str());
    }
  }
 
  if( vars->nTracesToRead <= 0 ) vars->nTracesToRead = -1;  // Do not bother how many traces, read in all
  //----------------------------------------------------
  // Check that all files exist
  //
 
  bool isError = false;
  for( int i = 0; i < vars->numFiles; i++ ) {
    FILE* file;
    file = fopen( vars->filenames[i], "rb" );
    if( file == NULL ) {
      isError = true;
      writer->line("Input file '%s' does not exist\n", vars->filenames[i] );
    }
    else {
      fclose( file );
    }
  }
  if( isError ) {
    writer->error("Could not find input file(s).\n" );
  }
  vars->currentFile = 0;
 
 
 
  //-------------------------------------------------------
  // Set non-standard SEGY trace headers
  //
 
  int nLines = param->getNumLines( "header" );
  cseis_geolib::csVector<std::string> valueList;
  for( int iline = 0; iline < nLines; iline++ ) {
    param->getAll( "header", &valueList, iline );
    int nValues = valueList.size();
    if( nValues < 3 ) {
      writer->error("Incorrect number of parameter values for parameter 'header'. Expected number of parameter values: At least 3, found: %d",
                    nValues);
    }
    headerName  = valueList.at(0);
    int byteLoc = atoi( valueList.at(1).c_str() ) - 1;  // -1 to reduce byte location to start at 0
    std::string typeNameIn  = valueList.at(2);
    std::string typeNameOut = "int";
    if( valueList.size() > 3 ) {
      typeNameOut = toLowerCase( valueList.at(3) );
    }
    std::string desc        = "";
    if( valueList.size() > 4 ) desc = valueList.at(4);
    if( byteLoc <= 0 || byteLoc > 240 ) {
      throw( csException("Inconsistent byte location for trace header '%s': %d. Must be number between 1-240",
                         headerName.c_str(), byteLoc) );
    }
 
    if( !vars->hdrMap->addHeader( byteLoc, typeNameIn, typeNameOut, headerName, desc ) ) {
      if( overrideHdrs && vars->hdrMap->removeHeader( headerName ) ) {
        if( !vars->hdrMap->addHeader( byteLoc, typeNameIn, typeNameOut, headerName, desc ) ) {
          writer->error( "Cannot add trace header '%s'. Header already exists.\n", headerName.c_str() );
        }
      }
      else if( !overrideHdrs ) {
        writer->line("Error: Cannot add user defined header '%s': Trace header is already defined in specified default header map.", headerName.c_str());
        writer->line("    Termination of job can be avoided by specifying parameter 'hdr_duplicate  override'.");
        env->addError();
      }
      else {
        writer->error( "Unknown error occurred while auto-removing previously defined trace header '%s'. Program bug? Does header already exist?\n", headerName.c_str() );
      }
    }
  }
 
  if( param->exists("dump_hdrmap") ) {
    param->getString( "dump_hdrmap", &yesno );
    yesno = toLowerCase( yesno );
    if( !yesno.compare( "yes" ) ) {
      writer->line(" *** Dump of SEGY trace header map, including user specified non-standard headers ***");
      vars->hdrMap->dump( writer->getFile() );
      writer->line("");
    }
    else if( !yesno.compare( "no" ) ) {
      // Nothing
    }
    else {
      writer->line("%s: Parameter value not recognized: %s", edef->moduleName().c_str(), yesno.c_str() );
      env->addError();
    }
  }
 
  if( param->exists("read_first_last_only") ) {
    std::string text;
    param->getString("read_first_last_only", &text);
    if( !text.compare("yes") ) {
      vars->readFirstLast = true;
    }
  } 

  //----------------------------------------------------
  //
  vars->hdrIndexSegy = NULL;
  vars->hdrTypeSegy  = NULL;
 
  vars->segyReader = NULL;
  vars->config.numTracesBuffer       = numTracesBuffer;
  vars->config.segyHeaderMapping     = vars->hdr_mapping;
  vars->config.reverseByteOrderHdr   = rev_byte_order_hdr;
  vars->config.reverseByteOrderData  = rev_byte_order_data;
  vars->config.autoscaleHdrs         = autoscale_hdrs;
  vars->config.overrideSampleFormat  = data_format;
  vars->config.isSUFormat            = isSUFormat;
  vars->config.enableRandomAccess    = vars->isHdrSelection || vars->readFirstLast;
  try {
    vars->segyReader = new csSegyReader( vars->filenames[vars->currentFile], vars->config, vars->hdrMap );
  }
  catch( csException& e ) {
    vars->segyReader = NULL;
    writer->error("Error when opening SEGY file '%s'.\nSystem message: %s", vars->filenames[vars->currentFile], e.getMessage() );
  }

  if( param->exists("var_trace_len") ) {
    std::string text;
    param->getString("var_trace_len", &text);
    if( !text.compare("yes") ) {
      vars->segyReader->setVariableTraceLength();
    }
  } 
 
  //  writer->line("Number of traces")  vars->segyReader->numTracesCapacity();
  vars->segyReader->setCharHdrFormat( isEBCDIC );
 
 
  // Create all input trace headers
  int nHeaders = vars->segyReader->numTraceHeaders();
  vars->hdrIndexSegy = new int[nHeaders];
  vars->hdrTypeSegy  = new type_t[nHeaders];
  for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
    csSegyHeaderInfo const* info = vars->segyReader->header(ihdr);
    if( !hdef->headerExists(info->name) ) {
      vars->hdrTypeSegy[ihdr]  = info->outType;
      if( vars->hdrTypeSegy[ihdr] != TYPE_STRING ) {
        vars->hdrIndexSegy[ihdr] = hdef->addHeader( vars->hdrTypeSegy[ihdr], info->name, info->description );
      }
      else { // TYPE_STRING
        vars->hdrIndexSegy[ihdr] = hdef->addHeader( vars->hdrTypeSegy[ihdr], info->name, info->description, info->byteSize );
      }
    }
    else {
      vars->hdrIndexSegy[ihdr] = hdef->headerIndex( info->name );
      vars->hdrTypeSegy[ihdr]  = hdef->headerType( info->name );
    }
  }
  // Add headers that will be filled from binary header
  vars->hdrId_jobID   = hdef->addHeader( TYPE_INT, "jobID", "Job ID");
  vars->hdrId_reelNum = hdef->addStandardHeader( "tapeno" );
  vars->hdrId_lineNum = hdef->addHeader( TYPE_INT, "linenum", "Line number");
  if( !hdef->headerExists(HDR_FILENO.name) ) {
    hdef->addStandardHeader( HDR_FILENO.name );
  }
  vars->hdrId_fileno = hdef->headerIndex(HDR_FILENO.name);
 
  if( param->exists("hdr_ens") ) {
    param->getString( "hdr_ens", &headerName );
    if( !hdef->headerExists( headerName ) ) {
      writer->error("Specified ensemble header does not exist: '%s'.", headerName.c_str() );
    }
    type_t type_hdr_ens = hdef->headerType( headerName );
    if( type_hdr_ens != TYPE_INT && type_hdr_ens != TYPE_FLOAT && type_hdr_ens != TYPE_DOUBLE ) {
      writer->error("Ensemble header can only be of number type, not a string or array type.");
    }
    shdr->setEnsembleKey( headerName, 0 );
  }
 
  try {
    if( overrideNumSamples ) {
      vars->segyReader->initialize(numSamplesOut);
    }
    else {
      vars->segyReader->initialize();
    }
  }
  catch( csException& e ) {
    writer->line("\n\n\n");
    mod_input_segy::dumpFileHeaders( writer, vars->segyReader, vars->hdr_mapping );
    writer->line("\n\n");
    vars->segyReader = NULL;
    writer->error("Error when initializing SEGY reader object.\nSystem message: %s", e.getMessage() );
  }
 
  //--------------------------------------------------------------------------------
  // Trace selection/sorting based on input trace header
  // ...important to do the following after all headers have been set for hdef and all other input files
  //
  if( vars->isHdrSelection ) {
    hdef->resetByteLocation();  // This will otherwise be done by base system AFTER init phase
    if( !hdef->headerExists( vars->selectionHdrName ) ) {
      writer->error("Selection trace header '%s' is not defined in input file '%s'", vars->selectionHdrName.c_str(), vars->filenames[vars->currentFile]);
    }
    type_t selectionHdrType = hdef->headerType( vars->selectionHdrName );
    if( !hdef->headerExists( vars->selectionHdrName ) || hdef->headerType( vars->selectionHdrName ) != selectionHdrType ) {
      writer->error("Selection trace header '%s' is not defined in input file, or has different type",
                    vars->selectionHdrName.c_str(), vars->filenames[vars->currentFile]);
    }
    bool success = vars->segyReader->setSelection( vars->selectionText, vars->selectionHdrName, vars->sortOrder, vars->sortMethod );
    if( !success ) {
      writer->error("Error occurred when intializing header selection for input file '%s'.\n --> No input traces found that match specified selection '%s' for header '%s'.\n",
                    vars->filenames[vars->currentFile], vars->selectionText.c_str(), vars->selectionHdrName.c_str() );
    }
  }
 
  //----------------------------------------
  // Log dump
  //
  // Dump char header and bin header info to log file...
 
  cseis_geolib::csSegyBinHeader const* binHdr = vars->segyReader->binHdr();
 
  if( param->exists("charhdr_extract") ) {
    vars->charInfo = new mod_input_segy::CharInfo();
    vars->charInfo->intValue    = 0;
    vars->charInfo->doubleValue = 0;

    param->getString( "charhdr_extract", &headerName, 0 );
    if( !hdef->headerExists( headerName ) ) {
      writer->error("Specified header does not exist: '%s'.", headerName.c_str() );
    }
    vars->charInfo->hdrId   = hdef->headerIndex( headerName );
    vars->charInfo->hdrType = hdef->headerType( headerName );
    param->getInt( "charhdr_extract", &vars->charInfo->line, 1 );
    param->getInt( "charhdr_extract", &vars->charInfo->byte1, 2 );
    param->getInt( "charhdr_extract", &vars->charInfo->byte2, 3 );
    int numLines = 40;
    int numLettersPerLine = 80;
    if( vars->charInfo->line < 1 || vars->charInfo->line > numLines ) writer->error("Wrong line number (=%d) provided to extract char header information. Must be in the range 1-%d", vars->charInfo->line, numLines);
    if( vars->charInfo->byte1 < 1 || vars->charInfo->byte1 > numLettersPerLine ) writer->error("Wrong byte number (=%d) provided to extract char header information. Must be in the range 1-%d", vars->charInfo->byte1, numLettersPerLine);
    if( vars->charInfo->byte2 < 1 || vars->charInfo->byte2 > numLettersPerLine ) writer->error("Wrong byte number (=%d) provided to extract char header information. Must be in the range 1-%d", vars->charInfo->byte2, numLettersPerLine);
    if( vars->charInfo->byte2 <= vars->charInfo->byte1 ) writer->error("Inconsistent byte numbers provided to extract char header information. byte1 %d >= byte2 %d", vars->charInfo->byte1, vars->charInfo->byte2);
    vars->charInfo->numBytes = vars->charInfo->byte2 - vars->charInfo->byte1 + 1;

    mod_input_segy::extractInfoFromCharHdr( vars->segyReader->charHdrBlock(), vars->charInfo );
    if( vars->charInfo->hdrType == TYPE_DOUBLE || vars->charInfo->hdrType == TYPE_FLOAT ) {
      writer->line("Value extracted from char hdr: %f, trace header '%s'", vars->charInfo->doubleValue, headerName.c_str());
    }
    else {
      writer->line("Value extracted from char hdr: %d, trace header '%s'", vars->charInfo->intValue, headerName.c_str());
    }
  }

  if( printEbcdic && !isSUFormat ) {
    char const* charHdr = vars->segyReader->charHdrBlock();
    writer->line( "Segy EBCDIC header:\n" );
    char line[83];
    line[80] = '+';
    line[81] = '+';
    line[82] = '\0';
    for( int i = 0; i < cseis_segy::SIZE_CHARHDR; i++) {
      line[i%80] = charHdr[i];
      if( ((i+1) % 80) == 0 ) {
        writer->line( line );
      }
    }
    writer->line( "Segy binary header:\n" );
    if( vars->hdr_mapping != csSegyHdrMap::SEGY_PSEGY ) {
      binHdr->dump( writer->getFile() );
    }
  }
  if( overrideSampleInt ) {
    shdr->sampleInt = sampleInt;
  }
  else {
    shdr->sampleInt = vars->segyReader->sampleIntMS();
  }
 
  writer->line("");
  writer->line( "       File name:            %s", vars->filenames[vars->currentFile]);
  writer->line( "       Sample interval [ms]: %f", shdr->sampleInt );
  writer->line( "       Number of samples:    %d", vars->segyReader->numSamples() );
  writer->write("       Sample data format:   %d   (Supported formats are: 1:IBM, 2:32bit INT, 3:16bit INT, 5:IEEE)\n", vars->segyReader->dataSampleFormat() );
  if( data_format != cseis_segy::AUTO ) {
    writer->line("                     NOTE: The sample format will be overridden by the user specified value: %d", data_format );
  }
  //  else if( edef->isDebug() && data_format != cseis_segy::AUTO ) {
  //    writer->write(" ...overridden by user specified value: %d", data_format );
  //  }
  writer->line("");
  writer->line("\nMaximum number of buffered traces:    %d", vars->segyReader->numTracesCapacity() );

  //--------------------------------------------------
  // Sanity check of binary header: Sample interval & number of samples OK?
  //
  if( binHdr->sampleIntUS <= 0 ) {
    if( overrideSampleInt ) {
      writer->warning("Incorrect sample interval in binary header: %d us", binHdr->sampleIntUS);
    }
    else {
      writer->error("Incorrect sample interval in binary header: %d us", binHdr->sampleIntUS);
    }
  }
  if( binHdr->numSamples <= 0 ) {
    if( overrideNumSamples ) {
      writer->warning("Incorrect number of samples in binary header: %d us", binHdr->numSamples);
    }
    else {
      writer->error("Incorrect number of samples in binary header: %d us", binHdr->numSamples);
    }
  }

  if( numSamplesOut == 0 ) {
    shdr->numSamples = vars->segyReader->numSamples();
  }
  else {
    shdr->numSamples = numSamplesOut;
  }
  vars->jobID      = binHdr->jobID;
  vars->reelNum    = binHdr->reelNum;
  vars->lineNum    = binHdr->lineNum;
 
  if( vars->hdr_mapping != csSegyHdrMap::SEGY_SU_ONLY ) {
    vars->hdrId_time_year = hdef->headerIndex("time_year");
    vars->hdrId_time_day  = hdef->headerIndex("time_day");
    vars->hdrId_time_hour = hdef->headerIndex("time_hour");
    vars->hdrId_time_min  = hdef->headerIndex("time_min");
    vars->hdrId_time_sec  = hdef->headerIndex("time_sec");
  }
  else {
    vars->hdrId_time_year = hdef->headerIndex("year");
    vars->hdrId_time_day  = hdef->headerIndex("day");
    vars->hdrId_time_hour = hdef->headerIndex("hour");
    vars->hdrId_time_min  = hdef->headerIndex("minute");
    vars->hdrId_time_sec  = hdef->headerIndex("sec");
  }
  vars->hdrId_time_samp1    = hdef->headerIndex( HDR_TIME_SAMP1.name );
  vars->hdrId_time_samp1_us = hdef->headerIndex( HDR_TIME_SAMP1_US.name );
  if( vars->hdr_mapping == csSegyHdrMap::SEGY_PSEGY ) {
    vars->hdrId_time_msec = hdef->headerIndex("ps_time_samp1_ms");
  }
 
  vars->traceCounter = 0;
 
  vars->segyReader->freeCharBinHdr();  // Don't need these headers anymore -> free memory
 
  if( vars->dumpTrchdr ) {
    writer->line( "Segy trace header dump:\n" );
    int nHeaders = vars->segyReader->numTraceHeaders();
    writer->line( "... %d trace headers\n", nHeaders );
    vars->segyReader->getTrcHdrMap()->dump( writer->getFile() );
  }
}
 
//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_input_segy_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csTrace* trace = traceGather->trace(0);
 
  csExecPhaseDef*         edef = env->execPhaseDef;
 
  if( vars->atEOF ) {
    traceGather->freeAllTraces();
    return;
  }
  
  csSuperHeader const*    shdr = env->superHeader;
  csTraceHeader* trcHdr = trace->getTraceHeader();
  float* samples = trace->getTraceSamples();
 
  int numSamplesActual = MIN( vars->segyReader->numSamples(), shdr->numSamples );
  for( int isamp = numSamplesActual; isamp < shdr->numSamples; isamp++ ) {
    samples[isamp] = 0.0;
  }
  if( edef->isDebug() ) writer->line("INPUT SEGY...");
  if( vars->readFirstLast ) {
    if( vars->traceCounter == 1 ) {
      vars->segyReader->moveToTrace( vars->segyReader->numTraces()-1, 1 );
    }
  }
  if( (vars->nTracesToRead > 0 && vars->nTracesToRead == vars->traceCounter) ||
      !vars->segyReader->getNextTrace( (byte_t*)samples, numSamplesActual ) ) {
    if( vars->traceCounter == 0 ) {
      writer->warning("SEGY file '%s' does not contain any data trace, only char & bin headers.", vars->segyReader->filename() );
    }
    vars->currentFile += 1;
    if( vars->currentFile == vars->numFiles ) {
      vars->atEOF = true;
      traceGather->freeAllTraces();
      return;
    }
    delete vars->segyReader;
    vars->segyReader = NULL;
 
    try {
      vars->segyReader = new csSegyReader( vars->filenames[vars->currentFile], vars->config, vars->hdrMap );
      if( vars->isHdrSelection ) {
        bool success = vars->segyReader->setSelection( vars->selectionText, vars->selectionHdrName, vars->sortOrder, vars->sortMethod );
        if( !success ) {
          writer->error("Error occurred when intializing header selection for input file '%s'.\n --> No input traces found that match specified selection '%s' for header '%s'.\n",
                        vars->filenames[vars->currentFile], vars->selectionText.c_str(), vars->selectionHdrName.c_str() );
        }
      }
    }
    catch( csException& e ) {
      vars->segyReader = NULL;
      writer->error("Error when opening SEGY file '%s'.\nSystem message: %s", vars->filenames[vars->currentFile], e.getMessage() );
    }
    try {
      vars->segyReader->initialize();
      if( edef->isDebug() ) {
        writer->line("");
        writer->line("  File name:                 %s", vars->filenames[vars->currentFile]);
        writer->line("  Sample interval [ms]: %f", vars->segyReader->sampleIntMS() );
        writer->line("  Number of samples:         %d", vars->segyReader->numSamples() );
        writer->line("  Sample data format:  %d   (Supported formats are: 1:IBM, 5:IEEE, 2:32bit INT)", vars->segyReader->dataSampleFormat() );
        writer->line("\n...Maximum number of traces buffered in reader: %d", vars->segyReader->numTracesCapacity() );
      }
    }
    catch( csException& e ) {
      vars->segyReader = NULL;
      writer->error("Error when initializing SEGY reader object.\nSystem message: %s", e.getMessage() );
    }
    cseis_geolib::csSegyBinHeader const* binHdr = vars->segyReader->binHdr();
 
    float diff = fabs( shdr->sampleInt - (float)(binHdr->sampleIntUS)/1000.0 );
    if( diff > 0.001 ) {
      writer->error("Input SEGY files have different sample intervals: %f (1)  !=       %f (2)\n", shdr->sampleInt, (float)(binHdr->sampleIntUS)/1000.0 );
    }
    vars->jobID      = binHdr->jobID;
    vars->reelNum    = binHdr->reelNum;
    vars->lineNum    = binHdr->lineNum;
    if( vars->charInfo != NULL ) {
      mod_input_segy::extractInfoFromCharHdr( vars->segyReader->charHdrBlock(), vars->charInfo );
    }
    vars->segyReader->freeCharBinHdr();  // Don't need these headers anymore -> free memory
 
    numSamplesActual = MIN( vars->segyReader->numSamples(), shdr->numSamples );
    for( int isamp = numSamplesActual; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] = 0.0;
    }
 
    if( (vars->nTracesToRead > 0 && vars->nTracesToRead == vars->traceCounter) ||
        !vars->segyReader->getNextTrace( (byte_t*)samples, numSamplesActual ) ) {
      if( vars->traceCounter == 0 ) {
        writer->warning("SEGY file '%s' does not contain any data trace, only char & bin headers.", vars->segyReader->filename() );
      }
      vars->atEOF = true;
      traceGather->freeAllTraces();
      return;
    }
  }
  vars->traceCounter++;
 
  csSegyTraceHeader const* segyTrcHdr = vars->segyReader->getTraceHeader();
 
  int nHeaders = segyTrcHdr->numHeaders();
  for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
    int hdrIdOut = vars->hdrIndexSegy[ihdr];
    switch( vars->hdrTypeSegy[ihdr] ) {
    case TYPE_FLOAT:
      trcHdr->setFloatValue( hdrIdOut, segyTrcHdr->floatValue(ihdr) );
      if( edef->isDebug() ) writer->line( "Segy header %2d: %d   '%s'  %f", ihdr, segyTrcHdr->intValue(ihdr), segyTrcHdr->headerName(ihdr), segyTrcHdr->floatValue(ihdr) );
      break;
    case TYPE_DOUBLE:
      trcHdr->setDoubleValue( hdrIdOut, segyTrcHdr->doubleValue(ihdr) );
      if( edef->isDebug() ) writer->line( "Segy header %2d: %d   '%s'  %f", ihdr, segyTrcHdr->intValue(ihdr), segyTrcHdr->headerName(ihdr), segyTrcHdr->doubleValue(ihdr) );
      break;
    case TYPE_INT:
      trcHdr->setIntValue( hdrIdOut, segyTrcHdr->intValue(ihdr) );
      if( edef->isDebug() ) writer->line( "Segy header %2d: %d   '%s'  %d", ihdr, segyTrcHdr->intValue(ihdr), segyTrcHdr->headerName(ihdr), segyTrcHdr->intValue(ihdr) );
      break;
    case TYPE_INT64:
      trcHdr->setInt64Value( hdrIdOut, (csInt64_t)segyTrcHdr->intValue(ihdr) );
      if( edef->isDebug() ) writer->line( "Segy header %2d: %d   '%s'  %d", ihdr, segyTrcHdr->intValue(ihdr), segyTrcHdr->headerName(ihdr), segyTrcHdr->intValue(ihdr) );
      break;
    case TYPE_STRING:
      trcHdr->setStringValue( hdrIdOut, segyTrcHdr->stringValue(ihdr) );
      if( edef->isDebug() ) writer->line( "Segy header %2d: %d   '%s'  %s", ihdr, segyTrcHdr->intValue(ihdr), segyTrcHdr->headerName(ihdr), segyTrcHdr->stringValue(ihdr).c_str() );
      break;
    }
  }
  trcHdr->setIntValue( vars->hdrId_jobID, vars->jobID );
  trcHdr->setIntValue( vars->hdrId_reelNum, vars->reelNum );
  trcHdr->setIntValue( vars->hdrId_lineNum, vars->lineNum );
  trcHdr->setIntValue( vars->hdrId_fileno, vars->currentFile+1 );
 
  int year = trcHdr->intValue(vars->hdrId_time_year);
  if( vars->override_year ) {
    year = vars->year;
  }
  int day  = trcHdr->intValue(vars->hdrId_time_day);
  int hour = trcHdr->intValue(vars->hdrId_time_hour);
  int min  = trcHdr->intValue(vars->hdrId_time_min);
  int sec  = trcHdr->intValue(vars->hdrId_time_sec);
  int msec = 0;
  if( vars->hdr_mapping == csSegyHdrMap::SEGY_PSEGY ) {
    msec = trcHdr->intValue(vars->hdrId_time_msec);
  }
  csInt64_t startTime = csGeolibUtils::date2UNIXmsec( year, day, hour, min, sec, msec );
  int time_samp1_s  = (int)(startTime/1000LL);
  int time_samp1_us = (int)(startTime%1000LL)*1000;
  trcHdr->setIntValue( vars->hdrId_time_samp1, time_samp1_s );
  trcHdr->setIntValue( vars->hdrId_time_samp1_us, time_samp1_us );
  if( edef->isDebug() ) {
    writer->line("Start time: %lldms  =  %ds (x1000)  +  %dus (/1000)", startTime, time_samp1_s, time_samp1_us);
  }
 
  if( vars->charInfo != NULL ) {
    switch( vars->charInfo->hdrType ) {
    case TYPE_DOUBLE:
    case TYPE_FLOAT:
      trcHdr->setDoubleValue( vars->charInfo->hdrId, vars->charInfo->doubleValue );
      break;
    default:
      trcHdr->setIntValue( vars->charInfo->hdrId, vars->charInfo->intValue );
    }
  }

  if( vars->dumpTrchdr ) {
    writer->line( "\n ****** Segy trace header dump, trace #%d ******", vars->traceCounter );
    vars->segyReader->dumpTrcHdr( writer->getFile() );
  }
  return;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_segy_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_SEGY", "Input SEGY/SU data", "Reads SEGY/SU data from disk file.");
 
  pdef->addDoc("This module reads in disk files in SEGY format or Seismic Unix (SU) format. For SEGY, only revision 0 is fully supported. For revision 1, only some but not all features are supported.");
  pdef->addDoc("For example, there is no support for extended textual header blocks, or varying trace length and sample intervals.");
  pdef->addDoc("SU format is basically read in as SEGY rev 0 except there is no textual and binary header block.");
  pdef->addDoc("SU files are assumed to be in big endian format; if this is not true, manually swap the endian byte order with user parameter 'reverse_byte_order'.");
  pdef->addDoc("SU file data sample format is assumed to be IEEE.");
  pdef->addDoc("");
  pdef->addDoc("To read input data in SU format regardless of the (first) file's extension, specify user parameter 'su_format yes'. To read in additional trace headers named the same as in SU, specify 'hdr_map su_both'. In order to read in trace headers solely with the SU naming convention, specify 'hdr_map su_standard'" );
 
  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );
 
  pdef->addParam( "directory", "Name of directory to search", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Name of directory", "All files in this directory and its subdirectories will be read" );
  pdef->addValue( "segy", VALTYPE_STRING, "File extension", "Only files with the given file extension will be read in" );
  pdef->addValue( "no", VALTYPE_OPTION, "Search subdirectories" );
  pdef->addOption( "no", "Do not search subdirectories" );
  pdef->addOption( "yes", "Also search subdirectories for files" );
 
  pdef->addParam( "nsamples", "Number of samples to read in", NUM_VALUES_VARIABLE,
                  "If number of samples in input data set is smaller, traces will be filled with zeros. Set 0 to set number of samples from input data set.");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of samples to read in" );
  pdef->addValue( "no", VALTYPE_OPTION, "Override number of samples in SEGY bin header?" );
  pdef->addOption( "no", "Do not override number of samples in SEGY bin header" );
  pdef->addOption( "yes", "Override number of samples in SEGY bin header" );
 
  pdef->addParam( "sample_int", "Sample interval", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Sample interval [ms]" );
 
  pdef->addParam( "ntraces", "Number of traces to read in", NUM_VALUES_FIXED,
                  "Input of traces will stop when all traces have been read in, or if the number of traces specified has been reached. Traces will not be filled up to the specified range");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to read in" );
 
  pdef->addParam( "ntraces_buffer", "Number of traces to read into buffer at once", NUM_VALUES_FIXED,
                  "Reading in a large number of traces at once may enhance performance, but requires more memory" );
  pdef->addValue( "20", VALTYPE_NUMBER, "Number of traces to buffer" );
 
  pdef->addParam( "print", "Print EBCDIC & binary header to log file", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Print EBCDIC & binary header" );
  pdef->addOption( "no", "Do not print EBCDIC & binary header" );
 
  pdef->addParam( "dump_trchdr", "Dump trace header byte block to log file", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Dump trace header" );
  pdef->addOption( "no", "Dump trace header" );
 
  pdef->addParam( "dump_hdrmap", "Dump specified trace header map (=byte locations) to log file", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Dump trace header map" );
  pdef->addOption( "no", "Dump trace header map" );
 
  pdef->addParam( "hdr_ens", "Ensemble trace header name", NUM_VALUES_FIXED,
                  "Ensembles are defined by this key header" );
  pdef->addValue( "", VALTYPE_STRING, "Ensemble trace header name" );
 
 
  pdef->addParam( "hdr_map", "Pre-set SEGY trace header mapping", NUM_VALUES_FIXED, "Use option 'none' to avoid decoding any pre-set trace headers. Specify user parameter 'dump_hdrmap  yes' to make a printout of the specified pre-set map");
  pdef->addValue( "standard", VALTYPE_OPTION );
  pdef->addOption( "standard", "Standard SEGY header mapping, SEGY revision 1" );
  //  pdef->addOption( "ascii_file", "Trace header map definition is given in external ASCII file, see user parameter 'filename_hdrmap'" );
  pdef->addOption( "none", "Do not use any pre-set trace header map. All trace headers that shall be extracted are given either by user parameter 'header' (non-standard trace header) or 'filename_hdrmap' (external ASCII file)" );
  pdef->addOption( "obc", "OBC 4C header mapping", "In principle SEGY revision 1, plus additional headers" );
  pdef->addOption( "send", "SEND header mapping" );
  pdef->addOption( "armss", "ARMSS header mapping" );
  pdef->addOption( "psegy", "PSEGY header mapping" );
  pdef->addOption( "node", "Node header mapping" );
  pdef->addOption( "node_old", "Node header mapping (obsolete)" );
  pdef->addOption( "su", "Read in Seismic Unix (SU) trace headers in Seaseis standard naming when possible. All trace headers only defined in SU will be read in with standard SU naming");
  pdef->addOption( "su_standard", "Same as option 'su', but read in all trace headers with SU standard names only", "NOTE: This may cause problems with other Seaseis modules which assume standard Seaseis header naming");
  pdef->addOption( "su_both", "Read in trace headers in Seismic Unix (SU) standard naming AND in Seaseis standard naming", "This means some byte locations are mapped to two identical trace headers but with different names" );
 
  pdef->addParam( "filename_hdrmap", "File name containing trace header map setup", NUM_VALUES_FIXED,
                  "All trace headers already defined in the specified pre-set trace header map at the same byte locations will be overriden. If no pre-set trace headers are required, specify user parameter 'hdr_map  none'");
  pdef->addValue( "", VALTYPE_STRING, "Input file name containing trace header map definition",
                  "File format:  NAME  BYTELOC  TYPE_SEGY  TYPE_CEIS   DESC, where\n"\
                  "BYTELOC:    Byte location, starting at 1 for first byte,\n"\
                  "TYPE_SEGY:  Header type as stored in SEGY file, see user parameter 'header' for options,\n"\
                  "TYPE_CSEIS: Header type of Seaseis trace header, see user parameter 'header' for options,\n"\
                  "NAME:       Trace header name,\n"\
                  "DESC:       Description of Seaseis trace header, enclosed in double-quotes.");
 
  pdef->addParam( "header", "Set user specified header byte location etc", NUM_VALUES_VARIABLE);
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Byte location (starting at 1)" );
  pdef->addValue( "int", VALTYPE_OPTION, "Trace header input type", "This is the number type in which the trace header is stored/formatted in the input file's trace header");
  pdef->addOption( "int", "Integer (4 byte)" );
  pdef->addOption( "short", "Short (2 byte)" );
  pdef->addOption( "ushort", "Unsigned short (2 byte)" );
  pdef->addOption( "float", "Float (4 byte) IEEE" );
  pdef->addOption( "float_ibm", "Float (4 byte) (IBM)" );
  pdef->addOption( "4", "Integer (4 byte)", "...for backward compatibility" );
  pdef->addOption( "2", "Short (2 byte)", "...for backward compatibility" );
  pdef->addValue( "int", VALTYPE_OPTION, "Trace header output type", "This is the number type of the Seaseis output trace header" );
  pdef->addOption( "int", "Integer" );
  pdef->addOption( "float", "Single precision floating point" );
  pdef->addOption( "double", "Double precision floating point" );
  pdef->addValue( "", VALTYPE_STRING, "Header description", "Descriptive text for Seaseis trace header" );
 
  pdef->addParam( "hdr_duplicate", "How to treat user specified trace headers (parameter 'header') that are already defined in the chosen trace header map", NUM_VALUES_FIXED );
  pdef->addValue( "override", VALTYPE_OPTION );
  pdef->addOption( "override", "Write header to user specified byte location only, i.e. override standard byte location" );
  pdef->addOption( "abort", "Abort if user specified header is already defined in trace header map" );
 
  pdef->addParam( "auto_scale", "Automatically scale SEGY headers", NUM_VALUES_FIXED, "Automatically scale SEGY headers using SEGY scalar values found at standard byte positions, e.g. src xy and rcv xy header values");
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Automatically scale header values" );
  pdef->addOption( "no", "Do not apply scaler." );
 
  pdef->addParam( "reverse_byte_order", "Reverse byte order of input file (endian byte swapping)", NUM_VALUES_VARIABLE, "Setting this to 'yes' means that the input SEGY file is in the wrong byte order" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Reverse byte order of input file for binary and trace headers" );
  pdef->addOption( "no", "Do not reverse byte order for binary & trace headers" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Reverse byte order of input file for data samples" );
  pdef->addOption( "no", "Do not reverse byte order for data samples" );
 
  pdef->addParam( "data_format", "Override SEGY data sample format (given in binary header)", NUM_VALUES_FIXED);
  pdef->addValue( "auto", VALTYPE_OPTION );
  pdef->addOption( "auto", "Automatically read data sample format from binary header" );
  pdef->addOption( "ibm", "4 byte IBM floating point" );
  pdef->addOption( "ieee", "4 byte IEEE floating point" );
 
  pdef->addParam( "charhdr_format", "Format of char header", NUM_VALUES_FIXED);
  pdef->addValue( "ebcdic", VALTYPE_OPTION );
  pdef->addOption( "ebcdic", "EBCDIC" );
  pdef->addOption( "ascii", "ASCII" );
  //  pdef->addOption( "ebcdic_su", "EBCDIC code format used in Seismic Unix", "Differs slightly from EBCDIC code defined in SEGY rev1" );
 
  pdef->addParam( "su_format", "Is input data in SU (Seismic Unix) format?", NUM_VALUES_FIXED,
                  "Use this parameter to override the default format which is determined by the (first) file's extension: *.su or *.SU == SU file" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Force input data to be read in as SU format, independent of file extension" );
  pdef->addOption( "no", "Force input data to be read in as SEGY format, independent of file extension" );
 
  pdef->addParam( "header_select", "Name of trace header used for trace selection", NUM_VALUES_FIXED, "Use in combination with user parameter 'select'. NOTE: With the current Seaseis disk data format, selecting traces on input is typically slower than reading in all traces and making the trace selection later on, e.g. by using module 'SELECT'" );
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
 
  pdef->addParam( "year", "Override year found in SEGY trace header", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Year" );
 
  pdef->addParam( "var_trace_len", "Does input data have variable trace length?", NUM_VALUES_FIXED, "Requires the number of samples to be defined in each trace's header bytes 115-117" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "No" );
  pdef->addOption( "yes", "Yes" );

  pdef->addParam( "read_first_last_only", "read first and last trace only", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "No" );
  pdef->addOption( "yes", "Yes" );
}
 

//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_input_segy_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_input_segy::VariableStruct* vars = reinterpret_cast<mod_input_segy::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_input_segy_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_input_segy::VariableStruct* vars = reinterpret_cast<mod_input_segy::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->segyReader != NULL ) {
    delete vars->segyReader;
    vars->segyReader = NULL;
  }
  if( vars->hdrMap != NULL ) {
    delete vars->hdrMap;
    vars->hdrMap = NULL;
  }
  if( vars->filenames != NULL ) {
    for( int i = 0; i < vars->numFiles; i++ ) {
      if( vars->filenames[i] != NULL ) {
        delete [] vars->filenames[i];
      }
    }
    delete [] vars->filenames;
    vars->filenames = NULL;
  }
  if( vars->charInfo != NULL ) {
    delete vars->charInfo;
    vars->charInfo = NULL;
  }
  delete [] vars->hdrIndexSegy; vars->hdrIndexSegy = NULL;
  delete [] vars->hdrTypeSegy; vars->hdrTypeSegy = NULL;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_input_segy_( csParamDef* pdef ) {
  params_mod_input_segy_( pdef );
}
extern "C" void _init_mod_input_segy_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_segy_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_segy_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_segy_( env, writer );
}
extern "C" void _exec_mod_input_segy_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_segy_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_segy_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_segy_( env, writer );
}

void mod_input_segy::dumpFileHeaders( csLogWriter* writer, csSegyReader* reader, int hdr_mapping ) {
  cseis_geolib::csSegyBinHeader const* binHdr = reader->binHdr();
 
  char const* charHdr = reader->charHdrBlock();
  writer->line( "Segy EBCDIC header:\n" );
  char line[83];
  line[80] = '+';
  line[81] = '+';
  line[82] = '\0';
  for( int i = 0; i < cseis_segy::SIZE_CHARHDR; i++) {
    line[i%80] = charHdr[i];
    if( ((i+1) % 80) == 0 ) {
      writer->line( line );
    }
  }
  writer->line( "Segy binary header:\n" );
  if( hdr_mapping != csSegyHdrMap::SEGY_PSEGY ) {
    binHdr->dump( writer->getFile() );
  }
 
}

void mod_input_segy::extractInfoFromCharHdr( char const* charHdr, CharInfo* charInfo ) {
  int numLettersPerLine = 80;
  char* text = new char[charInfo->numBytes+1];
  memcpy( text, &charHdr[(charInfo->line-1)*numLettersPerLine+charInfo->byte1-1], charInfo->numBytes*sizeof(char) );
  text[charInfo->numBytes] = '\0';
  if( charInfo->hdrType == TYPE_DOUBLE || charInfo->hdrType == TYPE_FLOAT ) {
    charInfo->doubleValue = atof( text );
  }
  else {
    charInfo->intValue = atoi( text );
  }
  delete []text;
}

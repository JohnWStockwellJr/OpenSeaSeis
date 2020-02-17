/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

#include "csStandardSegdHeader.h"
#include "csStandardHeaders.h"
#include "csSegdReader.h"
#include "csSegdHeader_SEAL.h"
#include "csFileUtils.h"
#include "csGeolibUtils.h"
#include "csSort.h"
#include <string>
#include <cstring>
#include <cmath>

char ebcdic2char( short c );

using namespace cseis_system;
using namespace cseis_segd;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: INPUT_SEGD
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_input_segd {
  struct VariableStruct {
    csSegdReader* segdReader;
    int recSystemID;
    cseis_segd::commonFileHeaderStruct const* comFileHdr;
    cseis_segd::commonRecordHeaderStruct comRecordHdr;
    int traceCounter;
    int recordCounter;
    long totalTraceCounter;
    bool isAtEOF;

    int dumpHdrFlag;
    std::ofstream* file_dumpHeaders;
    bool dumpEssFileInfo;
    bool dumpTrcHeaders;
    int numTracesToRead;
    int numRecordsToRead;

    int hdrSet;
    int* hdrIndexSegd;
    type_t* hdrTypeSegd;
    
    int numFiles;
    int currentFile;
    char** filenames;

    int hdrId_ffid;
    int hdrId_trcno;
    int hdrId_chan;
    int hdrId_time_year;
    int hdrId_time_day;
    int hdrId_time_hour;
    int hdrId_time_min;
    int hdrId_time_sec;
    int hdrId_time_nano;
    int hdrId_source;
    int hdrId_sline;
    int hdrId_rline;
    int hdrId_rcv;
    int hdrId_serial;
    int hdrId_seq;
    int hdrId_sou_z;
    int hdrId_rec_wdep;
    int hdrId_sou_elev;
    int hdrId_rec_elev;
    int hdrId_rec_x;
    int hdrId_rec_y;
    int hdrId_sou_x;
    int hdrId_sou_y;
    int hdrId_trc_type;
    int hdrId_chan_set;
    int hdrId_sampint_us;
    int hdrId_nsamp;
    int hdrId_time_samp1;
    int hdrId_time_samp1_us;

    int hdrId_fileno;

    int* hdrId_extra;
    int numExtraHdrs;
    int chanSetIndexToRead;
    std::string filenameDumpHeaders;
    bool isFirstCall;

    bool hasCharHdr;
    int numCharBytes;
  };
}
using mod_input_segd::VariableStruct;

//void dumpEssentialHeaders_record( essentialSegdHeaders const* essHdrs, csLogWriter* writer );
//void dumpEssentialHeaders_trace( essentialSegdHeaders const* essHdrs, csLogWriter* writer );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_segd_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  //---------------------------------------------------------------
  // Initialise
  //

  vars->segdReader  = NULL;
  vars->recSystemID = cseis_segd::UNKNOWN;
  vars->comFileHdr  = NULL;
//    vars->comRecordHdr;
  vars->traceCounter   = 0;
  vars->recordCounter  = 0;
  vars->totalTraceCounter = 0;
  vars->isAtEOF = false;

  vars->dumpHdrFlag = cseis_segd::DUMP_NONE;
  vars->file_dumpHeaders = NULL;
  vars->dumpTrcHeaders   = false;
  vars->dumpEssFileInfo  = true;
  vars->numRecordsToRead = -1;
  vars->numTracesToRead  = -1;
  
  vars->hdrSet       = csStandardSegdHeader::HDR_MAPPING_STANDARD;
  vars->hdrIndexSegd = NULL;
  vars->hdrTypeSegd  = NULL;

  vars->numFiles    = 0;
  vars->currentFile = 0;
  vars->filenames   = NULL;

  vars->hdrId_ffid     = -1;
  vars->hdrId_chan     = -1;
  vars->hdrId_time_year= -1;
  vars->hdrId_time_day = -1;
  vars->hdrId_time_hour= -1;
  vars->hdrId_time_min = -1;
  vars->hdrId_time_sec = -1;
  vars->hdrId_time_nano= -1;
  vars->hdrId_source   = -1;
  vars->hdrId_sline    = -1;
  vars->hdrId_rline    = -1;
  vars->hdrId_rcv      = -1;
  vars->hdrId_serial   = -1;
  vars->hdrId_seq      = -1;
  vars->hdrId_sou_z    = -1;
  vars->hdrId_rec_wdep = -1;
  vars->hdrId_sou_elev = -1;
  vars->hdrId_rec_elev = -1;
  vars->hdrId_rec_x    = -1;
  vars->hdrId_rec_y    = -1;
  vars->hdrId_sou_x    = -1;
  vars->hdrId_sou_y    = -1;
  vars->hdrId_trc_type = -1;
  vars->hdrId_fileno   = -1;
  vars->hdrId_chan_set = -1;
  vars->hdrId_sampint_us = -1;
  vars->hdrId_time_samp1 = -1;
  vars->hdrId_time_samp1_us = -1;
  vars->hdrId_nsamp = -1;
  vars->chanSetIndexToRead = -1;

  vars->hdrId_extra  = NULL;
  vars->numExtraHdrs = 0;
  vars->filenameDumpHeaders = "";
  vars->isFirstCall = true;

  vars->hasCharHdr = false;
  vars->numCharBytes = 0;

  //---------------------------------------------------------------
  csSegdReader::configuration config;
  config.readAuxTraces     = false;
  config.isDebug           = edef->isDebug();
  config.navInterfaceID    = cseis_segd::UNKNOWN;
  config.numSamplesAddOne  = false;
  config.thisIsRev0        = false;
  config.navSystemID       = cseis_segd::NAV_HEADER_NONE;
  config.numBytes_genhdr2_extHdrBlocks = 2;

  if( edef->isDebug() ) fprintf(stdout,"Starting init phase of INPUT SEGD...\n");

  string hdrSetName;
  string text;
  string yesno;

  //********************************************************************************

  if( param->exists("chan_set") ) {
    param->getInt( "chan_set", &vars->chanSetIndexToRead );
    vars->chanSetIndexToRead -= 1; // Convert to C-style index
  }

  bool searchSubDirs = false;
  int numSingleFiles = param->getNumLines( "filename" );
  csVector<string> fileList;
  if( param->exists("directory") ) {
    int numValues = param->getNumValues("directory");
    string directory;
    param->getString("directory", &directory, 0);
    string extension = "segd";
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
    int length = s1.size();
    vars->filenames[i] = new char[length+1];
    memcpy(vars->filenames[i],s1.c_str(),length);
    vars->filenames[i][length] = '\0';
  }

  for( int i = 0; i < numSingleFiles; i++ ) {
    string filename;
    param->getStringAtLine( "filename", &filename, i );
    int length = filename.length();
    vars->filenames[i] = new char[length+1];
    memcpy(vars->filenames[i],filename.c_str(),length);
    vars->filenames[i][length] = '\0';
  }

  for( int i = 0; i < vars->numFiles; i++ ) {
    writer->line("Input file #%4d: %s", i, vars->filenames[i]);
  }
  writer->line("");

  //********************************************************************************

  if( param->exists("ntraces") ) {
    param->getInt( "ntraces", &vars->numTracesToRead );
    if( vars->numTracesToRead <= 0 ) vars->numTracesToRead = -1;  // Do not bother how many traces, read in all
  }
  if( param->exists("nrecords") ) {
    param->getInt( "nrecords", &vars->numRecordsToRead );
    if( vars->numRecordsToRead <= 0 ) vars->numRecordsToRead = -1;  // Do not bother how many records, read in all
  }

  if( param->exists("hdr_map") ) {
    string hdrSetName;
    param->getString( "hdr_map", &hdrSetName );
    if( !hdrSetName.compare("obc") ) {
      vars->hdrSet = csStandardSegdHeader::HDR_MAPPING_OBC;
    }
    else if( !hdrSetName.compare("standard") ) {
      vars->hdrSet = csStandardSegdHeader::HDR_MAPPING_STANDARD;
    }
    else {
      writer->error("Header set name not recognised: %s", hdrSetName.c_str());
    }
  }

  if( param->exists("nsamples_plus_one") ) {
    string text;
    param->getString( "nsamples_plus_one", &text );
    if( !text.compare("yes") ) {
      config.numSamplesAddOne = true;
    }
    else if( !text.compare("no") ) {
      config.numSamplesAddOne = false;
    }
    else {
      writer->error("Unknow option: %s", text.c_str());
    }
  }

  if( param->exists("nbytes_genhdr2_exthdr") ) {
    param->getInt( "nbytes_genhdr2_exthdr", &config.numBytes_genhdr2_extHdrBlocks );
  }

  if( param->exists("override_trchdr_ext") ) {
    string text;
    param->getString( "override_trchdr_ext", &text );
    if( !text.compare("yes") ) {
      config.overrideTraceHdrExtensions = true;
      param->getInt( "override_trchdr_ext", &config.numTraceHdrExtensions, 1 );
    }
    else if( !text.compare("no") ) {
      config.overrideTraceHdrExtensions = false;
    }
    else {
      writer->error("Unknow option: %s", text.c_str());
    }
  }

//----------------------------------------------------
  if( param->exists("rec_system") ) {
    param->getString( "rec_system", &text );
    text = toLowerCase(text);
    if( !text.compare("seal") ) {
      vars->recSystemID = cseis_segd::RECORDING_SYSTEM_SEAL;
    }
    else if( !text.compare("geores") ) {
      vars->recSystemID = cseis_segd::RECORDING_SYSTEM_GEORES;
    }
    else if( !text.compare("fairfield") ) {
      vars->recSystemID = cseis_segd::RECORDING_SYSTEM_FAIRFIELD;
    }
   else if( !text.compare("unknown") ) {
      vars->recSystemID = cseis_segd::UNKNOWN;
    }
    else {
      writer->line("Option not recognized: '%s'", text.c_str() );
      env->addError();
    }
  }
  config.recordingSystemID = vars->recSystemID;

//----------------------------------------------------
  if( param->exists("nav_header") ) {
    param->getString( "nav_header", &text );
    text = toLowerCase(text);
    if( !text.compare("hydronav") ) {
      config.navSystemID = cseis_segd::NAV_HEADER_HYDRONAV_VER6;
    }
    else if( !text.compare("labo") ) {
      config.navSystemID = cseis_segd::NAV_HEADER_LABO;
    }
    else if( !text.compare("psi") ) {
      config.navSystemID = cseis_segd::NAV_HEADER_PSI;
    }
    else if( !text.compare("none") ) {
      config.navSystemID = cseis_segd::NAV_HEADER_NONE;
    }
    else {
      writer->line("Option not recognized: '%s'", text.c_str() );
      env->addError();
    }
  }

  //----------------------------------------------------
  if( param->exists("revision_0") ) {
    param->getString( "revision_0", &text );
    text = toLowerCase(text);
    if( !text.compare("yes") ) {
      config.thisIsRev0 = true;
    }
    else if( !text.compare("no") ) {
      config.thisIsRev0 = false;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  if( param->exists("read_aux") ) {
    param->getString( "read_aux", &text );
    text = toLowerCase(text);
    if( !text.compare("yes") ) {
      config.readAuxTraces = true;
    }
    else if( !text.compare("no") ) {
      config.readAuxTraces = false;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  /// Nav interface IDs
//  static int const CM_DIGI_COMP_A = 101;  // ..used on Venturer
//  static int const CM_DIGI_PSIB   = 102;  // ..used on Search
//  static int const CM_DIGI_TS     = 103;  //    -"-
//  static int const CM_DIGI_BIRD_DEPTHS = 104;


  //------------------------------------------------
  // Dump headers?
  //
  if( param->exists("dump_essential") ) {
    param->getString( "dump_essential", &yesno );
    yesno = toLowerCase( yesno );
    if( !yesno.compare( "yes" ) ) {
      vars->dumpEssFileInfo = true;
    }
    else if( !yesno.compare( "no" ) ) {
      vars->dumpEssFileInfo = false;
    }
    else {
      writer->line("Parameter value not recognized: %s", yesno.c_str() );
      env->addError();
    }
  }

  if( param->exists("dump_hdr") ) {
    csVector<std::string> valueList;
    param->getAll( "dump_hdr", &valueList );
    int nFields = valueList.size();
    for( int i = 0; i < nFields; i++ ) {
      std::string text = toLowerCase( valueList.at(i) );
      if( !text.compare("all") ) {
        vars->dumpHdrFlag |= cseis_segd::DUMP_ALL;
      }
      else if( !text.compare("general") ) {
        vars->dumpHdrFlag |= cseis_segd::DUMP_GENERAL;
      }
      else if( !text.compare("chanset") ) {
        vars->dumpHdrFlag |= cseis_segd::DUMP_CHANSET;
      }
      else if( !text.compare("extended") ) {
        vars->dumpHdrFlag |= cseis_segd::DUMP_EXTENDED;
      }
      else if( !text.compare("external") ) {
        vars->dumpHdrFlag |= cseis_segd::DUMP_EXTERNAL;
      }
      else if( !text.compare("birds") ) {
        vars->dumpHdrFlag |= cseis_segd::DUMP_EXTERNAL_BIRDS;
      }
      else if( !text.compare("raw") ) {
        vars->dumpHdrFlag |= cseis_segd::DUMP_RAW;
      }
      else if( !text.compare("none") ) {
        vars->dumpHdrFlag = cseis_segd::DUMP_NONE;
        break;
      }
      else if( !text.compare("trace") ) {
        vars->dumpTrcHeaders = true;      
      }
      else {
        writer->warning(" Option for user parameter 'dump_hdr' not recognised: '%s'", text.c_str());
      }
    }
  }

  if( param->exists("dump_filename") ) {
    param->getString( "dump_filename", &vars->filenameDumpHeaders );
  }
//  else if( vars->dumpHdrFlag != 0 ) {
//    vars->file_dumpHeaders = (std::ofstream*)(&std::cout);
//    printf("DUMP headers to standard output!\n");
//  }

  //-----------------------------------------------------
  // 
  //

  csVector<csHeaderInfo const*> stdHdrList;
  csStandardSegdHeader::setStandardHeaders( vars->hdrSet, &stdHdrList );

  int nHeaders = stdHdrList.size();
  vars->hdrIndexSegd = new int[nHeaders];
  vars->hdrTypeSegd  = new type_t[nHeaders];
  for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
    csHeaderInfo const* info = stdHdrList.at(ihdr);
    if( !hdef->headerExists(info->name) ) {
      vars->hdrTypeSegd[ihdr]  = info->type;
      vars->hdrIndexSegd[ihdr] = hdef->addHeader( info );
    }
    else {
      vars->hdrIndexSegd[ihdr] = hdef->headerIndex( info->name );
      vars->hdrTypeSegd[ihdr]  = hdef->headerType( info->name );
    }
    if( edef->isDebug() ) writer->line("Create new header %s", info->name.c_str());
  }
  stdHdrList.clear();

  //----------------------------------------------------
  std::string headerName;
  if( param->exists("header_ens") ) {
    param->getString( "header_ens", &headerName );
//    headerName = toLowerCase(headerName);
    writer->line("Specified ensemble header: '%s'.", headerName.c_str() );
    if( !hdef->headerExists( headerName ) ) {
      writer->error("Specified ensemble header does not exist: '%s'.", headerName.c_str() );
    }
    type_t type_hdr_ens = hdef->headerType( headerName );
    if( type_hdr_ens != TYPE_INT && type_hdr_ens != TYPE_FLOAT && type_hdr_ens != TYPE_DOUBLE ) {
      writer->error("Ensemble header can only be of number type, not string or array type.");
    }
    shdr->setEnsembleKey( headerName, 0 );
  }
  else {
    shdr->setEnsembleKey( "ffid", 0 );
  }

  if( !hdef->headerExists(HDR_FILENO.name) ) {
    hdef->addStandardHeader( HDR_FILENO.name );
  }
  vars->hdrId_fileno = hdef->headerIndex(HDR_FILENO.name);

  if( vars->dumpHdrFlag != 0 || vars->dumpTrcHeaders ) {
    if( vars->filenameDumpHeaders.size() > 0 ) {
      if( !csFileUtils::createDoNotOverwrite( vars->filenameDumpHeaders ) ) {
        writer->error("Error when opening header dump file '%s' for writing.", vars->filenameDumpHeaders.c_str() );
      }
    }
  }
  
  if( param->exists("char_hdr") ) {
    std::string text;
    param->getString("char_hdr",&text,0);
    if( !text.compare("yes") ) {
      vars->hasCharHdr = true;      
      param->getInt("char_hdr",&vars->numCharBytes,1);
    }
    else if( !text.compare("no") ) {
      vars->hasCharHdr = false;      
    }
    else {
      writer->warning(" Option for user parameter 'dump_hdr' not recognised: '%s'", text.c_str());
    }
  }
 
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
  
  //----------------------------------------------------
  // Open SEGD input file
  // Read in headers of first record
  //

  vars->segdReader = new csSegdReader();
  vars->segdReader->setConfiguration( config );

  try {
    vars->segdReader->open( vars->filenames[0] );
    if( vars->hasCharHdr ) {
      char* buffer_charHdr = new char[vars->numCharBytes+1];
      vars->segdReader->readCharHdr( buffer_charHdr, vars->numCharBytes );
      buffer_charHdr[vars->numCharBytes] = '\0';
      writer->write( "%d-byte ASCII header:\n%s", vars->numCharBytes, buffer_charHdr );
      delete [] buffer_charHdr;
    }
    vars->segdReader->readNewRecordHeaders();

    // TEMP
    //    if( vars->dumpHdrFlag != 0 ) {
    //  vars->segdReader->dumpFileHeaders( vars->dumpHdrFlag, vars->file_dumpHeaders );
    // }
    // exit(-1);
    // TEMP

    if( vars->chanSetIndexToRead >= 0 ) vars->segdReader->setChanSetToRead( vars->chanSetIndexToRead );
    bool isSuccess = vars->segdReader->readNextRecord( vars->comRecordHdr );
    if( !isSuccess ) writer->error("Cannot read in first record");
  }
  catch( csException& e ) {
    try {
      vars->segdReader->dumpEssentialFileInfo( writer->getFile() );
      vars->dumpHdrFlag = cseis_segd::DUMP_GENERAL;
      vars->dumpHdrFlag |= cseis_segd::DUMP_CHANSET;
      vars->dumpHdrFlag |= cseis_segd::DUMP_EXTENDED;
      vars->segdReader->dumpFileHeaders( vars->dumpHdrFlag, (std::ofstream*)(&std::cout) );
    } catch(...) {}
    writer->error("Error when opening SEGD file '%s'.\nSystem message: %s", vars->filenames[0], e.getMessage() );
  }
  catch( string text ) {
    writer->error("Error when opening SEGD file '%s'.\nSystem message: %s", vars->filenames[0], text.c_str() );
  }
  catch( ... ) {
    writer->error("Unknown error occurred when opening SEGD file '%s'.", vars->filenames[0] );
  }

  //------------------------------------
  // Set standard trace headers

  if( !hdef->headerExists("ffid") ) hdef->addStandardHeader("ffid");
  if( !hdef->headerExists("trcno") ) hdef->addStandardHeader("trcno");
  if( !hdef->headerExists("chan") ) hdef->addStandardHeader("chan");
  if( !hdef->headerExists("time_year") ) hdef->addStandardHeader("time_year");
  if( !hdef->headerExists("time_day") ) hdef->addStandardHeader("time_day");
  if( !hdef->headerExists("time_hour") ) hdef->addStandardHeader("time_hour");
  if( !hdef->headerExists("time_min") ) hdef->addStandardHeader("time_min");
  if( !hdef->headerExists("time_sec") ) hdef->addStandardHeader("time_sec");
  if( !hdef->headerExists("time_nano") ) hdef->addStandardHeader("time_nano");
  if( !hdef->headerExists("source") ) hdef->addStandardHeader("source");
  if( !hdef->headerExists("sou_line") ) hdef->addStandardHeader("sou_line");
  if( !hdef->headerExists("rec_line") ) hdef->addStandardHeader("rec_line");
  if( !hdef->headerExists("rcv") ) hdef->addStandardHeader("rcv");
  if( !hdef->headerExists("serial") ) hdef->addStandardHeader("serial");
  if( !hdef->headerExists("seq") ) hdef->addStandardHeader("seq");
  if( !hdef->headerExists("sou_z") ) hdef->addStandardHeader("sou_z");
  if( !hdef->headerExists("rec_wdep") ) hdef->addStandardHeader("rec_wdep");
  if( !hdef->headerExists("sou_elev") ) hdef->addStandardHeader("sou_elev");
  if( !hdef->headerExists("rec_elev") ) hdef->addStandardHeader("rec_elev");
  if( !hdef->headerExists("rec_x") ) hdef->addStandardHeader("rec_x");
  if( !hdef->headerExists("rec_y") ) hdef->addStandardHeader("rec_y");
  if( !hdef->headerExists("sou_x") ) hdef->addStandardHeader("sou_x");
  if( !hdef->headerExists("sou_y") ) hdef->addStandardHeader("sou_y");
  if( !hdef->headerExists("trc_type") ) hdef->addStandardHeader("trc_type");
  if( !hdef->headerExists("chan_set") ) hdef->addHeader(TYPE_INT,"chan_set","SEG-D channel set number");
  if( !hdef->headerExists("nsamp") ) hdef->addStandardHeader("nsamp");
  if( !hdef->headerExists("sampint_us") ) hdef->addStandardHeader("sampint_us");

  vars->hdrId_ffid      = hdef->headerIndex("ffid");
  vars->hdrId_trcno     = hdef->headerIndex("trcno");
  vars->hdrId_chan      = hdef->headerIndex("chan");
  vars->hdrId_time_year = hdef->headerIndex("time_year");
  vars->hdrId_time_day  = hdef->headerIndex("time_day");
  vars->hdrId_time_hour = hdef->headerIndex("time_hour");
  vars->hdrId_time_min  = hdef->headerIndex("time_min");
  vars->hdrId_time_sec  = hdef->headerIndex("time_sec");
  vars->hdrId_time_nano = hdef->headerIndex("time_nano");

  vars->hdrId_source    = hdef->headerIndex("source");
  vars->hdrId_sline    = hdef->headerIndex("sou_line");
  vars->hdrId_rline    = hdef->headerIndex("rec_line");
  vars->hdrId_rcv      = hdef->headerIndex("rcv");
  vars->hdrId_serial   = hdef->headerIndex("serial");
  vars->hdrId_seq      = hdef->headerIndex("seq");

  vars->hdrId_sou_z    = hdef->headerIndex("sou_z");
  vars->hdrId_rec_wdep = hdef->headerIndex("rec_wdep");

  vars->hdrId_sou_elev = hdef->headerIndex("sou_elev");
  vars->hdrId_rec_elev = hdef->headerIndex("rec_elev");
  vars->hdrId_rec_x    = hdef->headerIndex("rec_x");
  vars->hdrId_rec_y    = hdef->headerIndex("rec_y");
  vars->hdrId_sou_x    = hdef->headerIndex("sou_x");
  vars->hdrId_sou_y    = hdef->headerIndex("sou_y");
  vars->hdrId_trc_type = hdef->headerIndex("trc_type");
  vars->hdrId_chan_set = hdef->headerIndex("chan_set");
  vars->hdrId_sampint_us = hdef->headerIndex("sampint_us");
  vars->hdrId_nsamp = hdef->headerIndex("nsamp");
  vars->hdrId_time_samp1    = hdef->headerIndex( HDR_TIME_SAMP1.name );
  vars->hdrId_time_samp1_us = hdef->headerIndex( HDR_TIME_SAMP1_US.name );

  //-------------------------------------------------
  // Set non-standard trace headers
  if( vars->recSystemID == cseis_segd::RECORDING_SYSTEM_GEORES && vars->hdrSet == csStandardSegdHeader::HDR_MAPPING_OBC ) {
    vars->numExtraHdrs = 4;
    vars->hdrId_extra = new int[vars->numExtraHdrs];
    vars->hdrId_extra[0] = hdef->headerIndex("incl_i");
    vars->hdrId_extra[1] = hdef->headerIndex("incl_c");
    vars->hdrId_extra[2] = hdef->headerIndex("incl_v");
    vars->hdrId_extra[3] = hdef->headerIndex("sensor");
  }
  //-------------------------------------------------

  int traceByteSize = vars->segdReader->traceDataByteSize();
  if( edef->isDebug() ) fprintf(stderr,"Trace byte size: %d\n", traceByteSize);

  vars->comFileHdr = vars->segdReader->getCommonFileHeaders();

  if( vars->dumpEssFileInfo ) {
    vars->segdReader->dumpEssentialFileInfo( writer->getFile() );
  }

  commonChanSetStruct info;
  if( vars->chanSetIndexToRead < 0 ) {
    shdr->sampleInt = (float)vars->comFileHdr->sampleInt_us / 1000.0;
    shdr->numSamples = vars->comFileHdr->numSamples;
    for( int i = 0; i < vars->segdReader->numChanSets(); i++ ) {
      vars->segdReader->retrieveChanSetInfo( i, info );
      if( info.numChannels == 0 ) continue;
      if( info.numSamples > shdr->numSamples ) {
      vars->segdReader->dumpEssentialFileInfo( writer->getFile() );
      vars->dumpHdrFlag = cseis_segd::DUMP_GENERAL;
      vars->dumpHdrFlag |= cseis_segd::DUMP_CHANSET;
      vars->dumpHdrFlag |= cseis_segd::DUMP_EXTENDED;
      vars->segdReader->dumpFileHeaders( vars->dumpHdrFlag, (std::ofstream*)(&std::cout) );

        writer->error("Number of samples for channel set #%d (=%d) exceeds number of samples specified in SEG-D general header (=%d) or trace header extension (=%d).\nThis is not supported yet. Please read in a single channel set at a time.",
                   i+1, info.numSamples, (int)round(vars->segdReader->generalHdr1()->recordLenMS/shdr->sampleInt), shdr->numSamples);
      }
      else if( info.numSamples != shdr->numSamples ) {
        writer->warning("Number of samples for channel set #%d (=%d) differs from number of samples specified in SEG-D general header or trace header extension (=%d).\n  Suggestion: Dump all SEG-D header information to the terminal by setting parameter  'dump_hdr all'\n   ...otherwise try to set parameter  'nsamples_plus_one  yes/no'  and see if this helps (often it does)",
                     i+1, info.numSamples, shdr->numSamples);
      }
      if( vars->comFileHdr->sampleInt_us != info.sampleInt_us ) {
        writer->warning("Sample interval for channel set #%d (=%.4fms) differs from sample interval specified in SEG-D general header (=%.4fms).",
                     i+1, (float)(info.sampleInt_us/1000.0), shdr->sampleInt);
      }
    }
  }
  else {  // Set number of samples & sample interval to the values for this specific channel set
    if( vars->chanSetIndexToRead >= vars->segdReader->numChanSets() ) {
      writer->error("Cannot read in user specified channel set #%d: Input file only contains %d channel sets per scan type.",
                 vars->chanSetIndexToRead+1,
                 vars->segdReader->numChanSets() );
    }
    vars->segdReader->retrieveChanSetInfo( vars->chanSetIndexToRead, info );   // -1 to convert to C-style index
    shdr->sampleInt  = info.sampleInt_us / 1000.0;
    shdr->numSamples = info.numSamples;
  }

  writer->line("");
  writer->line("  File name:             %s", vars->filenames[0]);
  writer->line("  Sample interval [ms]:  %f", shdr->sampleInt);
  writer->line("  Number of samples:     %d", shdr->numSamples);
  writer->line("");

  vars->totalTraceCounter = 0;
  vars->traceCounter      = 0;
  vars->recordCounter     = 0;
  vars->isAtEOF           = false;
}

//==================================================================================
//
//
//==================================================================================

void exec_mod_input_segd_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  //  csSuperHeader const* shdr = env->superHeader;
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;

  csTrace* trace = traceGather->trace(0);



  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    if( vars->dumpHdrFlag != 0 || vars->dumpTrcHeaders ) {
      try {
        if( vars->filenameDumpHeaders.size() > 0 ) {
          vars->file_dumpHeaders = new std::ofstream ( vars->filenameDumpHeaders.c_str() );
          if( vars->file_dumpHeaders == NULL ) {
            writer->error("Could not open header dump file '%s'.", vars->filenameDumpHeaders.c_str() );
          }
        }
        else {
          vars->file_dumpHeaders = (std::ofstream*)(&std::cout);
        }
      }
      catch( string text ) {
        writer->error("Error when opening header dump file '%s'.\nSystem message: %s", vars->filenameDumpHeaders.c_str(), text.c_str() );
      }
      catch(...) {
        writer->error("Error when opening header dump file '%s'.", vars->filenameDumpHeaders.c_str() );
      }
      vars->segdReader->dumpFileHeaders( vars->dumpHdrFlag, vars->file_dumpHeaders );
    }
  }

  if( vars->isAtEOF ) {
    traceGather->freeAllTraces();
    return;
  }
//---------------------------------------------------------------
  cseis_segd::commonTraceHeaderStruct  comTrcHdr;

  float* samples = trace->getTraceSamples();
  
  // Retrieve next trace. If none is available for current record, try to read in the next record
  // (or the remaining part of the current record if it is the first one)
  //  float* buffer = new float[4*trace->numSamples()];


  while( !vars->segdReader->getNextTrace( samples, comTrcHdr ) ) {
    try {
      if( edef->isDebug() ) writer->write("Read next record... #%d", (vars->recordCounter+1) );
      cseis_segd::commonRecordHeaderStruct newComRecordHdr;
      bool isSuccess = vars->segdReader->readNextRecord( newComRecordHdr );
      while( !isSuccess && vars->currentFile < vars->numFiles-1 ) {
        vars->currentFile += 1;
        vars->segdReader->open( vars->filenames[vars->currentFile] );
        vars->segdReader->readNewRecordHeaders();
        isSuccess = true;
        if( vars->chanSetIndexToRead >= 0 ) vars->segdReader->setChanSetToRead( vars->chanSetIndexToRead );
        isSuccess = vars->segdReader->readNextRecord( newComRecordHdr );
        if( !isSuccess ) writer->warning("INPUT_SEGD: Cannot read in first record form file %s", vars->filenames[vars->currentFile]);
      }
      if( isSuccess ) {
        vars->comRecordHdr = newComRecordHdr;
        vars->recordCounter += 1;
        vars->traceCounter = 0;
        if( edef->isDebug() ) writer->line(", ...read file number (total traces so far): %d (%ld)", vars->comRecordHdr.fileNum, vars->totalTraceCounter );
        if( vars->dumpHdrFlag != 0 && vars->recordCounter != 0 ) {
          vars->segdReader->dumpFileHeaders( vars->dumpHdrFlag, vars->file_dumpHeaders );
        }
        if( vars->dumpEssFileInfo ) {
          vars->segdReader->dumpEssentialFileInfo( writer->getFile() );
        }
        if( vars->dumpHdrFlag != cseis_segd::DUMP_NONE ) {
          vars->segdReader->dumpRecordHeaders( writer->getFile() );
          vars->comRecordHdr.dump( writer->getFile() );
        }
      }
      // ...last record has been reached. Terminate with grace
      else {
        writer->line("");
        vars->isAtEOF = true;
        traceGather->freeAllTraces();
        return;
      }
    }
    catch( cseis_geolib::csException& exc ) {
      writer->error("\nException occurred: %s\nTerminated...\n", exc.getMessage());
    }
    catch( ... ) {
      writer->error("\nUnknown exception occurred during SEGD read operation...\nTerminated...\n");
    }
  }

  // Set any unassigned samples at end of trace to zero
  for( int isamp = comTrcHdr.numSamples; isamp < shdr->numSamples; isamp++ ) {
    samples[isamp] = 0.0;
  }

  // If specified number of traces or records have been read in --> terminate gracefully
  if( (vars->numTracesToRead > 0 && vars->totalTraceCounter > vars->numTracesToRead) ||
      ((vars->numRecordsToRead > 0 && vars->recordCounter > vars->numRecordsToRead) ) ) {
     vars->isAtEOF = true;
    traceGather->freeAllTraces();
    return;
  }
  
  vars->traceCounter += 1;
  vars->totalTraceCounter += 1;
  if( vars->dumpTrcHeaders ) {
    writer->line("Trace #%d (total #%d), traces to read: %d", vars->traceCounter, vars->totalTraceCounter, vars->numTracesToRead );
    // comTrcHdr.dump( writer->getFile() );
    vars->segdReader->dumpTrace( vars->file_dumpHeaders );
  }

  //-------------------------------------------------------------
  if( edef->isDebug() ) {
    writer->line( "Reading in trace #%d of %d", vars->traceCounter, vars->comFileHdr->totalNumChan );
    writer->flush();
  }
// Assume TOWED:
  cseis_system::csTraceHeader* trcHdr = trace->getTraceHeader();

  trcHdr->setIntValue(vars->hdrId_ffid, vars->comRecordHdr.fileNum );
  trcHdr->setIntValue(vars->hdrId_chan, comTrcHdr.chanNum );
  trcHdr->setIntValue(vars->hdrId_trcno, vars->traceCounter );
  int year = vars->comRecordHdr.shotTime.year;
  if( year < 1900 ) {
    if( year > 30 ) {
      year += 1900;
    }
    else {
      year += 2000;
    }
  }
  trcHdr->setIntValue(vars->hdrId_time_year, year );
  trcHdr->setIntValue(vars->hdrId_time_day, vars->comRecordHdr.shotTime.julianDay );
  trcHdr->setIntValue(vars->hdrId_time_hour, vars->comRecordHdr.shotTime.hour );
  trcHdr->setIntValue(vars->hdrId_time_min, vars->comRecordHdr.shotTime.minute );
  trcHdr->setIntValue(vars->hdrId_time_sec, vars->comRecordHdr.shotTime.second );

  int nanosec = vars->segdReader->traceHdrExtension()->nanoSeconds();
  //  fprintf(stderr,"Shottime nanosec:  %d\n", nanosec );
  if( nanosec == 0 ) {
    nanosec = vars->segdReader->extendedHdr()->nanoSeconds();
  }
  csInt64_t time_us = vars->segdReader->traceHdrExtension()->timeSamp1_us();
  int time_samp1_s  = 0;
  int time_samp1_us = 0;
  if( time_us != 0 ) {
    time_samp1_s  = (int)( time_us / 1000000ULL );
    time_samp1_us = (int)( time_us - (csInt64_t)time_samp1_s * 1000000ULL );
  }
  else {
    time_samp1_us = nanosec/1000;
    int extraSeconds  = time_samp1_us/1000/1000;
    time_samp1_us    -= extraSeconds*1000*1000;
    time_samp1_s  = csGeolibUtils::date2UNIXsec( year,
                                                 vars->comRecordHdr.shotTime.julianDay,
                                                 vars->comRecordHdr.shotTime.hour,
                                                 vars->comRecordHdr.shotTime.minute,
                                                 vars->comRecordHdr.shotTime.second + extraSeconds );
  }
  trcHdr->setIntValue(vars->hdrId_time_nano, nanosec );
  trcHdr->setIntValue( vars->hdrId_time_samp1, time_samp1_s );
  trcHdr->setIntValue( vars->hdrId_time_samp1_us, time_samp1_us );
  if( edef->isDebug() ) {
    writer->line("Start time: %lld.%06ds  =  %ds  +  %dus", time_samp1_s, nanosec/1000, time_samp1_s, time_samp1_us );
  }

  trcHdr->setIntValue(vars->hdrId_chan_set, comTrcHdr.chanSet );
  trcHdr->setIntValue(vars->hdrId_nsamp, comTrcHdr.numSamples );
  trcHdr->setIntValue(vars->hdrId_sampint_us, comTrcHdr.sampleInt_us );

  trcHdr->setIntValue(vars->hdrId_source, vars->comRecordHdr.shotNum );
  trcHdr->setIntValue(vars->hdrId_sline, 0 );
  trcHdr->setIntValue(vars->hdrId_rline, comTrcHdr.rcvLineNumber );
  trcHdr->setIntValue(vars->hdrId_rcv, comTrcHdr.rcvPointNumber );
  trcHdr->setIntValue(vars->hdrId_serial, comTrcHdr.serialNumber );
  trcHdr->setIntValue(vars->hdrId_seq, 0 );
//  trcHdr->setIntValue(vars->hdrId_traceEdit, comTrcHdr.traceEdit );

  trcHdr->setFloatValue(vars->hdrId_sou_z, 0 );
  trcHdr->setFloatValue(vars->hdrId_rec_wdep, 0 );

  trcHdr->setFloatValue(vars->hdrId_sou_elev, vars->comRecordHdr.srcElev );
  trcHdr->setFloatValue(vars->hdrId_rec_elev, comTrcHdr.rcvElevation );
  trcHdr->setDoubleValue(vars->hdrId_rec_x, comTrcHdr.rcvEasting );
  trcHdr->setDoubleValue(vars->hdrId_rec_y, comTrcHdr.rcvNorthing );
  trcHdr->setDoubleValue(vars->hdrId_sou_x, vars->comRecordHdr.srcEasting );
  trcHdr->setDoubleValue(vars->hdrId_sou_y, vars->comRecordHdr.srcNorthing );

  trcHdr->setIntValue(vars->hdrId_trc_type, comTrcHdr.chanTypeID  );

  if( vars->recSystemID == cseis_segd::RECORDING_SYSTEM_GEORES && vars->hdrSet == csStandardSegdHeader::HDR_MAPPING_OBC ) {
    trcHdr->setIntValue(vars->hdrId_extra[0], comTrcHdr.incl_i );
    trcHdr->setIntValue(vars->hdrId_extra[1], comTrcHdr.incl_c );
    trcHdr->setIntValue(vars->hdrId_extra[2], comTrcHdr.incl_v );
    trcHdr->setIntValue(vars->hdrId_extra[3], comTrcHdr.sensor );
  }
  trcHdr->setIntValue( vars->hdrId_fileno, vars->currentFile+1 );


  return;
  
}

//********************************************************************************
// Parameter definition
//
//
//********************************************************************************

void params_mod_input_segd_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_SEGD", "Input SEGD data", "Reads SEGD data from disk file" );

  pdef->addDoc("This module does not have production-type quality. The SEG-D format has been implemented fairly thoroughly, but with a few omissions regarding certain extensions, trailers and others.");
  pdef->addDoc("However, the complexity of the format, and its sometimes unclear definition leads to many variations in the implementation of the format by various vendors.");
  pdef->addDoc("One obvious problem is the lack of a clear definition for the number of samples in a trace. Hence, some vendors generate SEG-D files with N samples, while others create SEG-D data with N+1 samples (where N is the record length divided by the sample interval). Use the user parameter 'nsamples_plus_one' to add 1 to the number of samples in the general header when needed.");
  pdef->addDoc("<br>");
  pdef->addDoc("By default, INPUT_SEGD only reads in seismic data channels, not auxiliary channels. To read in auxiliary channels, specify user parameter 'read_aux yes'.");

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name including full path name" );

  pdef->addParam( "directory", "Name of directory to search", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Name of directory", "All files in this directory and its subdirectories will be read" );
  pdef->addValue( "segd", VALTYPE_STRING, "File extension", "Only files with the given file extension will be read in" );
  pdef->addValue( "no", VALTYPE_OPTION, "Search subdirectories" );
  pdef->addOption( "no", "Do not search subdirectories" );
  pdef->addOption( "yes", "Also search subdirectories for files" );

  pdef->addParam( "dump_filename", "Dump file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Dump file name including full path name" );

  pdef->addParam( "chan_set", "Channel set number to read in", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Specify channel set number to read in" );

  pdef->addParam( "nsamples_plus_one", "Add 1 to number of samples computed from general header record length?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Add 1 to number of samples from general header record length", "number of samples = record length / sample interval + 1" );
  pdef->addOption( "no", "Do not add 1 to number of samples" );
  //  pdef->addOption( "auto", "Attempt to automatically determine number of samples" );

  pdef->addParam( "ntraces", "Number of traces to read in", NUM_VALUES_FIXED,
                  "Input of traces will stop when all traces have been read in, or when the number of traces specified has been reached. Traces will not be filled up to the specified range");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to read in (0: read in all)" );

  pdef->addParam( "nrecords", "Number of records to read in", NUM_VALUES_FIXED,
                  "Input of data will stop when the number of records specified has been reached.");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of records to read in (0: read in all)" );

  pdef->addParam( "dump_hdr", "Dump header", NUM_VALUES_FIXED );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "none", "Do not dump SEGD headers" );
  pdef->addOption( "all", "Complete dump of all SEGD headers" );
  pdef->addOption( "general", "Dump SEGD general headers" );
  pdef->addOption( "chanset", "Dump SEGD chan set headers" );
  pdef->addOption( "extended", "Dump SEGD extended header" );
  pdef->addOption( "external", "Dump SEGD external header" );
  pdef->addOption( "external_bird", "Dump SEGD external bird header" );
  pdef->addOption( "trace", "Dump trace headers" );

  pdef->addParam( "header_ens", "Ensemble trace header name", NUM_VALUES_FIXED,
                  "Ensembles are defined by this key header" );
  pdef->addValue( "", VALTYPE_STRING, "Ensemble trace header name" );

  pdef->addParam( "dump_essential", "Dump essential headers", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Dump essential headers" );
  pdef->addOption( "no", "Do not dump essential headers" );

  pdef->addParam( "hdr_map", "Specify header set", NUM_VALUES_FIXED );
  pdef->addValue( "standard", VALTYPE_OPTION );
  pdef->addOption( "standard", "Standard header set" );
  pdef->addOption( "obc", "OBC 4C header set" );

  pdef->addParam( "rec_system", "Recording system ID name/product name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "seal", "Sercel SEAL recording system" );
  pdef->addOption( "geores", "Geospace GEORES recording system" );
  pdef->addOption( "fairfield", "Fairfield Nodal recording system" );

  pdef->addParam( "nav_header", "Navigation header type", NUM_VALUES_FIXED, "The navigation header forms the first part of the external header. If data does not read in correctly, try different nav headers, or 'none'" );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "hydronav", "Hydronav header" );
  pdef->addOption( "psi", "PSI header" );
  pdef->addOption( "labo", "LABO header" );
  pdef->addOption( "none", "SEGD file does not contain any nav header in the external header" );

  pdef->addParam( "revision_0", "Is this a revision 0 file?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "This is a revision 0 file" );
  pdef->addOption( "no", "Read revision number from general header 2" );

  pdef->addParam( "read_aux", "Read in auxiliary channels?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Read in auxiliary channels" );
  pdef->addOption( "no", "Do not read in auxiliary channels" );

  pdef->addParam( "override_trchdr_ext", "Override trace header extensions in channels set/trace header?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Override number of trace header extensions" );
  pdef->addOption( "no", "Do not override number of trace header extensions" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of trace header extensions" );

  pdef->addParam( "nbytes_genhdr2_exthdr", "Number of bytes of general header #2 field to read number of external header blocks from", NUM_VALUES_FIXED );
  pdef->addValue( "2", VALTYPE_NUMBER );

  pdef->addParam( "char_hdr", "Does SEGD file have ASCII character header? If yes, specify number of bytes in char header", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "No char header" );
  pdef->addOption( "yes", "Has char header" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of bytes in char header" );

//  pdef->addParam( "decode_hdr", "Decode additional header", NUM_VALUES_FIXED );
//  pdef->addValue( "", VALTYPE_STRING, "Name of additional SEGD header to decode", "Must be one of: nano_sec" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_input_segd_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_input_segd::VariableStruct* vars = reinterpret_cast<mod_input_segd::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_input_segd_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_input_segd::VariableStruct* vars = reinterpret_cast<mod_input_segd::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->hdrId_extra ) {
    delete [] vars->hdrId_extra;
    vars->hdrId_extra = NULL;
  }
  if( vars->segdReader != NULL ) {
    delete vars->segdReader;
  }
  if( vars->hdrIndexSegd ) {
    delete [] vars->hdrIndexSegd; vars->hdrIndexSegd = NULL;
  }
  if( vars->hdrTypeSegd ) {
    delete [] vars->hdrTypeSegd; vars->hdrTypeSegd = NULL;
  }
  if( vars->file_dumpHeaders ) {
    vars->file_dumpHeaders->close();
    //      delete vars->file_dumpHeaders;  Do not delete this
    vars->file_dumpHeaders = NULL;
  }
  if( vars->filenames ) {
    for( int i = 0; i < vars->numFiles; i++ ) {
      delete [] vars->filenames[i];
    }
    delete [] vars->filenames; vars->filenames = NULL;
  }
  delete vars;
  vars = NULL;
}

extern "C" void _params_mod_input_segd_( csParamDef* pdef ) {
  params_mod_input_segd_( pdef );
}
extern "C" void _init_mod_input_segd_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_segd_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_segd_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_segd_( env, writer );
}
extern "C" void _exec_mod_input_segd_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_segd_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_segd_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_segd_( env, writer );
}

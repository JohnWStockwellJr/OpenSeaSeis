/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSegyWriter.h"
#include "csSegyDefines.h"
#include "csSegyHeaderInfo.h"
#include "csSegyTraceHeader.h"
#include "csSegyHdrMap.h"
#include "csSegyBinHeader.h"
#include <cstring>
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: OUTOUT_SEGY
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_output_segy {
  struct VariableStruct {
    long traceCounter;
    csSegyWriter* segyWriter;
    int*  hdrIndexSegy;
    type_t* hdrTypeSegy;
    int   nsampMAX;
    int   nTracesRead;
    int   hdr_mapping;
    int   numSamplesOut;
    bool isFirstCall;
  };
  static int const HDR_MANDATORY_DEFAULT = 1;
  static int const HDR_MANDATORY_HEADER  = 2;
}
using namespace mod_output_segy;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_output_segy_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->traceCounter = 0;
  vars->segyWriter   = NULL;
  vars->hdrIndexSegy = NULL;
  vars->hdrTypeSegy  = NULL;
  vars->nsampMAX     = 0;
  vars->nTracesRead  = 0;
  vars->hdr_mapping  = csSegyHdrMap::SEGY_STANDARD;
  vars->numSamplesOut = shdr->numSamples;
  vars->isFirstCall = true;

  //--------------------------------------------------
  std::string headerName;
  std::string yesno;
  std::string filename;
  int numTracesBuffer;

  bool doForce = false;
  if( param->exists("force") ) {
    param->getString( "force", &yesno );
    yesno = toLowerCase( yesno );
    if( !yesno.compare( "yes" ) ) {
      doForce = true;
    }
    else if( !yesno.compare( "no" ) ) {
      doForce = false;
    }
    else {
      writer->line("%s: Parameter value not recognized: %s", edef->moduleName().c_str(), yesno.c_str() );
      env->addError();
    }
  }

  //---------------------------------------------------------
  float sampleInt = shdr->sampleInt;
  if( param->exists("sample_int") ) {
    param->getFloat( "sample_int", &sampleInt );
    if( sampleInt <= 0 ) {
      writer->error("Wrong sample interval specified: %f", sampleInt);
    }
  }

  unsigned int maxUnsignedShort = (unsigned int)pow(2.0,16);
  unsigned int sampleIntUS = (unsigned int)( sampleInt * 1000.0 );
 
  if( sampleIntUS <= 0 ) {
    if( doForce ) {
      writer->warning("Data sample interval (=%fus) is smaller than minimum possible value that can be saved to SEGHY format (=1us)\n", sampleInt * 1000.0 );
    }
    else {
      writer->error("Data sample interval (=%fus) is smaller than minimum possible value that can be saved to SEGHY format (=1us)\n", sampleInt * 1000.0 );
    }
  }
  if( sampleIntUS >= maxUnsignedShort ) {
    if( doForce ) {
      writer->warning("Data sample interval (=%uus) exceeds maximum possible value that can be saved to SEGY format (=%uus)\n",
                      sampleIntUS, maxUnsignedShort);
    }
    else {
      writer->error("Data sample interval (=%uus) exceeds maximum possible value that can be saved to SEGY format (=%uus)\n",
                    sampleIntUS, maxUnsignedShort);
    }
  }
  if( (unsigned int)shdr->numSamples >= maxUnsignedShort ) {
    if( doForce ) {
      writer->warning("Number of data samples (=%u) exceeds maximum possible value that can be saved to SEGY format (=%u)\n",
                      (unsigned int)shdr->numSamples, maxUnsignedShort);
    }
    else {
      writer->error("Number of data samples (=%u) exceeds maximum possible value that can be saved to SEGY format (=%u)\n",
                    (unsigned int)shdr->numSamples, maxUnsignedShort);
    }
  }

  if( edef->isDebug() ) writer->line("Starting init phase of INPUT_SEGY...");

  param->getString( "filename", &filename );
  if( param->exists("nsamples" ) ) {
    param->getInt( "nsamples", &vars->numSamplesOut );
    if( vars->numSamplesOut <= 0 || vars->numSamplesOut > shdr->numSamples ) {
      writer->error("Inconsistent number of samples specified: %d. Actual number of samples in trace: %d", vars->numSamplesOut, shdr->numSamples);
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

  int hdr_mandatory = HDR_MANDATORY_DEFAULT;
  if( param->exists("hdr_mandatory") ) {
    string text;
    param->getString( "hdr_mandatory", &text );
    if( !text.compare("default") ) {
      hdr_mandatory = HDR_MANDATORY_DEFAULT;
    }
    else if( !text.compare("header") ) {
      hdr_mandatory = HDR_MANDATORY_HEADER;
    }
    else {
      writer->error("Unknown option '%s'", text.c_str() );
    }
  }

  bool removeMissingHdrs = true;
  if( param->exists("hdr_missing") ) {
    string text;
    param->getString( "hdr_missing", &text );
    if( !text.compare("ignore") ) {
      removeMissingHdrs = true;
    }
    else if( !text.compare("abort") ) {
      removeMissingHdrs = false;
    }
    else {
      writer->error("Unknown option '%s'", text.c_str() );
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
  bool rev_byte_order = false;
  if( param->exists("reverse_byte_order") ) {
    string text;
    param->getString( "reverse_byte_order", &text );
    if( !text.compare("yes") ) {
      rev_byte_order = true;
    }
  }

  //---------------------------------------------------------
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
    int length = (int)filename.length();
    if( length >= 3 ) {
      std::string ext = filename.substr(length-3);
      if( !ext.compare(".su") || !ext.compare(".SU") || !ext.compare(".Su") || !ext.compare(".sU") ) {
        isSUFormat = true;
      }
      else {
        isSUFormat = false;
      }
    }
  }

  //------------------------------------------------------------
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
      vars->hdr_mapping = csSegyHdrMap::SEGY_SU;
    }
    else if( !hdr_mappingName.compare("su_standard") ) {
      vars->hdr_mapping = csSegyHdrMap::SEGY_SU_ONLY;
    }
    else if( !hdr_mappingName.compare("none") ) {
      vars->hdr_mapping = csSegyHdrMap::NONE;
    }
    else {
      writer->error("Option for user parameter 'hdr_map' not recognised: %s", hdr_mappingName.c_str());
    }
  }
  else if( isSUFormat ) { // Set pre-set map for SU input data
    vars->hdr_mapping   = csSegyHdrMap::SEGY_SU;
  }

  // Pre-set header map selected by the user (otherwise defaulted to STANDARD SEGY)
  csSegyHdrMap userHdrMap( vars->hdr_mapping, false );

  if( param->exists("filename_hdrmap") ) {
    string filename_hdrmap;
    param->getString("filename_hdrmap", &filename_hdrmap );
    userHdrMap.readHdrMapExternalFile( filename_hdrmap );
    //    vars->hdrMap = new csSegyHdrMap( vars->hdr_mapping, false, filename_hdrmap );
  }

  bool isHdrMapDump = false;
  if( param->exists("dump_hdrmap") ) {
    param->getString( "dump_hdrmap", &yesno );
    yesno = toLowerCase( yesno );
    if( !yesno.compare( "yes" ) ) {
      isHdrMapDump = true;
      //      writer->line(" *** Dump of SEGY trace header map (excluding user specified non-standard headers) ***");
      //   userHdrMap.dump( writer->getFile() );
      //  writer->line("");
    }
    else if( !yesno.compare( "no" ) ) {
      isHdrMapDump = false;
      // Nothing
    }
    else {
      writer->line("%s: Parameter value not recognized: %s", edef->moduleName().c_str(), yesno.c_str() );
      env->addError();
    }
  }

  //---------------------------------------------------------
  if( param->exists( "ntraces_buffer" ) ) {
    param->getInt( "ntraces_buffer", &numTracesBuffer );
    if( numTracesBuffer <= 0 || numTracesBuffer > 999999 ) {
      writer->warning("Number of buffered traces out of range (=%d). Changed to 20.", numTracesBuffer);
      numTracesBuffer = 20;
    }
  }
  else {
    numTracesBuffer = 20;
  }

  //----------------------------------------------------
  //
  try {
    vars->segyWriter = new csSegyWriter( filename, numTracesBuffer, rev_byte_order, autoscale_hdrs, isSUFormat );
  }
  catch( csException& e ) {
    vars->segyWriter = NULL;
    writer->error("Error when generating SEGY writer object.\nSystem message: %s", e.getMessage() );
  }

  //-------------------------------------------------------
  // Read in all parameter lines for non-standard headers
  //

  csSegyHdrMap finalHdrMap( vars->hdr_mapping, false );
  finalHdrMap.removeAll();

  int num_id = 0;
  int id_nsamp        = num_id++;
  int id_sampint_us   = num_id++;
  int id_scalar_coord = num_id++;
  int id_scalar_elev  = num_id++;
  int id_scalar_stat  = num_id++;
  int id_data_type    = num_id++;
  int id_unit_coord   = num_id++;
  int id_time_code    = num_id++;;
  int val_trcType  = 0;
  int val_fold     = 0;
  int val_foldVert = 0;
  int id_trc_type = 0;
  int id_fold     = 0;
  int id_fold_vert = 0;

  if( param->exists("trc_type") ) {
    param->getInt("trc_type", &val_trcType);
    //    mandatoryHdrValues.set( value, id_trc_type );
    id_trc_type     = num_id++;
  }
  if( param->exists("fold") ) {
    param->getInt("fold", &val_fold);
    //    mandatoryHdrValues.set( value, id_fold );
    id_fold         = num_id++;
  }
  if( param->exists("fold_vert") ) {
    param->getInt("fold_vert", &val_foldVert);
    //    mandatoryHdrValues.set( value, id_fold_vert );
    id_fold_vert    = num_id++;
  }

  // Add mandatory trace headers, if defined in user specified, pre-defined, trace header map
  csVector<string> mandatoryHdrNames(num_id);
  csVector<int>    mandatoryHdrValues(num_id);

  mandatoryHdrNames.allocate( num_id, "" );
  mandatoryHdrValues.allocate( num_id, 0 );

  // Set default values for mandatory SEGY output headers:
  mandatoryHdrNames.set("nsamp",       id_nsamp);        mandatoryHdrValues.set( vars->numSamplesOut, id_nsamp );
  mandatoryHdrNames.set("sampint_us",  id_sampint_us);   mandatoryHdrValues.set( (int)(sampleIntUS), id_sampint_us );
  mandatoryHdrNames.set("scalar_coord", id_scalar_coord); mandatoryHdrValues.set( 100, id_scalar_coord );
  mandatoryHdrNames.set("scalar_elev", id_scalar_elev);  mandatoryHdrValues.set( 10, id_scalar_elev );
  mandatoryHdrNames.set("scalar_stat", id_scalar_stat);  mandatoryHdrValues.set( 1, id_scalar_stat );
  mandatoryHdrNames.set("data_type",   id_data_type);    mandatoryHdrValues.set( 1, id_data_type );
  mandatoryHdrNames.set("unit_coord",  id_unit_coord);   mandatoryHdrValues.set( 1, id_unit_coord );
  mandatoryHdrNames.set("time_code",   id_time_code);    mandatoryHdrValues.set( 4, id_time_code );

  if( id_trc_type > 0 ) {
    mandatoryHdrNames.set("trc_type",    id_trc_type);     mandatoryHdrValues.set( val_trcType, id_trc_type );
  }
  if( id_fold > 0 ) {
    mandatoryHdrNames.set("fold",        id_fold);         mandatoryHdrValues.set( val_fold, id_fold );
  }
  if( id_fold_vert > 0 ) {
    mandatoryHdrNames.set("fold_vert",   id_fold_vert);    mandatoryHdrValues.set( val_foldVert, id_fold_vert );
  }

  int value;
  if( param->exists("scalar_coord") ) {
    param->getInt("scalar_coord", &value);
    mandatoryHdrValues.set( value, id_scalar_coord );
  }
  if( param->exists("scalar_elev") ) {
    param->getInt("scalar_elev", &value);
    mandatoryHdrValues.set( value, id_scalar_elev );
  }
  if( param->exists("scalar_stat") ) {
    param->getInt("scalar_stat", &value);
    mandatoryHdrValues.set( value, id_scalar_stat );
  }
  if( param->exists("data_type") ) {
    param->getInt("data_type", &value);
    mandatoryHdrValues.set( value, id_data_type );
  }
  
  // Check if mandatory headers are defined in user selected pre-set header map
  // Each header that is not defined will be removed from the list of mandatory headers
  int hdrIndex = -1;
  for( int ihdr = mandatoryHdrNames.size()-1; ihdr >= 0; ihdr-- ) {
    string name = mandatoryHdrNames.at(ihdr);
    if( userHdrMap.contains(name,&hdrIndex) ) {
      csSegyHeaderInfo const* info = userHdrMap.header(hdrIndex);
      finalHdrMap.addHeader( *info );
    }
    else {
      mandatoryHdrNames.remove( ihdr );
      mandatoryHdrValues.remove( ihdr );
    }
  }

  //-------------------------------------------------------
  // Add all trace headers that currently exist and that are defined in the SEGY header map
  //
  writer->line("The following CSEIS trace headers will not be written to output file, unless specified by the user\nin a parameter 'header': (the following list will be empty if all current\nCSEIS trace headers are defined in the chosen SEGY header map)");
  writer->line("START_OF_HEADER_LIST");
  for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
    string name = hdef->headerName(ihdr);
    int hdrIndex = -1;
    if( userHdrMap.contains( name, &hdrIndex ) ) {
      finalHdrMap.addHeader( *userHdrMap.header(hdrIndex) );
    }
    else {
      writer->line("     %s", name.c_str());
    }
  }
  writer->line("END_OF_HEADER_LIST\n");

  int counterMissingHdrs = 0;
  writer->line("The following trace headers are defined in the SEGY trace header map but do not exist as CSEIS trace headers in the input data:");
  writer->line("START_OF_HEADER_LIST");
  for( int ihdr = 0; ihdr < userHdrMap.numHeaders(); ihdr++ ) {
    string name = userHdrMap.header(ihdr)->name;
    if( !hdef->headerExists(name) ) {
      writer->write("     %s", name.c_str());
      counterMissingHdrs += 1;
      for( int i = 0; i < mandatoryHdrNames.size(); i++ ) {
        if( !mandatoryHdrNames.at(i).compare(name) ) {
          writer->write("    (default/manually specified value will be used: %d)", mandatoryHdrValues.at(i) );
          counterMissingHdrs -= 1;  // Do not count as missing header because of default value
          break;
        }
      }
      writer->write("\n");
    }
  }
  writer->line("END_OF_HEADER_LIST\n");
 
  if( !removeMissingHdrs && counterMissingHdrs > 0 ) {
    writer->error("%d trace headers (without default) missing in input file.\nSpecify parameter 'hdr_missing  ignore' to avoid flow termination\n", counterMissingHdrs);
  }

  int nLines = param->getNumLines( "header" );
  cseis_geolib::csVector<std::string> valueList;
  for( int iline = 0; iline < nLines; iline++ ) {
    param->getAll( "header", &valueList, iline );
    int nValues = valueList.size();
    if( nValues < 3 ) {
      writer->error("Incorrect number of parameter values for parameter 'header'. Expected: At least 3, found: %d", nValues);
    }
    headerName = valueList.at(0);
    if( !hdef->headerExists( headerName ) ) {
      writer->line("Error: Trace header '%s' does not exist in data. Cannot write to SEGY output file", headerName.c_str() );
      env->addError();
      continue;
    }
    std::string desc;
    int byteLoc = atoi( valueList.at(1).c_str() ) - 1;  // -1 to reduce byte location to start at 0

    type_t outType  = TYPE_INT;
    int byteSize = 0;
    std::string text = valueList.at(2);
    if( !text.compare("int") ) {
      outType   = TYPE_INT;
      byteSize = 4;
    }
    else if( !text.compare("short") ) {
      outType   = TYPE_SHORT;
      byteSize = 2;
    }
    else if( !text.compare("ushort") ) {
      outType   = TYPE_USHORT;
      byteSize = 2;
    }
    else if( !text.compare("float") ) {
      outType   = TYPE_FLOAT;
      byteSize = 4;
    }
    else if( !text.compare("2") ) {
      outType   = TYPE_SHORT;
      byteSize = 2;
    }
    else if( !text.compare("4") ) {
      outType   = TYPE_INT;
      byteSize = 4;
    }
    else {
      writer->error("Incorrect output number type given for header '%s': %s.", headerName.c_str(), valueList.at(2).c_str());
    }
    int hdrIndex = -1;
    if( finalHdrMap.contains( headerName, &hdrIndex ) ) {
      if( overrideHdrs ) {
        finalHdrMap.removeHeader( hdrIndex );
        finalHdrMap.addHeader( headerName, byteLoc, byteSize, outType );
      }
      else {
        writer->line("Error: Cannot add user defined header '%s': Trace header is already defined in specified pre-set header map.", headerName.c_str());
        writer->line("       Termination of job can be avoided by specifying parameter 'hdr_duplicate  override'.");
        env->addError();
      }
    }
    else {
      finalHdrMap.addHeader( headerName, byteLoc, byteSize, outType );
    }
  }
  //----------------------------------------------------
  // Set CHAR header
  //
  nLines = param->getNumLines( "charhdr" );
  if( nLines > 40 ) {
    writer->line( "Found more than 40 lines (%d) for Segy char header. Only the first 40 lines will be used.", nLines );
    nLines = 40;
  }
  else if( nLines < 40 ) {
    writer->line( "Found less than 40 lines for Segy char header. The remaining lines will be filled up with blanks." );
  }

  std::string lineCharHdr;   // One line of char hdr
  char* charHdr = new char[cseis_segy::SIZE_CHARHDR+1];
  int counter = 0;
  for( int iline = 0; iline < nLines; iline++ ) {
    param->getStringAtLine( "charhdr", &lineCharHdr, iline );
    writer->line("Line %d: %s", iline, lineCharHdr.c_str() );
    int length = MIN( (int)lineCharHdr.length(), 80 );
    for( int k = 0; k < length; k++ ) {
      charHdr[k+counter] = lineCharHdr[k];
    }
    for( int k = length; k < 80; k++ ) {
      charHdr[k+counter] = ' ';   // Set remaining places to blanks
    }
    counter += 80;
  }
  if( counter < cseis_segy::SIZE_CHARHDR ) {
    memset( &charHdr[counter], ' ', cseis_segy::SIZE_CHARHDR - counter );
  }
  charHdr[cseis_segy::SIZE_CHARHDR] = '\0';

  if( edef->isDebug() ) {
    writer->line( charHdr );
  }
  
  //----------------------------------------------------
  // Set BIN header
  //
  int jobID = 0;
  if( param->exists("bin_jobid") ) {
    param->getInt( "bin_jobid", &jobID );
  }
  int lineNum = 0;
  if( param->exists("bin_linenum") ) {
    param->getInt( "bin_linenum", &lineNum );
  }
  int reelNum = 0;
  if( param->exists("bin_reelnum") ) {
    param->getInt( "bin_reelnum", &reelNum );
  }
  int bin_numTraces = 1;
  if( param->exists("bin_ntraces") ) {
    param->getInt( "bin_ntraces", &bin_numTraces );
  }
  int dataSampleFormat = cseis_segy::DATA_FORMAT_IEEE;
  if( param->exists("data_format") ) {
    string text;
    param->getString( "data_format", &text );
    if( !text.compare("ieee") ) {
      dataSampleFormat = cseis_segy::DATA_FORMAT_IEEE;
    }
    else if( !text.compare("ibm") ) {
      dataSampleFormat = cseis_segy::DATA_FORMAT_IBM;
    }
    else {
      writer->error("Unknown sample format: '%s'", text.c_str());
    }
  }
  cseis_geolib::csSegyBinHeader* binHdr = vars->segyWriter->binHdr();
  binHdr->jobID            = jobID;
  binHdr->lineNum          = lineNum;
  binHdr->reelNum          = reelNum;
  binHdr->numTraces        = bin_numTraces;
  binHdr->numAuxTraces     = 0;
  binHdr->sampleIntUS     = (unsigned short)(sampleIntUS);
  binHdr->sampleIntOrigUS = (unsigned short)(sampleIntUS);
  binHdr->numSamples       = vars->numSamplesOut;
  binHdr->numSamplesOrig   = vars->numSamplesOut;
  binHdr->dataSampleFormat = (short unsigned int)dataSampleFormat;
  binHdr->fold             = 1;
  binHdr->sortCode         = 1;
  binHdr->vertSumCode      = 1;
  binHdr->unitSystem       = 1;
  binHdr->polarity         = 1;
  binHdr->revisionNum      = 0000;
  binHdr->fixedTraceLengthFlag = 0;
  binHdr->numExtendedBlocks    = 0;
  
  if( isHdrMapDump ) {
    writer->line(" *** Dump of SEGY trace header map (excluding user specified non-standard headers) ***");
    finalHdrMap.dump( writer->getFile() );
    writer->line("");
  }
  //----------------------------------------------------
  try {
    vars->segyWriter->initialize( &finalHdrMap, charHdr );
  }
  catch( csException& e ) {
    vars->segyWriter = NULL;
    writer->error("Error when initializing SEGY writer object.\nSystem message: %s", e.getMessage() );
  }
  delete [] charHdr;

  //----------------------------------------------------
  // Get all trace headers to be written to SEGY file, set fields...
  //
  //  csSegyTraceHeader const* trcHdr = vars->segyWriter->getTraceHeader();
  //  int nHeaders = trcHdr->numHeaders();
  csSegyTraceHeader const* segyTrcHdr = vars->segyWriter->getTraceHeader();

  int numHeaders = segyTrcHdr->numHeaders();
  writer->line("\nNumber of trace headers written to SEGY output file: %d", numHeaders);
  vars->hdrIndexSegy = new int[numHeaders];
  vars->hdrTypeSegy  = new type_t[numHeaders];
  for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
    string const& name = segyTrcHdr->headerName(ihdr);
    writer->line(" #%2d  %s", ihdr+1, name.c_str());
    if( hdef->headerExists(name) ) {
      vars->hdrIndexSegy[ihdr] = hdef->headerIndex( name );
      vars->hdrTypeSegy[ihdr]  = hdef->headerType( name );
    }
    else {  // Must be mandatory header...
      vars->hdrIndexSegy[ihdr] = -1;
    }
  }

  for( int ihdr = 0; ihdr < mandatoryHdrNames.size(); ihdr++ ) {
    string name = mandatoryHdrNames.at(ihdr);
    int hdrIndex = segyTrcHdr->headerIndex(name);
    if( hdrIndex < 0 ) writer->error("Program bug! Wrong mandatory header index....");

    if( hdef->headerExists(name) ) {
      if( hdr_mandatory == HDR_MANDATORY_DEFAULT ) {
        //      writer->line("Info:  Mandatory header '%s' exists in input trace header but a default or manually specified value (%d) will be used.", name.c_str(), mandatoryHdrValues.at(ihdr) );
        vars->hdrIndexSegy[hdrIndex] = -1;
      }
      else if( hdr_mandatory == HDR_MANDATORY_HEADER && param->exists(name.c_str()) ) {
        writer->line("Info:  Mandatory header '%s' exists in input trace header but was manually specified by user as %d.", name.c_str(), mandatoryHdrValues.at(ihdr) );
        vars->hdrIndexSegy[hdrIndex] = -1;
      }
    }
    vars->segyWriter->setIntValue( hdrIndex, mandatoryHdrValues.at(ihdr) );
  }

  writer->line("");
  writer->line("  File name:            %s", filename.c_str());
  writer->line("  Sample interval [ms]: %f", (float)(binHdr->sampleIntUS)/1000.0);
  writer->line("  Number of samples:    %d", binHdr->numSamples );
  writer->line("  Sample data format:   %d", binHdr->dataSampleFormat );
  writer->line("");

  //  shdr->sampleInt = (float)(binHdr->sampleIntUS)/1000.0; // Commented out: Do not update sample interval of Seaseis trace flows
  vars->traceCounter = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_output_segy_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;

  csTrace* trace = traceGather->trace(0);


  //----------------------------------------------------
  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    try {
      vars->segyWriter->openFile();
      vars->segyWriter->freeCharBinHdr();  // Don't need these headers anymore -> free memory
    }
    catch( csException& e ) {
      writer->error("Error when opening SEGY output file.\nSystem message: %s", e.getMessage() );
    }
  }


  csTraceHeader* trcHdr = trace->getTraceHeader();

  csSegyTraceHeader* segyTrcHdr = vars->segyWriter->getTraceHeader(); 
  int nHeadersSegy = segyTrcHdr->numHeaders();
  for( int ihdr = 0; ihdr < nHeadersSegy; ihdr++ ) {
    int hdrIdOut = vars->hdrIndexSegy[ihdr];
    if( hdrIdOut < 0 ) continue;  // Skip headers that have already been set (mandatory headers)
    if( vars->hdrTypeSegy[ihdr] == TYPE_FLOAT ) {
      segyTrcHdr->setFloatValue( ihdr, trcHdr->floatValue( hdrIdOut )  );
      //      if( edef->isDebug() ) writer->line( "Segy header %2d: %d   '%s'  %f", ihdr, segyWriter->hdrIntValue(ihdr), segyWriter->headerName(ihdr), segyWriter->hdrFloatValue(ihdr) );
    }
    else if( vars->hdrTypeSegy[ihdr] == TYPE_DOUBLE ) {
      segyTrcHdr->setDoubleValue( ihdr, trcHdr->doubleValue(hdrIdOut) );
      //      if( edef->isDebug() ) writer->line( "Segy header %2d: %d   '%s'  %f", ihdr, segyWriter->hdrIntValue(ihdr), segyWriter->headerName(ihdr), segyWriter->hdrDoubleValue(ihdr) );
    }
    else if( vars->hdrTypeSegy[ihdr] == TYPE_INT ) {
      segyTrcHdr->setIntValue( ihdr, trcHdr->intValue(hdrIdOut) );
      //      if( edef->isDebug() ) writer->line( "Segy header %2d: %d   '%s'  %d", ihdr, segyWriter->hdrIntValue(ihdr), segyWriter->headerName(ihdr), segyWriter->hdrIntValue(ihdr) );
    }
    else if( vars->hdrTypeSegy[ihdr] == TYPE_INT64 ) {
      segyTrcHdr->setIntValue( ihdr, (int)trcHdr->int64Value(hdrIdOut) );
    }
    else { // if( vars->hdrTypeSegy[ihdr] == TYPE_STRING ) {
      segyTrcHdr->setStringValue( ihdr, trcHdr->stringValue(hdrIdOut) );
    }
  }

  float* samples = trace->getTraceSamples();
  try {
    vars->segyWriter->writeNextTrace( (byte_t*)samples, vars->numSamplesOut );
  }
  catch( csException& exc ) {
    writer->error( "System message: %s" ,exc.getMessage() );
  }
  vars->traceCounter++;

  return;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_output_segy_( csParamDef* pdef ) {

  pdef->setModule( "OUTPUT_SEGY", "Output SEGY/SU data", "Writes SEGY or Seismic Unix (SU) data to disk file" );
  pdef->addParam( "filename", "Output file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );

  pdef->addParam( "auto_scale", "Automatically scale SEGY headers. This applies to standard coordinate, elevation, and statics trace headers", NUM_VALUES_FIXED,
                  "Automatically scale SEGY headers and save SEGY scalar values at standard byte positions.");
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Automatically scale trace header values before writing to output file" );
  pdef->addOption( "no", "Do not apply any scalars. Coordinates, elevations and statics values are written to output file as-is, scalars are set to 1." );

  pdef->addParam( "hdr_map", "Pre-set SEGY trace header mapping", NUM_VALUES_FIXED, "Use option 'none' to avoid encoding any pre-set trace headers. Specify user parameter 'dump_hdrmap  yes' to make a printout of the specified pre-set map");
  pdef->addValue( "standard", VALTYPE_OPTION );
  pdef->addOption( "standard", "Standard SEGY header mapping, SEGY revision 1" );
  pdef->addOption( "none", "Do not use any pre-set trace header map. All SEGY trace headers that shall be encoded are given either by user parameter 'header' (non-standard trace header) or 'filename_hdrmap' (external ASCII file)" );
  pdef->addOption( "obc", "OBC 4C header mapping", "In principle SEGY revision 1, plus additional headers" );
  pdef->addOption( "send", "SEND header mapping" );
  pdef->addOption( "armss", "ARMSS header mapping" );
  pdef->addOption( "psegy", "PSEGY header mapping" );
  pdef->addOption( "node", "Node header mapping" );
  pdef->addOption( "node_old", "Node header mapping (obsolete)" );
  pdef->addOption( "su", "Seismic Unix (SU) header mapping");
  pdef->addOption( "su_standard", "Same as option 'su', but assume all trace headers exist only with SU standard naming, not Seaseis standard naming");

  pdef->addParam( "filename_hdrmap", "File name containing trace header map setup", NUM_VALUES_FIXED,
                  "All trace headers already defined in the specified pre-set trace header map at the same byte locations will be overriden. If no pre-set trace headers are required, specify user parameter 'hdr_map  none'");
  pdef->addValue( "", VALTYPE_STRING, "Input file name containing trace header map definition",
                  "File format:  BYTELOC  TYPE_SEGY  TYPE_CEIS  NAME   DESC, where\n"\
                  "BYTELOC:    Byte location, starting at 1 for first byte,\n"\
                  "TYPE_SEGY:  Header type as stored in SEGY file, see user parameter 'header' for options,\n"\
                  "TYPE_CSEIS: Header type of Seaseis trace header, see user parameter 'header' for options,\n"\
                  "NAME:       Trace header name,\n"\
                  "DESC:       Description of Seaseis trace header, enclosed in double-quotes.");

  pdef->addParam( "header", "Set user specified SEGY trace headers", NUM_VALUES_VARIABLE);
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Byte location (starting at 1)" );
  pdef->addValue( "int", VALTYPE_OPTION, "Trace header output type", "This is the number type in which the trace header is stored/formatted in the output file's trace header");
  pdef->addOption( "int", "Integer (4 byte)" );
  pdef->addOption( "short", "Short (2 byte)" );
  pdef->addOption( "ushort", "Unsigned short (2 byte)" );
  pdef->addOption( "4", "Integer (4 byte)", "...for backward compatibility" );
  pdef->addOption( "2", "Short (2 byte)", "...for backward compatibility" );
  //  pdef->addOption( "float", "Float (4 byte)" );   // Not supported yet. Need to implement float conversion function first

  pdef->addParam( "data_format", "Floating point format for data samples", NUM_VALUES_FIXED);
  pdef->addValue( "ieee", VALTYPE_OPTION );
  pdef->addOption( "ieee", "Use 32-byte IEEE floating point format" );
  pdef->addOption( "ibm", "Use 32-byte IBM floating point format" );

  pdef->addParam( "hdr_mandatory", "How to treat mandatory SEGY trace headers", NUM_VALUES_FIXED,
                  "Mandatory SEGY trace headers are: data_type, trc_type, nsamp, sampint_us, scalar_coord, scalar_elev, scalar_stat, fold, fold_vert. Each of these can be set by a) Using a default value, b) Creating an CSEIS header with this name before OUTPUT_SEGY, c) Setting a header manually in one of the user parameters below" );
  pdef->addValue( "default", VALTYPE_OPTION );
  pdef->addOption( "default", "Set mandatory trace headers to default values, except for those that are manually specified in one of the according user parameters below", "Check log file to see which mandatory trace headers were defaulted or not" );
  pdef->addOption( "header", "If mandatory trace header exists in input data, use it to override default value. Otherwise use default. Manually specified value overrides all" );

  pdef->addParam( "scalar_coord", "SEGY 'coordinate scalar'. All SEGY standard coordinates will be scaled by this value", NUM_VALUES_FIXED, "Note that the given value will be stored as -value in the SEGY trace header, in accordance to the SEGY definition, in order to enable inverse scaling when reading this file" );
  pdef->addValue( "100", VALTYPE_NUMBER, "Scalar to be applied to all SEGY standard coordinate headers" );

  pdef->addParam( "scalar_elev", "SEGY 'elevation scalar'. All SEGY standard elevations will be scaled by this value", NUM_VALUES_FIXED, "See scalar_coord for further information." );
  pdef->addValue( "10", VALTYPE_NUMBER, "Scalar to be applied to all SEGY standard elevation headers" );

  pdef->addParam( "scalar_stat", "SEGY 'static scalar'. All SEGY standard statics will be scaled by this value", NUM_VALUES_FIXED, "See scalar_coord for further information." );
  pdef->addValue( "1", VALTYPE_NUMBER, "Scalar to be applied to all SEGY standard statics headers" );

  pdef->addParam( "fold", "Set SEGY 'fold' trace header for all traces", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "SEGY 'fold' header" );

  pdef->addParam( "fold_vert", "Set SEGY 'vertical fold' trace header for all traces", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "SEGY 'vertical fold' header" );

  pdef->addParam( "data_type", "Set SEGY 'data type' trace header for all traces", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "SEGY 'data type' header" );

  pdef->addParam( "trc_type", "Set SEGY 'trace type' trace header for all traces", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "SEGY 'trace type' header" );

  pdef->addParam( "hdr_duplicate", "How to treat user specified trace headers (parameter 'header') that are already defined in the chosen trace header map", NUM_VALUES_FIXED );
  pdef->addValue( "override", VALTYPE_OPTION );
  pdef->addOption( "override", "Write header to user specified byte location only, i.e. override standard byte location" );
  pdef->addOption( "abort", "Abort if user specified header is already defined in trace header map" );
  //  pdef->addOption( "both", "Write header to both byte locations: 1) The standard byte location defined in the trace header map, 2) The user specified byte location" );

  pdef->addParam( "hdr_missing", "How to deal with missing SEGY trace headers that are not defined in input data", NUM_VALUES_FIXED );
  pdef->addValue( "ignore", VALTYPE_OPTION );
  pdef->addOption( "ignore", "Ignore missing trace headers. Set to default value (usually 0)" );
  pdef->addOption( "abort", "Abort flow when trace headers required by the specified trace header map are missing in input data." );

  pdef->addParam( "charhdr", "80 characters for SEGY 3200-byte character header", NUM_VALUES_FIXED, "Specify up to 40 lines");
  pdef->addValue( "C00", VALTYPE_STRING, "Line 00 for SEGY 3200-byte character header" );

  //  pdef->addParam( "ntraces", "(Maximum) number of traces to output", NUM_VALUES_FIXED);
  //  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to output" );

  pdef->addParam( "bin_jobid", "Job ID, set in binary header", NUM_VALUES_FIXED);
  pdef->addValue( "0", VALTYPE_NUMBER, "Job ID" );
  pdef->addParam( "bin_linenum", "Line number, set in binary header", NUM_VALUES_FIXED);
  pdef->addValue( "0", VALTYPE_NUMBER, "Line number" );
  pdef->addParam( "bin_reelnum", "Reel number, set in binary header", NUM_VALUES_FIXED);
  pdef->addValue( "0", VALTYPE_NUMBER, "Reel number" );
  pdef->addParam( "bin_ntraces", "Number of traces, set in binary header", NUM_VALUES_FIXED);
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of traces" );

  pdef->addParam( "dump_hdrmap", "Dump specified trace header map (=byte locations) to log file", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Dump trace header map" );
  pdef->addOption( "no", "Dump trace header map" );

  pdef->addParam( "ntraces_buffer", "Number of traces to buffer before write operation", NUM_VALUES_FIXED,
                  "Writing a large number of traces at once enhances performance, but requires more memory" );
  pdef->addValue( "20", VALTYPE_NUMBER, "Number of traces to buffer before writing" );

  pdef->addParam( "sample_int", "Override output sample interval", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Sample interval [ms]" );

  pdef->addParam( "reverse_byte_order", "Reverse byte order in output file (endian byte swapping)", NUM_VALUES_FIXED, "Setting this to 'yes' means that the output SEGY file will be written in Little endian byte order" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Reverse byte order for output file" );
  pdef->addOption( "no", "Do not reverse byte order" );

  pdef->addParam( "force", "Force writing file", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Force writing SEGY output file" );
  pdef->addOption( "no", "Abort if output SEGY would become inconsistent" );

  pdef->addParam( "su_format", "Write output data in SU (Seismic Unix) format?", NUM_VALUES_FIXED,
                  "Use this parameter to override the default format which is determined by the (first) file's extension: *.su or *.SU == SU file" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Force output data to be written in SU format, independent of file extension" );
  pdef->addOption( "no", "Force output data to be written in SEGY format, independent of file extension" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_output_segy_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_output_segy::VariableStruct* vars = reinterpret_cast<mod_output_segy::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_output_segy_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_output_segy::VariableStruct* vars = reinterpret_cast<mod_output_segy::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->segyWriter != NULL ) {
    delete vars->segyWriter;
    vars->segyWriter = NULL;
  }
  delete [] vars->hdrIndexSegy; vars->hdrIndexSegy = NULL;
  delete [] vars->hdrTypeSegy; vars->hdrTypeSegy = NULL;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_output_segy_( csParamDef* pdef ) {
  params_mod_output_segy_( pdef );
}
extern "C" void _init_mod_output_segy_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_output_segy_( param, env, writer );
}
extern "C" bool _start_exec_mod_output_segy_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_output_segy_( env, writer );
}
extern "C" void _exec_mod_output_segy_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_output_segy_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_output_segy_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_output_segy_( env, writer );
}

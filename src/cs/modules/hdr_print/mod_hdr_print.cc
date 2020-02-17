/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFileUtils.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: HDR_PRINT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_hdr_print {
  struct VariableStruct {
    bool dumpAll;
    bool dumpShdr;
    std::string filename;
    int* indexHdr;
    type_t* typeHdr;
    std::string* hdrNames;
    std::string* titleFormats;
    std::string* printFormats;
    int nHeaders;
    int countTraces;
    int trcIncTitle;
    bool isExternalFile;
    FILE* f_out;
    bool auto_format;
    bool isFirstCall;
  };
}
using mod_hdr_print::VariableStruct;

#define CORRECT_CHAR 0
#define TRC_INC_TITLE 100
#define NONE -1
#define ONCE 0

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_hdr_print_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  std::string name;
  std::string text;
  csVector<std::string> valueList;


  vars->dumpAll     = false;
  vars->dumpShdr    = false;
  vars->nHeaders    = 0;
  vars->countTraces = 0;
  vars->trcIncTitle = NONE;
  vars->isExternalFile = false;
  vars->auto_format    = true;

  vars->indexHdr = NULL;
  vars->typeHdr  = NULL;
  vars->hdrNames = NULL;
  vars->printFormats = NULL;
  vars->titleFormats = NULL;
  vars->filename = "";
  vars->f_out    = NULL;
  vars->isFirstCall = true;

  //---------------------------------------------------------
  //
  if( param->exists( "dump_all" ) ) {
    param->getString( "dump_all", &text );
    if( !text.compare("yes") ) {
      vars->dumpAll = true;
      vars->trcIncTitle = 1;
    }
    else if( !text.compare("no") ) {
      vars->dumpAll = false;
    }
    else {
      writer->error("Unknown option for parameter 'dump_all': '%s'", text.c_str());
    }
  }
  if( param->exists( "dump_shdr" ) ) {
    string yesno;
    param->getString("dump_shdr",&yesno);
    yesno = toLowerCase(yesno);
    if( !yesno.compare("yes") ) {
      vars->dumpShdr = true;
      shdr->dump( writer->getFile() );
    }
  }

  if( param->exists( "title" ) ) {
    std::string text;
    param->getString( "title", &text );
    if( !text.compare("none") ) {
      vars->trcIncTitle = NONE;
    }
    else if( !text.compare("once") ) {
      vars->trcIncTitle = ONCE;
    }
    else if( !text.compare("trc_inc") ) {
      if( param->getNumValues("title") < 2 ) {
        writer->error("Missing input for option 'trc_inc': Specify trace increment (e.g. title trc_inc 100).");
      }
      param->getInt( "title", &vars->trcIncTitle, 1 );
      if( vars->trcIncTitle < 0 ) {
        writer->error("Incorrect trace increment for option 'trc_inc': Must be larger than 0 (e.g. title trc_inc 100)");
      }
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  if( vars->dumpAll ) {
    vars->nHeaders = hdef->numHeaders();
    for( int i = 0; i < vars->nHeaders; i++ ) {
      valueList.insertEnd( hdef->headerName(i) );
    }
  }
  else if( param->exists("header") ) {
    int nLines = param->getNumLines( "header" );
    if( nLines > 1 ) {
      writer->line( "More than one line encountered for user parameter '%s'. Only one line is supported.", "header" );
      //      throw( csException("More than one line encountered for user parameter '%s'. Only one line is supported.", "header" ) ); // TEMP
      env->addError();
    }

    valueList.clear();
    param->getAll( "header", &valueList );
    if( valueList.size() < 1 ) {
      writer->error("Wrong number of arguments for user parameter '%s'. Expected: 1, found: %d.", "header", valueList.size());
    }
    vars->nHeaders = valueList.size();
  }
  else if( vars->dumpShdr ) {
    return;
  }
  else {
    writer->error("None of the mandatory parameters specified");
  }

  vars->indexHdr  = new int[vars->nHeaders];
  vars->typeHdr   = new type_t[vars->nHeaders];
  vars->hdrNames  = new std::string[vars->nHeaders];
  for( int k = 0; k < vars->nHeaders; k++ ) {
    name = valueList.at(k);
    vars->hdrNames[k] = name;
    if( hdef->headerExists( name.c_str() ) ) {
      vars->indexHdr[k]  = hdef->headerIndex( name.c_str() );
      vars->typeHdr[k]   = hdef->headerType( name.c_str() );
    }
    else {
      writer->line("Unknown trace header: %s", name.c_str());
      //      throw( csException("Unknown trace header: %s", name.c_str()) ); // TEMP
      env->addError();
    }
  }

  //-----------------------------------------------------------
  // Set print format
  //
  vars->printFormats = new std::string[vars->nHeaders];
  vars->titleFormats = new std::string[vars->nHeaders];
  valueList.clear();

  bool auto_correct = false;
  if( param->exists("auto_correct") ) {
    param->getString("auto_correct", &text);
    if( !text.compare("yes") ) {
      auto_correct = true;
    }
    else if( !text.compare("no") ) {
      auto_correct = false;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  else {
    auto_correct = false;
  }

  // Dump all trace headers --> Automatically set print formats
  if( vars->dumpAll ) {
    for( int i = 0; i < vars->nHeaders; i++ ) {
      type_t type = hdef->headerType(i);
      switch( type ) {
        case TYPE_INT:
          text = "%-15d";
          break;
        case TYPE_INT64:
          text = "%-15lld";
          break;
        case TYPE_DOUBLE:
        case TYPE_FLOAT:
          text = "%-15e";
          break;
        case TYPE_VECTOR:
          text = "%-15e";
          break;
        case TYPE_VECTOR_X:
        case TYPE_VECTOR_Y:
        case TYPE_VECTOR_Z:
          text = "%-15e";
          break;
      case TYPE_STRING:
          text = "%-15s";
          break;
        case TYPE_CHAR:
          text = "%-15c";
          break;
      }
      if( (i+1) % 10 == 0 ) text = "\n" + text;
      vars->printFormats[i] = text;
      //      valueList.insertEnd( text );
    }
  }
  // Header list was specified --> Read in print formats, or set automatically
  else {
    // Read in user specified print formats
    if( param->exists( "format" ) ) {
      vars->auto_format = false;
      param->getAll( "format", &valueList );
      if( valueList.size() != vars->nHeaders ) {
        writer->error("Wrong number of parameters for option '%s'. Expected: %d, found: %d.", "format", vars->nHeaders, valueList.size());
      }
      for( int k = 0; k < vars->nHeaders; k++ ) {
        vars->printFormats[k] = valueList.at(k);
      }
    }

    if( param->exists("auto_format") ) {
      if( param->exists("format") ) {
        writer->warning("Specified automatic formatting ('auto_format'=yes), and at the same manual formatting (user parameter 'format').\nManual formatting will be ignored.");
      }
      param->getAll("auto_format", &valueList );
      text = valueList.at(0);
      if( !text.compare("yes") ) {
        vars->auto_format = true;
      }
      else if( !text.compare("no") ) {
        vars->auto_format = false;
      }
      else {
        writer->error("Unknown option: %s", text.c_str());
      }
    }
    // Automatically set print formats
    if( vars->auto_format ) {
      string text_format_length("12");
      if( param->exists("auto_format") && param->getNumValues("auto_format") > 1 ) {
        text_format_length = valueList.at(1);
      }
      int num_format_length = atoi( text_format_length.c_str() );
      if( fabs((float)num_format_length) > 99 ) writer->warning("User parameter 'auto_format': Very long format length encountered: %d", num_format_length);

      int minlen = 3 + num_format_length;
      char* tmpString = new char[2+minlen];
      int position = 0;
      for( int k = 0; k < vars->nHeaders; k++ ) {
        switch( vars->typeHdr[k] ) {
          case TYPE_FLOAT:
          case TYPE_DOUBLE:
          case TYPE_VECTOR_X:
          case TYPE_VECTOR_Y:
          case TYPE_VECTOR_Z:
            sprintf( tmpString, "%%%sf", text_format_length.c_str() );
            position = minlen;
            break;
          case TYPE_VECTOR:
            sprintf( tmpString, "%%%sf", text_format_length.c_str() );
            position = minlen;
            break;
          case TYPE_INT:
            sprintf( tmpString, "%%%dd", num_format_length );
            position = minlen;
            break;
          case TYPE_INT64:
            sprintf( tmpString, "%%%dlld", num_format_length );
            position = minlen+2;
            break;
          default:  // TYPE_CHAR, TYPE_STRING
            sprintf( tmpString, "%%%ss", text_format_length.c_str() );
            position = minlen;
        }
        tmpString[position-1] = '\0';
        vars->printFormats[k] = tmpString;
      }
      delete [] tmpString;
    }
  }

  text == "";
  int ic;
  int counter;
  for( int k = 0; k < vars->nHeaders; k++ ) {
// Convert format stringlets into text format ( %d --> %s, %f --> %s, %-12.3f --> %-12s )
// Replace all other characters with white space
    text = vars->printFormats[k];
    counter = 0;
    int j = 0;
    while( j < (int)text.length() ) {
      if( text[j] == '%' ) {  // Start of format description: '%'
        do {
          j++;
          counter++;
          ic = (int)text[j];
          if( text[j] == '.' ) {  // Remove all letters behind and including '.' in formats such as "%12.2f"
            do {
              j++;
              ic = (int)text[j];
            } while( j < (int)text.length() && (ic >= 48 && ic <= 57) );
          }
        } while( j < (int)text.length() && ( ic < 65 || ic > 122 || ( ic > 90 && ic < 97 ) ) );  // Stop when trailing letter is found
        char correct_format_char = CORRECT_CHAR;
        switch( vars->typeHdr[k] ) {
        case TYPE_FLOAT:
        case TYPE_DOUBLE:
        case TYPE_VECTOR_X:
        case TYPE_VECTOR_Y:
        case TYPE_VECTOR_Z:
          if( text[j] != 'f' && text[j] != 'e' ) {
            correct_format_char = 'f';
          }
          break;
        case TYPE_VECTOR:
          if( text[j] != 'f' && text[j] != 'e' ) {
            correct_format_char = 'f';
          }
          break;
        case TYPE_INT:
          if( text[j] != 'd' ) {
            correct_format_char = 'd';
          }
          break;
        case TYPE_CHAR:
        case TYPE_STRING:
          if( text[j] != 's' ) {
            correct_format_char = 's';
          }
          break;
        case TYPE_INT64:
          if( text[j] != 'l' || text[j+1] != 'l' || text[j+2] != 'd') {
            correct_format_char = 'l';
          }
          j += 2;  // Add 2 to running index for full format 'lld'
          break;
        default:
          break;
        }
        if( correct_format_char != CORRECT_CHAR ) {
          if( vars->typeHdr[k] == TYPE_INT64 ) {
            writer->error("Incorrect print format for trace header #%d: Expected '%%lld', user specified '%c...'",
                       k+1, correct_format_char );
          }
          else if( auto_correct ) {
            writer->line( "Changed incorrect print format for trace header #%d: Expected '%c', user specified '%c'",
                       k+1, correct_format_char, text[j] );
            vars->printFormats[k][j] = correct_format_char;
          }
          else {
            writer->warning("Incorrect print format for trace header #%d: Expected '%c', user specified '%c'",
                         k+1, correct_format_char, text[j] );
            env->addError();
          }
        }
        text[counter++] = 's';
        j++;
      }  // if '%'
      else if( text[j] == '\0' ) {
        j++;  // Nothing
      }
      else {
        if( text[j] == '\n' ) { // && j+1 < text.length() && text[j+1] == 'n' ) {
          text[counter++] = '\n';
        }
        else {
          text[counter++] = ' ';  // replace other characters with white space
        }
        j++;
      }
    }  // END: outer while loop
    vars->titleFormats[k] = text.substr( 0, counter );
    if( edef->isDebug() ) writer->line("Title format #%d: '%s'", k, vars->titleFormats[k].c_str() );
    if( edef->isDebug() ) writer->line("Print format #%d: '%s'", k, vars->printFormats[k].c_str() );
  }

  writer->write("\nPrint format string (C-style): '");
  for( int k = 0; k < vars->nHeaders; k++ ) {
    writer->write("%s ", vars->printFormats[k].c_str());
  }
  writer->write("'\n");
  writer->write("Title format string (C-style): '");
  for( int k = 0; k < vars->nHeaders; k++ ) {
    writer->write("%s ", vars->titleFormats[k].c_str());
  }
  writer->write("'\n\n");

  if( param->exists( "filename" ) ) {
    param->getString( "filename", &vars->filename );
    if( !csFileUtils::createDoNotOverwrite( vars->filename ) ) {
      writer->error("Unable to open output file %s. Wrong path name..?", vars->filename.c_str() );
    }
    vars->isExternalFile = true;
  }
  else {
    vars->f_out = writer->getFile();
    vars->isExternalFile = false;
  }

  vars->countTraces = vars->trcIncTitle;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_hdr_print_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csTrace* trace = traceGather->trace(0);

//----------------------------
  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    if( vars->isExternalFile ) {
      if( (vars->f_out = fopen( vars->filename.c_str(), "w" )) == (FILE*) NULL ) {
        writer->error("Could not open file: '%s'", vars->filename.c_str());
      }
    }
  }


  csTraceHeader* trcHdr = trace->getTraceHeader();

  if( vars->trcIncTitle < 0 ) {
    vars->countTraces = 1;
  }
  else if( vars->trcIncTitle > 0 || vars->countTraces <= 0 ) {
    if( vars->countTraces >= vars->trcIncTitle ) {
      std::string text;
      text = "";
      for( int i = 0; i < vars->nHeaders; i++ ) {
        fprintf( vars->f_out, vars->titleFormats[i].c_str(), vars->hdrNames[i].c_str() );
        fprintf( vars->f_out, " " );
        text += vars->titleFormats[i];
      }
      fprintf( vars->f_out, "\n" );
      vars->countTraces = 1;
    }
    else {
      vars->countTraces += 1;
    }
  }

  for( int i = 0; i < vars->nHeaders; i++ ) {
    //    char const* str = vars->printFormats[i].c_str();
    type_t type = vars->typeHdr[i];
    int index   = vars->indexHdr[i];
    if( type == TYPE_FLOAT ) {
      fprintf( vars->f_out, vars->printFormats[i].c_str(), trcHdr->floatValue(index) );
    }
    else if( type == TYPE_DOUBLE ) {
      fprintf( vars->f_out, vars->printFormats[i].c_str(), trcHdr->doubleValue(index) );
    }
    else if( type == TYPE_VECTOR_X || type == TYPE_VECTOR_Y || type == TYPE_VECTOR_Z ) {
      fprintf( vars->f_out, vars->printFormats[i].c_str(), trcHdr->vectorValue(index,type) );
    }
    else if( type == TYPE_VECTOR ) {
      csPoint3D p = trcHdr->vectorValue(index);
      fprintf( vars->f_out, vars->printFormats[i].c_str(), p.x );
      fprintf( vars->f_out, " " );
      fprintf( vars->f_out, vars->printFormats[i].c_str(), p.y );
      fprintf( vars->f_out, " " );
      fprintf( vars->f_out, vars->printFormats[i].c_str(), p.z );
    }
    else if( type == TYPE_INT ) {
      fprintf( vars->f_out, vars->printFormats[i].c_str(), trcHdr->intValue(index) );
    }
    else if( type == TYPE_INT64 ) {
      fprintf( vars->f_out, vars->printFormats[i].c_str(), trcHdr->int64Value(index) );
    }
    else if( type == TYPE_STRING ) {
      fprintf( vars->f_out, vars->printFormats[i].c_str(), trcHdr->stringValue(index).c_str() );
    }
    else {
      fprintf(stderr,"Unknown trace header: %d\n", type);
      //writer->error("Encountered unknown trace header type, code: %d", type);
    }
    fprintf( vars->f_out, " " );
  }

  if( vars->nHeaders > 0 ) fprintf( vars->f_out, "\n" );
  if( vars->dumpAll ) fprintf( vars->f_out, "\n\n" );

  return;
}
//********************************************************************************
//
//
void params_mod_hdr_print_( csParamDef* pdef ) {
  pdef->setModule( "HDR_PRINT", "Print trace header values", "Print trace header values to log or external file" );

  pdef->addParam( "filename", "Output file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );

  pdef->addParam( "header", "List of trace header names to print", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "List of trace header names" );

  pdef->addParam( "format", "Format codes", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "List of C-style format codes to use for trace header printout" );

  pdef->addParam( "title", "Print title line showing trace header names", NUM_VALUES_VARIABLE );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "once", "Print title line once at the beginning" );
  pdef->addOption( "none", "Do not print title line" );
  pdef->addOption( "trc_inc", "Print title line with header names every n'th trace.", "Specify trace increment in next user parameter" );
  pdef->addValue( "100", VALTYPE_NUMBER, "Trace increment in number of traces (only for option 'trc_inc')" );

  pdef->addParam( "auto_correct", "Automatically correct mismatches between user specified print format and header type", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Automatically correct print format" );
  pdef->addOption( "no", "Abort when mismatch has been found" );

  pdef->addParam( "auto_format", "Automatically format output", NUM_VALUES_VARIABLE );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Automatically format output" );
  pdef->addOption( "no", "Manually format output" );
  pdef->addValue( "12", VALTYPE_NUMBER, "Number of columns that each header spans. Examples: 12 -> '%12d'/'%12f', 0 -> '%d'/'%f', 12.2 --> '%12f'/'%12.2f', depending on type" );

  pdef->addParam( "dump_all", "Dump all headers", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Dump all headers. This overrides any specification of individual headers" );
  pdef->addOption( "no", "Do not create a dump of all headers" );

  pdef->addParam( "dump_shdr", "Dump super headers", NUM_VALUES_FIXED, "Print out super header during module init phase" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Print super header" );
  pdef->addOption( "no", "Do not print super header" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_hdr_print_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_hdr_print::VariableStruct* vars = reinterpret_cast<mod_hdr_print::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_hdr_print_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_hdr_print::VariableStruct* vars = reinterpret_cast<mod_hdr_print::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->isExternalFile && vars->f_out != NULL ) {
    fclose( vars->f_out );
  }
  if( vars->indexHdr ) {
    delete [] vars->indexHdr; vars->indexHdr = NULL;
  }
  if( vars->typeHdr ) {
    delete [] vars->typeHdr;
    vars->typeHdr = NULL;
  }
  if( vars->hdrNames ) {
    delete [] vars->hdrNames;
    vars->hdrNames = NULL;
  }
  if( vars->printFormats ) {
    delete [] vars->printFormats;
    vars->printFormats = NULL;
  }
  if( vars->titleFormats ) {
    delete [] vars->titleFormats;
    vars->titleFormats = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_hdr_print_( csParamDef* pdef ) {
  params_mod_hdr_print_( pdef );
}
extern "C" void _init_mod_hdr_print_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_hdr_print_( param, env, writer );
}
extern "C" bool _start_exec_mod_hdr_print_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_hdr_print_( env, writer );
}
extern "C" void _exec_mod_hdr_print_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_hdr_print_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_hdr_print_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_hdr_print_( env, writer );
}

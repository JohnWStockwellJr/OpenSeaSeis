/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdlib>
#include <cstring>
#include <stdarg.h>
#include <iostream>
#include <string>
#include "cseis_defines.h"
#include "csTrace.h"
#include "csLogWriter.h"
#include "csTraceHeaderData.h"
#include "csMethodRetriever.h"
#include "csUserConstant.h"
#include "csRunManager.h"
#include "csParamDef.h"
#include "csMemoryPoolManager.h"
#include "csHelp.h"

// From geolib:
#include "csVector.h"
#include "csHeaderInfo.h"
#include "csCompareVector.h"
#include "geolib_string_utils.h"
#include "csException.h"
#include "csStandardHeaders.h"
#include "csGeolibUtils.h"

#include <sys/timeb.h>
#include <ctime>
#include <sys/time.h>

#ifdef USE_MPI
 #include <mpi.h>
#endif

using namespace cseis_system;

//extern int cseis_help( std::string moduleName );
//extern int cseis_help_html();

namespace cseis_system {
  extern void cseis_read_spreadSheet( char const* filenameSheet, FILE* f_spreadSheet,
                                      cseis_geolib::csCompareVector<csUserConstant>* masterConstants, 
                                      cseis_geolib::csVector< cseis_geolib::csVector<std::string> >* spreadSheetConstantList );
  extern void cseis_read_globalConst( FILE* f_globalConst, cseis_geolib::csCompareVector<csUserConstant>* globalConstList );
  extern int cseis_create_flow( FILE* f_master_flow, FILE* f_flow, cseis_geolib::csVector<csUserConstant> const* masterConstants );
  extern std::string replaceUserConstants( char const* line, cseis_geolib::csVector<cseis_system::csUserConstant> const* list );
}
void check_all_modules_for_bugs();

/// Error output stream
int exitOnError( char const* text, ... );
FILE* gl_error_stream;


/**
 * Main method
 *
 * Command line submission tool for CSEIS jobs
 *
 * @author Bjorn Olofsson
 */
int main( int argc, char** argv ) {
  char const* filenameFlow = NULL;
  char* filenameSheet = NULL;
  char* filenameLog   = NULL;
  bool isLogStdout = false;
  char* filenameGlobalConst = NULL;
  char* dirLog = NULL;
  int lengthDir = 0;
//  string dirLog;
  FILE* f_flow = NULL;
  FILE* f_master_flow;
  
  csLogWriter* f_log = NULL;
  std::string moduleName;
  cseis_geolib::csVector<std::string> filenameList;

  bool isSpreadSheet  = false;
  bool isRunMaster    = false;
  bool isRunFlow      = true;
  bool isVerbose      = true;
  bool isRunExec      = true;
  bool isGlobalConst  = false;
  bool isUserConstant = false;
  bool isOutputFlow   = false;
  bool isDebug        = false;
  //  char* flowOutputDir = NULL;
  char* flowOutputName= NULL;
  int memoryPolicy    = csMemoryPoolManager::POLICY_SPEED;
  cseis_geolib::csCompareVector<csUserConstant> globalConstList;
  int mpiNumProc = 1;
  int mpiProcID  = 0;
  
  gl_error_stream = stderr;

  //--------------------------------------------------------------------------------
  // Start MPI
  //
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiProcID);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiNumProc);
  fprintf(stderr,"MPI proc ID: %d, num MPI procs: %d\n", mpiProcID, mpiNumProc);
#endif
  
  //--------------------------------------------------------------------------------
  // Command line options
  //
//  for( int i = 0; i < argc; i++ ) {
//    printf("%d: %s\n", i, argv[i]);
//  }
//  return(-1);

  if( argc == 1 ) {
    fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
    return(-1) ;
  }
  int iArg = 1;   // Start with 1 to skip name of program executable
  while( iArg < argc ) {
    if( argv[iArg][0] == '-' ) {
      if( strlen(argv[iArg]) < 2 ) {
        fprintf(stderr," Syntax error in command line: Unknown option: '%s'\n", argv[iArg]);
        fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
        return(-1);
      }
      char option = argv[iArg][1];
      if ( option == 'h' ) {
        if( !strcmp( argv[iArg], "-html" ) ) {
          csHelp help( stdout );
          help.fullHtmlHelp();
          return(-1);
        }
        fprintf( stderr, " SeaSeis job flow submission tool.\n");
        fprintf( stderr, " Usage:  %s -f <jobflow> [-o <joblog> | -d <joblog_dir>] [-h] [-m <name>] [-v] [-c] [-std] [-p {speed|memory} ] [-g <const_file>] [-s <spreadsheet>]\n", argv[0] );
        fprintf( stderr, " -f <flow1> <flow2> ... : File name(s) of job flow(s) to run\n");
        fprintf( stderr, " -o [<log>|stdout]      : File name of job log (defaulted to flowname.log if not specified)\n");
        fprintf( stderr, "                        : Use 'stdout' to redirect all log file output to standard output\n");
        fprintf( stderr, " -d <log_dir>           : Name of directory where job log shall be saved\n");
        fprintf( stderr, " -s <spreadsheet file>  : Name of file containing spread sheet for building of master flow\n");
        fprintf( stderr, " -g <global const. file>: Name of file containing global constants (%cdefine statements)\n", csUserConstant::LETTER_DEFINE);
        fprintf( stderr, " -D <name1>[=,:]<val1> <name2>[=,:]<val2> ... : Define user constant <name>, value <val>. No spaces allowed around equal sign.\n");
        fprintf( stderr, "                                        NOTE: The first user constant will be used to create the log file name.\n");
        //        fprintf( stderr, " -fd <flow_dir>         : Name of directory where job flow shall be saved (only works when -D option is specified). \n");
        fprintf( stderr, " -ff <flow_out>         : Name of file where job flow shall be saved (only works when -D option is specified). \n");
        fprintf( stderr, " -m <name>              : Print help for module <name>.\n");
        fprintf( stderr, "                        : If no module is found, print help of all modules starting with <name>.\n");
        fprintf( stderr, "                        : If <name> = empty: Print full help for all available modules.\n");
        fprintf( stderr, "                        : If <name> = .: Print short help for all available modules.\n");
        fprintf( stderr, " -run_master            : Immediately run all flows created from master flow(s)\n");
        fprintf( stderr, "                        : If not specified, the program will exit after the master flows have been created.\n");
        fprintf( stderr, " -h                     : Print this page\n");
        fprintf( stderr, " -html                  : Print full help for all modules & standard trace headers in HTML.\n");
        fprintf( stderr, " -v                     : Print out version info\n");
        fprintf( stderr, " -std                   : Dump all standard trace headers\n");
        fprintf( stderr, " -c                     : Check for link problems and consistency of all modules' params(=help) methods.\n");
        fprintf( stderr, " -p [speed | memory]    : Set memory policy: Optimised for speed or memory.\n");
        fprintf( stderr, " -no_run                : Do not run flow. This option is useful if an individual flow file is generated using option -ff\n");
        fprintf( stderr, " -init_only             : Run init phase only.\n");
        fprintf( stderr, " -no_verbose            : Do not output information messages.\n");
        fprintf( stderr, " -debug                 : Output extensive DEBUG information for trace flow preparation and execution.\n");
        return(-1);
      }
      else if ( option == 'v' ) {
        fprintf( stderr, "SeaSeis version %s, build date %s\n", CSEIS_VERSION, __DATE__ );
        fprintf( stderr, "Copyright (c) Colorado School of Mines, 2013\n" );
        fprintf( stderr, "This program comes with ABSOLUTELY NO WARRANTY;\n");
        fprintf( stderr, "This is free software. Redistribution and use in source and binary forms, with or without\n");
        fprintf( stderr, "modification, are permitted provided that the conditions given in the license agreement are met.\n");
        return(-1);
      }
      else if ( option == 'f' ) {
        if( strlen(argv[iArg]) > 2 ) {
          if( argv[iArg][2] == 'f' ) {
            isOutputFlow  = true;
            flowOutputName = argv[++iArg];
          }
          //          else if( argv[iArg][2] == 'd' ) {
            //   isOutputFlow  = true;
            // flowOutputDir = new char [lengthDir+1];
            // memcpy( flowOutputDir, argv[iArg], length );
            // flowOutputDir[length] = FORWARD_SLASH;  // Doesn't matter whether this is Unix or Windows system, '/' should work
            // flowOutputDir[length+1] = '\0';
            //            flowOutputDir = argv[++iArg];
          //  }
          else {
            fprintf(stderr,"Unknown command line option : %s\n", argv[iArg]);
            exit(-1);
          }
          iArg++;
        }
        else {
          while( ++iArg < argc && argv[iArg][0] != '-' ) {
//          std::string fnameFlow = std::string();    // ???
            std::string fnameFlow(argv[iArg]);     
            filenameList.insertEnd(fnameFlow);
          }
        }
      }
      else if ( option == 'D' ) {
        while( ++iArg < argc && argv[iArg][0] != '-' ) {
//        std::string text = std::string(argv[iArg]);        // ???
          std::string text(argv[iArg]);
          int pos = text.find_first_of('=');
          if( pos == (int)text.npos ) {
            pos = text.find_first_of(':');
          }
          if( pos >= (int)text.length() ) {
            fprintf(stderr,"Error in command line option -D %s...", text.c_str());
          }
          globalConstList.insertEnd( csUserConstant(text.substr(0,pos),text.substr(pos+1)) );
          isUserConstant = true;

          //                                      userConstNames.insertEnd();
          //      userConstValues.insertEnd(text.substr(pos+1));
        }
      }
      else if ( option == 'r' ) {
        if( !strcmp( argv[iArg], "-run_master" ) ) {
          isRunMaster = true;
          ++iArg;
        }
        else {
          fprintf(stderr,"Unknown option '%s'\n", argv[iArg]);
          return(-1);
        }
      }
      else if ( option == 'o' ) {
        iArg++;
        if( iArg == argc ) {
          fprintf(stderr,"Missing argument for option -%c\n", option);
          return(-1);
        }
        int len_name = strlen(argv[iArg]);
        int len_name_all = len_name + 1;
        if( !strcmp(argv[iArg],"stdout") ) {
          isLogStdout = true;
        }
        else {
          // In case of MPI jobs (more than one job): Add MPI Proc ID to log name
          if( mpiProcID > 0 ) len_name_all += 5;
          filenameLog = new char [len_name_all];
          memcpy( filenameLog, argv[iArg], len_name );
          if( mpiProcID > 0 ) sprintf( &filenameLog[len_name], "%05d", mpiProcID );
        }
        filenameLog[len_name_all-1] = '\0';
        ++iArg;
      }
      else if ( option == 's' ) {
        if( !strcmp( argv[iArg], "-std" ) ) {
          csHelp help( stdout );
          help.standardHeaderHelp();
          //          dump_all_standard_headers();
          return(-1);
        }
        else {
          iArg++;
          if( iArg == argc ) {
            return exitOnError("Missing argument for option -%c\n", option);
          }
          filenameSheet = new char [(strlen(argv[iArg])+1)];
          memcpy( filenameSheet, argv[iArg], strlen(argv[iArg]) );
          filenameSheet[strlen(argv[iArg])] = '\0';
          isSpreadSheet = true;
        }
        ++iArg;
      }
      else if ( option == 'g' ) {
        iArg++;
        if( iArg == argc ) {
          return exitOnError("Missing argument for option -%c\n", option);
        }
        filenameGlobalConst = (char*)malloc( (strlen(argv[iArg])+1)*sizeof(char) );
        memcpy( filenameGlobalConst, argv[iArg], strlen(argv[iArg]) );
        filenameGlobalConst[strlen(argv[iArg])] = '\0';
        isGlobalConst = true;
        ++iArg;
      }
      else if ( option == 'n' ) {
        if( !strcmp( argv[iArg], "-no_run" ) ) {
          isRunFlow = false;
        }
        else if( !strcmp( argv[iArg], "-no_verbose" ) ) {
          isVerbose = false;
        }
        else {
          fprintf(stderr,"Unknown option '%s'\n", argv[iArg]);
          return(-1);
        }
        ++iArg;
      }
      else if ( option == 'i' ) {
        if( !strcmp( argv[iArg], "-init_only" ) ) {
          isRunExec = false;
        }
        else {
          fprintf(stderr,"Unknown option '%s'\n", argv[iArg]);
          return(-1);
        }
        ++iArg;
      }
      else if ( option == 'm' ) {
        ++iArg;
        std::string versionString = "";
        if( iArg == argc ) {
          moduleName = ",";
        }
        else {
          moduleName = argv[iArg];
          if( ++iArg < argc ) {
            versionString = argv[iArg];
          }
        }
        csHelp help( stdout );
        fprintf(stderr,"'%s' '%s'\n", moduleName.c_str(), versionString.c_str() );
        help.moduleTextHelp( moduleName, versionString );
        //cseis_help( moduleName );
        return(-1);
      }
      else if ( option == 'p' ) {
        ++iArg;
        if( iArg == argc ) {
          return exitOnError("Missing argument for option -%c\n", option);
        }
        std::string text = argv[iArg];
        if( !text.compare("speed") ) {
          memoryPolicy = csMemoryPoolManager::POLICY_SPEED;
        }
        else if( !text.compare("memory") ) {
          memoryPolicy = csMemoryPoolManager::POLICY_MEMORY;
        }
        else {
          fprintf(stderr,"Unknown memory policy: '%s'. Valid options are: speed, memory\n", argv[iArg] );
          return(-1);
        }
        ++iArg;
      }
      else if ( option == 'd' ) {
        if( !strcmp( argv[iArg], "-debug" ) ) {
          isDebug = true;
        }
        else {
          ++iArg;
          if( iArg >= argc ) {
            fprintf(stderr," Syntax error in command line: Option '-d' has no argument\n");
            return(-1);
          }
          int length = strlen(argv[iArg]);
          lengthDir = length+1;
          dirLog = new char [lengthDir+1];
          memcpy( dirLog, argv[iArg], length );
          dirLog[length] = FORWARD_SLASH;  // Doesn't matter whether this is Unix or Windows system, '/' should work
          dirLog[length+1] = '\0';
        }
        ++iArg;
      }
      else if ( option == 'c' ) {
        check_all_modules_for_bugs();
        return(-1);
      }
      else {
        fprintf(stderr," Syntax error in command line: Unknown option: '%s'\n", argv[iArg]);
        fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
        return(-1) ;
      }
//      iArg++;
    }
    else {
      fprintf(stderr," Syntax error in command line:\n Expected command line option (starting with '-'), found '%s'\n", argv[iArg]);
      fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
      return(-1) ;
    }
  }
  //--------------------------------------------------------------------------------
  // Open flow/log files, start log
  //
  if( isOutputFlow && !isUserConstant ) {
    fprintf(stderr,"Output flow name/dir was specified, but no user constants were given (option -D) to support this option.\n");
    fprintf(stderr,"If input flow shall be saved as separate file, at least one user constant must be provided by option -D.\n");
    exit(-1);
  }
  if( filenameList.size() == 0 ) {
    fprintf(stderr," No job flow specified. Type '%s -h' for help.\n", argv[0]);
    return(-1);
  }
  for( int i = 0; i < filenameList.size(); i++ ) {
    int length = filenameList.at(i).length();
    if( length < 6 || filenameList.at(i).substr(length-4,4).compare("flow") ) {
      fprintf(stderr,"Input file not recognized as SeaSeis flow. File name needs '.flow' extension.\n");
      fprintf(stderr,"Specified file name: '%s'\n", filenameList.at(i).c_str());
      return(-1);
    }
  }
  if( filenameList.size() > 1 && filenameLog != NULL && !isLogStdout ) {
    fprintf(stderr,"Warning: When specifying more than one input files, the option -o cannot be used.\n");
    fprintf(stderr,"         The specified name of the output log file will be ignored.\n");
    delete filenameLog;
    filenameLog = NULL;
  }

//---------------------------------------------------------------
// For master jobs, create all individual flow files which have all user constants filled in
// from the spreadsheet
//
  if( isSpreadSheet ) {
    FILE* f_spreadSheet = NULL;
    if( (f_spreadSheet = fopen(filenameSheet, "r")) == (FILE *) NULL ) {
      return exitOnError("Could not open spread sheet '%s'\n", filenameSheet );
    }
    cseis_geolib::csCompareVector<csUserConstant>* masterConstants = new cseis_geolib::csCompareVector<csUserConstant>;
    cseis_geolib::csVector< cseis_geolib::csVector<std::string> >*  spreadSheetConstantList = new cseis_geolib::csVector< cseis_geolib::csVector<std::string> >;
    // Read entire spreadsheet, store information in the given fields
    cseis_read_spreadSheet( filenameSheet, f_spreadSheet, masterConstants, spreadSheetConstantList );
    fclose( f_spreadSheet );
    // Create individual flows
    cseis_geolib::csVector<std::string> individualFlows;
    std::string keyName = masterConstants->at(0).pureName();
    for( int i = 0; i < filenameList.size(); i++ ) {
      std::string masterName = filenameList.at(i);
      if( (f_master_flow = fopen(masterName.c_str(), "r")) == (FILE *) NULL ) {
        return exitOnError("Could not open master flow '%s'\n", masterName.c_str() );
      }
      for( int is = 0; is < spreadSheetConstantList->size(); is++ ) {
        if( is > 0 ) rewind( f_master_flow );
        cseis_geolib::csVector<std::string> const* constantListPtr = &(spreadSheetConstantList->at(is));
        std::string flowName = masterName.substr( 0, masterName.length()-5 ) + "_" + keyName + constantListPtr->at(0) + ".flow";
        individualFlows.insertEnd( flowName );
//        fprintf(stderr,"Individual flow file '%s'\n", flowName.c_str());
        // Set replacement values for user constants found in spreadsheet
        for( int i = 0; i < masterConstants->size(); i++ ) {
          masterConstants->at(i).setValue( constantListPtr->at(i) );
        }
        if( (f_flow = fopen(flowName.c_str(), "w")) == (FILE *) NULL ) {
          return exitOnError( "Could not open spawned job flow '%s'\n", flowName.c_str() );
        }
        cseis_create_flow( f_master_flow, f_flow, masterConstants );
        fclose( f_flow );
      }
      fclose( f_master_flow );
    }
    filenameList.clear();
    filenameList = individualFlows;
    delete masterConstants;
    delete spreadSheetConstantList;
  }
  if( isSpreadSheet && !isRunMaster ) {
    fprintf(stderr,"Successfully created all individual master flows. To run, submit individual flows.\n");
    exit(0);
  }
//---------------------------------------------------------------------------------
// Read in global define statements
// Globally defined user constants will be replaced in all flows that are run, and take precedence over
// Locally defined constants with the same name
//
  if( isGlobalConst ) {
    FILE* f_globalConst = NULL;
    if( (f_globalConst = fopen(filenameGlobalConst, "r")) == (FILE *) NULL ) {
      return exitOnError("Could not open global defines file '%s'\n", filenameGlobalConst );
    }
    cseis_read_globalConst( f_globalConst, &globalConstList );
    fclose( f_globalConst );
  }

//---------------------------------------------------------------------------------
// Run flow(s)
//
  int returnFlag = 0;

  for( int i = 0; i < filenameList.size(); i++ ) {
    filenameFlow = filenameList.at(i).c_str();
    if ((f_flow = fopen(filenameFlow, "r")) == (FILE *) NULL) {
      fprintf(stderr,"Could not open file '%s'\n", filenameFlow );
      fprintf(stderr,"Flow cannot not be submitted...\n\n");
      continue;  // Try to submit next flow instead...
    }
    int lengthFlow       = filenameList.at(i).length();
    int lengthLog        = lengthFlow;
    int lengthLog_mpiAdd = ( mpiProcID == 0 ) ? 0 : 4; // In case of MPI jobs (more than one job): Add MPI Proc ID to log name
    if( filenameLog == NULL ) {
      if( dirLog == NULL ) {
        filenameLog = new char [lengthLog+lengthLog_mpiAdd];
        memcpy( filenameLog, filenameFlow, lengthFlow );
      }
      else {  // Log pathname was specified. Do not copy flow pathname (if there is any)
        int counter = lengthFlow-1;
        while( counter >= 0 && filenameFlow[counter] != FORWARD_SLASH && filenameFlow[counter] != BACK_SLASH ) counter--;
        lengthLog -= counter+1;
        if( isUserConstant ) {
          csUserConstant us = globalConstList.at(0);
          int lengthValue = us.value.length();
          filenameLog = new char [(lengthLog+lengthDir+lengthValue+1+lengthLog_mpiAdd)*sizeof(char)];  // +1 to add underscore '_'
          memcpy( filenameLog, dirLog, lengthDir );
          memcpy( filenameLog+lengthDir, filenameFlow+counter+1, lengthLog );
          filenameLog[lengthDir+lengthLog-5] = '_';
          memcpy( &filenameLog[lengthDir+lengthLog-4], us.value.c_str(), lengthValue );
          lengthLog = lengthLog + lengthDir + lengthValue + 1;  // +1: Underscore
          filenameLog[lengthLog-5] = '.';
        }
        else {
          filenameLog = new char [(lengthLog+lengthDir+lengthLog_mpiAdd)*sizeof(char)];
          memcpy( filenameLog, dirLog, lengthDir );
          memcpy( filenameLog+lengthDir, filenameFlow+counter+1, lengthLog );
          lengthLog += lengthDir;
        }
      }
      filenameLog[lengthLog-4] = 'l';
      filenameLog[lengthLog-3] = 'o';
      filenameLog[lengthLog-2] = 'g';
    } // END: if filenameLog == NULL
    else if( dirLog != NULL ) {  // log name was specified, and log dir was specified
      if( isLogStdout ) {
        fprintf(stderr,"When specifying option -d stdout (standard output), no log directory should be given using option -d <logDir>\n");
        return(-1);
      }
      lengthLog     = strlen(filenameLog);
      int counter   = lengthLog-1;
      int lengthDir = strlen(dirLog);
      while( counter >= 0 && filenameLog[counter] != FORWARD_SLASH && filenameLog[counter] != BACK_SLASH ) counter--;
      if( counter >= 0 ) {
        fprintf(stderr,"When specifying both options -d <logDir> and -o <logName>, the log name\n");
        fprintf(stderr,"must be given without any preceding path.\n");
        return(-1);
      }
      char* tmp = new char[(lengthLog+lengthDir+lengthLog_mpiAdd+1)];
      memcpy( tmp, dirLog, lengthDir );
      memcpy( tmp+lengthDir, filenameLog, lengthLog );
      delete [] filenameLog;
      filenameLog = tmp;
      lengthLog += lengthDir + 1;
    }
    if( mpiProcID != 0 ) {
      sprintf( &filenameLog[lengthLog-1], "%04d", mpiProcID );
    }
    filenameLog[lengthLog+lengthLog_mpiAdd-1] = '\0';

    if( isVerbose ) {
      fprintf(stderr,"Job flow:   %s\n", filenameFlow);
      fprintf(stderr,"Job log :   %s\n", isLogStdout ? "Redirected to standard output" : filenameLog );
    }
    try {
      if( isLogStdout ) {
        f_log = new csLogWriter();
      }
      else {
        f_log = new csLogWriter( filenameLog );
      }
    }
    catch( cseis_geolib::csException& exc ) {
      fprintf(stderr,"Could not open output log file '%s'\n", filenameLog );
      return -1;
    }

    //--------------------------------------------------------------------------------
    // Main processing block: RunManager
    //
    try {
      csRunManager runManager( f_log, memoryPolicy, mpiProcID, mpiNumProc, isDebug );
      if( isOutputFlow ) {
        FILE* f_flow_in;
        FILE* f_flow_out;
        if( (f_flow_in = fopen(filenameFlow, "r")) == NULL ) {
          fprintf(stderr,"Could not open file '%s'\n", filenameFlow );
          fprintf(stderr,"Flow cannot not be submitted...\n\n");
          continue;  // Try to submit next flow instead...
        }
        if( (f_flow_out = fopen(flowOutputName, "w")) == (FILE *) NULL ) {
          return exitOnError("Could not open output flow file '%s'\n", flowOutputName );
        }
        cseis_create_flow( f_flow_in, f_flow_out, &globalConstList );
        fclose(f_flow_out);
        fclose(f_flow_in);
      }
      if( isRunFlow ) {
        returnFlag = runManager.runInitPhase( filenameFlow, f_flow, &globalConstList );
        fclose(f_flow);
        if( returnFlag == 0 && isRunExec ) {
          returnFlag = runManager.runExecStartPhase();    
        }
        if( returnFlag != 0 ) {
          fprintf(stderr,"Exec START phase process returned error code: %d\n", returnFlag);
        }
        else {
          if( mpiProcID == 0 ) {
            returnFlag = runManager.runExecPhase();
          }
          else {
            returnFlag = runManager.mpi_runExecPhase();
          }
          if( returnFlag != 0 ) {
            fprintf(stderr,"Exec phase process returned error code: %d\n", returnFlag);
          }
          else if( isVerbose ) {
            fprintf(stderr,"Successful completion...\n");
          }
        }
      }
    }
    catch( cseis_geolib::csException& exc ) {
      fprintf(stderr,"Program was terminated. Message: %s\n", exc.getMessage());
      returnFlag = 99;
    }

    // Clean up...
    if( f_log ) {
      delete f_log;
      f_log = NULL;
    }

    if( filenameLog ) {
      delete [] filenameLog;
      filenameLog = NULL;
    }

  }  // END for all input files

  if( dirLog ) {
    delete [] dirLog;
    dirLog = NULL;
  }

  return returnFlag;
} // END main()

int exitOnError( char const* text, ... ) {
  va_list argList;
  va_start( argList, text );
  fprintf( gl_error_stream, "FATAL ERROR:\n" );
  vfprintf( gl_error_stream, text, argList );
  fprintf( gl_error_stream, "\n" );
  return(-1);
}

void check_all_modules_for_bugs() {
  /*
  FILE* fin = fopen( filename, "r" );
  if( fin == NULL ) return;

  std::ifstream fin();
  fin.open( filename, ios::in );
  if( fin.fail() ) {
    fprintf(stderr,"Error occurred when opening file containing list of available modules: %s\n", filename);
    return;
  }
  char buffer[1024];
  csVector<string> moduleNameList();

  while( !fin.getline( buffer, 1024 ).eof() ) {
    int len = strlen(buffer);
    if( len == 0 || buffer[0] == '#' ) continue;
    string line(buffer);
    moduleNameList.insertEnd( line );
  }
  fin.close();
  */

  int numModules = csMethodRetriever::getNumStandardModules();
  std::string const* names = csMethodRetriever::getStandardModuleNames();

  csParamDef pdef;

  fprintf(stdout,"\n#...checking all %d modules for link problems and bugs in params(=help) method\n\n", numModules);
  if( numModules == 0 ) {
    return;
  }

  try {
    for( int i = 0; i < numModules; i++ ) {
      pdef.clear();
      fprintf( stdout, "Checking module '%s'...\n", names[i].c_str() );
      try {
        csMethodRetriever::getParamMethod( names[i] )( &pdef );
      }
      catch( cseis_geolib::csException& exc ) {
        fprintf(stdout,"%s", exc.getMessage() );
        continue;
      }
      int numParams = pdef.numParameters();
      for( int ip = 0; ip < numParams; ip++ ) {
        int numValues = pdef.numValues(ip);
        if( numValues == 0 ) {
          fprintf( stdout, " Module '%s', Param '%s': No values found for this parameter.\n",
                   pdef.module()->name(), pdef.param(ip)->name() );
        }
        for( int iv = 0; iv < numValues; iv++ ) {
          int numOptions = pdef.numOptions(ip,iv);
          csParamDescription const* valueIdent = pdef.value(ip,iv);
          if( valueIdent->type() == VALTYPE_OPTION ) {
            if( numOptions == 0 ) {
              fprintf( stdout, " Module '%s', Param '%s', Value '%s': Specified as VALTYPE_OPTION but no options found\n",
                       pdef.module()->name(), pdef.param(ip)->name(), valueIdent->desc() );
            }
          }
          else {
            if( numOptions != 0 ) {
              fprintf( stdout, " Module '%s', Param '%s', Value '%s': Not specified as VALTYPE_OPTION but found %d options\n",
                       pdef.module()->name(), pdef.param(ip)->name(), valueIdent->desc(), numOptions );
            }
          }
        }

      }
    }
//    throw csException( message );
  }
  catch( cseis_geolib::csException& exc ) {
    fprintf(stdout,"%s", exc.getMessage() );
  }

}

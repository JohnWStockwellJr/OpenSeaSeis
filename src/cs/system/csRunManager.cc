/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csRunManager.h"
#include "csModule.h"
#include "csExecPhaseDef.h"
#include "csUserConstant.h"
#include "csUserParam.h"
#include "csParamDef.h"
#include "csParamManager.h"
#include "csModuleIndexStack.h"
#include "csTraceHeaderDef.h"
#include "csSuperHeader.h"
#include "csMemoryPoolManager.h"
#include "csLogWriter.h"
#include "csTimer.h"
#include "csMethodRetriever.h"

#include <stdarg.h>
#include <ctime>
#include <cstring>

// cseis : Geolib
#include "csException.h"
#include "csVector.h"
#include "csTimer.h"
#include "csTable.h"
#include "csCompareVector.h"
#include "geolib_string_utils.h"
#include "csEquationSolver.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace cseis_system;

//#define return_data_tag 2002

namespace cseis_system {
  extern std::string replaceUserConstants( char const* line, cseis_geolib::csVector<cseis_system::csUserConstant> const* list );
}

csRunManager::csRunManager( csLogWriter* log, int memoryPolicy, int mpiProcID, int mpiNumProc, bool isDebug ) {
  myIsDebug = isDebug;
  myLog     = log;
  myModules = NULL;
  myNumModules = 0;
  myNextModuleID = NULL;
  myPrevModuleID = NULL;
  myMemoryPoolManager = new csMemoryPoolManager( memoryPolicy );
  myTimerCPU = new cseis_geolib::csTimer();
  myTables = NULL;
  myNumTables = 0;
  myMPIProcID  = mpiProcID;
  myMPINumProc = mpiNumProc;
}
csRunManager::~csRunManager() {
  if( myTables != NULL ) {
    for( int i = 0; i < myNumTables; i++ ) {
      delete myTables[i];
    }
    delete [] myTables;
    myTables = NULL;
  }
  if( myModules != NULL ) {
    for( int imodule = 0; imodule < myNumModules; imodule++ ) {
      if( myModules[imodule] != NULL ) {
        delete myModules[imodule];
      }
    }
    delete [] myModules;
    myModules = NULL;
  }
  if( myNextModuleID != NULL ) {
    for( int imodule = 0; imodule < myNumModules; imodule++ ) {
      if( myNextModuleID[imodule] != NULL ) delete myNextModuleID[imodule];
    }
    delete [] myNextModuleID;
    myNextModuleID = NULL;
  }
  if( myPrevModuleID != NULL ) {
    for( int imodule = 0; imodule < myNumModules; imodule++ ) {
      if( myPrevModuleID[imodule] != NULL ) delete myPrevModuleID[imodule];
    }
    delete [] myPrevModuleID;
    myPrevModuleID = NULL;
  }
  if( myMemoryPoolManager != NULL ) {
    delete myMemoryPoolManager;
    myMemoryPoolManager = NULL;
  }
  if( myTimerCPU != NULL ) {
    delete myTimerCPU;
    myTimerCPU = NULL;
  }
}
//*********************************************************************************
// Init phase
//
int csRunManager::runInitPhase( char const* filenameFlow, FILE* f_flow, cseis_geolib::csCompareVector<cseis_system::csUserConstant>* globalConstList ) {
  myTimerCPU->start();
  if( myModules != NULL ) {
    for( int imodule = 0; imodule < myNumModules; imodule++ ) {
      if( myModules[imodule] != NULL ) {
        delete myModules[imodule];
      }
    }
    delete [] myModules;
    myModules = NULL;
    myNumModules = 0;
  }

  char line[MAX_LINE_LENGTH];
  cseis_geolib::csVector<std::string> tokenList;

  time_t timer = time(NULL);
  myLog->line( "================================================================================\n" );
  myLog->line( " SeaSeis Version %s   (built  %s %s)\n"\
               " Date: %s\n",
               CSEIS_VERSION, __DATE__, __TIME__, asctime(localtime(&timer)) );
  myLog->line( "================================================================================\n" );

  //********************************************************************************
  //
  // Read user constants defined in '&define' statements
  // Count modules
  //

  myLog->line( "\nPre-parse input file '%s'...\n", filenameFlow );

  //--------------------------------------------------------------------------------
  int counterModules = 0;  // Counts modules in input flow
  int counterLines   = 0;  // Counts lines in input flow
  int counterDefines = 0;  // Counts number of define statements
  
  cseis_geolib::csCompareVector<csUserConstant>* userConstantList = new cseis_geolib::csCompareVector<csUserConstant>(5);
  cseis_geolib::csCompareVector<std::string> tableNameList;
  cseis_geolib::csVector<cseis_geolib::csTable*> tableList;

  char* pos;
  bool isCarriageReturnWarning = false;

  while( fgets( line, MAX_LINE_LENGTH, f_flow ) != NULL ) {
    if( (pos = strrchr(line,CARRIAGE_RETURN)) != NULL ) {
      pos[0] = ' '; 
      if( (pos = strrchr(line,CARRIAGE_RETURN)) != NULL ) {
        if( !isCarriageReturnWarning ) {
          myLog->warning("Found more than one carriage return in input flow. This may indicate a missing newline.\nPlease check your input flow for hidden characters, especially if the flow fails\n.");
          isCarriageReturnWarning = true;
        }
        pos[0] = ' '; 
      }
    }
    counterLines++;
    char c;
    if( cseis_geolib::firstNonBlankChar( line, c ) ) {
      if( c == LETTER_MODULE ) {
        counterModules++;
      }
      else if( c == csUserConstant::LETTER_DEFINE ) {
        tokenList.clear();
        tokenize( line, tokenList );

        if( !tokenList.at(0).compare(csUserConstant::defineWord()) ) {
          counterDefines += 1;
          if( tokenList.size() != 3 ) {
            myLog->error( NULL, 0, "Syntax error in input flow, line %d: '%s'. &define statement requires two arguments: &define <name> <value>.",
                          counterLines, line );
          }
          csUserConstant uc(tokenList.at(1), tokenList.at(2));
          if( globalConstList != NULL && globalConstList->contains(uc) ) {
            myLog->warning("Line %d: '%s': Locally defined user constant '%s = %s' is masked by global definition.",
                           counterLines, line, uc.pureName().c_str(), uc.value.c_str() );
            // BUGFIX 090329: Wrong parameter list in warning message printout (first two were missing)
          }
          else if( userConstantList->contains( uc ) ) {
            myLog->error(NULL,0,"Syntax error in input flow, line %d: %s Duplicate definition of user constant '%s'.",
                         counterLines, line, uc.pureName().c_str() );
          }
          else {
            userConstantList->insertEnd( uc );
          }
        }
        // Read in table
        else if( !tokenList.at(0).compare(csUserConstant::tableWord()) ) {
          if( tokenList.size() < 3 ) {
            myLog->error( NULL, 0, "Syntax error in input flow, line %d: '%s'. &table statement requires two arguments: &table <tablename> <filename> {methodInterpolation}.",
                          counterLines, line );
          }
          cseis_geolib::csTable* table = new cseis_geolib::csTable();
          std::string tableName = tokenList.at(1);
          int methodInterpolation = 0;
          if( tokenList.size() > 3 ) {
            methodInterpolation = atoi( tokenList.at(4).c_str() );
          }
          try {
            table->initialize( tableName, tokenList.at(2), methodInterpolation );
          }
          catch( cseis_geolib::csException& e1 ) {
            myLog->error("Error occurred when reading in table '%s'. Error message: %s", tableName.c_str(), e1.getMessage() );
          }
          if( tableNameList.contains( tableName ) ) {
            myLog->error(NULL,0,"Syntax error in input flow, line %d: '%s'. Duplicate table '%s'.", counterLines, line, tableName.c_str() );
          }
          else {
            tableNameList.insertEnd( tableName );
            tableList.insertEnd( table );
          }

        }
      }
    }
  }
  // Copy all local and global constants into one single list
  if( globalConstList != NULL ) {
    for( int i = 0; i < globalConstList->size(); i++ ) {
      userConstantList->insertEnd( globalConstList->at(i) );
    }
  }
  myNumModules = counterModules;
  if( myNumModules == 0 ) {
    throw( cseis_geolib::csException("No modules found in flow file") );
  }

  myModules = new csModule*[myNumModules];
  myNextModuleID = new cseis_geolib::csVector<int>*[myNumModules];
  myPrevModuleID = new cseis_geolib::csVector<int>*[myNumModules];
  for( int imodule = 0; imodule < myNumModules; imodule++ ) {
    myModules[imodule] = NULL;
    myNextModuleID[imodule] = NULL;
    myPrevModuleID[imodule] = NULL;
  }

  myNumTables = tableList.size();
  if( myNumTables > 0 ) {
    myTables    = new cseis_geolib::csTable const*[myNumTables];
    for( int i = 0; i < myNumTables; i++ ) {
      myTables[i] = tableList.at(i);
      //    printf("Set table %d %s \n", i, myTables[i]->tableName().c_str() );
    }
  }

  csModule** modules = myModules;
  cseis_geolib::csVector<csUserParam*>** userParamList = new cseis_geolib::csVector<csUserParam*>*[myNumModules];

  myLog->line( "Pre-parser summary:\n");
  myLog->line( "  Parsed %3d lines,", counterLines );
  myLog->line( "   found %3d 'define' statement(s),",counterDefines );
  myLog->line( "         %3d table(s), and", myNumTables );
  myLog->line( "         %3d module(s).\n", myNumModules );

  counterLines   = 0;

  //********************************************************************************
  //
  // Parse job flow and...
  //
  //  - Create all modules
  //  - Store user parameters
  //  - Print job flow to log
  //

  myLog->line( "Input flow '%s':\n", filenameFlow );

  csUserParam* userParam = NULL;
  int moduleIndexCurr = -1;

  rewind( f_flow );
  while( fgets( line, MAX_LINE_LENGTH, f_flow ) != NULL ) {
    tokenList.clear();
    counterLines++;
    int counterChars = strlen(line)-1;
    while( counterChars >= 0 && (line[counterChars] == NEW_LINE || line[counterChars] == CARRIAGE_RETURN) ) {
      line[counterChars] = '\0';
      counterChars -= 1;
    }
    char c;
    if( !cseis_geolib::firstNonBlankChar( line, c ) || c == cseis_system::LETTER_COMMENT ) {   // This line is empty or a comment -> Skip to next line
      myLog->line( "%s", line );
      continue;
    }
    std::string tmpStr = replaceUserConstants( line, userConstantList );
    myLog->line( "%s", tmpStr.c_str() );  // Print input line to log file
    tokenize( tmpStr.c_str(), tokenList );
    if( tokenList.size() == 0 ) continue;   // Comment or empty line (this case should have been caught already)
    // DEFINE & TABLE statements have already been processed -> Skip to next line
    if( c == csUserConstant::LETTER_DEFINE && (!tokenList.at(0).compare(csUserConstant::defineWord()) || !tokenList.at(0).compare(csUserConstant::tableWord())) ) continue;
    if( tokenList.at(0)[0] == LETTER_MODULE ) {  // This token is a module because it starts with the 'module letter'
      moduleIndexCurr += 1;
      if( moduleIndexCurr == myNumModules ) {
        myLog->error( "Program bug, or non-ASCII input flow file? Inconsistent number of modules found.\nMaybe there are hidden non-ASCII characters in the input file?");
      }
      std::string const moduleName = tokenList.at(0).substr(1,tokenList.at(0).length()-1);   // substr() in order to skip leading 'module letter'

      //--------------- CREATE MODULE ---------------------
      modules[moduleIndexCurr] = new csModule( moduleName, moduleIndexCurr, myMemoryPoolManager, myMPIProcID, myMPINumProc );
      myNextModuleID[moduleIndexCurr] = new cseis_geolib::csVector<int>(1);
      myPrevModuleID[moduleIndexCurr] = new cseis_geolib::csVector<int>(1);
      userParamList[moduleIndexCurr]  = new cseis_geolib::csVector<csUserParam*>(1);
    }
    else if( moduleIndexCurr < 0 ) {
      myLog->error( NULL, 0, "Syntax error in input flow, line %d: '%s'\nName of first module missing.\n", counterLines, line );
    }
    else {  // This line contains user parameters. Save to user parameter list of the current module
      userParam = new csUserParam( tokenList );
      userParamList[moduleIndexCurr]->insertEnd( userParam );
    }
  }

  //  fclose( f_flow );  // DO NOT CLOSE FILE

  myLog->line( "\nEnd of flow." );
  myLog->line( "================================================================================\n" );

  delete userConstantList;
  userConstantList = NULL;

  //********************************************************************************
  // Omit all modules above the last INPUT module, i.e. run the flow from the last INPUT module onwards.
  //
  for( int imodule = myNumModules-1; imodule > 0; imodule-- ) {
    if( !strcmp( modules[imodule]->getName(), "INPUT" ) ) {
      if( imodule == 0 ) {
        break;
      }
      if( imodule == myNumModules-1 ) {
        myLog->warning("Last module in flow is INPUT module. Flow will have no effect.\n");
      }
      if( myNumModules > 1 ) {
        myLog->warning("Module #%d in flow is INPUT module. Flow will run from this module onwards, previous modules will not be processed.\n", imodule+1 );
        // Remove earlier modules from flow. Reset various arrays.
        for( int imodule2 = 0; imodule2 < imodule; imodule2++ ) {
          delete modules[imodule2];
          delete myNextModuleID[imodule2];
          delete myPrevModuleID[imodule2];
          for( int iparam = 0; iparam < userParamList[imodule2]->size(); iparam++ ) {
            delete userParamList[imodule2]->at(iparam);
          }
          delete userParamList[imodule2];
        }
        // Move objects
        for( int imodule2 = imodule; imodule2 < myNumModules; imodule2++ ) {
          int index = imodule2-imodule;
          modules[index]        = modules[imodule2];
          myNextModuleID[index] = myNextModuleID[imodule2];
          myPrevModuleID[index] = myPrevModuleID[imodule2];
          userParamList[index]  = userParamList[imodule2];
        }
        myNumModules = myNumModules - imodule;
        break;
      }
    }
  }

  //********************************************************************************
  // Check module names, user parameters, options etc..
  // Convert all user inputs to lower case. Exception: String variables are kept as they are
  //
  //---------------- SET NEXT MODULE(S) & IF/SPLIT BLOCKS --------------
  //
  //  - Set flow succession, special attention goes to IF and SPLIT blocks
  //
  cseis_geolib::csStack<modInfoStruct> modIFStack;     // Stack that is used to determine flow succession in IF blocks
  cseis_geolib::csStack<modInfoStruct> modSPLITStack;  // Stack that is used to determine flow succession in SPLIT blocks
  modInfoStruct item;
  int moduleIndexNext;

  // !CHANGE! myNextModuleID has not been inserted a single ID when there is only one module in the flow. Possibly this also happens for last module...?

  for( int imodule = 0; imodule < myNumModules; imodule++ ) {
    csModule* module = modules[imodule];
    if( imodule < myNumModules-1 ) {
      myNextModuleID[imodule]->insertEnd( imodule+1 );
    }
    int type = module->getType();
    if( imodule == 0 ) continue;  // First module in flow. This MUST be an INPUT module. Skip to next module

    int modIndex_else_elseif;
    switch( type ) {
    case MODTYPE_UNKNOWN:  // Nothing to be done...
      break;
    case MODTYPE_IF:
      modIFStack.push( modInfoStruct( imodule, type ) );
      break;
    case MODTYPE_SPLIT:
      modSPLITStack.push( modInfoStruct( imodule, type ) );
      break;
    case MODTYPE_ELSEIF:
    case MODTYPE_ELSE:
      if( modIFStack.isEmpty() ) myLog->error( module->getName(), imodule, "Inconsistent IF block(s): No matching IF module found.");
      modIFStack.push( modInfoStruct( imodule, type ) );
      break;
    case MODTYPE_ENDIF:
      if( modIFStack.isEmpty() ) myLog->error( module->getName(), imodule, "Inconsistent IF block(s): No matching IF module found.");
      item = modIFStack.pop();
      modIndex_else_elseif = item.modIndex;
      moduleIndexNext = imodule;
      while( !modIFStack.isEmpty() && item.modType != MODTYPE_IF ) {  // Step backwards through IF block until starting IF module is found
        if( item.modType == MODTYPE_ELSEIF ) {
          myNextModuleID[item.modIndex]->insertEnd( moduleIndexNext ); // Set next module index to next ELSEIF/ELSE/ENDIF module
        }
        myNextModuleID[modIndex_else_elseif-1]->set( imodule, 0 );  // For module inside IF block that is preceding ELSE or ELSEIF, set next module to ENDIF module
        item = modIFStack.pop();
        moduleIndexNext = modIndex_else_elseif;
        modIndex_else_elseif = item.modIndex;
      }
      if( item.modType != MODTYPE_IF ) myLog->error( module->getName(), imodule, "Inconsistent IF block(s): No matching IF module found.");
      myNextModuleID[item.modIndex]->insertEnd( moduleIndexNext );  // For IF module, set next module index to next ELSEIF/ELSE/ENDIF module
      break;
    case MODTYPE_ENDSPLIT:
      if( modSPLITStack.isEmpty() ) myLog->error( module->getName(), imodule, "Inconsistent SPLIT block(s): No matching SPLIT module found.");
      item = modSPLITStack.pop();
      if( item.modIndex == imodule ) {
        myLog->error("Empty SPLIT-ENDSPLIT block. Remove or comment out $SPLIT/$ENDSPLITs");
      }
      // If ENDSPLIT is last module, set NextModule pointer of SPLIT module to ENDSPLIT. Otherwise, set it to next module after ENDSPLIT (imodule+1)
      //      myNextModuleID[item.modIndex]->insertEnd( std::min(imodule+1,myNumModules-1) );
      // TEMP!!!
      myNextModuleID[item.modIndex]->insertEnd( imodule+1 );
      break;
    default:
      break;
    }
  }

  // Last check of IF and SPLIT blocks
  if( !modIFStack.isEmpty() ) {
    item = modIFStack.pop();
    myLog->error( modules[item.modIndex]->getName(), item.modIndex, "Inconsistent IF block(s): No matching ENDIF found.");
  }
  if( !modSPLITStack.isEmpty() ) {
    item = modSPLITStack.pop();
    myLog->error( modules[item.modIndex]->getName(), item.modIndex, "Inconsistent SPLIT block(s): No matching ENDSPLIT found.");
  }
  modIFStack.dispose();
  modSPLITStack.dispose();

  //---------------- SET PREVIOUS MODULE(S) --------------
  // Each module has to know its previous module in the flow succession, in order to know which headers are defined and with what values.
  // Step through each module, retrieve its successive module(s), and add it to their list of previous modules (may be more than one for MERGE or ENDIF or...)
  //
  for( int imodule = 0; imodule < myNumModules; imodule++ ) {
    for( int iport = 0; iport < myNextModuleID[imodule]->size(); iport++ ) {
      int moduleIndexNext = myNextModuleID[imodule]->at( iport );
      if( moduleIndexNext > 0 && moduleIndexNext < myNumModules ) {
        myPrevModuleID[moduleIndexNext]->insertEnd( imodule );
      }
    }
  }
  // Check if all modules have a 'previous' modules: If not it cannot be reached
  for( int imodule = 1; imodule < myNumModules; imodule++ ) {
    if( myPrevModuleID[imodule]->size() == 0 ) {
      myLog->line("Warning: Module #%d '%s' cannot be reached.", imodule+1, modules[imodule]->getName() );
    }
  }

  //----------------------------
  if( myIsDebug ) {
    for( int imodule = 0; imodule < myNumModules; imodule++ ) {
      csModule* module = modules[imodule];
      if( myNextModuleID[imodule]->size() == 1 ) {
        myLog->write("Module %14s:  this(%2d),  next(%2d),   prev: ",
                     module->getName(), imodule+1, myNextModuleID[imodule]->at(0)+1 );
      }
      else if( myNextModuleID[imodule]->size() == 2 ) {  // Max 2 output ports
        myLog->write("\nModule %14s:  this(%2d),  next1(%2d)  next2(%d), prev: ",
                     module->getName(), imodule+1, myNextModuleID[imodule]->at(0)+1, myNextModuleID[imodule]->at(1)+1 );
      }
      else {
        myLog->write("Module %14s:  this(%2d):     prev: ",
                     module->getName(), imodule+1 );
      }
      for( int i = 0; i < myPrevModuleID[imodule]->size(); i++ ) {
        myLog->write(" (%2d)", myPrevModuleID[imodule]->at(i)+1 );
      }
      //for( int i = 0; i < module->getNumPrevModules(); i++ ) {
      //  myLog->write(" (%2d)", module->getPrevModuleIndex(i)+1 );
      //}
      myLog->line("");
    }
  }
  
  //***********************************************************************
  //
  // Init phase
  //
  //
  myLog->line( "Run init phase...\n" );

  cseis_geolib::csVector<csModule const*> prevModuleList;
  cseis_geolib::csVector<csModule const*> uniqModuleList;  // List of unique modules in current flow

  for( int imodule = 0; imodule < myNumModules; imodule++ ) {
    csModule* module = modules[imodule];
    //    if( !csMethodRetriever::isValidModule( module->getName() ) ) {
    //      myLog->line("Module '%s': No module with this name found. Name mis-spelled..?\n", module->getName() );
    //    throw( csException("Module not found") );
    //  }
    myLog->line("\nRun INIT phase for module  #%-5d %s...", imodule+1, module->getName() );
    if( myIsDebug ) fprintf(stdout,"Module %2d %s...\n", imodule, module->getName());
    csParamManager paramManager( userParamList[imodule], myLog );
    if( paramManager.exists( "debug" ) ) {
      module->setDebugFlag( true );
    }
    if( paramManager.exists( "version" ) ) {
      std::string versionString;
      paramManager.getString( "version", &versionString );
      int major = 0;
      int minor = 0;
      if( !parseVersionString( versionString, major, minor ) ) {
        myLog->error("Incorrect version string given: '%s'.  Version string must match the format  X.X (X=0-99)", versionString.c_str() );
      }
      myLog->line("User specified version: %s", versionString.c_str() );
      module->setVersion( major, minor );
    }
    try {
      // Set input ports (sets header definitions and other information from previous modules)
      int nInputPorts = myPrevModuleID[imodule]->size();
      prevModuleList.clear();
      for( int inPort = 0; inPort < nInputPorts; inPort++ ) {
        if( modules[myPrevModuleID[imodule]->at(inPort)]->getType() != MODTYPE_ENDSPLIT ) {  // Do not consider ENDSPLITs as 'previous' modules
          prevModuleList.insertEnd( modules[myPrevModuleID[imodule]->at(inPort)] );
        }
      }
      module->setInputPorts( &prevModuleList );

      csParamDef paramDef;
      paramDef.clear();
      // Check module name, parameters, number of values & options
      module->getParamDef( &paramDef );
      if( !csRunManager::checkParameters( module->getName(), &paramDef, userParamList[imodule], myLog ) ) {
        throw( cseis_geolib::csException("Error occurred while checking module parameters. See log file for details.") );
      }
      // MPI master process: Submit init phase, check if module supports MPI. Send MPI support flag toall pocesses. If MPI is supported, run init phase for all MPI processes
      if( myMPIProcID == 0 ) {
        // (a) First, run init phase for master process
        module->submitInitPhase( &paramManager, myLog, myTables, myNumTables );
#ifdef USE_MPI
        // START of MPI code
        bool isMPISupported = module->isMPISupported();
        int tag = 99;
        // (b) Send MPI supportflag to all other processes
        for( int mpiTargetProc = 1; mpiTargetProc < myMPINumProc; mpiTargetProc++ ) {
          MPI_Send( &isMPISupported, 1, MPI::BOOL, mpiTargetProc, tag, MPI_COMM_WORLD );
        }
        // If MPI is not supported, need to copy SuperHeader and TraceHeaderDef to all processes since init phase is not run in these.
        // Updated SuperHeader and TraceHeaderDef fields are needed in all MPI processes for subsequent modules which may support MPI
        if( myMPINumProc > 1 && !isMPISupported ) {
          int byteSize_shdr = module->getSuperHeader()->mpi_getByteSize();
          int byteSize_hdef = module->getHeaderDef()->mpi_getByteSize();
          int totalByteSize = byteSize_shdr + byteSize_hdef;
          int execType = module->getExecType();
          char* buffer = new char[totalByteSize];
          module->getSuperHeader()->mpi_compress( &buffer[0] );
          module->getHeaderDef()->mpi_compress( &buffer[byteSize_shdr] );
          for( int mpiTargetProc = 1; mpiTargetProc < myMPINumProc; mpiTargetProc++ ) {
            MPI_Send( &execType, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
            MPI_Send( &byteSize_shdr, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
            MPI_Send( &byteSize_hdef, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
            MPI_Send( buffer, totalByteSize, MPI::CHAR, mpiTargetProc, tag, MPI_COMM_WORLD );
          }
          delete [] buffer;
        }
        // END of MPI code
#endif
        if( module->getExecType() == EXEC_TYPE_INPUT && imodule != 0 ) {
          myLog->error( module->getName(), imodule, "Only the first module in the flow can be an INPUT module." );
        }
        if( imodule == 0 && module->getExecType() != EXEC_TYPE_INPUT ) {
          myLog->warning( "Module %s: First module in flow must be an INPUT module. If not, no exec phase will be run", module->getName() );
        }
      } // END if MPIProc == 0 : Main process, or single process
      else { // Other MPI processes (other than main)
#ifdef USE_MPI
        // START of MPI code
        bool isMPISupported;
        int tag = 99;
        // (c) Receive MPI support flag from master process
        MPI_Status stat;
        MPI_Recv( &isMPISupported, 1, MPI::BOOL, 0, tag, MPI_COMM_WORLD, &stat );
        // MPIDEBUG     if( myIsDebug ) fprintf(stderr,"Module #%2d, MPI proc %d. MPI supported: %s\n", imodule, myMPIProcID, isMPISupported ? "yes" : "no");
        // (d) If MPI is supported for this module, run init phase for all MPI processes
        if( isMPISupported ) {
          module->submitInitPhase( &paramManager, myLog, myTables, myNumTables );
        }
        // (e) If MPI is NOT supported, need to update superheader and headerdef from proc #0 for potential later MPI-supported modules
        else {
          int execType;
          int byteSize_shdr;
          int byteSize_hdef;
          MPI_Recv( &execType, 1, MPI::INT, 0, tag, MPI_COMM_WORLD, &stat );
          module->setExecType( execType );
          MPI_Recv( &byteSize_shdr, 1, MPI::INT, 0, tag, MPI_COMM_WORLD, &stat );
          MPI_Recv( &byteSize_hdef, 1, MPI::INT, 0, tag, MPI_COMM_WORLD, &stat );
          int totalByteSize = byteSize_shdr + byteSize_hdef;
          char* buffer = new char[totalByteSize];
          MPI_Recv( buffer, totalByteSize, MPI_BYTE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &stat );
          module->mpi_decompressInitPhaseObjects( &buffer[0], &buffer[byteSize_shdr] );
          delete [] buffer;
        }
        // END of MPI code
#endif
      } // END if/else other MPI processes
      //      module->getHeaderDef()->dump( myLog->getFile() );
      //      module->getSuperHeader()->dump( myLog->getFile() );
      //####################################################################################################################################
    }
    catch( cseis_geolib::csException& e ) {
      myLog->line("\n");
      myLog->line("Module #%d %s:", imodule+1, module->getName());
      myLog->line(" FATAL ERROR occurred in INIT phase.");
      myLog->line(" %s", e.getMessage() );
      myLog->line(" Flow terminated.\n");

      // Clean up...
      for( int imod2 = 0; imod2 <= imodule; imod2++ ) {
        module = modules[imod2];
        if( myMPIProcID == 0 || module->isMPISupported() ) {
          module->submitCleanupPhase( myLog, false );
        }
      }
      for( int imod2 = 0; imod2 < myNumModules; imod2++ ) {
        for( int i = 0; i < userParamList[imod2]->size(); i++ ) {
          delete userParamList[imod2]->at(i);
        }
        delete userParamList[imod2];
      }
      delete [] userParamList;
      throw( e );  // Rethrow to caller
    }
    catch(...) {
      if( myMPIProcID == 0 ) {
        myLog->line("\n");
        myLog->line("Module #%d %s:", imodule+1, module->getName());
        myLog->line(" FATAL ERROR occurred in INIT phase.");
        myLog->line(" Flow terminated.\n");
      }
      throw( cseis_geolib::csException("Program error. Flow terminated.") );  // Rethrow to caller
    }

    if( myMPIProcID == 0 ) {
      for( int i = 0; i < userParamList[imodule]->size(); i++ ) {
        if( paramManager.getNumValueCalls(i) == 0 ) {
          myLog->line( " Warning: Unused user parameter in module setup: '%s'", userParamList[imodule]->at(i)->getName().c_str() );
          int numLines = paramManager.getNumLines(userParamList[imodule]->at(i)->getName().c_str());
          if( numLines > 1 ) {
            myLog->error( modules[imodule]->getName(), imodule+1, "Duplicate user parameter '%s' (occurs %d times)",
                          userParamList[imodule]->at(i)->getName().c_str(), numLines );
          }
        }
      }
    }
  }
  prevModuleList.dispose();

  if( myIsDebug && myMPIProcID == 0 ) {
    fprintf(stdout,"--------------------------------------\n\n");
    for( int imodule = 0; imodule < myNumModules; imodule++ ) {
      fprintf(stdout,"   Module %d %s header dump:\n", imodule, modules[imodule]->getName() );
      modules[imodule]->getHeaderDef()->dump(stdout);
    }
  }
  //------------------------------------------------------
  // Check that all input ports have same sample interval and number of samples
  // Check if all modules have a 'previous' modules: If not it cannot be reached
  //
  if( myMPIProcID == 0 ) {
    for( int imodule = 1; imodule < myNumModules; imodule++ ) {
      if( myPrevModuleID[imodule]->size() == 0 ) continue;
      int firstModID = myPrevModuleID[imodule]->at(0);
      csSuperHeader const* shdrFirst = modules[firstModID]->getSuperHeader();
      for( int iport = 1; iport < myPrevModuleID[imodule]->size(); iport++ ) {
        if( modules[myPrevModuleID[imodule]->at(iport)]->getType() == MODTYPE_ENDSPLIT ) continue; // Do not consider ENDSPLITs as 'previous' modules
        csSuperHeader const* shdr = modules[myPrevModuleID[imodule]->at(iport)]->getSuperHeader();
        if( shdrFirst->numSamples != shdr->numSamples ) {
          myLog->error( modules[imodule]->getName(), imodule+1, "Unequal number of samples for input port %d", iport+1 );
        }
        else if( shdrFirst->sampleInt != shdr->sampleInt ) {
          myLog->error( modules[imodule]->getName(), imodule+1, "Unequal sample interval for input port %d", iport+1 );
        }
      }
    }
  }

  //***********************************************************************
  // Output module version in use
  //
  /*
   * Reproducable input flow:
   * Functionality is not implemented yet to use this, so it doesn't make sense to output it yet to the log file
   * The reproducable input flow shall enable running a flow by supplying the log file.
   * The reproducable input flow makes sure that it produces the exact same result as last time when the flow was run,
   * by specifying all parameters including the module version

   cseis_geolib::csVector<std::string> valueList;
   std::string flowText = "";
   char lineText[1024];
   //
   for( int imodule = 0; imodule < myNumModules; imodule++ ) {
   csModule* module = modules[imodule];
   sprintf(lineText,"\n%c%s\n", LETTER_MODULE, module->getName() );
   flowText.append(lineText);
   csParamManager paramManager( userParamList[imodule], myLog );
   if( !paramManager.exists( "version" ) ) {
   sprintf(lineText," %-14s  %-5s\n", "version", module->versionString().c_str() );
   flowText.append(lineText);
   }
   for( int ip = 0; ip < userParamList[imodule]->size(); ip++ ) {
   csUserParam* param = userParamList[imodule]->at(ip);
   param->getValues(&valueList);
   sprintf(lineText," %-14s  ",param->getName().c_str());
   flowText.append(lineText);
   for( int iv = 0; iv < valueList.size(); iv++ ) {
   // !CHANGE! Only put in quotation marks if value contains spaces:
   sprintf(lineText,"\"%s\" ", valueList.at(iv).c_str());
   flowText.append(lineText);
   }
   flowText.append("\n");
   }
   }

   myLog->line( "\n********************************************************************************" );
   myLog->line( "********************************************************************************" );
   myLog->line( "Reproducable input flow:" );
   myLog->line( "%s", flowText.c_str() );
   // !CHANGE! Output this flowText to SeaSeis output file (whenever one is generated in this flow)
   myLog->line( "\n********************************************************************************" );
   myLog->line( "********************************************************************************" );
  */

  for( int imodule = 0; imodule < myNumModules; imodule++ ) {
    for( int i = 0; i < userParamList[imodule]->size(); i++ ) {
      delete userParamList[imodule]->at(i);
    }
    delete userParamList[imodule];
  }
  delete [] userParamList;
  userParamList = NULL;

  //----------------------------------------------------------------------
  // Dump init phase info
  //
  if( myMPIProcID == 0 ) {
    myLog->line( "\n------------------------------------------------------------" );
    myLog->line( "Init phase summary\n" );
    myLog->line( "%3s  %-19s %11s %11s %9s", "#", "Module", "sampleInt", "numSamples", "#headers" );
    for(int iModule = 0; iModule < myNumModules; iModule++ ) {
      csSuperHeader const* h = modules[iModule]->getSuperHeader();
      myLog->line( "%3d  %-19s %11.5f %11d %9d", iModule+1, modules[iModule]->getName(), h->sampleInt, h->numSamples,
                   modules[iModule]->getHeaderDef()->numHeaders());
    }
    double timeCPUAll = myTimerCPU->getElapsedTime();
    myLog->line( "\nInit phase processing time:  %12.6f seconds\n\n", timeCPUAll );
  }

  return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------
//
int csRunManager::runExecStartPhase() {
  int returnFlag = 0;

  myLog->line( "\n================================================================================\n" );
  myLog->line( "Run 'start exec' phase...\n" );

  // Run start exec phase.
  for( int iModule = 0; iModule < myNumModules; iModule++ ) {
    //    if( myModules[iModule]->isMPISupported() ) {
    if( !myModules[iModule]->submitExecStartPhase( myLog ) ) {
      myLog->line("Error occurred during start exec phase of module #%d %s\n", iModule+1, myModules[iModule]->getName() );
      returnFlag = 22;
    }
    //    }
  }
  return returnFlag;
}

//***********************************************************************
//
// Exec phase processing loop - Other processes (not #0)
//
int csRunManager::mpi_runExecPhase() {
  int returnFlag = 0;
#ifdef USE_MPI
  // START of MPI code
  csModule** modules = myModules;
  bool atEOF = false;
  int tag = 99;
  int moduleIndex = -1;

  MPI_Status stat;
  try {
    do {
      // RECV: Wait for instructions from process #0
      MPI_Recv( &moduleIndex, 1, MPI::INT, 0, tag, MPI_COMM_WORLD, &stat );
      if( moduleIndex >= 0 ) {
        int outPort = 0;
        // bool success = modules[moduleIndex]->mpi_submitExecPhase( myLog, outPort );
        modules[moduleIndex]->mpi_submitExecPhase( myLog, outPort );
      }
      else { // if( mpiInstructionFlag == csRunManager::MPI_TERMINATE ) {
        atEOF = true;
      }
    } while( !atEOF );
  }
  catch( cseis_geolib::csException& e ) {
    if( moduleIndex >= myNumModules || moduleIndex < 0 ) {
      moduleIndex = 0;
    }
    myLog->line("\ncsRunManager::mpi_runExecPhase(): Exception caught while running exec phase, module #%2d %s. System message: \n%s",
                moduleIndex+1, myModules[moduleIndex]->getName(), e.getMessage());
    exit( -1 );
  }

  // Run cleanup phase. Cleaning up allocated memory in modules
  for( int iModule = 0; iModule < myNumModules; iModule++ ) {
    if( myModules[iModule]->isMPISupported() ) {
      myModules[iModule]->submitCleanupPhase( myLog );
    }
  }
  execPhaseCompletionTasks();
  // END of MPI code
#endif
  return returnFlag;
}

//*************************************************
// Exec phase processing loop - MPI process #0
//
int csRunManager::runExecPhase() {
  int returnFlag = 0;
  int indexLastModule  = myNumModules-1;
  cseis_geolib::csModuleIndexStack stack_moduleIndex( myNumModules+1 );

  myLog->line( "Run exec phase...\n" );
  if( myIsDebug ) fprintf(stdout,"--------------------------------------\n\n");

  myLog->flush();
  int iModule = 0;

  //MPIDEBUG   myLog->line("csRunManager: MPI %d: RunManager exec phase loop\n", myMPIProcID);
  
  try {
    bool isInputFinished = false;
    if( myModules[0]->getExecType() != EXEC_TYPE_INPUT ) isInputFinished = true;
    //---------------------------------------------------------------------
    // BIG LOOP FOR ALL INPUT TRACES
    //
    while( !isInputFinished ) {
      // STEP (1) Read in input trace. If none is read in, set isInputFinished = true. No trace will be passed on.
      int outPort = 0;
      iModule = myNumModules;

      if( myModules[0]->submitExecPhase( isInputFinished, myLog, outPort ) ) {  // Successful submission of exec phase
        if( myNextModuleID[0]->size() > 0 ) {  // Only move on traces if Input module is not the only (=last) module in flow
          iModule = myNextModuleID[0]->at(0);
          myModules[iModule]->moveTracesFrom( myModules[0], outPort );
        }
        else {  // Input module is only module in flow --> clean up traces
          myModules[0]->lastModuleTraceCleanup();
        }
      }
      else {
        isInputFinished = true;
        if( myNextModuleID[0]->size() > 0 ) { 
          iModule = myNextModuleID[0]->at(0);
        }
        else {
          iModule = myNumModules;  // ...for flows containing single module
        }
      }
      //MPIDEBUG myLog->line("csRunManager: MPI %d: Is Input finished: %d", myMPIProcID, isInputFinished);
      //      if( isInputFinished ) printf("Input is finished...\n");
      //---------------------------------------------------------------------
      // INNER LOOP FOR ALL MODULES EXCEPT INPUT MODULE
      //
      while( iModule < myNumModules ) {
        csModule* module = myModules[iModule];
        if( myIsDebug ) fprintf(stdout,"Module %2d %s... MPI process: %d\n", iModule, module->getName(), myMPIProcID);
        //MPIDEBUG myLog->line("MPI %d: Module %2d %s...\n", myMPIProcID, iModule, module->getName());
        outPort = 0;
        // Force module to run when...
        // a) last trace has been read in, and
        // b) previous modules have finished processing (stackModuleIndex.isEmpty()), and
        // c) current module has not finished processing yet  (this is checked inside method isReadyToSubmitExecPhase)
        bool forceToRun = isInputFinished && stack_moduleIndex.isEmpty();
        if( myIsDebug ) fprintf(stdout,"  Forced to run?  %d (module has %d stored traces)\n", forceToRun, module->tempNumTraces() );
        //MPIDEBUG myLog->line("MPI %d:  Forced to run?  %d (module has %d stored traces)\n", myMPIProcID, forceToRun, module->tempNumTraces() );
        // STEP (2) Check if module is ready for submission
        if( !module->isReadyToSubmitExecPhase( forceToRun ) ) {
          if( myIsDebug ) fprintf(stdout,"  ...is NOT ready to submit\n");
          if( !forceToRun ) {
            iModule = stack_moduleIndex.pop();  // returns myNumModules if empty
          }
          else {
            iModule += 1;  // If forced, step to next module in list
          }
        }
        // STEP (3) If module is ready for submission (or when forced), submit exec phase
        else {
#ifdef USE_MPI
          // START of MPI code
          static int tag = 99;
          if( module->isMPISupported() ) {
            for( int mpiTargetProc = 1; mpiTargetProc < myMPINumProc; mpiTargetProc++ ) {
              MPI_Send( &iModule, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
            }
          }
          // END of MPI code
#endif
          bool success = module->submitExecPhase(forceToRun,myLog,outPort);
          if( !success ) {
            //MPIDEBUG myLog->line("MPI %d:  STEP 3 ...did NOT run correctly\n", myMPIProcID );
            if( myIsDebug ) fprintf(stdout,"  STEP 3 ...did NOT run correctly\n");
            // Exec phase failed. This indicates that no trace was output
            if( module->finishedProcessing() ) {
              if( !forceToRun ) {
                iModule = stack_moduleIndex.pop();  // returns myNumModules if empty
              }
              else {
                iModule += 1;
              }
            }
            // else: nothing to be done. Continue processing this module, keep iModule
          }
          // STEP (4) Module has been processed, current module is last module in flow
          else if( iModule == indexLastModule ) {
            //MPIDEBUG      myLog->line("MPI %d:  STEP 4 ...is last module\n", myMPIProcID );
            if( myIsDebug ) fprintf(stdout,"  STEP 4 ...is last module\n");
            module->lastModuleTraceCleanup();
            if( module->finishedProcessing() ) {
              iModule = stack_moduleIndex.pop();  // returns myNumModules if empty
            } // else keep iModule, process next trace(s)
          }
          // STEP (5) Module has been processed...
          else {
            //MPIDEBUG myLog->line("MPI %d:  STEP 5 ...else\n", myMPIProcID );
            if( myIsDebug ) fprintf(stdout,"  STEP 5 ...else\n");
            if( outPort >= myNextModuleID[iModule]->size() ) {
              throw( cseis_geolib::csException("ERROR in csRunManager:runExecPhase(): output port number exceeds number of output ports: (module #%d), %d %d %s",
                                               iModule, outPort, myNextModuleID[iModule]->size(), myModules[iModule]->getName()) );
            }
            // STEP (6) Find correct combination: Output port of current module  <-->  Input port of next module
            int nextModuleID = myNextModuleID[iModule]->at(outPort);
            if( nextModuleID >= myNumModules ) { // Next module 'pointer' of given output port points beyond last module
              if( myIsDebug ) fprintf(stdout,"  STEP 7a ...is last module\n");
              module->lastModuleTraceCleanup();
              if( module->finishedProcessing() ) {
                iModule = stack_moduleIndex.pop();  // returns myNumModules if empty
              } // else keep iModule, process next trace(s)
            }
            else {
              int inPort = 0;
              while( myPrevModuleID[nextModuleID]->at(inPort) != iModule ) {  // Search matching input port
                inPort++;
              }
              // STEP (7) Move processed traces from current to next module, using correct output/input ports
              myModules[nextModuleID]->moveTracesFrom( myModules[iModule], inPort );
              
              if( !module->finishedProcessing() ) {
                if( myIsDebug ) fprintf(stdout,"  STEP 7b ...not finished processing yet\n");
                stack_moduleIndex.push( iModule );
                iModule = myNextModuleID[iModule]->at( outPort );
              }
              else if( !forceToRun ) {
                iModule = myNextModuleID[iModule]->at( outPort );
              }
              else {
                iModule += 1;  // If input has finished reading in traces, make sure no module is missed on the way down during the last pass
              }
            }
            //MPIDEBUG myLog->line("MPI %d:  END OF STEP 6, module #%d\n", myMPIProcID, iModule );
          } // END else (Step 5)
        } // END else (Step 3)
      }  // while( iModule < myNumModules ) {
      //MPIDEBUG myLog->line("MPI %d: LAST STEP, IsFinished %d\n", myMPIProcID, isInputFinished);
    }  // while( !isInputFinished ) {
  }
  catch( cseis_geolib::csException& e ) {
    if( iModule >= myNumModules || iModule < 0 ) {
      iModule = 0;
    }
    fprintf(stderr,"\nException caught while running exec phase, module #%2d %s. System message: \n%s\n",
            iModule+1, myModules[iModule]->getName(), e.getMessage());
    myLog->line("\nException caught while running exec phase, module #%2d %s. System message: \n%s",
                iModule+1, myModules[iModule]->getName(), e.getMessage());
    mpi_terminateExecPhase();
    exit( -1 );
  }
  //  catch(...) {
  //  printf("Unknown exception occurred\n");
  // }
  //MPIDEBUG myLog->line("MPI %d: EXITED exec phase\n", myMPIProcID);
  
  // Run cleanup phase. Cleaning up allocated memory in modules
  for( int iModule = 0; iModule < myNumModules; iModule++ ) {
    myModules[iModule]->submitCleanupPhase( myLog );
  }

  execPhaseCompletionTasks();
  //
  //
  // END Exec processing loop
  //
  //***********************************************************************
  
  mpi_terminateExecPhase();
  return returnFlag;
}

void csRunManager::execPhaseCompletionTasks() {
  //----------------------------------------------------------------------
  // Dump exec phase info
  //
  myLog->line( "\n------------------------------------------------------------" );
  myLog->line( "Exec phase summary\n" );
  myLog->line( "  #  Module                Traces in  Traces out     CPU time   CPU time all" );

  double timeCPUExecAll = 0.0;
  for( int iModule = 0; iModule < myNumModules; iModule++ ) {
    csModule* module = myModules[iModule];
    double timeCPU  = module->getExecPhaseCPUTime();
    timeCPUExecAll += timeCPU;
    myLog->line( "%3d  %-19s %11d %11d %12.3f   %12.3f", iModule+1, module->getName(),
                 module->numIncomingTraces(), module->numProcessedTraces(), timeCPU, timeCPUExecAll );
  }
  myLog->line( "\n------------------------------------------------------------\n" );
  
  time_t timer = time(NULL);
  myLog->line( " Total processing time:  %14.6f seconds / %10.3f hours\n", myTimerCPU->getElapsedTime(), myTimerCPU->getElapsedTime()/3600.0 );
  myLog->line( " Date: %s\n", asctime(localtime(&timer)) );
  myLog->line( " Trace allocation summary:");
  myMemoryPoolManager->dumpSummary( myLog->getFile() );
  myLog->line( "\n");
  myLog->line( "\n End of log." );  
  myLog->line( "================================================================================" );
}

void csRunManager::mpi_terminateExecPhase() {
#ifdef USE_MPI
  // START of MPI code
  if( myMPIProcID == 0 ) {
    int terminate = -1;
    int tag = 99;
    for( int mpiTargetProc = 1; mpiTargetProc < myMPINumProc; mpiTargetProc++ ) {
      MPI_Send( &terminate, 1, MPI::INT, mpiTargetProc, tag, MPI_COMM_WORLD );
    }
  }
  // END of MPI code
#endif
}

//**********************************************************************
//
// Helper methods
//
//**********************************************************************

/**
 * Check all user defined parameters/values etc
 *
 */
bool csRunManager::checkParameters( char const* moduleName, csParamDef const* paramDef, cseis_geolib::csVector<csUserParam*>* userParams, csLogWriter* log ) {
  bool returnFlag = true;
  if( paramDef->module() == NULL ) {
    log->line("Program bug: Module name not defined in parameter definition. Should be '%s'", moduleName );
    return false;
  }
  if( strcmp( moduleName, paramDef->module()->name() ) ) {
    log->line("Program bug: Inconsistent module names between parameter definition: '%s', and method definition: '%s'\n", paramDef->module()->name(), moduleName );
    returnFlag = false;
  }

  cseis_geolib::csVector<csParamDescription const*> paramDefList(5);  // User parameter definitions for current module, e.g. name & description of parameter, such as 'filename' 'Input file name'  (set in module 'params' method)
  cseis_geolib::csVector<csParamDescription const*> valueDefList(2);  // Definitions of all values for current user parameter (set in module 'params' method)
  cseis_geolib::csVector<csParamDescription const*> optionList(3);
  cseis_geolib::csVector<std::string> userValueList(2);         // All user specified values for current parameter

  paramDef->getParameters( &paramDefList );
  int nDefinedParams = paramDefList.size();

  int nUserParams = userParams->size();
  cseis_geolib::csEquationSolver equationSolver;
     
  // Go through all parameters specified by user
  for( int iUserParam = 0; iUserParam < nUserParams; iUserParam++ ) {
    csUserParam* userParam = userParams->at(iUserParam);
    int ip = -1;
    // Go through all parameters defined for this module, try to find matching parameter name
    for( int i = 0; i < nDefinedParams; i++ ) {
      if( userParam->getName().compare( paramDefList.at(i)->name() ) == 0 ) {
        ip = i;
        break;
      }
    }
    if( ip < 0 ) {
      if( userParam->getName().compare("debug") == 0 || userParam->getName().compare("version") == 0 ) {
        continue;
      }
      log->line("Error: Unknown user specified parameter: '%s'", userParam->getName().c_str());
      returnFlag = false;
      continue;  // Check other parameters
    }
    // fprintf(stdout,"\nUser param: '%s'  '%s'\n", userParamName, paramDefList.at(ip)->desc() );

    int nUserValues = userParam->getNumValues();  // Number of parameter values SPECIFIED BY USER
    valueDefList.clear();
    paramDef->getValues( ip, &valueDefList );
    int nDefinedValues = valueDefList.size();   // Number of values DEFINED for this parameter

    if( paramDefList.at(ip)->type() == NUM_VALUES_FIXED && nUserValues < nDefinedValues ) {
      log->line("Error: Too few user specified values for parameter '%s'. Required: %d, found: %d", userParam->getName().c_str(), nDefinedValues, nUserValues );
      returnFlag = false;
      continue;
    }
    else if( nUserValues == 0 ) {
      log->line("Error: No value specified for parameter '%s'.", userParam->getName().c_str() );
      returnFlag = false;
      continue;
    }
    // Convert all user specified parameter values to lower case (except for TYPE_STRING values)
    // Also, check all number values
    int minNumValues = std::min( nDefinedValues, nUserValues );
    cseis_geolib::csFlexNumber valueTmp;
    for( int i = 0; i < minNumValues; i++ ) {
      int valueType = valueDefList.at(i)->type();
      if( valueType == VALTYPE_NUMBER || valueType == VALTYPE_HEADER_NUMBER ) {
        std::string text = userParam->getValue( i );
        if( !valueTmp.convertToNumber( userParam->getValue( i ) ) ) {
          if( !equationSolver.prepare( text ) ) {
            if( valueType == VALTYPE_NUMBER ) {
              log->line("Error: User parameter '%s': Value is not recognised as valid number: '%s'", userParam->getName().c_str(), text.c_str() );
              returnFlag = false;
            }
          }
          else {
            userParam->setValue( i, equationSolver.solve() );
          }
        }
      }
      else if( valueType == VALTYPE_OPTION ) {
        //        userParam->setValueType( i, valueType );
        userParam->convertToLowerCase( i );
      }
    }
    if( nUserValues > nDefinedValues ) {
      for( int i = nDefinedValues; i < nUserValues; i++ ) {
        int valueType = valueDefList.at(nDefinedValues-1)->type();
        if( valueType == VALTYPE_OPTION ) {
          userParam->convertToLowerCase( i );
        }
      }
    }
    // For parameters taking a variable number of values, convert the last values (all same type)
    //    if( paramDefList.at(ip)->type() == NUM_VALUES_VARIABLE && valueDefList.at(nDefinedValues-1)->type() == VALTYPE_OPTION ) {
    //      userParam->setValueTypeVariable( nDefinedValues, valueDefList.at(nDefinedValues-1)->type() );
    //    }
    userValueList.clear();
    userParam->getValues( &userValueList );

    // Go through all parameter values, check correctness of OPTION values
    for( int iv = 0; iv < minNumValues; iv++ ) {
      if( valueDefList.at(iv)->type() == VALTYPE_OPTION ) {
        //        char const* userOptionName = userValueList.at(iv).c_str();
        // Check options:
        optionList.clear();
        paramDef->getOptions( ip, iv, &optionList );
        int optionFound = false;
        int nOptions = optionList.size();
        for( int io = 0; io < nOptions; io++ ) {
          if( userValueList.at(iv).compare( optionList.at(io)->name() ) == 0 ) {
            optionFound = true;
            break;
          }
        }
        if( !optionFound ) {
          log->line("Error: Unknown user specified option for parameter '%s':  %s", userParam->getName().c_str(), userValueList.at(iv).c_str() );
          log->write("Valid options are: ");
          for( int io = 0; io < nOptions; io++ ) {
            log->write("'%s' (%s)", optionList.at(io)->name(), optionList.at(io)->desc());
            if( io < nOptions-1 ) log->write(" / ");
          }
          log->line("");
          returnFlag = false;
        }
      } // END if
      else {
        if( valueDefList.at(iv)->type() != VALTYPE_STRING ) {
        }
      }
      //      fprintf(stdout," Name: %s, valueType: %d\n", valueDefList.at(iv)->name(), valueDefList.at(iv)->type());
    }
    
  }

  return returnFlag;
}

bool csRunManager::parseVersionString( std::string& text, int& major, int& minor ) {
  int length   = text.length();
  int indexDot = text.find_first_of('.');
  if( indexDot <= 0 || indexDot == length-1 || indexDot == (int)std::string::npos ) return false;
  //  fprintf(stdout,"Version string: '%s' = '%s'.'%s'\n", text.c_str(), text.substr(0,indexDot).c_str(), text.substr(indexDot+1,length-indexDot-1).c_str() );
  major = atoi( text.substr(0,indexDot).c_str() );
  minor = atoi( text.substr(indexDot+1,length-indexDot-1).c_str() );
  if( minor >= 0 && minor <= 99 && major >= 0 && major <= 99 ) return true;
  return false;
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_RUN_MANAGER_H
#define CS_RUN_MANAGER_H

#include <string>
#include <cstdio>

namespace cseis_geolib {
template<typename T> class csVector;
template<typename T> class csCompareVector;
class csTimer;
class csTable;
}

namespace cseis_system {

class csModule;
class csLogWriter;
class csUserConstant;
class csParamDef;
class csUserParam;
class csMemoryPoolManager;

 struct modInfoStruct {
    modInfoStruct() {
      modIndex = -1; modType = -1;
    }
    modInfoStruct( modInfoStruct const& obj ) {
      modIndex = obj.modIndex; modType = obj.modType;
    }
    modInfoStruct( int m, int t ) {
      modIndex = m; modType = t;
    }
    int modIndex;
    int modType;
  };

/**
 * Run Manager
 *
 * Manages the execution of one Cseis flow.
 * Parses input flow file, runs init and exec phase.
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
class csRunManager {
public:
  /**
  * Constructor
  *
  * @param log   Log writer
  * @param memoryPolicy   ...as defined in csMemoryPoolManager: POLICY_SPEED or POLICY_MEMORY
  * @param isDebug true if extended debug information shall be printed
  */
  csRunManager( csLogWriter* log, int memoryPolicy, int mpiProcID, int mpiNumProc, bool isDebug = false );
  ~csRunManager();
  /**
  * Run initialisation phase for all modules
  * 1) Parse input file, replace user constants etc
  * 2) Read modules and user parameters form input file
  * 3) Set flow sucession, i.e. for IF or SPLIT blocks
  * 4) Run init phase
  *
  * @param filenameFlow (i) The file name of the current flow
  * @param f_flow       (i) Pointer to the current flow file
  * @param globalConstList (i) List of global constants that shall be replaced in input flow
  */
  int runInitPhase( char const* filenameFlow, std::FILE* f_flow, cseis_geolib::csCompareVector<csUserConstant>* globalConstList );
  int runExecStartPhase();
  /**
  * Run execution phase for all modules
  */
  int runExecPhase();
  int mpi_runExecPhase();
  void mpi_terminateExecPhase();
  void execPhaseCompletionTasks();
  
  static bool checkParameters( char const* moduleName, csParamDef const* paramDef, cseis_geolib::csVector<csUserParam*>* userParams, csLogWriter* log );
private:
  csRunManager(csRunManager const& obj);
  csLogWriter* myLog;
  cseis_geolib::csTable const** myTables;
  int myNumTables;
  /// Array of pointers to all modules in one flow
  csModule** myModules;
  /// Array of lists that hold module ID of the next successive module(s)
  cseis_geolib::csVector<int>** myNextModuleID;
  cseis_geolib::csVector<int>** myPrevModuleID;
  int myNumModules;
  bool myIsDebug;
  int myMPIProcID;
  int myMPINumProc;
  /**
   * Parse version string.
   * @return false if string does not contain valid version string
   */
  static bool parseVersionString( std::string& text, int& major, int& minor );
  /**
  * Check all user parameters against parameter definition specified in module parameter method "void params_mod_MODULENAME_( csParamDef* pdef )"
  */
  csMemoryPoolManager* myMemoryPoolManager;
  cseis_geolib::csTimer* myTimerCPU;
};

} // namespace

#endif


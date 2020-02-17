/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_MODULE_H
#define CS_MODULE_H

#include <string>
#include "cseis_defines.h"
#include "csSuperHeader.h"

namespace cseis_geolib {
  template<typename T> class csVector;
  template<typename T> class csQueue;
  class csFlexNumber;
  class csTable;
}

namespace cseis_system {

class csTraceHeaderDef;
class csExecPhaseDef;
class csTraceHeaderDef;
class csTraceGather;
class csTrace;
class csUserParam;
class csParamManager;
class csParamDef;
class csExecPhaseEnv;
class csMemoryPoolManager;

/**
* Central Cseis class
* An csModule object holds all information regarding one Cseis processing module, during the flow execution
*
*  - Holds all information for one module, such as trace header definition, exec phase definition etc
*  - Manages access to module methods (param, init, exec phase)
*  - Manages complex trace flow between this and previous/following modules
*
* @author Bjorn Olofsson
* @date   2007
*/
class csModule {
public:
  csModule( std::string const& name, int id_number, csMemoryPoolManager* memManager, int mpiProcID, int mpiNumProc );
  ~csModule();
  /// Get parameter definition. These are option names and values that can be specified by the user in job flow
  void getParamDef( csParamDef* paramDef );
  /// Submit init phase
  void submitInitPhase( csParamManager* paramManager, csLogWriter* writer, cseis_geolib::csTable const** tables, int numTables );

  /// Submit start exec phase. Return false if error occurred.
  bool submitExecStartPhase( csLogWriter* writer );
  /// Return true if module is ready to submit exec phase (depends on state of trace buffer)
  bool isReadyToSubmitExecPhase( bool forceToRun ) const;
  /// Submit exec phase. Return false if no trace was processed.
  bool submitExecPhase( bool forceToProcess, csLogWriter* writer, int& port );

  void mpi_decompressInitPhaseObjects( char const* dataSuperHdr, char const* dataHdrDef );
  bool mpi_submitExecPhase( csLogWriter* writer, int& port );
  /// Submit clean-up phase. Returning false means there was some unspecified error, probably a program bug.
  void submitCleanupPhase( csLogWriter* writer, bool reportMissing = true );
  
  /// Move seismic traces from module 'module' that is connected to this module at input port 'inPort'
  void moveTracesFrom( csModule* module, int inPort );
  /// Call to clean up last module in flow. Call after each exec phase submission
  void lastModuleTraceCleanup();  //...
  /// @return true if there are still unprocessed traces waiting in the trace gather
  bool finishedProcessing() const;
  /// @return total number of incoming traces
  inline long numIncomingTraces() const { return( myTotalNumIncomingTraces ); }
  /// @return total number of processed traces
  inline long numProcessedTraces() const { return( myTotalNumProcessedTraces ); }
  /// @return module name
  inline char const* getName() const { return myName.c_str(); }
  /// @return module name
  inline std::string const* getNameString() const { return &myName; }
  /// @return module type
  inline int getType() const { return myModuleType; }
  /// @return module type
  void setExecType( int execType );
  /// @return module type
  int getExecType() const;
  /// @return trace header definition object
  csTraceHeaderDef const* getHeaderDef() const;
  /// @return exec phase definition object
  csExecPhaseDef const* getExecPhaseDef() const;
  /// @return super header
  csSuperHeader const* getSuperHeader() const { return mySuperHeader; }
  /**
  * Set input port modules, merge trace headers and super headers etc...
  * @param moduleList: List of all input port modules
  */
  void setInputPorts( cseis_geolib::csVector<csModule const*> const* moduleList );

  /// Set DEBUG flag for this module
  void setDebugFlag( bool doDebug );
  /// ...used for debugging purposes
  int tempNumTraces();
  /// @return CPU time used during module's exec phase
  inline double getExecPhaseCPUTime() const { return myTimeExecPhaseCPU; }
  /**
  * Set module version number
  * @param major: Major version (1-99)
  * @param minor: Minor version (0-99)
  */
  void setVersion( int major, int minor ) { myVersion[MAJOR] = major; myVersion[MINOR] = minor; }
  /// @return version string
  std::string versionString() const;
  /**
   * Set number of MPI processors
   */
  void mpiSetNumProc( int mpiNumProc );
  /**
   * Set MPI support ON/OFF
   */
  void setMPISupport( bool isMPISupported );
  /**
   * Does this module support MPI?
   * @return true: Yes, false: No
   */
  bool isMPISupported() const;
 private:
  /// Disabled copy constructor
  csModule( csModule const& module );  
  /// Initialize. Method is called by constructors
  void init();
  /// Retrieve module methods: Parameter definition method and init method
  void retrieveParamInitMethods();
  /// Maximum number of output ports
  static int const MAX_NUM_PORTS = 2;

  /// Manages some of the memory intensive allocation centrally
  csMemoryPoolManager* myMemoryPoolManager;

  /// Definitions of trace headers that are available to this module
  csTraceHeaderDef* myHeaderDef;
  /// Definitions for operation of exec phase
  csExecPhaseDef* myExecPhaseDef;
  /// Module name
  std::string myName;
  /// Super header (see class for description)
  csSuperHeader* mySuperHeader;
  /// Number of traces that have been processed already and are ready to be passed on to another module
  int myNumTracesToBePassed;

  /// Trace gather containing traces that are still to be processed by this module. This object is passed to the module exec method. (For single trace modules: the single trace inside the gather)
  csTraceGather* myTraceGather;
  /// Trace queue: Traces that are waiting to be moved to the actual trace gather.
  /// These are traces that were recently moved to this module or shall not be processed yet,
  /// for example because they belong to a different ensemble than the traces in myTraceGather.
  cseis_geolib::csQueue<csTrace*>* myTraceQueue;
  /// Contains pointers to objects used in module exec phase
  csExecPhaseEnv* myExecEnvPtr;
  /// Accumulated CPU time taken by module's exec phase
  double myTimeExecPhaseCPU;

  //-----------------------------------------------------------------------------------------
  // Fields relating to ensemble breaks -- maybe these could be wrapped up in another class? Maybe csExecPhaseDef?

  void mpi_distribute_traces_before_exec_phase( csLogWriter* writer );
  void mpi_distribute_traces_after_exec_phase( csLogWriter* writer );
  /**
  * Helper method: Add new trace to trace gather object
  * @param trace: Trace to add
  * @param inPort: Input port
  */
  void addNewTraceToGather( csTrace* trace, int inPort );
  /**
  * Helper method: Add new trace to trace queue object
  * @param trace: Trace to add
  * @param inPort: Input port
  */
  void addNewTraceToQueue( csTrace* trace, int inPort );

  /// Special method to update trace gathers for 'ensemble' modules
  void updateTracesEnsembleModule();
  /**
  * Is this trace part of the current ensemble?
  * @return true if trace is part of currently collected ensemble
  *    If this is the first trace in a new ensemble, the ensemble key value is set
  * @return false if trace is not part of current ensemble. If this is the case, ensemble is set to 'full' state,
  *   and next ensemble key header value is set
  */
  bool traceIsPartOfCurrentEnsemble( csTrace* trace );

  /// Ensemble keys: These specify ensemble breaks and are defined in super header
  /// Array of ensemble key header indexes (in trace header), for quick access of ensemble keys
  int* myEnsembleKeyHeaderIndex;
  /// Array of ensemble key values for currently processed trace
  cseis_geolib::csFlexNumber* myEnsembleKeyValue;
  /// Array of ensemble key values for next processed trace
  cseis_geolib::csFlexNumber* myNextEnsembleKeyValue;
  ///
  cseis_geolib::csFlexNumber* myHelperHdrValues;

  /// true if current ensemble is full, ready to be processed
  bool myIsEnsembleFull;
  /// true if processing of current trace gather is finished
  bool myIsFinishedProcessing;

  /// Parameter method (function pointer)
  MParamPtr myMethodParam;
  /// INIT phase method (function pointer)
  MInitPtr myMethodInit;
  /// EXEC START phase method (function pointer)
  MExecStartPtr myMethodExecStart;
  /// EXEC phase method (function pointer)
  MExecPtr myMethodExec;
  /// CLEANUP phase method (function pointer)
  MCleanupPtr myMethodCleanup;
  //--------------------
  // Module structure  
  //
  /// Unique module ID
  int const myUniqueID;
  /// Number of output ports. Only 2 ports supported in the current version
  int myNumOutputPorts;
  /// Number of input ports
  int myNumInputPorts;
  /// Module type (IF, ENDIF...)
  int myModuleType;
  static int const MAJOR = 0;
  static int const MINOR = 1;
  /// Module version
  int myVersion[2];
//---------------------------------------------------------------
// Members that may be put into their own object
  /// Total number of traces processed
  long myTotalNumProcessedTraces;
  /// Total number of incoming traces
  long myTotalNumIncomingTraces;
//---------------------------------------------------------------
// MPI variables
  /// Number of MPI processors
  int myMPINumProc;
  /// MPI process number
  int myMPIProcID;
  /// For FIXED-trace modules: Number of traces to collect before running exec phase
  int myNumFixedTracesToCollect;
};

} // namespace
#endif

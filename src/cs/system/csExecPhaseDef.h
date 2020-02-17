/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_EXEC_PHASE_DEF_H
#define CS_EXEC_PHASE_DEF_H

#include <string>
#include "csFlexNumber.h"

namespace cseis_system {
  class csTraceGather;
  class csSuperHeader;
  class csTraceHeaderDef;

/**
 * Exec phase definition.
 * This class 'defines' how to run the exec phase of one module.
 * An object of this class is passed to each module's init and exec phase subroutine for communication
 * between base system and module.
 * In the init phase, the module's 'exec type' and 'trace selection mode' (see cseis_defines.h) have to be defined.
 * Also, any memory which need to be accessed during both the init and exec phase must be set ('variables' field).
 * In the exec phase, 
 * If the method isCleanup() returns true this indicates to the module's exec phase is completed. No traces are passed
 * in or out when this call is made. The module's exec phase subroutine should clean up any allocated memory and return.
 *
 * @author Bjorn Olofsson
 * @date   2007
*/
class csExecPhaseDef {
public:
  csExecPhaseDef( std::string& moduleName, int mpiProc );
  ~csExecPhaseDef();
  /**
   * Set trace selection mode.
   * This is required for multi-trace modules. TRCMODE_ENSEMBLE means that the base system collects and passes
   * a full ensemble to the module's exec phase subroutine, as defined by the user-set ensemble trace header.
   * For TRCMODE_FIXED, a fixed number of traces are collected by the base system and passed on to the module's
   * exec phase subroutine.
   * @param mode    TRCMODE_FIXED or TRCMODE_ENSEMBLE
   */
  void setTraceSelectionMode( int mode );
  /**
   * Set trace selection mode.
   * This is required for multi-trace modules. TRCMODE_ENSEMBLE means that the base system collects and passes
   * a full ensemble to the module's exec phase subroutine, as defined by the user-set ensemble trace header.
   * For TRCMODE_FIXED, a fixed number of traces are collected by the base system and passed on to the module's
   * exec phase subroutine.
   * @param mode    TRCMODE_FIXED or TRCMODE_ENSEMBLE
   * @param nTraces Number of traces to pass to module each time exec phase is called (for TRCMODE_FIXED)
   */
  void setTraceSelectionMode( int mode, int nTraces );
  /**
   * Set exec phase type
   * EXEC_TYPE_INPUT is reserved for 'input' modules, i.e. modules which can act as first modules in a flow.
   * All other modules have type=EXEC_TYPE_NORMAL.
   * Modules with EXEC_TYPE_NORMAL receive a trace gather from the base system:
   * Either an ensemble or a fixed number of traces, based on the setting for the 'trace selection mode'.
   *
   * @param theExecType EXEC_TYPE_INPUT, EXEC_TYPE_NORMAL
   */
  void setExecType( int theExecType );
  /**
   * Save variables pointer that stores fields which need to be available to both init and exec phase
   * The pointer can point to any allocated block of memory (allocated during init phase) which shall
   * be accessed during the exec phase.
   * The module programmer is responsible for allocating and freeing all memory that this pointer points to.
   *
   * @param variables Variables pointer
   */
  void setVariables( void* variables ) { varPtr = variables; }
  /**
   * @return Access to variables pointer
   */
  inline void* variables() { return varPtr; }
  /**
   * @return true if this call to the module is the clean-up phase
   */
  inline bool isCleanup() const { return myIsCleanup; }
  /**
   * @return true if this module supports MPI
   */
  inline bool isMPISupported() const { return myIsMPISupported; }
  inline int getMPIProcID() const { return myMPIProcID; }
  /**
   * Debug flag is set/not set.
   * The module programmer can use this flag to print out special diagnostic messages that will help
   * the user to debug problems occurring during the job when calling the current module.
   * The user specifies 'debug' in the flow as one module parameter to set this flag.
   * @return true if the user has specified the debug flag in the seismic flow
   */
  inline bool isDebug() const { return myIsDebug; }
  /**
   * @return The module's name
   */
  inline std::string moduleName() const { return myModuleName; }
  /**
   * @return exec phase 'type': EXEC_TYPE_INPUT, EXEC_TYPE_NORMAL
   */
  inline int execType() const { return myExecType; }
  /**
   * In the module's Exec phase, call this method before returning in order to instruct
   * The Seaseis base system that the module still keeps traces in local memory that need to be passed
   * on to the next module.
   * This method has to be called each time before returning to give this instruction.
   * It makes sure that the base system knows it cannot "clean up" this module yet, even if
   * all earlier modules have finished processing traces.
   *
   * Note that if this method is ALWAYS called at the end one a module's Exec phase, the flow will
   * never finish an infinite loop.
   */
  void setTracesAreWaiting();
  void setTracesAreWaiting( bool areWaiting );
  /**
   * Has the 'tracesAreWaiting' flag been set?
   * Each time the exec phase of one module is called, this method return false (because it is reset each time).
   * @return true if module has traces stored internally that still need to be passed on to the base system
   */
  bool tracesAreWaiting() const;
  /**
   * @return true if this is the last call to the module's exec phase
   */
  bool isLastCall() const;
  int getTraceMode() const;

  void setTraceMode( int traceMode );
  void setLastCall( bool isLastCall );
  void setCleanUp( bool isCleanUp );
  void setMPISupport( bool isSupported );
  int getNumTraces() const { return numTraces; }

  /// class csModule needs access to the private fields of this class
  friend class csModule;

private:
  /// Module name
  std::string myModuleName;

  /// Trace selection mode defining which traces shall be passed to a module: TRCMODE_FIXED or TRCMODE_ENSEMBLE
  int traceMode;
  /// 'Type' of module execution method: EXEC_TYPE_INPUT, EXEC_TYPE_NORMAL
  int myExecType;
  /// Number of traces to pass to module (trace selection mode = TRCMODE_FIXED)
  int numTraces;
  /// true when all traces have been processed and modules shall 'clean up' allocated memory
  bool myIsCleanup;
  /// Debug flag. Is set if user has specified 'debug' parameters in flow
  bool myIsDebug;
  /// Pointer to parameters managed by module source code. Do not free this pointer within this class!
  void* varPtr;
  /// true if module has traces stored that still need to be passed on to base system (and hence the subsequent modules)
  bool myTracesAreWaiting;
  /// true if this is the last call from the base system to this module's exec phase (different from cleanup phase)
  bool myIsLastCall;
  /// true if module supports mpi
  bool myIsMPISupported;
  /// MPI process number
  int myMPIProcID;
};

} // namespace

#endif


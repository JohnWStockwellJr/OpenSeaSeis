/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CSEIS_DEFINES_H
#define CSEIS_DEFINES_H

#define CS_DEBUG

#include "geolib_defines.h"

namespace cseis_system {

class csParamManager;
class csTrace;
class csTraceGather;
class csParamDef;
class csInitPhaseEnv;
class csExecPhaseEnv;
class csLogWriter;
 
typedef void (*MParamPtr) ( csParamDef* paramDef );
typedef void (*MInitPtr) ( csParamManager* param, csInitPhaseEnv* initPhaseEnv, csLogWriter* writer );
/**
* Start of exec phase method, called once per module
*/
typedef bool (*MExecStartPtr) ( csExecPhaseEnv* env, csLogWriter* writer );
/**
* Exec phase method
* @param traceGather: Gather containing a number of traces (for example one 'ensemble')
* @param port:        Output port number
* @param numTrcToKeep: Number of traces to keep until next pass of exec method
*        The last numTrcToKeep traces in trace gather are not moved immediately to the next module
*        Instead, they stay until the next pass of the exec method
* @param env:     Exec phase environment.
* @param writer:  Output log writer
*/
typedef void (*MExecPtr) ( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer );
/**
* Clean-up phase method
*/
typedef void (*MCleanupPtr) ( csExecPhaseEnv* env, csLogWriter* writer );

// Execution type for modules
static int const EXEC_TYPE_INPUT      = 111;   // Input module
static int const EXEC_TYPE_NORMAL     = 333;   // Multi trace module
static int const EXEC_TYPE_UNKNOWN    = 0;

/// Trace mode types
/// This mode defines how to determine which traces shall be passed to the module exec phase method
static int const TRCMODE_FIXED    = -11;  // Fixed number of traces (1 or more) are passed to module exec phase method
static int const TRCMODE_ENSEMBLE = -22;  // Full ensemble/gather is passed to module exec phase method
static int const TRCMODE_UNKNOWN  = 0;

/// Module types
// Keep relative numbering!
static int const MODTYPE_UNKNOWN   = -1;
static int const MODTYPE_IF        = 0;
static int const MODTYPE_ELSEIF    = 1;
static int const MODTYPE_ELSE      = 2;
static int const MODTYPE_ENDIF     = 3;
static int const MODTYPE_SPLIT     = 4;
static int const MODTYPE_ENDSPLIT  = 5;
static int const MODTYPE_1PORT     = 10;
static int const MODTYPE_2PORTS    = 11;

/// Time domain unit types
static int const T_UNIT_UNKNOWN = -1;
static int const T_UNIT_TIME    = 220;
static int const T_UNIT_SAMPLE  = 221;

static int const MUTE_FRONT = 330;
static int const MUTE_END   = 331;

/// Maximum length of parameter names (this is just a recommendation)
static int const MAX_LENGTH_PARAMETER = 20;
static int const MAX_LINE_LENGTH = 1024;
static char const LETTER_MODULE = '$';   // Letter preceding module name in Cseis flow
static char const LETTER_COMMENT = '#';  // Letter indicating comment line in Cseis flow
 
#define CSEIS_VERSION "3.00"

}

#endif


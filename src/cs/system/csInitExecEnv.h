/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_INIT_EXEC_ENV_H
#define CS_INIT_EXEC_ENV_H

#include <string>

namespace cseis_geolib {
  class csTable;
}

namespace cseis_system {

class csTraceHeaderDef;
class csExecPhaseDef;
class csSuperHeader;

/**
 * Exec phase environment.
 *
 * Helper class carrying information between Cseis base system and processing modules
 * This helps to keep the interface between base system and the module subroutines stable
 * even if more information needs to be passed between the two in a future version.
 * @author Bjorn Olofsson
 * @date   2007
 */
class csExecPhaseEnv {
 public:
  csExecPhaseEnv(
    csTraceHeaderDef const* headerDef,
    csExecPhaseDef*         execPhaseDef,
    csSuperHeader const*    superHeader );
  ~csExecPhaseEnv();

  csTraceHeaderDef const* headerDef;
  csExecPhaseDef*         execPhaseDef;
  csSuperHeader const*    superHeader;
  void addError() { errorCounter++; }
  int errorCount() { return errorCounter; }
 private:
  int errorCounter;
  csExecPhaseEnv();
};

/**
 * Init phase environment
 * Helper class carrying information between Cseis base system and processing modules 
 * This helps to keep the interface between base system and the module subroutines stable
 * even if more information needs to be passed between the two in a future version.
*/
class csInitPhaseEnv {
 public:
  csInitPhaseEnv(
    csTraceHeaderDef* headerDef,
    csExecPhaseDef*   execPhaseDef,
    csSuperHeader*    superHeader,
    cseis_geolib::csTable const** tables,
    int               numTables );
  ~csInitPhaseEnv();
  csTraceHeaderDef* headerDef;
  csExecPhaseDef*   execPhaseDef;
  csSuperHeader*    superHeader;
  void addError() { errorCounter++; }
  int errorCount() { return errorCounter; }
  cseis_geolib::csTable const* getTable( std::string const& tableName ) const;
 private:
  int errorCounter;
  int myNumTables;
  cseis_geolib::csTable const** myTables;
  csInitPhaseEnv();
};

} // namespace

#endif


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SU_ARGUMENTS_H
#define CS_SU_ARGUMENTS_H

#include <string>

namespace cseis_geolib {
  class csSegyTraceHeader;
  class csSegyHdrMap;
  template <typename T> class csVector;
}

namespace cseis_su {

class csSUTraceManager;

/**
 * Carries SU module's command line arguments from CSEIS to SU
 */
class csSUArguments {
 public:
  csSUArguments();
  ~csSUArguments();
  /**
   * Set arguments. Populate member field 'argv'
   * @param argvList  List of arguments
   */
  void setArgv( cseis_geolib::csVector<std::string>* argvList );
 public:
  /// Trace manager for data transfer from CSEIS to SU
  csSUTraceManager* cs2su;
  /// Trace manager for data transfer from SU to CSEIS
  csSUTraceManager* su2cs;
  /// Command line arguments as they would appear when running SU form command line
  char** argv;
  /// Number of command line arguments
  int argc;
  /// Is debug mode?
  int debugFlag;
};

} // END namespace

#endif

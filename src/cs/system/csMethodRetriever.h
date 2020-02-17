/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_METHOD_RETRIEVER_H
#define CS_METHOD_RETRIEVER_H

#include "cseis_defines.h"
#include <string>

namespace cseis_geolib {
  template<typename T> class csVector;
}

namespace cseis_system {

/**
* Retrieval of module methods
* Use static member function of this class to retrieve function pointers to the standard methods of Cseis processing modules
*
* @author Bjorn Olofsson
* @date   2007
*/
class csMethodRetriever {
public:
  /**
  * Retrieve function pointer to parameter method
  * @param name (i) Module name whose method shall be retrieved
  */
  static MParamPtr getParamMethod( std::string const& name, int verMajor, int verMinor );
  static MParamPtr getParamMethod( std::string const& name, std::string versionString );
  static MParamPtr getParamMethod( std::string const& name );
  /**
  * Retrieve function pointers to parameter method and initialisation method
  * @param name (i) Module name whose method shall be retrieved
  */
  static void getParamInitMethod( std::string const& name, int verMajor, int verMinor, MParamPtr& param, MInitPtr& init );
  /**
  * Retrieve function pointer to execution method
  * @param name (i) Module name whose method shall be retrieved
  */
  static void getExecMethod( std::string const& name, int verMajor, int verMinor, MExecStartPtr& execStart, MExecPtr& exec, MCleanupPtr& cleanup );

  /**
  * NOTE: This method is only used when statically linking Cseis modules
  * @return number of standard modules defined
  */
  static int getNumStandardModules();
  /**
  * NOTE: This method is only used when statically linking Cseis modules
  * @return pointer to list of standard modules
  */
  static std::string const* getStandardModuleNames();

 private:
  static MParamPtr     getParamMethod( std::string const& name, void* handle );
  static MInitPtr      getInitMethod( std::string const& name, void* handle );
  static MExecPtr      getExecMethod( std::string const& name, void* handle );
  static MExecStartPtr getExecStartMethod( std::string const& name, void* handle );
  static MCleanupPtr   getCleanupMethod( std::string const& name, void* handle );

  static int getMethodIndex( std::string const& name );
  static int const METHOD_NOT_FOUND = -33;

  static MParamPtr getParamMethod( std::string const& name, char const* soName );

  csMethodRetriever();
  ~csMethodRetriever();  
};

} // namespace
#endif


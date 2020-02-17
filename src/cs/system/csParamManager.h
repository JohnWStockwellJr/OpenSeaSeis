/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_PARAM_MANAGER_H
#define CS_PARAM_MANAGER_H

#include <string>

namespace cseis_geolib {
  template <typename T> class csVector;
}

namespace cseis_system {

class csUserParam;
class csLogWriter;

/**
 * User Parameter manager.
 * Manages access to user parameters.
 * Encapsulates some of the functionality required for user parameters, for use in calls to Cseis module methods.
 *  What are the parameter's 'values' and 'lines':
 * Each defined user parameter in a module has a name and a list of values that the user needs to specify in the flow.
 * If the parameter requires more than one value, all values are specified in one line, separated by spaces.
 * Some modules may allow the parameter to be specified several times, i.e. in several 'lines'.
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
class csParamManager {

static int const NOT_FOUND = -1;

public:
  csParamManager( cseis_geolib::csVector<csUserParam*> const* userParamList, csLogWriter* log );
  ~csParamManager();

  /**
   * @return Total number of specified user parameters
   */
  int  numParameters() const;
  /**
   * How many calls were made during the module's init phase to a specific parameter value?
   * This can be used at run-time to find out if a parameter value defined by the user is actually used.
   *
   * @param index Index of parameter value (0 for first value, >0 if parameter has several values)
   * @return Number of calls that were made
   */
  int  getNumValueCalls( int index ) const;
  /**
  * Get number of lines in module setup that specify the given parameter (for multi-line parameters)
  * @return Number of lines
  */
  int getNumLines( char const* name );
  /**
  * Get number of values specified for this user parameter
  * @param name      (i)  Parameter name
  * @param indexLine (i)  Line index for this parameter
  * @return Number of values defined for this parameter
  */
  int getNumValues( char const* name, int indexLine = 0 );

  /**
   * Retrieve integer value specified for this user parameter
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  inline void getInt( char const* name, int* value, int indexValue = 0 ) {
    getIntAtLine( name, value, 0, indexValue );
  }
  /**
   * Retrieve float value specified for this user parameter
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  inline void getFloat( char const* name, float* value, int indexValue = 0 ) {
    getFloatAtLine( name, value, 0, indexValue );
  }
  /**
   * Retrieve double value specified for this user parameter
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  inline void getDouble( char const* name, double* value, int indexValue = 0 ) {
    getDoubleAtLine( name, value, 0, indexValue );
  }
  /**
   * Retrieve string value specified for this user parameter
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  inline void getString( char const* name, std::string* value, int indexValue = 0 ) {
    getStringAtLine( name, value, 0, indexValue );
  }

  /**
   * Retrieve integer value specified for this user parameter.
   * Retrieve it from the parameter at a specific line.
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexLine  (i) Index of parameter line.
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  void getIntAtLine( char const* name, int* value, int indexLine, int indexValue = 0 );
  /**
   * Retrieve float value specified for this user parameter.
   * Retrieve it from the parameter at a specific line.
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexLine  (i) Index of parameter line.
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  void getFloatAtLine( char const* name, float* value, int indexLine, int indexValue = 0 );
  /**
   * Retrieve double value specified for this user parameter.
   * Retrieve it from the parameter at a specific line.
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexLine  (i) Index of parameter line.
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  void getDoubleAtLine( char const* name, double* value, int indexLine, int indexValue = 0 );
  /**
   * Retrieve text value specified for this user parameter.
   * Retrieve it from the parameter at a specific line.
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexLine  (i) Index of parameter line.
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  void getText( char const* name, char** value, int indexLine = 0, int indexValue = 0 );
  /**
   * Retrieve string value specified for this user parameter.
   * Retrieve it from the parameter at a specific line.
   * @param name       (i) Parameter name
   * @param value      (o) Value to be retrieved
   * @param indexLine  (i) Index of parameter line.
   * @param indexValue (i) Index of parameter value. Retrieve (i+1)'th parameter value
   */
  void getStringAtLine( char const* name, std::string* value, int indexLine = 0, int indexValue = 0 );

  /**
   * Retrieve all user parameter values from one line
   * @param name (i)   Parameter name
   * @valueList  (o)   List of all parameter values as text strings
   * @indexLine  (i)   Line index for this parameter
   */
  void getAll( char const* name, cseis_geolib::csVector<std::string>* valueList, int indexLine = 0 );
  /**
   * @param name Parameter name
   * @return true if this parameter exists in module setup
   */
  bool exists( char const* name );
private:
  csParamManager( csParamManager const& obj );
  bool getParamValues( char const* name, int index, cseis_geolib::csVector<std::string>* valueList, bool exOnFail = true );
  std::string const getParamValue( char const* name, int indexLine, int indexValue, bool exOnFail = true );
  int getParamIndex( char const* name, int indexLine, bool exOnFail = true );

  int myIndexCurrentParam;
  char const* myNameCurrentParam;
  int myLineIndexCurrentParam;
  
  csLogWriter* myLog;
  /// Number of user calls to this
  int* myNumValueCalls;
  /// Do not free!
  cseis_geolib::csVector<csUserParam*> const* myUserParamList;
};

}
#endif


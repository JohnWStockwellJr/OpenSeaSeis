/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_PARAM_DEF_H
#define CS_PARAM_DEF_H

#include "csParamDescription.h"
#include "csVector.h"
#include "geolib_string_utils.h"
#include <cstdio>

namespace cseis_geolib {
  template <typename T> class csVector;
}

namespace cseis_system {
  /// Parameter type:
  static const int NUM_VALUES_FIXED    = 11;
  static const int NUM_VALUES_VARIABLE = 12;
  /// Value type:
  static const int VALTYPE_OPTION = 21;
  static const int VALTYPE_NUMBER = 22;
  static const int VALTYPE_STRING = 23;
  static const int VALTYPE_HEADER_NUMBER = 24;

/**
 * Parameter definition
 *
 * Defines all parameters for one module, including a short description/help
 * - Module name
 * - Parameter names & types
 * - Value types
 * - Option names (for option-type values)
 * - Default parameter values
 *
 * Definitions
 *   Parameter: One line in job flow. A parameter consists of a parameter name and associated 'values'
 *   Values:    Each parameter has a set of 'values'. Each value can be free-form text, or a specific option
 *   Option:    Some values 
 * NOTE:
 *  All methods in this class ensure that module names are UPPER CASE and parameter & option names are LOWER CASE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
class csParamDef {
public:
  static const int LOWER_CASE = 111;
  static const int UPPER_CASE = 222;
  class Param;
  class Value;

 private:
  // Version number indices:
  static int const MAJOR  = 0;
  static int const MINOR  = 1;

public:
  csParamDef();
  ~csParamDef();
  /**
   * Set module name and general description
   *
   * @param moduleName Module name (capital letters)
   * @param desc       Short description of module's purpose
   * @param descExtra  Additional short description of module's purpose
   */
  void setModule( char const* moduleName, char const* desc, char const* descExtra = NULL );
  /**
   * Set module's version number
   *
   * @param major Major version number
   * @param minor Minor version number
   */
  void setVersion( int major, int minor );
  /**
   * Add line of general module documentation
   * The documentation added here is the full description of this module's purpose, including for example:
   * - Theoretical background information
   * - Additional information necessary to run this module in CSEIS
   * - ...
   *
   * @param text  Line of text to be appended to self-documentation. Use of HTML instructions permitted.
   */
  void addDoc( char const* text );
  /**
   * Add user parameter to module's parameter list
   * An arbitrary number of parameters can be added to each module
   *
   * @param name       Parameter name
   * @param desc       Short description of this parameter's purpose
   * @param parameterType  NUM_VALUES_FIXED for parameters with fixed number of parameter values, NUM_VALUES_VARIABLE otherwise
   * @param descExtra  Additional description of this parameter's purpose
   */
  int addParam( char const* name, char const* desc, int parameterType = NUM_VALUES_FIXED, char const* descExtra = NULL );
  /**
   * Add parameter value to previously added user parameter
   * An arbitrary number of values can be added to each parameter
   *
   * @param defaultVal 
   * @param valueType  VALTYPE_STRING for strings,
   *                   VALTYPE_NUMBER for numbers,
   *                   VALTYPE_OPTION for options,
   *                   VALTYPE_HEADER_NUMBER for alternatively a trace header name, or a number (only works with number headers)
   * @param desc       Short description of this parameter value
   * @param descExtra  Additional description of this parameter value
   */
  int addValue( char const* defaultVal, int valueType, char const* desc = NULL, char const* descExtra = NULL );
  /**
   * Add parameter option to previously added value
   * An arbitrary number of options can be added to each parameter
   *
   * @param name       Option name
   * @param desc       Short description of this parameter option
   * @param descExtra  Additional description of this parameter option
   */
  csParamDescription addOption( char const* name, char const* desc, char const* descExtra = NULL );

  inline csParamDescription const* param( int ip ) const { return myParamList->at(ip)->Descriptor(); }
  inline csParamDescription const* value( int ip, int iv ) const { return myParamList->at(ip)->value(iv)->Descriptor(); }
  inline csParamDescription const* option( int ip, int iv, int io ) const { return myParamList->at(ip)->option(iv,io); }

  /// Return all 'parameters' defined in this Parameter Definition
  bool getParameters( cseis_geolib::csVector<csParamDescription const*>* paramList ) const;
  /// Return all 'values' defined for the parameter (at index ip)
  bool getValues( int ip, cseis_geolib::csVector<csParamDescription const*>* valueList ) const;
  /// Return all 'values' defined for the specified parameter
  bool getValues( char const* paramName, cseis_geolib::csVector<csParamDescription const*>* valueList ) const;
  /// Return all 'options' defined for the specified parameter and value (specified by index ip and iv)
  bool getOptions( int ip, int iv, cseis_geolib::csVector<csParamDescription const*>* optionList ) const;
  /// Return all 'options' defined for the specified parameter and value (specified by parameter name and index iv)
  bool getOptions( char const* paramName, int iv, cseis_geolib::csVector<csParamDescription const*>* optionList ) const;
  /// Retrieve module version number(s)
  void getVersion( int& major, int& minor ) const;
  /// Retrieve module version number(s)
  std::string versionString() const;
  /// Retrieve short self documentation
  std::string docString() const;

  /// Number of parameter defined in this parameter Definition
  inline int numParameters() const { return myParamList->size(); }
  inline int numValues( int ip ) const { return myParamList->at(ip)->numValues(); }
  inline int numOptions( int ip, int iv ) const { return myParamList->at(ip)->value(iv)->numOptions(); }
  csParamDescription const* module() const { return myDescriptorModule; }

  void dump() const;

  void clear();
private:
  csParamDef( csParamDef const& obj );
  //  inline Param const* param( int paramIndex ) const { return myParamList->at(paramIndex); }
  cseis_geolib::csVector<Param*>* myParamList;
  /// Descriptorifier for the module itself
  csParamDescription* myDescriptorModule;
  /// Module version, 2 numbers: 'NUM.NUM'
  int myVersion[2];
  std::string mySelfdoc;

  /// Storage of allocated strings. To be kept until destructor
  cseis_geolib::csVector<char*>* myJunkList;
  char const* createNewString( char const* text, int option = 0 );
  int getParamIndex( char const* name ) const;

  //********************************************************************************
  // Internal classes Param & Value
  //
public:

  /**
   * Parameter from cseis job flow
   * A parameter has a defined name, and a list of 'values'
   * In a job flow, the user has to supply the parameter name, and a list of entries for the parameter 'values'
   */
  class Param {
  public:
    Param( char const* name, char const* desc, char const* descExtra, int parameterType ) :
      myDescriptor( name, desc, descExtra, parameterType ) {
      myValueList = new cseis_geolib::csVector<Value*>(1);
    }
    ~Param() {
      for( int i = 0; i < myValueList->size(); i++ ) {
        delete myValueList->at(i);
      }
      delete myValueList;
      myValueList = NULL;
    }
    int addValue( char const* defaultVal, char const* desc, char const* descExtra, int valueType ) {
      return myValueList->insertEnd( new Value( defaultVal, desc, descExtra, valueType ) );
    }
    csParamDescription const* addOption( char const* name, char const* desc, char const* descExtra ) {
      csParamDescription* Descriptor = myValueList->at(myValueList->size()-1)->addOption( name, desc, descExtra );
      return Descriptor;
    }
    void getOptions( char const* argName, cseis_geolib::csVector<csParamDescription const*>* optionList ) const;
    inline int numValues() const { return myValueList->size(); }

    inline Value const* value( int iv ) const { return myValueList->at(iv); }
    inline csParamDescription const* option( int iv, int io ) const { return myValueList->at(iv)->option(io); }
    inline csParamDescription const* Descriptor() const { return &myDescriptor; }
  private:
    Param( Param const& obj );
    csParamDescription myDescriptor;
    cseis_geolib::csVector<Value*>* myValueList;
  };

  /**
   * Parameter value
   * A parameter value is defined by a 'value type', not by a name (it doesn't have one)
   */
  class Value {
  public:
    Value( char const* defaultVal, char const* desc, char const* descExtra, int valueType ) :
      myDescriptor( defaultVal, desc, descExtra, valueType )
    {
      myOptionList = new cseis_geolib::csVector<csParamDescription*>(1);
    }
    ~Value() {
      for( int i = 0; i < myOptionList->size(); i++ ) {
        delete myOptionList->at(i);
      }
      delete myOptionList;
      myOptionList = NULL;
    }
    csParamDescription* addOption( char const* name, char const* desc, char const* descExtra ) {
      int index = myOptionList->insertEnd( new csParamDescription(name, desc, descExtra, -1) );
      return myOptionList->at(index);
    }
    inline int numOptions() const { return myOptionList->size(); }
//    inline int valueType() const { return myDescriptor.type(); }
    inline csParamDescription const* option( int io ) const { return myOptionList->at(io); }
    inline csParamDescription const* Descriptor() const { return &myDescriptor; }
  private:
    cseis_geolib::csVector<csParamDescription*>* myOptionList;
    csParamDescription myDescriptor;
  };

};

} // namespace

#endif


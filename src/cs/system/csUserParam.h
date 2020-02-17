/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_USER_PARAM_H
#define CS_USER_PARAM_H

#include <string>

namespace cseis_geolib {
  // Forward declarations:
  template <typename T> class csVector;
}

namespace cseis_system {

/**
* CSeis user parameter
* Each cseis module has any number of parameters, which are each stored in one 'osUserParam' object.
* A 'osUserParam' object consists of the user parameter NAME and a list of its VALUES (all kept as strings).
*
* @author Bjorn Olofsson
* @date   2007
*/
class csUserParam {
public:
  csUserParam();
  csUserParam( cseis_geolib::csVector<std::string> const& tokenList );
  csUserParam( csUserParam const& obj );
  csUserParam& operator=( csUserParam const& obj );
  ~csUserParam();
  void set( char const* name, char const* value );
  void set( cseis_geolib::csVector<std::string> const& tokenList );
  void convertToLowerCase( int valIndex );
  void setValue( int valIndex, double value );
//  void setValueTypeVariable( int valIndexStart, int valueType );
  /**
  * @return name of user parameter
  */
  std::string const getName() const;

  /// Return true if lower case equivalent of name equals user parameter name
  bool equals( std::string name ) const;

  /// NOTE: Counts number of accesses
  std::string const getValue( int index = 0 ) const;
  /// NOTE: Counts number of accesses
  void getValues( cseis_geolib::csVector<std::string>* valueList) const;
  int getNumValues() const;
  void dump() const;
private:
  void init();
  void set( csUserParam const& );
  /// Name of user parameter. NOTE: The class methods ensure that the parameter name is a LOWER CASE STRING
  std::string myName;
  /// List of values for this user parameter
  cseis_geolib::csVector<std::string>* myValueList;
};

}

#endif



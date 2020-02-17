/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SELECTION_H
#define CS_SELECTION_H

#include <string>
#include "geolib_defines.h"
#include "csException.h"

namespace cseis_geolib {
  class csFlexNumber;
  class csSelectionField;
  template <typename T> class csVector;

/**
* Selection of values using selection string that follows a special syntax.
*
* A selection string contains:
*  - a list of '/'-separated selections
*  - for each individual header one ':'-separated expression
*  - in each such expression, a list of comma-separated 'selection fields' for each header
*
* Commas ',' mean logical OR, i.e. any of the comma-separated selection criteria that evaluates to true will be used.
* Slashes '/' mean logical OR.
* Colons ':' between individual header selections mean logical AND, i.e. criteria must evaluate to true to all headers
* for a selection to take place.
*
* Syntax of one selection field:
*  SELECTION_SINGLE     100         Number 100
*                       100+2       Number 100 +/-2, i.e. 98,99,100,101,102
*  SELECTION_RANGE      100-200     All numbers from 100 to 200
*                       100-200+2   As above, including each selected number +/-2
*  SELECTION_RANGE_INC  100-200(5)  Every 5th number from 100 to 200, starting with 100
*                       100-200(5)+2 As above, including each selected number +/-2
* Selections with other operators:
*                       <100        All numbers smaller than 100
*                       <=100       All numbers smaller than or equal to 100
*                       >100        All numbers greater than 100
*                       >=100       All numbers greater than or equal to 100
*                       !<any_selection>   All numbers NOT contained in the specified selection
* See @class csSelectionField
*
* @author Bjorn Olofsson
* @date 2005
*/
class csSelection {

public:
/// Selection types
  static int const SELECTION_SINGLE    = 101;
  static int const SELECTION_RANGE     = 103;
  static int const SELECTION_RANGE_INC = 104;
  static int const SELECTION_RANGE_INC_WIDTH = 105;
  static int const SELECTION_OPERATOR  = 106;
  static int const SELECTION_ALL       = 107;

/// Operator types
  static int const OPERATOR_SMALLER = 201;
  static int const OPERATOR_GREATER = 202;
  static int const OPERATOR_SMALLER_EQUAL = 203;
  static int const OPERATOR_GREATER_EQUAL = 204;
  static int const OPERATOR_INCREMENT = 205;
  static int const OPERATOR_RANGE_TO  = 206;

  static int const NONE    = -1;

public:
  csSelection( int numHeaders, type_t const* hdrTypes );
  ~csSelection();
  /**
  * Add selection item
  * Adding two selection strings subsequently is equivalent to first combining both strings into one string with a '/' in between them
  * @param text (i) Selection string
  */
  void add( std::string const& selectionText );
  /**
  * @params values (i) Value of each header in selection
  * @return true if selection contains the specified header values
  */
  bool contains( csFlexNumber const* const values );
  /**
  * Clear selection
  */
  void clear();
  /**
  * Reset headers: Number of headers and header types
  */
  void resetHeaders( int numHeaders, type_t const* hdrTypes );
  void dump();

protected:
  void tokenize( std::string const& text, char separator, csVector<std::string>& fieldTokenList );
  void parse( std::string const& text, csSelectionField* );

private:
  csVector<csSelectionField const*>* mySelectionList;
  csVector<int>* myNumFieldsList;
  int myNumHeaders;
  type_t* myHdrTypes;

  csSelection();
  csSelection( csSelection const& obj );
  csSelection& operator=( const csSelection& obj );
};

class csSelectionException : public csException {
public:
  csSelectionException( std::string const& message, std::string const& expression );
  ~csSelectionException();
private:
  csSelectionException();
};

} // namespace

#endif


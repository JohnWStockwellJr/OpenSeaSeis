/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TOKEN_H
#define CS_TOKEN_H

#include <string>
#include <cstdio>
#include <iostream>

namespace cseis_geolib {

class csMathFunction;

// Value types, can be used for QC output:
static int const STRING = 101;
static int const DOUBLE = 102;
static int const CHAR = 104;
static int const FUNC = 105;
static int const UNKNOWN = -1;

// csToken types:
static int const NUMBER             = 0;
static int const USER_CONSTANT      = 1;
static int const FUNCTION           = 2;
static int const BRACKET_OPEN       = 3;
static int const BRACKET_CLOSE      = 4;
static int const OPERATOR_PLUS_MINUS= 5;
static int const OPERATOR_MULT_DIV  = 6;
static int const COMMA_OPERATOR     = 7;
static int const INTERNAL_VAR       = 8;
static int const OPERATOR_PREFIX_MIN= 9;

// Function types:
static int const MATH_FUNCTION  = 51;


// Basic char types:
static int const OPERATOR = 1;
static int const BRACKET = 2;
static int const DIGIT   = 3;
static int const LETTER  = 4;
static int const COMMA   = 5;
static int const DOT     = 6;

/**
* Expression token
*
* Example usage:
* Mathematical expressions can be tokenized into csTokens, which are operators, numbers, functions etc..
*
* @author Bjorn Olofsson
* @date 2005
*/
class csToken {
public:
  csToken();
  csToken( const csToken& token );
  csToken( double value, int type, int posStart, int posEnd );
  csToken( std::string value, int type, int posStart, int posEnd );
  // User constant:
  csToken( std::string value, int type, int index, int posStart, int posEnd );
  csToken( char value, int type, int posStart );
  csToken( csMathFunction const* value, int posStart, int posEnd );
  csToken& operator=( const csToken& t );
  ~csToken();
  void dump();
  inline bool isNumber() {
    return( type == NUMBER || type == USER_CONSTANT || type == INTERNAL_VAR );
  }
  friend std::ostream& operator<<( std::ostream& out,  csToken const& t );

public:
  csMathFunction* function;
  int type;
  int typeValue;
  std::string valString;
  char   valChar;
  int    valInt;
  double valDouble;
  // The following two variables are useful for providing information to the user
  // for debugging errors in the equation:
  int posStart;  // csToken start position in original expression
  int posEnd;  // csToken end position in original expression

};

std::ostream& operator<<( std::ostream& out,  csToken const& t );

} // namespace

#endif


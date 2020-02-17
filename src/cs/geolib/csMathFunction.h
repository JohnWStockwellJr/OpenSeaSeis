/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_MATH_FUNCTION_H
#define CS_MATH_FUNCTION_H

#include <string>
#include <cmath>
#include <cstdlib>  // For srand() and rand()
#include "geolib_math.h"

namespace cseis_geolib {

typedef double (*Math1ArgPtr) ( double );
typedef double (*Math2ArgPtr) ( double, double );

/**
* Mathematical functions.
*
* Common math functions that can be called in a dynamic way by clients of this class
*
* @author Bjorn Olofsson
* @date 2005
*/
class csMathFunction {
public:
  std::string name;        // Name of function
  std::string description; // Short description of function
  int nArguments;     // Number of function arguments
  int key;            // Unique function key
  Math1ArgPtr method1ArgPtr;
  Math2ArgPtr method2ArgPtr;

  ~csMathFunction() {
  }
  csMathFunction( const csMathFunction& function ) {
    this->key = function.key;
    this->name = std::string(function.name);
    this->nArguments = function.nArguments;
    this->description = std::string(function.description);
    this->method1ArgPtr = function.method1ArgPtr;
    this->method2ArgPtr = function.method2ArgPtr;
    //    printf("Copied new math function: %s %d\n", this->name.c_str(), this->key);
  }

  csMathFunction( int key, std::string name, int nArguments ) {
    this->key = key;
    this->name = std::string(name);
    this->nArguments = nArguments;
    this->description = std::string(name);
    this->method1ArgPtr = NULL;
    this->method2ArgPtr = NULL;
  }
  csMathFunction& operator=( const csMathFunction& function ) {
    this->key = function.key;
    this->name = std::string(function.name);
    this->nArguments = function.nArguments;
    this->description = std::string(function.description);
    this->method1ArgPtr = function.method1ArgPtr;
    this->method2ArgPtr = function.method2ArgPtr;
    //    printf("Equalled new math function: %s %d\n", this->name.c_str(), this->key);
    return *this;
  }

  csMathFunction( int key, std::string name, Math1ArgPtr fPtr ) {
    this->key = key;
    this->name = std::string(name);
    this->nArguments = 1;
    this->description = std::string(name);
    this->method1ArgPtr = fPtr;
    this->method2ArgPtr = NULL;
  }

  csMathFunction( int key, std::string name, Math2ArgPtr fPtr ) {
    this->key = key;
    this->name = std::string(name);
    this->nArguments = 2;
    this->description = std::string(name);
    this->method1ArgPtr = NULL;
    this->method2ArgPtr = fPtr;
  }

  static double INT( double argument ) {
    return double( int(argument) );
  }
  static double ROUND( double argument ) {
    return( (double)( (argument >= 0) ? (int)( argument + 0.55555555 ) : ( (int)( argument - 0.55555555 ) ) ) );
//    return double( int(argument+0.5555555) );
  }
  static double RANDOM( double argument ) {
    return( fmod(rand(),argument) );
  }

  static int const NUM_PRESETS = 27;
  static csMathFunction const* presets() {
    static csMathFunction const presets[NUM_PRESETS] = {
      csMathFunction( 0, "abs",   fabs ),
      csMathFunction( 1, "acos",  acos ),
      csMathFunction( 2, "asin",  asin ),
      csMathFunction( 3, "atan",  atan ),
      csMathFunction( 4, "atan2", atan2 ),
      csMathFunction( 5, "ceil",  ceil ),
      csMathFunction( 6, "cos",   cos ),
      csMathFunction( 7, "cosh",  cosh ),
      csMathFunction( 8, "exp",   exp ),
      csMathFunction( 9, "floor", floor ),
      csMathFunction( 10, "log",   log ),
      csMathFunction( 11, "log10",  log10 ),
      csMathFunction( 12, "max",  cseis_geolib::MAX ),
      csMathFunction( 13, "min",  cseis_geolib::MIN ),
      csMathFunction( 14, "mod",  fmod ),
      csMathFunction( 15, "pow",  pow ),
      csMathFunction( 16, "int",  INT),
      csMathFunction( 17, "round", ROUND ),
      csMathFunction( 18, "sin",  sin ),
      csMathFunction( 19, "sinh",  sinh ),
      csMathFunction( 20, "sqrt", sqrt ),
      csMathFunction( 21, "tan",  tan ),
      csMathFunction( 22, "tanh",  tanh ),
      csMathFunction( 23, "todegrees", RAD2DEG ),
      csMathFunction( 24, "toradians", DEG2RAD ),
      csMathFunction( 25, "sign",  SIGN ),
      csMathFunction( 26, "random",  RANDOM )
    };
    return presets;
  }


private:
  csMathFunction() {
    name = "NONE";
    key = -1;
    nArguments = 0;
    description = "NONE";
    method1ArgPtr = NULL;
    method2ArgPtr = NULL;
  }
};

} // END namespace

#endif


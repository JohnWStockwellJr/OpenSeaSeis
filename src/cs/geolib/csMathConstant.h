/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_MATH_CONSTANT_H
#define CS_MATH_CONSTANT_H

#include <string>
#include <cmath>

namespace cseis_geolib {

/**
 * Mathematical constants
 *
 * @author Bjorn Olofsson
 * @date 2005
 */
class csMathConstant {
public:
  static int const NUM_PRESETS = 2;
  std::string name;
  int key;
  double value;

  csMathConstant( int theKey, std::string theName, double theValue ) {
    this->key   = theKey;
    this->name  = std::string(theName);
    this->value = theValue;
  }
  // Use static function instead of member to avoid static 'initialization fiasco' and other problems
  // g++ optimization created wrong code when using static member
  static csMathConstant const* presets() {
    static csMathConstant const presets[NUM_PRESETS] = {
      csMathConstant( 0, "e", M_E ),
      csMathConstant( 1, "pi", M_PI )
    };
    return presets;
  }
};

} // END namespace

#endif


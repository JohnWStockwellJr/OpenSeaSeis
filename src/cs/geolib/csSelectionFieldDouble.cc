/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <string>
#include <cstdio>
#include <cmath>
#include "csSelection.h"
#include "csSelectionField.h"
#include "csSelectionFieldDouble.h"
#include "csFlexNumber.h"

using namespace cseis_geolib;

csSelectionFieldDouble::csSelectionFieldDouble() : csSelectionField() {
  myDoInvert   = false;
  mySelectType = csSelection::NONE;
  myOperator   = csSelection::NONE;

  myValue      = 0;
  myRangeMin   = 0;
  myRangeMax   = 0;
  myRangeInc   = 0;
  myWidth      = 0;
}
csSelectionFieldDouble::~csSelectionFieldDouble() {
}
void csSelectionFieldDouble::setValue( csFlexNumber const& value ) {
  myValue = value.doubleValue();
}
void csSelectionFieldDouble::setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax ) {
  myRangeMin = rangeMin.doubleValue();
  myRangeMax = rangeMax.doubleValue();
}
void csSelectionFieldDouble::setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax, csFlexNumber const& rangeInc ) {
  setRange( rangeMin, rangeMax );
  myRangeInc = rangeInc.doubleValue();
}
void csSelectionFieldDouble::setWidth( csFlexNumber const& width ) {
  myWidth = width.doubleValue();
}
//--------------------------------------------------------------------------------
bool csSelectionFieldDouble::contains( csFlexNumber const& value_in ) const {
  double value = value_in.doubleValue();
  bool ret;

  switch( mySelectType ) {
  case csSelection::SELECTION_SINGLE:
    ret = ( value == myValue );
    break;
  case csSelection::SELECTION_RANGE:
    ret = ( value >= myRangeMin && value <= myRangeMax );
    break;
  case csSelection::SELECTION_RANGE_INC:
    ret = ( value >= myRangeMin && value <= myRangeMax && fmod( value-myRangeMin, myRangeInc ) == 0 );
    break;
  case csSelection::SELECTION_RANGE_INC_WIDTH:
//    printf("---- %f %f %f %f %f (%f) -----", value, myRangeMin, myRangeMax, myRangeInc, myWidth, fmod(value - myRangeMin + myWidth, myRangeInc));
    ret = ( value >= myRangeMin-myWidth && value <= myRangeMax+myWidth && fmod(value - myRangeMin + myWidth, myRangeInc) <= 2*myWidth );
    break;
  case csSelection::SELECTION_OPERATOR:
    {
      switch( myOperator ) {
      case csSelection::OPERATOR_SMALLER:
        ret = ( value < myValue );
        break;
      case csSelection::OPERATOR_GREATER:
        ret = ( value > myValue );
        break;
      case csSelection::OPERATOR_SMALLER_EQUAL:
        ret = ( value <= myValue );
        break;
      case csSelection::OPERATOR_GREATER_EQUAL:
        ret = ( value >= myValue );
        break;
      default:
        ret = false;
      }  // END inner switch
    }
    break;
  case csSelection::SELECTION_ALL:
    ret = true;
    break;
  default:
    ret = false;
  }  // END switch
//  printf(" Contains (%d,%f,%f,%f,%f) %d %d.. ", mySelectType, myValue, myRangeMin, myRangeMax, myRangeInc, ret, myDoInvert );

  return( !myDoInvert ? ret : !ret );
}

void csSelectionFieldDouble::dump() const {

  if( mySelectType == csSelection::SELECTION_SINGLE ) {
    printf(" Type: %d, operator: %d, value: %f\n", mySelectType, myOperator, myValue );
  }
  else if( mySelectType == csSelection::SELECTION_RANGE ) {
    printf(" Type: %d, operator: %d, min: %f, max: %f\n", mySelectType, myOperator, myRangeMin, myRangeMax );
  }
  else if( mySelectType == csSelection::SELECTION_RANGE_INC ) {
    printf(" Type: %d, operator: %d, min: %f, max: %f, inc: %f\n", mySelectType, myOperator, myRangeMin, myRangeMax, myRangeInc );
  }
  else if( mySelectType == csSelection::SELECTION_RANGE_INC_WIDTH ) {
    printf(" Type: %d, operator: %d, min: %f, max: %f, inc: %f, width: %f\n", mySelectType, myOperator, myRangeMin, myRangeMax, myRangeInc, myWidth );
  }
  else if ( mySelectType == csSelection::SELECTION_ALL ) {
    printf(" Type: %d, ALL\n", mySelectType );
  }
  else {
    printf(" Operator type: %d, operator: %d, value: %f\n", mySelectType, myOperator, myValue );
  }
}


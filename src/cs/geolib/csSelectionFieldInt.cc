/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <string>
#include <cstdio>
#include <cmath>
#include "csSelection.h"
#include "csSelectionFieldInt.h"
#include "csFlexNumber.h"

using namespace cseis_geolib;

csSelectionFieldInt::csSelectionFieldInt() : csSelectionField() {
  myDoInvert   = false;
  mySelectType = csSelection::NONE;
  myOperator   = csSelection::NONE;

  myValue      = 0;
  myRangeMin   = 0;
  myRangeMax   = 0;
  myRangeInc   = 0;
  myWidth      = 0;
}
csSelectionFieldInt::~csSelectionFieldInt() {
}
void csSelectionFieldInt::setValue( csFlexNumber const& value ) {
  myValue = value.intValue();
}
void csSelectionFieldInt::setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax ) {
  myRangeMin = rangeMin.intValue();
  myRangeMax = rangeMax.intValue();
}
void csSelectionFieldInt::setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax, csFlexNumber const& rangeInc ) {
  setRange( rangeMin, rangeMax );
  myRangeInc = rangeInc.intValue();
}
void csSelectionFieldInt::setWidth( csFlexNumber const& width ) {
  myWidth = width.intValue();
}

//--------------------------------------------------------------------------------
bool csSelectionFieldInt::contains( csFlexNumber const& value_in ) const {
  int value = value_in.intValue();
  bool ret;

  switch( mySelectType ) {
  case csSelection::SELECTION_SINGLE:
    ret = ( value == myValue );
    break;
  case csSelection::SELECTION_RANGE:
    ret = ( value >= myRangeMin && value <= myRangeMax );
    break;
  case csSelection::SELECTION_RANGE_INC:
    ret = ( value >= myRangeMin && value <= myRangeMax && ( (value-myRangeMin) % myRangeInc ) == 0 );
    break;
  case csSelection::SELECTION_RANGE_INC_WIDTH:
    ret = ( value >= myRangeMin-myWidth && value <= myRangeMax+myWidth && ((value - myRangeMin + myWidth) % myRangeInc) <= 2*myWidth );
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
  return( !myDoInvert ? ret : !ret );
}

void csSelectionFieldInt::dump() const {
  if( mySelectType == csSelection::SELECTION_SINGLE ) {
    printf(" Type: %d, operator: %d, value: %d\n", mySelectType, myOperator, myValue );
  }
  else if( mySelectType == csSelection::SELECTION_RANGE ) {
    printf(" Type: %d, operator: %d, min: %d, max: %d\n", mySelectType, myOperator, myRangeMin, myRangeMax );
  }
  else if( mySelectType == csSelection::SELECTION_RANGE_INC ) {
    printf(" Type: %d, operator: %d, min: %d, max: %d, inc: %d\n", mySelectType, myOperator, myRangeMin, myRangeMax, myRangeInc );
  }
  else if( mySelectType == csSelection::SELECTION_RANGE_INC_WIDTH ) {
    printf(" Type: %d, operator: %d, min: %d, max: %d, inc: %d, width: %d\n", mySelectType, myOperator, myRangeMin, myRangeMax, myRangeInc, myWidth );
  }
  else if ( mySelectType == csSelection::SELECTION_ALL ) {
    printf(" Type: %d, ALL\n", mySelectType );
  }
  else {
    printf(" Operator type: %d, operator: %d, value: %d\n", mySelectType, myOperator, myValue );
  }
}


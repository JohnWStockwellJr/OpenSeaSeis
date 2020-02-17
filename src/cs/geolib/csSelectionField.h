/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SELECTION_FIELD_H
#define CS_SELECTION_FIELD_H

#include <string>

namespace cseis_geolib {
class csFlexNumber;

/**
* Selection field.
*
* A 'selection field' selects certain values.
* Values can be selected by specifying a single number, a range of numbers, or any other selection
* of numbers indicated by a certain syntax involving operators.
*
* @author Bjorn Olofsson
* @date 2005
*/

class csSelectionField {
public:
  csSelectionField() {}
  virtual ~csSelectionField() {}
  virtual bool contains( csFlexNumber const& value ) const = 0;
  virtual void dump() const = 0;
  virtual void set( int selectType, int theOperator, bool doInvert ) {
    myDoInvert = doInvert;
    mySelectType = selectType;
    myOperator = theOperator;
  }
  virtual void setValue( csFlexNumber const& value ) {}
  virtual void setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax ) {}
  virtual void setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax, csFlexNumber const& rangeInc ) {}
  virtual void setWidth( csFlexNumber const& width ) {}
protected:
  /// false if inverse of selection shall be taken
  bool myDoInvert;
  int  mySelectType;
  int  myOperator;
};

} // namespace

#endif


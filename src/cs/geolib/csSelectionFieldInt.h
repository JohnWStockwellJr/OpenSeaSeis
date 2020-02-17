/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SELECTION_FIELD_INT_H
#define CS_SELECTION_FIELD_INT_H

#include <string>
#include "csSelectionField.h"

namespace cseis_geolib {

/**
* Integer selection field
*
* @author Bjorn Olofsson
*/
class csSelectionFieldInt : public csSelectionField {
public:
  csSelectionFieldInt();
  virtual ~csSelectionFieldInt();
  virtual bool contains( csFlexNumber const& value ) const;
  virtual void dump() const;
  virtual void setValue( csFlexNumber const& value );
  virtual void setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax );
  virtual void setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax, csFlexNumber const& rangeInc );
  virtual void setWidth( csFlexNumber const& width );

private:
  int myValue;
  /// Inclusive minimum of range
  int myRangeMin;
  /// Inclusive maximum of range
  int myRangeMax;
  /// Increment in range
  int myRangeInc;
  /// Width ....
  int myWidth;
};

} // namespace

#endif


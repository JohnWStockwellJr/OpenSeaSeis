/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SELECTION_FIELD_DOUBLE_H
#define CS_SELECTION_FIELD_DOUBLE_H

#include <string>
#include "csSelectionField.h"

namespace cseis_geolib {

/**
* Double selection field
*
* @author Bjorn Olofsson
*/
class csSelectionFieldDouble : public csSelectionField {
public:
  csSelectionFieldDouble();
  virtual ~csSelectionFieldDouble();
  virtual bool contains( csFlexNumber const& value ) const;
  virtual void setValue( csFlexNumber const& value );
  virtual void setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax );
  virtual void setRange( csFlexNumber const& rangeMin, csFlexNumber const& rangeMax, csFlexNumber const& rangeInc );
  virtual void setWidth( csFlexNumber const& width );
  virtual void dump() const;
private:
  double myValue;
  /// Inclusive minimum of range
  double myRangeMin;
  /// Inclusive maximum of range
  double myRangeMax;
  /// Increment in range
  double myRangeInc;
  /// Width ....
  double myWidth;
};

} // namespace

#endif


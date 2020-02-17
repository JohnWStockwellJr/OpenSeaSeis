/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csStack.h"

namespace cseis_geolib {

/**
* Simple stack realisation, for specific use of module indices.
* This class is used by csRunManager.
*
* @author Bjorn Olofsson
* @date   2007
*/
class csModuleIndexStack : public csStack<int> {
public:
  csModuleIndexStack( int none_value );
  csModuleIndexStack( csModuleIndexStack const& obj );
  virtual ~csModuleIndexStack();
  inline virtual int pop();
private:
  int myNoneValue;
};

csModuleIndexStack::csModuleIndexStack( int none_value ) :
  csStack<int>(0), myNoneValue(none_value) {
}
//-------------------------------------------------
csModuleIndexStack::csModuleIndexStack( csModuleIndexStack const& obj ) : 
  csStack<int>( obj ) {
}
//-------------------------------------------------
csModuleIndexStack::~csModuleIndexStack() { 
}
//-------------------------------------------------
inline int csModuleIndexStack::pop() {
  if( csCollection<int>::mySize > 0 ) {
    return csStack<int>::pop();
  }
  else {
    return myNoneValue;
  }
}

}


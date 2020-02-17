/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csUserConstant.h"

using namespace cseis_system;

csUserConstant::csUserConstant() {
  name  = "";
  nameLength = 0;
  value = "";
}
csUserConstant::csUserConstant( std::string const& theName, std::string const& theValue ) {
  name = LETTER_DEFINE + theName + LETTER_DEFINE;
  nameLength = name.length();
  value = theValue;
}
csUserConstant::csUserConstant( csUserConstant const& obj ) {
  name = obj.name;
  nameLength = obj.nameLength;
  value = obj.value;
}
csUserConstant::~csUserConstant() {
}
csUserConstant& csUserConstant::operator=( csUserConstant const& obj ) {
  name = obj.name;
  nameLength = obj.nameLength;
  value = obj.value;
  return *this;
}
bool csUserConstant::operator==( csUserConstant const& obj ) {
  return( !name.compare(obj.name) );
}
std::string csUserConstant::pureName() {
  if( nameLength < 2 ) return name;
//  return string("11111");
  return name.substr( 1, nameLength-2 );
}
void csUserConstant::setValue( std::string theValue ) {
  value = theValue;
}


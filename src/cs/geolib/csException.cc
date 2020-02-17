/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include "csException.h"
#include <stdarg.h>
#include <cstring>

using namespace cseis_geolib;

csException::csException() { myMessage = ""; }
csException::csException( std::string const& message ) { myMessage = message; }
csException::csException( int dummy, char const* text ) {
  myMessage = text;
}
csException::csException( char const* text, ... ) {
  char tmpString[1000];
  va_list argList;
  va_start( argList, text );
  vsprintf( tmpString, text, argList );
  myMessage = tmpString;
}
csException::csException( csException const& obj ) {
  myMessage = obj.myMessage;
}
csException::~csException() {}
const char* csException::getMessage() { return myMessage.c_str(); }


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "csFlexHeader.h"
#include "csException.h"

using namespace cseis_geolib;

const int csFlexHeader::MIN_BYTE_SIZE = 8;

csFlexHeader::csFlexHeader() {
  init();
  myByteSize    = MIN_BYTE_SIZE;
  myValue       = new char[myByteSize];
  myType        = TYPE_INT;
}
void csFlexHeader::init() {  
  myValue    = NULL;
  myByteSize = 0;
}
void csFlexHeader::init( int byteSize ) {  
  myByteSize    = byteSize;
  myValue       = new char[myByteSize];
}
csFlexHeader::csFlexHeader( cseis_geolib::type_t type, std::string value ) {
  myByteSize    = MIN_BYTE_SIZE;
  myValue       = new char[myByteSize];
  myType        = type;
  setValueFromString( type, value );
}
csFlexHeader::csFlexHeader( csFlexHeader const& obj ) {
  init( obj.myByteSize );
  myType     = obj.myType;
  memcpy( myValue, obj.myValue, myByteSize );
}
csFlexHeader::csFlexHeader( double d ) {
  init( MIN_BYTE_SIZE );
//  fprintf(stdout,"Constructor double %x %f\n", myValue, d);
  myType        = TYPE_DOUBLE;
  *((double*)myValue) = d;
}
csFlexHeader::csFlexHeader( float f ) {
  init( MIN_BYTE_SIZE );
//  fprintf(stdout,"Constructor float %x %f\n", myValue, f);
  myType        = TYPE_FLOAT;
  *((float*)myValue) = f;
}
csFlexHeader::csFlexHeader( int i ) {
  init( MIN_BYTE_SIZE );
//  fprintf(stdout,"Constructor int   %x %d\n", myValue, i );
  myType        = TYPE_INT;
  *((int*)myValue) = i;
}
csFlexHeader::csFlexHeader( char c ) {
  init( MIN_BYTE_SIZE );
//  fprintf(stdout,"Constructor int   %x %d\n", myValue, i );
  myType        = TYPE_CHAR;
  myValue[0] = c;
  myValue[1] = '\0';
}
csFlexHeader::csFlexHeader( csInt64_t i ) {
  init( MIN_BYTE_SIZE );
  myType        = TYPE_INT64;
  *((csInt64_t*)myValue) = i;
}
csFlexHeader::csFlexHeader( std::string& value ) {
//  fprintf(stdout,"Constructor int   %x %d\n", myValue, i );
  int newByteSize  = value.length()+1;
  init( std::max( newByteSize, MIN_BYTE_SIZE ) );
  myType = TYPE_STRING;

  memcpy( myValue, value.c_str(), newByteSize-1 );
  myValue[newByteSize-1] = '\0';
}
csFlexHeader::~csFlexHeader() {
//  fprintf(stdout,"Destructor %x\n", myValue);
  myType = TYPE_UNKNOWN;
  if( myValue ) {
    delete [] myValue;
    myValue = NULL;
  }
  myByteSize = 0;
}
//-------------------------------------------------------------
int csFlexHeader::stringSize() {
  if( myType != TYPE_STRING ) {
    return 0;
  }
  else {
    return strlen(myValue);
  }
}

//-------------------------------------------------------------
void csFlexHeader::setValueFromString( cseis_geolib::type_t type, std::string value ) {
  if( type == TYPE_INT ) {
    setIntValue( atoi(value.c_str()) );
  }
  else if( type == TYPE_FLOAT ) {
    setFloatValue( atof(value.c_str()) );
  }
  else if( type == TYPE_DOUBLE ) {
    setDoubleValue( atof(value.c_str()) );
  }
  else if( type == TYPE_INT64 ) {
    setInt64Value( atoi(value.c_str()) );
  }
  else if( type == TYPE_STRING ) {
    setStringValue( value );
  }
  else if( type == TYPE_CHAR ) {
    setCharValue( value[0] );
  }
  else {
    throw( csException("csFlexHeader::setValueFromString: Unknown data type...") );
  }
}
//-------------------------------------------------------------
csFlexHeader& csFlexHeader::operator=( const csFlexHeader& obj ) {
//  fprintf(stdout,"Operator=            %x\n", myValue );
  if( obj.myType == TYPE_STRING && obj.myByteSize > myByteSize ) {
    freeMemory();
    myByteSize = obj.myByteSize;
    myValue    = new char[myByteSize];
  }
  myType = obj.myType;
  memcpy( myValue, obj.myValue, obj.myByteSize );
  return *this;
}
csFlexHeader& csFlexHeader::operator=( int i ) {
//  fprintf(stdout,"Operator=(int)       %x\n", myValue );
  setIntValue( i );
  return *this;
}
csFlexHeader& csFlexHeader::operator=( csInt64_t i ) {
//  fprintf(stdout,"Operator=(int)       %x\n", myValue );
  setInt64Value( i );
  return *this;
}
csFlexHeader& csFlexHeader::operator=( char c ) {
//  fprintf(stdout,"Operator=(int)       %x\n", myValue );
  setCharValue( c );
  return *this;
}
csFlexHeader& csFlexHeader::operator=( double d ) {
//  fprintf(stdout,"Operator=(double)    %x\n", myValue );
  setDoubleValue( d );
  return *this;
}
csFlexHeader& csFlexHeader::operator=( std::string& text ) {
  setStringValue( text );
  return *this;
}
//---------------------------------------------------------------------
//
bool csFlexHeader::operator==( const csFlexHeader& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    if( obj.myType != TYPE_STRING ) {
      return( doubleValue_internal() == obj.doubleValue() );
    }
    break;
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() == obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT || obj.myType == TYPE_FLOAT ) {
      return( floatValue_internal() == obj.floatValue() );
    }
    break;
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() == obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() == obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return( intValue_internal() == obj.intValue() );
    }
    break;
  case TYPE_STRING:
    if( obj.myType == TYPE_STRING ) {
      return( !strcmp( myValue, obj.myValue)  );
    }
    break;
  }
  return false;
}
//---------------------------------------------------------------------
bool csFlexHeader::operator!=( const csFlexHeader& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() != obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() != obj.doubleValue_internal() );
    }
    else {
      return( floatValue_internal() != obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() != obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() != obj.floatValue_internal() );
    }
    else {
      return( intValue_internal() != obj.intValue() );
    }
  case TYPE_STRING:
    if( obj.myType == TYPE_STRING ) {
      return( strcmp( myValue, obj.myValue)  );
    }
    break;
  }
  return false;
}


//---------------------------------------------------------------------
bool csFlexHeader::operator>( const csFlexHeader& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() > obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() > obj.doubleValue_internal() );
    }
    else {
      return( floatValue_internal() > obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() > obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() > obj.floatValue_internal() );
    }
    else {
      return( intValue_internal() > obj.intValue() );
    }
  case TYPE_STRING:
    if( obj.myType == TYPE_STRING ) {
      return( strcmp( myValue, obj.myValue) > 0 );
    }
    break;
  }
  return false;
}

//---------------------------------------------------------------------
bool csFlexHeader::operator<( const csFlexHeader& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() < obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() < obj.doubleValue_internal() );
    }
    else {
      return( floatValue_internal() < obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() < obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() < obj.floatValue_internal() );
    }
    else {
      return( intValue_internal() < obj.intValue() );
    }
  case TYPE_STRING:
    if( obj.myType == TYPE_STRING ) {
      return( strcmp( myValue, obj.myValue) < 0 );
    }
    break;
  }
  return false;
}


//---------------------------------------------------------------------
bool csFlexHeader::operator<=( const csFlexHeader& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() <= obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() <= obj.doubleValue_internal() );
    }
    else {
      return( floatValue_internal() <= obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() <= obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() <= obj.floatValue_internal() );
    }
    else {
      return( intValue_internal() <= obj.intValue() );
    }
  case TYPE_STRING:
    if( obj.myType == TYPE_STRING ) {
      return( strcmp( myValue, obj.myValue) <= 0 );
    }
    break;
  }
  return false;
}


//---------------------------------------------------------------------
bool csFlexHeader::operator>=( const csFlexHeader& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() >= obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() >= obj.doubleValue_internal() );
    }
    else {
      return( floatValue_internal() >= obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() >= obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() >= obj.floatValue_internal() );
    }
    else {
      return( intValue_internal() >= obj.intValue() );
    }
  case TYPE_STRING:
    if( obj.myType == TYPE_STRING ) {
      return( strcmp( myValue, obj.myValue) >= 0 );
    }
    break;
  }
  return false;
}


/*
bool csFlexNumber::operator>( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() > obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() > obj.doubleValue_internal() );
    }
    else {
      return( floatValue_internal() > obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() > obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() > obj.floatValue_internal() );
    }
    else {
      return( intValue_internal() > obj.intValue() );
    }
  }
  return false;
}
*/
//---------------------------------------------------------------------
//

void csFlexHeader::dump() const {

  fprintf(stdout," Type: %2d: ", myType );
  if( myType == TYPE_INT ) {
    fprintf(stdout," valueInt:     %d\n", *((int*)myValue) );
  }
  else if( myType == TYPE_FLOAT ) {
    fprintf(stdout," valueFloat:   %f\n", *((float*)myValue) );
  }
  else if( myType == TYPE_DOUBLE ) {
    fprintf(stdout," valueDouble:  %f\n", *((double*)myValue) );
  }
  else if( myType == TYPE_INT64 ) {
    fprintf(stdout," valueInt64:   %lld\n", *((csInt64_t*)myValue) );
  }
  else if( myType == TYPE_CHAR ) {
    fprintf(stdout," valueChar:   '%c'\n", myValue[0] );
  }
  else { // if( myType == TYPE_STRING ) {
    fprintf(stdout," valueString: '%s'\n", myValue);
  }
}

std::string csFlexHeader::toString() const {
  if( myType == TYPE_STRING || myType == TYPE_CHAR ) {
    std::string text2 = myValue;
    return text2;
  }
  char text[80];
  if( myType == TYPE_INT ) {
    sprintf(text,"%d", *((int*)myValue) );
  }
  else if( myType == TYPE_FLOAT ) {
    sprintf(text,"%f", *((float*)myValue) );
  }
  else if( myType == TYPE_DOUBLE ) {
    sprintf(text,"%f", *((double*)myValue) );
  }
  else if( myType == TYPE_INT64 ) {
    sprintf(text,"%lld", *((csInt64_t*)myValue) );
  }
  else { // 
    sprintf(text," ???");
  }
  return text;
}


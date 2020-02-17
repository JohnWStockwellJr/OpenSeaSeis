/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstring>
#include <string>
#include <cstdio>
#include <cstdlib>
#include "csFlexNumber.h"
#include "csFlexHeader.h"
#include "csException.h"
#include "geolib_string_utils.h"


using namespace cseis_geolib;

csFlexNumber::csFlexNumber() {
  init();
  myType = TYPE_INT;
//  fprintf(stdout,"NEW Constructor %x\n", myValue );
}
void csFlexNumber::init() {
  myValue = new unsigned char[MAX_SIZE];
  memset( myValue, 0, MAX_SIZE );
}
csFlexNumber::csFlexNumber( csFlexNumber const& obj ) {
  init();
//  fprintf(stdout,"COPY Constructor %x\n", myValue );
  myType = obj.myType;
  memcpy( myValue, obj.myValue, MAX_SIZE );
}

csFlexNumber::csFlexNumber( csFlexNumber const* flexNum, bool invertPolarity ) {
  init();
  myType = flexNum->type();
  int sign = invertPolarity ? -1 : 1;
  switch( myType ) {
  case TYPE_INT:
    setIntValue( sign*flexNum->intValue() );
    break;
  case TYPE_INT64:
    setInt64Value( sign*flexNum->int64Value() );
    break;
  case TYPE_FLOAT:
    setFloatValue( (float)sign*flexNum->floatValue() );
    break;
  case TYPE_DOUBLE:
    setDoubleValue( (double)sign*flexNum->doubleValue() );
    break;
  }
}


csFlexNumber::csFlexNumber( csFlexHeader const* flexHdr, bool invertPolarity ) {
  init();
  myType = flexHdr->type();
  int sign = invertPolarity ? -1 : 1;
  switch( myType ) {
  case TYPE_CHAR:
  case TYPE_INT:
    setIntValue( sign*flexHdr->intValue() );
    break;
  case TYPE_INT64:
    setInt64Value( sign*flexHdr->int64Value() );
    break;
  case TYPE_FLOAT:
    setFloatValue( (float)sign*flexHdr->floatValue() );
    break;
  case TYPE_DOUBLE:
    setDoubleValue( sign*flexHdr->doubleValue() );
    break;
  case TYPE_STRING:
    convertToNumber( flexHdr->stringValue() );
    break;
  default:
    convertToNumber( flexHdr->stringValue() );
    break;
  }
}

csFlexNumber::csFlexNumber( double d ) {
  init();
//  fprintf(stdout,"Constructor double %x %f\n", myValue, d);
  myType        = TYPE_DOUBLE;
  *((double*)myValue) = d;
}
csFlexNumber::csFlexNumber( float f ) {
  init();
//  fprintf(stdout,"Constructor float %x %f\n", myValue, f);
  myType        = TYPE_FLOAT;
  *((float*)myValue) = f;
}csFlexNumber::csFlexNumber( int i ) {
  init();
//  fprintf(stdout,"Constructor int   %x %d\n", myValue, i );
  myType        = TYPE_INT;
  *((int*)myValue) = i;
}
csFlexNumber::csFlexNumber( csInt64_t i ) {
  init();
//  fprintf(stdout,"Constructor int   %x %d\n", myValue, i );
  myType        = TYPE_INT64;
  *((csInt64_t*)myValue) = i;
}
csFlexNumber::csFlexNumber( type_t type, std::string value ) {
  init();
  myType = type;
  setValueFromString( type, value );
}
csFlexNumber::~csFlexNumber() {
//  fprintf(stdout,"Destructor %x\n", myValue);
  myType = TYPE_UNKNOWN;
  if( myValue ) {
    delete [] myValue;
    myValue = NULL;
  }
}
//-------------------------------------------------------------
void csFlexNumber::setValueFromString( type_t type, std::string value ) {
  if( type == TYPE_INT ) {
    setIntValue( atoi(value.c_str()) );
  }
  else if( type == TYPE_FLOAT ) {
    setFloatValue( (float)atof(value.c_str()) );
  }
  else if( type == TYPE_DOUBLE ) {
    setDoubleValue( atof(value.c_str()) );
  }
  else {
    throw( csException("csFlexHeader::setValueFromString: Unknown data type...") );
  }
}
//-------------------------------------------------------------
csFlexNumber& csFlexNumber::operator=( csFlexNumber const& obj ) {
  //  fprintf(stdout,"Operator=            %x\n", obj.myValue );
  myType = obj.myType;
  memcpy( myValue, obj.myValue, MAX_SIZE );
  return *this;
}
csFlexNumber& csFlexNumber::operator=( int i ) {
//  fprintf(stdout,"Operator=(int)       %x\n", myValue );
  setIntValue( i );
  return *this;
}
csFlexNumber& csFlexNumber::operator=( csInt64_t i ) {
//  fprintf(stdout,"Operator=(csInt64)       %x\n", myValue );
  setInt64Value( i );
  return *this;
}
csFlexNumber& csFlexNumber::operator=( double d ) {
//  fprintf(stdout,"Operator=(double)    %x\n", myValue );
  setDoubleValue( d );
  return *this;
}
//-------------------------------------------------------------
/*
csFlexNumber::operator double() {
  return doubleValue();
}
csFlexNumber::operator float() {
  return floatValue();
}
csFlexNumber::operator int() {
  return intValue();
}
*/
//-------------------------------------------------------------
//
bool csFlexNumber::convertToNumber( std::string expression ) {
  int counter = 0;
  int length = (int)expression.length();
  while( counter < length && expression[counter] == ' ' ) counter++;
  if( counter == length ) return false;
  myType = TYPE_UNKNOWN;
  int beginIndex = counter;

  // 0. Number starts with a + or - sign
  if( expression[counter] == '-' || expression[counter] == '+' ) {
    counter++;
  }

  // 1. Number starts with a series of digits:
  while( counter < length && isDigit( expression[ counter ] ) ) {
    counter++;
  }
     
  // 2. optional dot for floating point number
  if( counter < length && expression[ counter ] == '.' ) {
    myType = TYPE_DOUBLE;
    counter++;
    while( counter < length && isDigit( expression[ counter ] ) ) {
      counter++;
    }
  }
     
  // 3. Optional 'e' notation
  if( counter < length && expression[ counter ] == 'e' ) {
    myType = TYPE_DOUBLE;
    counter++;
    if( counter >= length || ( expression[ counter ] != '+' && expression[ counter ] != '-' ) ) {
      //throw(  ExpressionException( "Wrong number format (missing sign in 'e' notation)",
      //         csToken( 0.0, NUMBER, beginIndex, counter) ) );
      return false;
    }
    counter++;
    int checkIndex = counter;
    while( counter < length && isDigit( expression[ counter ] ) ) {
      counter++;
    }
    if( counter == checkIndex ) {
      //      throw( ExpressionException( "Wrong number format (no number in 'e' notation)",
      //        csToken(0.0, NUMBER, beginIndex, counter) ) );
      return false;
    }
  }
  else if( myType != TYPE_DOUBLE ) {
    myType = TYPE_INT;
  }

  if( myType == TYPE_DOUBLE ) {
    *((double*)myValue) = atof( expression.substr( beginIndex, counter ).c_str() );
  }
  else {
    *((int*)myValue) = atoi( expression.substr( beginIndex, counter ).c_str() );
  }
  //  fprintf(stdout,"Counter/length %d/%d\n", counter, length);
  return( counter == length );
}
//---------------------------------------------------------------------
//
csFlexNumber csFlexNumber::operator+( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return csFlexNumber( doubleValue_internal() + obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)floatValue_internal() + obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( (double)floatValue_internal() + (double)obj.int64Value_internal() );
    }
    else {
      return csFlexNumber( floatValue_internal() + obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_INT ) {
      return csFlexNumber( intValue_internal() + obj.intValue_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)intValue_internal() + obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( (csInt64_t)intValue_internal() + obj.int64Value_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (float)intValue_internal() + obj.floatValue_internal() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( int64Value_internal() + obj.int64Value_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)int64Value_internal() + obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (double)int64Value_internal() + (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return csFlexNumber( int64Value_internal() + (csInt64_t)obj.intValue_internal() );
    }
  }
  return *this;
}
//---------------------------------------------------------------------
csFlexNumber csFlexNumber::operator-( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return csFlexNumber( doubleValue_internal() - obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)floatValue_internal() - obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( (double)floatValue_internal() - (double)obj.int64Value_internal() );
    }
    else {
      return csFlexNumber( floatValue_internal() - obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_INT ) {
      return csFlexNumber( intValue_internal() - obj.intValue_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)intValue_internal() - obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( (csInt64_t)intValue_internal() - obj.int64Value_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (float)intValue_internal() - obj.floatValue_internal() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( int64Value_internal() - obj.int64Value_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)int64Value_internal() - obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (double)int64Value_internal() - (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return csFlexNumber( int64Value_internal() - (csInt64_t)obj.intValue_internal() );
    }
  }
/*
  switch( myType ) {
  case TYPE_DOUBLE:
    return csFlexNumber( doubleValue_internal() - obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)floatValue_internal() - obj.doubleValue_internal() );
    }
    else {
      return csFlexNumber( floatValue_internal() - obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)intValue_internal() - obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (float)intValue_internal() - obj.floatValue_internal() );
    }
    else {
      return csFlexNumber( intValue_internal() - obj.intValue() );
    }
  }
*/
  return *this;
}
//---------------------------------------------------------------------
csFlexNumber csFlexNumber::operator/( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return csFlexNumber( doubleValue_internal() / obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)floatValue_internal() / obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( (double)floatValue_internal() / (double)obj.int64Value_internal() );
    }
    else {
      return csFlexNumber( floatValue_internal() / obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_INT ) {
      return csFlexNumber( intValue_internal() / obj.intValue_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)intValue_internal() / obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( (csInt64_t)intValue_internal() / obj.int64Value_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (float)intValue_internal() / obj.floatValue_internal() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( int64Value_internal() / obj.int64Value_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)int64Value_internal() / obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (double)int64Value_internal() / (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return csFlexNumber( int64Value_internal() / (csInt64_t)obj.intValue_internal() );
    }
  }
  return *this;
}
//---------------------------------------------------------------------
csFlexNumber csFlexNumber::operator*( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return csFlexNumber( doubleValue_internal() * obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)floatValue_internal() * obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( (double)floatValue_internal() * (double)obj.int64Value_internal() );
    }
    else {
      return csFlexNumber( floatValue_internal() * obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_INT ) {
      return csFlexNumber( intValue_internal() * obj.intValue_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)intValue_internal() * obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( (csInt64_t)intValue_internal() * obj.int64Value_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (float)intValue_internal() * obj.floatValue_internal() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_INT64 ) {
      return csFlexNumber( int64Value_internal() * obj.int64Value_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)int64Value_internal() * obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (double)int64Value_internal() * (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return csFlexNumber( int64Value_internal() * (csInt64_t)obj.intValue_internal() );
    }
  }
  return *this;
}

//---------------------------------------------------------------------
csFlexNumber& csFlexNumber::operator*=( const csFlexNumber& obj ) {
  switch( myType ) {
  case TYPE_DOUBLE:
    setDoubleValue( doubleValue_internal() * obj.doubleValue() );
    break;
  case TYPE_FLOAT:
    setFloatValue( floatValue_internal() * obj.floatValue() );
    break;
  case TYPE_INT:
    setIntValue( intValue_internal() * obj.intValue() );
    break;
  case TYPE_INT64:
    setInt64Value( int64Value_internal() * obj.int64Value() );
    break;
  }
  return *this;
}

//---------------------------------------------------------------------
csFlexNumber csFlexNumber::operator%( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_INT:
    if( obj.myType == TYPE_INT ) {
      return csFlexNumber( intValue_internal() % obj.intValue_internal() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)intValue_internal() - (double)(intValue_internal()/obj.intValue())*obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (float)intValue_internal() - (float)(intValue_internal()/obj.intValue())*obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      throw( csException("csFlexNumber::operator: Type INT64 not supported for operator '%%'") );
    }
  case TYPE_DOUBLE:
    if( obj.myType == TYPE_INT ) {
      return csFlexNumber( doubleValue_internal() - (double)(intValue()/obj.intValue_internal())*obj.doubleValue() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( doubleValue_internal() - (double)(intValue()/obj.intValue())*obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( doubleValue_internal() - (double)(intValue()/obj.intValue())*obj.doubleValue() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      throw( csException("csFlexNumber::operator: Type INT64 not supported for operator '%%'") );
    }
  case TYPE_FLOAT:
    if( obj.myType == TYPE_INT ) {
      return csFlexNumber( (double)floatValue_internal() - (double)(intValue()/obj.intValue_internal())*obj.doubleValue() );
    }
    else if( obj.myType == TYPE_DOUBLE ) {
      return csFlexNumber( (double)floatValue_internal() - (double)(intValue()/obj.intValue())*obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return csFlexNumber( (double)floatValue_internal() - (double)(intValue()/obj.intValue())*obj.doubleValue() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      throw( csException("csFlexNumber::operator: Type INT64 not supported for operator '%%'") );
    }
  case TYPE_INT64:
    throw( csException("csFlexNumber::operator: Type INT64 not supported for operator '%%'") );
  }
  return *this;

}

//---------------------------------------------------------------------
bool csFlexNumber::operator==( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() == obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() == obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return( (double)floatValue_internal() == (double)obj.int64Value_internal() );
    }
    else {
      return( floatValue_internal() == obj.floatValue() );
    }
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
    else {
      return( (csInt64_t)intValue_internal() == obj.int64Value() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)int64Value_internal() == obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (double)int64Value_internal() == (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return( int64Value_internal() == (csInt64_t)obj.intValue_internal() );
    }
    else {
      return( int64Value_internal() == obj.int64Value() );
    }
  }
  return false;
}
//---------------------------------------------------------------------
bool csFlexNumber::operator!=( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() != obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() != obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return( (double)floatValue_internal() != (double)obj.int64Value_internal() );
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
    else if( obj.myType == TYPE_INT ) {
      return( intValue_internal() != obj.intValue() );
    }
    else {
      return( (csInt64_t)intValue_internal() != obj.int64Value() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)int64Value_internal() != obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (double)int64Value_internal() != (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return( int64Value_internal() != (csInt64_t)obj.intValue_internal() );
    }
    else {
      return( int64Value_internal() != obj.int64Value() );
    }
  }
  return false;
}
//---------------------------------------------------------------------
bool csFlexNumber::operator>=( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() >= obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() >= obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return( (double)floatValue_internal() >= (double)obj.int64Value_internal() );
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
    else if( obj.myType == TYPE_INT ) {
      return( intValue_internal() >= obj.intValue() );
    }
    else {
      return( (csInt64_t)intValue_internal() >= obj.int64Value() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)int64Value_internal() >= obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (double)int64Value_internal() >= (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return( int64Value_internal() >= (csInt64_t)obj.intValue_internal() );
    }
    else {
      return( int64Value_internal() >= obj.int64Value() );
    }
  }
  return false;
}
//---------------------------------------------------------------------
bool csFlexNumber::operator<=( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() <= obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() <= obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return( (double)floatValue_internal() <= (double)obj.int64Value_internal() );
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
    else if( obj.myType == TYPE_INT ) {
      return( intValue_internal() <= obj.intValue() );
    }
    else {
      return( (csInt64_t)intValue_internal() <= obj.int64Value() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)int64Value_internal() <= obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (double)int64Value_internal() <= (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return( int64Value_internal() <= (csInt64_t)obj.intValue_internal() );
    }
    else {
      return( int64Value_internal() <= obj.int64Value() );
    }
  }
  return false;
}
//---------------------------------------------------------------------
bool csFlexNumber::operator>( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() > obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() > obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return( (double)floatValue_internal() > (double)obj.int64Value_internal() );
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
    else if( obj.myType == TYPE_INT ) {
      return( intValue_internal() > obj.intValue() );
    }
    else {
      return( (csInt64_t)intValue_internal() > obj.int64Value() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)int64Value_internal() > obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (double)int64Value_internal() > (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return( int64Value_internal() > (csInt64_t)obj.intValue_internal() );
    }
    else {
      return( int64Value_internal() > obj.int64Value() );
    }
  }
  return false;
}
//---------------------------------------------------------------------
bool csFlexNumber::operator<( const csFlexNumber& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() < obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() < obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT64 ) {
      return( (double)floatValue_internal() < (double)obj.int64Value_internal() );
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
    else if( obj.myType == TYPE_INT ) {
      return( intValue_internal() < obj.intValue() );
    }
    else {
      return( (csInt64_t)intValue_internal() < obj.int64Value() );
    }
  case TYPE_INT64:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)int64Value_internal() < obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (double)int64Value_internal() < (double)obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return( int64Value_internal() < (csInt64_t)obj.intValue_internal() );
    }
    else {
      return( int64Value_internal() < obj.int64Value() );
    }
  }
  return false;
}
//---------------------------------------------------------------------
//
void csFlexNumber::dump() const {

  printf(" Type: %d: ", myType );
  if( myType == TYPE_INT ) {
    printf(" valueInt: %d\n", *((int*)myValue) );
  }
  else if( myType == TYPE_FLOAT ) {
    printf(" valueFloat: %f\n", *((float*)myValue) );
  }
  else if( myType == TYPE_DOUBLE ) {
    printf(" valueDouble: %f\n", *((double*)myValue) );
  }
  else if( myType == TYPE_INT64 ) {
    printf(" valueInt: %lld\n", *((csInt64_t*)myValue) );
  }
}

csFlexNumberException::csFlexNumberException( std::string expression ) {
  myMessage = std::string("Expression does not evaluate into number: '") + expression + std::string("':\n");
}
csFlexNumberException::~csFlexNumberException() {
}


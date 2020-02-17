/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <string>
#include <cstdio>
#include <iostream>
#include "csToken.h"
#include "csMathFunction.h"

using namespace cseis_geolib;

csToken::csToken() :
    function( NULL ),
    type( UNKNOWN ),
    typeValue( UNKNOWN ),
    valString(""),
    valChar(' '),
    valInt( 0 ),
    valDouble( 0.0 ),
    posStart( 0 ),
    posEnd( 0 )
{
}

csToken::csToken( const csToken& token ) :
    function( NULL ),
    type( token.type ),
    typeValue( token.typeValue ),
    valString( token.valString ),
    valChar( token.valChar ),
    valInt( token.valInt ),
    valDouble( token.valDouble ),
    posStart( token.posStart ),
    posEnd( token.posEnd )
{
  if( token.function != NULL ) {
    function = new csMathFunction( *(token.function) );
  }
}

csToken::csToken( double value, int type, int posStart, int posEnd ) :
  function( NULL ),
  type( type ),
  typeValue( DOUBLE ),
  valString(""),
  valChar(' '),
  valInt( 0 ),
  valDouble( value ),
  posStart( posStart ),
  posEnd( posEnd )
{
}

csToken::csToken( std::string value, int type, int posStart, int posEnd ) :
  function( NULL ),
  type( type ),
  typeValue( STRING ),
  valString( value ),
  valChar(' '),
  valInt( 0 ),
  valDouble( 0.0 ),
  posStart( posStart ),
  posEnd( posEnd )
{
}

// User constant
csToken::csToken( std::string value, int type, int index, int posStart, int posEnd ) :
  function( NULL ),
  type( type ),
  typeValue( STRING ),
  valString( value ),
  valChar(' '),
  valInt( index ),
  valDouble( 0.0 ),
  posStart( posStart ),
  posEnd( posEnd )
{
}

csToken::csToken( char value, int type, int posStart ) :
  function( NULL ),
  type( type ),
  typeValue( CHAR ),
  valString(""),
  valChar( value ),
  valInt( 0 ),
  valDouble( 0.0 ),
  posStart( posStart ),
  posEnd( posStart+1 )
{
}

csToken::csToken( csMathFunction const* value, int posStart, int posEnd ) :
  function( NULL ),
  type( FUNCTION ),
  typeValue( FUNC ),
  valString(""),
  valChar(' '),
  valInt( 0 ),
  valDouble( 0.0 ),
  posStart( posStart ),
  posEnd( posEnd )
{
  if( value != NULL ) {
    function = new csMathFunction( *value );
  }
}

csToken& csToken::operator=( const csToken& t ) {
  if( this->function != NULL ) {
    delete this->function;
    this->function = NULL;
  }
  if( t.function != NULL ) {
    this->function = new csMathFunction( *(t.function) );
  }
  this->type      = t.type;
  this->typeValue = t.typeValue;
  this->valString = t.valString;
  this->valChar   = t.valChar;
  this->valDouble = t.valDouble;
  this->valInt   = t.valInt;
  this->posStart = t.posStart;
  this->posEnd   = t.posEnd;
  return *this;
}

csToken::~csToken() {
  if( function != NULL ) {
    delete function;
    function = NULL;
  }
}

std::ostream& cseis_geolib::operator<<( std::ostream& out,  csToken const& t ) {
//std::ostream& operator<<( std::ostream& out,  const csToken& t ) {
  switch( t.type ) {
  case NUMBER:
    out << t.valDouble;
    break;
  case OPERATOR_PLUS_MINUS:
  case OPERATOR_MULT_DIV:
  case COMMA_OPERATOR:
  case BRACKET_OPEN:
  case BRACKET_CLOSE:
    out << t.valChar;
    break;
  case OPERATOR_PREFIX_MIN:
    out << "prefix-";
    break;
  case FUNCTION:
    out << t.function->name;
    break;
  case USER_CONSTANT:
    out << t.valString;
    break;
  case INTERNAL_VAR:
    out << "STK#" << t.valInt;
    break;
  default:
    out << "???";
  }
  return out;
}

void csToken::dump() {
  if( typeValue == DOUBLE ) {
    printf( "Double: %f\n", valDouble );
  }
  else if( typeValue == STRING ) {
    printf( "String: %s\n", valString.c_str() );
  }
  else if( typeValue == CHAR ) {
    printf( "Char:   '%c'\n", valChar );
  }
  else if( typeValue == FUNC ) {
    printf( "Func:   '%s'  %d\n", function->name.c_str(), function->key );
  }
  else {
    printf( "????:   %f    %s    %c\n", valDouble, valString.c_str(), valChar );
  }
}



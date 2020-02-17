/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include "csUserParam.h"
#include "csParamDef.h"
// From geolib
#include "csVector.h"
#include "csException.h"
#include "geolib_string_utils.h"

using namespace cseis_system;

csUserParam::csUserParam() {
  init();
}
void csUserParam::init() {
  myName = "_no_name_";
  myValueList = new cseis_geolib::csVector<std::string>(3);
}
csUserParam::csUserParam( cseis_geolib::csVector<std::string> const& tokenList ) {
  myValueList = new cseis_geolib::csVector<std::string>(3);
  set( tokenList );
}
csUserParam::~csUserParam() {
  if( myValueList != NULL ) {
    delete myValueList;
    myValueList = NULL;
  }
  //    fprintf(stdout,"UserParam: Destruct %x\n", this);
  //    fflush(stdout);
}

csUserParam::csUserParam( csUserParam const& obj ) {
  myValueList = new cseis_geolib::csVector<std::string>(1);
  set( obj );
}
void csUserParam::set( csUserParam const& obj ) {
  myName = obj.myName;
  //  myValueList = obj.myValueList;
  for( int i = 0; i < obj.myValueList->size(); i++ ) {
    myValueList->insertEnd( obj.myValueList->at(i) );
  }
}
csUserParam& csUserParam::operator=( csUserParam const& obj ) {
  myValueList->clear();
  set( obj );
  return *this;
}

bool csUserParam::equals( std::string name ) const {
  return( !myName.compare( cseis_geolib::toLowerCase(name) ) );
}

void csUserParam::set( char const* name, char const* value ) {
  myName = cseis_geolib::toLowerCase(name);
  myValueList->clear();
  myValueList->insertEnd( std::string(value) );
}

void csUserParam::set( cseis_geolib::csVector<std::string> const& tokenList ) {
  myName = cseis_geolib::toLowerCase(tokenList.at(0));
  myValueList->clear();
  for( int i = 1; i < tokenList.size(); i++ ) {
    myValueList->insertEnd( tokenList.at(i) );
  }
}

//  const char* getName() const { return myName.c_str(); }
//  const char* getValue() const { return myValueList.at(0).c_str(); }
std::string const csUserParam::getName() const {
  return myName;
}
//--------------------------------------------------------------------------------
//
std::string const csUserParam::getValue( int index ) const {
  if( myValueList->size() > index ) {
    return myValueList->at(index);
  }
  else {
    throw( cseis_geolib::csException("csUserParam::getValue: Wrong index provided (=%d) for user parameter '%s'. Number of elements = %d. This is a program bug in the calling function.",
      index, myName.c_str(), myValueList->size()) );
//    return std::string("_none_");
  }
}
//--------------------------------------------------------------------------------
//
void csUserParam::getValues( cseis_geolib::csVector<std::string>* valueList) const {
  valueList->clear();
  for( int i = 0; i < myValueList->size(); i++ ) {
    valueList->insertEnd( myValueList->at(i) );
  }
}
int csUserParam::getNumValues() const {
  return myValueList->size();
}
void csUserParam::csUserParam::dump() const {
  std::cout << "*** PARAMETER " << myName << ", numValues: " << myValueList->size() << " ***" << std::endl;
  myValueList->dump();
  fflush(stdout);
}
//--------------------------------------------------------------------------------
//
void csUserParam::convertToLowerCase( int valIndex ) {
  myValueList->set( cseis_geolib::toLowerCase( myValueList->at(valIndex) ), valIndex );
}
void csUserParam::setValue( int valIndex, double value ) {
  if( valIndex < getNumValues() ) {
    double diff = fabs( value - (double)(int)value );
    std::ostringstream sstream;
    if( value == 0 || diff < 1e-7 ) {
      sstream << (int)value;
    }
    else {
      sstream << value;
    }
    std::string valueAsText = sstream.str();
    myValueList->set( valueAsText, valIndex );
  }
}
/*
void csUserParam::setValueTypeVariable( int valIndexStart, int valueType ) {
  if( valueType == VALTYPE_OPTION ) {
    int nValues = getNumValues();
    for( int i = valIndexStart; i < nValues; i++ ) {
      myValueList->set( toLowerCase( myValueList->at(i) ), i );
    }
  }
}
*/


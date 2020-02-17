/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csParamManager.h"
#include "csUserParam.h"
#include "csVector.h"
#include "csLogWriter.h"
#include <cstring>

using namespace cseis_system;

csParamManager::csParamManager( cseis_geolib::csVector<csUserParam*> const* userParamList, csLogWriter* log ) {
  myUserParamList = userParamList;
  myLog = log;
  myNumValueCalls = NULL;
  if( myUserParamList->size() > 0 ) {
    myNumValueCalls = new int[myUserParamList->size()];
    for( int i = 0; i < myUserParamList->size(); i++ ) {
      myNumValueCalls[i] = 0;
    }
  }
  myIndexCurrentParam = NOT_FOUND;
  myLineIndexCurrentParam = NOT_FOUND;
  myNameCurrentParam = NULL;
}
csParamManager::~csParamManager() {
  if( myNumValueCalls != NULL ) {
    delete [] myNumValueCalls;
    myNumValueCalls = NULL;
  }
}
int csParamManager::numParameters() const {
  return myUserParamList->size();
}
int csParamManager::getNumValueCalls( int index ) const {
  return myNumValueCalls[index];
}

int csParamManager::getNumLines( char const* name ) {
  int counter = 0;
  std::string idName = std::string(name);
  for( int i = 0; i < myUserParamList->size(); i++ ) {
    if( myUserParamList->at(i)->equals(idName) ) {
      counter++;
    }
  }
  return counter;
}
//---------------------------------------------------------------------------
int csParamManager::getNumValues( char const* name, int indexLine ) {
  int paramIndex = getParamIndex( name, indexLine );
  if( paramIndex != NOT_FOUND ) {
//   !CHANGE! CHECK csParamDef instead
    return myUserParamList->at(paramIndex)->getNumValues();
  }
  else {
    return 0;
  }
}
//---------------------------------------------------------------------------
//
void csParamManager::getIntAtLine( char const* name, int* value, int indexLine, int indexValue ) {
  if( myLineIndexCurrentParam == indexLine && myNameCurrentParam != NULL && !strcmp( myNameCurrentParam, name ) ) {
    *value = atoi( myUserParamList->at(myLineIndexCurrentParam)->getValue( indexValue ).c_str() );
  }
  else {
    std::string const valueStr = getParamValue( name, indexLine, indexValue );
    *value = atoi( valueStr.c_str() );
  }
}
void csParamManager::getFloatAtLine( char const* name, float* value, int indexLine, int indexValue ) {
  double valueDouble;
  getDoubleAtLine( name, &valueDouble, indexLine, indexValue );
  *value = (float)valueDouble;
}
void csParamManager::getDoubleAtLine( char const* name, double* value, int indexLine, int indexValue ) { 
 if( myLineIndexCurrentParam == indexLine && myNameCurrentParam != NULL && !strcmp( myNameCurrentParam, name ) ) {
    *value = atof( myUserParamList->at(myLineIndexCurrentParam)->getValue( indexValue ).c_str() );
  }
  else {
    std::string const valueStr = getParamValue( name, indexLine, indexValue );
    *value = atof( valueStr.c_str() );
  }
}
//---------------------------------------------------------------------------
void csParamManager::getText( char const* name, char** value, int indexLine, int indexValue ) {
  std::string const valueStr = getParamValue( name, indexLine, indexValue );
  int length = valueStr.length();
  *value = new char[length+1];
  memcpy( *value, valueStr.c_str(), length );
  (*value)[length] = '\0';
}

void csParamManager::getStringAtLine( char const* name, std::string* value, int indexLine, int indexValue ) {
  if( myLineIndexCurrentParam == indexLine && myNameCurrentParam != NULL && !strcmp( myNameCurrentParam, name ) ) {
    *value = myUserParamList->at(myLineIndexCurrentParam)->getValue( indexValue );
  }
  else {
    *value = getParamValue( name, indexLine, indexValue );
  }
}
void csParamManager::getAll( char const* name, cseis_geolib::csVector<std::string>* valueList, int indexLine ) {
  getParamValues( name, indexLine, valueList );
}
bool csParamManager::exists( char const* name ) {
  std::string idName = std::string(name);
  for( int i = 0; i < myUserParamList->size(); i++ ) {
    if( myUserParamList->at(i)->equals(idName) ) {
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------
// Internal method
int csParamManager::getParamIndex( char const* name, int indexLine, bool exOnFail ) {
  int counter = 0;
  std::string idName = std::string(name);
  for( int i = 0; i < myUserParamList->size(); i++ ) {
    if( myUserParamList->at(i)->equals(idName) ) {
      if( counter == indexLine ) {
        myNumValueCalls[i] += 1;
        return i;
      }
      counter++;
    }
  }
  if( exOnFail ) myLog->error("User parameter '%s' not specified.\n", name );
  return NOT_FOUND;
}

std::string const csParamManager::getParamValue( char const* name, int indexLine, int indexValue, bool exOnFail ) {
  int counter = 0;
  std::string idName = std::string(name);
  for( int i = 0; i < myUserParamList->size(); i++ ) {
    if( myUserParamList->at(i)->equals(idName) ) {
      if( counter == indexLine ) {
        myNumValueCalls[i] += 1;
        return myUserParamList->at(i)->getValue( indexValue );
      }
      counter++;
    }
  }
  if( exOnFail ) myLog->error("User parameter '%s' not specified.\n", name );
  return std::string("");
}
bool csParamManager::getParamValues( char const* name, int indexLine, cseis_geolib::csVector<std::string>* valueList, bool exOnFail ) {
  int counter = 0;
  std::string idName = std::string(name);
  for( int i = 0; i < myUserParamList->size(); i++ ) {
    //      if( !strcmp(myUserParamList[i].getName(),name) ) {
    if( myUserParamList->at(i)->equals(idName) ) {
      if( counter == indexLine ) {
        myNumValueCalls[i] += 1;
        myUserParamList->at(i)->getValues( valueList );
        return true;
      }
      counter++;
    }
  }
  if( exOnFail ) myLog->error("User parameter '%s' not specified.\n", name );
  return false;
}



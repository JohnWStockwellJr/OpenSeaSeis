/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csParamDef.h"
#include "cseis_defines.h"
#include "csException.h"
#include "geolib_math.h"
#include <cstring>

using namespace cseis_system;

csParamDef::csParamDef() : myDescriptorModule( NULL ) {
  myParamList = new cseis_geolib::csVector<Param*>(1);
  myVersion[MAJOR]  = 1;
  myVersion[MINOR]  = 0;
  mySelfdoc         = "";
  myJunkList = new cseis_geolib::csVector<char*>;
}

csParamDef::~csParamDef() {
  if( myParamList != NULL ) {
    for( int i = 0; i < myParamList->size(); i++ ) {
      delete myParamList->at(i);
    }
    delete myParamList;
    myParamList = NULL;
  }
  if( myDescriptorModule != NULL ) {
    delete myDescriptorModule;
    myDescriptorModule = NULL;
  }
  if( myJunkList != NULL ) {
    for( int i = 0; i < myJunkList->size(); i++ ) {
      delete [] myJunkList->at(i);
    }
    delete myJunkList;
    myJunkList = NULL;
  }
}

void csParamDef::setModule( char const* moduleName, char const* desc, char const* descExtra ) {
  if( myDescriptorModule ) delete myDescriptorModule;
  myDescriptorModule = new csParamDescription( createNewString(moduleName,UPPER_CASE), desc, descExtra, -1 );
}
void csParamDef::setVersion( int major, int minor ) {
  myVersion[MAJOR]  = major;
  myVersion[MINOR]  = minor;
}
void csParamDef::addDoc( char const* text ) {
  mySelfdoc.append( text );
  mySelfdoc.append( "\n" );
}
/// Add parameter
int csParamDef::addParam( char const* name, char const* desc, int parameterType, char const* descExtra ) {
  int index = myParamList->insertEnd( new Param( createNewString(name,LOWER_CASE), desc, descExtra, parameterType) );
  return index;
}
/// Add value to previously added parameter
int csParamDef::addValue( char const* defaultVal, int valueType, char const* desc, char const* descExtra ) {
  if( myParamList->size() == 0 ) {
    throw( cseis_geolib::csException("csParamDef::addValue: Can't add value before adding the associated parameter") );
  }
  if( valueType != VALTYPE_NUMBER && valueType != VALTYPE_STRING && valueType != VALTYPE_OPTION && valueType != VALTYPE_HEADER_NUMBER ) {
    throw( cseis_geolib::csException("csParamDef::addValue: Unknown value type: %d", valueType) );
  }
  return myParamList->at(myParamList->size()-1)->addValue( defaultVal, desc, descExtra, valueType );
}
// Add option to previously added value
csParamDescription csParamDef::addOption( char const* name, char const* desc, char const* descExtra ) {
  if( myParamList->size() == 0 ) {
    throw( cseis_geolib::csException("csParamDef::addOption: Can't add option before adding the associated parameter and value") );
  }
  else if( numValues(myParamList->size()-1) <= 0 ) {
    throw( cseis_geolib::csException("csParamDef::addOption: Can't add option before adding the associated value") );
  }

  return *myParamList->at(myParamList->size()-1)->addOption( createNewString(name,LOWER_CASE), desc, descExtra );
}

//--------------------------------------------------------------------------------
//
bool csParamDef::getParameters( cseis_geolib::csVector<csParamDescription const*>* paramList ) const {
  for( int i = 0; i < myParamList->size(); i++ ) {
    paramList->insertEnd( myParamList->at(i)->Descriptor() );
  }
  return true;
}
bool csParamDef::getValues( int ip, cseis_geolib::csVector<csParamDescription const*>* argList ) const {
  Param* param = myParamList->at(ip);
  int nValues = param->numValues();
  for( int iv = 0; iv < nValues; iv++ ) {
    argList->insertEnd( param->value(iv)->Descriptor() );
  }
  return true;
}
bool csParamDef::getValues( char const* paramName, cseis_geolib::csVector<csParamDescription const*>* argList ) const {
  int paramIndex = getParamIndex( paramName );
  if( paramIndex < 0 ) return false;
  return getValues( paramIndex, argList );
}
std::string csParamDef::docString() const {
  return mySelfdoc;
}
//--------------------------------------------------------------------------------
//
int csParamDef::getParamIndex( char const* name ) const {
  int paramIndex = -1;
  for( int i = 0; i < myParamList->size(); i++ ) {
    if( !strcmp( name, myParamList->at(i)->Descriptor()->name() ) ) {
      paramIndex = i;
      break;
    }
  }
  return paramIndex;
}
//--------------------------------------------------------------------------------
//
bool csParamDef::getOptions( int ip, int iv, cseis_geolib::csVector<csParamDescription const*>* optionList ) const {
  Value const* value = myParamList->at(ip)->value(iv);
  int nOptions = value->numOptions();
  for( int io = 0; io < nOptions; io++ ) {
    optionList->insertEnd( value->option(io) );
  }
  return true;
}
bool csParamDef::getOptions( char const* paramName, int iv, cseis_geolib::csVector<csParamDescription const*>* optionList ) const {
  int paramIndex = getParamIndex( paramName );
  if( paramIndex < 0 ) return false;
  return getOptions( paramIndex, iv, optionList );
}

void csParamDef::clear() {
  if( myParamList != NULL ) {
    for( int i = 0; i < myParamList->size(); i++ ) {
      delete myParamList->at(i);
    }
    myParamList->clear();
  }
  if( myDescriptorModule != NULL ) {
    delete myDescriptorModule;
    myDescriptorModule = NULL;
  }
}
//--------------------------------------------------------------------------------
//
//
char const* csParamDef::createNewString( char const* text, int option ) {
  int length = strlen(text);
  char* textNew = new char[length+1];
  memcpy( textNew, text, length );
  textNew[length] = '\0';
  if( option == LOWER_CASE ) cseis_geolib::toLowerCase( textNew );
  else if( option == UPPER_CASE ) cseis_geolib::toUpperCase( textNew );
  myJunkList->insertEnd( textNew );
  return textNew;
}
//--------------------------------------------------------------------------------
void csParamDef::getVersion( int& major, int& minor ) const {
  major = myVersion[MAJOR];
  minor = myVersion[MINOR];
}
std::string csParamDef::versionString() const {
  char ver[40];
  sprintf( ver, "%d.%d%c", myVersion[MAJOR], myVersion[MINOR], '\0' );
  return( std::string(ver) );
}
//--------------------------------------------------------------------------------
// Dump methods
//
void csParamDef::dump() const {
  csParamDescription const* valueDescriptor;
  csParamDescription const* optionDescriptor;

  std::fprintf(stdout,"\n#****************************************************\n");

  if( myDescriptorModule != NULL ) {
    fprintf(stdout,"  Module: '%s': %s\n", myDescriptorModule->name(), myDescriptorModule->desc());
    fprintf(stdout,"   Extra: %s\n", myDescriptorModule->descExtra());
  }
  else {
    fprintf(stdout,"  ---Module name not defined--- \n");
  }
  int nParams = numParameters();
  for( int ip = 0; ip < nParams; ip++ ) {
    fprintf(stdout,"  Parameter #%d: ", ip);
    fprintf(stdout, " Name: %s, desc: %s, descExtra: %s, type: %d\n",
            param(ip)->name(), param(ip)->desc(), param(ip)->descExtra(), param(ip)->type() );
    
    int nValues = numValues(ip);
    for( int iv = 0; iv < nValues; iv++ ) {
      fprintf(stdout,"   Value #%d: ", iv);
      valueDescriptor = value(ip,iv);
      fprintf(stdout, "  Name: %s, desc: %s, descExtra: %s, type: %d\n",
              valueDescriptor->name(), valueDescriptor->desc(), valueDescriptor->descExtra(), valueDescriptor->type() );
      int nOptions = numOptions(ip,iv);
      for( int io = 0; io < nOptions; io++ ) {
        fprintf(stdout,"     Option #%d: ", io);
        optionDescriptor = option(ip,iv,io);
        fprintf(stdout, "  Name: %s, desc: %s, descExtra: %s, type: %d\n",
                optionDescriptor->name(), optionDescriptor->desc(), optionDescriptor->descExtra(), optionDescriptor->type() );
      }
    }
  }
}

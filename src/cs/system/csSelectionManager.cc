/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSelectionManager.h"
#include "geolib_defines.h"
#include "csFlexNumber.h"
#include "csFlexHeader.h"
#include "csVector.h"
#include "csSelection.h"
#include "csException.h"
#include "csTraceHeaderDef.h"
#include "csTraceHeader.h"
#include "csException.h"
#include <string>

using namespace cseis_system;

//--------------------------------------------------------------
//
csSelectionManager::csSelectionManager() {
  myNumHeaders  = 0;
  myHeaderIndex = NULL;
  myHeaderType  = NULL;
  myHeaderNames = NULL;
  myValues      = NULL;
  mySelection   = NULL;
  myStringSelectionText = "";
}
//--------------------------------------------------------------
//
csSelectionManager::~csSelectionManager() {
  if( myHeaderIndex ) { delete [] myHeaderIndex; myHeaderIndex = NULL; }
  if( myHeaderType ) { delete [] myHeaderType; myHeaderType = NULL; }
  if( myHeaderNames ) { delete [] myHeaderNames; myHeaderNames = NULL; }
  if( myValues ) { delete [] myValues; myValues = NULL; }
  if( mySelection ) { delete mySelection; mySelection = NULL; }
  myNumHeaders   = 0;
}
//--------------------------------------------------------------
//
void csSelectionManager::set( cseis_geolib::csVector<std::string> const* headerList, std::string const* selectionText,
    cseis_system::csTraceHeaderDef const* hdef ) {
  if( headerList->size() == 0 ) {
    throw cseis_geolib::csException( "Empty trace header list" );
  }
  myNumHeaders   = headerList->size();
  myHeaderIndex  = new int[myNumHeaders];
  myHeaderType   = new cseis_geolib::type_t[myNumHeaders];
  myHeaderNames  = new std::string[myNumHeaders];
  myValues       = new cseis_geolib::csFlexNumber[myNumHeaders];
  for( int k = 0; k < myNumHeaders; k++ ) {
    std::string name = headerList->at(k);
    myHeaderNames[k] = name;
    if( hdef->headerExists( name.c_str() ) ) {
      myHeaderIndex[k] = hdef->headerIndex( name.c_str() );
      myHeaderType[k]  = hdef->headerType( name.c_str() );
      if( myHeaderType[k] != cseis_geolib::TYPE_INT && myHeaderType[k] != cseis_geolib::TYPE_FLOAT && myHeaderType[k] != cseis_geolib::TYPE_DOUBLE ) {
        if( myNumHeaders == 1 && myHeaderType[k] == cseis_geolib::TYPE_STRING ) {
          myStringSelectionText= *selectionText;
          return;
        }
        throw cseis_geolib::csException( "Trace header '%s': Selection using multiple headers is only supported for number headers.", name.c_str() );
      }
    }
    else {
      throw cseis_geolib::csException( "Unknown trace header: %s", name.c_str() );
    }
  }
  mySelection = new cseis_geolib::csSelection( myNumHeaders, myHeaderType );
  mySelection->add( *selectionText );
}
void csSelectionManager::dump() {
  mySelection->dump();
}
//----------------------------------------------------------
//
bool csSelectionManager::contains( cseis_system::csTraceHeader const* trcHeader ) {
  for( int i = 0; i < myNumHeaders; i++ ) {
    if( myHeaderType[i] == cseis_geolib::TYPE_FLOAT ) {
      myValues[i].setDoubleValue( trcHeader->floatValue(myHeaderIndex[i]) );
    }
    else if( myHeaderType[i] == cseis_geolib::TYPE_DOUBLE ) {
      myValues[i].setDoubleValue( trcHeader->doubleValue(myHeaderIndex[i]) );
    }
    else if( myHeaderType[i] == cseis_geolib::TYPE_INT ) {
      myValues[i].setIntValue( trcHeader->intValue( myHeaderIndex[i] ) );
    }
    else { // TYPE_STRING
      return( strcmp( myStringSelectionText.c_str(), trcHeader->stringValue(myHeaderIndex[i]).c_str() ) == 0 );
    }
  }
  if( mySelection != NULL ) {
    return mySelection->contains( myValues );
  }
  else {
    return true;
  }
}
bool csSelectionManager::contains( cseis_geolib::csFlexHeader const* hdrValue ) {
  for( int i = 0; i < myNumHeaders; i++ ) {
    if( myHeaderType[i] == cseis_geolib::TYPE_FLOAT ) {
      myValues[i].setDoubleValue( hdrValue->floatValue() );
    }
    else if( myHeaderType[i] == cseis_geolib::TYPE_DOUBLE ) {
      myValues[i].setDoubleValue( hdrValue->doubleValue());
    }
    else if( myHeaderType[i] == cseis_geolib::TYPE_INT ) {
      myValues[i].setIntValue( hdrValue->intValue() );
    }
    else { // TYPE_STRING
      return false; // Not sure what to do here...
    }
  }
  if( mySelection != NULL ) {
    return mySelection->contains( myValues );
  }
  else {
    return true;
  }
}
std::string csSelectionManager::headerName( int index ) const {
  if( index >= myNumHeaders ) throw( cseis_geolib::csException("csSelectionManager::headerName: Wrong index passed (%d). This is program bug in the calling function", index) );
  return myHeaderNames[index];
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TRACE_HEADER_H
#define CS_TRACE_HEADER_H

#include "csTraceHeaderData.h"
#include "csTraceHeaderDef.h"

#include <string>
#include <cstdio>
#include "csException.h"
#include "geolib_defines.h"
#include "csPoint3D.h"

namespace cseis_system {

/**
 * Trace header for one seismic trace
 *
 * This class gives access to the trace header values & trace header definition (const pointer object)
 * Users should see objects of this class only, not the trace header data object which is encapsulated in this class
 *
 *
 * @author Bjorn Olofsson
 * @date   2007
 */

class csTraceHeader {
public:
  csTraceHeader();
//  csTraceHeader( csTraceHeaderData* const traceHeaderData, csTraceHeaderDef const* traceHeaderDef );
  ~csTraceHeader();
  //--------------------------------------------------
  /**
  * Sets trace header definition and data
  * @param traceHeaderDef  Trace header definition object
  * @param inPort  Input port index
  */
  void setHeaders( csTraceHeaderDef const* traceHeaderDef, int inPort = 0 );
  /**
  * Delete trace headers from trace header block
  */
  void deleteHeaders( csTraceHeaderDef const* hdef ) {
    myTraceHeaderData->deleteHeaders( hdef );
  }
  /**
  * Copy trace header from another instance
  */
  void copyFrom( csTraceHeader const* hdr );
  //----------------------------------------------------------------
  inline int numHeaders() const {
    return myTraceHeaderData->numHeaders();
  }
  inline int getIndex( std::string const& name ) const {
    return myHeaderDefPtr->headerIndex( name );
  }
  inline cseis_geolib::type_t type( int index ) const {
    return myHeaderDefPtr->headerType( index );
  }
  inline cseis_geolib::type_t type( std::string const& name ) const {
    return myHeaderDefPtr->headerType( name );
  }
  inline std::string const name( int index ) const {
    return myHeaderDefPtr->headerName( index );
  }
  //  inline void setAbsoluteTime( int days, int sec, int us ) {
  //   myTraceHeaderData->setAbsoluteTime( days, sec, us );
  // }
  // inline void getAbsoluteTime( int& days, int& sec, int& us ) const {
  //  myTraceHeaderData->getAbsoluteTime( days, sec, us );
  // }
//------------------------------------------------------
  char const* getTraceHeaderValueBlock() const;
  char* getTraceHeaderValueBlockHandle();
  void setTraceHeaderValueBlock( char const* hdrValueBlock, int byteSize );
  void writeTraceHeaderValueBlock( char const* hdrValueBlock_in, int const* byteMap_in, int const* hdrMap_out, int numHeaders_in );
  void readTraceHeaderValueBlock( char* hdrValueBlock_out, int const* byteMap_out, int const* hdrMap_out, int numHeaders_out ) const;
//------------------------------------------------------
  void clear();
  void clearMemory();
  csTraceHeaderDef const* getHeaderDefPtr() const { return myHeaderDefPtr; }
  /**
   * Dump all trace header names and values
   */
  void dump( std::FILE* stream = NULL ) const;
private:
  /// The actual trace header values:
  csTraceHeaderData* myTraceHeaderData;
  /// Constant pointer to trace header definition (managed somewhere else)
  csTraceHeaderDef const* myHeaderDefPtr;

public:
  //---------------------------------------------------------------------
inline void setIntValue( int index, int value ) {
#ifdef CS_DEBUG
  if( index < numHeaders() && index >= 0) {
#endif
    cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
    switch( type ) {
    case cseis_geolib::TYPE_INT:
      myTraceHeaderData->setIntValue( index, value );
      break;
    case cseis_geolib::TYPE_FLOAT:
      myTraceHeaderData->setFloatValue( index, (float)value );
      break;
    case cseis_geolib::TYPE_DOUBLE:
      myTraceHeaderData->setDoubleValue( index, (double)value );
      break;
    case cseis_geolib::TYPE_INT64:
      myTraceHeaderData->setInt64Value( index, (csInt64_t)value );
      break;
#ifdef CS_DEBUG
    default:
      throw( cseis_geolib::csException("csTraceHeader::setIntValue(): Program bug. Header does not have a number type.") );
#endif
    }
#ifdef CS_DEBUG
  }
  else {
    throw( cseis_geolib::csException("csTraceHeader::setIntValue(): Program bug. Passed wrong header index.") );
  }
#endif
}

inline void setFloatValue( int index, float value ) {
#ifdef CS_DEBUG
    if( index < numHeaders() && index >= 0) {
#endif
      cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
      switch( type ) {
        case cseis_geolib::TYPE_FLOAT:
          myTraceHeaderData->setFloatValue( index, value );
          break;
        case cseis_geolib::TYPE_DOUBLE:
          myTraceHeaderData->setDoubleValue( index, (double)value );
          break;
        case cseis_geolib::TYPE_INT:
          myTraceHeaderData->setIntValue( index, (int)value );
          break;
        case cseis_geolib::TYPE_INT64:
          myTraceHeaderData->setInt64Value( index, (csInt64_t)value );
          break;
#ifdef CS_DEBUG
        default:
          throw( cseis_geolib::csException("csTraceHeader::setFloatValue(): Program bug. Header does not have a number type.") );
#endif
      }
#ifdef CS_DEBUG
    }
    else {
      throw( cseis_geolib::csException("csTraceHeader::setFloatValue(): Program bug. Passed wrong header index.") );
    }
#endif
  }
  //----------------------------------------------------------------
inline   void setDoubleValue( int index, double value ) {
#ifdef CS_DEBUG
    if( index < numHeaders() && index >= 0) {
#endif
      cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
      switch( type ) {
        case cseis_geolib::TYPE_DOUBLE:
          myTraceHeaderData->setDoubleValue( index, value );
          break;
        case cseis_geolib::TYPE_FLOAT:
          myTraceHeaderData->setFloatValue( index, (float)value );
          break;
        case cseis_geolib::TYPE_INT:
          myTraceHeaderData->setIntValue( index, (int)value );
          break;
        case cseis_geolib::TYPE_INT64:
          myTraceHeaderData->setInt64Value( index, (csInt64_t)value );
          break;
#ifdef CS_DEBUG
        default:
          throw( cseis_geolib::csException("csTraceHeader::setDoubleValue(): Program bug. Header does not have a number type.") );
#endif
      }
#ifdef CS_DEBUG
    }
    else {
      throw( cseis_geolib::csException("csTraceHeader::setDoubleValue(): Program bug. Passed wrong header index.") );
    }
#endif
 }
 inline void setInt64Value( int index, csInt64_t value ) {
#ifdef CS_DEBUG
   if( index < numHeaders() && index >= 0) {
#endif
     cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
     switch( type ) {
     case cseis_geolib::TYPE_INT64:
       myTraceHeaderData->setInt64Value( index, value );
       break;
     case cseis_geolib::TYPE_INT:
       myTraceHeaderData->setIntValue( index, (int)value );
       break;
     case cseis_geolib::TYPE_FLOAT:
       myTraceHeaderData->setFloatValue( index, (float)value );
       break;
     case cseis_geolib::TYPE_DOUBLE:
       myTraceHeaderData->setDoubleValue( index, (double)value );
       break;
#ifdef CS_DEBUG
     default:
       throw( cseis_geolib::csException("csTraceHeader::setInt64Value(): Program bug. Header does not have a number type.") );
#endif
     }
#ifdef CS_DEBUG
   }
   else {
     throw( cseis_geolib::csException("csTraceHeader::setInt64Value(): Program bug. Passed wrong header index.") );
   }
#endif
 }
 void setStringValue( int index, std::string value ) {
   myTraceHeaderData->setStringValue( index, value );
 }
 void setVectorValue( int index, cseis_geolib::csPoint3D value ) {
   myTraceHeaderData->setVectorValue( index, value );
 }
 void setVectorValue( int index, double value, cseis_geolib::type_t type ) {
   myTraceHeaderData->setVectorValue( index, value, type );
 }
 void setVectorValueX( int index, double value ) {
   myTraceHeaderData->setVectorValueX( index, value );
 }
 void setVectorValueY( int index, double value ) {
   myTraceHeaderData->setVectorValueY( index, value );
 }
 void setVectorValueZ( int index, double value ) {
   myTraceHeaderData->setVectorValueZ( index, value );
 }
//------------------------------------------------------
 inline   double doubleValue( int index ) const {
#ifdef CS_DEBUG
    if( index < numHeaders() && index >= 0) {
#endif
      cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
      switch( type ) {
        case cseis_geolib::TYPE_DOUBLE:
          return myTraceHeaderData->doubleValue( index );
          break;
        case cseis_geolib::TYPE_FLOAT:
          return (double)myTraceHeaderData->floatValue( index );
          break;
        case cseis_geolib::TYPE_INT:
          return (double)myTraceHeaderData->intValue( index );
          break;
        case cseis_geolib::TYPE_INT64:
          return (double)myTraceHeaderData->int64Value( index );
          break;
        default:
#ifdef CS_DEBUG
          throw( cseis_geolib::csException("csTraceHeader::doubleValue(): Program bug. Header does not have a number type.") );
#else
          return 0.0;
#endif
      }
#ifdef CS_DEBUG
    }
    else {
      throw( cseis_geolib::csException("csTraceHeader::doubleValue(): Program bug. Passed wrong header index.") );
    }
#endif
  }
//------------------------------------------------------
 inline  float floatValue( int index ) const {
#ifdef CS_DEBUG
    if( index < numHeaders() && index >= 0 ) {
#endif
      cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
      switch( type ) {
        case cseis_geolib::TYPE_FLOAT:
          return myTraceHeaderData->floatValue( index );
          break;
        case cseis_geolib::TYPE_DOUBLE:
          return (float)myTraceHeaderData->doubleValue( index );
          break;
        case cseis_geolib::TYPE_INT:
          return (float)myTraceHeaderData->intValue( index );
          break;
        case cseis_geolib::TYPE_INT64:
          return (float)myTraceHeaderData->int64Value( index );
          break;
        default:
#ifdef CS_DEBUG
          throw( cseis_geolib::csException("csTraceHeader::floatValue(): Program bug. Header does not have a number type.") );
#else
          return 0.0;
#endif
      }
#ifdef CS_DEBUG
    }
    else {
      throw( cseis_geolib::csException("csTraceHeader::floatValue(): Program bug. Passed wrong header index.") );
    }
#endif
  }
//------------------------------------------------------
inline int intValue( int index ) const {
#ifdef CS_DEBUG
    if( index < numHeaders() && index >= 0 ) {
#endif
      cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
      switch( type ) {
        case cseis_geolib::TYPE_INT:
          return myTraceHeaderData->intValue( index );
          break;
        case cseis_geolib::TYPE_FLOAT:
          return (int)myTraceHeaderData->floatValue( index );
          break;
        case cseis_geolib::TYPE_DOUBLE:
          return (int)myTraceHeaderData->doubleValue( index );
          break;
        case cseis_geolib::TYPE_INT64:
          return (int)myTraceHeaderData->int64Value( index );
          break;
        default:
#ifdef CS_DEBUG
          throw( cseis_geolib::csException("csTraceHeader::intValue(): Program bug. Header does not have a number type.") );
#else
          return 0;
#endif
      }
#ifdef CS_DEBUG
    }
    else {
      throw( cseis_geolib::csException("csTraceHeader::intValue(): Program bug. Passed wrong header index.") );
    }
#endif
  }
//------------------------------------------------------
inline csInt64_t int64Value( int index ) const {
#ifdef CS_DEBUG
    if( index < numHeaders() && index >= 0 ) {
#endif
      cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
      switch( type ) {
        case cseis_geolib::TYPE_INT64:
          return myTraceHeaderData->int64Value( index );
          break;
        case cseis_geolib::TYPE_INT:
          return (csInt64_t)myTraceHeaderData->intValue( index );
          break;
        case cseis_geolib::TYPE_FLOAT:
          return (csInt64_t)myTraceHeaderData->floatValue( index );
          break;
        case cseis_geolib::TYPE_DOUBLE:
          return (csInt64_t)myTraceHeaderData->doubleValue( index );
          break;
        default:
#ifdef CS_DEBUG
          throw( cseis_geolib::csException("csTraceHeader::int64Value(): Program bug. Header does not have a number type.") );
#else
          return 0;
#endif
      }
#ifdef CS_DEBUG
    }
    else {
      throw( cseis_geolib::csException("csTraceHeader::int64Value(): Program bug. Passed wrong header index.") );
    }
#endif
  }
//------------------------------------------------------
inline std::string const stringValue( int index ) const {
#ifdef CS_DEBUG
  if( index < numHeaders() && index >= 0 ) {
    cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
    if( type != cseis_geolib::TYPE_STRING ) {
      throw( cseis_geolib::csException("csTraceHeader::stringValue(): Program bug. Header is not of string type.") );
    }
  }
  else {
    throw( cseis_geolib::csException("csTraceHeader::stringValue(): Program bug. Passed wrong header index.") );
  }
#endif
  return myTraceHeaderData->stringValue( index );
}
//------------------------------------------------------
 inline cseis_geolib::csPoint3D const vectorValue( int index ) const {
#ifdef CS_DEBUG
  if( index < numHeaders() && index >= 0 ) {
    cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
    if( type != cseis_geolib::TYPE_VECTOR ) {
      if( index > numHeaders()-2 || 
          (myHeaderDefPtr->headerType(index) != cseis_geolib::TYPE_DOUBLE || myHeaderDefPtr->headerType(index+1) != cseis_geolib::TYPE_DOUBLE || myHeaderDefPtr->headerType(index+2) != cseis_geolib::TYPE_DOUBLE ) ) {
        throw( cseis_geolib::csException("csTraceHeader::vectorValue(): Program bug. Header is not of vector type.") );
      }
    }
  }
  else {
    throw( cseis_geolib::csException("csTraceHeader::vectorValue(): Program bug. Passed wrong header index.") );
  }
#endif
  return myTraceHeaderData->vectorValue( index );
}
//------------------------------------------------------
 inline double const vectorValue( int index, cseis_geolib::type_t type ) const {
#ifdef CS_DEBUG
  if( index < numHeaders() && index >= 0 ) {
    cseis_geolib::type_t type = myHeaderDefPtr->headerType(index);
    if( type != cseis_geolib::TYPE_VECTOR ) {
      throw( cseis_geolib::csException("csTraceHeader::vectorValue(): Program bug. Header is not of vector type.") );
    }
  }
  else {
    throw( cseis_geolib::csException("csTraceHeader::vectorValue(): Program bug. Passed wrong header index.") );
  }
#endif
  return myTraceHeaderData->vectorValue( index, type );
}

};

} // namespace

#endif


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
#include "cseis_defines.h"
#include "csTraceHeaderDef.h"
#include "csTraceHeaderInfo.h"
#include "csHeaderInfo.h"
#include "csGeolibUtils.h"
#include "csStandardHeaders.h"
#include "csMemoryPoolManager.h"
#include "csVector.h"
#include "csException.h"
#include <string>
#include <cstring>

using namespace cseis_system;

std::string csTraceHeaderDef::HEADER_NAME_UNKNOWN = std::string("UNKNOWN_HEADER");

csTraceHeaderDef::csTraceHeaderDef( csTraceHeaderDef* refHdrDefPtr ) {
  init( 1, refHdrDefPtr->myMemoryManager );
  // Always set absolute time headers. These must ALWAYS exist. Added 080724
  addHeader_internal( cseis_geolib::HDR_TIME_SAMP1.type, cseis_geolib::HDR_TIME_SAMP1.name, cseis_geolib::HDR_TIME_SAMP1.description, 1 );
  addHeader_internal( cseis_geolib::HDR_TIME_SAMP1_US.type, cseis_geolib::HDR_TIME_SAMP1_US.name, cseis_geolib::HDR_TIME_SAMP1_US.description, 1 );
}
csTraceHeaderDef::csTraceHeaderDef( csMemoryPoolManager* memManager ) {
  init( 1, memManager );
  // Always set absolute time headers. These must ALWAYS exist.
  addHeader_internal( cseis_geolib::HDR_TIME_SAMP1.type, cseis_geolib::HDR_TIME_SAMP1.name, cseis_geolib::HDR_TIME_SAMP1.description, 1 );
  addHeader_internal( cseis_geolib::HDR_TIME_SAMP1_US.type, cseis_geolib::HDR_TIME_SAMP1_US.name, cseis_geolib::HDR_TIME_SAMP1_US.description, 1 );
}
csTraceHeaderDef::csTraceHeaderDef( int numInputPorts, csTraceHeaderDef const** hdef, csMemoryPoolManager* memManager ) {
  if( numInputPorts > 0 ) {
    init( numInputPorts, memManager );
    initInputPorts( hdef );
  }
  else {
    init( 1, memManager );
    // Always set absolute time headers. These must ALWAYS exist. Added 080724
    addHeader_internal( cseis_geolib::HDR_TIME_SAMP1.type, cseis_geolib::HDR_TIME_SAMP1.name, cseis_geolib::HDR_TIME_SAMP1.description, 1 );
    addHeader_internal( cseis_geolib::HDR_TIME_SAMP1_US.type, cseis_geolib::HDR_TIME_SAMP1_US.name, cseis_geolib::HDR_TIME_SAMP1_US.description, 1 );
  }
}

int csTraceHeaderDef::mpi_getByteSize() const {
  char* dataBuffer = new char[10000];
  int byteSize = mpi_compress( dataBuffer );
  delete [] dataBuffer;
  return byteSize;
}
int csTraceHeaderDef::mpi_compress( char* data ) const {
  int byteLoc = 0;
  int nHeaders = numHeaders();  
  memcpy( &data[byteLoc], &nHeaders, sizeof(int) );
  byteLoc += sizeof(int);
  memcpy( &data[byteLoc], &myTotalNumBytes, sizeof(int) );
  byteLoc += sizeof(int);

  //  fprintf(stderr,"SEND numHeaders: %d, numToAdd[0]: %d\n", numHeaders(), myTraceHeadersToAdd[0].size() );
  for( int inPort = 0; inPort < myNumInputPorts; inPort++ ) {
    int numToAdd = myTraceHeadersToAdd[inPort].size();
    memcpy( &data[byteLoc], &numToAdd, sizeof(int) );
    byteLoc += sizeof(int);
    for( int ihdr = 0; ihdr < numToAdd; ihdr++ ) {
      csTraceHeaderInfo const* info = myTraceHeadersToAdd[inPort].at(ihdr);
      memcpy( &data[byteLoc], &info->type, sizeof(char) );
      byteLoc += sizeof(char);
      memcpy( &data[byteLoc], &info->nElements, sizeof(short) );
      byteLoc += sizeof(short);
      int nameLen = info->name.length();
      int descLen = info->description.length();
      //      fprintf(stderr,"Send name/desc length: %d %d\n", nameLen, descLen );
      memcpy( &data[byteLoc], &nameLen, sizeof(int) );
      byteLoc += sizeof(int);
      memcpy( &data[byteLoc], &descLen, sizeof(int) );
      byteLoc += sizeof(int);
      if( nameLen > 0 ) {
        memcpy( &data[byteLoc], info->name.c_str(), nameLen*sizeof(char) );
        byteLoc += nameLen*sizeof(char);
      }
      if( descLen > 0 ) {
        memcpy( &data[byteLoc], info->description.c_str(), descLen*sizeof(char) );
        byteLoc += descLen*sizeof(char);
      }
    }
  }

  int numDel = myIndexOfHeadersToDel->size();
  memcpy( &data[byteLoc], &numDel, sizeof(int) );
  byteLoc += sizeof(int);
  for( int i = 0; i < numDel; i++ ) {
    memcpy( &data[byteLoc], &myIndexOfHeadersToDel->at(i), sizeof(int) );
    byteLoc += sizeof(int);
  }

  memcpy( &data[byteLoc], myByteLocation, nHeaders*sizeof(int) );
  byteLoc += nHeaders*sizeof(int);
  return byteLoc;
}

void csTraceHeaderDef::mpi_decompress( char const* data ) {
  int byteLoc = 0;
  int nHeadersCurrent = numHeaders();
  int nHeaders = 0;
  memcpy( &nHeaders, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  int nHeadersDiff = nHeaders - nHeadersCurrent;
  if( nHeadersDiff < 0 ) throw( cseis_geolib::csException("mpi_decompress(): More trace headers in recipient than sending MPI proc: %d > %d", nHeadersCurrent, nHeaders) );
  
  memcpy( &myTotalNumBytes, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);

  for( int inPort = 0; inPort < myNumInputPorts; inPort++ ) {
    myTraceHeadersToAdd[inPort].clear();
    int numToAdd;
    memcpy( &numToAdd, &data[byteLoc], sizeof(int) );
    byteLoc += sizeof(int);
    for( int ihdr = 0; ihdr < numToAdd; ihdr++ ) {
      csTraceHeaderInfo info(0,"","");
      memcpy( &info.type, &data[byteLoc], sizeof(char) );
      byteLoc += sizeof(char);
      memcpy( &info.nElements, &data[byteLoc], sizeof(short) );
      byteLoc += sizeof(short);
      int nameLen = 0;
      int descLen = 0;
      memcpy( &nameLen, &data[byteLoc], sizeof(int) );
      byteLoc += sizeof(int);
      memcpy( &descLen, &data[byteLoc], sizeof(int) );
      byteLoc += sizeof(int);
      //      fprintf(stderr,"RECV name/desc length: %d %d\n", nameLen, descLen );
      char* textStr = new char[std::max(nameLen,descLen)+1];
      if( nameLen > 0 ) {
        memcpy( textStr, &data[byteLoc], nameLen*sizeof(char) );
        textStr[nameLen] = '\0';
        info.name = textStr;
        byteLoc += nameLen*sizeof(char);
      }
      if( descLen > 0 ) {
        memcpy( textStr, &data[byteLoc], descLen*sizeof(char) );
        textStr[descLen] = '\0';
        info.description = textStr;
        byteLoc += descLen*sizeof(char);
      }
      delete [] textStr;
      csTraceHeaderInfo const* infoPtr = myMemoryManager->getNewTraceHeaderInfo( info.type, info.nElements, info.name, info.description );
      myTraceHeadersToAdd[inPort].insertEnd( infoPtr );
      // Add trace headers to main list only if not existing yet. For example, INPUT modules have time_samp1 and time_samp1_us already set by default.
      if( numToAdd-ihdr <= nHeadersDiff ) myTraceHeaderInfoList->insertEnd( infoPtr );
    }
  }

  int numDel = 0;
  memcpy( &numDel, &data[byteLoc], sizeof(int) );
  byteLoc += sizeof(int);
  myIndexOfHeadersToDel->clear();
  for( int i = 0; i < numDel; i++ ) {
    int value;
    memcpy( &value, &data[byteLoc], sizeof(int) );
    myIndexOfHeadersToDel->insertEnd(value);
    byteLoc += sizeof(int);
  }

  if( myByteLocation != NULL ) delete [] myByteLocation;
  myByteLocation = new int[nHeaders];
  memcpy( myByteLocation, &data[byteLoc], nHeaders*sizeof(int) );
  byteLoc += nHeaders*sizeof(int);

  if( nHeaders != numHeaders() ) {
    for( int ihdr = 0; ihdr < numHeaders(); ihdr++ ) {
      csTraceHeaderInfo const* info = myTraceHeaderInfoList->at(ihdr);
      fprintf(stderr," Header #%d: %d  %s  %s\n", ihdr, info->type, info->name.c_str(), info->description.c_str());
    }
    throw( cseis_geolib::csException("mpi_decompress(): Inconsistent number of trace headers: %d != %d", nHeaders, numHeaders()) );
  }
}

void csTraceHeaderDef::init( int numInputPorts, csMemoryPoolManager* memManager ) {
  myMemoryManager = memManager;
  myNumInputPorts = numInputPorts;

  myTraceHeaderInfoList = new cseis_geolib::csVector<csTraceHeaderInfo const*>;
  myTotalNumBytes = 0;
  myByteLocation = NULL;
  myIndexOfHeadersToDel = new cseis_geolib::csVector<int>(0);
  myTraceHeadersToAdd = new cseis_geolib::csVector<csTraceHeaderInfo const*>[numInputPorts];
}
csMemoryPoolManager* csTraceHeaderDef::getMemoryManager() {
  return myMemoryManager;
}

//----------------------------------------------
//
void csTraceHeaderDef::initInputPorts( csTraceHeaderDef const** hdefPrev ) {
  // Copy all headers from first input port
  csTraceHeaderDef const* hdef = hdefPrev[0];
  int nHeaders = hdef->numHeaders();
  cseis_geolib::csVector<int> const* indexOfHeadersToDel = hdef->getIndexOfHeadersToDel();
  int nHeadersToDel = indexOfHeadersToDel->size();
  int counterDel = 0;
  for( int iHdr = 0; iHdr < nHeaders; iHdr++ ) {
    csTraceHeaderInfo const* info = hdef->myTraceHeaderInfoList->at(iHdr);
    // This header shall be deleted, not to be passed on...
    if( nHeadersToDel > counterDel && indexOfHeadersToDel->at(counterDel) == iHdr ) {
      counterDel += 1;
    }
    else {
      myTraceHeaderInfoList->insertEnd( info );
    }
  }
  myTotalNumBytes = hdefPrev[0]->myTotalNumBytes;

  if( myNumInputPorts > 1 ) {
    // Copy all new headers from input ports 1,2,..N
    for( int inPort = 1; inPort < myNumInputPorts; inPort++ ) {
      hdef = hdefPrev[inPort];
      nHeaders = hdef->numHeaders();
      indexOfHeadersToDel = hdef->getIndexOfHeadersToDel();
      nHeadersToDel = indexOfHeadersToDel->size();
      counterDel = 0;
      for( int iHdr = 0; iHdr < nHeaders; iHdr++ ) {
        csTraceHeaderInfo const* info = hdef->myTraceHeaderInfoList->at(iHdr);
        int index = 0;
        // This header shall be deleted, not to be passed on...
        if( nHeadersToDel > counterDel && indexOfHeadersToDel->at(counterDel) == iHdr ) {
          counterDel += 1;
          continue;
        }
        if( !getIndex(info->name,index) ) {
          myTraceHeaderInfoList->insertEnd( info );
          myTotalNumBytes += cseis_geolib::csGeolibUtils::numBytes( info->type )*info->nElements;
        }
        else { // Header already exists -> check if type is the same
          if( myTraceHeaderInfoList->at( index )->type != info->type ) {
            throw( cseis_geolib::csException("Trace header from different input ports already exists but has different type...") );
          }
          else if( myTraceHeaderInfoList->at( index )->nElements != info->nElements ) {
            throw( cseis_geolib::csException("Trace header from different input ports already exists but has different number of elements...") );
          }
        }
      }
    }
    // Check every input port; if trace header does not exist yet, add it to this port. If header exists and is at a different location --> Exception
    nHeaders = numHeaders();
    for( int inPort = 0; inPort < myNumInputPorts; inPort++ ) {
      hdef = hdefPrev[inPort];
      myTraceHeadersToAdd[inPort].clear();
      int ihead = 0;
      while( ihead < nHeaders && hdef->headerName(ihead) == headerName(ihead) ) { ihead++; }
      for( ; ihead < nHeaders; ihead++ ) {
        csTraceHeaderInfo const* info = myTraceHeaderInfoList->at(ihead);
        int indexDummy = 0;
        if( !hdef->getIndex( info->name, indexDummy ) ) {
          myTraceHeadersToAdd[inPort].insertEnd( info );
        }
        else {
//          bool comp = (hdef->headerName(ihead) == headerName(ihead));
//          bool comp2 = hdef->headerName(ihead).compare(headerName(ihead));
//          fprintf(stdout,"CHECK-X Header %d %d: %s  .. %d .. %d ..\n", inPort, ihead, hdef->headerName(ihead).c_str(), comp, comp2);
          // !CHANGE! Add capability to merge input ports with chaotic differences in trace headers, by reshuffling trace headers when merging
          throw( cseis_geolib::csException("Trace headers from different input ports of this module do not match, or have different order. Trace header not found: %s. Make sure this trace header exists in all input streams. Merge not possible.", info->name.c_str()) );
        }
      } // END for ihead
    } // END for inPort
  } // END   if( myNumInputPorts > 1 ) {
}
//----------------------------------------------------------------
//
csTraceHeaderDef::~csTraceHeaderDef() {
  if( myTraceHeaderInfoList != NULL ) {
    delete myTraceHeaderInfoList;
    myTraceHeaderInfoList = NULL;
  }
  if( myTraceHeadersToAdd != NULL ) {
    delete [] myTraceHeadersToAdd;
    myTraceHeadersToAdd = NULL;
  }
  if( myIndexOfHeadersToDel != NULL ) {
    delete myIndexOfHeadersToDel;
    myIndexOfHeadersToDel = NULL;
  }
  if( myByteLocation != NULL ) {
    delete [] myByteLocation;
    myByteLocation = NULL;
  }
}
int csTraceHeaderDef::numHeaders() const {
  return myTraceHeaderInfoList->size();
}
//----------------------------------------------------------------
//
bool csTraceHeaderDef::equals( csTraceHeaderDef const* hdef ) const {
  int nHeaders = numHeaders();
  if( nHeaders != hdef->numHeaders() || getTotalNumBytes() != hdef->getTotalNumBytes() ) return false;
  for( int i = 0; i < nHeaders; i++ ) {
    csTraceHeaderInfo const* info1 = myTraceHeaderInfoList->at( i );
    csTraceHeaderInfo const* info2 = hdef->myTraceHeaderInfoList->at( i );
    if( info1->name.compare( info2->name ) || (info1->type != info2->type ) ) {
      return false;
    }
  }
  return true;
}

/*
  int index = 0;
  if( !getIndex( name, index ) ) {
    throw( cseis_geolib::csException("Trace header not found: '%s'\n", name.c_str() ) );
  }
  return index;
*/
//----------------------------------------------------------------
int csTraceHeaderDef::headerIndex( std::string const& name ) const {
  int index = 0;
  if( !getIndex( name, index ) ) {
    // Test if this is a vector header:
    int pos = name.find_first_of('.');
    if( pos != (int)std::string::npos ) {
      std::string nameVec = name.substr(0,pos);
      //      fprintf(stderr,"Vec name %s\n", nameVec.c_str());
      if( !getIndex( nameVec, index ) ) {
        throw( cseis_geolib::csException("Trace header not found: '%s'\n", name.c_str() ) );
      }
    }
    else {
    // Test if this is a vector header:
    //    std::string nameVecX = name + ".x";
    // if( !getIndex( nameVecX, index ) ) {
      throw( cseis_geolib::csException("Trace header not found: '%s'\n", name.c_str() ) );
    }
  }
  return index;
}
//----------------------------------------------------------------
//
int csTraceHeaderDef::numElements( int index ) const {
  if( index >= 0 && index < myTraceHeaderInfoList->size() ) {
    return myTraceHeaderInfoList->at(index)->nElements;
  }
  else {
    return 0;
  }
}
//----------------------------------------------------------------
//
cseis_geolib::type_t csTraceHeaderDef::headerType( int index ) const {
  if( myTraceHeaderInfoList == NULL ) return 0;
  if( index >= 0 && index < myTraceHeaderInfoList->size() ) {
    return myTraceHeaderInfoList->at(index)->type;
  }
  else {
    return cseis_geolib::TYPE_UNKNOWN;
  }
}
//----------------------------------------------------------------
//
cseis_geolib::type_t csTraceHeaderDef::headerType( std::string const& name ) const {
  int index = 0;
  if( !getIndex( name, index ) ) {
    // Test if this is a vector header:
    int pos = name.find_first_of('.');
    if( pos != (int)std::string::npos ) {
      std::string nameVec = name.substr(0,pos);
      //      fprintf(stderr,"Vec name %s\n", nameVec.c_str());
      if( getIndex( nameVec, index ) ) {
        std::string letter = name.substr(pos+1);
        //        fprintf(stderr,"Vec letter %s\n", letter.c_str());
        if( letter.compare("x") == 0 ) return cseis_geolib::TYPE_VECTOR_X;
        else if( letter.compare("y") == 0 ) return cseis_geolib::TYPE_VECTOR_Y;
        else if( letter.compare("z") == 0 ) return cseis_geolib::TYPE_VECTOR_Z;
      }
    }
    return cseis_geolib::TYPE_UNKNOWN;
  }
  return headerType( index );
}
//----------------------------------------------------------------
//
std::string csTraceHeaderDef::headerName( int index ) const {
  if( index >= 0 && index < myTraceHeaderInfoList->size() ) {
    return myTraceHeaderInfoList->at(index)->name;
  }
  else {
    return HEADER_NAME_UNKNOWN;
  }
}
std::string csTraceHeaderDef::headerDesc( int index ) const {
  if( index >= 0 && index < myTraceHeaderInfoList->size() ) {
    return myTraceHeaderInfoList->at(index)->description;
  }
  else {
    return HEADER_NAME_UNKNOWN;
  }
}
//----------------------------------------------------------------
//
bool csTraceHeaderDef::headerExists( std::string const& name ) const {
  int index = 0;
  if( getIndex( name, index ) ) {
    return true;
  }
  else {
    //    fprintf(stderr,"headerexists %s\n", name.c_str() );
    // Test if this is a vector header XYZ:
    int pos = name.find_first_of('.');
    if( pos != (int)std::string::npos ) {
      std::string nameVec = name.substr(0,pos);
      //      fprintf(stderr,"Vec name %s\n", nameVec.c_str());
      return( getIndex( nameVec, index ) );
    }
    return false;
  }
}
//----------------------------------------------------------------
//
bool csTraceHeaderDef::getIndex( std::string const& name, int& index ) const {
  for( index = 0; index < myTraceHeaderInfoList->size(); index++ ) {
    //    if( name.compare("plane_p1") == 0 ) fprintf(stderr,"%d %s\n", index, myTraceHeaderInfoList->at(index)->name.c_str() );
    if( myTraceHeaderInfoList->at(index)->name == name ) {
      return true;
    }
  }
  index = HEADER_NOT_FOUND;
  return false;
}
//----------------------------------------------------------------
//
csTraceHeaderInfo const* csTraceHeaderDef::headerInfo( int index ) const {
  return myTraceHeaderInfoList->at(index);
}
//----------------------------------------------------------------
//
//int csTraceHeaderDef::addStandardHeader( cseis_geolib::csHeaderInfo const* info ) {
//  return addStandardHeader( info->name );
//}
int csTraceHeaderDef::addStandardHeader( std::string const& name ) {
  cseis_geolib::csHeaderInfo const* info = cseis_geolib::csStandardHeaders::get( name );
  if( info == NULL ) {
    throw( cseis_geolib::csException("Cannot add standard trace header '%s': Standard trace header with this name does not exist. Specify name and type if new non-standard header shall be added.", name.c_str() ) );
  }
  else {
    return addHeader_internal( info->type, info->name, info->description, 1 );
  }
}
int csTraceHeaderDef::addHeader( cseis_geolib::type_t type, std::string const& name, int nElements ) {
  cseis_geolib::csHeaderInfo const* info = cseis_geolib::csStandardHeaders::get( name );
  if( info == NULL ) {
    return addHeader_internal( type, name, "NONE", nElements );
  }
  else if( info->type == type ) {
    return addHeader_internal( type, name, info->description, nElements );
  }
  else {
    throw( cseis_geolib::csException("Cannot add new trace header '%s': Standard header exists with same name but with different type (%s).",
           name.c_str(), cseis_geolib::csGeolibUtils::typeText(info->type)) );
  }
}
int csTraceHeaderDef::addHeader( cseis_geolib::csHeaderInfo const* info, int nElements ) {
  return addHeader( info->type, info->name, info->description, nElements );
}
int csTraceHeaderDef::addHeader( cseis_geolib::type_t type, std::string const& name, std::string const& description, int nElements ) {
  if( name.length() < 1 ) throw( cseis_geolib::csException("csTraceHeaderDef::addHeader: Header name empty. This is a program bug in the calling function") );
  cseis_geolib::csHeaderInfo const* info = cseis_geolib::csStandardHeaders::get( name );
  if( info == NULL || (info->type == type ) ) {
    //  if( info == NULL || (info->type == type && (!info->description.compare(description)) ) ) { // 100510: COMMENTED out to avoid breaking existing module code
    return addHeader_internal( type, name, description, nElements );
  }
  else if( info->type != type ) {
    throw( cseis_geolib::csException("Cannot add new trace header '%s': Standard header with same name exists, but with different type (%s).",
           name.c_str(), cseis_geolib::csGeolibUtils::typeText(info->type)) );
  }
  else { //if( info->description ) {
    throw( cseis_geolib::csException("Cannot add new trace header '%s': Standard header with same name exists, but with different description (%s).",
           name.c_str(), info->description.c_str() ) );
  }
}
int csTraceHeaderDef::addHeader_internal( cseis_geolib::type_t type, std::string const& name, std::string const& description, int nElements ) {
#ifdef CS_DEBUG
  if( type == cseis_geolib::TYPE_STRING && nElements <= 1 ) {
    throw cseis_geolib::csException("csTraceHeaderDef::addHeader: Attempted to add STRING header of length = %d. Must specify length of string (argument nElements). This is a program bug in the calling function.", nElements);
  }
#endif
  int index = 0;
  if( !getIndex( name, index ) ) {
    csTraceHeaderInfo const* info = myMemoryManager->getNewTraceHeaderInfo( type, nElements, name, description );
    myTraceHeaderInfoList->insertEnd( info );
    index = numHeaders()-1;
  for( int iport = 0; iport < myNumInputPorts; iport++ ) {
      myTraceHeadersToAdd[iport].insertEnd( info );
    }
    myTotalNumBytes += cseis_geolib::csGeolibUtils::numBytes( type ) * nElements;
    // fprintf(stderr,"New total num bytes: %d (header %s) %d %d\n", myTotalNumBytes, info->name.c_str(), type, nElements);
    return( index );
  }
  else if( type == myTraceHeaderInfoList->at(index)->type ) {
    if( nElements == myTraceHeaderInfoList->at(index)->nElements ) {
      return( index );
    }
    else {
      throw( cseis_geolib::csException("csTraceHeaderDef::addHeader(): Trace header already exists but has different number of elements") );
    }
  }
  else {
    throw( cseis_geolib::csException("csTraceHeaderDef::addHeader(): Trace header already exists but has different type") );
  }
}
//----------------------------------------------------------------
//
int csTraceHeaderDef::getByteLocation( int index ) const {
#ifdef CS_DEBUG
  if( index < 0 || index >= myTraceHeaderInfoList->size() ) {
    throw( cseis_geolib::csException("csTraceHeaderDef::getByteLocation: Wrong header index passed to function") );
  }
#endif
  return myByteLocation[index];  
}
//----------------------------------------------------------------
//
int csTraceHeaderDef::deleteHeader( std::string const& name ) {
  try {
    int index = 0;
    if( getIndex( name, index ) ) {
      if( !isSystemTraceHeader( index ) ) {
        deleteHeader( index );
//      fprintf(stdout,"Delete header: %s index %d  (%d)\n", name.c_str(), index, insertAt);
        return index;
      }
    }
  }
  catch( ... ) {
    //
  }
  return HEADER_NOT_FOUND;
}
void csTraceHeaderDef::deleteAllHeaders() {
  try {
    for( int index = 0; index < numHeaders(); index++ ) {
      if( !isSystemTraceHeader( index ) ) {
        deleteHeader( index );
      }
    }
  }
  catch( ... ) {
    //
  }
}
void csTraceHeaderDef::deleteHeader( int index ) {
  if( index < 0 || index >= numHeaders() ) throw cseis_geolib::csException("csTraceHeaderDef::deleteHeader(): Wrong trace header index passed to method.");
  int insertAt = 0;
  int size = myIndexOfHeadersToDel->size();
  while( insertAt < size && myIndexOfHeadersToDel->at(insertAt) < index ) {
    insertAt++;
  }
  myIndexOfHeadersToDel->insert( index, insertAt );
}
cseis_geolib::csVector<int> const* csTraceHeaderDef::getIndexOfHeadersToDel() const {
  return myIndexOfHeadersToDel;
}
void csTraceHeaderDef::dump( FILE* fout ) const {
  fprintf(fout,"********* csTraceHeaderDef::dump(), total num bytes: %d *********\n", myTotalNumBytes);
  for( int i = 0; i < myTraceHeaderInfoList->size(); i++ ) {
    csTraceHeaderInfo const* info = myTraceHeaderInfoList->at( i );
    fprintf(fout,"Info %2d, type %8s, byte %d: '%-20s', Desc: '%s'\n", i, cseis_geolib::csGeolibUtils::typeText(info->type),
            getByteLocation(i), info->name.c_str(), info->description.c_str() );
  }
}
void csTraceHeaderDef::resetByteLocation() {
  if( myByteLocation != NULL ) {
    delete [] myByteLocation;
    myByteLocation = NULL;
  }
  int nHeaders = numHeaders();
  myByteLocation = new int[nHeaders];
  myByteLocation[0] = 0;
  int numBytes = 0;
  for( int ihdr = 1; ihdr < nHeaders; ihdr++ ) {
    csTraceHeaderInfo const* info = myTraceHeaderInfoList->at(ihdr-1);
    if( info->type != cseis_geolib::TYPE_STRING ) {
      numBytes += cseis_geolib::csGeolibUtils::numBytes( info->type );
    }
    else { // if( type == cseis_geolib::TYPE_STRING ) {
      // BUGFIX 080704: Number of bytes for string headers were previously computed as length of description string. This was preliminary code
      numBytes += info->nElements;
    }
    myByteLocation[ihdr] = numBytes;
  }
}

bool csTraceHeaderDef::isSystemTraceHeader( std::string const& name ) const {
  int index = 0;
  getIndex( name, index );
  return( isSystemTraceHeader( index ) );
}

bool csTraceHeaderDef::isSystemTraceHeader( int index ) {
  return( index == csTraceHeaderDef::HDR_ID_TIME_SAMP1_S || index == csTraceHeaderDef::HDR_ID_TIME_SAMP1_US );
}


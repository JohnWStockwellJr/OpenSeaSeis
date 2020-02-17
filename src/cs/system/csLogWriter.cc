/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csLogWriter.h"
#include "csException.h"
#include <stdarg.h>
#include <cstring>

using namespace cseis_system;

csLogWriter::csLogWriter( char const* filename ) : myFilename( filename ), myIsOpenedLocally(false) {
  myLogFile = NULL;
  if( myFilename != NULL ) {
    myLogFile = fopen( myFilename, "w" );
    if( myLogFile == NULL ) {
      throw cseis_geolib::csException("Can't open file: ");
    }
    myIsOpenedLocally = true;
  }
}
csLogWriter::csLogWriter( FILE* file ) : myFilename( "Unknown" ), myIsOpenedLocally( false ) {
  myLogFile = file;
}
csLogWriter::csLogWriter() {
  myFilename = NULL;
  myIsOpenedLocally = false;
  myLogFile = stdout;
}
//-----------------------------------------
csLogWriter::~csLogWriter() {
  if( myIsOpenedLocally && myLogFile ) {
    fclose( myLogFile );
    myLogFile = NULL;
  }
}
//-----------------------------------------
void csLogWriter::line( char const* text, ... ) {
  if( myLogFile ) {
    va_list argList;
    va_start( argList, text );
    vfprintf( myLogFile, text, argList );
    fprintf( myLogFile, "\n" );
#ifdef CS_DEBUG
    fflush(myLogFile);
#endif
  }
}
//-----------------------------------------
void csLogWriter::write( char const* text, ... ) {
  if( myLogFile ) {
    va_list argList;
    va_start( argList, text );
    vfprintf( myLogFile, text, argList );
  #ifdef CS_DEBUG
    fflush(myLogFile);
  #endif
  }
}
//-----------------------------------------
void csLogWriter::error( char const* text, ... ) {
  if( myLogFile ) {
    va_list argList;
    va_start( argList, text );
    fprintf( myLogFile, "FATAL ERROR:\n" );
    vfprintf( myLogFile, text, argList );
    fprintf( myLogFile, "\n" );
  }
  throw( cseis_geolib::csException("Fatal error occurred. See log file for details.") );  
}
void csLogWriter::error( char const* moduleName, int moduleIndex, char const* text, ... ) {
  if( myLogFile ) {
    va_list argList;
    va_start( argList, text );
    if( moduleName != NULL ) {
      if( moduleIndex > 0 ) {
        fprintf( myLogFile, "\nFATAL ERROR in module %s (%d):\n", moduleName, moduleIndex+1 );
      }
      else {
        fprintf( myLogFile, "\nFATAL ERROR in module %s:\n", moduleName );
      }
    }
    else {
      fprintf( myLogFile, "\nFATAL ERROR:\n" );
    }
    vfprintf( myLogFile, text, argList );
    fprintf( myLogFile, "\n" );
  }
  throw( cseis_geolib::csException("Fatal error occurred. See log file for details.") );
}
//-----------------------------------------
void csLogWriter::warning( char const* text, ... ) {
  if( myLogFile ) {
    va_list argList;
    va_start( argList, text );
    fprintf( myLogFile, "WARNING:\n" );
    vfprintf( myLogFile, text, argList );
    fprintf( myLogFile, "\n" );
  }
}
//-----------------------------------------
void csLogWriter::flush() {
  if( myLogFile ) {
    fflush(myLogFile);
  }
}

char const* csLogWriter::getFilename() const {
  return myFilename;
}



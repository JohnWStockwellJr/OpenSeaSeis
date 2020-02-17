/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csFileUtils.h"
#include "csVector.h"
#include "csException.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstring>
#include <limits>
#include "geolib_platform_dependent.h"
#include "geolib_defines.h"

extern "C" {
#include <sys/stat.h>
#include <dirent.h>
#include <limits.h>
}

using namespace cseis_geolib;

csFileUtils::csFileUtils() {
//  myFileSize = FILESIZE_UNKNOWN;
}
csFileUtils::~csFileUtils() {
}

csInt64_t csFileUtils::retrieveFileSize( std::string const& filename ) {
  return determineFileSize( filename );
}

bool csFileUtils::retrieveTimeStamp( std::string const& filename, int* timeStamp_s ) {
  csInt64_t fileSize;
  return csFileUtils::retrieveFileInfo( filename, &fileSize, timeStamp_s );
}

bool csFileUtils::retrieveFileInfo( std::string const& filename, csInt64_t* fileSize, int* timeStamp_s ) {
#ifdef PLATFORM_WINDOWS
  struct _stati64 statField;
  _stati64( filename.c_str(), &statField );
  csInt64_t fileSize64 = statField.st_size;

  if( fileSize64 >= 0 ) {
    *fileSize = (csInt64_t)fileSize64;
  }
  else {
    *fileSize = FILESIZE_UNKNOWN;
  }
  struct stat entryInfo;
  if( stat( filename.c_str(), &entryInfo ) == 0 ) {
    *timeStamp_s = (int)entryInfo.st_mtime; // Time of last modification. !CHANGE! Is st_mtime in [s], [ms] or [us]???
  }
  else {
    *timeStamp_s = 0;
  }
#elif PLATFORM_APPLE
  struct stat64 statField;
  FILE* fileC = fopen( filename.c_str(), "rb" );
  if( fileC == NULL ) {
    return false;
    //    throw( csException("Could not open input file '%s'\n", filename.c_str() ) );    
  }
  stat64( filename.c_str(), &statField );
  if( statField.st_size >= 0 ) {
    *fileSize = (csInt64_t)statField.st_size;
  }
  else {
    if( fileC == NULL ) {
      return false;
      //      throw( csException("Could not open input file '%s'\n", filename.c_str() ) );
    }
    off64_t fileSize64 = 0;
    fseeko64( fileC, 0, SEEK_END );
    fileSize64 = ftello64( fileC );
    if( fileSize64 >= 0 ) {
      *fileSize = (csInt64_t)fileSize64;
    }
    else {
      *fileSize = FILESIZE_UNKNOWN;
    }
  }
  fclose(fileC);

  struct stat entryInfo;
  if( stat( filename.c_str(), &entryInfo ) == 0 ) {
    *timeStamp_s = entryInfo.st_mtime; // Time of last modification.
  }
  else {
    *timeStamp_s = 0;
  }
#else
  struct stat64 statField;
  FILE* fileC = fopen( filename.c_str(), "rb" );
  if( fileC == NULL ) {
    return false;
    //    throw( csException("Could not open input file '%s'\n", filename.c_str() ) );    
  }
  stat64( filename.c_str(), &statField );
  if( statField.st_size >= 0 ) {
    *fileSize = (csInt64_t)statField.st_size;
  }
  else {
    if( fileC == NULL ) {
      return false;
      //      throw( csException("Could not open input file '%s'\n", filename.c_str() ) );
    }
    off64_t fileSize64 = 0;
    fseeko64( fileC, 0, SEEK_END );
    fileSize64 = ftello64( fileC );
    if( fileSize64 >= 0 ) {
      *fileSize = (csInt64_t)fileSize64;
    }
    else {
      *fileSize = FILESIZE_UNKNOWN;
    }
  }
  fclose(fileC);

  struct stat entryInfo;
  if( stat( filename.c_str(), &entryInfo ) == 0 ) {
    *timeStamp_s = entryInfo.st_mtim.tv_sec; // Time of last modification.
  }
  else {
    *timeStamp_s = 0;
  }
#endif
  return true;
}

csInt64_t csFileUtils::determineFileSize( std::string const& filename ) {
  csInt64_t fileSize = 0;
#ifdef PLATFORM_WINDOWS
  struct _stati64 statField;
  _stati64( filename.c_str(), &statField );
  csInt64_t fileSize64 = statField.st_size;

  if( fileSize64 >= 0 ) {
    fileSize = (csInt64_t)fileSize64;
  }
  else {
    fileSize = FILESIZE_UNKNOWN;
  }
#else
  struct stat64 statField;
  FILE* fileC = fopen( filename.c_str(), "rb" );
  if( fileC == NULL ) {
    throw( csException("Could not open input file '%s'\n", filename.c_str() ) );    
  }
  stat64( filename.c_str(), &statField );
  if( statField.st_size >= 0 ) {
    fileSize = (csInt64_t)statField.st_size;
  }
  else {
    if( fileC == NULL ) {
      throw( csException("Could not open input file '%s'\n", filename.c_str() ) );
    }
    off64_t fileSize64 = 0;
    fseeko64( fileC, 0, SEEK_END );
    fileSize64 = ftello64( fileC );
    if( fileSize64 >= 0 ) {
      fileSize = (csInt64_t)fileSize64;
    }
    else {
      fileSize = FILESIZE_UNKNOWN;
    }
  }
  fclose(fileC);

#endif
  return fileSize;
}
//----------------------------------------------------------------------------
//
bool csFileUtils::retrieveFilesStartingWith( std::string const& directory, std::string const& startString, cseis_geolib::csVector<std::string>* fileList_out, bool searchSubdirs, FILE* msgStream ) {
  return( csFileUtils::retrieveFiles_internal( directory, startString, fileList_out, csFileUtils::OPTION_START_STRING, searchSubdirs, msgStream ) );
}
bool csFileUtils::retrieveFiles( std::string const& directory, std::string const& extension, cseis_geolib::csVector<std::string>* fileList_out, bool searchSubdirs, FILE* msgStream ) {
  return( csFileUtils::retrieveFiles_internal( directory, extension, fileList_out, csFileUtils::OPTION_EXTENSION, searchSubdirs, msgStream ) );
}
bool csFileUtils::retrieveFiles_internal( std::string const& directory, std::string const& text, cseis_geolib::csVector<std::string>* fileList_out, int selectOption, bool searchSubdirs, FILE* msgStream ) {
  csVector<std::string> fileList;
  unsigned count = 0;
  int textLength = text.size();

  DIR* dir = opendir( directory.c_str() );
  if( dir == NULL ) {
    if( msgStream ) {
      fprintf( msgStream, "Error opening %s\n", directory.c_str() );
    }
    return false;
  }
        
  struct dirent *entryPtr = NULL;
        
  entryPtr = readdir( dir );
  while( entryPtr != NULL ) {
    struct stat entryInfo;
    if( ( strncmp( entryPtr->d_name, ".", PATH_MAX ) == 0 ) || ( strncmp( entryPtr->d_name, "..", PATH_MAX ) == 0 ) ) {
      entryPtr = readdir( dir );
      continue;
    }
    std::string fullPath = directory;
    fullPath.append( "/" );
    fullPath.append( entryPtr->d_name );
 
    if( stat( fullPath.c_str(), &entryInfo ) == 0 ) {
      count++;
      if( S_ISDIR( entryInfo.st_mode ) ) { // Directory
        //        printf( "Directory %s/\n", fullPath.c_str() );
        if( searchSubdirs ) {
          if( !csFileUtils::retrieveFiles_internal( fullPath, text, fileList_out, selectOption, searchSubdirs, msgStream ) ) {
            return false;
          }
        }
      }
      else if( S_ISREG( entryInfo.st_mode ) ) { // regular file
        //        printf( "\tFile %s has %lld bytes\n", fullPath.c_str(), (csInt64_t)entryInfo.st_size );
        int fileLength = strlen(entryPtr->d_name);
        if( selectOption == csFileUtils::OPTION_EXTENSION ) {
          if( textLength == 0 ||
              (textLength <= fileLength &&
              !text.compare( fullPath.substr(fullPath.size()-textLength,textLength) ) ) ) {
            fileList_out->insertEnd( fullPath );
          }
        }
        else {
          if( textLength == 0 ||
              (textLength <= fileLength &&
               !text.compare( fullPath.substr(fullPath.size()-fileLength,textLength) ) ) ) {
            fileList_out->insertEnd( fullPath );
          }
        }
      }
      //      else if( S_ISLNK( entryInfo.st_mode ) ) { // symbolic link
      //  if( msgStream ) {
      //    fprintf( msgStream, "File is a symbolic link: %s\n", fullPath.c_str() );
      //  }
        /*
          char targetName[PATH_MAX + 1];
          if( readlink( fullPath.c_str(), targetName, PATH_MAX ) != -1 ) {
          printf( "\tLink %s -> %s\n", fullPath.c_str(), targetName );
          }
          else {
          printf( "\tLink %s -> invalid link\n", fullPath.c_str() );
          return false;
          }
        */
      //      }
    }
    else {
      if( msgStream ) {
        fprintf( msgStream, "Error retrieving status information on file %s\n", fullPath.c_str() );
      }
    }
    entryPtr = readdir( dir );
  }
    
  closedir( dir );
  return true;
}

//--------------------------------------------------------------------------------
//
//

bool csFileUtils::retrieveDirectories( std::string const& directory,
                                       cseis_geolib::csVector<std::string>* dirList_out,
                                       FILE* msgStream )
{
  bool searchSubdirs = false;
  int option = RETRIEVE_DIRS;
  std::string dummyText("");
  return csFileUtils::retrieveFileDir_internal( directory, dummyText, dirList_out, option, searchSubdirs, msgStream );
}

bool csFileUtils::retrieveFileDir_internal( std::string const& directory,
              std::string const& text,
              cseis_geolib::csVector<std::string>* fileList_out,
              int option,
              bool searchSubdirs,
              FILE* msgStream )
{

  csVector<std::string> fileList;
  unsigned count = 0;
  int textLength = text.size();

  DIR* dir = opendir( directory.c_str() );
  if( dir == NULL ) {
    if( msgStream ) {
      fprintf( msgStream, "Error opening %s\n", directory.c_str() );
    }
    return false;
  }

  struct dirent *entryPtr = NULL;

  entryPtr = readdir( dir );
  while( entryPtr != NULL ) {
    struct stat entryInfo;
    if( ( strncmp( entryPtr->d_name, ".", PATH_MAX ) == 0 ) || ( strncmp( entryPtr->d_name, "..", PATH_MAX ) == 0 ) ) {
      entryPtr = readdir( dir );
      continue;
    }
    std::string fullPath = directory;
    fullPath.append( "/" );
    fullPath.append( entryPtr->d_name );

    if( stat( fullPath.c_str(), &entryInfo ) == 0 ) {
      count++;
      int fileLength = strlen(entryPtr->d_name);
      if( textLength == 0 || textLength >= fileLength ) {
        if( S_ISDIR( entryInfo.st_mode ) ) { // Directory
          if( option == RETRIEVE_DIRS ) {
            fileList_out->insertEnd( fullPath );
          }
          else if( searchSubdirs ) {
            if( !csFileUtils::retrieveFileDir_internal( fullPath, text, fileList_out, option, searchSubdirs, msgStream ) ) {
              return false;
            }
          }
        }
        else if( S_ISREG( entryInfo.st_mode ) ) { // regular file
          if( option == RETRIEVE_DIRS ) {
          }
          else if( option == csFileUtils::OPTION_EXTENSION ) {
            if( !text.compare( fullPath.substr(fullPath.size()-textLength,textLength) ) ) {
              fileList_out->insertEnd( fullPath );
            }
          }
          else {
            if( !text.compare( fullPath.substr(fullPath.size()-fileLength,textLength) ) ) {
              fileList_out->insertEnd( fullPath );
            }
          }
        }
        /*        else if( S_ISLNK( entryInfo.st_mode ) ) { // symbolic link
          if( msgStream ) {
            fprintf( msgStream, "File is a symbolic link: %s\n", fullPath.c_str() );
          }
        }
        */
      }
      else {
        if( msgStream ) {
          fprintf( msgStream, "Error retrieving status information on file %s\n", fullPath.c_str() );
        }
      }
    }
    entryPtr = readdir( dir );
  }
    
  closedir( dir );
  return true;

}

bool csFileUtils::createDirectory( std::string const& directory ) {
  struct stat entryInfo;
  if( stat( directory.c_str(), &entryInfo ) == 0 ) {
    if( S_ISDIR( entryInfo.st_mode ) ) { // Directory
      return true;
    }
    else {
      return false;
    }
  }
  else {
    std::string command( "mkdir " );
    command.append( directory );
    int returnFlag = system( command.c_str() );
    if( returnFlag != 0 ) {
      return false;
    }
    return true;
  }
}
bool csFileUtils::rename( std::string const& filenameOld, std::string const& filenameNew ) {
  bool ret = false;
  struct stat entryInfo;
  if( stat( filenameOld.c_str(), &entryInfo ) == 0 ) {
    if( S_ISREG( entryInfo.st_mode ) || S_ISDIR( entryInfo.st_mode ) ) { // regular file or Directory
#ifdef PLATFORM_WINDOWS
      std::string command( "move " );
#else
      std::string command( "mv -f " );
#endif
      command.append( filenameOld + " " + filenameNew );
      int returnFlag = system( command.c_str() );
      if( returnFlag == 0 ) {
        ret = true;
      }
    }
  }
  return ret;
}
bool csFileUtils::removeFile( std::string const& filename ) {
  bool ret = false;
  struct stat entryInfo;
  if( stat( filename.c_str(), &entryInfo ) == 0 ) {
    if( S_ISREG( entryInfo.st_mode ) ) { // regular file
#ifdef PLATFORM_WINDOWS
      std::string command( "del " );
#else
      std::string command( "rm -f " );
#endif
      command.append( filename );
      int returnFlag = system( command.c_str() );
      if( returnFlag == 0 ) {
        ret = true;
      }
    }
  }
  return ret;
}
bool csFileUtils::seekg_relative( csInt64_t bytePosRelative, std::ifstream* file ) {

  // Complex algorithm to be able to make step that is larger than 2Gb : Make several smaller steps instead
  
  int maxInt = std::numeric_limits<int>::max() - 1;   // -1 to be on the safe side, also for negative byte positions
  if( bytePosRelative < 0 ) maxInt *= -1;
  
  csInt64_t numSteps  = bytePosRelative / (csInt64_t)maxInt + 1LL;
  int bytePosResidual = (int)(bytePosRelative % (csInt64_t)maxInt);
  
  for( csInt64_t istep = 0; istep < numSteps - 1; ++istep ) { 
    file->clear(); // Clear all flags
    file->seekg( maxInt, std::ios_base::cur );
    if( file->fail() ) return false;
  }
  
  file->seekg( bytePosResidual, std::ios_base::cur );
  
  return true;
}

bool csFileUtils::fileExists( std::string const& filename ) {
  FILE* ftest = fopen( filename.c_str(), "r" );
  if( ftest == NULL ) { // File does not exist
    return false;
  }
  fclose(ftest);
  return true;
}
bool csFileUtils::createDoNotOverwrite( std::string const& filename ) {
  FILE* ftest = fopen( filename.c_str(), "r" );
  if( ftest != NULL ) { // File already exists
    fclose(ftest);
    return true;
  }
  ftest = fopen( filename.c_str(), "w" ); // Create file
  if( ftest == NULL ) return false;
  fclose(ftest);
  return true;
}

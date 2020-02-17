/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_FILE_UTILS_H
#define CS_FILE_UTILS_H

#include <string>
#include <cstdio>
#include "geolib_defines.h"

namespace cseis_geolib {

  template <typename T> class csVector;

  /**
   * File utilities.
   *
   * Accessing, listing, finding files etc
   * @author Bjorn Olofsson
   * @date 2007
   */
  class csFileUtils {
  public:
    static csInt64_t const FILESIZE_UNKNOWN = -111LL;
    static int const OPTION_START_STRING = 1;
    static int const OPTION_EXTENSION    = 2;

  public:
    csFileUtils();
    ~csFileUtils();
    static csInt64_t retrieveFileSize( std::string const& filename );
    static bool retrieveFileInfo( std::string const& filename, csInt64_t* fileSize, int* timeStamp_s );
    /**
     * Retrieve file's time stamp in UNIX seconds (seconds since 1-1-1970)
     */
    static bool retrieveTimeStamp( std::string const& filename, int* timeStamp_s );

    static bool retrieveDirectories( std::string const& directory,
             cseis_geolib::csVector<std::string>* dirList_out,
             FILE* msgStream = NULL );

    /**
     * Retrieve files in directory
     * Retrieve all files in the specified directory and its subdirectories with matching file extension.
     * @param directory     (i) Directory to scan
     * @param extension     (o) File 'extension'. Last N characters in filename. Only files with matching extension will be returned.
     *                          Pass empty string "" if all file names shall be returned.
     * @param fileList_out  (o) List of file names including full path
     * @param searchDir     (i) True if subdirectories shall be searched, false if not
     * @return false if some problem occurred during the operation (Proper exception handling not implemented yet)
     */
    static bool retrieveFiles( std::string const& directory,
             std::string const& extension,
             cseis_geolib::csVector<std::string>* fileList_out,
             bool searchSubdirs,
             FILE* msgStream = NULL );
    static bool retrieveFilesStartingWith( std::string const& directory,
             std::string const& startString,
             cseis_geolib::csVector<std::string>* fileList_out,
             bool searchSubdirs,
             FILE* msgStream = NULL );
    static bool createDirectory( std::string const& directory );
    static bool removeFile( std::string const& filename );
    static bool rename( std::string const& filenameOld, std::string const& filenameNew );

    static bool seekg_relative( csInt64_t bytePosRelative, std::ifstream* file );

    /**
     * Create file
     * @return true if file exists or was created
     */
    static bool createDoNotOverwrite( std::string const& filename );
    static bool fileExists( std::string const& filename );
  private:
    static int const RETRIEVE_FILES_EXTENSION = 10;
    static int const RETRIEVE_FILES_FIRST     = 11;
    static int const RETRIEVE_DIRS            = 12;


//    std::string myFilename;
//    csInt64_t myFileSize;
    int myRetrieveOption;

    static csInt64_t determineFileSize( std::string const& filename );
    static csInt64_t determineFileSize_bruteForceMethod( std::string const& filename );
    static bool retrieveFiles_internal( std::string const& directory,
          std::string const& text,
          cseis_geolib::csVector<std::string>* fileList_out,
          int selectOption,
          bool searchSubdirs,
          FILE* msgStream );

    static bool retrieveFileDir_internal( std::string const& directory,
            std::string const& text,
            cseis_geolib::csVector<std::string>* list_out,
            int option,
            bool searchSubdirs,
            FILE* msgStream );
  };

} // namespace
#endif


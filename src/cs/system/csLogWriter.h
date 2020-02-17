/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_LOG_WRITER_H
#define CS_LOG_WRITER_H

#include <cstdio>
#include <string>

namespace cseis_system {

/**
* Log file writer
*
* THE log writer for Cseis processing flows.
*
* If no log file is required, pass NULL to constructor.
* The methods will not do anything if no filename was given, or the file did not open correctly
* This is to ensure that Cseis flows will not fail even if no log file is available.
*
* @author Bjorn Olofsson
* @date   2007
*/
class csLogWriter {
public:
  /**
   * Constructor
   * Pass NULL if no log information shall be output.
   * @param filename File name with full path where log information shall be written to
   */
  csLogWriter( char const* filename );
  /**
   * Constructor
   * Pass NULL if no log information shall be output.
   * @param file File stream where log information shall be written to
   */
  csLogWriter( std::FILE* file );
  /**
   * Default constructor
   * Use this constructor to output to standard output
   */
  csLogWriter();
  ~csLogWriter();
  /**
  * Print message line to log file, with terminating newline, similar to printf("...\n")
  */
  void line( char const* text, ... );
  /**
  * Print message line to log file, without terminating newline, similar to printf("...")
  */
  void write( char const* text, ... );
  /**
  * Write error message to log file
  * The message is formatted in a standard way to inform users that an error occurred
  */
  void error( char const* text, ... );
  /**
  * Write error message to log file
  * The message is formatted in a standard way to inform users that an error occurred,
  * including a message about the module specified in the argument list
  * @param moduleName
  * @param moduleIndex
  */
  void error( char const* moduleName, int moduleIndex, char const* text, ... );
  /**
  * Write warning message to log file
  * The message is formatted in a standard way to inform users that a warning occurred
  */
  void warning( char const* text, ... );
  /**
  * Flush output
  */
  void flush();
  /**
  * @return file object
  */
  FILE* getFile() {
    return myLogFile;
  }
  /**
   * @return File name
   */
  char const* getFilename() const;
 private:
  csLogWriter( csLogWriter const& obj );
  std::FILE* myLogFile;
  char const* myFilename;
  /// true if file is opened within class. If yes then the file srtream will be closed during object destruction
  bool myIsOpenedLocally;
};
  
} // namespace
#endif


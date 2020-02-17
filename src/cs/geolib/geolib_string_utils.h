/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

/**
 * String/char utilities
 *
 * @author Bjorn Olofsson
 * @date 2007
 */

#ifndef CSEIS_STRING_UTILS_H
#define CSEIS_STRING_UTILS_H

#include <string>   // C++ 'string' class

namespace {
  //-------------------------------------
  // System methods
  static const char DOUBLE_QUOTE = '"';
  static const char TAB = 9;
  static const char WHITE_SPACE = ' ';
//  static const char BACK_SLASH  = '\\';
  static const char BACK_SLASH  = 92;
  static const char FORWARD_SLASH  = '/';
  static const char LETTER_COMMENT = '#';
  static const char NEW_LINE = 10;
  static const char CARRIAGE_RETURN = 13;
  static const char NULL_TERMINATOR = '\0';
  static const char COMMA = ',';
  static const char SEMICOLON = ';';
  static const char COLON = ':';

  static char const GT_CHAR = '>';
  static char const LT_CHAR = '<';
  static char const EQ_CHAR = '=';
  static char const DASH_CHAR = '-';
  static char const EXCL_CHAR = '!';
  static char const PLUS_CHAR = '+';
  static char const BRACKET_OPEN_CHAR = '(';
  static char const BRACKET_CLOSE_CHAR = ')';
  static char const BRACKET_SQ_OPEN_CHAR = '[';
  static char const BRACKET_SQ_CLOSE_CHAR = ']';
  static char const STAR_CHAR = '*';

}

namespace cseis_geolib {
  // Forward declarations:
  template <typename T> class csVector;

  /**
   * Remove outer white spaces and tabs
   */
  char* removeOuterSpaces( char* str );
  /**
   * Determine first non-blank character in given string
   * @return false if none is found (empty or fully blank string)
   */
  bool firstNonBlankChar( char const* str, char& c );
  /**
   * Tokenize input string, separators: white space, tab, "" enclosed (except /")
   */
  void tokenize( char const* line, csVector<std::string>& tokenList, bool removeComments = true );
  /**
   * Convert to upper case
   */
  std::string toUpperCase( std::string const& text );
  /**
   * Convert to upper case
   */
  void toUpperCase( char* text );
  /**
   * Convert to lower case
   */
  std::string toLowerCase( std::string const& text );
  /**
   * Convert to lower case
   */
  void toLowerCase( char* text );
  /**
   * @return true if given char is a letter
   */
  bool isLetter( char c );
  /**
   * @return true if given char is a digit
   */
  bool isDigit( char c );
  /**
   * Remove given character from given string
   */
  std::string removeChar( std::string const& str, char c );
  /**
   * Replace 'old' substring with 'new' substring in given text string
   */
  std::string replaceStr( std::string const& text, std::string const& strOld, std::string const& strNew );
  /**
   * Remove leading and trailing whitespace (including tabs and newlines)
   */
  std::string trim( std::string const& str );

} // namespace

#endif



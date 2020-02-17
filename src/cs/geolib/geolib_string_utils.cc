/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "geolib_string_utils.h"
#include "csVector.h"
#include <cstring>

namespace cseis_geolib {

void toUpperCase( char* text ) {
  int length = (int)strlen(text);
  for( int i = 0; i < length; i++ ) {
    if( text[i] >= 97 && text[i] <= 122 ) {
      text[i] = (char)(text[i] - 32);
    }
  }
}
void toLowerCase( char* text ) {
  int length = (int)strlen(text);
  for( int i = 0; i < length; i++ ) {
    if( text[i] >=65 && text[i] <= 90 ) {
      text[i] = (char)(text[i] + 32);
    }
  }
}

std::string toUpperCase( std::string const& text ) {
  std::string out = text;
  int length = (int)text.length();
  for( int i = 0; i < length; i++ ) {
    if( text[i] >= 97 && text[i] <= 122 ) {
      out[i] = (char)(out[i] - 32);
    }
  }
  return std::string( out );
}

std::string toLowerCase( std::string const& text ) {
  std::string out = text;
  int length = (int)text.length();
  for( int i = 0; i < length; i++ ) {
    if( text[i] >=65 && text[i] <= 90 ) {
      out[i] = (char)(out[i] + 32);
    }
  }
  return std::string( out );
}

bool isLetter( char c ) {
  return ( ( c >= 65 && c <= 90 ) || ( c >= 97 && c <= 122 ) || c == '_' );
}

bool isDigit( char c ) {
  return ( c >= 48 && c <= 57 );
}

std::string removeChar( std::string const& str, char c ) {
  char* out = new char[str.length()+1];
  int pos = 0;
  int length = (int)str.length();
  for( int i = 0; i < length; i++ ) {
    if( str[i] != c ) {
      out[pos++] = str[i];
    }
  }
  out[pos] = '\0';
  std::string outstr(out);
  delete [] out;
  return( outstr );
}

//-----------------------------------------------------------------
void tokenize( char const* line, csVector<std::string>& tokenList, bool removeComments ) {
  std::string str(line);
  int sizeIn = (int)str.length();  // strlen(line);
  int pos1, ipos;
  std::string tmp;
  ipos = sizeIn-1;
  // Reduce input string by ending newline or carriage returns
  while( ipos >= 0 && (line[ipos] == NEW_LINE || line[ipos] == CARRIAGE_RETURN) ) ipos--;
  sizeIn = ipos + 1;
  ipos = 0;
  do {
    while( ipos < sizeIn && (line[ipos] == WHITE_SPACE || line[ipos] == TAB || line[ipos] == COMMA || line[ipos] == SEMICOLON ) ) ipos++;
    pos1 = ipos;
    if( ipos == sizeIn ) break;
    if( line[ipos] == DOUBLE_QUOTE ) {
      ipos++;
      pos1 = ipos;
      while( ipos < sizeIn && (line[ipos] != DOUBLE_QUOTE || line[ipos-1] == BACK_SLASH) ) ipos++;
      tmp = str.substr( pos1, ipos-pos1 );
      ipos++;
      tokenList.insertEnd( tmp );
    }
    else {
      if( removeComments && line[ipos] == LETTER_COMMENT ) return;  // Remaining text are comments
      pos1 = ipos;
      while( ipos < sizeIn && 
             ( (line[ipos] != WHITE_SPACE && line[ipos] != TAB) || line[ipos] == COMMA || line[ipos] == SEMICOLON ) ) ipos++;
      tmp = str.substr( pos1, ipos-pos1 );
      tokenList.insertEnd( tmp );
      // tokenList.push( toLowerCase(tmp) );   // Convert to lower case unless in double quotes..
    }
  } while( ipos < sizeIn );
}

//-----------------------------------------------------------------
// Input string must be terminated by '/0'
//
char* removeOuterSpaces( char* str ) {
  int pos;
  char* strPtr = str;
  if( !str ) return str;

  while( strPtr && (*strPtr == ' ' || *strPtr == TAB) ) strPtr++;
  pos = (int)strlen(strPtr);
  while( pos > 0 && (strPtr[pos-1] == ' ' || strPtr[pos-1] == TAB || strPtr[pos-1] == '\n' ) ) pos--;
  strPtr[pos] = '\0';
  return strPtr;
}
//---------------------------------------------------------------------
bool firstNonBlankChar( char const* str, char& c ) {
  if( !str ) return false;

  char const* strPtr = str;
  while( strPtr && (*strPtr == WHITE_SPACE || *strPtr == TAB || *strPtr == NEW_LINE || *strPtr == CARRIAGE_RETURN ) ) strPtr++;
  if( strPtr == NULL || *strPtr == NULL_TERMINATOR ) return false;

  c = *strPtr;
  return true;
}

//---------------------------------------------------------------------
//

std::string replaceStr( std::string const& textIn, std::string const& strOld, std::string const& strNew ) {
  std::string textOut;
  int length = (int)textIn.length();
  int lenOld = (int)strOld.length();
  int counterIn = 0;
  int pos1 = 0;

  if( lenOld == 0 ) return std::string( textIn );
  //  if( lenOld < lenNew ) fprintf(stderr,"ERROR in replaceStr: Unsupported action, lenNew > lenOld\n");
  //...for some reason, this case does not work...? But it should

  char first = strOld[0];
  while( counterIn < length ) {
    if( textIn[counterIn] == first && length-counterIn >= lenOld ) {
      if( textIn.compare( counterIn, lenOld, strOld, 0, lenOld ) == 0 ) {
        if( pos1 < counterIn ) textOut.append( textIn.substr( pos1, counterIn-pos1) );
        textOut.append( strNew );
        counterIn += lenOld;
        pos1 = counterIn;
      }
      else {
        counterIn += 1;
      }        
    }
    else {
      counterIn += 1;
    }
  }
  if( pos1 < counterIn ) textOut.append( textIn.substr( pos1, counterIn-pos1) );

  return textOut;
}

std::string trim( std::string const& str ) {
  std::string const& whitespace = " \t\n";
  unsigned int strBegin = str.find_first_not_of(whitespace);
  if( strBegin == std::string::npos ) return "";
  unsigned int strEnd = str.find_last_not_of(whitespace);
  unsigned int strRange = strEnd - strBegin + 1;
  return str.substr( strBegin, strRange );
}

} // namespace


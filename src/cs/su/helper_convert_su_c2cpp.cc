/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include <string>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include "csVector.h"
#include "csFileUtils.h"
 
using namespace std;
using namespace cseis_geolib;

void writeCSEISIncludes( FILE* fout );
void writeCInclude( FILE* fout, char const* buffer );
void writeNamespace( FILE* fout, char const* buffer, char const* moduleName );
void writeMain( FILE* fout, char const* moduleName );
void writeInit( FILE* fout, char const* buffer, char const* moduleName );
void writeEnd( FILE* fout, char const* buffer );
void writeCatch( FILE* fout, char const* buffer );

void findReplace( std::string& str, char const* textIn, char const* textOut );
//void findReplaceMany( std::string& str, char const* textIn, char const* textOut );
void findReplaceStart( std::string& str, char const* textIn, char const* textOut, char doNotStartChar );
void findReplaceEnd( std::string& str, char const* textIn, char const* textOut, char doNotEndChar );
void findReplace( std::string& str, char const* textIn, char const* textOut, char doNotStartChar, char doNotEndChar );

void findReplaceErrorMessage( std::string& str );
// Replace ambiguous 'pow(x,y)' statements
void findCorrectPow( std::string& str );
// Remove leading and trailing white spaces
std::string trim( std::string const& str );
void getWhiteSpace( char* whiteSpace, int numSpaces, char const* buffer );
void countCurlyBrackets( int* numCurlyBrackets, char const* buffer );
void findReplaceDecodeReflectors( std::string& str );
void findReplaceGetParVal( std::string& str );

static int const STATE_00_START    = 0;
static int const STATE_01_INCLUDE  = 1;
static int const STATE_02_SDOC_START = 2;
static int const STATE_03_SDOC     = 3;
static int const STATE_04_SDOC_END = 4;
static int const STATE_05_MAIN     = 5;
static int const STATE_06_INIT     = 6;
static int const STATE_07_BODY     = 7;
static int const STATE_08_END      = 8;


/*
  Helper program for CSEIS 'SU' module.
  Converts SU module C source code into CSEIS-compatible C++ source code

  How to compile:
export SRCDIR=/disk/sources/seaseis/seaseis_v2.04
export LIBDIR=/opt/seaseis/lib
cd $SRCDIR/src/cs/su
g++ helper_convert_su_c2cpp.cc -o helper_convert_su_c2cpp -I../geolib -Wl,-rpath,${LIBDIR} -L${LIBDIR} -lgeolib

  How to use:
./helper_convert_su_c2cpp  $CWPROOT/src/su/main/filters  $SRCDIR/src/cs/su/main

Next, enter the new modules into the list of available SU modules in file $SRCDIR/src/su/su_modules.txt.

 */


int main( int argc, char** argv ) {
  int numCurlyBrackets = 0;

  if( argc < 3 ) {
    fprintf(stderr,"Error. Missing arguments: input directory, output directory\n");
    exit(-1);
  }
  std::string dirIn(argv[1]);
  std::string dirOut(argv[2]);
  cseis_geolib::csVector<std::string> fileList;
  bool searchSubDirs = true;
  bool success = csFileUtils::retrieveFiles( dirIn, ".c", &fileList, searchSubDirs, NULL );

  csVector<std::string> moduleNameList;

  fprintf(stderr,"Number of SU module source files: %d\n", fileList.size());
  for( int ifile = 0; ifile < fileList.size(); ifile++ ) {
    std::string filenameIn = fileList.at(ifile);
    int firstChar = -1;
    int lastChar  = -1;
    for( int i = filenameIn.length()-1; i > 0; i-- ) {
      if( lastChar < 0 && filenameIn[i] == '.' ) {
        lastChar = i-1;
      }
      else if( filenameIn[i] == '/' ) {
        firstChar = i+1;
        break;
      }
    }
    if( firstChar < 0 ) firstChar = 0;
    if( lastChar < 0 ) {
      fprintf(stderr,"Cannot find C source file extension '.c' in file name '%s'\n", filenameIn.c_str());
      continue;
    }
    std::string moduleNameStr = filenameIn.substr(firstChar,lastChar-firstChar+1);
    moduleNameList.insertEnd(moduleNameStr);
    char const* moduleName = moduleNameStr.c_str();
    string filenameOut(dirOut);
    filenameOut.append("/");
    filenameOut.append(moduleName);
    filenameOut.append(".cc");

    fprintf(stderr,"File #%-2d Convert SU module '%s' from %s to %s\n",
            ifile+1, moduleNameStr.c_str(), filenameIn.c_str(), filenameOut.c_str() );

  int state = STATE_00_START;
  char buffer[512];
  int lineNum = 0;

  FILE* fin  = fopen(filenameIn.c_str(),"r");
  FILE* fout = fopen(filenameOut.c_str(),"w");

  bool mainBracketFound = false;

  //  fprintf(stderr,"Input file:  %s\n", filenameIn.c_str());
  //  fprintf(stderr,"Output file: %s\n", filenameOut.c_str());

  while( fgets( buffer, 512, fin ) != NULL ) {
    string text(buffer);
    bool writeBuffer = true;
    int found = 0;

    if( state > STATE_05_MAIN && state < STATE_08_END ) {
      countCurlyBrackets( &numCurlyBrackets, buffer );
      if( numCurlyBrackets == 0 ) {  // End of main method reached --> Close try-catch block
        if( mainBracketFound ) {
          writeCatch( fout, buffer );
          state += 1;
        }
      }
      else {
        mainBracketFound  = true;
      }
    }
    //--------------------------------------------------
    if( state == STATE_07_BODY ) {
      found = text.find( "checkpars" );
      if( found >= 0 ) {
        findReplace(text,"checkpars","parObj.checkpars");
      }
      else {
        found = text.find( "requestdoc" );
        if( found >= 0 ) {
          writeBuffer = false;
        }
        else {
          found = text.find( "CWP_Exit" );
          int found2 = text.find( "EXIT_SUCCESS" );
          if( found >= 0 || found2 >= 0 ) {
            writeEnd( fout, buffer );
            writeBuffer = false;
          }

          findReplace(text,"fgettr(stdin,","cs2su->getTrace(");
          findReplace(text,"fgettr( stdin,","cs2su->getTrace(");
          findReplace(text,"fgettr( stdin ,","cs2su->getTrace(");
          findReplace(text,"fputtr(stdout,","cs2su->putTrace(");
          findReplace(text,"fputtr( stdout,","cs2su->putTrace(");
          findReplace(text,"fputtr( stdout ,","cs2su->putTrace(");

          findReplace(text,"gettr","cs2su->getTrace",'f','a');
          findReplaceStart(text,"puttr","su2cs->putTrace",'f');

          // Replace any getparval function
          findReplaceGetParVal( text );

          findReplace(text,"MUSTGETPAR","CSMUSTGETPAR");
          findReplace(text,"getpar","parObj.getpar");
          findReplace(text,"getnpar","parObj.getnpar");
          findReplace(text,"countpar","parObj.countpar");
          findReplace(text,"countnpar","parObj.countnpar");

          // Replace any SU error messages
          findReplaceErrorMessage( text );

          // Replace any decodeReflectors function
          findReplaceDecodeReflectors( text );

          // Correct any calls to pow()
          findCorrectPow( text );
        }
      }
    }
    else if( state == STATE_08_END ) {
      // Replace any SU error messages
      findReplaceErrorMessage( text );

      // Replace any decodeReflectors function
      findReplaceDecodeReflectors( text );

      // Replace any getparval function
      findReplaceGetParVal( text );

      // Correct any calls to pow()
      findCorrectPow( text );
    }

    //--------------------------------------------------
    if( state == STATE_00_START ) {
      found = text.find("#include");
      if( found >= 0 ) {
        writeCSEISIncludes(fout);
        writeCInclude(fout,"#include <pthread.h>");
        state += 1;
      }
    }
    //--------------------------------------------------
    if( state == STATE_01_INCLUDE ) {
      found = text.find("<");
      if( found >= 0 ) {
        writeCInclude(fout,buffer);
        writeBuffer = false;
      }
      else {
        found = text.find("self documentation");
        if( found >= 0 ) {
          state += 1;
        }
      }
    }
    //--------------------------------------------------
    else if( state == STATE_02_SDOC_START ) {
      found = text.find("*sdoc");
      if( found >= 0 ) {
        fprintf(fout,"std::string sdoc_%s =\n", moduleName);
        writeBuffer = false;
        state += 1;
      }
    }
    //--------------------------------------------------
    else if( state == STATE_03_SDOC ) {
      found = text.find("*sdoc");
      if( found >= 0 ) {
        fprintf(fout,"std::string sdoc_%s =\n", moduleName);
        writeBuffer = false;
        state += 1;
      }
      else {
        found = text.find( "NULL" );
        if( found >= 0 ) {
          found = text.find( "}" );
        }
        if( found >= 0 ) {
          fprintf(fout,";\n");
          writeNamespace( fout, buffer, moduleName );
          writeBuffer = false;
          state += 2;  // Skip STATE 04
        }
        else {
          found = text.find("\",");
          while( found >= 0 ) {
            findReplace(text,"\",","\"");
            found = text.find("\",");
          }
          found = text.find("\" ,");
          while( found >= 0 ) {
            findReplace(text,"\" ,","\"");
            found = text.find("\" ,");
          }
        }
      }
    }
    //--------------------------------------------------
    else if( state == STATE_04_SDOC_END ) {
      found = text.find("end self doc");
      if( found >= 0 ) {
        writeNamespace( fout, buffer, moduleName );
        state += 1;
        writeBuffer = false;
      }
    }
    //--------------------------------------------------
    else if( state == STATE_05_MAIN ) {
      found = text.find("main");
      if( found >= 0 ) found = text.find("argc");
      if( found >= 0 ) {
        writeMain( fout, moduleName );
        countCurlyBrackets( &numCurlyBrackets, buffer );
        if( numCurlyBrackets > 0 ) { // Opening bracket for main method is in same line
          fprintf(fout,"{\n");
        }
        state += 1;
        writeBuffer = false;
      }
      else {
        found = text.find("int");
        if( found >= 0 ) {
          std::string text2 = trim( text );
          if( text2.length() < 4 ) {
            writeBuffer = false;
          }
        }
      }
    }
    //--------------------------------------------------
    else if( state == STATE_06_INIT ) {
      found = text.find("initargs");
      if( found >= 0 ) {
        writeInit( fout, buffer, moduleName );
        state += 1;
        writeBuffer = false;
      }
      else {
        found = text.find( "requestdoc" );
        if( found >= 0 ) {
          writeBuffer = false;
        }
      }
    }

    if( writeBuffer ) {
      fprintf(fout,"%s", text.c_str());
    }
    lineNum += 1;
  }

  fprintf(fout,"\n} // END namespace\n");
  fclose(fout);
  fclose(fin);
  } // END for fileList

  for( int imod = 0; imod < moduleNameList.size(); imod++ ) {
    fprintf(stderr,"%s\n", moduleNameList.at(imod).c_str() );
  }
  return 0;
}

void writeCSEISIncludes( FILE* fout ) {
  fprintf(fout,"#include \"csException.h\"\n");
  fprintf(fout,"#include \"csSUTraceManager.h\"\n");
  fprintf(fout,"#include \"csSUArguments.h\"\n");
  fprintf(fout,"#include \"csSUGetPars.h\"\n");
  fprintf(fout,"#include \"su_complex_declarations.h\"\n");
  fprintf(fout,"#include \"cseis_sulib.h\"\n");
  fprintf(fout,"#include <string>\n");
  fprintf(fout,"\n");
}
void writeCInclude( FILE* fout, char const* buffer ) {
  fprintf(fout,"extern \"C\" {\n");
  fprintf(fout,"  %s\n", buffer);
  fprintf(fout,"}\n");
}
void writeNamespace( FILE* fout, char const* buffer, char const* moduleName ) {
  time_t timer = time(NULL);
  fprintf(fout,"\n/*  %s  Automatically modified for usage in SeaSeis  */\n", asctime(localtime(&timer)));

  //  fprintf(fout,"%s\n", buffer);
  fprintf(fout,"namespace %s {\n\n", moduleName);
}
void writeMain( FILE* fout, char const* moduleName ) {
  fprintf(fout,"void* main_%s( void* args )\n", moduleName);
}

void getWhiteSpace( char* whiteSpace, int numSpaces, char const* buffer ) {
  int counter = 0;
  int length = strlen(buffer);
  while( (counter < length && counter < numSpaces-1) && (buffer[counter] == ' ' || buffer[counter] == '\t' ) ) {
    whiteSpace[counter] = buffer[counter];
    counter += 1;
  }
  whiteSpace[counter] = '\0';
}

void writeInit( FILE* fout, char const* buffer, char const* moduleName ) {
  char whiteSpace[100];
  getWhiteSpace( whiteSpace, 100, buffer );

  fprintf(fout,"%scseis_su::csSUArguments* suArgs = (cseis_su::csSUArguments*)args;\n", whiteSpace);
  fprintf(fout,"%scseis_su::csSUTraceManager* cs2su = suArgs->cs2su;\n", whiteSpace);
  fprintf(fout,"%scseis_su::csSUTraceManager* su2cs = suArgs->su2cs;\n", whiteSpace);
  fprintf(fout,"%sint argc = suArgs->argc;\n", whiteSpace);
  fprintf(fout,"%schar **argv = suArgs->argv;\n", whiteSpace);
  fprintf(fout,"%scseis_su::csSUGetPars parObj;\n", whiteSpace);
  fprintf(fout,"\n");
  fprintf(fout,"%svoid* retPtr = NULL;  /*  Dummy pointer for return statement  */\n", whiteSpace);
  fprintf(fout,"%ssu2cs->setSUDoc( sdoc_%s );\n", whiteSpace, moduleName);
  fprintf(fout,"%sif( su2cs->isDocRequestOnly() ) return retPtr;\n",whiteSpace);
  fprintf(fout,"%sparObj.initargs(argc, argv);\n", whiteSpace);
  fprintf(fout,"\n");
  fprintf(fout,"%stry {  /* Try-catch block encompassing the main function body */\n", whiteSpace);
  fprintf(fout,"\n");
}

void writeEnd( FILE* fout, char const* buffer ) {
  char whiteSpace[100];
  getWhiteSpace( whiteSpace, 100, buffer );

  fprintf(fout,"%ssu2cs->setEOF();\n",whiteSpace);
  fprintf(fout,"%spthread_exit(NULL);\n",whiteSpace);
  fprintf(fout,"%sreturn retPtr;\n",whiteSpace);
}

void writeCatch( FILE* fout, char const* buffer ) {
  char whiteSpace[100];
  getWhiteSpace( whiteSpace, 100, buffer );

  fprintf(fout,"%s}\n",whiteSpace);
  fprintf(fout,"%scatch( cseis_geolib::csException& exc ) {\n",whiteSpace);
  fprintf(fout,"%s  su2cs->setError(\"%%s\",exc.getMessage());\n",whiteSpace);
  fprintf(fout,"%s  pthread_exit(NULL);\n",whiteSpace);
  fprintf(fout,"%s  return retPtr;\n",whiteSpace);
  fprintf(fout,"%s}\n",whiteSpace);
}

void findReplace( std::string& str, char const* textIn, char const* textOut, char doNotStartChar, char doNotEndChar ) {
  int found  = str.find( textIn );
  int lengthIn  = (int)strlen(textIn);
  int lengthOut = (int)strlen(textOut);
  while( found >= 0 && found < str.length() ) {
    if( found > 0 && str[found-1] == doNotStartChar ) {
      // nothing
    }
    else if( found+lengthIn < ((int)str.length()-2) && str[found+lengthIn] == doNotEndChar ) {
      // nothing
    }
    else {
      str.erase(found,lengthIn);
      str.insert(found,textOut);
    }
    found = str.find( textIn, found+lengthOut );
  }
}
void findReplaceStart( std::string& str, char const* textIn, char const* textOut, char doNotStartChar ) {
  findReplace( str, textIn, textOut, doNotStartChar, '@' );
}
void findReplaceEnd( std::string& str, char const* textIn, char const* textOut, char doNotEndChar ) {
  findReplace( str, textIn, textOut, '@', doNotEndChar );
}
//void findReplaceMany( std::string& str, char const* textIn, char const* textOut ) {
//  findReplace( str, textIn, textOut, '@, '@' );
//}
void findReplace( std::string& str, char const* textIn, char const* textOut ) {
  //  findReplace( str, textIn, textOut, '@', '@' );
  int found  = str.find( textIn );
  if( found >= 0 ) {
    str.erase(found,strlen(textIn));
    str.insert(found,textOut);
  }
}

std::string trim( std::string const& str ) {
  std::string const& whitespace = " \t\n";
  int strBegin = str.find_first_not_of(whitespace);
  if( strBegin == std::string::npos ) return "";
  int strEnd = str.find_last_not_of(whitespace);
  int strRange = strEnd - strBegin + 1;
  return str.substr( strBegin, strRange );
}

void findReplaceErrorMessage( std::string& text ) {
  string textErr("err(");
  int found = text.find( textErr.c_str() );
  if( found < 0 ) {
    found = text.find( "err (" );
    if( found < 0 ) return;
    textErr = "err (";
  }
  // Make sure this is not a call to 'clearerr(':
  if( found > 4 && !( (text.substr(found-5,8)).compare("clearerr")) ) {
    fprintf(stderr,"CLEAR: %s\n", text.c_str() );
    return;
  }
  findReplace( text, textErr.c_str(), "throw cseis_geolib::csException(");
}

void findCorrectPow( std::string& text ) {
  int startPos = 0;
  // Continue searching until end of string reached, or no function pow found
  while( startPos < text.length()-2 ) {
    int found = text.find( "pow(", startPos );
    if( found < 0 ) return;

    int length = text.length();
    int counter = found+4;  // Skip 4 letters = pow(
    int counterOpenBrackets = 0;
    // Step 1: Search for ending double quotes
    while( counter < length ) {
      if( text[counter] == '(' ) {
        counterOpenBrackets += 1;
      }
      else if( text[counter] == ')' ) {
        if( counterOpenBrackets == 0 ) break;
        counterOpenBrackets -= 1;
      }
      counter += 1;
    }
    if( counter == length ) return; // Ending brackets of pow function not found. Why??
    
    if( text[counter-1] == '2' || text[counter-1] == '4' ) {
      text.insert(counter,".0");
    }
    if( text[counter-1] == ' ' && text[counter-2] == '2' ) {
      text.insert(counter-1,".0");
    }
    startPos = counter;
  }
  
}
void countCurlyBrackets( int* numCurlyBrackets, char const* buffer ) {
  int length = strlen(buffer);
  for( int i = 0; i < length; i++ ) {
    if( buffer[i] == '{' ) *numCurlyBrackets += 1;
    else if( buffer[i] == '}' ) *numCurlyBrackets -= 1;
  }
}

void findReplaceDecodeReflectors( std::string& text ) {
  string textFind("decodeReflectors(");
  int found = text.find( textFind.c_str() );
  findReplace( text, textFind.c_str(), "cseis_decodeReflectors(&parObj,");
}

void findReplaceGetParVal( std::string& text ) {
  string textFind("getparval(");
  int found = text.find( textFind.c_str() );
  findReplace( text, textFind.c_str(), "cseis_getParVal(&parObj,");
}

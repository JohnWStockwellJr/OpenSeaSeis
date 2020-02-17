/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <iostream>
#include <string>
#include <cstring>
#include "csHelp.h"
#include "cseis_defines.h"

/// !CHANGE! these methods, maybe put in different file or class...

using namespace cseis_system;
using namespace std;

/**
 * CSEIS Help
 *
 * Command line tool for CSEIS help
 *
 * @author Bjorn Olofsson
 */
int main( int argc, char** argv ) {
  char** argPtr;
  string moduleName;

  //--------------------------------------------------------------------------------
  // Command line options
  //

  argPtr = argv;
  if( argc == 1 ) {
    fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
    return(-1) ;
  }
  int iArg = 1;   // Start with 1 to skip name of program executable
  while( iArg < argc ) {
    if( argv[iArg][0] == '-' ) {
      if( strlen(argv[iArg]) < 2 ) {
        fprintf(stderr," Syntax error in command line: Unknown option: '%s'\n", argv[iArg]);
        fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
        return(-1);
      }
      char option = argv[iArg][1];
      if ( option == 'h' ) {
        if( !strcmp( argv[iArg], "-html" ) ) {
          csHelp help( stdout );
          help.fullHtmlHelp();
          return(-1);
        }
        fprintf( stderr, " SeaSeis help tool.\n");
        fprintf( stderr, " Usage:  %s [-h] [-m <name>] [-v] [-std]\n", argv[0] );
        fprintf( stderr, " -m <name>              : Print help for module <name>.\n");
        fprintf( stderr, "                        : If no module is found, print help of all modules starting with <name>.\n");
        fprintf( stderr, "                        : If <name> = empty: Print full help for all available modules.\n");
        fprintf( stderr, "                        : If <name> = .: Print short help for all available modules.\n");
        fprintf( stderr, " -h                     : Print this page\n");
        fprintf( stderr, " -html                  : Print full help for all modules & standard trace headers in HTML\n");
        fprintf( stderr, " -v                     : Print out version info\n");
        fprintf( stderr, " -std                   : Dump all standard trace headers\n");
        return(-1);
      }
      else if ( option == 'v' ) {
        fprintf( stderr, " SeaSeis version %s, build date %s\n", CSEIS_VERSION, __DATE__);
        fprintf( stderr, " Bjorn Olofsson, 2006-2012\n");
        return(-1);
      }
      else if ( option == 's' ) {
        if( !strcmp( argv[iArg], "-std" ) ) {
          csHelp help( stdout );
          help.standardHeaderHelp();
          return(-1);
        }
        else {
          fprintf(stderr," Syntax error in command line: Unknown option: '%s'\n", argv[iArg]);
          fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
          return(-1) ;
        }
        ++iArg;
      }
      else if ( option == 'm' ) {
        ++iArg;
        if( iArg == argc ) {
          moduleName = ",";
        }
        else {
          moduleName = argv[iArg];
        }
        csHelp help( stdout );
        help.moduleTextHelp( moduleName );
        return(-1);
      }
      else {
        fprintf(stderr," Syntax error in command line: Unknown option: '%s'\n", argv[iArg]);
        fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
        return(-1) ;
      }
    }
    else {
      fprintf(stderr," Syntax error in command line:\n Expected command line option (starting with '-'), found '%s'\n", argv[iArg]);
      fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
      return(-1) ;
    }
  }
}



/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSUArguments.h"
#include "csSUTraceManager.h"
#include "csVector.h"
#include <cstdio>
#include <cstring>

using namespace cseis_su;

csSUArguments::csSUArguments() {
  argv = NULL;
  argc = 0;
  debugFlag = 0;
}

csSUArguments::~csSUArguments() {
  if( argv != NULL ) {
    for( int i = 0; i < argc; i++ ) {
      if( argv[i] != NULL ) {
        delete [] argv[i];
      }
    }
    delete [] argv;
    argv = NULL;
  }
}

void csSUArguments::setArgv( cseis_geolib::csVector<std::string>* argvList ) {
  argc = argvList->size();
  argv = new char*[argc];
  for( int i = 0; i < argvList->size(); i++ ) {
    std::string text = argvList->at(i);
    int length = (int)text.length();
    argv[i] = new char[length+1];
    memcpy(argv[i],text.c_str(),length);
    argv[i][length] = '\0';
  }
}

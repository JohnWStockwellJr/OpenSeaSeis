/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include <string>
#include <cstdlib>
#include <cstring>
#include "csVector.h"
 
using namespace std;
using namespace cseis_geolib;

void writeInclude( FILE* fout, cseis_geolib::csVector<std::string>* moduleList );
void writeMake( FILE* fout, cseis_geolib::csVector<std::string>* moduleList );

/*
  Helper program for CSEIS 'SU' module.
  Creates a make file (Makefile_su) and header file (su_modules.h), both required to compile CSEIS 'SU' module (mod_su.cc).
  Reads in a file called "su_modules.txt" containing a list of available SU modules.

  How to compile and use:

export SRCDIR=/disk/sources/seaseis/seaseis_v2.04
export LIBDIR=/opt/seaseis/lib
cd $SRCDIR/src/cs/su
g++ helper_create_su_make_include.cc -o helper_create_su_make_include -I../geolib -Wl,-rpath,${LIBDIR} -L${LIBDIR} -lgeolib
./helper_create_su_make_include < su_modules.txt


  To build Seaseis, run the make script make_seaseis.sh from the directory $SRCDIR.
*/


int main( int argc, char** argv ) {
  char buffer[512];
  int lineNum = 0;

  string filenameIncl("su_modules.h");
  string filenameMake("Makefile_su");

  fprintf(stderr,"Include file:  %s\n", filenameIncl.c_str());
  fprintf(stderr,"Make file:     %s\n", filenameMake.c_str());

  csVector<string> moduleList;

  while( fgets( buffer, 512, stdin ) != NULL ) {
    lineNum += 1;
    if( buffer[0] == '#' ) continue;
    if( strlen(buffer) < 3 ) continue;
    if( lineNum == 1 ) continue;
    string text(buffer);
    text = text.substr(0,text.length()-1);
    moduleList.insertEnd( text );
  }

  FILE* fincl = fopen(filenameIncl.c_str(),"w");
  writeInclude(fincl, &moduleList);
  fclose(fincl);

  FILE* fmake = fopen(filenameMake.c_str(),"w");
  writeMake(fmake, &moduleList);
  fclose(fmake);

  return 0;
}

void writeInclude( FILE* fout, cseis_geolib::csVector<std::string>* moduleList ) {
  int numModules = moduleList->size();
  fprintf(fout,"#ifndef SU_MODULES_H\n");
  fprintf(fout,"#define SU_MODULES_H\n");
  fprintf(fout,"\n");
  fprintf(fout,"typedef void* (*SUMainMethod) ( void* );\n");
  fprintf(fout,"\n");
  fprintf(fout,"char *sdoc[] = {\n");
  fprintf(fout,"\"                                                       \",\n");
  fprintf(fout,"\" GENERIC SU self-documentation                         \",\n");
  fprintf(fout,"\"                                                       \",\n");
  fprintf(fout,"NULL};\n");
  fprintf(fout,"\n");
  fprintf(fout,"const int NUM_SU_MODULES = %d;\n", numModules);
  fprintf(fout,"\n");
  fprintf(fout,"\n");

  fprintf(fout,"std::string SU_MODULE_NAMES[NUM_SU_MODULES] = {\n");
  for( int imod = 0; imod < numModules; imod++ ) {
    if( imod > 0 ) fprintf(fout,",\n");
    fprintf(fout, "  std::string(\"%s\")", moduleList->at(imod).c_str() );
  }
  fprintf(fout,"\n};\n");

  for( int imod = 0; imod < numModules; imod++ ) {
    fprintf(fout,"namespace %s {\n", moduleList->at(imod).c_str() );
    fprintf(fout,"  void* main_%s( void* args );\n", moduleList->at(imod).c_str() );
    fprintf(fout,"}\n");
  }
  fprintf(fout,"\n");

  fprintf(fout,"void*(*SU_MAIN_METHODS[NUM_SU_MODULES])( void* args ) = {\n");
  for( int imod = 0; imod < numModules; imod++ ) {
    if( imod > 0 ) fprintf(fout,",\n");
    fprintf(fout, "  %s::main_%s", moduleList->at(imod).c_str(), moduleList->at(imod).c_str() );
  }
  fprintf(fout,"\n};\n");
  fprintf(fout,"#endif\n");
}


void writeMake( FILE* fout, cseis_geolib::csVector<std::string>* moduleList ) {
  int numModules = moduleList->size();
  fprintf(fout,"MODDIR    = $(SRCDIR)/cs/su\n");
  fprintf(fout,"LIB_SU    = libsumodules.so\n");

  fprintf(fout,"OBJS_SU =");
  for( int imod = 0; imod < numModules; imod++ ) {
    fprintf(fout, " $(OBJDIR)/%s.o", moduleList->at(imod).c_str() );
  }
  fprintf(fout,"\n");

  fprintf(fout,"\n");
  fprintf(fout,"INCS =  -I\"$(SRCDIR)/cs/geolib\" -I\"$(MODDIR)\" -I\"$(CWPROOT)/include\" -I\"$(CWPROOT)/src/Complex/include\"\n");
  fprintf(fout,"\n");
  fprintf(fout,"# -Wno-write-strings: Suppress warnings about deprecated string assignments\n");
  fprintf(fout,"# -Wno-sign-compare:  Suppress warnings about comparison between signed and unsigned integer expressions\n");
  fprintf(fout,"# -fpermissive:       Downgrade errors to warnings for certain language structures\n");
  fprintf(fout,"ALL_FLAGS = $(INCS) $(COMMON_FLAGS) -fPIC -g -Wno-write-strings -Wno-sign-compare -fpermissive\n");
  fprintf(fout,"\n");
  fprintf(fout,"default: $(OBJS_SU) $(LIBDIR)/$(LIB_SU)\n");
  fprintf(fout,"\n");
  fprintf(fout,"clean:\n");
  fprintf(fout,"	rm -f $(OBJS_SU)\n");
  fprintf(fout,"\n");
  fprintf(fout,"bleach: clean\n");
  fprintf(fout,"	rm -f $(LIBDIR)/$(LIB_SU)\n");
  fprintf(fout,"\n");
  fprintf(fout,"$(LIBDIR)/$(LIB_SU): $(OBJS_SU)\n");
  fprintf(fout,"	$(CPP) -shared -Wl,-soname,$(LIB_SU) -o $(LIBDIR)/$(LIB_SU) $(OBJS_SU) -L$(LIBDIR) -L$(CWPROOT)/lib -lc -lgeolib -lsu -lpar -lcwp -lm -lpthread\n");
  fprintf(fout,"\n");
  fprintf(fout,"##################################################\n");
  fprintf(fout,"###### SU modules ##############\n");
  fprintf(fout,"\n");

  for( int imod = 0; imod < numModules; imod++ ) {
    char const* name = moduleList->at(imod).c_str() ;
    fprintf(fout,"$(OBJDIR)/%s.o: $(MODDIR)/main/%s.cc\n", name, name);
    fprintf(fout,"	$(CPP) -c $(ALL_FLAGS) $(MODDIR)/main/%s.cc -I$(CWPROOT)/include -Wall -pedantic -L$(CWPROOT)/lib -lsu -lpar -lcwp -lm -o $(OBJDIR)/%s.o\n", name, name);
    fprintf(fout,"\n");
  }
}

/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdlib>
#include <cstdio>
#include <string>
#include "cseis_defines.h"
#include "geolib_string_utils.h"
#include "csVector.h"
#include "csSort.h"

using namespace std;
using namespace cseis_geolib;
using namespace cseis_system;


namespace cseis_system {
  class csTrace;
  class csParamManager;
  class csTraceGather;
  class csParamDef;
  class csInitPhaseEnv;
  class csExecPhaseEnv;
  class csLogWriter;
}

namespace cseis {
  static string const filenameModuleList = "/include/cseis_modules_win.txt";

/**
 * CSEIS_BUILD for Windows
 *
 * @author Bjorn Olofsson
 * @date 2006
 */

void cseis_build( string const& srcdir, string const& bindir,
                  cseis_geolib::csVector<std::string> const& moduleFileList )
{
  string filenameOut( bindir + "/include/cseis_modules_win.h" );
  string filenameOutAll( bindir + "/include/cseis_modules_win_all.h" );
  FILE* fout_names = fopen( filenameOut.c_str(), "w" );
  if( fout_names == NULL ) {
    fprintf(stderr,"Could not open module header output file '%s'\n", filenameOut.c_str());
    exit(-1);
  }
  FILE* fout_all = fopen( filenameOutAll.c_str(), "w" );
  if( fout_all == NULL ) {
    fprintf(stderr,"Could not open module header output file '%s'\n", filenameOutAll.c_str());
    exit(-1);
  }
  char line[MAX_LINE_LENGTH];
  csVector<string> tokenList;
  csVector<string> singleModuleList;
  csVector<string> multiModuleList;

  for( int ifile = 0; ifile < moduleFileList.size(); ifile++ ) {
    string filenameIn = moduleFileList.at(ifile);
    fprintf(stderr,"Generate module header file from input file  %s\n", filenameIn.c_str());
    FILE* fin = fopen( filenameIn.c_str(), "r" );
    if( fin == NULL ) {
      fprintf(stderr,"Could not open module listing input file '%s'\n", filenameIn.c_str());
      exit(-1);
    }
    int counterLines = 0;

    while( fgets( line, MAX_LINE_LENGTH, fin ) != NULL ) {
      counterLines++;
      tokenList.clear();
      tokenize( line, tokenList );
      if( tokenList.size() >= 2 ) {
        string type = tokenList.at(1);
        string name = toLowerCase(tokenList.at(0));
        if( !type.compare("SINGLE") ||
            !type.compare("INPUT") ) {
          singleModuleList.insertEnd( name );
        }
        else if( !type.compare("MULTI") ) {
          multiModuleList.insertEnd( name );
        }
        else {
          fprintf(stderr,"Module type not recognized: '%s'. Input line: '%s'.\n", type.c_str(), line );
          exit(-1);
        }
      }
      else if( tokenList.size() != 0 ) {
        fprintf(stderr,"Syntax error in module listing input file '%s', line %d:\nExpected two columns, found %d.\n",
                filenameIn.c_str(), counterLines, tokenList.size() );


        exit(-1);
      }
    }
    fclose(fin);
  }

  if( singleModuleList.size() == 0 && multiModuleList.size() == 0 ) {
    fprintf(stderr,"No modules specified in module listing input file '%s'.\n", moduleFileList.at(0).c_str() );
    exit(-1);
  }
  int numSingleModules = singleModuleList.size();
  int numMultiModules  = multiModuleList.size();
  int numModules       = singleModuleList.size() + multiModuleList.size();

  csVector<string> allModuleList(numModules);
  for( int i = 0; i < numSingleModules; i++ ) {
    allModuleList.insertEnd( singleModuleList.at(i) );
  }
  for( int i = 0; i < numMultiModules; i++ ) {
    allModuleList.insertEnd( multiModuleList.at(i) );
  }


  // Sort module lists
  std::string* list = new std::string[allModuleList.size()];
  for( int i = 0; i < allModuleList.size(); i++ ) {
    list[i] = allModuleList.at(i);
  }
  csSort<string> sort;
  sort.treeSort( list, allModuleList.size() );
  for( int i = 0; i < allModuleList.size(); i++ ) {
    allModuleList.at(i) = list[i];
  }
  delete [] list;
  
  list = new std::string[singleModuleList.size()];
  for( int i = 0; i < singleModuleList.size(); i++ ) {
    list[i] = singleModuleList.at(i);
  }
  sort.treeSort( list, singleModuleList.size() );
  for( int i = 0; i < singleModuleList.size(); i++ ) {
    singleModuleList.at(i) = list[i];
  }
  delete [] list;
  
  list = new std::string[multiModuleList.size()];
  for( int i = 0; i < multiModuleList.size(); i++ ) {
    list[i] = multiModuleList.at(i);
  }
  sort.treeSort( list, multiModuleList.size() );
  for( int i = 0; i < multiModuleList.size(); i++ ) {
    multiModuleList.at(i) = list[i];
  }
  delete [] list;


  //-------------------------------------------------------
  // Create module header file
  //
  fprintf( fout_all, "#ifndef CS_MODULES_ALL_H\n" );
  fprintf( fout_all, "#define CS_MODULES_ALL_H\n" );
  fprintf( fout_all, "\n" );
  fprintf( fout_all, "namespace cseis_system {\n");
  fprintf( fout_all, "  class csTrace;\n");
  fprintf( fout_all, "  class csParamManager;\n");
  fprintf( fout_all, "  class csTraceGather;\n");
  fprintf( fout_all, "  class csParamDef;\n");
  fprintf( fout_all, "  class csInitPhaseEnv;\n");
  fprintf( fout_all, "  class csExecPhaseEnv;\n");
  fprintf( fout_all, "  class csLogWriter;\n");
  fprintf( fout_all, "}\n");
  fprintf( fout_all, "\n" );
  fprintf( fout_all, "const int N_METHODS_SINGLE = %d;\n", numSingleModules);
  fprintf( fout_all, "const int N_METHODS_MULTI  = %d;\n", numMultiModules);
  fprintf( fout_all, "const int N_METHODS = N_METHODS_MULTI + N_METHODS_SINGLE;\n", numModules);
  fprintf( fout_all, "\n" );


  for( int i = 0; i < numModules; i++ ) {
    fprintf( fout_all, "void init_mod_%s_(cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter*);\n",
      allModuleList.at(i).c_str() );
  }
  fprintf( fout_all, "\n" );


  fprintf( fout_all, "bool exec_mod_dummy_single_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) { return true; }\n");
  for( int i = 0; i < numSingleModules; i++ ) {
    fprintf( fout_all, "bool exec_mod_%s_(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);\n",
      singleModuleList.at(i).c_str() );
  }
  fprintf( fout_all, "\n" );
  fprintf( fout_all, "void exec_mod_dummy_multi_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) {}\n");
  for( int i = 0; i < numMultiModules; i++ ) {
    fprintf( fout_all, "void exec_mod_%s_(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*);\n",
      multiModuleList.at(i).c_str() );
  }
  fprintf( fout_all, "\n" );

  fprintf( fout_all, "void(*METHODS_INIT[N_METHODS])( cseis_system::csParamManager*, cseis_system::csInitPhaseEnv*, cseis_system::csLogWriter* ) = {\n");
  for( int i = 0; i < numSingleModules; i++ ) {
    if( i != 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  init_mod_%s_", singleModuleList.at(i).c_str() );
  }
  for( int i = 0; i < numMultiModules; i++ ) {
    if( i != 0 || numSingleModules > 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  init_mod_%s_", multiModuleList.at(i).c_str() );
  }
  fprintf( fout_all, "\n};\n" );
  fprintf( fout_all, "\n" );

  fprintf( fout_all, "bool(*METHODS_EXEC_SINGLE[N_METHODS])(cseis_system::csTrace*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) = {\n");
  for( int i = 0; i < numSingleModules; i++ ) {
    if( i != 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  exec_mod_%s_", singleModuleList.at(i).c_str() );
  }
  for( int i = 0; i < numMultiModules; i++ ) {
    if( i != 0 || numSingleModules > 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  exec_mod_dummy_single_" );
  }
  fprintf( fout_all, "\n};\n" );
  fprintf( fout_all, "\n" );

  fprintf( fout_all, "void(*METHODS_EXEC_MULTI[N_METHODS])(cseis_system::csTraceGather*, int*, int*, cseis_system::csExecPhaseEnv*, cseis_system::csLogWriter*) = {\n");
  for( int i = 0; i < numSingleModules; i++ ) {
    if( i != 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  exec_mod_dummy_multi_" );
  }
  for( int i = 0; i < numMultiModules; i++ ) {
    if( i != 0 || numSingleModules > 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  exec_mod_%s_", multiModuleList.at(i).c_str() );
  }
  fprintf( fout_all, "\n};\n" );
  fprintf( fout_all, "\n" );
  fprintf( fout_all, "std::string NAMES[N_METHODS] = {\n");


  for( int i = 0; i < numSingleModules; i++ ) {
    if( i != 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  std::string(\"%s\")", toUpperCase(singleModuleList.at(i)).c_str() );
  }
  for( int i = 0; i < numMultiModules; i++ ) {
    fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  std::string(\"%s\")", toUpperCase(multiModuleList.at(i)).c_str() );
  }

  fprintf( fout_all, "\n};\n" );
  fprintf( fout_all, "\n" );

  for( int i = 0; i < numModules; i++ ) {
    fprintf( fout_all, "void params_mod_%s_(cseis_system::csParamDef*);\n", allModuleList.at(i).c_str() );
  }
  fprintf( fout_all, "\n" );

  fprintf( fout_all, "void (*METHODS_PARAM[N_METHODS])(cseis_system::csParamDef*) = {\n");

  for( int i = 0; i < numSingleModules; i++ ) {
    if( i != 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  params_mod_%s_", singleModuleList.at(i).c_str() );
  }
  for( int i = 0; i < numMultiModules; i++ ) {
    if( i != 0 || numSingleModules > 0 ) fprintf( fout_all, ",\n" );
    fprintf( fout_all, "  params_mod_%s_", multiModuleList.at(i).c_str() );
  }
  fprintf( fout_all, "\n};\n" );
  fprintf( fout_all, "\n" );
  fprintf( fout_all, "#endif\n" );

  //-------------------------------------------------------------
  //
  fprintf( fout_names, "#ifndef CS_MODULES_H\n" );
  fprintf( fout_names, "#define CS_MODULES_H\n" );
  fprintf( fout_names, "\n" );
  fprintf( fout_names, "\n" );
  fprintf( fout_names, "const int N_METHODS_SINGLE = %d;\n", numSingleModules);
  fprintf( fout_names, "const int N_METHODS_MULTI  = %d;\n", numMultiModules);
  fprintf( fout_names, "const int N_METHODS = N_METHODS_MULTI + N_METHODS_SINGLE;\n", numModules);
  fprintf( fout_names, "\n" );
  fprintf( fout_names, "std::string NAMES[N_METHODS] = {\n");
  for( int i = 0; i < numModules; i++ ) {
    if( i != 0 ) fprintf( fout_names, ",\n" );
    fprintf( fout_names, "  std::string(\"%s\")", toUpperCase(allModuleList.at(i)).c_str() );
  }
  fprintf( fout_names, "\n};\n" );
  fprintf( fout_names, "\n" );
  fprintf( fout_names, "#endif\n" );

}

} // namespace

using namespace cseis;

int main(int argc, char *argv[])
{
  if( argc < 3 ) {
    fprintf(stderr,"Generate C++ header files containing all available Cseis modules.\n");
    fprintf(stderr,"This is required before compilation whenever the module list file is changed.\n");
    fprintf(stderr,"Usage: %s <source_path_name> [module_list_file...]\n", argv[0]);
    fprintf(stderr," <source_path_name>    : Path name to Cseis source code.\n");
    fprintf(stderr," <bin_path_name>       : Path name to Cseis binary output.\n");
    fprintf(stderr," [module_list_file...] : Space separated list of files containing a list of all modules. Syntax:\n");
    fprintf(stderr,"                         MODULE_NAME   <SINGLE|MULTI|INPUT>\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Default module list file: <source_path_name>%s\n\n", filenameModuleList.c_str());
    exit(-1);
  }
  string srcdir = argv[1];
  string bindir = argv[2];
  bool doWriteModuleHeaderFile = true;
  csVector<string> moduleFileList;
  int counter = 3;
  while( counter < argc ) {
    string text( argv[counter] );
    moduleFileList.insertEnd( text );
    counter += 1;
  }
  if( moduleFileList.size() == 0 ) {
    string text( srcdir + filenameModuleList );
    moduleFileList.insertEnd( text );
  }
  cseis_build( srcdir, bindir, moduleFileList );
  fprintf(stderr,"Build successful.\n");
  return 0;
}


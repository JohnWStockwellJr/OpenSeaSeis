/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include <string>
#include <cstring>
#include <stdarg.h>

#include "cseis_defines.h"
#include "geolib_string_utils.h"
#include "csUserConstant.h"
#include "csCompareVector.h"
#include "csKey.h"
#include "csVector.h"
#include "csUserParam.h"
#include "csException.h"
#include "csParamDescription.h"
#include "csParamDef.h"
#include "csLogWriter.h"
#include "csFlexNumber.h"


/// !CHANGE! these methods, maybe put in different file or class...

namespace cseis_system {


void exitOnError( char const* text, ... ) {
  va_list argList;
  va_start( argList, text );
  vfprintf( stderr, text, argList );
  fprintf( stderr, "\n" );
  exit(-1);
}

  std::string replaceUserConstants( char const* line, cseis_geolib::csVector<cseis_system::csUserConstant> const* list ) {
  int length = strlen(line);
  std::string strNew = line;
  int nFound = 0;
  int counter = 0;
  while( counter < length ) {
    if( line[counter] == csUserConstant::LETTER_DEFINE ) {
      int smallLength = strlen(line+counter);
      csUserConstant const* defPtr;
      for( int i = 0; i < list->size(); i++ ) {
        defPtr = &(list->at(i));
        if( smallLength >= defPtr->nameLength && !defPtr->name.compare(0,defPtr->nameLength,line+counter,0,defPtr->nameLength) ) {
          counter += defPtr->nameLength-1;
          strNew = cseis_geolib::replaceStr( strNew, defPtr->name, defPtr->value );
          nFound++;
          break;
        }
      }
    }
    counter++;
  }
  return strNew;
}

/**
* Create individual flow from master flow
*
* @return Number of modules
*/
int cseis_create_flow( FILE* f_master_flow, FILE* f_flow, cseis_geolib::csVector<csUserConstant> const* masterConstants ) {
  int counterModules = 0;
  int counterLines   = 0;
  char line[MAX_LINE_LENGTH];
  cseis_geolib::csVector<std::string> tokenList;

  while( fgets( line, MAX_LINE_LENGTH, f_master_flow ) != NULL ) {
    counterLines++;
    char c;
    if( cseis_geolib::firstNonBlankChar( line, c ) ) {
      std::string tmpStr = replaceUserConstants( line, masterConstants );
      fprintf(f_flow,"%s",tmpStr.c_str());
    }
    else {
      fprintf(f_flow,"%s",line);
    }
  }
  return counterModules;
}
//----------------------------------------------------------------------------------
/**
* Read spreadsheet file that contains user constant table (sort of survey database)
*
*/
void cseis_read_spreadSheet( char const* filenameSheet, FILE* f_spreadSheet,
                 cseis_geolib::csCompareVector<csUserConstant>* masterConstants, 
                 cseis_geolib::csVector< cseis_geolib::csVector<std::string> >* spreadSheetConstantList ) {  
    
  cseis_geolib::csVector<std::string> tokenList;
  int counterLines = 0;
  int nConstants = 0;
  char line[MAX_LINE_LENGTH];
  // Read in first line. Expected to find list of constant names in first line
  // Constant names must start with a letter
  // First column must contain unique key
  if( fgets( line, MAX_LINE_LENGTH, f_spreadSheet ) != NULL ) {
    counterLines++;
    tokenize( line, tokenList );
    nConstants = tokenList.size();
  }
  if( nConstants == 0 ) {
    throw( cseis_geolib::csException("Error in spread sheet. Empty first line, expected constant names.") );
  }
  // Tokenize constant name list
  for( int i = 0; i < tokenList.size(); i++ ) {
    csUserConstant def( tokenList.at(i), "" );
    if( def.nameLength == 2 ) {
      exitOnError("Error in spread sheet '%s'. Empty first line, expected constant names.\n", filenameSheet );
    }
    else if( !cseis_geolib::isLetter(def.name[1]) ) {
      exitOnError("Error in spread sheet '%s', line %d. Constant names must start with a letter: '%s'\n",
        filenameSheet, counterLines, def.pureName().c_str() );
    }
    if( masterConstants->contains(def) ) {
      exitOnError("Error in spread sheet '%s', line %d. Duplicate constant name: '%s'\n",
        filenameSheet, counterLines, def.pureName().c_str() );
    }
    masterConstants->insertEnd( def );
  } // END: for loop, tokenize all constant names
  cseis_geolib::csCompareVector<cseis_geolib::csKey> keyList;
  while( fgets( line, MAX_LINE_LENGTH, f_spreadSheet ) != NULL ) {
    counterLines++;
    tokenList.clear();
    tokenize( line, tokenList );
    if( tokenList.size() != masterConstants->size() ) {
      exitOnError("Error in spread sheet '%s', line %d: %s Inconsistent number of columns. Expected: %d, found: %d.\n",
        filenameSheet, counterLines, line, masterConstants->size(), tokenList.size() );
    }
    if( keyList.contains( cseis_geolib::csKey(tokenList.at(0)) ) ) {
      exitOnError("Error in spread sheet '%s', line %d: %s Duplicate key found. First column must contain a unique key.\n",
        filenameSheet, counterLines, line );
    }
    keyList.insertEnd( cseis_geolib::csKey(tokenList.at(0)) );
    spreadSheetConstantList->insertEnd( tokenList );
  }
  if( counterLines < 2 ) {
    exitOnError("Error in spread sheet '%s'. File does not contain any values.\n", filenameSheet );
  }
/*
  for( int i = 0; i < masterConstants->size(); i++ ) {
    fprintf(stderr,"Master constant %d %s...\n", i, masterConstants->at(i).pureName().c_str() );
    for( int k = 0; k < spreadSheetConstantList->size(); k++ ) {
      fprintf(stderr,"Master constant value %d %s...\n", k, spreadSheetConstantList->at(k).at(i).c_str() );
    }
  }
*/
}
//---------------------------------------------------------------
/**
*
*
*
*/
void cseis_read_globalConst( FILE* f_globalConst, cseis_geolib::csCompareVector<csUserConstant>* globalConstList ) {
  cseis_geolib::csVector<std::string> tokenList(3);
  int counterLines = 0;
  char line[MAX_LINE_LENGTH];
  while( fgets( line, MAX_LINE_LENGTH, f_globalConst ) != NULL ) {
    counterLines++;
    char c;
    if( cseis_geolib::firstNonBlankChar( line, c ) && c != LETTER_COMMENT ) {
      if( c == csUserConstant::LETTER_DEFINE ) {
        tokenList.clear();
        tokenize( line, tokenList );
        if( !tokenList.at(0).compare(csUserConstant::defineWord()) ) {
          if( tokenList.size() != 3 ) {
            exitOnError("Syntax error in global constant input file, line %d: %s &define statement requires two arguments: &define <name> <value>.",
                  counterLines, line );
          }
          csUserConstant uc(tokenList.at(1), tokenList.at(2));
          if( globalConstList->contains( uc ) ) {
            exitOnError("Syntax error in global constant input file, line %d: %s Duplicate definition of global constant '%s'.",
                  counterLines, line, uc.pureName().c_str() );
          }
          globalConstList->insertEnd( uc );
        }
      }
    }
  }
  
}

bool checkParameters( char const* moduleName, csParamDef const* paramDef, cseis_geolib::csVector<csUserParam*>* userParams, csLogWriter* log ) {
  bool returnFlag = true;
  if( paramDef->module() == NULL ) {
    log->line("Program bug: Module name not defined in parameter definition. Should be '%s'", moduleName );
    return false;
  }
  if( strcmp( moduleName, paramDef->module()->name() ) ) {
    log->line("Program bug: Inconsistent module names between parameter definition: '%s', and method definition: '%s'\n", paramDef->module()->name(), moduleName );
    returnFlag = false;
  }

  cseis_geolib::csVector<csParamDescription const*> paramDefList(5);  // User parameter definitions for current module, e.g. name & description of parameter, such as 'filename' 'Input file name'  (set in module 'params' method)
  cseis_geolib::csVector<csParamDescription const*> valueDefList(2);  // Definitions of all values for current user parameter (set in module 'params' method)
  cseis_geolib::csVector<csParamDescription const*> optionList(3);
  cseis_geolib::csVector<std::string> userValueList(2);         // All user specified values for current parameter

  paramDef->getParameters( &paramDefList );
  int nDefinedParams = paramDefList.size();

  int nUserParams = userParams->size();

  // Go through all parameters specified by user
  for( int iUserParam = 0; iUserParam < nUserParams; iUserParam++ ) {
    csUserParam* userParam = userParams->at(iUserParam);
    //    char const* userParamName = userParam->getName().c_str();
    int ip = -1;
    // Go through all parameters defined for this module, try to find matching parameter name
    for( int i = 0; i < nDefinedParams; i++ ) {
      //      if( !strcmp( userParamName, paramDefList.at(i)->name() ) ) {
      if( userParam->getName().compare( paramDefList.at(i)->name() ) == 0 ) {
        ip = i;
        break;
      }
    }
    if( ip < 0 ) {
      //      if( !strcmp(userParamName,"debug") || !strcmp(userParamName,"version") ) {
      if( userParam->getName().compare("debug") == 0 || userParam->getName().compare("version") == 0 ) {
        continue;
      }
      log->line("Error: Unknown user specified parameter: '%s'", userParam->getName().c_str() );
      returnFlag = false;
      continue;  // Check other parameters
    }
    // fprintf(stdout,"\nUser param: '%s'  '%s'\n", userParamName, paramDefList.at(ip)->desc() );

    int nUserValues = userParam->getNumValues();  // Number of parameter values SPECIFIED BY USER
    valueDefList.clear();
    paramDef->getValues( ip, &valueDefList );
    int nDefinedValues = valueDefList.size();   // Number of values DEFINED for this parameter

    if( paramDefList.at(ip)->type() == NUM_VALUES_FIXED && nUserValues < nDefinedValues ) {
      log->line("Error: Too few user specified values for parameter '%s'. Required: %d, found: %d", userParam->getName().c_str(), nDefinedValues, nUserValues );
      returnFlag = false;
      continue;
    }
    else if( nUserValues == 0 ) {
      log->line("Error: No value specified for parameter '%s'.", userParam->getName().c_str() );
      returnFlag = false;
      continue;
    }
    // Convert all user specified parameter values to lower case (except for TYPE_STRING values)
    // Also, check all number values
    int minNumValues = std::min( nDefinedValues, nUserValues );
    cseis_geolib::csFlexNumber valueTmp;
    for( int i = 0; i < minNumValues; i++ ) {
      int valueType = valueDefList.at(i)->type();
      if( valueType == VALTYPE_NUMBER ) {
        std::string text = userParam->getValue( i );
        if( !valueTmp.convertToNumber( userParam->getValue( i ) ) ) {
          log->line("Error: User parameter '%s': Value is not recognised as valid number: '%s'", userParam->getName().c_str(), text.c_str() );
          returnFlag = false;
        }
      }
      else if( valueType == VALTYPE_OPTION ) {
        //        userParam->setValueType( i, valueType );
        userParam->convertToLowerCase( i );
      }
    }
    if( nUserValues > nDefinedValues ) {
      for( int i = nDefinedValues; i < nUserValues; i++ ) {
        int valueType = valueDefList.at(nDefinedValues-1)->type();
        if( valueType == VALTYPE_OPTION ) {
          userParam->convertToLowerCase( i );
        }
      }
    }
    // For parameters taking a variable number of values, convert the last values (all same type)
    //    if( paramDefList.at(ip)->type() == NUM_VALUES_VARIABLE && valueDefList.at(nDefinedValues-1)->type() == VALTYPE_OPTION ) {
    //      userParam->setValueTypeVariable( nDefinedValues, valueDefList.at(nDefinedValues-1)->type() );
    //    }
    userValueList.clear();
    userParam->getValues( &userValueList );

    // Go through all parameter values, check correctness of OPTION values
    for( int iv = 0; iv < minNumValues; iv++ ) {
      if( valueDefList.at(iv)->type() == VALTYPE_OPTION ) {
        //        char const* userOptionName = userValueList.at(iv).c_str();
        // Check options:
        optionList.clear();
        paramDef->getOptions( ip, iv, &optionList );
        int optionFound = false;
        int nOptions = optionList.size();
        for( int io = 0; io < nOptions; io++ ) {
          if( userValueList.at(iv).compare( optionList.at(io)->name() ) == 0 ) {
            optionFound = true;
            break;
          }
        }
        if( !optionFound ) {
          log->line("Error: Unknown user specified option for parameter '%s':  %s", userParam->getName().c_str(), userValueList.at(iv).c_str() );
          log->write("Valid options are: ");
          for( int io = 0; io < nOptions; io++ ) {
            log->write("'%s' (%s)", optionList.at(io)->name(), optionList.at(io)->desc());
            if( io < nOptions-1 ) log->write(" / ");
          }
          log->line("");
          returnFlag = false;
        }
      } // END if
      else {
        if( valueDefList.at(iv)->type() != VALTYPE_STRING ) {
        }
      }
      //      fprintf(stdout," Name: %s, valueType: %d\n", valueDefList.at(iv)->name(), valueDefList.at(iv)->type());
    }
    
  }
  
  //int nParam = ;
  //  for( int 


  return returnFlag;
}

} // namespace


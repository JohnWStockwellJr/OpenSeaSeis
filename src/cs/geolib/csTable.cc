/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csTable.h"
#include "csException.h"
#include "csTimeFunction.h"
#include "csVector.h"
#include "csGeolibUtils.h"
#include "geolib_string_utils.h"
#include <string>
#include <cstring>

using namespace cseis_geolib;

csTable::csTable() {
  myValues1D   = NULL;
  myKeyValues  = NULL;
  myKeyNames   = NULL;
  myValueNames = NULL;
  myTimeFunctions2D = NULL;
  myColumnTypes = NULL;
  myFile       = NULL;
  myNumKeys    = 0;
  myNumValues  = 0;
  myTableDimension = TABLE_DIM_UNKNOWN;

  myNumCols        = 0;
  myNumLocations   = 0;
  myIndexTimeCol   = 0;
  myIndexFirstValueCol = 0;
  myHasReadTableContents = false;
  myHasBeenInitialized   = false;

  myCurrentTimeFunction = new csTimeFunction<double>();
}
csTable::~csTable() {
  if( myFile != NULL ) {
    fclose( myFile );
    myFile = NULL;
  }
  if( myValueNames ) {
    delete [] myValueNames;
    myValueNames = NULL;
  }
  if( myKeyNames ) {
    delete [] myKeyNames;
    myKeyNames = NULL;
  }
  if( myTimeFunctions2D != NULL ) {
    for( int i = 0; i < myNumLocations; i++ ) {
      if( myTimeFunctions2D[i] != NULL ) {
        delete myTimeFunctions2D[i];
      }
    }
    delete [] myTimeFunctions2D;
    myTimeFunctions2D = NULL;
  }
  if( myCurrentTimeFunction ) {
    delete myCurrentTimeFunction;
    myCurrentTimeFunction = NULL;
  }
  if( myColumnTypes != NULL ) {
    delete [] myColumnTypes;
    myColumnTypes = NULL;
  }
  clearBuffers();
}
void csTable::clearBuffers() {
  if( myValues1D != NULL ) {
    for( int ival = 0; ival < myNumValues; ival++ ) {
      delete myValues1D[ival];
    }
    delete [] myValues1D;
    myValues1D = NULL;
  }
  if( myKeyValues != NULL ) {
    for( int i = 0; i < myNumLocations; i++ ) {
      delete [] myKeyValues[i];
    }
    delete [] myKeyValues;
    myKeyValues = NULL;
  }
}
//-----------------------------------------------------------------------------------
//
//
void csTable::initialize( std::string const& tablename, std::string const& filename, int methodInterpolation ) {
  if( myHasBeenInitialized ) {
    throw( csException("Table has been initialized before. Previous input file name: '%s', new file name: %s", myFilename.c_str(), filename.c_str()) );
  }
  myTablename = tablename;
  myFilename  = filename;
  myMethodInterpolation = methodInterpolation;

  clearBuffers();
  char buffer[1024];

  if( (myFile = fopen( filename.c_str(), "r" )) == (FILE*) NULL ) {
    throw csException("Could not open file: '%s'", filename.c_str());
  }
  
  if( fgets( buffer, 1024, myFile ) == NULL ) {
    throw( csException("Table input file '%s' contains no lines", myFilename.c_str()) );
  }

  // Convert names of keys/time/values to lower case
  toLowerCase( buffer );

  // 1) Set key names, and determine number of keys, values and time axis
  csVector<std::string> tokenList;
  tokenize( buffer, tokenList, false );
  myNumCols = tokenList.size();
  if( myNumCols == 0 ) {
    throw( csException("Table input file '%s' contains no lines", myFilename.c_str()) );
  }
  else if( myNumCols < 2 ) {
    throw( csException("Table input file '%s' contains no key or value column", myFilename.c_str()) );
  }

  myNumKeys = 0;
  while( myNumKeys < myNumCols && tokenList.at(myNumKeys)[0] == KEY_CHAR ) {
    myNumKeys += 1;
  }
  //if( myNumKeys == 0 ) {
  //  throw( csException("Table input file '%s' contains no key column. Key column names are preceded by the character '%c', e.g. %csource",
  //                     myFilename.c_str(), KEY_CHAR, KEY_CHAR) );    
  // }
  if( myNumKeys == myNumCols ) {
    throw( csException("Table input file '%s' contains no value column. Names of value columns must NOT be preceded by the key character '%c'",
                       myFilename.c_str(), KEY_CHAR) );    
  }

  std::string token = tokenList.at(myNumKeys);
  if( !token.compare("time") ) {
    myTableDimension  = TABLE_DIM_2D;
    myNumValues       = myNumCols-myNumKeys-1;
    myIndexTimeCol    = myNumKeys;
    if( myNumValues > 1 ) {
      throw( csException("Table input file '%s' contains %d value columns. Only one value column is currently supported for 2D tables",
                         myFilename.c_str(), myNumValues ) );    
    }
  }
  else {
    myTableDimension = TABLE_DIM_1D;
    myNumValues      = myNumCols-myNumKeys;
  }
  //timeFunctionList
  if( myNumValues == myNumCols ) {
    throw( csException("Table input file '%s' contains no value column.", myFilename.c_str()) );    
  }

  myValueNames = new std::string[myNumValues];

  if( myNumKeys > 0 ) {
    myKeyNames = new std::string[myNumKeys];
    for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
      int length = tokenList.at(ikey).length();
      myKeyNames[ikey] = tokenList.at(ikey).substr(1,length-1);
      printf("Key #%d: '%s'\n", ikey+1, myKeyNames[ikey].c_str() );
    }
  }
  myIndexFirstValueCol = myNumCols - myNumValues;
  myColumnTypes = new type_t[myNumCols];
  for( int icol = 0; icol < myNumCols; icol++ ) {
    myColumnTypes[icol] = TYPE_UNKNOWN;
  }  
  myHasBeenInitialized = true;
}
//------------------------------------------------------------------------
//
//
void csTable::readTableContents( type_t const* types, int numCols ) {
  if( !myHasBeenInitialized ) {
    throw( csException("csTable::readTableContents(): Need to call method 'initialize' first before calling this method. This is a program bug in the calling function") );
  }
  else if( numCols != myNumCols ) {
    throw( csException("csTable::setColumnTypes: Wrong number of columns supplied: %d, actual number of columns: %d", numCols, myNumCols) );
  }
  if( myHasReadTableContents ) return;
  myHasReadTableContents = true;

  for( int icol = 0; icol < myNumCols; icol++ ) {
    if( types[icol] == TYPE_DOUBLE || types[icol] == TYPE_FLOAT || types[icol] == TYPE_INT ) {
      myColumnTypes[icol] = types[icol];
    }
    else throw( csException("Wrong column type specified. Table columns must have number types. Specified column type: '%s'", csGeolibUtils::typeText(types[icol])) );
  }

  // 2) Read in all values, key values, time values...
  csVector<double> timeList;
  csVector<double> valueListTime;

  csVector<double>* valueList = NULL;
  if( myTableDimension == TABLE_DIM_1D ) {
    valueList = new csVector<double>[myNumValues];
  }

  csVector<csTimeFunction<double>*> timeFunctionList;
  csVector<double*> keyValueList;
  double* keysCurrent = NULL;
  double* keysNew     = NULL;
  if( myNumKeys > 0 ) {
    keysCurrent = new double[myNumKeys];
    keysNew     = new double[myNumKeys];
    for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
      keysCurrent[ikey] = -99999.2;
    }
  }

  char buffer[1024];
  int counterLines = 1;
  csVector<std::string> tokenList;

  //-------------------------------
  // (1) Loop through all input lines
  //
  while( fgets( buffer, 1024, myFile ) != NULL ) {
    tokenList.clear();
    tokenize( buffer, tokenList );
    if( tokenList.size() == 0 ) continue;  // Assume blank line, do nothing

    // Extract key values
    if( myNumKeys > 0 ) {
      bool isSame = true;
      for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
        keysNew[ikey] = atof(tokenList.at(ikey).c_str());
        if( keysNew[ikey] != keysCurrent[ikey] ) isSame = false;
//        printf(" Key value #%d (%s):  %f\n", ikey+1, myKeyNames[ikey].c_str(), keysNew[ikey] );
      }
      if( !isSame ) {  // Key value in current line differs from key value in previous line
        if( myTableDimension == TABLE_DIM_1D ) {
          double* keysTMP = new double[myNumKeys];
          memcpy( keysTMP, keysNew, myNumKeys*sizeof(double) );
          keyValueList.insertEnd( keysTMP );
          memcpy( keysCurrent, keysNew, myNumKeys*sizeof(double) );
          //          keysNew = new double[myNumKeys];
        }
        else if( myTableDimension == TABLE_DIM_2D && timeList.size() != 0 ) {
          csTimeFunction<double>* timeFunc = new csTimeFunction<double>();
          timeFunc->set( &valueListTime, &timeList );
          timeFunctionList.insertEnd( timeFunc );
          valueListTime.clear();
          timeList.clear();

          double* keysTMP = new double[myNumKeys];
          memcpy( keysTMP, keysCurrent, myNumKeys*sizeof(double) );
          keyValueList.insertEnd( keysTMP );
          //          keysCurrent = new double[myNumKeys];
          memcpy( keysCurrent, keysNew, myNumKeys*sizeof(double) );
        }
        else {
          memcpy( keysCurrent, keysNew, myNumKeys*sizeof(double) );
        }
      }
      else if( myTableDimension == TABLE_DIM_1D ) {
        throw( csException("Input table (1D) file '%s' has two lines with same key values, at line #%d.",
                           myFilename.c_str(), counterLines ) );
      }
    }

    if( myTableDimension == TABLE_DIM_1D ) {
      for( int ival = 0; ival < myNumValues; ival++ ) {
        //        csFlexNumber* value = new csFlexNumber( myColumnTypes[ival+myIndexFirstValueCol],  );
        valueList[ival].insertEnd( atof(tokenList.at(myIndexFirstValueCol+ival).c_str()) );
      }
    }
    else { // if( myTableDimension == TABLE_DIM_2D ) {
      timeList.insertEnd( atof(tokenList.at(myIndexTimeCol).c_str()) );
      valueListTime.insertEnd( atof(tokenList.at(myIndexFirstValueCol).c_str()) );
    }
    counterLines++;
  }  // end while reading lines from input file
  fclose( myFile );
  myFile = NULL;

  if( myTableDimension == TABLE_DIM_2D && timeList.size() != 0 ) {
    double* keysTMP = new double[myNumKeys];
    memcpy( keysTMP, keysCurrent, myNumKeys*sizeof(double) );
    keyValueList.insertEnd( keysTMP );
    csTimeFunction<double>* timeFunc = new csTimeFunction<double>();
    timeFunc->set( &valueListTime, &timeList );
    timeFunctionList.insertEnd( timeFunc );
  }

  myNumLocations = keyValueList.size();
  if( myNumLocations > 0 ) {
    myKeyValues = new double const*[myNumLocations];
    for( int i = 0; i < myNumLocations; i++ ) {
      myKeyValues[i] = keyValueList.at(i);
    }
  }
  else { // No keys specified --> there is only a single location
    myNumLocations = 1;
  }
  if( myTableDimension == TABLE_DIM_1D ) {
    myValues1D = new double*[myNumValues];
    for( int ival = 0; ival < myNumValues; ival++ ) {
      myValues1D[ival] = new double[myNumLocations];
      for( int i = 0; i < myNumLocations; i++ ) {
        myValues1D[ival][i] = valueList[ival].at(i);
        //        delete valueList[ival].at(i);  // Free csFlexNumber pointer
      }
    }
  }
  else if( myTableDimension == TABLE_DIM_2D ) {
    myTimeFunctions2D = new csTimeFunction<double>* [myNumLocations];
    for( int i = 0; i < myNumLocations; i++ ) {
      myTimeFunctions2D[i] = timeFunctionList.at(i);
    }
  }

  if( valueList != NULL ) delete [] valueList;
  if( keysCurrent ) {
    delete [] keysCurrent;
  }
  if( keysNew ) {
    delete [] keysNew;
  }

}

//-----------------------------------------------------------------------------------
//
void csTable::initialize_simpleTimeTable( int numLocations ) {
  myNumKeys   = 0;
  myNumValues = 1;
  myNumLocations = numLocations;
  myKeyValues = new double const*[myNumLocations];

  myTableDimension = TABLE_DIM_2D;
  myTimeFunctions2D = new csTimeFunction<double>*[myNumLocations];
  for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
    myTimeFunctions2D[iloc] = new csTimeFunction<double>();
  }
}
void csTable::setValue_simpleTimeTable( int locationIndex, csTimeFunction<double> const* timeFunction ) {
  myTimeFunctions2D[locationIndex]->set( timeFunction );
}
/*void setValues( int locationIndex, double const* keyValues, double const* values ) {
  myKeyValues[locationIndex] = new double[myNumKeys];
  for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
    myKeyValues[locationIndex][ikey] = keyValues[ikey];
  }
  if( myTableDimension == TABLE_DIM_1D ) {
    myValues1D = new double[myNumLocations];
    for( int i = 0; i < myNumLocations; i++ ) {
      myValues1D[i] = valueList.at(i);
    }
  }
  else if( myTableDimension == TABLE_DIM_2D ) {
    myTimeFunctions2D = new csTimeFunction<double> const* [myNumLocations];
    for( int i = 0; i < myNumLocations; i++ ) {
      myTimeFunctions2D[i] = timeFunctionList.at(i);
    }
  }
}
*/
//-----------------------------------------------------------------------------------
//

double const* csTable::getKeyValues( int indexLocation ) const {
  if( indexLocation >= 0 && indexLocation < myNumLocations && myNumKeys > 0 ) {
    return myKeyValues[indexLocation];
  }
  else {
    return NULL;
  }
}

//-----------------------------------------------------------------------------------
//
csTimeFunction<double> const* csTable::getFunction( double const* keyValues_in ) const {
  int locLeft;
  int locRight;
  double weightLoc;

  findKeyLocation( keyValues_in, locLeft, locRight, weightLoc );

  csTimeFunction<double> const* timeFuncLeft  = myTimeFunctions2D[locLeft];
  csTimeFunction<double> const* timeFuncRight = myTimeFunctions2D[locRight];

//  fprintf(stdout,"Time set location left %d,  right %d, weight: %f\n", locLeft, locRight, weightLoc );

  csVector<double> timeList;
  int timeIndexLeft  = 0;
  int timeIndexRight = 0;
  int numTimesLeft   = timeFuncLeft->numValues();
  int numTimesRight  = timeFuncRight->numValues();


  while( timeIndexLeft < numTimesLeft || timeIndexRight < numTimesRight ) {
    if( timeIndexRight == numTimesRight ) {
      timeList.insertEnd(timeFuncLeft->timeAtIndex(timeIndexLeft));
      timeIndexLeft += 1;
    }
    else if( timeIndexLeft == numTimesLeft ) {
      timeList.insertEnd(timeFuncRight->timeAtIndex(timeIndexRight));
      timeIndexRight += 1;
    }
    else {
      double timeLeft  = timeFuncLeft->timeAtIndex(timeIndexLeft);
      double timeRight = timeFuncRight->timeAtIndex(timeIndexRight);
      if( timeLeft == timeRight ) {
        timeList.insertEnd(timeLeft);
        timeIndexLeft  += 1;
        timeIndexRight += 1;
      }
      else if( timeLeft < timeRight ) {
        timeList.insertEnd(timeLeft);
        timeIndexLeft  += 1;
      }
      else {
        timeList.insertEnd(timeRight);
        timeIndexRight += 1;
      }
    }
  } // END while loop through all times

  int numTimes = timeList.size();
  csVector<double> valueList;
  for( int itime = 0; itime < numTimes; itime++ ) {
    double time = timeList.at(itime);
    double valueLeft = timeFuncLeft->valueAt( time );
    double valueRight= timeFuncRight->valueAt( time );
//    fprintf(stdout,"Time #%d: %f, values: %f %f\n", itime, time, valueLeft, valueRight);
    valueList.insertEnd( valueLeft + weightLoc * ( valueRight - valueLeft ) );
  }

  myCurrentTimeFunction->set( &valueList, &timeList );
  return myCurrentTimeFunction;
}
//-----------------------------------------------------------------------------------
//
csTimeFunction<double> const* csTable::getFunction( int indexLocation ) const {
  return NULL;
}
//-----------------------------------------------------------------------------------
//
std::string const csTable::keyName( int indexKey ) const {
  if( indexKey >= 0 && indexKey < myNumKeys ) {
    return myKeyNames[indexKey];
  }
  else {
    throw( csException("csTable::keyName(): Incorrect key index passed") );
  }
}
//-----------------------------------------------------------------------------------
//
std::string const csTable::valueName( int indexValue ) const {
  if( indexValue >= 0 && indexValue < myNumValues ) {
    return myValueNames[indexValue];
  }
  else {
    throw( csException("csTable::valueName(): Incorrect value index passed") );
  }
}

//-----------------------------------------------------------------------------------
//
double csTable::getKeyValue( int indexLocation, int indexKey ) const {
  if( indexLocation >= 0 && indexLocation < myNumLocations && indexKey >= 0 && indexKey < myNumKeys ) {
    return myKeyValues[indexLocation][indexKey];
  }
  else {
    throw( csException("csTable::getKeyValue(): Incorrect location or key index passed") );
  }
}
//-----------------------------------------------------------------------------------
//
double csTable::getValue( int indexLocation, int indexValue ) const {
  if( indexValue >= 0 && indexValue < myNumValues && indexLocation >= 0 && indexLocation < myNumLocations && myTableDimension == TABLE_DIM_1D ) {
    return myValues1D[indexValue][indexLocation];
  }
  else {
    throw( csException("csTable::getValue(): Incorrect value or location index passed, or incorrect use of getValue for non-1D table") );
  }
}
//-----------------------------------------------------------------------------------
//
double csTable::getValue( int indexLocation, double time ) const {
  if( indexLocation >= 0 && indexLocation < myNumLocations && myTableDimension == TABLE_DIM_2D ) {
    return 0.0;
  }
  else {
    throw( csException("csTable::getValue(): Incorrect location index passed, or incorrect use of getValue for non-2D table") );
  }
}
//-----------------------------------------------------------------------------------
//
double csTable::getValue( double const* keyValues, int indexValue ) const {
  if( myNumKeys > 0 || myTableDimension == TABLE_DIM_1D ) {
    return interpolate1D( indexValue, keyValues );
  }
  else {
    throw( csException("csTable::getValue(): Table does not contain any keys, or incorrect use for non-1D table") );
  }
}
//-----------------------------------------------------------------------------------
//
double csTable::getValue( double const* keyValues, double time ) const {
  if( myNumKeys > 0 || myTableDimension == TABLE_DIM_2D ) {
    return interpolate2D( keyValues, time );
  }
  else {
    throw( csException("csTable::getValue(): Table does not contain any keys, or incorrect use for non-2D table") );
  }
}

//**************************************************************************************
//
// Internal interpolation methods
//
//**************************************************************************************

double csTable::interpolate2D( double const* keyValues_in, double time ) const {
  int locLeft;
  int locRight;
  double weightLoc;

  findKeyLocation( keyValues_in, locLeft, locRight, weightLoc );

  double valueLeft  = myTimeFunctions2D[locLeft]->valueAt( time );
  double valueRight = myTimeFunctions2D[locRight]->valueAt( time );

  return( valueLeft + weightLoc * ( valueRight - valueLeft ) );
}

double csTable::interpolate1D( int valueIndex, double const* keyValues_in ) const {
  int locLeft;
  int locRight;
  double weight;

  findKeyLocation( keyValues_in, locLeft, locRight, weight );

  double valLeft  = myValues1D[valueIndex][locLeft];
  double valRight = myValues1D[valueIndex][locRight];

  return( valLeft + weight * ( valRight - valLeft ) );
}

// Assumes that keys are all sorted with increasing order
// Not implemented fully! Deals with single key only
void csTable::findKeyLocation( double const* keyValues_in, int& locLeft, int& locRight, double& weight ) const {
  if( myNumKeys > 1 ) {
    throw( csException("csTable::findKeyLocation(): Table contains more than 1 key. This is not supported yet.") );
  }
  else if( myNumKeys == 0 ) {
    locLeft  = 0;
    locRight = 0;
    weight   = 1.0;
    return;
  }
//  int location = myNumLocations-1;
//  for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
  int ikey = 0;
  double currentKeyValue = keyValues_in[ikey];
  for( int location = myNumLocations-1; location >= 0; location-- ) {
    if( currentKeyValue >= myKeyValues[location][ikey] ) {
      if( location == myNumLocations-1 ) {
        locLeft  = location;
        locRight = location;
        weight   = 1.0;
      }
      else if( currentKeyValue != myKeyValues[location][ikey] ) {
        locLeft  = location;
        locRight = location+1;
        // Interpolate linearly
        weight = (currentKeyValue-myKeyValues[location][ikey])/(myKeyValues[location+1][ikey]-myKeyValues[location][ikey]);
      }
      else { // currentKeyValue == myKeyValues[location][ikey]
        locLeft  = location;
        locRight = location;
        weight = 1.0;
      }
      return;
    }
  }
  locLeft  = 0;
  locRight = 0;
  weight   = 1.0;
  return;
}


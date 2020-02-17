/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csTableNew.h"
#include "csTimeFunction.h"
#include "csException.h"
#include "csVector.h"
#include "csGeolibUtils.h"
#include "geolib_string_utils.h"
#include "csSort.h"
#include "csSortManager.h"
#include "csFlexNumber.h"
#include <string>
#include <cstring>
#include <algorithm>

using namespace cseis_geolib;

csTableNew::csTableNew( int tableType, int timeColumnIndex ) {
  if( tableType != TABLE_TYPE_TIME_FUNCTION ) {
    throw(csException("csTable::csTable: Program bug in calling function: For tables other than 'time function tables', use different constructor"));
  }
  init(tableType);
  myIndexTimeCol = timeColumnIndex;
}
csTableNew::csTableNew( int tableType ) {
  if( tableType == TABLE_TYPE_TIME_FUNCTION ) {
    throw(csException("csTable::csTable: Program bug in calling function: For 'time function tables', use different constructor"));
  }
  init(tableType);
}
void csTableNew::init( int tableType ) {
  myValues     = NULL;
  myKeyValues  = NULL;
  myFile       = NULL;
  myNumKeys    = 0;
  myNumInterpKeys    = 0;
  myNumValues  = 0;
  myTableType  = tableType;
  myKeyCols = NULL;
  myKeyAllCols = NULL;
  myKeyInterpCols = NULL;
  //  myKeyInterpolate = NULL;
  myValueColumns   = NULL;

  myNumLocations   = 0;
  myHasReadTableContents = false;
  myHasBeenInitialized   = false;

  myTimeFunctions2D = NULL;
  myIndexTimeCol    = 0;
  myCurrentTimeFunction = NULL;
  if( myTableType == TABLE_TYPE_TIME_FUNCTION ) {
    myCurrentTimeFunction = new csTimeFunction<double>();
  }
  else {
    myValues = new csVector<csTableValueList*>();
  }
}
csTableNew::~csTableNew() {
  clearBuffers();
  if( myKeyCols != NULL ) {
    delete [] myKeyCols;
    myKeyCols = NULL;
  }
  if( myKeyAllCols != NULL ) {
    delete [] myKeyAllCols;
    myKeyAllCols = NULL;
  }
  if( myKeyInterpCols != NULL ) {
    delete [] myKeyInterpCols;
    myKeyInterpCols = NULL;
  }
  if( myValueColumns != NULL ) {
    delete [] myValueColumns;
    myValueColumns = NULL;
  }
  if( myValues != NULL ) {
    for( int ie = 0; ie < myValues->size(); ie++ ) {
      delete myValues->at(ie);
    }
    delete myValues;
    myValues = NULL;
  }
  if( myFile != NULL ) {
    fclose( myFile );
    myFile = NULL;
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
}
void csTableNew::addKey( int columnIndex, bool doInterpolate ) {
  if( doInterpolate ) {
    if( myNumInterpKeys == 2 ) {
      throw(csException("csTableNew::addKey(): Maximum 2 interpolation keys supported."));
    }
    int newNumKeys = myNumInterpKeys + 1;
    int* keyCols = new int[newNumKeys];
    if( myNumInterpKeys > 0 ) {
      for( int ikey = 0; ikey < myNumInterpKeys; ikey++ ) {
        keyCols[ikey]     = myKeyInterpCols[ikey];
      }
      delete [] myKeyInterpCols;
    }
    myKeyInterpCols = keyCols;
    myNumInterpKeys = newNumKeys;
    myKeyInterpCols[myNumInterpKeys-1] = columnIndex;
    if( myNumInterpKeys > 1 ) {
      if( myKeyInterpCols[myNumInterpKeys-1] <= myKeyInterpCols[myNumInterpKeys-2] ) {
        throw( csException("csTableNew::addKey: Key columns must be specified in increasing order") );
      }
    }
  }
  else {
    int newNumKeys = myNumKeys + 1;
    int* keyCols = new int[newNumKeys];
    if( myNumKeys > 0 ) {
      for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
        keyCols[ikey]     = myKeyCols[ikey];
      }
      delete [] myKeyCols;
    }
    myKeyCols = keyCols;
    myNumKeys = newNumKeys;
    myKeyCols[myNumKeys-1] = columnIndex;
    if( myNumKeys > 1 ) {
      if( myKeyCols[myNumKeys-1] <= myKeyCols[myNumKeys-2] ) {
        throw( csException("csTableNew::addKey: Key columns must be specified in increasing order") );
      }
    }
  }
  myNumAllKeys = myNumKeys + myNumInterpKeys;
}


void csTableNew::addValue( int columnIndex ) {
  int newNumValues = myNumValues + 1;
  int* valueColumns = new int[newNumValues];
  if( myNumValues > 0 ) {
    for( int i = 0; i < myNumValues; i++ ) {
      valueColumns[i]     = myValueColumns[i];
    }
    delete [] myValueColumns;
  }
  myValueColumns   = valueColumns;
  myNumValues      = newNumValues;
  myValueColumns[myNumValues-1] = columnIndex;
  if( myTableType == TABLE_TYPE_TIME_FUNCTION ) {
    //    if( myNumValues > 1 ) {
    //   throw( csException("csTableNew::addValue: Only one value column is currently supported for time function tables") );
    //  }
  }
}


void csTableNew::clearBuffers() {
  if( myKeyValues != NULL ) {
    for( int i = 0; i < myNumAllKeys; i++ ) {
      delete [] myKeyValues[i];
    }
    delete [] myKeyValues;
    myKeyValues = NULL;
  }
}

//-----------------------------------------------------------------------------------
//
//
void csTableNew::initialize( std::string const& filename, bool doSort ) {
  if( myHasBeenInitialized ) {
    throw( csException("csTableNew::initialize: Table has been initialized before. Previous input file name: '%s', new file name: %s", myFilename.c_str(), filename.c_str()) );
  }
  if( myNumValues == 0 ) {
    throw( csException("csTableNew::initialize: No value column specified. Use addValue() to add value column to table.") );
  }
  myFilename  = filename;

  clearBuffers();

  if( (myFile = fopen( filename.c_str(), "r" )) == (FILE*) NULL ) {
    throw csException("Could not open file: '%s'", filename.c_str());
  }
  myHasBeenInitialized = true;

  myNumAllKeys = myNumKeys + myNumInterpKeys;
  if( myNumAllKeys > 0 ) {
    // Check that all interpolated keys are located at higher column indexes than non-interpolated ones
    //    if( myNumKeys > 0 && myNumInterpKeys > 0 ) {
    //  if( myKeyCols[myNumKeys-1] >= myKeyInterpCols[myNumInterpKeys-1] ) {
    //    throw( csException("Interpolated keys must be located 'behind' non-interpolated keys. In other words, keys that shall be interpolated must be located (in the input ASCII file) at higher column locations (to the right) of non-interpolated keys.") );
    //  }
    // }
    myKeyAllCols = new int[myNumAllKeys];
    for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
      myKeyAllCols[ikey] = myKeyCols[ikey];
    }
    for( int ikey = 0; ikey < myNumInterpKeys; ikey++ ) {
      myKeyAllCols[ikey+myNumKeys] = myKeyInterpCols[ikey];
    }
  }
  readTableContents( doSort );
}



//------------------------------------------------------------------------
//
//
void csTableNew::readTableContents( bool doSort ) {

  // 2) Read in all values, key values, time values...
  double* keysCurrent = NULL;
  double* keysNew     = NULL;
  int maxColumnIndex = 0;
  if( myNumAllKeys > 0 ) {
    keysCurrent = new double[myNumAllKeys];
    keysNew     = new double[myNumAllKeys];
    for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
      keysCurrent[ikey] = -99999.2;
    }
  }
  if( myNumKeys > 0 ) {
    maxColumnIndex = myKeyCols[myNumKeys-1];
  }
  if( myNumInterpKeys > 0 ) {
    maxColumnIndex = std::max( maxColumnIndex, myKeyInterpCols[myNumInterpKeys-1] );
  }
  maxColumnIndex = std::max( maxColumnIndex, myValueColumns[myNumValues-1] );

  char buffer[1024];
  int counterLines = 1;
  // Helper field for storing tokens extracted from input line
  csVector<std::string> tokenList;
  // List of values for each 'value' column in input file
  csVector<double>* valueList = new csVector<double>[myNumValues];
  // List of key values for each 'location'. Each list item is an array of key values (one value for each key)
  csVector<double*> keyValueList;

  csVector<double> timeList;
  //  csVector<double> valueListTime;
  csVector<csTimeFunction<double>*> timeFunctionList;

  //-------------------------------
  // (1) Loop through all input lines
  //
  int maxTableColumns = 0;
  while( fgets( buffer, 1024, myFile ) != NULL ) {
    tokenList.clear();
    tokenize( buffer, tokenList );
    int numColumns = tokenList.size();
    if( numColumns > maxTableColumns ) maxTableColumns = numColumns;
    if( numColumns == 0 ) continue;  // Assume blank line, do nothing
    else if( tokenList.at(0)[0] == '#') continue; // Comment line, do nothing
    else if( numColumns <= maxColumnIndex ) continue;

    // Extract key values
    if( myNumAllKeys > 0 ) {
      bool isSame = true;
      for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
        keysNew[ikey] = atof(tokenList.at( myKeyAllCols[ikey] ).c_str());
        if( keysNew[ikey] != keysCurrent[ikey] ) isSame = false;
      }
      if( !isSame ) {  // Key value in current line differs from key value in previous line --> Store values up to this point under previous key
        if( myTableType != TABLE_TYPE_TIME_FUNCTION ) {
          double* keysTMP = new double[myNumAllKeys];
          memcpy( keysTMP, keysNew, myNumAllKeys*sizeof(double) );
          keyValueList.insertEnd( keysTMP );
          memcpy( keysCurrent, keysNew, myNumAllKeys*sizeof(double) );
          if( myTableType == TABLE_TYPE_DUPLICATE_KEYS && valueList[0].size() != 0 ) addData_internal( valueList );
        }
        else if ( timeList.size() != 0 ) {
          csTimeFunction<double>* timeFunc = new csTimeFunction<double>( myNumValues );  // myNumValues: spatial values in input table
          for( int ival = 0; ival < myNumValues; ival++ ) {
            //  csVector<double>* valueList = new csVector<double>[myNumValues];
            //  csVector<double> valueListTime;
            timeFunc->set( &(valueList[ival]), &timeList, ival );
            valueList[ival].clear();
          }
          /*
no matching function for call to 'cseis_geolib::csTimeFunction<double>::set(cseis_geolib::csVector<double>&, cseis_geolib::csVector<double>*, 
no matching function for call to 'cseis_geolib::csTimeFunction<double>::set(cseis_geolib::csVector<double>&, 
                                                                            cseis_geolib::csVector<T>*, const cseis_geolib::csVector<double>*
           */

          //          valueListTime.clear();
          timeFunctionList.insertEnd( timeFunc );
          timeList.clear();

          double* keysTMP = new double[myNumAllKeys];
          memcpy( keysTMP, keysCurrent, myNumAllKeys*sizeof(double) );
          keyValueList.insertEnd( keysTMP );
          memcpy( keysCurrent, keysNew, myNumAllKeys*sizeof(double) );          
        }
        else {
          memcpy( keysCurrent, keysNew, myNumAllKeys*sizeof(double) );
        }
      }
      else if( myTableType == TABLE_TYPE_UNIQUE_KEYS ) {
        throw( csException("Input table file '%s' has two lines with same key values, at line #%d.",
                           myFilename.c_str(), counterLines ) );
      }
    } // End: Extract key values
    
    for( int ival = 0; ival < myNumValues; ival++ ) {
      valueList[ival].insertEnd( atof(tokenList.at( myValueColumns[ival] ).c_str()) );
    }
    if( myTableType == TABLE_TYPE_TIME_FUNCTION ) {
      timeList.insertEnd( atof(tokenList.at(myIndexTimeCol).c_str()) );
      //      valueListTime.insertEnd( atof(tokenList.at( myValueColumns[0] ).c_str()) );
    }

    counterLines++;
  }  // end while reading lines from input file

  //  if( myTableType != TABLE_TYPE_TIME_FUNCTION ) {
    if( valueList[0].size() == 0 ) {
      throw( csException("No valid line found in input table.\nWrong column numbers specified?\nMax number of columns found in table = %d.\nMax column number specified for key(s)/value(s) = %d", maxTableColumns, maxColumnIndex ) );
    }
    //  }
    // else {
    //  if( valueListTime.size() == 0 ) {
    //    throw( csException("No valid line found in input table.\nWrong column numbers specified?\nMax number of columns found in table = %d.\nMax column number specified for key(s)/value(s) = %d", maxTableColumns, maxColumnIndex ) );
    //  }
    //  }

  fclose( myFile );
  myFile = NULL;

  //-------------------------------------------------------
  //
  csSortManager sortManager( myNumAllKeys, csSortManager::TREE_SORT );

  myNumLocations = keyValueList.size();
  if( myNumAllKeys > 0 ) {
    if( myTableType == TABLE_TYPE_TIME_FUNCTION && timeList.size() != 0 ) {
      double* keysTMP = new double[myNumAllKeys];
      memcpy( keysTMP, keysCurrent, myNumAllKeys*sizeof(double) );
      keyValueList.insertEnd( keysTMP );
      myNumLocations += 1;
    }
    myKeyValues = new double*[myNumAllKeys];
    for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
      myKeyValues[ikey] = new double[myNumLocations];
    }
    if( !doSort ) {
      for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
        for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
          myKeyValues[ikey][iloc] = keyValueList.at(iloc)[ikey];
        }
      }
    }
    else {
      sortManager.resetValues( myNumLocations );
      for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
        for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
          double value =  keyValueList.at(iloc)[ikey];
          sortManager.setValue( iloc, myNumAllKeys-ikey-1, csFlexNumber(value) );
        }
      }
      sortManager.sort();

      for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
        for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
          int locIndex = sortManager.sortedIndex(iloc);
          myKeyValues[ikey][iloc] = keyValueList.at(locIndex)[ikey];;
        }
      }
    }
  }
  else {
    myNumLocations = 1;
  }

  if( myTableType == TABLE_TYPE_DUPLICATE_KEYS ) {
    addData_internal( valueList );
  }
  else if( myTableType == TABLE_TYPE_UNIQUE_KEYS ) {
    for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
      int locIndex = iloc;
      if( doSort ) locIndex = sortManager.sortedIndex(iloc);
      csTableValueList* tvl = new csTableValueList( myNumValues, 1 );
      for( int ival = 0; ival < myNumValues; ival++ ) {
        tvl->set( ival, 0, valueList[ival].at(locIndex) );
      }
      myValues->insertEnd( tvl );
    }
  }
  else {
    if( timeList.size() != 0 ) {
      csTimeFunction<double>* timeFunc = new csTimeFunction<double>( myNumValues );
      for( int ival = 0; ival < myNumValues; ival++ ) {
        timeFunc->set( &(valueList[ival]), &timeList, ival );
      }
      timeFunctionList.insertEnd( timeFunc );
    }

    myTimeFunctions2D = new csTimeFunction<double>* [myNumLocations];
    for( int i = 0; i < myNumLocations; i++ ) {
      myTimeFunctions2D[i] = timeFunctionList.at(i);
    }
    myCurrentTimeFunction->resetNumSpatialValues( myNumValues );
  }

  if( valueList != NULL ) delete [] valueList;
  if( keysCurrent ) {
    delete [] keysCurrent;
  }
  if( keysNew ) {
    delete [] keysNew;
  }

  if( myNumAllKeys > 0 ) {
    for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
      if( keyValueList.at(iloc) != NULL ) {
        delete [] keyValueList.at(iloc);
      }
    }
  }
}

//-----------------------------------------------------------------------------------
//
double csTableNew::getKeyValue( int indexLocation, int indexKey ) const {
  if( indexLocation >= 0 && indexLocation < myNumLocations && indexKey >= 0 && indexKey < myNumAllKeys ) {
    return myKeyValues[indexKey][indexLocation];
  }
  else {
    throw( csException("csTableNew::getKeyValue(): Incorrect location or key index passed") );
  }
}

//-----------------------------------------------------------------------------------
// Find exact key values (non-interpolated keys only)
// Assumes that keys are all sorted with increasing order
// Reminder: Non-interpolated keys come first, in fields keyValues_in and myKeyValues 
//
bool csTableNew::findExactKeyLocation( double const* keyValues_in, int& locStart, int& locEnd ) const {
  locStart = -1;
  locEnd   = -1;
  int startLocationIndex = 0;
  for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
    double currentKeyValue = keyValues_in[ikey];
    for( int location = startLocationIndex; location < myNumLocations; location++ ) {
      if( currentKeyValue <= myKeyValues[ikey][location] ) {
        if( currentKeyValue == myKeyValues[ikey][location] ) {      
          locStart = location;
          break;
        }
        else {
          return false;
        }
      }
    } // END for location
    if( locStart == -1 ) return false;
    startLocationIndex = locStart;
    // Check that this location still matches all previous keys:
    for( int ik2 = 0; ik2 < ikey; ik2++ ) {
      if( keyValues_in[ik2] != myKeyValues[ik2][locStart] ) {
        locStart = -1;
        return false;
      }
    }
  } // END for ikey

  // Find last location with matching keys:
  locEnd = locStart;
  for( int location = locStart+1; location < myNumLocations; location++ ) {
    for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
      if( keyValues_in[ikey] != myKeyValues[ikey][location] ) {
        return true;
      }
    }
    locEnd = location;
  }

  return true;
}
//-----------------------------------------------------------------------------------

bool csTableNew::findInterpKeyLocation( double keyValue_in, int keyIndex, int& locLeft, int& locRight, double& weight ) const {
  int counter = locRight - locLeft + 1;

  weight = 1.0;

  double valueLeft  = myKeyValues[keyIndex][locLeft];
  double valueRight = myKeyValues[keyIndex][locRight];

  int direction = valueRight > valueLeft ? 1 : -1;

  double diffLeft  = (double)direction * (keyValue_in - valueLeft);
  double diffRight = (double)direction * (keyValue_in - valueRight);

  // If key value is outside of given range, extrapolate with closest, constant value
  if( diffLeft <= 0 ) {
    locRight = locLeft;
    return false;
  }
  else if( diffRight >= 0 ) {
    locLeft = locRight;
    return false;
  }

  do {
    int locMid  = (locRight+locLeft)/2;
    double valueMid  = myKeyValues[keyIndex][locMid];
    double diff = (double)direction * (keyValue_in - valueMid);
    if( diff > 0 ) {
      locLeft  = locMid;
    }
    else {
      locRight = locMid;
    }
    if( locRight-locLeft == 1 ) {
      weight =
        ( keyValue_in - myKeyValues[keyIndex][locLeft] ) /
        ( myKeyValues[keyIndex][locRight] - myKeyValues[keyIndex][locLeft] );
      if( weight != 1.0 ) return false;
      else return true;
    }
    counter -= 1;
  } while( counter >= 0 );

  return false;
}

void csTableNew::findInterpKeyLocation2D( double const* keyValues_in,
                                          int locStart, int locEnd,
                                          int& locLeftUp, int& locLeftDown,
                                          int& locRightUp, int& locRightDown,
                                          double& weight1, double& weight2 ) const
{
  findInterpKeyLocation2D( keyValues_in,
                           locStart, locEnd,
                           locLeftUp, locLeftDown,
                           locRightUp, locRightDown,
                           weight1, weight2,
                           myNumKeys, myNumKeys+1 );
}

void csTableNew::findInterpKeyLocation2D( double const* keyValues_in,
                                          int locStart, int locEnd,
                                          int& locLeftUp, int& locLeftDown,
                                          int& locRightUp, int& locRightDown,
                                          double& weight1, double& weight2,
                                          int keyIndex1, int keyIndex2 ) const
{
  locLeftUp    = locStart;
  locRightDown = locEnd;
  findInterpKeyLocation( keyValues_in[keyIndex1], keyIndex1, locLeftUp, locRightDown, weight1 );

  // Save locations for later
  int loc1 = locLeftUp;
  int loc2 = locRightDown;
  // Widen the field of search for the second key: Extend left/right key locations to locations with same first key value
  while( locLeftUp > locStart && myKeyValues[keyIndex1][locLeftUp] == myKeyValues[keyIndex1][locLeftUp-1] ) {
    locLeftUp -= 1;
  }
  while( locRightDown < locEnd && myKeyValues[keyIndex1][locRightDown] == myKeyValues[keyIndex1][locRightDown+1] ) {
    locRightDown += 1;
  }

  // First interpolation key value is constant over the range of possible locations
  // --> Simply search for second key value
  if( myKeyValues[keyIndex1][locLeftUp] == myKeyValues[keyIndex1][locRightDown] ) {
    findInterpKeyLocation( keyValues_in[keyIndex2], keyIndex2, locLeftUp, locRightDown, weight2 );
    locRightUp   = locLeftUp;
    locLeftDown  = locRightDown;
  }
  else {
    locLeftDown = loc1;
    findInterpKeyLocation( keyValues_in[keyIndex2], keyIndex2, locLeftUp, locLeftDown, weight2 );

    locRightUp = loc2;
    findInterpKeyLocation( keyValues_in[keyIndex2], keyIndex2, locRightUp, locRightDown, weight2 );
  }
}


//-----------------------------------------------------------------------------------
//
double csTableNew::getValue( int indexLocation, int indexValue ) const {
  if( indexValue >= 0 && indexValue < myNumValues && indexLocation >= 0 && indexLocation < myNumLocations ) {
    return myValues->at(indexLocation)->get(indexValue);  // Get first value
  }
  else {
    throw( csException("csTableNew::getValue(): Incorrect value or location index passed") );
  }
}
//-----------------------------------------------------------------------------------
//
double csTableNew::getValue( double const* keyValues_in, int indexValue ) const {
  if( myNumAllKeys > 0 ) {
    if( myTableType == TABLE_TYPE_UNIQUE_KEYS ) {
      return interpolate( indexValue, keyValues_in );
    }
    else {
      int indexLocation = 0;
      int indexEnd = myNumLocations-1;
      if( findExactKeyLocation( keyValues_in, indexLocation, indexEnd ) ) {
        return myValues->at(indexLocation)->get( indexValue );  // Return first value
      }
    }
  }
  else {
    throw( csException("csTableNew::getValue(): Table does not contain any keys") );
  }
  return 0.0;
}
//-----------------------------------------------------------------------------------
csTableValueList const* csTableNew::getValues( int indexLocation ) const {
  if( myTableType == TABLE_TYPE_DUPLICATE_KEYS ) {
    return myValues->at(indexLocation);
  }
  else {
    throw( csException("csTableNew::getValues(): Incorrect use of this function for current table type. This is a program bug in the calling function") );
  }
}
//-----------------------------------------------------------------------------------
//
csTableValueList const* csTableNew::getValues( double const* keyValues_in ) const {
  if( myTableType == TABLE_TYPE_DUPLICATE_KEYS ) {
    int indexLocation = 0;
    int indexEnd = myNumLocations-1;
    if( findExactKeyLocation( keyValues_in, indexLocation, indexEnd ) ) {
      return getValues( indexLocation );
    }
    return NULL;
  }
  else {
    throw( csException("csTableNew::getValues(): Incorrect use of this function for current table type. This is a program bug in the calling function") );
  }
}
//-----------------------------------------------------------------------------------
//
void csTableNew::addData_internal( csVector<double>* valueList ) {
  csTableValueList* tvl = new csTableValueList( myNumValues, valueList[0].size() );
  for( int ival = 0; ival < myNumValues; ival++ ) {
    for( int ie = 0; ie < tvl->numLines(); ie++ ) {
      tvl->set( ival, ie, valueList[ival].at(ie) );
      //      printf(" Value #%d/%d:  %f\n", ival+1, ie+1, valueList[ival].at(ie) );
    }
    valueList[ival].clear();
  }
  myValues->insertEnd( tvl );
}
//-----------------------------------------------------------------------------------
//
double csTableNew::interpolate( int valueIndex, double const* keyValues_in ) const {
  // If there are keys that shall not be interpolated, search for these first and narrow down the locations:
  int locStart = 0;
  int locEnd   = myNumLocations-1;
  if( myNumKeys > 0 ) {
    if( !findExactKeyLocation( keyValues_in, locStart, locEnd ) ) {
      throw( csException("Exact key location not found. First key value: %f, number of non-interpolated keys: %d. (Possible solution: Interpolate based on key)",
                         keyValues_in[0], myNumKeys) );
    }
  }

  if( myNumInterpKeys < 2 ) {
    int locLeft   = locStart;
    int locRight  = locEnd;
    double weight = 1.0;
    // Table has one key that needs to be interpolated
    if( myNumInterpKeys > 0 ) findInterpKeyLocation( keyValues_in[myNumKeys], myNumKeys, locLeft, locRight, weight );

    double valLeft  = myValues->at(locLeft)->get( valueIndex );
    double valRight = myValues->at(locRight)->get( valueIndex );
    return( valLeft + weight * ( valRight - valLeft ) );
  }
  else {
    int locLeftUp    = 0;
    int locLeftDown  = 0;
    int locRightUp   = 0;
    int locRightDown = 0;
    double weightLeft;
    double weightRight;

    findInterpKeyLocation2D( keyValues_in, locStart, locEnd, locLeftUp, locLeftDown, locRightUp, locRightDown, weightLeft, weightRight );

    int keyIndex1 = myNumKeys;
    int keyIndex2 = myNumKeys+1;

    double result = interpolateStep2( keyValues_in, keyIndex1, keyIndex2, valueIndex,
                                      locLeftUp, locLeftDown, locRightUp, locRightDown );

    return result;
  }
}

double csTableNew::interpolateStep2( double const* keyValues_in, int keyIndex1, int keyIndex2, int valueIndex,
                                     int locLeftUp, int locLeftDown, int locRightUp, int locRightDown ) const
{
  double x1 = std::min( myKeyValues[keyIndex1][locLeftUp], myKeyValues[keyIndex1][locLeftDown] );
  double x2 = std::max( myKeyValues[keyIndex1][locRightUp], myKeyValues[keyIndex1][locRightDown] );

  double y1 = std::min( myKeyValues[keyIndex2][locLeftUp], myKeyValues[keyIndex2][locRightUp] );
  double y2 = std::max( myKeyValues[keyIndex2][locLeftDown], myKeyValues[keyIndex2][locRightDown] );

  double dx12 = x2 - x1;
  double dy12 = y2 - y1;
  double dx1 = keyValues_in[keyIndex1] - x1;
  double dx2 = x2 - keyValues_in[keyIndex1];
  double dy1 = keyValues_in[keyIndex2] - y1;
  double dy2 = y2 - keyValues_in[keyIndex2];
  
  double result = 0.0;
  if( dx12 == 0.0 ) {
    if( dy12 == 0.0 ) return myValues->at(locLeftUp)->get( valueIndex );
    result = myValues->at(locLeftUp)->get( valueIndex ) + dy1/dy12 *
      ( myValues->at(locLeftDown)->get( valueIndex ) - myValues->at(locLeftUp)->get( valueIndex ) );
  }
  else if( dy12 == 0.0 ) {
    result = myValues->at(locLeftUp)->get( valueIndex ) + dx1/dx12 *
      ( myValues->at(locRightUp)->get( valueIndex ) - myValues->at(locLeftUp)->get( valueIndex ) );
  }
  else { 
    result += dx2 * dy2 * myValues->at(locLeftUp)->get( valueIndex );
    result += dx1 * dy2 * myValues->at(locRightUp)->get( valueIndex );
    result += dx2 * dy1 * myValues->at(locLeftDown)->get( valueIndex );
    result += dx1 * dy1 * myValues->at(locRightDown)->get( valueIndex );
    result *= 1.0 / (dx12 * dy12);
  }
  return result;
}
//-----------------------------------------------------------------------------------
//
void csTableNew::dump() const {
  fprintf(stdout,"Table file name:  %s\n", myFilename.c_str() );
  fprintf(stdout,"Number of keys:   %d  (= %d + %d )\n", myNumAllKeys, myNumKeys, myNumInterpKeys );
  fprintf(stdout,"Number of values: %d\n", myNumValues );
  fprintf(stdout,"Num locations:    %d\n", myNumLocations);

  if( myTableType == TABLE_TYPE_DUPLICATE_KEYS ) {
    for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
      csTableValueList* tvl = myValues->at(iloc);
      fprintf(stdout,"Key %12f:\n", getKeyValue(iloc,0) );
      for( int iline = 0; iline < tvl->numLines(); iline++ ) {
        fprintf(stdout," Line %3d: ", iline+1 );
        for( int icol = 0; icol < tvl->numColumns(); icol++ ) {
          fprintf(stdout," %12f", tvl->get(icol,iline) );
        }
        fprintf(stdout,"\n");
      }
    }
  }
  else if( myTableType == TABLE_TYPE_UNIQUE_KEYS ) {
    for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
      csTableValueList* tvl = myValues->at(iloc);
      for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
        fprintf(stdout,"%12f ", getKeyValue(iloc,ikey) );
      }
      if( tvl->numLines() == 1 ) {
        for( int icol = 0; icol < tvl->numColumns(); icol++ ) {
          fprintf(stdout," %12f", tvl->get(icol) );
        }
      }
      else {
        fprintf(stdout,":\n" );
        for( int iline = 0; iline < tvl->numLines(); iline++ ) {
          fprintf(stdout," Line %3d: ", iline+1 );
          for( int icol = 0; icol < tvl->numColumns(); icol++ ) {
            fprintf(stdout," %12f", tvl->get(icol,iline) );
          }
        }
      }
      fprintf(stdout,"\n");
    }
  }
  else if( myTableType == TABLE_TYPE_TIME_FUNCTION ) {
    double* keyValues = new double[myNumAllKeys];

    for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
      fprintf(stdout,"Keys at loc #%d: ", iloc);
      for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
        keyValues[ikey] = getKeyValue( iloc, ikey );
        fprintf(stdout," %12f", keyValues[ikey] );
      }
      fprintf(stdout,"\n");
      csTimeFunction<double> const* timeFunction = getFunction( keyValues );
      int numValues = timeFunction->numValues();
      int numSpatialValues = timeFunction->numSpatialValues();
      fprintf(stdout,"Number of spatial values/values in function: %d / %d\n", numSpatialValues, numValues );
      for( int ival = 0; ival < numValues; ival++ ) {
        double time  = timeFunction->timeAtIndex( ival );
        fprintf(stdout," Time/Values #%-3d  %12f ", ival+1, time );
        for( int is = 0; is < numSpatialValues; is++ ) {
          double value = timeFunction->valueAtIndex( ival, is );
          fprintf(stdout,"%12f ", value );
        }
        fprintf(stdout,"\n" );
      }
    }
    delete [] keyValues;
  }
}


//********************************************************************************
//********************************************************************************
// Time function
//

//-----------------------------------------------------------------------------------
//
//
csTimeFunction<double> const* csTableNew::getFunction( double const* keyValues_in, bool dump ) const {
  int locLeft    = 0;
  int locRight   = myNumLocations-1;
  double weight1 = 1.0;

  csTimeFunction<double> const* timeFuncLeft  = myTimeFunctions2D[locLeft];  // Time function on 'left' side of requested location (key1 is smaller)
  csTimeFunction<double> const* timeFuncRight = myTimeFunctions2D[locRight]; // Time function on 'right' side of requested location (key1 is larger)
  csTimeFunction<double> tempTimeFuncLeft( timeFuncLeft->numSpatialValues() );
  csTimeFunction<double> tempTimeFuncRight( timeFuncLeft->numSpatialValues() );

  if( myNumAllKeys > 0 ) {
    if( myNumInterpKeys == 0 ) { // No interpolation, need to find excact key
      if( !findInterpKeyLocation( keyValues_in[myNumKeys-1], 0, locLeft, locRight, weight1 ) ) {
        throw( csException("Key value not found: %f. Suggest to specify to interpolate key value.\n", keyValues_in[myNumKeys-1]) );
      }
      timeFuncLeft  = myTimeFunctions2D[locLeft];
      timeFuncRight = myTimeFunctions2D[locRight];
      //      fprintf(stderr,"Locations: %f  %d %d   %f\n", keyValues_in[0], locLeft, locRight, weight1 );
    }
    else if( myNumInterpKeys == 1 ) {
      findInterpKeyLocation( keyValues_in[myNumKeys], 0, locLeft, locRight, weight1 );
      timeFuncLeft  = myTimeFunctions2D[locLeft];
      timeFuncRight = myTimeFunctions2D[locRight];
      //      fprintf(stderr,"Locations: %f %f  %d %d\n", keyValues_in[0], keyValues_in[1], locLeft, locRight );
    }
    else { // More than one interpolated key --> Need 2D interpolation
      int locLeftUp;
      int locLeftDown;
      int locRightUp;
      int locRightDown;
      double weight2;

      // Step 1: Find four time functions surrounding the requested location
      // Left: key1 smaller, Right: key1 larger, Up: key2 smaller, Down: key2 larger
      findInterpKeyLocation2D( keyValues_in, locLeft, locRight, locLeftUp, locLeftDown, locRightUp, locRightDown, weight1, weight2 );

      timeFuncLeft  = &tempTimeFuncLeft;  // Set link to internal, temporary time function, for last interpolation step
      timeFuncRight = &tempTimeFuncRight;

      // Step 2a: Interpolate between up/down locations on left hand side, along key2 axis
      interpolateTimeFunction( myTimeFunctions2D[locLeftUp], myTimeFunctions2D[locLeftDown], weight2, &tempTimeFuncLeft );
      // Step 2b: Interpolate between up/down locations on right hand side, along key2 axis
      interpolateTimeFunction( myTimeFunctions2D[locRightUp], myTimeFunctions2D[locRightDown], weight2, &tempTimeFuncRight );
      
    } // END more than one interpolated key
  } // END if myNumAllKeys > 0

  // Step X: Interpolate between left and right sides, along key1 axis
  interpolateTimeFunction( timeFuncLeft, timeFuncRight, weight1, myCurrentTimeFunction );

  if( dump ) {
    for(int i = 0; i < myCurrentTimeFunction->numValues(); i++ ) {
      for( int ikey = 0; ikey < myNumAllKeys; ikey++ ) {
        fprintf(stdout,"%f ", keyValues_in[ikey] );
      }
      fprintf(stdout,"%f %f\n",myCurrentTimeFunction->timeAtIndex( i ), myCurrentTimeFunction->valueAtIndex( i ) );
    }
  }

  return myCurrentTimeFunction;
}

void csTableNew::interpolateTimeFunction( csTimeFunction<double> const* timeFuncLeft, csTimeFunction<double> const* timeFuncRight,
                                          double weightLoc, csTimeFunction<double>* newTimeFunction ) const
{
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
  for( int ival = 0; ival < timeFuncLeft->numSpatialValues(); ival++ ) {
    for( int itime = 0; itime < numTimes; itime++ ) {
      double time = timeList.at(itime);
      double valueLeft = timeFuncLeft->valueAt( time, ival );
      double valueRight= timeFuncRight->valueAt( time, ival );
      
      valueList.insertEnd( valueLeft + weightLoc * ( valueRight - valueLeft ) );
    }
    newTimeFunction->set( &valueList, &timeList, ival );
    valueList.clear();
  }

}

//-----------------------------------------------------------------------------------
//
double csTableNew::getTimeValue( double const* keyValues, double time ) const {
  if( myNumKeys > 0 || myTableType == TABLE_TYPE_TIME_FUNCTION ) {
    return interpolateTimeFunction( keyValues, time );
  }
  else {
    throw( csException("csTableNew::getTimeValue(): Table does not contain any keys, or incorrect use of non-timefunction table") );
  }
}

//-----------------------------------------------------------------------------------
//
double csTableNew::interpolateTimeFunction( double const* keyValues_in, double time ) const {
  int locLeft = 0;
  int locRight = myNumLocations-1;
  double weightLoc;

  findInterpKeyLocation( keyValues_in[myNumKeys], 0, locLeft, locRight, weightLoc );
  //  findExactKeyLocation( keyValues_in, locLeft, locRight, weightLoc );

  double valueLeft  = myTimeFunctions2D[locLeft]->valueAt( time );
  double valueRight = myTimeFunctions2D[locRight]->valueAt( time );

  return( valueLeft + weightLoc * ( valueRight - valueLeft ) );
}



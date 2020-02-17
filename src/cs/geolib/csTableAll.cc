/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csTableAll.h"
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

csTableAll::csTableAll( int tableType ) {
  myValues     = NULL;
  myKeyValues  = NULL;
  myKeyNames   = NULL;
  myValueNames = NULL;
  myColumnTypes = NULL;
  myFile       = NULL;
  myNumKeys    = 0;
  myNumValues  = 0;
  myTableType  = tableType;

  myNumCols        = 0;
  myNumLocations   = 0;
  myIndexFirstValueCol = 0;
  myHasReadTableContents = false;
  myHasBeenInitialized   = false;

  myTimeFunctions2D = NULL;
  myIndexTimeCol    = 0;
  myCurrentTimeFunction = NULL;
  if( myTableType == TABLE_TYPE_TIME_FUNCTION ) {
    myCurrentTimeFunction = new csTimeFunction<double>();
  }
  else {
    myValues = new csVector<TableValueList*>();
  }
}
csTableAll::~csTableAll() {
  clearBuffers();
  if( myValues != NULL ) {
    for( int ie = 0; ie < myValues->size(); ie++ ) {
      delete myValues->at(ie);
    }
    delete myValues;
    myValues = NULL;
  }
  if( myColumnTypes != NULL ) {
    delete [] myColumnTypes;
    myColumnTypes = NULL;
  }
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
}
void csTableAll::clearBuffers() {
  if( myKeyValues != NULL ) {
    for( int i = 0; i < myNumKeys; i++ ) {
      delete [] myKeyValues[i];
    }
    delete [] myKeyValues;
    myKeyValues = NULL;
  }
}

//-----------------------------------------------------------------------------------
//
//
void csTableAll::initialize( std::string const& filename, int methodInterpolation ) {
  if( myHasBeenInitialized ) {
    throw( csException("Table has been initialized before. Previous input file name: '%s', new file name: %s", myFilename.c_str(), filename.c_str()) );
  }
  myTablename = "Table";
  myFilename  = filename;
  myMethodInterpolation = methodInterpolation;

  clearBuffers();
  char buffer[1024];

  if( (myFile = fopen( filename.c_str(), "r" )) == (FILE*) NULL ) {
    throw csException("Could not open file: '%s'", filename.c_str());
  }
  
  if( fgets( buffer, 1024, myFile ) == NULL ) {
    throw( csException("Table input file '%s' contains no valid lines. First line: %s", myFilename.c_str(), buffer ) );
  }
  // Convert names of keys/time/values to lower case
  toLowerCase( buffer );

  // 1) Set key names, and determine number of keys, values and time axis
  csVector<std::string> tokenList;
  tokenize( buffer, tokenList );
  myNumCols = tokenList.size();
  if( myNumCols == 0 ) {
    throw( csException("Table input file '%s' contains no valid lines. First line: %s", myFilename.c_str(), buffer ) );
  }
  else if( myNumCols < 2 ) {
    throw( csException("Table input file '%s' contains no key or value column", myFilename.c_str()) );
  }

  myNumKeys = 0;
  while( myNumKeys < myNumCols && tokenList.at(myNumKeys)[0] == KEY_CHAR ) {
    myNumKeys += 1;
  }
  //  if( myNumKeys == 0 ) {
  //  throw( csException("Table input file '%s' contains no key column. Key column names are preceded by the character '%c', e.g. %csource",
  //                     myFilename.c_str(), KEY_CHAR, KEY_CHAR) );    
  // }
  if( myNumKeys == myNumCols ) {
    throw( csException("Table input file '%s' contains no value column. Names of value columns must NOT be preceded by the key character '%c'",
                       myFilename.c_str(), KEY_CHAR) );    
  }

  std::string token = tokenList.at(myNumKeys);
  int addTimeColumn = 0;
  if( myTableType == TABLE_TYPE_TIME_FUNCTION ) {
    if( token.compare("time") ) {
      throw( csException("Table input file '%s' contains no 'time' column ('time' NOT preceded by the key character '%c'). This is required for the requested table type (time function table)",
                         myFilename.c_str(), KEY_CHAR ) );
    }
    myNumValues    = myNumCols-myNumKeys-1;
    addTimeColumn  = 1;
    myIndexTimeCol = myNumKeys;
    if( myNumValues > 1 ) {
      throw( csException("Table input file '%s' contains %d value columns. Only one value column is currently supported for time function tables",
                         myFilename.c_str(), myNumValues ) );    
    }
  }
  else {
    myNumValues = myNumCols-myNumKeys;
  }

  if( myNumValues == myNumCols ) {
    throw( csException("Table input file '%s' contains no key column.", myFilename.c_str()) );    
  }

  if( myNumKeys > 0 ) {
    myKeyNames = new std::string[myNumKeys];
    for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
      int length = tokenList.at(ikey).length();
      myKeyNames[ikey] = tokenList.at(ikey).substr(1,length-1);
    }
  }
  myValueNames = new std::string[myNumValues];
  for( int ival = 0; ival < myNumValues; ival++ ) {
    //    int length = tokenList.at(myNumKeys+ival).length();
    myValueNames[ival] = tokenList.at(myNumKeys+ival+addTimeColumn);
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
void csTableAll::readTableContents( bool doSort ) {
//void csTableAll::readTableContents( type_t const* columnTypes, int numCols ) {

  // 2) Read in all values, key values, time values...
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
  // Helper field for storing tokens extracted from input line
  csVector<std::string> tokenList;
  // List of values for each 'value' column in input file
  csVector<double>* valueList = new csVector<double>[myNumValues];
  // List of key values for each 'location'. Each list item is an array of key values (one value for each key)
  csVector<double*> keyValueList;

  csVector<double> timeList;
  csVector<double> valueListTime;
  csVector<csTimeFunction<double>*> timeFunctionList;

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
      }
      if( !isSame ) {  // Key value in current line differs from key value in previous line --> Store values up to this point under previous key
        if( myTableType != TABLE_TYPE_TIME_FUNCTION ) {
          double* keysTMP = new double[myNumKeys];
          memcpy( keysTMP, keysNew, myNumKeys*sizeof(double) );
          keyValueList.insertEnd( keysTMP );
          memcpy( keysCurrent, keysNew, myNumKeys*sizeof(double) );
          if( myTableType == TABLE_TYPE_DUPLICATE_KEYS && valueList[0].size() != 0 ) addData_internal( valueList );
        }
        else if ( timeList.size() != 0 ) {
          csTimeFunction<double>* timeFunc = new csTimeFunction<double>();
          timeFunc->set( &valueListTime, &timeList );
          timeFunctionList.insertEnd( timeFunc );
          valueListTime.clear();
          timeList.clear();

          double* keysTMP = new double[myNumKeys];
          memcpy( keysTMP, keysCurrent, myNumKeys*sizeof(double) );
          keyValueList.insertEnd( keysTMP );
          memcpy( keysCurrent, keysNew, myNumKeys*sizeof(double) );          
        }
        else {
          memcpy( keysCurrent, keysNew, myNumKeys*sizeof(double) );
        }
      }
      else if( myTableType == TABLE_TYPE_UNIQUE_KEYS ) {
        throw( csException("Input table file '%s' has two lines with same key values, at line #%d.",
                           myFilename.c_str(), counterLines ) );
      }
    } // End: Extract key values
    
    if( myTableType != TABLE_TYPE_TIME_FUNCTION ) {
      for( int ival = 0; ival < myNumValues; ival++ ) {
        valueList[ival].insertEnd( atof(tokenList.at(myIndexFirstValueCol+ival).c_str()) );
      }
    }
    else {
      timeList.insertEnd( atof(tokenList.at(myIndexTimeCol).c_str()) );
      valueListTime.insertEnd( atof(tokenList.at(myIndexFirstValueCol).c_str()) );
    }

    counterLines++;
  }  // end while reading lines from input file

  fclose( myFile );
  myFile = NULL;

  //-------------------------------------------------------
  //
  csSortManager sortManager( myNumKeys, csSortManager::TREE_SORT );

  myNumLocations = keyValueList.size();
  if( myNumKeys > 0 ) {
    if( myTableType == TABLE_TYPE_TIME_FUNCTION && timeList.size() != 0 ) {
      double* keysTMP = new double[myNumKeys];
      memcpy( keysTMP, keysCurrent, myNumKeys*sizeof(double) );
      keyValueList.insertEnd( keysTMP );
      myNumLocations += 1;
    }
    myKeyValues = new double*[myNumKeys];
    for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
      myKeyValues[ikey] = new double[myNumLocations];
    }
    if( !doSort ) {
      for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
        for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
          myKeyValues[ikey][iloc] = keyValueList.at(iloc)[ikey];
        }
      }
    }
    else {
      sortManager.resetValues( myNumLocations );
      for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
        for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
          double value =  keyValueList.at(iloc)[ikey];
          sortManager.setValue( iloc, myNumKeys-ikey-1, csFlexNumber(value) );
        }
      }
      sortManager.sort();

      for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
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
      TableValueList* tvl = new TableValueList( myNumValues, 1 );
      for( int ival = 0; ival < myNumValues; ival++ ) {
        tvl->set( ival, 0, valueList[ival].at(locIndex) );
      }
      myValues->insertEnd( tvl );
    }
  }
  else {
    if( timeList.size() != 0 ) {
      csTimeFunction<double>* timeFunc = new csTimeFunction<double>();
      timeFunc->set( &valueListTime, &timeList );
      timeFunctionList.insertEnd( timeFunc );
    }

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

  for( int iloc = 0; iloc < keyValueList.size(); iloc++ ) {
    delete [] keyValueList.at(iloc);
  }
}
  

//-----------------------------------------------------------------------------------
//

/*double const* csTableAll::getKeyValues( int indexLocation ) const {
  if( indexLocation >= 0 && indexLocation < myNumLocations && myNumKeys > 0 ) {
    return myKeyValues[indexLocation];
  }
  else {
    return NULL;
  }
}
*/
//-----------------------------------------------------------------------------------
//
std::string const csTableAll::keyName( int indexKey ) const {
  if( indexKey >= 0 && indexKey < myNumKeys ) {
    return myKeyNames[indexKey];
  }
  else {
    throw( csException("csTableAll::keyName(): Incorrect key index passed") );
  }
}
//-----------------------------------------------------------------------------------
//
std::string const csTableAll::valueName( int indexValue ) const {
  if( indexValue >= 0 && indexValue < myNumValues ) {
    return myValueNames[indexValue];
  }
  else {
    throw( csException("csTableAll::valueName(): Incorrect value index passed") );
  }
}

//-----------------------------------------------------------------------------------
//
double csTableAll::getKeyValue( int indexLocation, int indexKey ) const {
  if( indexLocation >= 0 && indexLocation < myNumLocations && indexKey >= 0 && indexKey < myNumKeys ) {
    return myKeyValues[indexKey][indexLocation];
  }
  else {
    throw( csException("csTableAll::getKeyValue(): Incorrect location or key index passed") );
  }
}

//-----------------------------------------------------------------------------------
// Assumes that keys are all sorted with increasing order
// Not implemented fully! Deals with single key only
bool csTableAll::findKeyLocation( double keyValue_in, int& locationIndex ) const {
  if( myNumKeys > 1 ) {
    throw( csException("csTableAll::...(): Table contains more than 1 key. This is not supported yet.") );
  }
  locationIndex = -1;
  int ikey = 0;
  double currentKeyValue = keyValue_in;
  for( int location = myNumLocations-1; location >= 0; location-- ) {
    if( currentKeyValue >= myKeyValues[ikey][location] ) {
      if( currentKeyValue == myKeyValues[ikey][location] ) {      
        locationIndex = location;
        return true;
      }
      else {
        return false;
      }
    }
  }
  return false;
}
//-----------------------------------------------------------------------------------
void csTableAll::findKeyLocation( double keyValue_in, int& locLeft, int& locRight, double& weight, int keyIndex ) const {
  if( keyIndex < 0 && myNumKeys > 1 ) {
    throw( csException("csTable::findKeyLocation: Table contains more than 1 key. This is not supported yet.") );
  }
  int ikey = keyIndex;
  if( ikey < 0 ) ikey = 0;
  double currentKeyValue = keyValue_in;

  weight = 1.0;
  locLeft  = 0;
  locRight = myNumLocations-1;
  double valueLeft  = myKeyValues[ikey][locLeft];
  double valueRight = myKeyValues[ikey][locRight];

  int direction = valueRight > valueLeft ? 1 : -1;

  double diffLeft  = (double)direction * (currentKeyValue - valueLeft);
  double diffRight = (double)direction * (currentKeyValue - valueRight);

  if( diffLeft <= 0 ) {
    locRight = locLeft;
    return;
  }
  else if( diffRight >= 0 ) {
    locLeft = locRight;
    return;
  }

  int counter = myNumLocations;
  do {
    int locMid  = (locRight+locLeft)/2;
    double valueMid  = myKeyValues[ikey][locMid];
    double diff = (double)direction * (currentKeyValue - valueMid);
    if( diff > 0 ) {
      locLeft  = locMid;
    }
    else {
      locRight = locMid;
    }
    if( locRight-locLeft == 1 ) {
      weight = (currentKeyValue-myKeyValues[ikey][locLeft])/(myKeyValues[ikey][locRight]-myKeyValues[ikey][locLeft]);
      return;
    }
    counter -= 1;
  } while( counter >= 0 );

  throw(csException("osTableAll::...() Error occurred: Unable to determine key values..."));
}

//-----------------------------------------------------------------------------------
//
void csTableAll::findKeyLocation2D( double const* keyValues_in, int& locLeftUp, int& locLeftDown, int& locRightUp, int& locRightDown ) const {
  if( myNumKeys > 2 ) {
    throw( csException("csTable::findKeyLocation2D: Table contains more than 2 keys. This is not supported yet.") );
  }
  double weightDummy;
  int loc1;
  int loc2;
  findKeyLocation( keyValues_in[0], loc1, loc2, weightDummy, 0 );

  //  fprintf(stderr,"Location1D  %d %d   %f %f\n", loc1, loc2, keyValues_in[0], keyValues_in[1] );
  
  if( myKeyValues[0][loc1] == myKeyValues[0][loc2] ) {
    while( loc2 > 0 && myKeyValues[0][loc2] == myKeyValues[0][loc2-1] ) {
      loc2 -= 1;
    }
    while( loc2 < myNumLocations-1 && myKeyValues[1][loc2] < keyValues_in[1] ) {
      if( myKeyValues[0][loc2] != myKeyValues[0][loc2+1] ) break;
      loc2 += 1;
    }
    loc1 = loc2;
    while( loc1 > 0 && myKeyValues[1][loc1] > keyValues_in[1] ) {
      if( myKeyValues[0][loc1] != myKeyValues[0][loc1-1] ) break;
      loc1 -= 1;
    }

    locLeftUp    = loc1;
    locRightUp   = loc1;
    locLeftDown  = loc2;
    locRightDown = loc2;
  }
  else {
    int locIndex = loc1;
    locLeftDown  = loc1;
    while( locIndex > 0 && myKeyValues[1][locIndex] > keyValues_in[1] ) {
      if( myKeyValues[0][locIndex] != myKeyValues[0][locIndex-1] ) break;
      locIndex    -= 1;
    }
    locLeftUp = locIndex;
    if( locIndex != loc1 ) locLeftDown = locIndex+1;
    if( myKeyValues[1][locLeftUp] == keyValues_in[1] ) locLeftDown = locLeftUp;

    locIndex    = loc2;
    locRightUp  = loc2;

    while( locIndex < myNumLocations-1 && myKeyValues[1][locIndex] < keyValues_in[1] ) {
      if( myKeyValues[0][locIndex] != myKeyValues[0][locIndex+1] ) break;
      locIndex    += 1;
      //    fprintf(stderr,"Increased to %d   %f\n", locIndex, keyValues_in[1]);
    }
    locRightDown = locIndex;
    if( locIndex != loc2 ) locRightUp = locIndex-1;
    if( myKeyValues[1][locRightDown] == keyValues_in[1] ) locRightUp = locRightDown;
  }
  
  //  fprintf(stderr,"Location  %d %d  %d %d\n", locLeftUp, locLeftDown, locRightUp, locRightDown );
  //  for( int ikey = 0; ikey < myNumKeys; ikey++ ) { 
  //  fprintf(stderr,"Location key #%d:  %f %f  %f %f\n", ikey, myKeyValues[ikey][locLeftUp], myKeyValues[ikey][locLeftDown],
  //          myKeyValues[ikey][locRightUp], myKeyValues[ikey][locRightDown] );
  // }
}

//-----------------------------------------------------------------------------------
//
double csTableAll::getValue( int indexLocation, int indexValue ) const {
  if( indexValue >= 0 && indexValue < myNumValues && indexLocation >= 0 && indexLocation < myNumLocations ) {
    return myValues->at(indexLocation)->get(indexValue);  // Get first value
  }
  else {
    throw( csException("csTableAll::getValue(): Incorrect value or location index passed") );
  }
}
//-----------------------------------------------------------------------------------
//
double csTableAll::getValue( double const* keyValues, int indexValue ) const {
  if( myNumKeys > 0 ) {
    if( myTableType == TABLE_TYPE_UNIQUE_KEYS ) {
      return interpolate( indexValue, keyValues );
    }
    else {
      int indexLocation = 0;
      if( findKeyLocation( keyValues[0], indexLocation ) ) {
        return myValues->at(indexLocation)->get( indexValue );  // Return first value
      }
    }
  }
  else {
    throw( csException("csTableAll::getValue(): Table does not contain any keys") );
  }
  return 0.0;
}
//-----------------------------------------------------------------------------------
TableValueList const* csTableAll::getValues( int indexLocation ) const {
  if( myTableType == TABLE_TYPE_DUPLICATE_KEYS ) {
    return myValues->at(indexLocation);
  }
  else {
    throw( csException("csTableAll::getValues(): Incorrect use of this function for current table type. This is a program bug in the calling function") );
  }
}
//-----------------------------------------------------------------------------------
//
TableValueList const* csTableAll::getValues( double const* keyValues ) const {
  if( myTableType == TABLE_TYPE_DUPLICATE_KEYS ) {
    int indexLocation = 0;
    if( findKeyLocation( keyValues[0], indexLocation ) ) {
      return getValues( indexLocation );
    }
    return NULL;
  }
  else {
    throw( csException("csTableAll::getValues(): Incorrect use of this function for current table type. This is a program bug in the calling function") );
  }
}
//-----------------------------------------------------------------------------------
//
void csTableAll::addData_internal( csVector<double>* valueList ) {
  TableValueList* tvl = new TableValueList( myNumValues, valueList[0].size() );
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
double csTableAll::interpolate( int valueIndex, double const* keyValues_in ) const {
  if( myNumKeys < 2 ) {
    int locLeft   = 0;
    int locRight  = 0;
    double weight = 1.0;
    if( myNumKeys > 0 ) findKeyLocation( keyValues_in[0], locLeft, locRight, weight );

    double valLeft  = myValues->at(locLeft)->get( valueIndex );
    double valRight = myValues->at(locRight)->get( valueIndex );
    return( valLeft + weight * ( valRight - valLeft ) );
  }
  else {
    int locLeftUp    = 0;
    int locLeftDown  = 0;
    int locRightUp   = 0;
    int locRightDown = 0;

    findKeyLocation2D( keyValues_in, locLeftUp, locLeftDown, locRightUp, locRightDown );

    double x1 = std::min( myKeyValues[0][locLeftUp], myKeyValues[0][locLeftDown] );
    double x2 = std::max( myKeyValues[0][locRightUp], myKeyValues[0][locRightDown] );

    double y1 = std::min( myKeyValues[1][locLeftUp], myKeyValues[1][locRightUp] );
    double y2 = std::max( myKeyValues[1][locLeftDown], myKeyValues[1][locRightDown] );

    double dx12 = x2 - x1;
    double dy12 = y2 - y1;
    double dx1 = keyValues_in[0] - x1;
    double dx2 = x2 - keyValues_in[0];
    double dy1 = keyValues_in[1] - y1;
    double dy2 = y2 - keyValues_in[1];

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
}
//-----------------------------------------------------------------------------------
//
void csTableAll::dump() const {
  fprintf(stdout,"Key names  (%d): ", myNumKeys );
  for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
    fprintf(stdout,"@%-12s ", keyName(ikey).c_str());
  }
  fprintf(stdout,"\n");
  fprintf(stdout,"Value names (%d): ", myNumValues );
  for( int ival = 0; ival < myNumValues; ival++ ) {
    fprintf(stdout,"%-12s ", valueName(ival).c_str());
  }
  fprintf(stdout,"\nNum locations: %d\n", myNumLocations);

  if( myTableType == TABLE_TYPE_DUPLICATE_KEYS ) {
    for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
      TableValueList* tvl = myValues->at(iloc);
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
      TableValueList* tvl = myValues->at(iloc);
      fprintf(stdout,"Key " );
      for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
        fprintf(stdout," %12f", getKeyValue(iloc,ikey) );
      }
      fprintf(stdout,":\n" );
      for( int iline = 0; iline < tvl->numLines(); iline++ ) {
        fprintf(stdout," Line %3d: ", iline+1 );
        for( int icol = 0; icol < tvl->numColumns(); icol++ ) {
          fprintf(stdout," %12f", tvl->get(icol,iline) );
        }
        fprintf(stdout,"\n");
      }
    }
  }
  else if( myTableType == TABLE_TYPE_TIME_FUNCTION ) {
    for( int i = 0; i < myNumLocations; i++ ) {
      fprintf(stdout,"++++++++++++ %d %f\n", i, myTimeFunctions2D[i]->valueAt(0));
    }
    double* keyValues = new double[myNumKeys];

    for( int iloc = 0; iloc < myNumLocations; iloc++ ) {
      fprintf(stdout,"Key ");
      for( int ikey = 0; ikey < myNumKeys; ikey++ ) {
        keyValues[ikey] = getKeyValue( ikey, iloc );
        fprintf(stdout," %12f", keyValues[ikey] );
      }
      fprintf(stdout,"\n");
      csTimeFunction<double> const* timeFunction = getFunction( keyValues );
      int numValues = timeFunction->numValues();
      for( int ival = 0; ival < numValues; ival++ ) {
        double value = timeFunction->valueAtIndex( ival );
        double time  = timeFunction->timeAtIndex( ival );
        fprintf(stdout," Time/Value #%-3d  %12f %12f\n", ival+1, time, value );
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
csTimeFunction<double> const* csTableAll::getFunction( double const* keyValues_in ) const {
  int locLeft    = 0;
  int locRight   = 0;
  double weightLoc = 1.0;

  if( myNumKeys > 0 ) findKeyLocation( keyValues_in[0], locLeft, locRight, weightLoc );

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
double csTableAll::getTimeValue( double const* keyValues, double time ) const {
  if( myNumKeys > 0 || myTableType == TABLE_TYPE_TIME_FUNCTION ) {
    return interpolateTimeFunction( keyValues, time );
  }
  else {
    throw( csException("csTableAll::getTimeValue(): Table does not contain any keys, or incorrect use for non-timefunction table") );
  }
}

//-----------------------------------------------------------------------------------
//
double csTableAll::interpolateTimeFunction( double const* keyValues_in, double time ) const {
  int locLeft;
  int locRight;
  double weightLoc;

  findKeyLocation( keyValues_in[0], locLeft, locRight, weightLoc );

  double valueLeft  = myTimeFunctions2D[locLeft]->valueAt( time );
  double valueRight = myTimeFunctions2D[locRight]->valueAt( time );

  return( valueLeft + weightLoc * ( valueRight - valueLeft ) );
}



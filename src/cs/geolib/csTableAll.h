/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TABLE_ALL_H
#define CS_TABLE_ALL_H

#include <cstdio>
#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {

  template<typename T> class csTimeFunction;
  template <typename T> class csVector;

  class TableValueList {
  public:
    TableValueList( int theNumCols, int theNumLines ) {
      myNumLines   = theNumLines;
      myNumCols    = theNumCols;
      values       = new double[myNumCols*myNumLines];
    }
    ~TableValueList() {
      if( values != NULL ) {
        delete [] values;
        values = NULL;
      }
    }
    void set( int indexCol, int indexLine, double value ) {
      values[indexCol*myNumLines+indexLine] = value;
    }
    double get( int indexCol, int indexLine = 0 ) const {
      return values[indexCol*myNumLines+indexLine];
    }
    double const* getAll( int indexCol ) const {
      return &values[indexCol*myNumLines];
    }
    int numColumns() const { return myNumCols; }
    int numLines() const { return myNumLines; }
  private:
    double* values;
    int myNumCols;
    int myNumLines;
  };


/**
* CSEIS ASCII Table - Base class
*
* CSEIS table given in ASCII file
*
---------------------------------
* This class provides functionality to...
*  - ...read in an ASCII file containing a table
*  - ...interpolate table values to any key location
* Supported types of tables are:
*  (a) Tables with 'unique key values'
*      - Each combination of key values (N keys) exists only once in table
*      - Values can be interpolated for key value combinations that are not explicitely given in table
*  (b) Tables with 'duplicate key values'
*      - Key values may occur more than once in table, i.e.
*        a list of values may be given for each key value combination
*      - Interpolation is not possible with such tables, due to ambiguity of which values to use
*  (c) Time function tables
*      - 
*
* Format of first line in input table:
* @name_key1  @name_key2 .. @name_keyN  [time]  name_val1  name_val2 .. name_valM
*  
* name_key?: Names of table keys
* name_val?: Names of table values
* time:      Optional time key for tables containing a time axis (e.g. velocity profiles)
--------------------------------
*
* @author Bjorn Olofsson
* @date 2007
*/

class csTableAll {
public:
  static const int TABLE_TYPE_UNKNOWN        = -1;
  static const int TABLE_TYPE_TIME_FUNCTION  = 1;
  static const int TABLE_TYPE_DUPLICATE_KEYS = 2;
  static const int TABLE_TYPE_UNIQUE_KEYS    = 3;
  static const int TABLE_TYPE_MIXED_KEYS     = 4;

  static const char KEY_CHAR = '@';

public:
  csTableAll( int tableType );
  virtual ~csTableAll();

  virtual void initialize( std::string const& filename, int methodInterpolation = 0 );
  virtual void readTableContents( bool doSort = false );
  //  virtual void readTableContents( type_t const* columnTypes, int numCols );

  virtual double getKeyValue( int indexLocation, int indexKey ) const;
  //  virtual double const* getKeyValues( int indexLocation ) const;

  virtual inline int numLocations() const { return myNumLocations; }
  virtual inline int numKeys() const { return myNumKeys; }
  virtual inline int numValues() const { return myNumValues; }
  virtual inline std::string tableName() const { return myTablename; }

  virtual std::string const keyName( int indexKey ) const;
  virtual std::string const valueName( int indexValue = 0 ) const;
  virtual inline int type() const { return myTableType; }

  virtual double getValue( int indexLocation, int indexValue ) const;
  virtual double getValue( double const* keyValues, int indexValue = 0 ) const;
  virtual TableValueList const* getValues( int indexLocation ) const;
  virtual TableValueList const* getValues( double const* keyValues ) const;

  virtual void dump() const;

  //---------------------------------------------------------
  // Time function
  virtual csTimeFunction<double> const* getFunction( double const* keyValues ) const;
  virtual double getTimeValue( double const* keyValues, double time ) const;

protected:
  /// Table name
  std::string myTablename;
  /// Table input file name
  std::string myFilename;
  /// Table input file stream
  std::FILE* myFile;
  /// Interpolation method
  int myMethodInterpolation;
  /// Type of table
  int myTableType;
  /// Data type of each column (for example TYPE_DOUBLE, TYPE_INT, TYPE_STRING...)
  type_t* myColumnTypes;

  /// Number of keys describing each location
  int myNumKeys;
  /// Names of keys
  std::string* myKeyNames;

  /// Number of values at each location
  int myNumValues;
  /// Names of values
  std::string* myValueNames;

  /// Number of locations
  int myNumLocations;

  /// Key values at each location (for 1D & 2D tables)
  double ** myKeyValues;

  // Index values, column locations in input file
  /// Number of columns in table
  int myNumCols;
  /// Index of first value column
  int myIndexFirstValueCol;
  bool myHasReadTableContents;
  bool myHasBeenInitialized;

  csVector<TableValueList*>* myValues;

 protected:
  csTableAll();
  csTableAll( csTableAll const& obj );
  bool findKeyLocation( double keyValue_in, int& locationIndex ) const;
  void findKeyLocation( double keyValue_in, int& locLeft, int& locRight, double& weight, int keyIndex = -1 ) const;
  void findKeyLocation2D( double const* keyValues_in, int& locLeftUp, int& locLeftDown, int& locRightUp, int& locRightDown ) const;
  double interpolate( int valueIndex, double const* keyValues_in ) const;
  void clearBuffers();
  void addData_internal( cseis_geolib::csVector<double>* valueList );

  //------------------------------------------------------------------------------------
  // Time function
  //
 protected:
  /// Values at each knee point of each location, including time values
  csTimeFunction<double>** myTimeFunctions2D;
  csTimeFunction<double>* myCurrentTimeFunction;
  /// Index of time column
  int myIndexTimeCol;
  double interpolateTimeFunction( double const* keyValues_in, double time ) const;
};

} // end namespace
#endif


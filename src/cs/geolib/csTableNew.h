/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TABLE_NEW_H
#define CS_TABLE_NEW_H

#include <cstdio>
#include <string>
#include "geolib_defines.h"
#include "csTableValueList.h"

namespace cseis_geolib {

  template<typename T> class csTimeFunction;
  template <typename T> class csVector;

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
* Format of input table:
* .. key1 ..  key2 .. keyN ..  [time] .. val1 ..  val2 .. valM
*  
* key?: Table keys, columns specified by user
* val?: Table values, columns specified by user
* time: Time key for tables containing a time axis (e.g. velocity profiles)
--------------------------------
*
* @author Bjorn Olofsson
*/

class csTableNew {
public:
  static const int TABLE_TYPE_UNKNOWN        = -1;
  static const int TABLE_TYPE_TIME_FUNCTION  = 1;
  static const int TABLE_TYPE_DUPLICATE_KEYS = 2;
  static const int TABLE_TYPE_UNIQUE_KEYS    = 3;

  static const char KEY_CHAR = '@';

public:
  csTableNew( int tableType );
  csTableNew( int tableType, int timeColumnIndex );
  virtual ~csTableNew();

  /**
   * Define table key 
   * @param column        Column index where key is located in table
   * @param doInterpolate true: Values in this table shall be interpolated according to this key. false: Do not interpolate values based to this key (key serves a pure look-up for table value)
   */
  void addKey( int column, bool doInterpolate );
  /**
   * Define table value
   * @param column        Column index where value is located in table
   */
  void addValue( int column );

  virtual void initialize( std::string const& filename, bool doSort = false );

  /**
   * @param indexLocation  Location index
   * @param indexKey       Key index (0 for first key, etc)
   * @return Key value
   */
  virtual double getKeyValue( int indexLocation, int indexKey ) const;

  /**
   * @return Number of 'locations' in table. Each location is defined by a unique set of key values
   */
  virtual inline int numLocations() const { return myNumLocations; }
  /**
   * @return Number of keys specified for this table 
   */
  virtual inline int numKeys() const { return myNumAllKeys; }
  /**
   * @return Number of values specified for this table 
   */
  virtual inline int numValues() const { return myNumValues; }

  virtual inline int type() const { return myTableType; }

  /**
   * @param indexLocation  Location index
   * @param indexValue     Value index (0 for first value, etc)
   * @return Value
   */
  virtual double getValue( int indexLocation, int indexValue ) const;
  /**
   * @param keyValues  Array containing key values
   * @param indexValue Value index (0 for first value, etc)
   * @return Value
   */
  virtual double getValue( double const* keyValues, int indexValue = 0 ) const;
  /**
   * @param indexLocation  Location index
   * @return List of values found at specified location
   */
  virtual csTableValueList const* getValues( int indexLocation ) const;
  /**
   * @param keyValues  Array containing key values
   * @return List of values found at specified location
   */
  virtual csTableValueList const* getValues( double const* keyValues ) const;

  virtual void dump() const;

  //---------------------------------------------------------
  // Time function
  virtual csTimeFunction<double> const* getFunction( double const* keyValues, bool dump = false ) const;
  virtual double getTimeValue( double const* keyValues, double time ) const;
  void interpolateTimeFunction( csTimeFunction<double> const* timeFuncLeft, csTimeFunction<double> const* timeFuncRight,
                                double weightLoc, csTimeFunction<double>* newTimeFunction ) const;

  /**
   * @return true if exact key was found
   */
  bool findInterpKeyLocation( double keyValue_in, int keyIndex, int& locLeft, int& locRight, double& weight ) const;
  void findInterpKeyLocation2D( double const* keyValues_in,
                                int locStart, int locEnd,
                                int& locLeftUp, int& locLeftDown,
                                int& locRightUp, int& locRightDown,
                                double& weight1, double& weight2 ) const;
  void findInterpKeyLocation2D( double const* keyValues_in,
                                int locStart, int locEnd,
                                int& locLeftUp, int& locLeftDown,
                                int& locRightUp, int& locRightDown,
                                double& weight1, double& weight2,
                                int keyIndex1, int keyIndex2 ) const;
  bool findExactKeyLocation( double const* keyValues_in, int& locStart, int& locEnd ) const;
  double interpolateStep2( double const* keyValues_in, int keyIndex1, int keyIndex2, int valueIndex,
                           int locLeftUp, int locLeftDown, int locRightUp, int locRightDown ) const;

protected:
  /// Table input file name
  std::string myFilename;
  /// Table input file stream
  std::FILE* myFile;
  /// Type of table
  int myTableType;
  /// Data type of each column (for example TYPE_DOUBLE, TYPE_INT, TYPE_STRING...)
  type_t* myColumnTypes;

  /// Number of keys describing each location
  int myNumKeys;
  int myNumAllKeys;
  int myNumInterpKeys;
  /// Column numbers containing keys
  int* myKeyCols;
  int* myKeyInterpCols;
  int* myKeyAllCols;
  /// Interpolate according to this key?
  //  bool* myKeyInterpolate;

  /// Number of values at each location
  int myNumValues;
  /// Column numbers containing values
  int* myValueColumns;

  /// Number of locations
  int myNumLocations;

  /// Key values at each location (for 1D & 2D tables)
  double ** myKeyValues;

  bool myHasReadTableContents;
  bool myHasBeenInitialized;

  csVector<csTableValueList*>* myValues;

 protected:
  csTableNew();
  csTableNew( csTableNew const& obj );
  void init( int tableType );
  void readTableContents( bool doSort = false );
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


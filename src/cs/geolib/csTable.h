/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_INPUT_TABLE_H
#define CS_INPUT_TABLE_H

#include <cstdio>
#include <string>
#include "geolib_defines.h"

namespace cseis_geolib {
  template<typename T> class csTimeFunction;
  template<typename T> class csVector;

/**
* Table of input time functions, specified at key locations
*
* This class shall provide functionality to...
*  - ...read in an ASCII file containing a table
*  - ...interpolate table values to any key location
*
* A table must be specified within the CSEIS flow file, using directive &table <tablename> <filename> {methodInterpolation}*
*
* Special key is 'time', which gives the vertical dimension
*
* NOTE: Implementation of this class is far from finished.
* Currently provided is
*  - reading in function at single key location
*
* Format of first line in input table:
* #name_key1  #name_key2 .. #name_keyN  [time]  name_val1  name_val2 .. name_valM
*  
* name_key?: Names of table keys
* name_val?: Names of table values
* time:      Optional time key for tables containing a time axis (e.g. velocity profiles)
*
* @author Bjorn Olofsson
* @date 2007
*/
class csTable {
public:
  static const int TABLE_DIM_UNKNOWN = -1;
  static const int TABLE_DIM_1D = 1;
  static const int TABLE_DIM_2D = 2;
  static const int TABLE_DIM_3D = 3;

  static const char KEY_CHAR = '#';

public:
  csTable();
  virtual ~csTable();

  virtual void initialize( std::string const& tablename, std::string const& filename, int methodInterpolation = 0 );
  void readTableContents( type_t const* types, int numCols );

  virtual void initialize_simpleTimeTable( int numLocations );
  virtual void setValue_simpleTimeTable( int locationIndex, csTimeFunction<double> const* timeFunction );

  virtual csTimeFunction<double> const* getFunction( double const* keyValues ) const;
  virtual csTimeFunction<double> const* getFunction( int indexLocation ) const;

  virtual double getValue( int indexLocation, int indexValue = 0 ) const;
  virtual double getValue( double const* keyValues, int indexValue = 0 ) const;
  virtual double getValue( int indexLocation, double time ) const;
  virtual double getValue( double const* keyValues, double time ) const;

  virtual double getKeyValue( int indexLocation, int indexKey ) const;

  virtual double const* getKeyValues( int indexLocation ) const;

  virtual inline int numLocations() const { return myNumLocations; }
  virtual inline int numKeys() const { return myNumKeys; }
  virtual inline std::string tableName() const { return myTablename; }

  virtual std::string const keyName( int indexKey ) const;
  virtual std::string const valueName( int indexValue = 0 ) const;
  virtual inline int dimension() const { return myTableDimension; }

private:
  /// Table name
  std::string myTablename;
  /// Table input file name
  std::string myFilename;
  /// Table input file stream
  FILE* myFile;  
  /// Interpolation method
  int myMethodInterpolation;
  /// Dimension of table
  int myTableDimension;
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

  /**
   * Values at each 1D location (for 1D tables)
   * There may be more than one value at one location, i.e. a list of values
   */
  double** myValues1D;
//  csVector<double>** myValues1D;
  /// Key values at each location (for 1D & 2D tables)
  double const** myKeyValues;
  /// Values at each knee point of each location, including time values
  csTimeFunction<double>** myTimeFunctions2D;
  csTimeFunction<double>* myCurrentTimeFunction;

  // Index values, column locations in input file
  /// Number of columns in table
  int myNumCols;
  /// Index of time column
  int myIndexTimeCol;
  /// Index of first value column
  int myIndexFirstValueCol;
  bool myHasReadTableContents;
  bool myHasBeenInitialized;

private:
  void clearBuffers();
  double interpolate2D( double const* keyValues_in, double time ) const;
  double interpolate1D( int indexValue, double const* keyValues_in ) const;
  void findKeyLocation( double const* keyValues_in, int& locLeft, int& locRight, double& weight ) const;
};

} // end namespace
#endif


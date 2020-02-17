/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TABLE_VALUE_LIST_H
#define CS_TABLE_VALUE_LIST_H

namespace cseis_geolib {

class csTableValueList {
 public:
  csTableValueList( int numCols, int numLines ) {
    myNumLines   = numLines;
    myNumCols    = numCols;
    myValues       = new double[myNumCols*myNumLines];
  }
  ~csTableValueList() {
    if( myValues != NULL ) {
      delete [] myValues;
      myValues = NULL;
    }
  }
  void set( int indexCol, int indexLine, double value ) {
    myValues[indexCol*myNumLines+indexLine] = value;
  }
  double get( int indexCol, int indexLine = 0 ) const {
    return myValues[indexCol*myNumLines+indexLine];
  }
  double const* getAll( int indexCol ) const {
    return &myValues[indexCol*myNumLines];
  }
  int numColumns() const { return myNumCols; }
  int numLines() const { return myNumLines; }
 private:
  double* myValues;
  int myNumCols;
  int myNumLines;
};

} // end namespace
#endif


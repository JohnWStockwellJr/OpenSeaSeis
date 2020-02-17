/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_TIME_FUNCTION_H
#define CS_TIME_FUNCTION_H

#include "csException.h"
#include "csVector.h"
#include "cstring"

namespace cseis_geolib {

/**
 * Time function, used in Table
 *
 * @author Bjorn Olofsson
 * @date 2008
 */
template <typename T>
class csTimeFunction {
public:
  csTimeFunction();
  csTimeFunction( int numSpatialValues );
  ~csTimeFunction();
  inline void resetNumSpatialValues( int numSpatialValues );
  void set( csTimeFunction<T> const* timeFunc );
  void set( T const* values, double const* times, int numValues ) { set( values, times, numValues, 0 ); }
  void set( T const* values, double const* times, int numValues, int spatialIndex );

  void set( csVector<T> const* valueList, csVector<double> const* timeList ) { set( valueList, timeList, 0 ); }
  void set( csVector<T> const* valueList, csVector<double> const* timeList, int spatialIndex );

  inline T valueAt( double time ) const { return valueAt( time, 0 ); }
  inline T valueAt( double time, int spatialIndex ) const;
  inline T valueAtIndex( int index ) const { return valueAtIndex( index, 0 ); }
  inline T valueAtIndex( int index, int spatialIndex ) const;

  inline double timeAtIndex( int timeIndex ) const;
  inline int numValues() const;
  inline int numSpatialValues() const;
  inline void applyPercentageChange( float percentage ) { applyPercentageChange(percentage,0); }
  inline void applyPercentageChange( float percentage, int spatialIndex );

  inline void dump() const;

private:
  void resetBuffer( int numValues );
  float mySampleInt;
  double* myTimeBuffer;
  int myNumValues;
  int myNumSpatialValues;
  cseis_geolib::csVector<T*>* myValueList;
};

template<typename T>csTimeFunction<T>::csTimeFunction() {
  myTimeBuffer  = NULL;
  myNumValues   = 0;
  myNumSpatialValues = 1;
  myValueList = new csVector<T*>();
  myValueList->insertEnd( NULL );
}

template<typename T>csTimeFunction<T>::csTimeFunction( int numSpatialValues ) {
  myTimeBuffer  = NULL;
  myNumValues   = 0;
  myNumSpatialValues = numSpatialValues;
  myValueList = new csVector<T*>();
  for( int i = 0; i < numSpatialValues; i++ ) {
    myValueList->insertEnd( NULL );
  }
}
template<typename T>csTimeFunction<T>::~csTimeFunction() {
  if( myTimeBuffer != NULL ) {
    delete [] myTimeBuffer;
    myTimeBuffer = NULL;
  }
  if( myValueList != NULL ) {
    for( int i = 0; i < myValueList->size(); i++ ) {
      if( myValueList->at(i) != NULL ) {
        delete [] myValueList->at(i);
      }
    }
    delete myValueList;
    myValueList = NULL;
  }
}

template<typename T> inline double csTimeFunction<T>::timeAtIndex( int timeIndex ) const {
  return myTimeBuffer[timeIndex];
}

template<typename T> void csTimeFunction<T>::set( csTimeFunction<T> const* timeFunc ) {
  resetNumSpatialValues( timeFunc->numSpatialValues() );
  resetBuffer( timeFunc->numValues() );
  std::memcpy( myTimeBuffer, timeFunc->myTimeBuffer, sizeof(T)*myNumValues );
  for( int i = 0; i < myNumSpatialValues; i++ ) {
    T const* valuesIn = timeFunc->myValueList->at(i);
    T* valuesOut = myValueList->at(i);
    std::memcpy( valuesOut, valuesIn, sizeof(T)*myNumValues );
    //    set( timeFunc->myValueList->at(i), timeFunc->myTimeBuffer, timeFunc->myNumValues, i );
  }
}

template<typename T> void csTimeFunction<T>::set( csVector<T> const* valueList, csVector<double> const* timeList, int spatialIndex ) {
  int numValues = valueList->size();
  resetBuffer( numValues );
  T* valuesOut = myValueList->at(spatialIndex);
  for( int i = 0; i < myNumValues; i++ ) {
    myTimeBuffer[i]  = timeList->at(i);
    valuesOut[i]     = valueList->at(i);
  }
}
 template<typename T> void csTimeFunction<T>::set( T const* valuesIn, double const* timesIn, int numValues, int spatialIndex ) {
  resetBuffer( numValues );
  T* valuesOut = myValueList->at(spatialIndex);
  for( int i = 0; i < myNumValues; i++ ) {
    myTimeBuffer[i]  = timesIn[i];
    valuesOut[i]     = valuesIn[i];
  }
}

 template<typename T> inline void csTimeFunction<T>::applyPercentageChange( float percentage, int spatialIndex ) {
  float scalar = (percentage+100) / 100.0;
  T* valuesOut = myValueList->at(spatialIndex);
  for( int i = 0; i < myNumValues; i++ ) {
    valuesOut[i] *= scalar;
  }
}

template<typename T> inline int csTimeFunction<T>::numValues() const {
  return myNumValues;
}
template<typename T> inline int csTimeFunction<T>::numSpatialValues() const {
  return myNumSpatialValues;
}
template<typename T> inline T csTimeFunction<T>::valueAtIndex( int index, int spatialIndex ) const {
  return myValueList->at(spatialIndex)[index];
}
template<typename T> inline T csTimeFunction<T>::valueAt( double time, int spatialIndex ) const {
  T const* values = myValueList->at(spatialIndex);
  if( myNumValues == 0 ) {
    throw( csException("Called csTimeFunction 'value' functon, but object is not set yet") );
  }
  else if( time <= myTimeBuffer[0] ) {
    return values[0];
  }
  else if( time >= myTimeBuffer[myNumValues-1] ) {
    return values[myNumValues-1];
  }
  // Linear interpolation
  else {
    int topIndex = -1;
    // Temp, slow 'find' method for correct index... !CHANGE!
    for( int i = 1; i < myNumValues; i++ ) {
      if( myTimeBuffer[i] >= time ) {
        topIndex = i-1;
        break;
      }
    }
    int bottomIndex = topIndex+1;
    double weight = (time-myTimeBuffer[topIndex]) / (myTimeBuffer[bottomIndex]-myTimeBuffer[topIndex]);
    return( (T)( (double)values[topIndex] + weight*( values[bottomIndex] - values[topIndex] ) ) );
  }
}

template<typename T> void csTimeFunction<T>::resetNumSpatialValues( int numSpatialValues ) {
  for( int i = myValueList->size(); i < numSpatialValues; i++ ) {
    myValueList->insertEnd( NULL );
    //myValueList->insertEnd( new T[myNumValues] );
  }
  myNumSpatialValues = numSpatialValues;
}
template<typename T> void csTimeFunction<T>::resetBuffer( int numValues ) {
  if( numValues > myNumValues ) {
    for( int i = 0; i < myNumSpatialValues; i++ ) {
      T* values = myValueList->at(i);
      if( values != NULL ) delete [] values;
      myValueList->set( new T[numValues], i );
    }
    if( myTimeBuffer != NULL ) {
      delete [] myTimeBuffer;
      myTimeBuffer = NULL;
    }
    myTimeBuffer = new double[numValues];
  }
  myNumValues = numValues;
}
template<typename T> void csTimeFunction<T>::dump() const {
  for( int index = 0; index < myNumSpatialValues; index++ ) {
    for( int itime = 0; itime < myNumValues; itime++ ) {
      fprintf(stdout,"%d %d   %f %f\n", index, itime, myTimeBuffer[itime], myValueList->at(index)[itime]);
    }
  }
}

} // namespace
#endif


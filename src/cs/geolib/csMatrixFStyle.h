/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_MATRIX_FSTYLE_H
#define CS_MATRIX_FSTYLE_H

#include <iostream>

namespace cseis_geolib {

/**
 * csMatrix class, FORTRAN style.
 *
 * Up to 4 dimensions.
 * If all matrix elements are accessed in a double loop, the outer loop should
 * be for the second dimension, and the inner loop for the first dimension.
 *
 * Example:
 *   for( int j = 0; j < myDim2; ++j ) {
 *      for( int i = 0; i < myDim1; ++i ) {
 *         matrix(i,j) = .....;
 *      }
 *   }
 *
 * @author Bjorn Olofsson
 * @date 2005
 */
template <typename T>
class csMatrixFStyle
{
public:
    /**
     * Constructor
     * @param d1 : size in first dimension
     * @param d2 : size in second dimension
     */
    csMatrixFStyle( int d1, int d2 );
    /**
     * Constructor
     * @param d1 : size in first dimension
     * @param d2 : size in second dimension
     * @param d3 : size in third dimension
     */
    csMatrixFStyle( int d1, int d2, int d3 );
    /**
     * Constructor
     * @param d1 : size in first dimension
     * @param d2 : size in second dimension
     * @param d3 : size in third dimension
     * @param d4 : size in fourth dimension
     */
    csMatrixFStyle( int d1, int d2, int d3, int d4 );

    /// Default constructor.
    csMatrixFStyle();

    /**
     * Copy constructor
     * @param mat : matrix to copy from
     */
    csMatrixFStyle( const csMatrixFStyle<T>& mat );
  
    /// Destructor.
    ~csMatrixFStyle();

    /** Access operator
     * @param d1 : index in first dimension
     * @param d2 : index in second dimension
     * @param d3 : index in third dimension
     * @param d4 : index in fourth dimension
     * @return element (d1,d2,d3,d4) by reference
     */
    inline T& operator()(int d1, int d2=0, int d3=0, int d4=0);
  
    /** Access operator
     * @param d1 : index in first dimension
     * @param d2 : index in second dimension
     * @param d3 : index in third dimension
     * @param d4 : index in fourth dimension
     * @return element (d1,d2,d3,d4) by reference
     */
    inline const T& operator()( int d1, int d2=0, int d3=0, int d4=0 ) const;
  
    /** Value access function
     * @param d1 : index in first dimension
     * @param d2 : index in second dimension
     * @param d3 : index in third dimension
     * @param d4 : index in fourth dimension
     * @return element (d1,d2,d3,d4) by value
     */
    inline T getValue( int d1, int d2=0, int d3=0, int d4=0 ) const;

    /**
     * Set matrix to single value
     */
    inline void set( T value );
  
    /** Standard Assignement operator.
     * @param mat : matrix to assign from
     * @return reference to current instance
     */    
    inline csMatrixFStyle<T>& operator=  (const csMatrixFStyle<T>& mat);

    /** Assignement operator.
     * @param aScalar : scalar to assign from
     * @return reference to current instance
     */    
    inline csMatrixFStyle<T>& operator=  (const T& aScalar); 

    /** Assignement and add  operator.
     * @param mat : matrix to add
     * @return reference to current instance
     */
    inline csMatrixFStyle<T>& operator+= (const csMatrixFStyle<T>& mat);

    /** Assignement and add operator.
     * @param aScalar : scalar to add
     * @return reference to current instance
     */
    inline csMatrixFStyle<T>& operator+= (const T& aScalar); 

    /** Assignement and substract operator.
     * @param mat : matrix to substract
     * @return reference to current instance
     */
    inline csMatrixFStyle<T>& operator-= (const csMatrixFStyle<T>& mat);

    /** Assignement and substract operator.
     * @param aScalar : scalar to substract
     * @return reference to current instance
     */
    inline csMatrixFStyle<T>& operator-= (const T& aScalar); 

    /** Assignement and divide operator.
     * @param aScalar : scalar to divide by
     * @return reference to current instance
     */
    inline csMatrixFStyle<T>& operator/= (const T& aScalar); 

    /** Assignement and multiply operator.
     * @param aScalar : scalar to multiply by
     * @return reference to current instance
     */
    inline csMatrixFStyle<T>& operator*= (const T& aScalar); 
    
    /**
     * @return size of n'th dimension, n=1,2,3,4
     */
    inline int getDim( int n ) const;
    
    /** Pointer getter
     * @param d1 : index in first dimension
     * @param d2 : index in second dimension
     * @param d3 : index in third dimension
     * @param d4 : index in fourth dimension
     * @return pointer to element (d1,d2,d3,d4)
     */
    inline T* getPointer(int d1=0, int d2=0, int d3=0, int d4=0) const;

    /**
     * Resize matrix
     * @param s1 : size in first dimension
     * @param s2 : size in second dimension
     * @param s3 : size in second dimension
     * @param s4 : size in second dimension
     * @note Previous content of matrix is lost
     */
    inline void resize(int s1, int s2, int s3=1, int s4=1);

    /** Minimum function
     * For each element the minimum of the current instance
     * and the parameter matrix is taken.
     * @param mat : matrix to compare with
     */
    inline void min(const csMatrixFStyle<T>& mat);

    /** Maximum function
     * For each element the maximum the current instance
     * and the parameter matrix is taken.
     * @param mat : matrix to compare with
     */
    inline void max(const csMatrixFStyle<T>& mat);

    /** Element by element product
     * @param mat : matrix to multiply by
     */
    inline void dotProduct(const csMatrixFStyle<T>& mat);

    /** Non null function
     * Elements are set to 1 for non null elements
     * of a matrix, 0 otherwise.
     * @param mat : matrix to examine
     */
    inline void nonNull(const csMatrixFStyle<T>& mat);

    /// Shift the 4 quadrants
    void fftShift();

    /// Transpose the matrix
    void transpose();
 
    /// Formatted print
    inline void print() const;
 
    /// Dump elements as a column vector
    inline void dump() const;

    /** Check if same dimension
     * @return true if same dimension else otherwise
     */
    inline bool sameDim(const csMatrixFStyle<T>& mat) const;

    /** Check if same size
     * @return true if product of dimensions is the same
     */
    inline bool sameSize(const csMatrixFStyle<T>& mat) const;

    /** Total size
     * @return product of first and second dimension
     */
    int inline size() const;

private:
    void init( int d1, int d2, int d3, int d4 );

    /// Number of dimensions
    int myNumDimensions;
    /// First dimension
    int myDim1;
    /// Second dimension
    int myDim2;
    /// Third dimension
    int myDim3;
    /// Fourth dimension
    int myDim4;
    /// Total Size
    int mySize;
    /// Pointer to data
    T* myData;  
};

 template <typename T> csMatrixFStyle<T>::csMatrixFStyle(int d1, int d2) {
   myNumDimensions = 2;
   init(d1,d2,1,1);
 }

 template <typename T> csMatrixFStyle<T>::csMatrixFStyle(int d1, int d2, int d3) {
   myNumDimensions = 3;
   init(d1,d2,d3,1);
 }

 template <typename T> csMatrixFStyle<T>::csMatrixFStyle(int d1, int d2, int d3, int d4) {
   myNumDimensions = 4;
   init(d1,d2,d3,d4);
 }

 template <typename T> void csMatrixFStyle<T>::init(int d1, int d2, int d3, int d4) {
   myDim1 = d1;
   myDim2 = d2;
   myDim3 = d3;
   myDim4 = d4;
   mySize = myDim1 * myDim2 * myDim3 * myDim4;
   myData = new T[mySize];
 }

 template <typename T> csMatrixFStyle<T>::csMatrixFStyle() :
 myDim1(0),  myDim2(0), myDim3(0), myDim4(0), mySize(0), myData(0)
   {}

 template <typename T> csMatrixFStyle<T>::csMatrixFStyle( const csMatrixFStyle<T>& mat ) {
   myDim1 = mat.getFirstDim();
   myDim2 = mat.getSecondDim();
   mySize = mat.size();
   myData = new T[mySize];
   for( int i = 0; i < mySize; i++ ) {
     myData[i] = mat(i);
   }
 }

 template <typename T> csMatrixFStyle<T>::~csMatrixFStyle() {
   delete[] myData;
 }

 /*
s1=3 s2=2
  1 1  <-- a(1,1)
  2 1
  3 1
  1 2
  2 2
  3 2  <-- a(3,2) = b[ (d1-1) + (d2-1)*s1 ] =  b[ 2 + 1*s1 ]

s1=3 s2=2 s3=4
b[ d1 + d2*s1 + d3*s1*s2 ]
b[ d1 + d2*s1 + d3*s1*s2 + d4*s1*s2*s3 ]
b[ ((d4*s3 + d3)*s2 + d2)*s1 + d1 ]

  1 1 1  <-- a(1,1,1)
  2 1 1
  3 1 1
  1 2 1
  2 2 1
  3 2 1  <-- a(3,2,1) = b[ (d1-1) + (d2-1)*s1 + (d3-1)*s1*s2 ] =  b[ 2 + 1*s1 + 0*s1*s2 ] = b[5]

  1 1 2  <-- a(1,1,2) = b[ (d1-1) + (d2-1)*s1 + (d3-1)*s1*s2 ] =  b[ 0 + 0*s1 + 1*s1*s2 ] = b[6]
  2 1 2
  3 1 2
  1 2 2
  2 2 2
  3 2 2

  1 1 3
  2 1 3
  3 1 3
  1 2 3
  2 2 3
  3 2 3

  1 1 4  <-- a(1,1,2) = b[ (d1-1) + (d2-1)*s1 + (d3-1)*s1*s2 ] =  b[ 0 + 0*s1 + 3*s1*s2 ] = b[18]
  2 1 4
  3 1 4
  1 2 4
  2 2 4
  3 2 4  <-- a(1,1,2) = b[ (d1-1) + (d2-1)*s1 + (d3-1)*s1*s2 ] =  b[ 2 + 1*s1 + 3*s1*s2 ] = b[23]

  */

 template <typename T> inline const T& csMatrixFStyle<T>::operator()(int d1, int d2, int d3, int d4) const {
   return myData[ ((d4*myDim3 + d3)*myDim2 + d2)*myDim1 + d1 ];
 }
 template <typename T> inline T& csMatrixFStyle<T>::operator()(int d1, int d2, int d3, int d4) {
   //   fprintf(stderr,"Retrieve element #%-3d  (%d)\n", (((d4*myDim3 + d3)*myDim2 + d2)*myDim1 + d1), mySize );
   return myData[ ((d4*myDim3 + d3)*myDim2 + d2)*myDim1 + d1 ];
 }

 template <typename T> inline T csMatrixFStyle<T>::getValue(int d1, int d2, int d3, int d4) const {
   return myData[ ((d4*myDim3 + d3)*myDim2 + d2)*myDim1 + d1 ];
 }

 template <typename T> inline void csMatrixFStyle<T>::set( T value ) {
   memset(myData,int(value),mySize*sizeof(T));
 }

 template <typename T> inline int csMatrixFStyle<T>::getDim( int n ) const {
   switch( n ) {
   case 1:
     return myDim1;
   case 2:
     return myDim2;
   case 3:
     return myDim3;
   case 4:
     return myDim4;
   }
   return myDim1;
 }

 template <typename T> inline T* csMatrixFStyle<T>::getPointer(int d1, int d2, int d3, int d4) const {
   return( myData + ((d4*myDim3 + d3)*myDim2 + d2)*myDim1 + d1 );
 }

 template <typename T> inline csMatrixFStyle<T>& csMatrixFStyle<T>::operator=( const csMatrixFStyle<T>& mat ) {
   if( &mat != this ) {
     delete[] myData;
     init( mat.getDim(1), mat.getDim(2), mat.getDim(3), mat.getDim(4) );
     myNumDimensions = mat.myNumDimenstions;
     for( int i = 0; mySize < 0; i++ ) {
       myData[i] = mat(i);
     } 
   }
   return *this;
 }

 template <typename T> inline csMatrixFStyle<T>& csMatrixFStyle<T>::operator=(const T& value) {
   for(int i = 0; i < mySize; i++ ) {
     myData[i] = value;
   }
   return *this;
 }

 template <typename T> inline csMatrixFStyle<T>& csMatrixFStyle<T>::operator+=( const csMatrixFStyle<T>& mat ) {
   for(int i = 0; i < mySize; i++ ) {
     myData[i] += mat(i);
   }
   return *this;
 }

 template <typename T> inline csMatrixFStyle<T>& csMatrixFStyle<T>::operator+=( const T& value ) {
   for(int i = 0; i < mySize; i++ ) {
     myData[i] += value;
   }
   return *this;
 }

 template <typename T> inline csMatrixFStyle<T>& csMatrixFStyle<T>::operator-=( const csMatrixFStyle<T>& mat ) {
   for(int i = 0; i < mySize; i++ ) {
     myData[i] -= mat(i);
   }
   return *this;
 }

 template <typename T> inline csMatrixFStyle<T>& csMatrixFStyle<T>::operator-=(const T& value) {
   for(int i = 0; i < mySize; i++ ) {
     myData[i] -= value;
   }
   return *this;
 }

 template <typename T> inline csMatrixFStyle<T>& csMatrixFStyle<T>::operator/=( const T& value ) {
   if( value !=0 ) {
     T val2 = 1.0/value;
     for(int i = 0; i < mySize; i++ ) {
       myData[i] *= val2;
     }
   }
   return *this; 
 }

 template <typename T> inline csMatrixFStyle<T>& csMatrixFStyle<T>::operator*=( const T& value ) {
   for(int i = 0; i < mySize; i++ ) {
     myData[i] *= value;
   }
   return *this;
 }

 template <typename T> inline void csMatrixFStyle<T>::resize(int d1, int d2, int d3, int d4) {
   if( myData ) delete [] myData;
   init(d1,d2,d3,d4);
   if( d3 == 1 && d4 == 1 ) {
     myNumDimensions = 2;
   }
   else if( d4 == 1 ) {
     myNumDimensions = 3;
   }
   else {
     myNumDimensions = 4;
   }
 }

 template <typename T> inline void csMatrixFStyle<T>::min( const csMatrixFStyle<T>& mat ) {
   for(int i = 0; i < mySize; i++ ) {
     if( mat(i) < myData[i] ) {
       myData[i] = mat(i);
     }
   }
 }

 template <typename T> inline void csMatrixFStyle<T>::max( const csMatrixFStyle<T>& mat ) {
   for( int i = 0; i < mySize; i++ ) {
     if( mat(i) > myData[i] ) { 
       myData[i] = mat(i);
     }
   }
 }

 template <typename T> inline void csMatrixFStyle<T>::dotProduct( const csMatrixFStyle<T>& mat ) {
   for( int i = 0; i < mySize; i++ ) {
     myData[i] *= mat(i);
   }
 }

 template <typename T> void csMatrixFStyle<T>::transpose() {
   T* tmp = new T[mySize];
   for(int j=myDim2-1;j>=0;--j) {
     for(int i=myDim1-1;i>=0;--i)
       { tmp[j + i*myDim2]= myData[i + j*myDim1];} }
         
   int tmpDim = myDim2; myDim2 = myDim1;  myDim1 =tmpDim;        
   delete[] myData;
   myData=tmp;
 }

 template <typename T> inline void csMatrixFStyle<T>::print() const {
   for( int i1 = 0; i1 < myDim1; i1++ ) {
     for( int i2 = 0; i2 < myDim2; i2++ ) {
       for( int i3 = 0; i3 < myDim3; i3++ ) {
         for( int i4 = 0; i4 < myDim3; i4++ ) {
           std::cout << this(i1,i2,i3,i4) << "  ";
         }
         std::cout << std::endl;
       }
     }
   }
 }

 template <typename T> inline void csMatrixFStyle<T>::dump() const {
   for(int i=0; i < mySize; i++ ) { 
     std::cout <<  myData[i] <<  std::endl;
   }
 }

 template <typename T> inline void csMatrixFStyle<T>::nonNull( const csMatrixFStyle<T>& mat ) {
   for( int i = 0; i < mySize; i++ ) {
     myData[i] = (mat(i)!=0) ? 1 : 0;
   }
 }

 template <typename T> void csMatrixFStyle<T>::fftShift() {
   int d1 = myDim1/2, d2 = myDim2/2;
   T* tmp = new T[mySize];
   for(int i=d1-1; i>=0; --i) {
     for(int j=d2-1; j>=0; --j)
       {
         tmp[i + j*myDim1] = myData[i+d1 + (j+d2)*myDim1];
         tmp[i+d1 + (j+d2)*myDim1] = myData[i + j*myDim1];
         tmp[i+d1 + j*myDim1] = myData[i + (j+d2)*myDim1];
         tmp[i + (j+d2)*myDim1] = myData[i+d1 + j*myDim1];
       } }
   
   delete[] myData;
   myData=tmp;
 }

 template <typename T> inline bool  csMatrixFStyle<T>::sameDim( const csMatrixFStyle<T>& mat ) const {
   return ( (mat.getDim(1)==myDim1) &&  
            (mat.getDim(2)==myDim2) && 
            (mat.getDim(3)==myDim3) && 
            (mat.getDim(4)==myDim4));
 }

 template <typename T> inline int csMatrixFStyle<T>::size() const {
   return(mySize);
 }

 template <typename T> inline bool csMatrixFStyle<T>::sameSize(const csMatrixFStyle<T>& mat) const {
   return ( mat.size() == mySize );
 }

} // END namespace

#endif


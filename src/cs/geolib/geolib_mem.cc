/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "geolib_mem.h"
#include <cstring>

namespace cseis_geolib {

//-----------------------------------------------------------------
// Memory allocation
//
void deleteString( char* str ) {
  if( str ) { delete [] str; str = NULL; }
}

//-----------------------------------------------------------------
char* createString( char const* strIn ) {
  char* strOut;
  int length = strlen( strIn );
  if( length <= 0 ) {
    return NULL;
  }
  strOut = new char[length+1];
  memcpy( strOut, strIn, length );
  strOut[length] = '\0';
  return strOut;
}

//-----------------------------------------------------------------
char** allocate_string_list( int numStrings, int strLength ) {
  int i;
  //  char** strList = (char **) malloc((size_t)((numStrings)*sizeof(char*)));
  char** strList = new char* [numStrings];
  if (!strList) {
    return(0);
  }
  //  strList[0] = (char *) malloc((size_t)((numStrings*strLength)*sizeof(char)));
  strList[0] = new char[numStrings*strLength];
  if (!strList[0]) {
    return(0);
  }
  for( i = 1; i < numStrings; i++ ) {
    strList[i] = strList[i-1] + strLength;
  }
  return strList;
}

//-----------------------------------------------------------------
float** allocate_matrix( int nrows, int ncols ) {
  int i;
  float** matrix = new float*[nrows];
  //float** matrix = (float **) malloc((size_t)((nrows)*sizeof(float*)));
  if (!matrix) {
    return(0);
  }
  //matrix[0]=(float *) malloc((size_t)((nrows*ncols)*sizeof(float)));
  matrix[0] = new float[nrows*ncols];
  if (!matrix[0]) {
    return(0);
  }
  for( i = 1; i < nrows; i++ ) {
    matrix[i] = matrix[i-1] + ncols;
  }
  return matrix;
}

//-----------------------------------------------------------------
float* allocate_vector( int nrows ) {
  float* vec;
  //  vec = (float *)malloc((size_t) ((nrows)*sizeof(float)));
  vec = new float [nrows];
  if (!vec) {
    return 0;
  }
  return vec;
}

//-----------------------------------------------------------------
int* allocate_vector_int( int nrows ) {
  int* vec;
  //  vec = (int *)malloc((size_t) ((nrows)*sizeof(int)));
  vec = new int[nrows];
  if (!vec) {
    return 0;
  }
  return vec;
}

//-----------------------------------------------------------------
void free_vector( float*& vec ) {
  if( vec ) {
    delete [] vec;
    vec = NULL;
  }
}

//-----------------------------------------------------------------
void free_vector_int( int*& vec ) {
  if( vec ) {
    delete [] vec;
    vec = NULL;
  }
}

//-----------------------------------------------------------------
void free_matrix( float**& matrix ) {
  if( matrix ) {
    delete [] matrix[0];
    delete [] matrix;
    matrix = NULL;
  }
}

//-----------------------------------------------------------------
void free_string_list( char** strList ) {
  if( strList ) {
    delete [] strList[0];
    delete [] strList;
    strList = NULL;
  }
}

} // namespace


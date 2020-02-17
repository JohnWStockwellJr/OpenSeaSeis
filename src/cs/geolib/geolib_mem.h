/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

/*********************************************************************************
 *
 * Memory 'library' routines
 *
 */

#ifndef GEOLIB_MEM_H
#define GEOLIB_MEM_H

namespace cseis_geolib {

  char** allocate_string_list( int numStrings, int strLength );
  float** allocate_matrix( int nrows, int ncols );
  float* allocate_vector( int nrows );
  int* allocate_vector_int( int nrows );

  void free_string_list( char**& strList );
  void free_matrix( float**& matrix );
  void free_vector( float*& vec );
  void free_vector_int( int*& vec );

  char* createString( char const* strIn );
  void deleteString( char* str );
} // namespace

#endif


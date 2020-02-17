/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_RSF_HEADER
#define CS_RSF_HEADER

#include <cstdio>
#include <string>
#include <cstdlib>

using namespace std;

namespace cseis_io {

class csRSFHeader {
 public:
  static const int DATA_FORMAT_UNKNOWN = 0;
  static const int DATA_FORMAT_FLOAT = 1;

  // Don't change index of these fields:
  static int const HDR_NONE  = -1;
  static int const HDR_DIM2  = 0;
  static int const HDR_DIM3  = 1;
  static int const HDR_DIM4  = 2;
  static int const HDR_DIM5  = 3;
  static int const HDR_TRCNO = 4;


  /// Sample data unit = Unit of vertical axis = Unit of sample interval, specified in superheader
  static int const SAMPLE_UNIT_UNKNOWN = -1;
  static int const SAMPLE_UNIT_MS = 401;
  static int const SAMPLE_UNIT_M  = 402;
  static int const SAMPLE_UNIT_HZ = 403;
  static int const SAMPLE_UNIT_S  = 404;
  static int const SAMPLE_UNIT_KM = 405;

  csRSFHeader();
  void set( csRSFHeader const& hdr );
  ~csRSFHeader();
  bool setField( char const* name, char const* value );
  void dump();
  void dump( FILE* stream );
  void dump( FILE* stream, bool outputGrid );
  static bool readRSFHeaderFile( char const* filename_rsf, csRSFHeader& hdr );

  void setFilename( char const* theFilename_bin, char const* theFilename_bin_full_path );
  static int readSampleUnit( std::string text );
public:
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;

  double o1;
  double o2;
  double o3;
  double o4;
  double o5;

  double d1;
  double d2;
  double d3;
  double d4;
  double d5;

  double e1;
  double e2;
  double e3;
  double e4;
  double e5;

  double world_x1;
  double world_x2;
  double world_x3;

  double world_y1;
  double world_y2;
  double world_y3;

  double il1;
  double il2;
  double il3;

  double xl1;
  double xl2;
  double xl3;

  double ild;
  double xld;

  double sou_x;
  double sou_y;
  double sou_z;

  int unit1;
  int unit2;
  int unit3;
  int unit4;
  int unit5;

  int domain1;

  int data_format;
  int esize;
  std::string filename_bin;
  std::string filename_bin_full_path;

 private:
  csRSFHeader( csRSFHeader const& hdr );
};

} // end namespace

#endif


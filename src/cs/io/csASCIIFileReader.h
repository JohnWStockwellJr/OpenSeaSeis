/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_ASCII_FILE_READER_H
#define CS_ASCII_FILE_READER_H

#include <cstdio>

namespace cseis_geolib {
  template <typename T> class csVector;
}

namespace cseis_io {
  class csASCIIFileReader;

class ASCIIParam {
 public:
  ASCIIParam();
  ~ASCIIParam();
  float sample(int index) const;
  double time(int index) const;
  int numSamples() const;
  float const* getSamples() const;
  double const* getTimes() const;
  void clear();

  float sampleInt;
  float srcDepth;
  double timeFirstSamp;
  double timeLastSamp;
  int traceNumber;
  friend class csASCIIFileReader;
 private:
  int myNumSamples;
  cseis_geolib::csVector<float>* sampleList;
  cseis_geolib::csVector<double>* timeList;
};

class csASCIIFileReader {
 public:
  static int const FORMAT_NUCLEUS_SIGNATURE = 11;
  static int const FORMAT_COLUMNS  = 12;
  static int const FORMAT_ZMAP     = 13;
  static int const FORMAT_NUCLEUS_PLUS = 14;
  static int const FORMAT_SPIKOGRAM  = 15;

  /**
   * @param filename  ASCII file name.
   * @param format    ASCII file format, see constants above.
   */
  csASCIIFileReader( std::string const& filename, int format );
  ~csASCIIFileReader();
  /**
   * Read one trace from input file.
   * @param param  Contains all extracted data.
   */
  bool readNextTrace( ASCIIParam* param );
  bool initialize( ASCIIParam* param, int colIndexTime, int colIndexValue, int colIndexTrace );
  bool initialize( ASCIIParam* param );
  bool initializeZMap( ASCIIParam* param,
                       int& zmap_numTraces,
                       double& zmap_x1,
                       double& zmap_y1,
                       double& zmap_x2,
                       double& zmap_y2 );
  bool isAtEOF() const;

 private:
  //  bool readOneTraceColumnFormat( ASCIIParam* param );
  bool readOneTraceColumnFormat( cseis_geolib::csVector<double>* timeList,
                                 cseis_geolib::csVector<float>* sampleList,
                                 int maxSamplesToRead,
                                 int* newTraceNumber,
                                 bool bailOut = true );
  bool readOneTraceSpikogramFormat( cseis_geolib::csVector<double>* timeList,
                                    cseis_geolib::csVector<float>* sampleList );

  FILE* myFileASCII;
  int myFormat;
  char myBuffer[1024];
  int myNumColumns;
  int myCounterLines;
  int myCurrentTraceIndex;
  bool myIsAtEOF;
  int myNumSamples;
  int myNumMinColumns;
  int myColIndexTime;
  int myColIndexValue;
  int myColIndexTrace;

  bool myIsOneSampleWaiting;
  float mySampleValueWaiting;
  float myTimeValueWaiting;
  int myTraceIndexWaiting;

  double myZMap_noValue;
};


} // end namespace

#endif


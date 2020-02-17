/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_P190_READER_H
#define CS_P190_READER_H

#include <cstdio>
#include <string>
#include <fstream>

namespace cseis_geolib {
  template <typename T> class csVector;
}

namespace cseis_io {

class csDataSource {
 public:
  int id;
  int point;
  double wdep;
  int day;
  int hour;
  int min;
  int sec;
  double x;
  double y;
  double z;
  int lineNumber;
};

class csDataChan {
 public:
  int chan;
  int cable;
  double x;
  double y;
  double z;
};


class csP190Reader {
 public:
  csP190Reader( std::string const& filename );
  ~csP190Reader();

  void initialize();
  csDataSource const* getSource( int source );
  csDataChan const* getChan( int source, int chan, int cable );

  void dump() const;
  void dumpAllSourceInfo() const;
 private:
  bool readAllSources();
  bool readChannels( int source );
  int scanChan( std::string& line, int index );
  void scanSource( std::string& line, csDataSource* data );
  void setFilePointerToStartOfLine( int lineNum );
  int getSourceIndex( int source, int startIndex ) const;

  std::ifstream* myFile;
  //  FILE* myFile;
  char myBuffer[1024];
  int myCurrentSource;
  int myCurrentSourceIndex;
  int myCurrentChanIndex;
  int myCurrentLineNumber; // Pointer to current line in input file
  int myCurrentMinChan;
  int myCurrentMaxChan;
  //  csDataSource myCurrentSourceData;
  cseis_geolib::csVector<csDataChan*>* myCurrentChanData;
  cseis_geolib::csVector<csDataSource*>* mySourceData;
};


} // end namespace

#endif


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csP190Reader.h"
#include "csVector.h"
#include "limits"

using namespace cseis_io;

csP190Reader::csP190Reader( std::string const& filename ) {
  myFile = new std::ifstream( filename.c_str() );
  
  if( myFile == NULL || !myFile->is_open()) {
    throw( cseis_geolib::csException("Could not open file '%s'", filename.c_str() ) );
  }
  myCurrentChanData = new cseis_geolib::csVector<csDataChan*>();
  mySourceData      = new cseis_geolib::csVector<csDataSource*>();

  myCurrentSource      = -1;
  myCurrentSourceIndex = -1;

  myCurrentMinChan = -1;
  myCurrentMaxChan = -1;

  myCurrentLineNumber = 0; // Pointer to current line in input file
}

csP190Reader::~csP190Reader() {
  if( myCurrentChanData != NULL ) {
    for( int i = 0; i < myCurrentChanData->size(); i++ ) {
      delete myCurrentChanData->at(i);
    }
    delete myCurrentChanData;
    myCurrentChanData = NULL;
  }
  if( mySourceData->size() > 0 ) {
    for( int i = 0; i < mySourceData->size(); i++ ) {
      delete mySourceData->at(i);
    }
    delete mySourceData;
    mySourceData = NULL;
  }
  if( myFile != NULL ) {
    myFile->close();
    delete myFile;
    myFile = NULL;
  }
}

void csP190Reader::initialize() {
  if( mySourceData->size() > 0 ) {
    for( int i = 0; i < mySourceData->size(); i++ ) {
      delete mySourceData->at(i);
    }
  }
  mySourceData->clear();
  if( !readAllSources() ) {
    throw( cseis_geolib::csException("Unknown error occurred when trying to read in all source information from P190 file") );
  }

  myCurrentSource      = -1;
  myCurrentSourceIndex = -1;
  myCurrentChanIndex   = 0;
}
//----------------------------------------------------------------------
//
csDataSource const* csP190Reader::getSource( int source ) {
  if( source != myCurrentSource ) {
    int sourceIndexNew = getSourceIndex( source, myCurrentSourceIndex );
    if( sourceIndexNew < 0 ) return NULL;
    bool success = readChannels( source );
    if( !success ) throw( cseis_geolib::csException("Unknown problem occurred when trying to read in channel information for source %d", source) );
    myCurrentSource      = source;
    myCurrentSourceIndex = sourceIndexNew;
    myCurrentChanIndex = 0;
  }
  return mySourceData->at(myCurrentSourceIndex);
}

//----------------------------------------------------------------------
//
csDataChan const* csP190Reader::getChan( int source, int chan, int cable ) {
  if( myCurrentSource != source ) {
    csDataSource const* ds = getSource( source );
    if( ds == NULL ) return NULL;
  }
  if( myCurrentChanIndex < myCurrentChanData->size() ) {
    csDataChan* dc = myCurrentChanData->at( myCurrentChanIndex );
    if( chan == dc->chan && cable == dc->cable ) return dc;
  }
  else {
    myCurrentChanIndex = 0;
  } 
  // Search matching cable/channel starting at current position
  for( int ichan = myCurrentChanIndex; ichan < myCurrentChanData->size(); ichan++ ) {
    csDataChan* dc = myCurrentChanData->at( ichan );
    if( chan == dc->chan && cable == dc->cable ) {
      myCurrentChanIndex = ichan;
      return dc;
    }
  }
  // Search again starting at first position
  for( int ichan = 0; ichan < myCurrentChanIndex; ichan++ ) {
    csDataChan* dc = myCurrentChanData->at( ichan );
    if( chan == dc->chan && cable == dc->cable ) {
      myCurrentChanIndex = ichan;
      return dc;
    }
  }
 return NULL;
}
//----------------------------------------------------------------------
//
int csP190Reader::getSourceIndex( int source, int startIndex ) const {
  for( int i = std::max(startIndex,0); i < mySourceData->size(); i++ ) {
    if( mySourceData->at(i)->point == source ) {
      return i;
    }
  }
  for( int i = 0; i < startIndex; i++ ) {
    if( mySourceData->at(i)->point == source ) {
      return i;
    }
  }
  return -1;
}

//----------------------------------------------------------------------
// lineNum starts at 0 for first line
//
void csP190Reader::setFilePointerToStartOfLine( int lineNum ) {
  if( lineNum < myCurrentLineNumber ) {
    myFile->clear();                 // clear fail and eof bits
    myFile->seekg( 0, std::ios::beg); // back to the start
    myCurrentLineNumber = 0;
  }
  //  fprintf(stderr,"Line number: requested=%d, current=%d, EOF: %d\n", lineNum, myCurrentLineNumber, myFile->eof());
  for( int iline = myCurrentLineNumber; iline < lineNum; iline++ ) {
    if( myFile->eof() ) throw( cseis_geolib::csException("Unexpected end of file encountered in input P190 file. Line number = %d", iline+1) );
    if( !myFile->ignore( std::numeric_limits<std::streamsize>::max(), myFile->widen('\n') ) ) {
      throw( cseis_geolib::csException("Unexpected end of file encountered in input P190 file. Line number = %d", iline+1) );
    }
  }
  myCurrentLineNumber = lineNum;
}
//--------------------------------------------------
//
bool csP190Reader::readAllSources() {
  myFile->clear();                 // clear fail and eof bits
  myFile->seekg( 0, std::ios::beg); // back to the start
  int lineCounter = 0;

  while( !myFile->eof() ) {
    std::string line;
    std::getline( *myFile, line );
    if( line[0] == 'S' ) {
      csDataSource* dataSource = new csDataSource();
      scanSource( line, dataSource );
      dataSource->lineNumber = lineCounter;
      mySourceData->insertEnd( dataSource );
    }
    lineCounter += 1;
  }

  myFile->clear();                 // clear fail and eof bits
  myFile->seekg( 0, std::ios::beg); // back to the start
  return true;
}
//--------------------------------------------------
//
bool csP190Reader::readChannels( int source ) {
  if( myCurrentSource == source ) return true;

  int sourceIndexNew = getSourceIndex( source, myCurrentSourceIndex );
  if( sourceIndexNew < 0 ) return false; // Source not found

  csDataSource* dsNew = mySourceData->at(sourceIndexNew);
  setFilePointerToStartOfLine( dsNew->lineNumber );
  int endLine = (sourceIndexNew == mySourceData->size()-1) ? std::numeric_limits<int>::max() : mySourceData->at(sourceIndexNew+1)->lineNumber;

  int counterChannels = 0;
  for( int iline = myCurrentLineNumber; iline < endLine; iline++ ) {
    if( myFile->eof() ) break;
    std::string line;
    std::getline( *myFile, line );
    if( line[0] == 'R' ) {
      int numChan = scanChan( line, counterChannels );
      counterChannels += numChan;
    } // == 'R'
  } // END for iline

  myCurrentLineNumber = endLine;
  for( int ichan = counterChannels; ichan < myCurrentChanData->size(); ichan++ ) {
    delete myCurrentChanData->at(ichan);
  }
  myCurrentChanData->remove( counterChannels, myCurrentChanData->size()-counterChannels );

  //  myCurrentMinChan = myCurrentChanData->at(0)->chan;
  //  myCurrentMaxChan = myCurrentChanData->at( myCurrentChanData->size()-1 )->chan;

  return true;
}

//----------------------------------------------------------------------
//
void csP190Reader::scanSource( std::string& line, csDataSource* data ) {
  data->id    = atoi( line.substr(17,1).c_str() );
  data->point = atoi( line.substr(19,6).c_str() );
  data->x     = atof( line.substr(46,9).c_str() );
  data->y     = atof( line.substr(55,9).c_str() );
  data->wdep  = atof( line.substr(64,6).c_str() );
  data->day   = atoi( line.substr(70,3).c_str() );
  data->hour  = atoi( line.substr(73,2).c_str() );
  data->min   = atoi( line.substr(75,2).c_str() );
  data->sec   = atoi( line.substr(77,2).c_str() );
}

//----------------------------------------------------------------------
//
int csP190Reader::scanChan( std::string& line, int index ) {
  int numChan = 1;

  if( index >= myCurrentChanData->size() ) {
    myCurrentChanData->insertEnd( new csDataChan() );
  }
  csDataChan* data1 = myCurrentChanData->at(index+0);
  // Cable number is stored as hexadecimal, single letter
  data1->cable = (int)line[79];
  if( data1->cable <= 57 ) {
    data1->cable -= 48;
  }
  else {
    data1->cable -= 55;
  }

  data1->chan = atoi( line.substr(1,4).c_str() );
  data1->x    = atof( line.substr(5,9).c_str() );
  data1->y    = atof( line.substr(14,9).c_str() );
  data1->z    = atof( line.substr(23,4).c_str() );

  int chan = (int)line[30];
  if( chan >= 48 && chan <= 57 ) {
    if( index+1 >= myCurrentChanData->size() ) {
      myCurrentChanData->insertEnd( new csDataChan() );
    }
    csDataChan* data2 = myCurrentChanData->at(index+1);
    data2->chan = atoi( line.substr(27,4).c_str() );
    data2->x    = atof( line.substr(31,9).c_str() );
    data2->y    = atof( line.substr(40,9).c_str() );
    data2->z    = atof( line.substr(49,4).c_str() );
    data2->cable = data1->cable;
    numChan = 2;

    chan = (int)line[56];
    if( chan >= 48 && chan <= 57 ) {
      if( index+2 >= myCurrentChanData->size() ) {
        myCurrentChanData->insertEnd( new csDataChan() );
      }
      csDataChan* data3 = myCurrentChanData->at(index+2);
      data3->chan = atoi( line.substr(53,4).c_str() );
      data3->x    = atof( line.substr(57,9).c_str() );
      data3->y    = atof( line.substr(66,9).c_str() );
      data3->z    = atof( line.substr(75,4).c_str() );
      data3->cable = data1->cable;
      numChan = 3;
    }
  }

  return numChan;
}

void csP190Reader::dump() const {
  fprintf(stdout,"*** P1/90 source data dump ***\n" );
  for( int i = 0; i < mySourceData->size(); i++ ) {
    fprintf(stdout,"#%-4d %d\n", i, mySourceData->at(i)->point);
  }
}
void csP190Reader::dumpAllSourceInfo() const {
  fprintf(stdout,"*** P1/90 data dump for source %d (numChan=%d) ***\n", myCurrentSource, myCurrentChanData->size() );
  for( int i = 0; i < myCurrentChanData->size(); i++ ) {
    csDataChan* dc = myCurrentChanData->at(i);
    fprintf(stdout,"#%-5d: Cable/chan %d/%-5d\n", i+1, dc->cable, dc->chan);
  }
}

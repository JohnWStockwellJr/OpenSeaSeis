/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEISMIC_READER_H
#define CS_SEISMIC_READER_H

#include <cstdio>
#include <string>
#include "geolib_defines.h"
#include "csIReader.h"

namespace cseis_io {
  class csSeismicReader_ver;
}

namespace cseis_geolib {
  class csFlexHeader;
  class csIOSelection;
}
namespace cseis_system {

class csSuperHeader;
class csTraceHeaderDef;

/**
 * Seismic file reader interface, Cseis format
 * This class purely provides an interface for the seismic reader class defined in namespace cseis_io.
 * It also provides functionality to translate header information between the (generic) reader class and the
 * Cseis base system.
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
 class csSeismicReader : public cseis_geolib::csIReader {
public:
  /**
   * Constructor
   * @param filename  Name of the seismic file
   * @param numTraces Number of traces to be buffered at once
   */
  csSeismicReader( std::string filename, int numTraces );
  /**
   * Constructor
   * @param filename           Name of the seismic file
   * @param enableRandomAccess true if random access to the file shall be enabled
   * @param numTraces          Number of traces to be buffered at once
   */
  csSeismicReader( std::string filename, bool enableRandomAccess, int numTraces );
  ~csSeismicReader();
  /**
  * Read all header information.
  * The 'file' header is the one at the front of the seismic data file, giving information on...see arguments
  * Call once to retrieve general file information
  * @param shdr  (o) Cseis super header
  * @param hdef  (o) Cseis trace header definition object
  * @param hdrValueBlockSize  (o) Size in bytes of seismic trace header value block that sits on top of every seismic trace
  * @param stream (i) Error output stream
  */
  bool readFileHeader( csSuperHeader* shdr, csTraceHeaderDef* hdef, int* hdrValueBlockSize, std::FILE* stream = NULL );
  /**
  * Read next trace
  * @param samples  (o) Trace samples
  * @param hdrValueBlock (o) Buffer to hold all trace header values in the format defined in the trace header definition
   * @return false if something went wrong.
  */
  bool readTrace( float* samples, char* hdrValueBlock );
  bool readTrace( float* samples, char* hdrValueBlock, int numSamples );
  /**
   * Set header to peek
   * Initiates header 'peeking' operation. Next, the method 'peekheaderValue' may be called to retrieve the header value
   * for an arbitrary trace in the file.
   * @param @headerName Name of the header to peek
   * @return false if something went wrong.
   */
  bool setHeaderToPeek( std::string const& headerName );
  int getPeekByteOffset() const;
  /**
   * Set header to peek
   * Initiates header 'peeking' operation. Next, the method 'peekheaderValue' may be called to retrieve the header value
   * for an arbitrary trace in the file.
   * @param @headerName Name of the header to peek
   * @param @headerType Type of the header to peek
   * @return false if something went wrong.
   */
  bool setHeaderToPeek( std::string const& headerName, cseis_geolib::type_t& headerType );
  /**
   * Peek header value of arbitrary trace
   * @param hdrValue   Header value
   * @param traceIndex Trace index of trace to peek header
   * @return false     if something went wrong.
   */
  bool peekHeaderValue( cseis_geolib::csFlexHeader* hdrValue, int traceIndex = -1 );
  /** 
   * @return Number of traces in file
   */
  int numTraces() const { return myNumTraces; }
  /**
   * @return Number of traces that can be buffered at once
   */
  int numTracesCapacity() const;
  /**
   * Move file pointer to arbitrary trace.
   * @param traceIndex      Index of trace to move to
   * @param numTracesToRead Number of traces to read into buffer at once after move operation
   */
  bool moveToTrace( int traceIndex, int numTracesToRead );
  /**
   * Move file pointer to arbitrary trace.
   * @param traceIndex      Index of trace to move to
   */
  bool moveToTrace( int traceIndex );
  /**
   * @return Number of samples in file
   */
  int numSamples() const;

  /**
   * Set specific trace selection.
   * Only traces matching the selection will be read in and sorted, if requested.
   * @param hdrValueSelectionText  A string defining the selection of certain trace header values.
   * @param headerName  Header name to select/sort on
   * @param sortOrder   csIOSelection::SORT_NONE, SORT_INCREASING, or SORT_DECREASING
   * @param sortMethod  csSortManager::SIMPLE_SORT or csSortManager::TREE_SORT
   */
  bool setSelection( std::string const& hdrValueSelectionText, std::string const& headerName, int sortOrder, int sortMethod );
  int getCurrentTraceIndex() const;

private:
  void init();
  bool performIOSelection();

  cseis_io::csSeismicReader_ver* myReader;
  int myNumTraces;
  /// Header to check. 
  int  myHdrCheckByteOffset;
  cseis_geolib::type_t myHdrCheckType;
  int    myHdrCheckByteSize;
  char* myHdrCheckBuffer;
  cseis_system::csTraceHeaderDef const* myTrcHdrDef;  // Pointer only, do not free!
  cseis_geolib::csIOSelection* myIOSelection;
};

} // end namespace
#endif


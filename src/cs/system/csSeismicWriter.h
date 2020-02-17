/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_SEISMIC_WRITER_H
#define CS_SEIMSIC_WRITER_H

#include <cstdio>
#include <string>

namespace cseis_io {
  class csSeismicWriter_ver;
}

namespace cseis_system {

class csSuperHeader;
class csTraceHeaderDef;

/**
 * Seismic file writer, Cseis format
 * This class purely provides an interface for the seismic writer class defined in namespace cseis_io.
 * It also provides functionality to translate header information between the (generic) writer class and the
 * Cseis base system.
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
class csSeismicWriter {
public:
  /**
   * Constructor
   * @param filename Name of output file
   * @param numTracesBuffer Number of traces to buffer at once
   * @param sampleByteSize  Byte size of one sample. 4: No compression, 2: 2x compression, 1: 4x compression
   * @param overwrite       true to overwrite even if file exists. false if existing file shall not be overwritten.
   */
  csSeismicWriter( std::string filename, int numTracesBuffer, int sampleByteSize = 4, bool overwrite = true );
  ~csSeismicWriter();

  /**
   * Write file header
   * Call only once before starting to write trace data
   * @param shdr  (i) Cseis super header
   * @param hdef  (i) Cseis trace header definition object
   */
  bool writeFileHeader( csSuperHeader const* shdr, csTraceHeaderDef const* hdef );
  /**
   * Write next trace
   * @param samples  (i) Trace samples
   * @param hdrValueBlock (i) Buffer holding all trace header values in the format defined in the trace header definition
   */
  bool writeTrace( float* samples, char const* hdrValueBlock );

private:
  cseis_io::csSeismicWriter_ver* myWriter;
  char* myHdrTempBuffer;
  csTraceHeaderDef const* myHdef;
};

} // end namespace
#endif


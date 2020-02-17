/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexHeader.h"
#include "csSplineInterpolation2D.h"
#include "csRandomDataReader.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: INPUT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_input_interpolate {

  template <typename T> class Point {
  public:
    Point( T r, T c) {
      row = r;
      col = c;
    }
    Point() {
      row = 0;
      col = 0;
    }
    T row;
    T col;
  };

  template <typename T> class GridDef {
  public:
    GridDef( T n1, T n2, T inc ) {
      val1 = n1;
      val2 = n2;
      increment = inc;
      num = (int)round( (val2-val1)/inc ) + 1;
    }
    GridDef() {
      val1 = 0;
      val2 = 0;
      increment = 0;
      num = 0;
    }
    void computeNum() {
      num = (int)round( (val2-val1)/increment ) + 1;
    }
    T val1;
    T val2;
    T increment;
    int num;
  };

  struct VariableStruct {
    csRandomDataReader* reader;

    int hdrId_row;
    int hdrId_col;
    int hdrId_rowOut;
    int hdrId_colOut;
    type_t hdrType_row;
    type_t hdrType_col;
    std::string hdrName_row;
    std::string hdrName_col;

    int method;
    int methodLineDef;
    int numTracesOut;
    int numTracesIn;
    Point<int> incIn;
    Point<int> p1In;
    Point<int> p2In;

    int numRowsIn;
    int numColsIn;

    Point<double>* pointsToCompute;

    bool isColFastIndex;
    bool atEOF;
    long traceCounter;
    int nTracesToRead;

    int numValuesSplineBuffer;
    float** splineSampleBuffer;
    Point<int> p1Buffer;
    Point<int> p2Buffer;

    csFlexHeader* rowValuesIn;
    csFlexHeader* colValuesIn;
  };
  static int const DEFINE_LINE_ROWCOL  = 11;
  static int const DEFINE_GRID_ROWCOL  = 12;
  static int const METHOD_SPLINE        = 21;
  static int const METHOD_SPLINE_MINMAX = 22;
  static int const METHOD_NEAREST       = 23;
  static int const METHOD_BILINEAR      = 24;

  int getTraceIndex( mod_input_interpolate::VariableStruct* vars, int row, int col ) {
    // Compute 'reduced' row/col numbers which adhere to input data range in order to be able to work close to edges
    int rowRed = std::min( std::max( row, vars->p1In.row ), vars->p2In.row );
    int colRed = std::min( std::max( col, vars->p1In.col ), vars->p2In.col );
    int drow = (rowRed - vars->p1In.row)/vars->incIn.row;
    int dcol = (colRed - vars->p1In.col)/vars->incIn.col;
    //    int numRows = (vars->p2In.row - vars->p1In.row)/vars->incIn.row + 1;
    //    int numCols = (vars->p2In.col - vars->p1In.col)/vars->incIn.col + 1;
    //    if( row == 2120 ) {
    //  fprintf(stderr,"row/col: %d %d, isColFast: %d, drow/dcol: %d %d\n", row, col, vars->isColFastIndex, drow, dcol);
    // }
    if( vars->isColFastIndex ) {
      return( drow * vars->numColsIn + dcol );
    }
    else {
      return( dcol * vars->numRowsIn + drow );
    }
  } // END getTraceIndex
}
using mod_input_interpolate::VariableStruct;


//*********************************************************************************bool****************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_interpolate_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->pointsToCompute = NULL;

  vars->reader       = NULL;
  vars->atEOF         = false;

  vars->traceCounter      = 0;
  vars->nTracesToRead  = 10;

  vars->methodLineDef = mod_input_interpolate::DEFINE_LINE_ROWCOL;
  vars->method = mod_input_interpolate::METHOD_SPLINE;
  vars->hdrId_row = -1;
  vars->hdrId_col = -1;
  vars->hdrId_rowOut = -1;
  vars->hdrId_colOut = -1;
  vars->numTracesIn = 0;
  vars->atEOF = false;

  vars->numValuesSplineBuffer = 4;
  vars->incIn.row = 0;
  vars->incIn.col = 0;
  vars->numRowsIn = 0;
  vars->numColsIn = 0;

  vars->isColFastIndex = true;

//------------------------------------------------------------
  std::string text;
  if( param->exists( "define" ) ) {
    param->getString( "define", &text );
    if( !text.compare("line") ) {
      vars->methodLineDef = mod_input_interpolate::DEFINE_LINE_ROWCOL;
    }
    else if( !text.compare("grid") ) {
      vars->methodLineDef = mod_input_interpolate::DEFINE_GRID_ROWCOL;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists( "method" ) ) {
    param->getString( "method", &text );
    if( !text.compare("spline") ) {
      vars->method = mod_input_interpolate::METHOD_SPLINE;
      vars->numValuesSplineBuffer = 4;
    }
    else if( !text.compare("spline_minmax") ) {
      vars->method = mod_input_interpolate::METHOD_SPLINE_MINMAX;
      vars->numValuesSplineBuffer = 4;
    }
    else if( !text.compare("nearest") ) {
      vars->method = mod_input_interpolate::METHOD_NEAREST;
      vars->numValuesSplineBuffer = 2;
    }
    else if( !text.compare("bilinear") ) {
      vars->method = mod_input_interpolate::METHOD_BILINEAR;
      vars->numValuesSplineBuffer = 2;
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  std::string filename;
  param->getString( "filename", &filename );

  vars->hdrName_row = "row";
  vars->hdrName_col = "col";
  if( param->exists( "hdrin_rowcol" ) ) {
    param->getString( "hdrin_rowcol", &vars->hdrName_row, 0 );
    param->getString( "hdrin_rowcol", &vars->hdrName_col, 1 );
  }

  if( vars->methodLineDef == mod_input_interpolate::DEFINE_LINE_ROWCOL ) {
    mod_input_interpolate::Point<double> line_p1;
    mod_input_interpolate::Point<double> line_p2;
    mod_input_interpolate::Point<double> line_delta;
    double trc_inc_m;
    double grid_binsize_il;
    double grid_binsize_xl;

    param->getDouble( "line_p1", &line_p1.row, 0 );
    param->getDouble( "line_p1", &line_p1.col, 1 );
    param->getDouble( "line_p2", &line_p2.row, 0 );
    param->getDouble( "line_p2", &line_p2.col, 1 );
    param->getDouble( "trc_inc", &trc_inc_m, 0 );

    param->getDouble( "line_binsize", &grid_binsize_il, 0 );
    param->getDouble( "line_binsize", &grid_binsize_xl, 1 );
    if( grid_binsize_il <= 0 || grid_binsize_xl <= 0 ) writer->error("Inconsistent bin size specified. Must be larger than zero. IL size: %f, XL size: %f", grid_binsize_il, grid_binsize_xl);

    double drow = line_p2.row - line_p1.row;
    double dcol = line_p2.col - line_p1.col;

    double dist_row_m = drow * grid_binsize_xl;
    double dist_col_m = dcol * grid_binsize_il;

    double line_length = sqrt( pow(dist_row_m,2) + pow(dist_col_m,2) );
    vars->numTracesOut = (int)round( line_length / trc_inc_m ) + 1;

    // Recompute last point:
    double length_correction = ( (double)(vars->numTracesOut-1) * trc_inc_m ) / line_length;
    dist_row_m *= length_correction;
    dist_col_m *= length_correction;

    line_delta.row = ( dist_row_m / grid_binsize_xl ) / (double)(vars->numTracesOut-1);
    line_delta.col = ( dist_col_m / grid_binsize_il ) / (double)(vars->numTracesOut-1);
    line_p2.row = line_p1.row + line_delta.row * (double)(vars->numTracesOut-1);
    line_p2.col = line_p1.col + line_delta.col * (double)(vars->numTracesOut-1);

    vars->pointsToCompute = new mod_input_interpolate::Point<double>[vars->numTracesOut];
    for( int itrc = 0; itrc < vars->numTracesOut; itrc++ ) {
      vars->pointsToCompute[itrc].row = (double)itrc * line_delta.row + line_p1.row;
      vars->pointsToCompute[itrc].col = (double)itrc * line_delta.col + line_p1.col;
    }

    writer->line("Random line definition:");
    writer->line("  Point #1:  %.3f %.3f", line_p1.row, line_p1.col);
    writer->line("  Point #2:  %.3f %.3f", line_p2.row, line_p2.col );
    writer->line("  Increment: %.3f %.3f", line_delta.row, line_delta.col);
    writer->line("  Line length: %.2fm, numTraces: %d  (length correction term: %.4e)", line_length, vars->numTracesOut, length_correction);
  }
  else if( vars->methodLineDef == mod_input_interpolate::DEFINE_GRID_ROWCOL ) {
    mod_input_interpolate::GridDef<double> griddef_row;
    mod_input_interpolate::GridDef<double> griddef_col;

    param->getDouble( "grid_row", &griddef_row.val1, 0 );
    param->getDouble( "grid_row", &griddef_row.val2, 1 );
    param->getDouble( "grid_row", &griddef_row.increment, 2 );
    param->getDouble( "grid_col", &griddef_col.val1, 0 );
    param->getDouble( "grid_col", &griddef_col.val2, 1 );
    param->getDouble( "grid_col", &griddef_col.increment, 2 );

    bool isRow = true;
    if( param->exists( "grid_sort" ) ) {
      param->getString( "grid_sort", &text );
      if( !text.compare("row") ) {
        isRow = true;
      }
      else if( !text.compare("col") ) {
        isRow = false;
      }
      else {
        writer->error("Unknown option: %s", text.c_str());
      }
    }

    griddef_row.computeNum();
    griddef_col.computeNum();

    vars->numTracesOut = griddef_row.num * griddef_col.num;
    vars->pointsToCompute = new mod_input_interpolate::Point<double>[vars->numTracesOut];
    int counter = 0;
    if( isRow ) {
      for( int irow = 0; irow < griddef_row.num; irow++ ) {
        for( int icol = 0; icol < griddef_col.num; icol++ ) {
          vars->pointsToCompute[counter].row = griddef_row.val1 + griddef_row.increment * irow;
          vars->pointsToCompute[counter++].col = griddef_col.val1 + griddef_col.increment * icol;
        }
      }
    }
    else {
      for( int icol = 0; icol < griddef_col.num; icol++ ) {
        for( int irow = 0; irow < griddef_row.num; irow++ ) {
          vars->pointsToCompute[counter].row = griddef_row.val1 + griddef_row.increment * irow;
          vars->pointsToCompute[counter++].col = griddef_col.val1 + griddef_col.increment * icol;
        }
      }
    }
    writer->line("Grif definition:");
    writer->line("  Point #1:  %.3f %.3f", griddef_row.val1, griddef_col.val1 );
    writer->line("  Point #2:  %.3f %.3f", griddef_row.val2, griddef_col.val2 );
    writer->line("  Increment: %.3f %.3f", griddef_row.increment, griddef_col.increment );
    writer->line("  NumTraces: %d ", vars->numTracesOut);
  }
  else {
    writer->error("Method not implemented...");
  }

  if( edef->isDebug() ) {
    fprintf(stdout,"LIST OF POINTS TO COMPUTE:\n");
    for( int itrc = 0; itrc < vars->numTracesOut; itrc++ ) {
      fprintf(stdout,"Point %d: %f %f\n", itrc+1, vars->pointsToCompute[itrc].row,  vars->pointsToCompute[itrc].col );
    }
  }
  //-------------------------------------------------------------------------
  // Open input file
  //
  bool enableRandomAccess = true;
  int numTracesBuffer = 30;

  try {
    vars->reader = new cseis_system::csRandomDataReader( filename, enableRandomAccess, numTracesBuffer );
    bool success = vars->reader->readFileHeader( shdr, hdef, writer->getFile() );
    if( !success ) {
      writer->error("Unknown error occurred when reading file header from SeaSeis file '%s'.\n", filename.c_str() );
    }
    vars->numTracesIn = vars->reader->numTraces();
    if( vars->numTracesIn < vars->numValuesSplineBuffer*2 ) writer->error("Input data contains too few traces: %d", vars->numTracesIn);

    vars->rowValuesIn = new csFlexHeader[vars->numTracesIn];
    vars->colValuesIn = new csFlexHeader[vars->numTracesIn];
  }
  catch( csException& exc ) {
    writer->error("Error occurred when opening SeaSeis file '%s'. System message:\n%s", filename.c_str(), exc.getMessage() );
  }

  //--------------------------------------------------------------------------------

  std::string hdrName_rowOut("row_dbl");
  std::string hdrName_colOut("col_dbl");
  if( param->exists( "hdrout_rowcol" ) ) {
    param->getString( "hdrout_rowcol", &hdrName_rowOut, 0 );
    param->getString( "hdrout_rowcol", &hdrName_colOut, 1 );
  }

  if( !hdef->headerExists(hdrName_rowOut) ){
    hdef->addHeader( cseis_geolib::TYPE_DOUBLE, hdrName_rowOut );
  }
  else if( hdef->headerType(hdrName_rowOut) != cseis_geolib::TYPE_DOUBLE && hdef->headerType(hdrName_rowOut) != cseis_geolib::TYPE_FLOAT ) {
    writer->error("Output row trace header '%s' must be of floating point type. Choose different trace header existing in input file, or specify a new trace header", hdrName_rowOut.c_str());
  }
  if( !hdef->headerExists(hdrName_colOut) ){
    hdef->addHeader( cseis_geolib::TYPE_DOUBLE, hdrName_colOut );
  }
  else if( hdef->headerType(hdrName_rowOut) != cseis_geolib::TYPE_DOUBLE && hdef->headerType(hdrName_rowOut) != cseis_geolib::TYPE_FLOAT ) {
    writer->error("Output col trace header '%s' must be of floating point type. Choose different trace header existing in input file, or specify a new trace header", hdrName_colOut.c_str());
  }

  // ...important to do the following after all headers have been set for hdef and all other input files
  hdef->resetByteLocation();  // This will otherwise be done by base system AFTER init phase

  vars->hdrId_row = hdef->headerIndex( vars->hdrName_row );
  vars->hdrId_col = hdef->headerIndex( vars->hdrName_col );
  vars->hdrType_row = hdef->headerType( vars->hdrName_row );
  vars->hdrType_col = hdef->headerType( vars->hdrName_col );
  vars->hdrId_rowOut = hdef->headerIndex( hdrName_rowOut );
  vars->hdrId_colOut = hdef->headerIndex( hdrName_colOut );

  vars->splineSampleBuffer = new float*[vars->numValuesSplineBuffer * vars->numValuesSplineBuffer];
  int num = vars->numValuesSplineBuffer * vars->numValuesSplineBuffer;
  for( int i = 0; i < num; i++ ) {
    vars->splineSampleBuffer[i] = new float[shdr->numSamples];
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      vars->splineSampleBuffer[i][isamp] = 0;
    }
  }

  vars->traceCounter = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_input_interpolate_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const*  shdr = env->superHeader;

  if( vars->atEOF ) return false;
  if( vars->numTracesOut > 0 && vars->numTracesOut == vars->traceCounter ) {
    vars->atEOF = true;
    return false;
  }

  csTraceHeader* trcHdr = traceGather->trace(0)->getTraceHeader();
  float* samplesOut = traceGather->trace(0)->getTraceSamples();

  //********************************************************************************
  // For first trace, read in row/col trace header values from all input traces.
  // --> More efficiently, assume rectangular input cube: Read in only first, second and last trace and compute all row/col values including row/col increments.
  // Initialize all fields for retrieval of all necessary input traces.
  //
  if( vars->traceCounter == 0 ) {
    csFlexHeader rowVal_trc1;
    csFlexHeader rowVal_trc2;
    csFlexHeader rowVal_trcN;
    csFlexHeader colVal_trc1;
    csFlexHeader colVal_trc2;
    csFlexHeader colVal_trcN;

    bool success;
    success = vars->reader->setHeaderToPeek( vars->hdrName_row );
    if( success ) success = vars->reader->peekHeaderValue( &rowVal_trc1, 0 );
    if( success ) success = vars->reader->peekHeaderValue( &rowVal_trc2, 1 );
    if( success ) success = vars->reader->peekHeaderValue( &rowVal_trcN, vars->numTracesIn-1 );
    if( !vars->reader->moveToTrace( 0 ) ) writer->error("Error occurred when moving file pointer to first trace");
    if( success ) success = vars->reader->setHeaderToPeek( vars->hdrName_col );
    if( success ) success = vars->reader->peekHeaderValue( &colVal_trc1, 0 );
    if( success ) success = vars->reader->peekHeaderValue( &colVal_trc2, 1 );
    if( success ) success = vars->reader->peekHeaderValue( &colVal_trcN, vars->numTracesIn-1 );
    if( !success ) {
      writer->error("Error occurred when retrieving trace headers '%s' and '%s' for first trace",
                    vars->hdrName_row.c_str(), vars->hdrName_col.c_str() );
    }
  
    vars->p1In.row = rowVal_trc1.intValue();
    vars->p2In.row = rowVal_trcN.intValue();
    int row_trc2 = rowVal_trc2.intValue();

    vars->p1In.col = colVal_trc1.intValue();
    vars->p2In.col = colVal_trcN.intValue();
    int col_trc2 = colVal_trc2.intValue();

    if( vars->p1In.row == row_trc2 ) {
      if( vars->p1In.col == col_trc2 ) {
        writer->error("Unexpected input row/col sort order encountered. Expected either rows or columns to increase between trace 1 and 2.\nEncountered values: row(%d) col(%d)", vars->p1In.row, vars->p1In.col);
      }
      vars->isColFastIndex = true;
      vars->incIn.col = col_trc2 - vars->p1In.col;
      int numCols = (int)( ( vars->p2In.col - vars->p1In.col ) / vars->incIn.col ) + 1;
      int numRows = vars->numTracesIn / numCols;
      if( numRows * numCols != vars->numTracesIn ) writer->error("Unexpected input row/col sort order encountered. Expected rectangular volume with regular rows & columns.\nEncountered: First row/col: %d/%d, Last row/col: %d %d   (fast columns)",
                                                                 vars->p1In.row, vars->p1In.col, vars->p2In.row, vars->p2In.col);
      if( numRows > 1 ) {
        vars->incIn.row = ( vars->p2In.row - vars->p1In.row ) / ( numRows - 1 );
      }
      else {
        vars->incIn.row = 1;
      }
    }
    else {
      vars->isColFastIndex = false;
      vars->incIn.row = row_trc2 - vars->p1In.row;
      int numRows = (int)( ( vars->p2In.row - vars->p1In.row ) / vars->incIn.row ) + 1;
      int numCols = vars->numTracesIn / numRows;
      if( numCols * numRows != vars->numTracesIn ) writer->error("Unexpected input row/col sort order encountered. Expected rectangular volume with regular rows & columns.\nEncountered: First row/col: %d/%d, Last row/col: %d %d  (fast rows)",
                                                                 vars->p1In.row, vars->p1In.col, vars->p2In.row, vars->p2In.col);
      if( numCols > 1 ) {
        vars->incIn.col = ( vars->p2In.col - vars->p1In.col ) / ( numCols - 1 );
      }
      else {
        vars->incIn.col = 1;
      }
    }

    vars->numRowsIn = (vars->p2In.row - vars->p1In.row)/vars->incIn.row + 1;
    vars->numColsIn = (vars->p2In.col - vars->p1In.col)/vars->incIn.col + 1;

    writer->line("Row/col values:  Trace #1: %d/%d  #2: %d/%d  #%d: %d/%d\n", vars->p1In.row, vars->p1In.col, row_trc2, col_trc2, vars->numTracesIn, vars->p2In.row, vars->p2In.col );
    writer->line("Row/col increment: %d / %d,  #rows/cols: %d / %d\n", vars->incIn.row, vars->incIn.col, vars->numRowsIn, vars->numColsIn );

    // Read in header value block of first trace
    if( !vars->reader->moveToTrace( 0 ) ) writer->error("Error occurred when moving file pointer to first trace");
    if( !vars->reader->readTrace( vars->splineSampleBuffer[0], shdr->numSamples, trcHdr ) ) writer->error("Error occurred when reading in first trace");

    // Check user-specified output trace range. CHANGE: Remove traces which cannot be interpolated from input data, i.e. which fall outside input data row/col range
    double rowMin = vars->pointsToCompute[0].row;
    double rowMax = rowMin;
    double colMin = vars->pointsToCompute[0].col;
    double colMax = colMin;
    for( int itrc = 1; itrc < vars->numTracesOut; itrc++ ) {
      int row = vars->pointsToCompute[itrc].row;
      int col = vars->pointsToCompute[itrc].col;
      if( row > rowMax ) rowMax = row;
      if( row < rowMin ) rowMin = row;
      if( col > colMax ) colMax = col;
      if( col < colMin ) colMin = col;
    }
    if( rowMin < vars->p1In.row || rowMax > vars->p2In.row || colMin < vars->p1In.col || colMax > vars->p2In.col ) {
      writer->error("Specified row/col output range (rows %.1f-%.1f/cols %.1f-%.1f) exceeds...\n...row/col range in input file (rows %d-%d/cols %d-%d)", rowMin, rowMax, colMin, colMax, vars->p1In.row, vars->p2In.row, vars->p1In.col, vars->p2In.col );
    }
  }
  // END: Read in (peek) all required trace header values
  //********************************************************************************

  int totalNumTracesBuffer = vars->numValuesSplineBuffer*vars->numValuesSplineBuffer;
  double rowOut = 0;
  double colOut = 0;
  int col1 = 0;
  int row1 = 0;

  // Outer loop is repeated if current output trace cannot be computed since its contributing input traces are outside of area covered by input data cube
  bool traceIsOutsideArea = false;
  do {
    traceIsOutsideArea = false; // true if at least one trace required to interplate current output trace is outside input data area
    // Decide which traces to read in
    // rowOut & colOut: The actual row/col where we need a trace (row/col in binning grid of input file)
    if( vars->traceCounter >= vars->numTracesOut ) writer->error("Maximum trace counter reached: %d", vars->traceCounter+1 );
    rowOut = vars->pointsToCompute[vars->traceCounter].row;
    colOut = vars->pointsToCompute[vars->traceCounter].col;
    //    rowOut = (double)vars->traceCounter * vars->line_delta.row + vars->line_p1.row;
    //    colOut = (double)vars->traceCounter * vars->line_delta.col + vars->line_p1.col;
    if( edef->isDebug() ) fprintf(stdout,"...processing trace %ld:  %.2f %.2f\n", vars->traceCounter, rowOut, colOut);

    int numValHalf = (int)( (vars->numValuesSplineBuffer-1)/2 );

    // row1Prev & col1Prev etc: Patch of input traces used to interpolate previous output trace
    int row1Prev = vars->p1Buffer.row;
    int col1Prev = vars->p1Buffer.col;
    int row2Prev = vars->p2Buffer.row;
    int col2Prev = vars->p2Buffer.col;

    // row1/col1 & row2/col2: Min/max row/col numbers (NxN traces) that are needed for spline interpolation
    row1 = (int)( (rowOut-vars->p1In.row)/vars->incIn.row ) * vars->incIn.row + vars->p1In.row - numValHalf * vars->incIn.row;
    col1 = (int)( (colOut-vars->p1In.col)/vars->incIn.col ) * vars->incIn.col + vars->p1In.col - numValHalf * vars->incIn.col;
    int row2 = row1 + (vars->numValuesSplineBuffer-1) * vars->incIn.row;
    int col2 = col1 + (vars->numValuesSplineBuffer-1) * vars->incIn.col;

    vars->p1Buffer.row = row1;
    vars->p1Buffer.col = col1;
    vars->p2Buffer.row = row2;
    vars->p2Buffer.col = col2;

    //  fprintf(stderr,"row/col out: %f %f,  area:  %d %d  %d %d   arePrev: %d %d  %d %d\n", rowOut, colOut, row1, col1, row2, col2, row1Prev, col1Prev, row2Prev, col2Prev);
    if( row1 == row1Prev && col1 == col1Prev ) {
      // Nothing to be done. All traces have already been read into buffer. Current output trace is interpolated from the same patch of input traces
    }
    else { // Need to reshuffle buffer and read in some more traces from input file
      // Determine overlapping row & col interval
      int row1Overlap = ( row1 > row1Prev ) ? row1 : row1Prev;
      int row2Overlap = ( row2 > row2Prev ) ? row2Prev : row2;
      int col1Overlap = ( col1 > col1Prev ) ? col1 : col1Prev;
      int col2Overlap = ( col2 > col2Prev ) ? col2Prev : col2;
      
      //    fprintf(stderr,"Overlap: %d %d %d %d\n", row1Overlap, col1Overlap, row2Overlap, col2Overlap);
      if( vars->traceCounter == 0 ) { // For first input trace, jumble overlap to ensure all traces are read in, not copied from buffer
        row1Overlap = row2;
        row2Overlap = row1;
        col1Overlap = col2;
        col2Overlap = col1;
      }

      // Test if all required traces exist in input data: Compute required trace number (traceIndex) to see if input data actually has all necessary traces
      /*
      for( int irow = 0; irow < vars->numValuesSplineBuffer; irow++ ) {
        int row = irow*vars->incIn.row + row1;
        bool isRowOverlap = ( row >= row1Overlap && row <= row2Overlap );
        for( int icol = 0; icol < vars->numValuesSplineBuffer; icol++ ) {
          int col = icol*vars->incIn.col + col1;
          if( isRowOverlap && (col >= col1Overlap && col <= col2Overlap) ) continue;
          int traceIndex = getTraceIndex( vars, row, col );
          fprintf(stdout," === %d %d  %d %d\n", row, col, isRowOverlap, traceIndex);
          if( traceIndex < 0 || traceIndex >= vars->numTracesIn ) {
            if( traceIndex >= vars->numTracesIn ) {
              vars->traceCounter = vars->numTracesOut+1;
              fprintf(stdout,"END: Trace not found??? %d\n", traceIndex );
              return false;
            }
            fprintf(stdout," === %d %d  %d %d\n", row, col, isRowOverlap, traceIndex);
            traceIsOutsideArea = true;
            vars->p1Buffer.row = 0;
            vars->p1Buffer.col = 0;
            vars->p2Buffer.row = 0;
            vars->p2Buffer.col = 0;
            break;
          }
        }
        if( traceIsOutsideArea ) break;
      }
      if( traceIsOutsideArea ) fprintf(stdout," === Trace is outside area %.2f %.2f\n", rowOut, colOut);
      */
      traceIsOutsideArea = false;

      if( !traceIsOutsideArea ) {
        // Temp buffer: Buffer where trace samples are read into first. These are then later moved over to the permanent array vars->splineSampleBuffer
        float** sampleBufferTemp = new float*[totalNumTracesBuffer];
        for( int row = row1Overlap; row <= row2Overlap; row += vars->incIn.row ) {
          int rowIndexPrev = (row - row1Prev)/vars->incIn.row;
          int rowIndexNow  = (row - row1)/vars->incIn.row;
          //        fprintf(stdout," ROW prev/now index:  %d %d   - %d\n", rowIndexPrev, rowIndexNow, row);
          for( int col = col1Overlap; col <= col2Overlap; col += vars->incIn.col ) {
            int colIndexPrev = (col - col1Prev)/vars->incIn.col;
            int colIndexNow  = (col - col1)/vars->incIn.col;
            //        fprintf(stdout," COL prev/now index:  %d %d   - %d\n", colIndexPrev, colIndexNow, col);
            int fullIndexPrev = rowIndexPrev*vars->numValuesSplineBuffer + colIndexPrev;
            int fullIndexNow  = rowIndexNow*vars->numValuesSplineBuffer + colIndexNow;
            //        fprintf(stdout," prev/now index:      %d %d   %d %d %d   %d %d %d\n", fullIndexPrev, fullIndexNow,rowIndexPrev, rowIndexNow, row, colIndexPrev, colIndexNow, col );
            sampleBufferTemp[ fullIndexNow ] = vars->splineSampleBuffer[ fullIndexPrev ];
            vars->splineSampleBuffer[ fullIndexPrev ] = NULL;
          }
        }
        // Read in new traces
        int counter = 0;
        for( int irow = 0; irow < vars->numValuesSplineBuffer; irow++ ) {
          int row = irow*vars->incIn.row + row1;
          bool isRowOverlap = ( row >= row1Overlap && row <= row2Overlap );
          for( int icol = 0; icol < vars->numValuesSplineBuffer; icol++ ) {
            int col = icol*vars->incIn.col + col1;
            if( isRowOverlap && (col >= col1Overlap && col <= col2Overlap) ) continue;
            int fullIndex = irow * vars->numValuesSplineBuffer + icol;
            // Read in new trace and store in samples
            // Re-use already allocated samples buffer
            while( vars->splineSampleBuffer[ counter ] == NULL ) {
              counter++;
            }
            sampleBufferTemp[fullIndex] = vars->splineSampleBuffer[counter]; // Copy over pointer to allocated sample array. Next, overwrite sample array with newly read-in values
            vars->splineSampleBuffer[counter++] = NULL;
            int traceIndex = getTraceIndex( vars, row, col );
            if( traceIndex < 0 || traceIndex >= vars->numTracesIn ) {
              if( traceIndex >= vars->numTracesIn ) {
                vars->traceCounter = vars->numTracesOut+1;
                return false;
              }
              traceIsOutsideArea = true;
              vars->p1Buffer.row = 0;
              vars->p1Buffer.col = 0;
              vars->p2Buffer.row = 0;
              vars->p2Buffer.col = 0;
              break;
            }
            //        fprintf(stderr,"Trace index: %d,   %d %d   ~=  %d %d\n", traceIndex, row, col, vars->rowValuesIn[traceIndex].intValue(), vars->colValuesIn[traceIndex].intValue() );
            bool success = vars->reader->moveToTrace( traceIndex );
            if( !success ) writer->error("Unknown error occurred trying to move file pointer to trace %d", traceIndex+1);
            //        success = vars->reader->readTrace( sampleBufferTemp[fullIndex], vars->hdrValueBlock_dummy, shdr->numSamples );
            success = vars->reader->readTrace( sampleBufferTemp[fullIndex], shdr->numSamples );
            if( !success ) writer->error("Unknown error occurred reading trace from input file");
          } // END for icol
        } // END for irow
        for( int index = 0; index < totalNumTracesBuffer; index++ ) {
          vars->splineSampleBuffer[index] = sampleBufferTemp[index];
        }
        delete [] sampleBufferTemp;
      } // END if !traceIsOutsideArea
    } // END else (One or more new traces need to read in from input data)

    vars->traceCounter += 1;
    //    fprintf(stderr,"Loop end trace %d, outside: %d\n", vars->traceCounter, traceIsOutsideArea);
  } while( traceIsOutsideArea && vars->traceCounter < vars->numTracesOut );

  if( traceIsOutsideArea ) return false;

  vars->reader->setTraceHeaderValues( trcHdr );
  trcHdr->setDoubleValue( vars->hdrId_rowOut, rowOut );
  trcHdr->setDoubleValue( vars->hdrId_colOut, colOut );

  if( vars->method == mod_input_interpolate::METHOD_SPLINE || vars->method == mod_input_interpolate::METHOD_SPLINE_MINMAX ) {
    float* values2D = new float[totalNumTracesBuffer];
    csSplineInterpolation2D spline2D( vars->numValuesSplineBuffer, (float)col1, (float)vars->incIn.col, vars->numValuesSplineBuffer, (float)row1, (float)vars->incIn.row );

    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      for( int index = 0; index < totalNumTracesBuffer; index++ ) {
        values2D[index] = vars->splineSampleBuffer[index][isamp];
      }
      spline2D.prepareAll( values2D );
      samplesOut[isamp] = spline2D.compute( (float)colOut, (float)rowOut );
    }

    delete [] values2D;
  }

  if( vars->method != mod_input_interpolate::METHOD_SPLINE ) {
    int indexRow1 = (int)( ( rowOut - row1 ) / vars->incIn.row );
    int indexCol1 = (int)( ( colOut - col1 ) / vars->incIn.col );

    if( vars->method == mod_input_interpolate::METHOD_NEAREST ) {
      int row = indexRow1 * vars->incIn.row + row1;
      int col = indexCol1 * vars->incIn.col + col1;
      if( (rowOut-row) > (row+vars->incIn.row-rowOut) ) row = row+vars->incIn.row;
      if( (colOut-col) > (col+vars->incIn.col-colOut) ) col = col+vars->incIn.col;
      int indexRow = (int)( ( row - row1 ) / vars->incIn.row );
      int indexCol = (int)( ( col - col1 ) / vars->incIn.col );
      int index = indexRow * vars->numValuesSplineBuffer + indexCol;
      float* ptr = vars->splineSampleBuffer[index];
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samplesOut[isamp] = ptr[isamp];
      }
    }
    else if( vars->method == mod_input_interpolate::METHOD_BILINEAR ) {
      int index1 = indexRow1 * vars->numValuesSplineBuffer + indexCol1;
      int index2 = indexRow1 * vars->numValuesSplineBuffer + indexCol1+1;
      int index3 = (indexRow1+1) * vars->numValuesSplineBuffer + indexCol1;
      int index4 = (indexRow1+1) * vars->numValuesSplineBuffer + indexCol1+1;
      float* ptr1 = vars->splineSampleBuffer[index1];
      float* ptr2 = vars->splineSampleBuffer[index2];
      float* ptr3 = vars->splineSampleBuffer[index3];
      float* ptr4 = vars->splineSampleBuffer[index4];

      double tmp0 = 1.0 / ( ( vars->p2Buffer.col - vars->p1Buffer.col )*( vars->p2Buffer.row - vars->p1Buffer.row ) );
      double tmp_row1 = vars->p2Buffer.row - rowOut;
      double tmp_row2 = rowOut - vars->p1Buffer.row;
      double tmp_col1 = vars->p2Buffer.col - colOut;
      double tmp_col2 = colOut - vars->p1Buffer.col;
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samplesOut[isamp] = tmp0 * ( tmp_row1 * ( ptr1[isamp] * tmp_col1 + ptr2[isamp] * tmp_col2 ) + tmp_row2 * ( ptr3[isamp] * tmp_col1 + ptr4[isamp] * tmp_col2 ) );
      }
    }
    else { // SPLINE min/max method
      int index1 = indexRow1 * vars->numValuesSplineBuffer + indexCol1;
      int index2 = indexRow1 * vars->numValuesSplineBuffer + indexCol1+1;
      int index3 = (indexRow1+1) * vars->numValuesSplineBuffer + indexCol1;
      int index4 = (indexRow1+1) * vars->numValuesSplineBuffer + indexCol1+1;
      float* ptr1 = vars->splineSampleBuffer[index1];
      float* ptr2 = vars->splineSampleBuffer[index2];
      float* ptr3 = vars->splineSampleBuffer[index3];
      float* ptr4 = vars->splineSampleBuffer[index4];

      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        float value = samplesOut[isamp];
        float minValue = std::min( ptr1[isamp], std::min( ptr2[isamp], std::min( ptr3[isamp], ptr4[isamp] ) ) );
        float maxValue = std::max( ptr1[isamp], std::max( ptr2[isamp], std::max( ptr3[isamp], ptr4[isamp] ) ) );
      
        if( value < minValue ) samplesOut[isamp] = minValue;
        else if( value > maxValue ) samplesOut[isamp] = maxValue;
      }
    }
  } // END not a pure spline interpolation

  return true;

}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_interpolate_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_INTERPOLATE", "Read/interpolate random line or regular grid from regular, rectangular 3d cube" );

  pdef->addParam( "filename", "Input file name. Supported file formats: '.cseis', '.rsf' and SEGY (all other extensions)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING );

  pdef->addParam( "method", "Interpolation method", NUM_VALUES_FIXED );
  pdef->addValue( "spline", VALTYPE_OPTION );
  pdef->addOption( "spline", "Standard 2D spline interpolation" );
  pdef->addOption( "spline_minmax", "Spline interpolation but enforce min/max: Check interpolated value against nearest four neighbors. Replace with min/max if necessary" );
  pdef->addOption( "nearest", "Nearest (four) neighbor interpolation" );
  pdef->addOption( "bilinear", "Bilinear interpolation of nearest four neighbors" );

  pdef->addParam( "define", "Method of model positions definition", NUM_VALUES_FIXED);
  pdef->addValue( "line", VALTYPE_OPTION );
  pdef->addOption( "line", "Read in regular 2D line defined by row/col numbers. Specify user parameters 'line_p1', 'line_p2', 'trc_inc'" );
  pdef->addOption( "grid", "Read in regular 3D grid defined by row/col numbers. Specify user parameters 'grid_il', 'grid_xl'" );
  //  pdef->addOption( "line_xy", "Read in regular 2D line defined by X/Y coordinates. Specify user parameters 'line_p1', 'line_p2', 'trc_inc'. Requires to specify binning grid (unless input data file is in Seaseis format and already has a grid defined)" );
  //  pdef->addOption( "trace", "Read in one trace per input file, then repeat until first file has been fully read in. Stop then." );

  pdef->addParam( "line_p1", "Start point of straight line", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Point #1 row/X coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Point #1 col/Y coordinate" );

  pdef->addParam( "line_p2", "End point of straight line", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Point #2 row/X coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Point #2 col/Y coordinate" );

  pdef->addParam( "line_binsize", "Bin/cell size in input data", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Bin/cell size in inline/row direction [m]    = distance between two columns", "" );
  pdef->addValue( "", VALTYPE_NUMBER, "Bin/cell size in crossline/col direction [m] = distance between two rows", "" );

  pdef->addParam( "trc_inc", "Output trace increment [m]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER );

  pdef->addParam( "hdrin_rowcol", "Input trace header names containing row & column numbers", NUM_VALUES_FIXED );
  pdef->addValue( "row", VALTYPE_STRING, "Inline/row trace header name" );
  pdef->addValue( "col", VALTYPE_STRING, "Crossline/column trace header name" );

  pdef->addParam( "hdrout_rowcol", "Output trace header names where row & column numbers are written", NUM_VALUES_FIXED, "These are double floating point headers" );
  pdef->addValue( "row_dbl", VALTYPE_STRING, "Inline/row trace header name" );
  pdef->addValue( "col_dbl", VALTYPE_STRING, "Crossline/column trace header name" );

  pdef->addParam( "grid_row", "Grid definition in row/inline direction", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "First inline" );
  pdef->addValue( "", VALTYPE_NUMBER, "Last  inline" );
  pdef->addValue( "", VALTYPE_NUMBER, "Inline increment" );

  pdef->addParam( "grid_col", "Grid definition in col/xline direction", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "First xline" );
  pdef->addValue( "", VALTYPE_NUMBER, "Last  xline" );
  pdef->addValue( "", VALTYPE_NUMBER, "Xline increment" );

  pdef->addParam( "grid_sort", "Sort order for grid definition", NUM_VALUES_FIXED);
  pdef->addValue( "row", VALTYPE_OPTION );
  pdef->addOption( "row", "Sort in row order" );
  pdef->addOption( "col", "Sort in column order" );

  //  pdef->addParam( "sample_int", "Specify to resample input data to a different sample interval", NUM_VALUES_FIXED, "Spline function will be used in the X-T domain to interpolate or decimate samples" );
  //  pdef->addValue( "", VALTYPE_NUMBER );

  /*
  pdef->addParam( "grid_azim", "Output grid azimuth", NUM_VALUES_VARIABLE, "Inline/crossline (row/col) directions, clock-wise from North" );
  pdef->addValue( "", VALTYPE_NUMBER, "Inline/row direction [deg]", "Direction of increasing crossline/col numbers" );
  pdef->addValue( "+90", VALTYPE_OPTION, "Crossline/col direction, relative to inline direction [deg]", "Direction of increasing inline/row numbers" );
  pdef->addOption( "+90", "Crossline direction is 90deg clock-wise from inline direction" );
  pdef->addOption( "-90", "Crossline direction is 90deg anticlock-wise from inline direction" );

  pdef->addParam( "grid_orig_xy", "Output grid origin XY coordinates", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Grid origin X coordinate [m]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Grid origin Y coordinate [m]" );

  pdef->addParam( "grid_orig_rowcol", "Output grid origin row/col number", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Origin inline/row number" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Origin crossline/col number" );

  pdef->addParam( "grid_binsize", "Output grid bin/cell size", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Bin/cell size in inline/row direction [m]", "" );
  pdef->addValue( "", VALTYPE_NUMBER, "Bin/cell size in crossline/col direction [m]", "" );

  pdef->addParam( "hdrin_binxy", "Trace headers containing X/Y coordinates", NUM_VALUES_FIXED);
  pdef->addValue( "bin_x", VALTYPE_STRING );
  pdef->addValue( "bin_y", VALTYPE_STRING );
  */
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_input_interpolate_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_input_interpolate::VariableStruct* vars = reinterpret_cast<mod_input_interpolate::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_input_interpolate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_input_interpolate::VariableStruct* vars = reinterpret_cast<mod_input_interpolate::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;

  if( vars->reader != NULL ) {
    delete vars->reader;
    vars->reader = NULL;
  }
  if( vars->splineSampleBuffer != NULL ) {
    int num = vars->numValuesSplineBuffer * vars->numValuesSplineBuffer;
    for( int i = 0; i < num; i++ ) {
      if( vars->splineSampleBuffer[i] != NULL ) delete [] vars->splineSampleBuffer[i];
    }
    delete [] vars->splineSampleBuffer;
    vars->splineSampleBuffer = NULL;
  }
  if( vars->colValuesIn != NULL ) {
    delete [] vars->colValuesIn;
    vars->colValuesIn = NULL;
  }
  if( vars->rowValuesIn != NULL ) {
    delete [] vars->rowValuesIn;
    vars->rowValuesIn = NULL;
  }

  if( vars->pointsToCompute != NULL ) {
    delete [] vars->pointsToCompute;
    vars->pointsToCompute = NULL;
  }
  delete vars; vars = NULL;
}


extern "C" void _params_mod_input_interpolate_( csParamDef* pdef ) {
  params_mod_input_interpolate_( pdef );
}
extern "C" void _init_mod_input_interpolate_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_interpolate_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_interpolate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_interpolate_( env, writer );
}
extern "C" void _exec_mod_input_interpolate_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_interpolate_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_interpolate_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_interpolate_( env, writer );
}

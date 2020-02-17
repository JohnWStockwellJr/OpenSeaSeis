/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_methods.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: BIN
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_bin {
  struct VariableStruct {
    int method;
    int hdrId_row;
    type_t hdrType_row;
    type_t hdrType_col;
    int hdrId_col;
    int hdrId_bin_x;
    int hdrId_bin_y;
  };
  static int const METHOD_DEFINE = 1;
  static int const METHOD_ROWCOL_TO_BINXY = 2;
  static int const METHOD_BINXY_TO_ROWCOL = 3;
}
using namespace mod_bin;

//*************************************************************************************************
// Init phase
//
//
//*************************************************************************************************
void init_mod_bin_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->method        = -1;
  vars->hdrId_bin_x   = -1;
  vars->hdrId_bin_y   = -1;
  vars->hdrId_row     = -1;
  vars->hdrId_col     = -1;
  vars->hdrType_row   = TYPE_INT;
  vars->hdrType_col   = TYPE_INT;

  //-------------------------------------------------------------
  //
  std::string text;
  if( param->exists("method") ) {
    param->getString("method",&text);
    if( !text.compare("define") ) {
      vars->method = METHOD_DEFINE;
    }
    else if( !text.compare("rowcol_to_binxy") ) {
      vars->method = METHOD_ROWCOL_TO_BINXY;
    }
    else if( !text.compare("binxy_to_rowcol") ) {
      vars->method = METHOD_BINXY_TO_ROWCOL;
    }
    else {
      writer->error("Unknown option: '%s'", text.c_str());
    }
  }

  if( param->exists("grid_orig_xy") ) {
    param->getDouble( "grid_orig_xy", &shdr->grid_orig_x, 0 );
    param->getDouble( "grid_orig_xy", &shdr->grid_orig_y, 1 );
  }
  else if( vars->method == METHOD_DEFINE ) {
    writer->line("Required user parameter not specified: 'grid_orig_xy'");
    env->addError();
  }

  if( param->exists("grid_orig_rowcol") ) {
    param->getInt("grid_orig_rowcol", &shdr->grid_orig_il, 0 );
    param->getInt("grid_orig_rowcol", &shdr->grid_orig_xl, 1 );
  }
  else if( vars->method == METHOD_DEFINE ) {
    shdr->grid_orig_il = 1;
    shdr->grid_orig_xl = 1;
    //    writer->line("Required user parameter not specified: 'grid_orig_rowcol'");
    //    env->addError();
  }

  if( param->exists("grid_binsize") ) {
    param->getDouble( "grid_binsize", &shdr->grid_binsize_il, 0 );
    param->getDouble( "grid_binsize", &shdr->grid_binsize_xl, 1 );
  }
  else if( vars->method == METHOD_DEFINE ) {
    writer->line("Required user parameter not specified: 'grid_binsize'");
    env->addError();
  }

  if( param->exists("grid_azim") ) {
    param->getDouble( "grid_azim", &shdr->grid_azim_il, 0 );
    if( param->getNumValues("grid_azim") > 1 ) {
      param->getString( "grid_azim", &text, 1 );
      if( !text.compare("+90") ) {
        shdr->grid_azim_xl = shdr->grid_azim_il + 90.0;
      }
      else if( !text.compare("-90") ) {
        shdr->grid_azim_xl = shdr->grid_azim_il - 90.0;
      }
      else {
        writer->line("Unknown option: '%s'", text.c_str());
        env->addError();
      }
    }
    else {
      shdr->grid_azim_xl = shdr->grid_azim_il + 90.0;
    }
  }
  else if( vars->method == METHOD_DEFINE ) {
    writer->line("Required user parameter not specified: 'grid_binsize'");
    env->addError();
  }

  if( vars->method != METHOD_DEFINE ) {
    if( shdr->grid_binsize_il == 0.0 || shdr->grid_binsize_xl == 0.0 ) {
      writer->error("Bin size found to be zero. Before binning can be performed, a survey grid has to be defined using parameter 'method', set to 'define'.");
    }
    string headerName_row("row");
    string headerName_col("col");
    string headerName_binx("bin_x");
    string headerName_biny("bin_y");

    if( param->exists("hdr_rowcol") ) {
      param->getString("hdr_rowcol", &headerName_row, 0);
      param->getString("hdr_rowcol", &headerName_col, 1);
    }
    if( param->exists("hdr_binxy") ) {
      param->getString("hdr_binxy", &headerName_binx, 0);
      param->getString("hdr_binxy", &headerName_biny, 1);
    }

    if( vars->method == METHOD_BINXY_TO_ROWCOL ) {
      if( !hdef->headerExists( headerName_binx ) ) {
        writer->error("Trace header '%s' does not exist.", headerName_binx.c_str());
      }
      else if( hdef->headerType( headerName_binx ) != TYPE_DOUBLE ) {
        writer->error("Trace header '%s' exists but has wrong type. Should be of type DOUBLE.", headerName_binx.c_str());
      }
      if( !hdef->headerExists( headerName_biny ) ) {
        writer->error("Trace header '%s' does not exist.", headerName_biny.c_str());
      }
      else if( hdef->headerType( headerName_biny ) != TYPE_DOUBLE ) {
        writer->error("Trace header '%s' exists but has wrong type. Should be of type DOUBLE.", headerName_biny.c_str());
      }
      if( !hdef->headerExists( headerName_row ) ) {
        hdef->addHeader(TYPE_INT, headerName_row,"Row number, survey grid");
      }
      if( !hdef->headerExists( headerName_col ) ) {
        hdef->addHeader(TYPE_INT,headerName_col,"Column number survey grid");
      }
    }
    else if( vars->method == METHOD_ROWCOL_TO_BINXY ) {
      if( !hdef->headerExists( headerName_row ) ) {
        writer->error("Trace header '%s' does not exist.", headerName_row.c_str());
      }
      if( !hdef->headerExists( headerName_col ) ) {
        writer->error("Trace header '%s' does not exist.", headerName_col.c_str());
      }
      if( !hdef->headerExists( headerName_binx ) ) {
        hdef->addHeader(TYPE_DOUBLE,headerName_binx,"Bin X coordinate [m]");
      }
      if( !hdef->headerExists( headerName_biny ) ) {
        hdef->addHeader(TYPE_DOUBLE,headerName_biny,"Bin Y coordinate [m]");
      }
    }
    vars->hdrId_row   = hdef->headerIndex(headerName_row);
    vars->hdrId_col   = hdef->headerIndex(headerName_col);
    vars->hdrType_row = hdef->headerType(headerName_row);
    vars->hdrType_col = hdef->headerType(headerName_col);
    vars->hdrId_bin_x = hdef->headerIndex(headerName_binx);
    vars->hdrId_bin_y = hdef->headerIndex(headerName_biny);
  }

  shdr->dump( writer->getFile() );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_bin_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace = traceGather->trace(0);


  csTraceHeader* trcHdr = trace->getTraceHeader();

  //  double an_rad = ( shdr->grid_azim - 90.0 ) * M_PI / 180.0;
  double an_il_rad = ( shdr->grid_azim_il ) * M_PI / 180.0;
  double an_xl_rad = ( shdr->grid_azim_xl ) * M_PI / 180.0;
  if( vars->method == METHOD_BINXY_TO_ROWCOL ) {
    double bin_x = trcHdr->doubleValue( vars->hdrId_bin_x );
    double bin_y = trcHdr->doubleValue( vars->hdrId_bin_y );
    double dx1 = bin_x - shdr->grid_orig_x;
    double dy1 = bin_y - shdr->grid_orig_y;
    double dx2 = cos(an_il_rad)*dx1 - sin(an_il_rad)*dy1;
    double dy2 = sin(an_il_rad)*dx1 + cos(an_il_rad)*dy1;
    double row = 0;
    double col = 0;

    row = dx2 / shdr->grid_binsize_xl * sin(an_xl_rad-an_il_rad) + (double)shdr->grid_orig_il;
    col = dy2 / shdr->grid_binsize_il                            + (double)shdr->grid_orig_xl;

    switch( vars->hdrType_row ) {
    case TYPE_INT:
      trcHdr->setIntValue( vars->hdrId_row, (int)round(row) );
      break;
    case TYPE_INT64:
      trcHdr->setInt64Value( vars->hdrId_row, (csInt64_t)round(row) );
      break;
    case TYPE_DOUBLE:
      trcHdr->setDoubleValue( vars->hdrId_row, row );
      break;
    case TYPE_FLOAT:
      trcHdr->setFloatValue( vars->hdrId_row, (float)row );
      break;
    }
    switch( vars->hdrType_col ) {
    case TYPE_INT:
      trcHdr->setIntValue( vars->hdrId_col, (int)round(col) );
      break;
    case TYPE_INT64:
      trcHdr->setInt64Value( vars->hdrId_col, (csInt64_t)round(col) );
      break;
    case TYPE_DOUBLE:
      trcHdr->setDoubleValue( vars->hdrId_col, col );
      break;
    case TYPE_FLOAT:
      trcHdr->setFloatValue( vars->hdrId_col, (float)col );
      break;
    }
  }
  else if( vars->method == METHOD_ROWCOL_TO_BINXY ) {
    double row = 0.0;
    double col = 0.0;
    switch( vars->hdrType_row ) {
    case TYPE_INT:
      row = (double)trcHdr->intValue( vars->hdrId_row );
      break;
    case TYPE_INT64:
      row = (double)trcHdr->int64Value( vars->hdrId_row );
      break;
    case TYPE_DOUBLE:
      row = trcHdr->doubleValue( vars->hdrId_row );
      break;
    case TYPE_FLOAT:
      row = (double)trcHdr->floatValue( vars->hdrId_row );
      break;
    default:
      writer->error("Program bug: Row number trace header has wrong type");
    }
    switch( vars->hdrType_col ) {
    case TYPE_INT:
      col = (double)trcHdr->intValue( vars->hdrId_col );
      break;
    case TYPE_INT64:
      col = (double)trcHdr->int64Value( vars->hdrId_col );
      break;
    case TYPE_DOUBLE:
      col = trcHdr->doubleValue( vars->hdrId_col );
      break;
    case TYPE_FLOAT:
      col = (double)trcHdr->floatValue( vars->hdrId_col );
      break;
    default:
      writer->error("Program bug: Col number trace header has wrong type");
    }

    double dx1 = ( row - (double)shdr->grid_orig_il ) * shdr->grid_binsize_xl / sin(an_xl_rad-an_il_rad);
    double dy1 = ( col - (double)shdr->grid_orig_xl ) * shdr->grid_binsize_il;
    double bin_x =  cos(an_il_rad)*dx1 + sin(an_il_rad)*dy1 + shdr->grid_orig_x;
    double bin_y = -sin(an_il_rad)*dx1 + cos(an_il_rad)*dy1 + shdr->grid_orig_y;
    //fprintf(stderr,"%d %d    %f %f  %f %f   %f\n",dx1,dy1,bin_x,bin_y,sin(an_xl_rad-an_il_rad)*180.0/M_PI);
    //    fprintf(stderr,"%f %f    %f\n", shdr->grid_azim_il, shdr->grid_azim_xl, sin(an_xl_rad-an_il_rad));

    trcHdr->setDoubleValue( vars->hdrId_bin_x, bin_x );
    trcHdr->setDoubleValue( vars->hdrId_bin_y, bin_y );
  }

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_bin_( csParamDef* pdef ) {
  pdef->setModule( "BIN", "Perform binning" );

  pdef->addDoc("Note that binning is performed to bin centers, i.e. bin_x & bin_y after binning are bin centers.");

  pdef->addParam( "method", "Binning method", NUM_VALUES_FIXED,
                  "If new grid is defined (option 'define'), all parameters named 'grid_...' must be specified. For the other options, any or all grid parameter may be redefined. Grid parameters that are not redefined are taken from the existing input super header" );
  pdef->addValue( "rowcol_to_binxy", VALTYPE_OPTION );
  pdef->addOption( "rowcol_to_binxy", "Set bin XY coordinates from row/col numbers" );
  pdef->addOption( "binxy_to_rowcol", "Set row/col numbers from bin XY coordinates" );
  pdef->addOption( "define", "Define survey grid", "Requires to specify all user parameters named 'grid_...'" );

  pdef->addParam( "grid_azim", "Grid azimuth", NUM_VALUES_VARIABLE, "Inline/crossline (row/col) directions, clock-wise from North" );
  pdef->addValue( "", VALTYPE_NUMBER, "Inline/row direction [deg]", "Direction of increasing crossline/col numbers" );
  pdef->addValue( "+90", VALTYPE_OPTION, "Crossline/col direction, relative to inline direction [deg]", "Direction of increasing inline/row numbers" );
  pdef->addOption( "+90", "Crossline direction is 90deg clock-wise from inline direction" );
  pdef->addOption( "-90", "Crossline direction is 90deg anticlock-wise from inline direction" );

  pdef->addParam( "grid_orig_xy", "Grid origin XY coordinates", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Grid origin X coordinate [m]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Grid origin Y coordinate [m]" );

  pdef->addParam( "grid_orig_rowcol", "Grid origin row/col number", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Origin inline/row number" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Origin crossline/col number" );

  pdef->addParam( "grid_binsize", "Grid bin/cell size", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Bin/cell size in inline/row direction [m]    = distance between crosslines", "" );
  pdef->addValue( "", VALTYPE_NUMBER, "Bin/cell size in crossline/col direction [m] = distance between inlines", "" );

  pdef->addParam( "hdr_rowcol", "Trace header names for row & column output", NUM_VALUES_FIXED );
  pdef->addValue( "row", VALTYPE_STRING, "Inline/row output trace header name" );
  pdef->addValue( "col", VALTYPE_STRING, "Crossline/column output trace header name" );

  pdef->addParam( "hdr_binxy", "Trace header names for bin X & Y output", NUM_VALUES_FIXED );
  pdef->addValue( "bin_x", VALTYPE_STRING, "Bin X trace header name" );
  pdef->addValue( "bin_y", VALTYPE_STRING, "Bin Y output trace header name" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_bin_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_bin::VariableStruct* vars = reinterpret_cast<mod_bin::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_bin_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_bin::VariableStruct* vars = reinterpret_cast<mod_bin::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_bin_( csParamDef* pdef ) {
  params_mod_bin_( pdef );
}
extern "C" void _init_mod_bin_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_bin_( param, env, writer );
}
extern "C" bool _start_exec_mod_bin_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_bin_( env, writer );
}
extern "C" void _exec_mod_bin_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_bin_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_bin_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_bin_( env, writer );
}

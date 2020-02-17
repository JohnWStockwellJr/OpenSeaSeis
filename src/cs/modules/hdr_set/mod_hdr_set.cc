/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csTableNew.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: HDR_SET
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_hdr_set {
  struct VariableStruct {
    csTableNew* table;
    int tableValueIndex; // Index of value column in table containing given header name that shall be set
    int* hdrId_keys;
    int hdrId;
    int hdrType;
  };
}
using mod_hdr_set::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_hdr_set_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


//---------------------------------------------
//
  vars->table = NULL;
  vars->hdrId_keys     = NULL;
  vars->hdrId   = -1;
  vars->hdrType = 0;
  vars->tableValueIndex = -1;

  //-----------------------------------------------
  std::string headerName;
  
  param->getString("header", &headerName );

  //-----------------------------------------------
  std::string text;
  bool sortTable = false;
  if( param->exists("table_sort" ) ) {
    param->getString("table_sort", &text );
    if( !text.compare("yes") ) {
      sortTable = true;
    }
    else if( !text.compare("no") ) {
      sortTable = false;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }
  //-----------------------------------------------

  vars->table = new csTableNew( csTableNew::TABLE_TYPE_UNIQUE_KEYS );

  //-----------------------------------------------
  int numKeys = param->getNumLines("key");
  if( numKeys == 0 ) {
    writer->error("No table key(s) specified.");
  }
  vars->hdrId_keys     = new int[numKeys];
  for( int ikey = 0; ikey < numKeys; ikey++ ) {
    std::string headerName;
    int col;
    bool interpolate = false;
    param->getStringAtLine( "key", &headerName, ikey, 0 );
    param->getIntAtLine( "key", &col, ikey, 1 );
    if( param->getNumValues( "key", ikey ) > 2 ) {
      param->getStringAtLine( "key", &text, ikey, 2 );
      if( !text.compare("yes") ) {
        interpolate = true;
      }
      else if( !text.compare("no") ) {
        interpolate = false;
      }
      else {
        writer->error("Unknown option: %s", text.c_str() );
      }
    }
    vars->table->addKey( col-1, interpolate );  // -1 to convert from 'user' column to 'C' column
    if( !hdef->headerExists( headerName ) ) {
      writer->error("No matching trace header found for table key '%s'", headerName.c_str() );
    }
    vars->hdrId_keys[ikey] = hdef->headerIndex( headerName );
  } // END for ikey

  //-----------------------------------------------
  int numValues = param->getNumLines("header");
  if( numValues == 0 ) {
    writer->error("No table 'value' header specified.");
  }
  if( numValues > 1 ) {
    writer->error("Only one table 'value' header is supported right now.");
  }
  for( int ival = 0; ival < numValues; ival++ ) {
    std::string headerName;
    int col;
    param->getStringAtLine("header",&headerName,ival,0);
    param->getIntAtLine("header",&col,ival,1);
    vars->table->addValue( col-1 );  // -1 to convert from 'user' column to 'C' column

    if( !hdef->headerExists(headerName) ) {
      writer->error("Unknown trace header '%s'", headerName.c_str());
    }
    vars->hdrId   = hdef->headerIndex(headerName);
    vars->hdrType = hdef->headerType(headerName);
    if( vars->hdrType == TYPE_STRING ) {
      writer->error("String headers are not supported by this module.");
    }
    vars->tableValueIndex = ival;
  } // END for ival

  //--------------------------------------------------
  param->getString("table", &text );
  try {
    vars->table->initialize( text, sortTable );
  }
  catch( csException& exc ) {
    writer->error("Error when initializing input table '%s': %s\n", text.c_str(), exc.getMessage() );
  }

  if( edef->isDebug() ) {
    for( int i = 0; i < vars->table->numLocations(); i++ ) {
      int station = (int)vars->table->getKeyValue( i, 0 );
      double value = vars->table->getValue( i, 0 );
      fprintf(stdout,"%d  %f\n", station, value);
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_hdr_set_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;

  csTrace* trace = traceGather->trace(0);


  csTraceHeader* trcHdr = trace->getTraceHeader();

  //  if( edef->isDebug() ) fprintf(stdout,"MOD_HDR_SET: Value: %f\n", value );
  double* keyValueBuffer = new double[vars->table->numKeys()];
  for( int ikey = 0; ikey < vars->table->numKeys(); ikey++ ) {
    keyValueBuffer[ikey] = trace->getTraceHeader()->doubleValue( vars->hdrId_keys[ikey] );
  }
  double value;
  try {
    value = vars->table->getValue( keyValueBuffer, vars->tableValueIndex );
  }
  catch( csException& e ) {
    delete [] keyValueBuffer;
    writer->error("Error occurred in HDR_SET: %s", e.getMessage());
    throw(e);
  }

  if( vars->hdrType == TYPE_INT ) {
    trcHdr->setIntValue( vars->hdrId, (int)value );
  }
  else if( vars->hdrType == TYPE_FLOAT ) {
    trcHdr->setFloatValue( vars->hdrId, (float)value );
  }
  else if( vars->hdrType == TYPE_DOUBLE ) {
    trcHdr->setDoubleValue( vars->hdrId, value );
  }
  else if( vars->hdrType == TYPE_INT64 ) {
    trcHdr->setInt64Value( vars->hdrId, (csInt64_t)value );
  }

  delete [] keyValueBuffer;

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_hdr_set_( csParamDef* pdef ) {
  pdef->setModule( "HDR_SET", "Set trace header from table", "Trace header values are read from table and interpolated if necessary" );

  pdef->addParam( "table", "Table containing trace header values", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Name of table containing trace header values");

  pdef->addParam( "key", "Key trace header used to match values found in specified table columns", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name of key header" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table" );
  pdef->addValue( "no", VALTYPE_OPTION, "Interpolate based to this key?" );
  pdef->addOption( "yes", "Use this key for interpolation of value" );
  pdef->addOption( "no", "Do not use this key for interpolation", "The input table is expected to contain the exact key values for this trace header" );

  pdef->addParam( "header", "Trace header to be read in/interpolated from specified table column", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number in input table" );

  pdef->addParam( "table_sort", "Sort table on key columns prior to interpolation?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Sort input table" );
  pdef->addOption( "no", "Do not sort table on input. Assume input table is sorted according to its key columns" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_hdr_set_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_hdr_set::VariableStruct* vars = reinterpret_cast<mod_hdr_set::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_hdr_set_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_hdr_set::VariableStruct* vars = reinterpret_cast<mod_hdr_set::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->table != NULL ) {
    delete vars->table;
    vars->table = NULL;
  }
  if( vars->hdrId_keys != NULL ) {
    delete [] vars->hdrId_keys;
    vars->hdrId_keys = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_hdr_set_( csParamDef* pdef ) {
  params_mod_hdr_set_( pdef );
}
extern "C" void _init_mod_hdr_set_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_hdr_set_( param, env, writer );
}
extern "C" bool _start_exec_mod_hdr_set_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_hdr_set_( env, writer );
}
extern "C" void _exec_mod_hdr_set_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_hdr_set_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_hdr_set_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_hdr_set_( env, writer );
}

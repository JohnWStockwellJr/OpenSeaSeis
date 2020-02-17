/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csP190Reader.h"
#include "csTime.h"
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: P190
 *
 * @author Bjorn Olofsson
 */
namespace mod_p190 {
  struct VariableStruct {
    cseis_io::csP190Reader* reader;

    int hdrID_chan;
    int hdrID_source;
    int hdrID_cable;

    int hdrID_recx;
    int hdrID_recy;
    int hdrID_recz;

    int hdrID_sou_index;
    int hdrID_soux;
    int hdrID_souy;
    int hdrID_sou_wdep;
    int hdrID_time_hour;
    int hdrID_time_min;
    int hdrID_time_sec;

    bool dropUnmatchedTraces;
  };
}
using namespace mod_p190;

//*************************************************************************************************
// Init phase
//
//
//*************************************************************************************************
void init_mod_p190_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->reader = NULL;

  vars->hdrID_chan = -1;
  vars->hdrID_recx = -1;
  vars->hdrID_recy = -1;
  vars->hdrID_recz = -1;
  vars->hdrID_cable = -1;

  vars->hdrID_source = -1;
  vars->hdrID_sou_index = -1;
  vars->hdrID_soux = -1;
  vars->hdrID_souy = -1;
  vars->hdrID_sou_wdep = -1;
  vars->hdrID_time_hour = -1;
  vars->hdrID_time_min = -1;
  vars->hdrID_time_sec = -1;

  vars->dropUnmatchedTraces = false;

  //----------------------------------------------------

  if( param->exists("drop_traces") ) {
    std::string text;
    param->getString( "drop_traces", &text );
    if( !text.compare("yes") ) {
      vars->dropUnmatchedTraces = true;
    }
    else if( !text.compare("no") ) {
      vars->dropUnmatchedTraces = false;
    }
    else {
      writer->error("Unknown option '%s'", text.c_str());
    }
  }

  std::string filename;

  param->getString( "filename", &filename );

  std::string text;
  if( param->exists("hdr_source") ) {
    param->getString( "hdr_source", &text );
    vars->hdrID_source = hdef->headerIndex( text.c_str() );
  }
  else {
    vars->hdrID_source = hdef->headerIndex( "source" );
  }
  if( param->exists("hdr_chan") ) {
    param->getString( "hdr_chan", &text );
    vars->hdrID_chan = hdef->headerIndex( text.c_str() );
  }
  else {
    vars->hdrID_chan = hdef->headerIndex( "chan" );
  }
  if( param->exists("hdr_cable") ) {
    param->getString( "hdr_cable", &text );
    vars->hdrID_cable = hdef->headerIndex( text.c_str() );
  }
  else {
    vars->hdrID_cable = hdef->headerIndex( "cable" );
  }


  try {
    vars->reader = new cseis_io::csP190Reader( filename );
    vars->reader->initialize();
  }
  catch( csException e ) {
    writer->error("Error occurred when opening file %s. System message:\n%s", filename.c_str(), e.getMessage() );
  }

  if( !hdef->headerExists(HDR_REC_X.name) ) {
    hdef->addStandardHeader( HDR_REC_X.name );
  }
  if( !hdef->headerExists(HDR_REC_Y.name) ) {
    hdef->addStandardHeader( HDR_REC_Y.name );
  }
  if( !hdef->headerExists(HDR_REC_Z.name) ) {
    hdef->addStandardHeader( HDR_REC_Z.name );
  }
  if( !hdef->headerExists(HDR_SOU_X.name) ) {
    hdef->addStandardHeader( HDR_SOU_X.name );
  }
  if( !hdef->headerExists(HDR_SOU_Y.name) ) {
    hdef->addStandardHeader( HDR_SOU_Y.name );
  }
  if( !hdef->headerExists(HDR_SOU_WDEP.name) ) {
    hdef->addStandardHeader( HDR_SOU_WDEP.name );
  }
  if( !hdef->headerExists(HDR_SOU_INDEX.name) ) {
    hdef->addStandardHeader( HDR_SOU_INDEX.name );
  }
  if( !hdef->headerExists(HDR_TIME_HOUR.name) ) {
    hdef->addStandardHeader( HDR_TIME_HOUR.name );
  }
  if( !hdef->headerExists(HDR_TIME_MIN.name) ) {
    hdef->addStandardHeader( HDR_TIME_MIN.name );
  }
  if( !hdef->headerExists(HDR_TIME_SEC.name) ) {
    hdef->addStandardHeader( HDR_TIME_SEC.name );
  }
  vars->hdrID_recx   = hdef->headerIndex( HDR_REC_X.name );
  vars->hdrID_recy   = hdef->headerIndex( HDR_REC_Y.name );
  vars->hdrID_recz   = hdef->headerIndex( HDR_REC_Z.name );

  vars->hdrID_sou_index = hdef->headerIndex( HDR_SOU_INDEX.name );
  vars->hdrID_soux      = hdef->headerIndex( HDR_SOU_X.name );
  vars->hdrID_souy      = hdef->headerIndex( HDR_SOU_Y.name );
  vars->hdrID_sou_wdep  = hdef->headerIndex( HDR_SOU_WDEP.name );

  vars->hdrID_time_hour = hdef->headerIndex( HDR_TIME_HOUR.name );
  vars->hdrID_time_min  = hdef->headerIndex( HDR_TIME_MIN.name );
  vars->hdrID_time_sec  = hdef->headerIndex( HDR_TIME_SEC.name );

  if( edef->isDebug() ) {
    vars->reader->dump();
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_p190_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csTrace* trace = traceGather->trace(0);


  csTraceHeader* trcHdr = trace->getTraceHeader();
  int source = trcHdr->intValue( vars->hdrID_source );
  int chan   = trcHdr->intValue( vars->hdrID_chan );
  int cable  = trcHdr->intValue( vars->hdrID_cable );

  cseis_io::csDataSource const* dataSource = vars->reader->getSource( source );
  if( dataSource == NULL ) {
    if( !vars->dropUnmatchedTraces ) {
      fprintf(stderr,"Error: Source not found in P190: %d\n", source);
      writer->error("Error: Source not found in P190: %d", source);
    }
    traceGather->freeAllTraces();
    return;
  }
  cseis_io::csDataChan const* dataChan = vars->reader->getChan( source, chan, cable );
  if( dataChan == NULL ) {
    if( !vars->dropUnmatchedTraces ) {
      fprintf(stderr,"Error: Channel/cable not found in P190: %d/%d\n", chan, cable);
      writer->error("Error: Channel/cable not found in P190: %d/%d", chan, cable);
    }
    traceGather->freeAllTraces();
    return;
  }

  trcHdr->setIntValue( vars->hdrID_sou_index, dataSource->id );
  trcHdr->setIntValue( vars->hdrID_time_hour, dataSource->hour );
  trcHdr->setIntValue( vars->hdrID_time_min, dataSource->min );
  trcHdr->setIntValue( vars->hdrID_time_sec, dataSource->sec );
  trcHdr->setDoubleValue( vars->hdrID_soux, dataSource->x );
  trcHdr->setDoubleValue( vars->hdrID_souy, dataSource->y );
  trcHdr->setDoubleValue( vars->hdrID_sou_wdep, dataSource->wdep );

  trcHdr->setDoubleValue( vars->hdrID_recx, dataChan->x );
  trcHdr->setDoubleValue( vars->hdrID_recy, dataChan->y );
  trcHdr->setDoubleValue( vars->hdrID_recz, dataChan->z );

  return;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_p190_( csParamDef* pdef ) {
  pdef->setModule( "P190", "Merge with P190 navigation data" );

  pdef->setVersion( 1, 0 );

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "hdr_source", "Trace header containing shot point number", NUM_VALUES_FIXED, "This must match the shot point number in the P1/90 file" );
  pdef->addValue( "source", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "hdr_chan", "Trace header containing channel number", NUM_VALUES_FIXED, "This must match the channel number in the P1/90 file" );
  pdef->addValue( "chan", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "hdr_cable", "Trace header containing cable number", NUM_VALUES_FIXED, "This must match the cable number in the P1/90 file" );
  pdef->addValue( "cable", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "drop_traces", "Drop traces which could not be found in P190 file?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Drop traces for which no match could be found in P190 file");
  pdef->addOption( "no", "Do not drop unmatched traces" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_p190_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_p190::VariableStruct* vars = reinterpret_cast<mod_p190::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_p190_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_p190::VariableStruct* vars = reinterpret_cast<mod_p190::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->reader ) {
    delete vars->reader;
    vars->reader = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_p190_( csParamDef* pdef ) {
  params_mod_p190_( pdef );
}
extern "C" void _init_mod_p190_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_p190_( param, env, writer );
}
extern "C" bool _start_exec_mod_p190_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_p190_( env, writer );
}
extern "C" void _exec_mod_p190_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_p190_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_p190_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_p190_( env, writer );
}

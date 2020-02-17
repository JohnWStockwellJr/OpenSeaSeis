/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSelectionManager.h"
#include "csSelection.h"
#include "csPoint2D.h"
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SELECT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_select {
  struct VariableStruct {
    csSelectionManager* selectionManager;
    bool isArea;
    int hdrID_x;
    int hdrID_y;
    std::string filename;
    cseis_geolib::csPoint2D* polygon;
    int numPolygonPoints;
  };
}
using mod_select::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_select_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


  vars->selectionManager = NULL;
  vars->hdrID_x = -1;
  vars->hdrID_y = -1;
  vars->isArea = false;
  vars->filename = "";
  vars->polygon = NULL;
  vars->numPolygonPoints = 0;

  //---------------------------------------------------------
  if( param->exists("select_xy_area") ) {
    std::string name1;
    std::string name2;
    param->getString("select_xy_area", &name1, 0);
    param->getString("select_xy_area", &name2, 1 );
    param->getString("select_xy_area", &vars->filename, 2 );
    vars->hdrID_x = hdef->headerIndex(name1);
    vars->hdrID_y = hdef->headerIndex(name2);
    vars->isArea = true;
  }

  //---------------------------------------------------------
  // Create new headers
  if( !vars->isArea ) {
    int nLines = param->getNumLines( "header" );
    if( nLines > 1 ) {
      writer->line( "More than one line encountered for user parameter '%s'. Only one line is supported.", "header" );
      env->addError();
    }

    csVector<std::string> valueList;
    param->getAll( "header", &valueList );

    if( valueList.size() == 0 ) {
      writer->line("Wrong number of arguments for user parameter '%s'. Expected: >0, found: %d.", "header", valueList.size());
      env->addError();
    }
    else {
      std::string text;
      param->getString( "select", &text );
      
      vars->selectionManager = NULL;
      try {
        vars->selectionManager = new csSelectionManager();
        vars->selectionManager->set( &valueList, &text, hdef );
      }
      catch( csException& e ) {
        vars->selectionManager = NULL;
        writer->error( "%s", e.getMessage() );
      }
      if( edef->isDebug() ) vars->selectionManager->dump();
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_select_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csTrace* trace = traceGather->trace(0);

  if( !vars->isArea ) {
    if( vars->selectionManager->contains( trace->getTraceHeader() ) ) {
      return;
    }
    else {
      traceGather->freeAllTraces();
      return;
    }
  }
  else if( vars->numPolygonPoints == 0 ) {
    FILE* f_polygon = fopen( vars->filename.c_str(), "r" );
    if( f_polygon == NULL ) {
      writer->error("Polygon ASCII (X Y) file not found: %s\n", vars->filename.c_str());
    }
    char buffer[256];
    while( fgets(buffer,256,f_polygon) != NULL ) {
      if( strlen(buffer) < 5 ) continue;
      vars->numPolygonPoints++;
    }
    vars->polygon = new cseis_geolib::csPoint2D[vars->numPolygonPoints];
    rewind(f_polygon);
    int counter = 0;
    while( fgets(buffer,256,f_polygon) != NULL ) {
      if( strlen(buffer) < 5 ) continue;
      csPoint2D p;
      sscanf(buffer,"%lf %lf", &p.x, &p.z );
      vars->polygon[counter++] = p;
    }
    fclose(f_polygon);
  }
  csTraceHeader* trcHdr = trace->getTraceHeader();
  cseis_geolib::csPoint2D p;
  p.x = trcHdr->doubleValue( vars->hdrID_x );
  p.z = trcHdr->doubleValue( vars->hdrID_y );
  if( !p.isPointInPolygon( vars->polygon, vars->numPolygonPoints ) ) {
    traceGather->freeAllTraces();
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_select_( csParamDef* pdef ) {
  pdef->setModule( "SELECT", "Select traces", "Select traces that match specified header selection" );

  pdef->addParam( "header", "Names of trace headers used for trace selection", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING,
    "List of selection strings, one for each specified header. See documentation for more detailed description of selection syntax",
    "Examples: '10,13' '10-13'-->'10,11,12,13' '10-30(10)'-->'10,20,30'. Operators: <,<=,>=,>,!" );
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_select_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_select::VariableStruct* vars = reinterpret_cast<mod_select::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_select_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_select::VariableStruct* vars = reinterpret_cast<mod_select::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->selectionManager ) {
    if( vars->selectionManager ) {
      delete vars->selectionManager; vars->selectionManager = NULL;
    }
    if( vars->polygon != NULL ) {
      delete [] vars->polygon;
      vars->polygon = NULL;
    }
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_select_( csParamDef* pdef ) {
  params_mod_select_( pdef );
}
extern "C" void _init_mod_select_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_select_( param, env, writer );
}
extern "C" bool _start_exec_mod_select_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_select_( env, writer );
}
extern "C" void _exec_mod_select_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_select_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_select_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_select_( env, writer );
}
